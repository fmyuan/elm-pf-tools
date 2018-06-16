#!/usr/bin/env python

#----------------------
# Extract point(s) data needed for running off-line ELM
#    (1) INPUT data source: global datasets
#    (2) 5 datasets: land domain, land surfdata, and optionally
#                    datm domain, metdata, initial data
#

import netcdf_modules as nfmod
import os, sys, time, math
import numpy as np
from optparse import OptionParser
import re

#-------------------Local functions --------------------------------------------
def pointLocator(lons_all, lats_all, lons_pt, lats_pt):
    
    xindxpts = []
    yindxpts = []
    
    if ('str' in str(type(lons_pt))):#values with operators
        nondecimal = re.compile(r'[^\d.,:;]+')
        v_pt = nondecimal.sub("",lons_pt)
        v_pt = np.float64(v_pt.split("[,:;]"))
        if("~" in lons_pt):#nearest point(s)
            pt_val=[]
            for v in v_pt:
                pt_nearest=np.argmin(abs(lons_all-v))
                pt_nearest=np.unravel_index(pt_nearest,lons_all.shape)
                pt_val=np.append(pt_val,lons_all[pt_nearest])
        else:
            if("<=" in lons_pt): 
                pt_val=lons_all[lons_all<=max(v_pt)]
            elif("<" in lons_pt):
                pt_val=lons_all[lons_all<=max(v_pt)]
            
            if(">=" in lons_pt): 
                pt_val=pt_val[pt_val>=min(v_pt)]
            elif(">" in lons_pt):
                pt_val=pt_val[pt_val>min(v_pt)]
 
        if pt_val is None:
            print('Error:  invalid ranged value expression - '+lons_pt)
            sys.exit()
        else:
            lons_pt = np.copy(pt_val)

    else: # exactly point(s)

        if(lons_all.dtype=='float64'):
            lons_pt = np.float64(lons_pt)
        else:
            lons_pt = np.float(lons_pt)
        
    #   
    if ('str' in str(type(lats_pt))):#values with operators
        nondecimal = re.compile(r'[^\d.,:;]+')
        v_pt = nondecimal.sub("",lats_pt)
        v_pt = np.float64(v_pt.split("[,:;]"))
        if("~" in lats_pt):#nearest point(s)
            pt_val=[]
            for v in v_pt:
                pt_nearest=np.argmin(abs(lats_all-v))
                pt_nearest=np.unravel_index(pt_nearest,lats_all.shape)
                pt_val=np.append(pt_val,lats_all[pt_nearest])
        else:
            if("<=" in lats_pt): 
                pt_val=lats_all[lats_all<=max(v_pt)]
            elif("<" in lats_pt):
                pt_val=lats_all[lats_all<=max(v_pt)]
            
            if(">=" in lats_pt): 
                pt_val=pt_val[pt_val>=min(v_pt)]
            elif(">" in lats_pt):
                pt_val=pt_val[pt_val>min(v_pt)]
 
        if pt_val is None:
            print('Error:  invalid ranged value expression - '+lats_pt)
            sys.exit()
        else:
            lats_pt = np.copy(pt_val)

    else: # exactly point(s)
        if(lats_all.dtype=='float64'):
            lats_pt = np.float64(lats_pt)
        else:
            lats_pt = np.float(lats_pt)
        
                     
    # joined indices
    if np.size(lons_pt)>0:  
        xind=np.where(np.isin(lons_all,lons_pt))
    else:
        xind=[]
    if np.size(lats_pt)>0:
        yind=np.where(np.isin(lats_all,lats_pt))
    else:
        yind=[]

    if(np.size(xind)>0 and np.size(yind)>0):
        ij=np.isin(np.isin(xind[0],yind[0]),np.isin(xind[1],yind[1]))
        xstart = xind[0][ij]
        ystart = xind[1][ij]
    elif(np.size(xind)>0):
        xstart = xind[0]
        ystart = xind[1]        
    elif(np.size(yind)>0):
        xstart = yind[0]
        ystart = yind[1]
   
    #unique and countinuing indices count
    x= np.array(xstart,dtype=int)
    if(np.size(x)>1):
        x=np.unique(x)
        xdiff=np.insert(np.diff(x),0,-1)
        xi=x[np.where(xdiff!=1)] # unique and non-continuing index
        loc=np.where(np.isin(x,xi))[0] # index of xi in x, from which continuing indices can be counted
        loc=np.append(loc,np.size(x))
        pts=np.diff(loc)       
    elif(np.size(x)==1):
        xi=np.copy(x)
        pts=np.array([1],dtype=int)
    xindxpts = np.column_stack((xi,pts))
     
    # 
    x= np.array(ystart,dtype=int)
    if(np.size(x)>1):
        x=np.unique(x)
        xdiff=np.insert(np.diff(x),0,-1)
        xi=x[np.where(xdiff!=1)] # unique and non-continuing index
        loc=np.where(np.isin(x,xi))[0] # index of xi in x, from which continuing indices can be counted
        loc=np.append(loc,np.size(x))
        pts=np.diff(loc)       
    elif(np.size(x)==1):
        xi=np.copy(x)
        pts=np.array([1],dtype=int)
    yindxpts = np.column_stack((xi,pts))
    
    return xindxpts, yindxpts



#-------------------Parse options-----------------------------------------------
parser = OptionParser()

parser.add_option("--e3sminput", dest="ccsm_input", \
                  default='./e3sm_inputdata', \
                  help = "input data directory for E3SM (required)")
parser.add_option("--caserundir", dest="casedataroot", default="", \
                  help = "component set to use (required)")
parser.add_option("--finidat", dest="finidat", default='', \
                  help = "initial data file to use" \
                  +" (should be in your run directory)")
parser.add_option("--metdomain", dest="metdomain", default="", \
                  help = 'datm domain file for met data forcing')
parser.add_option("--metdir", dest="metdir", default="", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--shared_domain", dest="shared_domain", default="", \
                  help = 'land domain file either in fullpath or filename in default path')
parser.add_option("--surfdata", dest="surfdata", default="none", \
                  help = 'land surfdata file either in fullpath or filename in default path')
parser.add_option("--outdir", dest="outdir", default="", \
                  help = 'new dataset output directory')
parser.add_option("--newdata_affix", dest="newdata_affix", default="", \
                  help = 'land surfdata file either in fullpath or filename in default path')
parser.add_option("--lons", dest="ptlon", default="", \
                      help = "point or ranged(<,<=,>=,>,~) longitude(s), separated with',:;'")
parser.add_option("--lats", dest="ptlat", default="", \
                      help = "point or ranged (<,<=,>=,>,~) latitude(s), separated with',:;'")
parser.add_option("--ncobinpath", dest="ncobinpath", default="", \
                      help = "NCO bin path if not in $PATH")


(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------
currdir = os.getcwd()

#--------------- nco bin path -----------
if(not options.ncobinpath.strip()==''):
    ncopath = options.ncobinpath

#---------------case information -----------------
if(not options.casedataroot.strip()==''):
    casedatadir = os.path.abspath(options.casedataroot)

#------------ccsm input data directory ------------------
print(options.ccsm_input)
if (options.ccsm_input == '' or (os.path.exists(options.ccsm_input) \
                                 == False)):
    print('Error:  invalid input data directory')
    sys.exit()
else:
    ccsm_input = os.path.abspath(options.ccsm_input)

#-----------shared land/datm domain directory/file ------------------
if (options.shared_domain == '' ):
    print('Error:  invalid input data directory')
#    sys.exit()
else:
    ccsm_input = os.path.abspath(options.ccsm_input)


#------------met input data and domain directory ------------------
if(options.metdir.startswith('/') or \
    options.metdir.startswith('./')):
    met_dir = os.path.abspath(options.metdir)
else:
    met_dir   = os.path.abspath(ccsm_input+"/atm/datm7/"+options.metdir)

if(options.metdomain.startswith('/') or \
   options.metdomain.startswith('./')):
    metdomain = os.path.abspath(options.metdomain)
else:
    metdomain = os.path.abspath(ccsm_input+"/atm/datm7/"+options.metdomain)

#------------surfdata input data directory/file ------------------
if(options.surfdata.startswith('/') or \
   options.surfdata.startswith('./')):
    surfdata = os.path.abspath(options.surfdata)
else:
    surfdata = os.path.abspath(ccsm_input+"/lnd/clm2/surfdata_map/"+options.surfdata)
surfname = os.path.basename(surfdata)
surfname = os.path.splitext(surfname)[0]

###########################################################################
if(options.outdir.startswith('/')):
    outdir = os.path.abspath(options.outdir)
    if(not outdir.endswith('/')): outdir = outdir+'/'
else:
    outdir = './'

#------------new surfdata  ------------------
if (not options.surfdata.strip()==''):    
    surfdata_new = outdir + surfname + '_'+ \
               options.newdata_affix + '.nc'

    print('INFO: extracting surfdata - \n', surfdata, '\n -->', surfdata_new)

    #
    [alllats, vdim, vattr]=nfmod.getvar(surfdata,['LATIXY'])
    alllats = list(alllats.values())[0] # dict --> list
    [alllons, vdim, vattr]=nfmod.getvar(surfdata,['LONGXY'])
    alllons = list(alllons.values())[0]
    
    #indx/indy: 2-D array with paired [startindex, numpts]
    [indx,indy] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)
            
    for pt in range(len(indx)):
        pts_idx=[indx[pt][0],indy[pt][0]]
        pts_num=[indx[pt][1],indy[pt][1]]
        nfmod.nco_extract(surfdata, surfdata_new, \
            ['lsmlat','lsmlon'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')

