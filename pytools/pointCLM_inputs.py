#!/usr/bin/env python

#----------------------
# Extract point(s) data needed for running off-line ELM
#    (1) INPUT data source: global datasets
#    (2) 5 datasets: land domain, land surfdata, and optionally
#                    datm domain, metdata, initial data
#

import os, sys, time, math
import numpy as np
from optparse import OptionParser
import re
import netcdf_modules as nfmod

#-------------------Local functions --------------------------------------------
def pointLocator(lons_all, lats_all, lons_pt, lats_pt):
    
    xindxpts = []
    yindxpts = []
    
    
    if(lons_pt.strip()==''):
        lons_pt=[]
    else:                
        if ('str' in str(type(lons_pt))):#values with operators
            nondecimal = re.compile(r'[^\d.,:;-]+')
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
                    pt_val=lons_all[lons_all<max(v_pt)]
            
                if(">=" in lons_pt): 
                    pt_val=lons_all[lons_all>=min(v_pt)]
                elif(">" in lons_pt):
                    pt_val=lons_all[lons_all>min(v_pt)]
 
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
    if(lats_pt.strip()==''):
        lats_pt =[]
    else:
        if ('str' in str(type(lats_pt))):#values with operators
            nondecimal = re.compile(r'[^\d.,:;-]+')
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
                    pt_val=lats_all[lats_all<max(v_pt)]
            
                if(">=" in lats_pt): 
                    pt_val=lats_all[lats_all>=min(v_pt)]
                elif(">" in lats_pt):
                    pt_val=lats_all[lats_all>min(v_pt)]
 
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
        
                     
    if np.size(lons_pt)>0:  
        xind=np.where(np.isin(lons_all,lons_pt)) # paired indice
        xind=np.ravel_multi_index(xind,lons_all.shape) # flatten indice
    else:
        xind=[]
    if np.size(lats_pt)>0:
        yind=np.where(np.isin(lats_all,lats_pt))
        yind=np.ravel_multi_index(yind,lats_all.shape)
    else:
        yind=[]

    if(np.size(xind)>0 and np.size(yind)>0):
        # here we're picking points in 'xind' if it's in 'yind'. It shall be same if picking 'yind' if it's in 'xind'
        ij=np.where(np.isin(xind, yind))        
        ijindx=np.unravel_index(xind[ij], lons_all.shape)
    elif(np.size(xind)>0):
        ijindx=np.unravel_index(xind, lons_all.shape)
    elif(np.size(yind)>0):
        ijindx=np.unravel_index(yind, lats_all.shape)
    xstart = ijindx[0]
    ystart = ijindx[1]
   
    #unique and countinuing indices count
    x= np.array(xstart,dtype=int)
    if(np.size(x)>1):
        x=np.unique(x)
        xdiff=np.insert(np.diff(x),0,-1)
        xi=x[np.where(xdiff!=1)] # unique and non-continuing index
        xloc=np.where(np.isin(x,xi))[0] # index of xi in x, from which continuing indices can be counted
        xloc=np.append(xloc,np.size(x))
        xpts=np.diff(xloc)       
    elif(np.size(x)==1):
        xi=np.copy(x)
        xpts=np.array([1],dtype=int)
    xindxpts = np.column_stack((xi,xpts))
     
    # 
    y= np.array(ystart,dtype=int)
    if(np.size(y)>1):
        y=np.unique(y)
        ydiff=np.insert(np.diff(y),0,-1)
        yi=y[np.where(ydiff!=1)] # unique and non-continuing index
        yloc=np.where(np.isin(y,yi))[0] # index of yi in y, from which continuing indices can be counted
        yloc=np.append(yloc,np.size(y))
        ypts=np.diff(yloc)       
    elif(np.size(y)==1):
        yi=np.copy(y)
        ypts=np.array([1],dtype=int)
    yindxpts = np.column_stack((yi,ypts))
    
    return xindxpts, yindxpts



#-------------------Parse options-----------------------------------------------
parser = OptionParser()

parser.add_option("--e3sminput", dest="ccsm_input", \
                  default='./e3sm_inputdata', \
                  help = "input data directory for E3SM (required)")
parser.add_option("--caserundir", dest="casedataroot", default="", \
                  help = "component set to use")
parser.add_option("--finidat", dest="finidat", default='', \
                  help = "initial data file to use" \
                  +" (should be in your run directory)")
parser.add_option("--met_domain", dest="metdomain", default="", \
                  help = 'datm domain file for met data forcing')
parser.add_option("--met_dir", dest="metdir", default="", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--met_fileheader", dest="metfileheader", default="", \
                  help = 'headers of files for met data forcing')
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
if(not (options.shared_domain == '' or options.shared_domain == 'None')):
    if(options.shared_domain.startswith('/') or \
       options.shared_domain.startswith('./')):
        shareddomain = os.path.abspath(options.shared_domain)   #user-defined file with full path
    else: # file with the default path
        shareddomain = os.path.abspath(ccsm_input+"/share/domains/domain.clm/"+options.shared_domain)
    
    shareddomainname = os.path.basename(shareddomain)
    shareddomainname = os.path.splitext(shareddomainname)[0]


#------------met input data and domain directory ------------------
if(options.metdir.startswith('/') or \
    options.metdir.startswith('./')):
    metdir = os.path.abspath(options.metdir)
else:
    metdir   = os.path.abspath(ccsm_input+"/atm/datm7/"+options.metdir)

if(options.metdomain.startswith('/') or \
   options.metdomain.startswith('./')):
    metdomain = os.path.abspath(options.metdomain)
else:
    metdomain = os.path.abspath(ccsm_input+"/atm/datm7/"+options.metdomain)
    if(not os.path.isfile(metdomain) and os.path.isdir(metdir)):
        metdomain = os.path.abspath(metdir+"/"+options.metdomain)
metdomainname = os.path.basename(metdomain)
metdomainname = os.path.splitext(metdomainname)[0]


#------------surfdata input data directory/file ------------------
if(not (options.surfdata == '' or options.surfdata == 'none')):
    if(options.surfdata.startswith('/') or \
       options.surfdata.startswith('./')):
        surfdata = os.path.abspath(options.surfdata)   #user-defined file with full path
    else: # file with the default path
        surfdata = os.path.abspath(ccsm_input+"/lnd/clm2/surfdata_map/"+options.surfdata)
    surfname = os.path.basename(surfdata)
    surfname = os.path.splitext(surfname)[0]

###########################################################################
if(options.outdir.startswith('/')):
    outdir = os.path.abspath(options.outdir)
    if(not outdir.endswith('/')): outdir = outdir+'/'
else:
    outdir = './'

#------------new shared domain for land surfdata and datm -------------------------------------
if (not options.shared_domain.strip()==''):    
    shareddomain_new = outdir + shareddomainname + '_'+ \
               options.newdata_affix + '.nc'

    print('INFO: extracting shared_domain for land surfdata - \n', shareddomain, '\n -->', shareddomain_new)

    #
    [alllats, vdim, vattr]=nfmod.getvar(shareddomain,['yc'])
    alllats = list(alllats.values())[0] # dict --> list
    [alllons, vdim, vattr]=nfmod.getvar(shareddomain,['xc'])
    alllons = list(alllons.values())[0]
    
    #indx/indy: 2-D array with paired [startindex, numpts]
    [indx,indy] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)
            
    for pt in range(len(indx)):
        pts_idx=[indx[pt][0],indy[pt][0]]
        pts_num=[indx[pt][1],indy[pt][1]]
        nfmod.nco_extract(shareddomain, shareddomain_new, \
            ['nj','ni'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')

#------------new surfdata  -------------------------------------------------------------------
if (not (options.surfdata.strip()=='' or options.surfdata=='none')):
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

#------------new metdata domain -------------------------------------------------------------------
if (not options.metdomain.strip()==''):
    metdomain_new = outdir + metdomainname + '_'+ \
               options.newdata_affix + '.nc'

    print('INFO: extracting metdata domain - \n', metdomain, '\n -->', metdomain_new)

    #
    [alllats, vdim, vattr]=nfmod.getvar(metdomain,['yc'])
    alllats = list(alllats.values())[0] # dict --> list
    [alllons, vdim, vattr]=nfmod.getvar(metdomain,['xc'])
    alllons = list(alllons.values())[0]
    
    #indx/indy: 2-D array with paired [startindex, numpts]
    [indx,indy] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)
            
    for pt in range(len(indx)):
        pts_idx=[indx[pt][0],indy[pt][0]]
        pts_num=[indx[pt][1],indy[pt][1]]
        nfmod.nco_extract(metdomain, metdomain_new, \
            ['nj','ni'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')

#------------new metdata -------------------------------------------------------------------
if (not options.metdir.strip()==''):

    metdir_new = outdir+ '/atm/datm7/newmetdata_'+ \
               options.newdata_affix
    os.system('mkdir -p ' + metdir_new)

    #checking where is the original data
    dirfiles = os.listdir(metdir)
    for dirfile in dirfiles:        
        filehead = options.metfileheader
        filehead_new = filehead+'_'+options.newdata_affix
        indx=[]; indy=[]
        # in metdata directory
        if(os.path.isfile(metdir+'/'+dirfile)): # file names
            if(filehead in dirfile):
            
                print('\n file: '+dirfile)
    
                metfile_old = metdir+'/'+dirfile
                dirfile_new = dirfile.replace(filehead,filehead_new)
                metfile_new = metdir_new+'/'+ dirfile_new

                print('INFO: extracting met-data - \n', metfile_old, '\n -->', metfile_new)
                #
                if(len(indx)<=0):#only needed from the first metfile
                    [alllats, vdim, vattr]=nfmod.getvar(metfile_old,['LATIXY'])
                    alllats = list(alllats.values())[0] # dict --> list
                    [alllons, vdim, vattr]=nfmod.getvar(metfile_old,['LONGXY'])
                    alllons = list(alllons.values())[0]    
                    #indx/indy: 2-D array with paired [startindex, numpts]
                    [indx,indy] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)            
                
                for pt in range(len(indx)):
                    pts_idx=[indx[pt][0],indy[pt][0]]
                    pts_num=[indx[pt][1],indy[pt][1]]
                    nfmod.nco_extract(metfile_old, metfile_new, \
                            ['lat','lon'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')

        # in subdirectory of metdata directory
        elif(os.path.isdir(metdir+'/'+dirfile)): # subdirectories
            subfiles = os.listdir(metdir+'/'+dirfile)
            for subfile in subfiles:
    
                if(os.path.isfile(metdir+'/'+dirfile+'/'+subfile)):
                    if(filehead in subfile):
     
                        metfile_old = metdir+'/'+dirfile+'/'+subfile                        
                        subfile_new = subfile.replace(filehead,filehead_new)
                        metfile_new = metdir_new+'/'+subfile_new
                        
                        print('INFO: extracting met-data - \n', metfile_old, '\n -->', metfile_new)
                        #
                        if(len(indx)<=0):#only needed from the first metfile
                            [alllats, vdim, vattr]=nfmod.getvar(metfile_old,['LATIXY'])
                            alllats = list(alllats.values())[0] # dict --> list
                            [alllons, vdim, vattr]=nfmod.getvar(metfile_old,['LONGXY'])
                            alllons = list(alllons.values())[0]    
                            #indx/indy: 2-D array with paired [startindex, numpts]
                            [indx,indy] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)            
                        
                        for pt in range(len(indx)):
                            pts_idx=[indx[pt][0],indy[pt][0]]
                            pts_num=[indx[pt][1],indy[pt][1]]
                            nfmod.nco_extract(metfile_old, metfile_new, \
                                    ['lat','lon'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')

#------------END of pointCLM_data -----------------------------------------------------------------------
    

