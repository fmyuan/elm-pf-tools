#!/usr/bin/env python

#----------------------
# Extract point(s) data from ELM outputs
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
        if(len(xind)>1):
            ystart = xind[1]
    elif(np.size(yind)>0):
        xstart = yind[0]
        if(len(yind)>1):
            ystart = yind[1]
   
    #unique and continuing indices count
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
    if():
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
    else:
        yindxpts = []
    
    return xindxpts, yindxpts



#-------------------Parse options-----------------------------------------------
parser = OptionParser()

parser.add_option("--outdir", dest="outdir", default="", \
                  help = 'output directory')
parser.add_option("--ofheader", dest="header", default="none", \
                  help = 'header of output files either in fullpath or filename in default path')
parser.add_option("--newdata_affix", dest="newdata_affix", default="", \
                  help = 'new dataset output with affix')
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



###########################################################################
if(options.outdir.startswith('/')):
    outdir = os.path.abspath(options.outdir)
    if(not outdir.endswith('/')): outdir = outdir+'/'
else:
    outdir = './'


#------------extract all data as needed -----------------------------------
if (os.path.isdir(outdir.strip())):

    #checking where is the original data
    dirfiles = os.listdir(outdir)
    for dirfile in dirfiles:        
        filehead = options.header
        filehead_new = filehead+'_'+options.newdata_affix
        indx1=[]; indy1=[]
        indx2=[]; indy2=[]
        indx3=[]; indy3=[]
        indx4=[]; indy4=[]
        
        # in outdir directory
        if(os.path.isfile(outdir+'/'+dirfile)): # file names
            if(filehead in dirfile):
            
                print('\n file: '+dirfile)
    
                outfile_old = outdir+'/'+dirfile
                dirfile_new = dirfile.replace(filehead,filehead_new)
                outfile_new = outdir+'/'+dirfile_new

                print('INFO: extracting clm output - \n', outfile_old, '\n -->', outfile_new)
                #

                outfile_temp = outdir+'/temp.nc'
                # (1) for pft-level data: Doing this first may be faster (because the generated 'temp.nc' would be for rest ncks operation)
                  
                if(len(indx1)<=0):#only needed from the first output file                   
                    [alllats, vdim, vattr]=nfmod.getvar(outfile_old,['pfts1d_lat'])
                    alllats = list(alllats.values())[0] # dict --> list
                    [alllons, vdim, vattr]=nfmod.getvar(outfile_old,['pfts1d_lon'])
                    alllons = list(alllons.values())[0]    
                    #indx/indy: 2-D array with paired [startindex, numpts]
                    [indx1,indy1] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)                
                for pt in range(len(indx1)):
                    pts_idx=[indx1[pt][0]]
                    pts_num=[indx1[pt][1]]
                    nfmod.nco_extract(outfile_old,outfile_new, \
                            ['pft'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')

                # (2) for column-level data                   
                os.system('mv -f ' + outfile_new + ' '+ outfile_temp)   # do the extraction continuously

                if(len(indx2)<=0):#only needed from the first output file                   
                    [alllats, vdim, vattr]=nfmod.getvar(outfile_temp,['cols1d_lat'])
                    alllats = list(alllats.values())[0] # dict --> list
                    [alllons, vdim, vattr]=nfmod.getvar(outfile_temp,['cols1d_lon'])
                    alllons = list(alllons.values())[0]    
                    #indx/indy: 2-D array with paired [startindex, numpts]
                    [indx2,indy2] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)                
                for pt in range(len(indx2)):
                    pts_idx=[indx2[pt][0]]
                    pts_num=[indx2[pt][1]]
                    nfmod.nco_extract(outfile_temp,outfile_new, \
                            ['column'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')

                # (3) for landunit-level data                   
                os.system('mv -f ' + outfile_new + ' '+ outfile_temp)   # do the extraction continuously

                if(len(indx3)<=0):#only needed from the first output file                   
                    [alllats, vdim, vattr]=nfmod.getvar(outfile_temp,['land1d_lat'])
                    alllats = list(alllats.values())[0] # dict --> list
                    [alllons, vdim, vattr]=nfmod.getvar(outfile_temp,['land1d_lon'])
                    alllons = list(alllons.values())[0]    
                    #indx/indy: 2-D array with paired [startindex, numpts]
                    [indx3,indy3] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)                
                for pt in range(len(indx3)):
                    pts_idx=[indx3[pt][0]]
                    pts_num=[indx3[pt][1]]
                    nfmod.nco_extract(outfile_temp,outfile_new, \
                            ['landunit'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')


                # (4) for gridcell-level data                   
                os.system('mv -f ' + outfile_new + ' '+ outfile_temp)   # do the extraction continuously

                if(len(indx4)<=0):#only needed from the first output file                   
                    [alllats, vdim, vattr]=nfmod.getvar(outfile_temp,['grid1d_lat'])
                    alllats = list(alllats.values())[0] # dict --> list
                    [alllons, vdim, vattr]=nfmod.getvar(outfile_temp,['grid1d_lon'])
                    alllons = list(alllons.values())[0]    
                    #indx/indy: 2-D array with paired [startindex, numpts]
                    [indx4,indy4] = pointLocator(alllons, alllats, options.ptlon, options.ptlat)                
                for pt in range(len(indx4)):
                    pts_idx=[indx4[pt][0]]
                    pts_num=[indx4[pt][1]]
                    nfmod.nco_extract(outfile_temp,outfile_new, \
                            ['gridcell'], pts_idx, pts_num,ncksdir='/usr/local/nco/bin')

                os.system("rm -f " + outfile_temp)
            
#------------END of pointCLM_data -----------------------------------------------------------------------
    

