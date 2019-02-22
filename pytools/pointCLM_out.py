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
def pointLocator(lons_all, lats_all, lons_pt, lats_pt, is1D):
    
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

    #
    xstart = []
    if(len(xind)>0): xstart = xind[0]
    ystart = []
    if(len(yind)>0): ystart = yind[0]

    if(is1D):
        if(np.size(xind)>0 and np.size(yind)>0):
            ij=np.isin(np.isin(xind[0],yind[0]),np.isin(xind[1],yind[1]))
            xstart = xind[0][ij]
            ystart = xind[1][ij]
        
        elif(np.size(xind)>0):
            xstart = xind[0]
            if(len(xind)>1):
                ystart = xind[1]
    
        elif(np.size(yind)>0):
            ystart = yind[0]
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
    if(len(x)>0): xindxpts = np.column_stack((xi,pts))
     
    # 
    y= np.array(ystart,dtype=int)
    if(np.size(y)>1):
        y=np.unique(y)
        ydiff=np.insert(np.diff(y),0,-1)
        yi=y[np.where(ydiff!=1)] # unique and non-continuing index
        loc=np.where(np.isin(y,yi))[0] # index of xi in x, from which continuing indices can be counted
        loc=np.append(loc,np.size(y))
        pts=np.diff(loc)       
    elif(np.size(y)==1):
        yi=np.copy(y)
        pts=np.array([1],dtype=int)
    if(len(y)>0): yindxpts = np.column_stack((yi,pts))
    
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
else:
    ncopath = '/usr/local/nco/bin'

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
                os.system('cp -f ' + outfile_old + ' '+ outfile_temp)

                #loop through possible longitude/latitude pairs, with their dimension(s)
                coord_vars = [['lon','lat'], 
                               ['pfts1d_lon','pfts1d_lat'],
                               ['cols1d_lon','cols1d_lat'],
                               ['land1d_lon','land1d_lat'],
                               ['grid1d_lon','grid1d_lat'] ];
                coord_dims = [['lon','lat'], 
                               ['pft'],
                               ['column'],
                               ['landunit'],
                               ['gridcell'] ];
                
 
                # loop for all dimensions possibly in the nc file
                for idim in range(len(coord_dims)):
                    indx = []; indy =[];

                    dim = coord_dims[idim]
                    coord = coord_vars[idim]
                    is1D = True
                    if(len(dim)>1): is1D = False
                    
                    # all lon/lat pairs
                    [alllons, vdim, vattr]=nfmod.getvar(outfile_temp,[coord[0]])
                    if(len(alllons)>0):
                        alllons = list(alllons.values())[0] # dict --> list
                        
                        [alllats, vdim, vattr]=nfmod.getvar(outfile_temp,[coord[1]])
                        if(len(alllats)>0):    
                            alllats = list(alllats.values())[0]    
                            #indx/indy: 2-D array with paired [startindex, numpts]
                            [indx,indy] = pointLocator(alllons, alllats, options.ptlon, options.ptlat, is1D)                
                
                    vdim = dim[0]
                    if(len(indx)>0):
                        for pt in range(len(indx)):
                            pts_idx=[indx[pt][0]]
                            pts_num=[indx[pt][1]]
                            nfmod.nco_extract(outfile_temp,outfile_new, \
                                              [vdim], pts_idx, pts_num,ncksdir=ncopath)
                    
                            #newly-updated nc file for further processing, if any
                            os.system('mv -f ' + outfile_new + ' '+ outfile_temp)   # do the extraction continuously

                    if(is1D):
                        vdim = dim[0]
                    else:
                        vdim = dim[1]
                    if(len(indy)>0):
                        for pt in range(len(indy)):
                            pts_idx=[indy[pt][0]]
                            pts_num=[indy[pt][1]]
                            nfmod.nco_extract(outfile_temp,outfile_new, \
                                              [vdim], pts_idx, pts_num,ncksdir=ncopath)
                    
                            #newly-updated nc file for further processing, if any
                            os.system('mv -f ' + outfile_new + ' '+ outfile_temp)   # do the extraction continuously
                 
                # (END) loop for all dimensions possibly in the nc file
                if(os.path.isfile(outfile_temp)): 
                    os.system('mv -f ' + outfile_temp + ' '+ outfile_new)
                    
            # (END) if-file contains 'fileheader'
        # (END) if-it's a file
    #(END) for-loop of all files in outdir directory

           
#------------END of pointCLM_data -----------------------------------------------------------------------
    

