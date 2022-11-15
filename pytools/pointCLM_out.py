#!/usr/bin/env python

#----------------------
# Extract point(s) data from ELM outputs
#
import os, sys
import glob
import re
import math
from optparse import OptionParser
import numpy as np

import netcdf_modules as nfmod

#-------------------Local functions --------------------------------------------
def pointLocator(lons_all, lats_all, lons_pt, lats_pt, is1D):
    
    xindxpts = []
    yindxpts = []
        
    if(lons_pt.strip()==''):
        lons_pt=[]
    else:                
        if ('str' in str(type(lons_pt))):#values with operators
            nondecimal = re.compile(r'[^\d.,:;-]+')
            v_pt = nondecimal.sub("",lons_pt)
            
            v_pt = np.float64(v_pt.split(":"))
            if("~" in lons_pt):#nearest point(s)
                pt_val=[]
                for v in v_pt:
                    if v<0.0: v=v+360.0
                    pt_nearest=np.argmin(abs(lons_all-v))
                    pt_nearest=np.unravel_index(pt_nearest,lons_all.shape)
                    pt_val=np.append(pt_val,lons_all[pt_nearest])
                    # if lons_all are multiple points, i.e. max. distance can be calculated,
                    # then 'pt_nearest' cannot be beyond that max. distance 
                    # (This is for 'lons_all' are from multiple files or datasets)
                    if(len(lons_all)>1):
                        max_dist = abs(np.diff(lons_all))
                        max_dist = np.argmax(max_dist)
                    else:
                        max_dist = abs(lons_all[pt_nearest]-v)
                    if(abs(lons_all[pt_nearest]-v)<=max_dist):
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
            v_pt = np.float64(v_pt.split(":"))
            if("~" in lats_pt):#nearest point(s)
                pt_val=[]
                for v in v_pt:
                    pt_nearest=np.argmin(abs(lats_all-v))
                    pt_nearest=np.unravel_index(pt_nearest,lats_all.shape)
                    pt_val=np.append(pt_val,lats_all[pt_nearest])
                    # if lats_all are multiple points, i.e. max. distance can be calculated,
                    # then 'pt_nearest' cannot be beyond that max. distance 
                    # (This is for 'lats_all' are from multiple files or datasets)
                    if(len(lats_all)>1):
                        max_dist = abs(np.diff(lats_all))
                        max_dist = np.argmax(max_dist)
                    else:
                        max_dist = abs(lats_all[pt_nearest]-v)
                    if(abs(lats_all[pt_nearest]-v)<=max_dist):
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
    
    xind=[]
    yind=[]
    if np.size(lons_pt)>0:  
        #xind=np.argwhere(np.isin(lons_all,lons_pt)).ravel() # cannot get duplicated indices
        for x in lons_pt:
            xind.append(np.argwhere(np.isin(lons_all,x)).ravel()[0])
    if np.size(lats_pt)>0:
        #yind=np.argwhere(np.isin(lats_all,lats_pt)).ravel() # cannot get duplicated indices
        for y in lats_pt:
            yind.append(np.argwhere(np.isin(lats_all,y)).ravel()[0])
    
    #print(xind, '/n', yind)
    return xind, yind



#-------------------Parse options-----------------------------------------------
parser = OptionParser()

parser.add_option("--outdir", dest="outdir", default="", \
                  help = 'output directory')
parser.add_option("--ofheader", dest="header", default="none", \
                  help = 'header of output files either in fullpath or filename in default path')
parser.add_option("--newdata_affix", dest="newdata_affix", default="", \
                  help = 'new dataset output with affix')
parser.add_option("--lons", dest="ptlon", default="", \
                      help = "point or ranged(<,<=,>=,>,~) longitude(s), separated with':'")
parser.add_option("--lats", dest="ptlat", default="", \
                      help = "point or ranged (<,<=,>=,>,~) latitude(s), separated with':'")
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
                coord_vars = [['lon','lat']];
                coord_dims = [['lon','lat']];
                
 
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
                            
                            #indx/indy
                            [pts_indx,pts_indy] = pointLocator(alllons, alllats, options.ptlon, options.ptlat, is1D)                
                    
                    #print('Points found: ', len(pts_indx), pts_indx, pts_indy)
                    
                    if(len(pts_indx)>0):
                        for ipt in range(len(pts_indx)):
                            wdigit=int(math.log10(len(pts_indx)))
                            dirfile_new = dirfile.replace(filehead,filehead_new+str(ipt).zfill(wdigit))
                            outfile_new = outdir+'/'+dirfile_new
                            
                            dim_string = ' -d lon,'+str(pts_indx[ipt])+','+str(pts_indx[ipt]) \
                                        +' -d lat,'+str(pts_indy[ipt])+','+str(pts_indy[ipt])
                            os.system(ncopath+'./ncks --no_abc -O'+ dim_string+ \
                                  ' '+outfile_temp.strip()+' -o '+outfile_new.strip())
                            
                        #
                        os.system('rm -rf '+outfile_temp.strip())
                        
                        #merge file
                        ierr = os.system('find '+outdir+'/ -name "'+filehead_new+'*'+ \
                                         '" | xargs ls | sort | '+ncopath+'./ncecat -O -h -o ' \
                                         +outdir+'/temp0.nc')
                        if(ierr!=0): 
                            raise RuntimeError('Error: ncecat ')
                        else:
                            os.system('find '+outdir+'/ -name "'+filehead_new+'*'+ \
                                         '" -exec rm {} \; ')
                        
                        #point dim 'record' from 'ncecat' swap to 'geox' position
                        ierr = os.system(ncopath+'./ncpdq -h -O -a lon,record '+outdir+'/temp0.nc -o '+outdir+'/temp1.nc')
                        if(ierr!=0): raise RuntimeError('Error: ncpdq ')
                        
                        # remove 'geox/geoy' dims
                        ierr = os.system(ncopath+'ncwa -h -O -a lat '+outdir+'/temp1.nc -o '+outdir+'/temp2.nc')
                        if(ierr!=0): raise RuntimeError('Error: ncwa ')
                        ierr = os.system(ncopath+'ncwa -h -O -a lon '+outdir+'/temp2.nc -o '+outdir+'/temp1.nc')
                        if(ierr!=0): raise RuntimeError('Error: ncwa ')
                        # rename 'record' dim to 'gridcell'
                        ierr = os.system(ncopath+'ncrename -h -O -d record,gridcell '+outdir+'/temp1.nc '+outdir+'/temp2.nc')
                        if(ierr!=0): raise RuntimeError('Error: ncrename ')
                        #rename temp filename to new file name
                        ierr = os.system(ncopath+'./ncks --no_abc -O -x -v lon,lat '+outdir+'/temp2.nc -o '+outdir+'/multi_'+dirfile.replace(filehead,filehead_new))
                        if(ierr!=0): 
                            raise RuntimeError('Error: ncks ')
                        else:
                            os.system('rm -f '+outdir+'/temp*.nc')



                    
            # (END) if-file contains 'fileheader'
        # (END) if-it's a file
    #(END) for-loop of all files in outdir directory

           
#------------END of pointCLM_data -----------------------------------------------------------------------
    

