#!/usr/bin/env python

import os, sys, time, math
import glob
import re
import numpy as np
#from datetime import datetime, date
from optparse import OptionParser
from copy import deepcopy

import netcdf_modules as ncmod


#------------------------------------------------------------------------------------------------------------------
def Daymet_ELM_mapinfo(mapfile, redoxy=False):
    # read-in mapping file
    #mapfile = options.gridmap.strip()
    with open(mapfile, 'r') as f:
        # remove header txts
        next(f)
        try:
            data = [x.strip().split() for x in f]
            f.close()
        except Exception as e:
            print(e)
            print('Error in reading - '+mapfile)
            sys.exit(-1)
    data = np.asarray(data,np.float)
    lon=data[:,0]
    lat=data[:,1]
    geox=data[:,2]
    geoy=data[:,3]
    xidx=np.asanyarray(data[:,4],np.int)-1  # in mappings.txt, those index are 1-based
    yidx=np.asanyarray(data[:,5],np.int)-1
    gidx=np.asanyarray(data[:,6],np.int)-1
    
    # re-do yidx/xidx, which from orginal file are error due to a bug
    missing_idx = np.any(xidx<0) or np.any(yidx<0)
    if redoxy or missing_idx:
        xmin = np.min(geox)
        xmax = np.max(geox)
        xx = np.arange(xmin, xmax+1000.0, 1000.0)
        ymin = np.min(geoy)
        ymax = np.max(geoy)
        yy = np.arange(ymin, ymax+1000.0, 1000.0)
        #
        f = open(mapfile, 'w')
        fheader='   lon          lat            geox            geoy        i     j     g '
        f.write(fheader+'\n')
        for ig in range(len(gidx)):
            ii=np.argmin(abs(geox[ig]-xx))
            jj=np.argmin(abs(geoy[ig]-yy))
            xidx[ig] = ii
            yidx[ig] = jj
            gidx[ig] = ig
            
            # re-write daymet_elm_mapping.txt    
            #'(f12.5,1x,f12.6,1x, 2(f15.1, 1x),3(I5,1x))'
            f.write('%12.5f ' % lon[ig] )
            f.write('%12.6f ' % lat[ig] )
            f.write('%15.1f ' % geox[ig] )
            f.write('%15.1f ' % geoy[ig] )
            f.write('%5d ' % (xidx[ig]+1) )  #x/yidx were 0-based, but need to 1-based for the mapping file
            f.write('%5d ' % (yidx[ig]+1) )
            f.write('%5d ' % (gidx[ig]+1) )
            f.write('\n')
        f.close()
    #
    else:
        # xidx/yidx is really actual indices
        # geox/geoy need to sort by xidx/yidx and removal of duplicate
        # the following approach can guaranttee xidx/yidx match with original geox/geoy order 
        [idx, i] = np.unique(xidx, return_index=True)
        ii = np.argsort(idx)
        idx = i[ii]
        xx = geox[idx]

        [idx, i] = np.unique(yidx, return_index=True)
        ii = np.argsort(idx)
        idx = i[ii]
        yy = geoy[idx]

    # output 2-D grid net geox/geoy, mapping index of  1D gidx <==> (xidx,yidx)
    return lon, lat, geox, geoy, xx, yy, xidx, yidx, gidx


#------------------------------------------------------------------------------------------------------------------
parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
parser.add_option("--workdir2", dest="workdir2", default="", \
                  help = "data work directory2 (default = "", i.e., NO other directory for data-merge)")
parser.add_option("--daymet_elm_mapfile", dest="gridmap", default="daymet_elm_mapping.txt", \
                  help = "DAYMET tile 2D to 1D landmasked grid mapping file ")
parser.add_option("--mapfile_only", dest="mapfile_only", default=False, \
                  help = " ONLY merge mapfile ", action="store_true")
parser.add_option("--mapfile_redoxy", dest="redoxy", default=False, \
                  help = " redo x/y index in mapfile ", action="store_true")
parser.add_option("--subzone_gridnumber", dest="zone_gnumber", default=36000, \
                  help = " new tile zone grid numbers ")
parser.add_option("--ncobinpath", dest="ncobinpath", default="", \
                      help = "NCO bin path if not in $PATH")

(options, args) = parser.parse_args()

cwdir = './'

#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)

if (options.workdir2 != ''):
    print('outputs will be merged from: '+ options.workdir2+' into: ' +options.workdir)
    workdir2=options.workdir2.strip().split(',') # multiple directories, separated by ','

if (options.gridmap == ''):
    print('MUST input daymet_elm_mapping txt file, including fullpath, by " --daymet_elm_mapfile=???"')
    sys.exit()
# 

if not options.workdir.endswith('/'): options.workdir = options.workdir+'/'
if not cwdir.endswith('/'): cwdir = cwdir+'/'

#------------------------------------------------------------------------------

# 
# Note: the output data is in 1D of gridcell or landgridcell
if True:
    
    pathfileheader = options.workdir+'GSWP3_daymet4'
    
    # mapping file reading for the first (primary) workdir
    mapfile = options.gridmap.strip()  # for 1D <--> 2D
    [lon, lat, geox, geoy, xx, yy, xidx, yidx, gidx] = Daymet_ELM_mapinfo(options.workdir+mapfile, redoxy=options.redoxy)
    cumsum_gid = len(gidx)
    count_gid  = [np.asarray(cumsum_gid)]
    tile_gid = [np.asarray(cumsum_gid)]
    tile_gid[:] = gidx
    tile_dir = np.zeros(len(gidx), dtype=object)
    tile_dir[:] = options.workdir
    
    zline_file = cwdir+'zone_mappings.txt' # this is for CPL_BYPASS in ELM
    if(options.workdir+'zone_mappings.txt'==zline_file): # don't over-write the first one
        zline_file = cwdir+'merged-zone_mappings.txt'
    
    if (options.workdir2!=''):
        for dir2 in workdir2:
            # mapfile is NOT with path, but in different dir2
            [lon2, lat2, geox2, geoy2, xx2, yy2, xidx2, yidx2, gidx2] = Daymet_ELM_mapinfo(dir2.strip()+'/'+mapfile, redoxy=options.redoxy)
            
            # x/yidx of the 1st in original/individual tiles
            xidx1st_1 = xidx[0]
            xidx1st_2 = xidx2[0]
            yidx1st_1 = yidx[0]
            yidx1st_2 = yidx2[0]
            xx1st_1 = xx[xidx1st_1]
            xx1st_2 = xx2[xidx1st_2]
            yy1st_1 = yy[yidx1st_1]
            yy1st_2 = yy2[yidx1st_2]
            
            #appending 2D grid x/y lines and resort them (as dimensions)
            xx = np.unique(np.append(xx, xx2))
            yy = np.unique(np.append(yy, yy2))
            
            
            lon = np.append(lon, lon2)
            lat = np.append(lat, lat2)
            geox = np.append(geox, geox2)
            geoy = np.append(geoy, geoy2)
            if (options.mapfile_only): 
                # gidx are those in individual nc files, so keep it as original if not merge nc files
                gidx = np.append(gidx, gidx2)
            else:
                # gidx are those in individual nc files, need to add it up when merge nc files
                gidx = np.append(gidx, gidx[-1]+1+gidx2)
            cumsum_gid = cumsum_gid + len(gidx2)
            count_gid = np.append(count_gid, cumsum_gid)
            
            tmp_dir = np.zeros(len(gidx2), dtype=object)
            tmp_dir[:] = dir2
            tile_dir = np.append(tile_dir, tmp_dir)
            tile_gid = np.append(tile_gid, gidx2)
            
            # xidx/yidx are in to-be-merged nc file, so need to redo them by adding x/y offsets
            xidx1st_1new = np.where(xx==xx1st_1)[0][0]  #np.where() output is a tuple, but really here needed is an interger index
            xidx = xidx + (xidx1st_1new - xidx1st_1)
            yidx1st_1new = np.where(yy==yy1st_1)[0][0]
            yidx = yidx + (yidx1st_1new - yidx1st_1)
            
            xidx1st_2new = np.where(xx==xx1st_2)[0][0]
            xidx2 = xidx2 + (xidx1st_2new - xidx1st_2)
            yidx1st_2new = np.where(yy==yy1st_2)[0][0]
            yidx2 = yidx2 + (yidx1st_2new - yidx1st_2)
            
            xidx = np.append(xidx, xidx2)
            yidx = np.append(yidx, yidx2)
        
    # daymet_elm_mapping.txt, zone_mappings.txt
    
    mapfile_out = cwdir+options.gridmap.strip().split('/')[-1]
    if(options.workdir+mapfile==mapfile_out):
        mapfile_out = cwdir+'merged-'+options.gridmap.strip().split('/')[-1]
    f = open(mapfile_out, 'w')
    fheader='     lon          lat        geox            geoy       i     j     g '
    f.write(fheader+'\n')
    
    f2 = open(zline_file, 'w')
    
    # newly-creating subzone
    z_gsize = int(options.zone_gnumber)       # this is max.
    z_num = math.ceil(len(gidx)/z_gsize)
    z_gsize = math.floor(len(gidx)/z_num)     # this is the actual
    z_gres = len(gidx)%z_gsize                # this number of residual will distribute to first 'g_res' zones with 1 for each
    z_gid = gidx.copy()                       # original gidx for new subzone
    z_gid[:] = 0
    
    g_0 = 0
    g_accsum = z_gsize
    g_zno = 1
    if (g_zno<=z_gres): g_accsum = g_accsum + 1
    for ig in range(len(gidx)):
        #'(f12.5,1x,f12.6,1x, 2(f15.1, 1x),3(I5,1x))'
        f.write('%12.5f ' % lon[ig] )
        f.write('%12.6f ' % lat[ig] )
        f.write('%15.1f ' % geox[ig] )
        f.write('%15.1f ' % geoy[ig] )
        f.write('%5d ' % (xidx[ig]+1) )  #x/yidx were 0-based, but need to 1-based for the mapping file
        f.write('%5d ' % (yidx[ig]+1) )
        f.write('%5d ' % (gidx[ig]+1) )
        f.write('\n')
        
        if ig>=g_accsum: 
            g_0 = g_accsum
            g_accsum = g_accsum + z_gsize
            g_zno = g_zno + 1
            if (g_zno<=z_gres): g_accsum = g_accsum + 1
        z_gid[ig] = g_zno
        
        f2.write('%12.5f ' % lon[ig] )
        f2.write('%12.6f ' % lat[ig] )
        f2.write('%5d ' % g_zno )
        f2.write('%5d ' % (gidx[ig]+1-g_0) )
        f2.write('\n')
        
    
    f.close()
    f2.close()
    if options.mapfile_only: os.sys.exit()
    
    
    #--------------------------------------------------------------------------------------------
    # data
    ftype = 'nc'
    alldirfiles = sorted(glob.glob("%s*.%s" % (pathfileheader, ftype))) # first daymet tile
    if(len(alldirfiles)<=0):
        sys.exit("No file exists - %s*.%s IN %s" %(pathfileheader, ftype, options.workdir))
    else:
        print('Total Files of: '+str(len(alldirfiles)))
    
    # read-in datasets one by one and do merging/divding
    
    for zno in range(z_num):
        # output directories
        zno_name = str(int(zno+1)).zfill(3)
        ncdirout = "DOMAIN_SUB"+zno_name
        if (os.path.exists(ncdirout)): os.system('rm -rf '+ncdirout)
        os.system('mkdir '+ncdirout)
        
        for ncfile in alldirfiles:
            
            # output filenames
            ncfname = ncfile.split('/')[-1]
            if re.search("_z??", ncfname): 
                z0str = ncfname.split('_')[-1]
                z0str = z0str.split('.')[0]
                ncfileout = ncdirout+"/"+ncfname.replace(z0str,'z'+zno_name)
            else:
                ncfileout = ncdirout+"/"+ncfname
            
            zworkdir=tile_dir[np.where(z_gid==zno+1)]
            zworkdir=np.unique(zworkdir)
           
            for dir in zworkdir:
                file_matched = False
                ncfile2 = sorted(glob.glob("%s" % (dir.strip()+'/'+ncfname)))
                if len(ncfile2)==1: 
                    file_matched = True
                elif len(ncfile2)>1:
                    print('Error: multiple files matched found in: '+ dir.strip())
                    print(ncfile2)
                    sys.exit(-1)
                else:
                    print('Warning: NO file matched found in: '+ dir.strip())
                    print('Warning: skip file merging: '+ncfile)
                
                if(not file_matched): 
                    break # break 'for dir2 in workdir2:'
                else:
                    print ('Processing - ', ncfile2[0], '==> ', ncfileout)
                    
                    gidx=tile_gid[np.where((tile_dir==dir) & (z_gid==zno+1))]
                    if (options.ncobinpath==""):
                        ncmod.nco_extract(ncfile2[0], 'temp.nc', ['n'], 
                                [gidx[0]], [len(gidx)])
                    else:
                        ncmod.nco_extract(ncfile2[0], 'temp.nc', ['n'], 
                                [gidx[0]], [len(gidx)],ncksdir=options.ncobinpath)
                    
                    if (os.path.isfile(ncfileout)):
                        os.system('mv '+ncfileout+' temp0.nc')
                        ncmod.mergefilesby1dim('temp0.nc', 'temp.nc', ncfileout, 'n')
                    else:
                        os.system('mv temp.nc '+ncfileout)
                    
                    
            #end of for dir in zworkdir
        
        #end of for zno in zones
        os.system('rm -f temp.nc temp0.nc')
        
        print ('DONE with ncfile: ', ncfile)
    # end of 'for ncfile in alldirfiles:'
    
    
# end
#-------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------
