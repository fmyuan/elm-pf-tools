#!/usr/bin/env python

import os, sys, time, math
import glob
import re
import numpy as np
#from datetime import datetime, date
from optparse import OptionParser
from copy import deepcopy
from mpi4py import MPI

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
parser.add_option("--subdomain_gridnumber", dest="zone_grdnum", default=36000, \
                  help = " sud-domain grid numbers, by default 36,000 ")
parser.add_option("--subdomain_number", dest="zone_num", default=0, \
                  help = " sub-domain numbers , by default uppon to subdomain_gridnumber which can be overriden by giving a number over 0")
parser.add_option("--subdomain_cores", dest="zone_cores", default=6, \
                  help = " sub-domain on number of cores , by default 6 i.e. subdomain is for running at 6 cores ")
parser.add_option("--ncfile_header", dest="fileheader", default="GSWP3_daymet", \
                  help = "netcdf file name header, by default 'GSWP3_daymet*', i.e. high-res forcing data files only ")
parser.add_option("--ncobinpath", dest="ncobinpath", default="", \
                      help = "NCO bin path if not in $PATH")

(options, args) = parser.parse_args()

cwdir = './'
zdigit = 3

mycomm = MPI.COMM_WORLD
myrank = mycomm.Get_rank()
mysize = mycomm.Get_size()

#
if (options.workdir == './'):
    if myrank==0: print('data directory is the current')
else:
    if myrank==0: print('data directory: '+ options.workdir)

if (options.workdir2 != ''):
    workdir2=options.workdir2.strip().split(',') # multiple directories, separated by ','
    if myrank==0: print('outputs will be merged from: '+ workdir2[0]+' into: ' +options.workdir)

if (options.gridmap == ''):
    if myrank==0: print('MUST input daymet_elm_mapping txt file, including fullpath, by " --daymet_elm_mapfile=???"')
    sys.exit()
# 

if not options.workdir.endswith('/'): options.workdir = options.workdir+'/'
if not cwdir.endswith('/'): cwdir = cwdir+'/'

#------------------------------------------------------------------------------

#
pathfileheader = options.workdir+options.fileheader
mapfile = options.gridmap.strip()  # for 1D <--> 2D

mapfile_out = cwdir+options.gridmap.strip().split('/')[-1]
if(options.workdir+mapfile==mapfile_out):
    mapfile_out = cwdir+'merged-'+options.gridmap.strip().split('/')[-1]

zline_file = cwdir+'zone_mappings.txt' # this is for CPL_BYPASS in ELM
if(options.workdir+'zone_mappings.txt'==zline_file): # don't over-write the first one
    zline_file = cwdir+'merged-zone_mappings.txt'

gidx_len = None
z_num = None

#it's hard to mpi bcast np string array
wdir_all ={}
wdir_all[0] = options.workdir
if (options.workdir2!=''):
    for idir2 in range(len(workdir2)): 
        wdir_all[idir2+1]=workdir2[idir2]

# Note: the output data is in 1D of gridcell or landgridcell
if myrank==0:
    
    # mapping file reading for the first (primary) workdir
    [lon, lat, geox, geoy, xx, yy, xidx, yidx, gidx] = Daymet_ELM_mapinfo(options.workdir+mapfile, redoxy=options.redoxy)
    cumsum_gid = len(gidx)
    count_gid  = [np.asarray(cumsum_gid)]
    tile_gid = [np.asarray(cumsum_gid)]
    tile_gid[:] = gidx
    tile_idir = np.zeros(len(gidx), dtype=np.int)
    tile_idir[:] = 0
    
    
    if (options.workdir2!=''):
        for idir2 in range(len(workdir2)):
            dir2 = workdir2[idir2]
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
            
            tmp_idir = np.zeros(len(gidx2), dtype=np.int)
            tmp_idir[:] = idir2+1
            tile_idir = np.append(tile_idir, tmp_idir)
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
    f = open(mapfile_out, 'w')
    fheader='     lon          lat        geox            geoy       i     j     g '
    f.write(fheader+'\n')
    f2 = open(zline_file, 'w')
    
    # newly-creating subzone
    z_gsize = int(options.zone_grdnum)        # this is max.
    if (int(options.zone_num)>0):
        z_num = int(options.zone_num)
    else:
        z_num = math.ceil(len(gidx)/z_gsize)
    z_gsize = math.floor(len(gidx)/z_num)     # this is the actual
    z_gres = len(gidx)%z_gsize                # this number of residual will distribute to first 'g_res' zones with 6 for each and less than 6 for the last
    z_gid = gidx.copy()                       # original gidx for new subzone
    z_gid[:] = 0
    
    g_0 = 0
    g_zno = 1
    g_accsum = z_gsize
    z_cores = int(options.zone_cores)
    z_gres1 = math.floor(z_gres/z_cores)       # so z_gres = z_gres1*cores+z_gresx
    z_gresx = z_gres%z_cores
    if (g_zno<=z_gres1): 
        g_accsum = g_accsum + z_cores
    elif (g_zno==(z_gres1+1) and z_gresx>0):
        g_accsum = g_accsum + z_gresx
    
    fsub = open(mapfile_out+str(int(g_zno)).zfill(zdigit), 'w')
    fsub.write(fheader+'\n')
    fsub2 = open(zline_file+str(int(g_zno)).zfill(zdigit), 'w')
    
    gidx_len = len(gidx)
    for ig in range(gidx_len):
        if ig>=g_accsum:
            fsub.close()
            fsub2.close()
            
            g_0 = g_accsum
            g_zno = g_zno + 1
            g_accsum = g_accsum + z_gsize
            if (g_zno<=z_gres1): 
                g_accsum = g_accsum + z_cores
            elif (g_zno==(z_gres1+1) and z_gresx>0):
                g_accsum = g_accsum + z_gresx
                    
            fsub = open(mapfile_out+str(int(g_zno)).zfill(zdigit), 'w')
            fsub.write(fheader+'\n')
            fsub2 = open(zline_file+str(int(g_zno)).zfill(zdigit), 'w')
            
        z_gid[ig] = g_zno
        
        #merged mapping files
        #'(f12.5,1x,f12.6,1x, 2(f15.1, 1x),3(I5,1x))'
        f.write('%12.5f ' % lon[ig] )
        f.write('%12.6f ' % lat[ig] )
        f.write('%15.1f ' % geox[ig] )
        f.write('%15.1f ' % geoy[ig] )
        f.write('%5d ' % (xidx[ig]+1) )  #x/yidx were 0-based, but need to 1-based for the mapping file
        f.write('%5d ' % (yidx[ig]+1) )
        f.write('%5d ' % (gidx[ig]+1) )
        f.write('\n')
        f2.write('%12.5f ' % lon[ig] )
        f2.write('%12.6f ' % lat[ig] )
        f2.write('%5d ' % g_zno )
        f2.write('%5d ' % (gidx[ig]+1-g_0) )
        f2.write('\n')
        
        # sub mapping files
        fsub.write('%12.5f ' % lon[ig] )
        fsub.write('%12.6f ' % lat[ig] )
        fsub.write('%15.1f ' % geox[ig] )
        fsub.write('%15.1f ' % geoy[ig] )
        fsub.write('%5d ' % (xidx[ig]+1) )  #x/yidx were 0-based, but need to 1-based for the mapping file
        fsub.write('%5d ' % (yidx[ig]+1) )
        fsub.write('%5d ' % (gidx[ig]+1-g_0) )
        fsub.write('\n')
        fsub2.write('%12.5f ' % lon[ig] )
        fsub2.write('%12.6f ' % lat[ig] )
        fsub2.write('%5d ' % g_zno )
        fsub2.write('%5d ' % (gidx[ig]+1-g_0) )
        fsub2.write('\n')
        
    #
    f.close()
    f2.close()
    fsub.close()
    fsub2.close()
    
    
#--------------------------------------------------------------------------------------------
# the following need Bcast to each mpi rank
gidx_len = mycomm.bcast(gidx_len, root=0)
z_num = mycomm.bcast(z_num, root=0)

if myrank != 0:
    tile_idir = np.empty(gidx_len, dtype=np.int)
    tile_gid = np.empty(gidx_len, dtype=np.int)
    z_gid = np.empty(gidx_len, dtype=np.int)
else:
    print('Total Grids: ', gidx_len, 'Total Sub-Zones: ', z_num)
    
mycomm.Bcast([tile_idir, MPI.INT], root=0)
mycomm.Bcast([tile_gid, MPI.INT], root=0)
mycomm.Bcast([z_gid, MPI.INT], root=0)

# data
ftype = 'nc'
alldirfiles = sorted(glob.glob("%s*.%s" % (pathfileheader, ftype))) # first daymet tile
if(len(alldirfiles)<=0):
    sys.exit("No file exists - %s*.%s IN %s" %(pathfileheader, ftype, options.workdir))
else:
    if myrank==0: print('Total Files of: '+str(len(alldirfiles)))


mycomm.Barrier()

# read-in datasets one by one and do merging/divding
#if False: # to skip the following loop
if (not options.mapfile_only):
    
    # mpi implementation - simply round-robin 'z_num' over cpu_cores
    # (best if matching -np with known z_num) 
    nz_total = z_num
    nz = math.floor(nz_total/mysize)
    nz_rank = np.full([mysize], np.int(1))
    nz_rank = np.cumsum(nz_rank)*nz
    xz = int(math.fmod(nz_total, mysize))
    xz_rank = np.full([mysize], np.int(0))
    if xz>0: xz_rank[:xz]=1
    nz_rank = nz_rank + np.cumsum(xz_rank) - 1        # ending z index, starting 0, for each rank
    nz0_rank = np.hstack((0, nz_rank[0:mysize-1]+1))  # starting z index, starting 0, for each rank

    
    #for zno in range(z_num):
    for zno in range(nz0_rank[myrank], nz_rank[myrank]+1):
        # output directories
        zno_name = str(int(zno+1)).zfill(zdigit)
        ncdirout = "sub"+zno_name
        if (os.path.exists(ncdirout)): os.system('rm -rf '+ncdirout)
        os.system('mkdir '+ncdirout)
        
        fsub = mapfile_out+str(int(zno+1)).zfill(zdigit)
        fsub2 = zline_file+str(int(zno+1)).zfill(zdigit)
        os.system('mv '+fsub+' '+ncdirout+'/'+mapfile)
        os.system('mv '+fsub2+' '+ncdirout+'/zone_mappings.txt')
        #
        for ifile in range(len(alldirfiles)):
            ncfile = alldirfiles[ifile]
            
            # output filenames
            ncfname = ncfile.split('/')[-1]
            if re.search("_z??", ncfname): 
                z0str = ncfname.split('_')[-1]
                z0str = z0str.split('.')[0]
                ncfileout = ncdirout+"/"+ncfname.replace(z0str,'z'+zno_name)
            else:
                ncfileout = ncdirout+"/"+ncfname
            if (os.path.isfile(ncfileout) and zno==0): 
                os.system('rm -f '+ncfileout)
            
            zworkdirno=tile_idir[np.where(z_gid==zno+1)]
            zworkdirno=np.unique(zworkdirno)
            
            for idir in zworkdirno:
                dir = wdir_all[idir]
                file_matched = False
                ncfile2 = sorted(glob.glob("%s" % (dir.strip()+'/'+ncfname)))
                if len(ncfile2)==1: 
                    file_matched = True
                elif len(ncfile2)>1:
                    if myrank ==0: print('Error: multiple files matched found in: '+ dir.strip())
                    if myrank ==0: print(ncfile2)
                    sys.exit(-1)
                else:
                    if myrank ==0: print('Warning: NO file matched found in: '+ dir.strip())
                    if myrank ==0: print('Warning: skip file merging: '+ncfile)
                
                if(not file_matched): 
                    break # break 'for dir2 in workdir2:'
                else:
        
                    dim_name = ['n']
                    if ('domain' in ncfile):
                        dim_name = ['ni']
                    elif ('surfdata' in ncfile):
                        dim_name = ['gridcell']
                    
                    gidx=tile_gid[np.where((tile_idir==idir) & (z_gid==zno+1))]
                    tmpnc = 'temp_'+zno_name+'_'+str(ifile)+'.nc'
                    if (os.path.isfile(tmpnc)): os.system('rm -f '+tmpnc)
                    
                    if (options.ncobinpath==""):
                        ncmod.nco_extract(ncfile2[0], tmpnc, dim_name, 
                                [gidx[0]], [len(gidx)])
                    else:
                        ncmod.nco_extract(ncfile2[0], tmpnc, dim_name, 
                                [gidx[0]], [len(gidx)],ncksdir=options.ncobinpath)
                    
                    if (os.path.isfile(ncfileout)):
                       os.system('mv '+ncfileout+' '+tmpnc+'0')
                       ncmod.mergefilesby1dim(tmpnc+'0', tmpnc, ncfileout, dim_name[0])
                    else:
                       os.system('cp '+ tmpnc +' '+ncfileout)
                    
                
            #end of for idir in zworkdirno
            if (dir==wdir_all[len(wdir_all)-1] and zno==z_num-1): 
                print ('DONE with ncfile: ', ncfile.split('/')[-1], ' ON rank: ', myrank)
        
        #end of 'for ncfile in alldirfiles:'
        mycomm.Barrier()
        if (myrank==0): os.system('rm -f temp_*nc*')
        
    # end of 'for zno in range(z_num)'
    
# end if (not options.mapfile_only)
#-------------------------------------------------------------------------


MPI.Finalize()

#------------------------------------------------------------------------------------------------
