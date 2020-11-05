#!/usr/bin/env python

import os, sys, time, math
import glob
import re
import numpy as np
from datetime import datetime, date
from matplotlib.dates import date2num, num2date

from optparse import OptionParser

from numpy import long, int16
from netCDF4 import Dataset
from copy import deepcopy

import netCDF4

from pyproj import Proj
from pyproj import Transformer
from pyproj import CRS

#------------------------------------------------------------------------------------------------------------------
def Daymet_ELM_mapinfo(mapfile):
    # read-in mapping file
    #mapfile = options.gridmap.strip()
    with open(mapfile, 'r') as f:
        
        next(f) # if to remove header txts
        
        try:
            data = [x.strip().split() for x in f]
        except:
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
    
    #xidx/yidx may be missing (<0)
    resx = 1000.0 #daymet cell resolution in meters
    resy = 1000.0
    if any(xidx<=0) or any(yidx<=0):
        #
        xmin = np.min(geox)
        xmax = np.max(geox)
        x = np.arange(xmin, xmax+resx, resx)
        ymin = np.min(geoy)
        ymax = np.max(geoy)
        y = np.arange(ymin, ymax+resy, resy)
        
        for idx in range(len(gidx)):
            ii=np.argmin(abs(geox[idx]-x))
            jj=np.argmin(abs(geoy[idx]-y))
            xidx[idx] = ii
            yidx[idx] = jj
        geox = deepcopy(x)
        geoy = deepcopy(y)
        
    else:
        # xidx/yidx is really actual indices
        # geox/geoy need to sort by xidx/yidx and removal of duplicate
        # the following approach can guaranttee xidx/yidx match with original geox/geoy order 
        [idx, i] = np.unique(xidx, return_index=True)
        ii = np.argsort(idx)
        idx = i[ii]
        geox = geox[idx]

        [idx, i] = np.unique(yidx, return_index=True)
        ii = np.argsort(idx)
        idx = i[ii]
        geoy = geoy[idx]

    # output 2-D grid net geox/geoy, mapping index of  1D gidx <==> (xidx,yidx)
    return lon, lat, geox, geoy, xidx, yidx, gidx


#------------------------------------------------------------------------------------------------------------------
def Daymet_ELM_gridmatching(Grid1_Xdim, Grid1_Ydim, Grid2_x, Grid2_y, \
                         Grid1ifxy=False, Grid2ifxy=True, Grid1_cells=()):
    
    # Match Grid2 within each Grid1, so return Grid2-xy index for each Grid1 cell
    # by default, (1) Grid1 (parent) in lon/lat (aka ifxy=False), Grid2 in geox/y (aka ifxy=True)
    #             (2) Grid1 XY is grid-mesh-nodes, Grid2xy is grid-centroids. Then is good for searching Grid2 in Grid1
    #             (3) all cells in Grid1 are assigned Grid2's cell-index - maybe only those indiced in 'Grid1_cells'
    
    # it's supposed that: Grid1 X/Y are in 1-D regularly-intervaled (may not evenly) nodes along axis
    # while, Grid2 might be either like Grid2 or in 2-D mesh.
    if (len(Grid2_x.shape)<2): 
        # Grid2 must be converted to 2D paired x/y mesh, if not
        Grid2_xx, Grid2_yy = np.meshgrid(Grid2_x, Grid2_y) # mid-points of grid
    elif (len(Grid2_x.shape)==2):
        # Grid2 grid-centroids are in paired x/y for each grid
        Grid2_xx = Grid2_x
        Grid2_yy = Grid2_y
        
    if (len(Grid1_Xdim.shape)==1): #  Grid1 mesh in TWO 1-D dimensional nodes
        Grid1_x = Grid1_Xdim
        Grid1_y = Grid1_Ydim
        Grid1_xx, Grid1_yy = np.meshgrid(Grid1_Xdim, Grid1_Ydim) # nodes of grid-mesh
    else:
        #Grid1 mesh in 2-D for X/Y axis 
        print ('TODO - matching range-Grid1 in 2D mesh')
        sys.exit()
    
    # For projection conversion
    #     short lambert_conformal_conic ;
    #    lambert_conformal_conic:grid_mapping_name = "lambert_conformal_conic" ;
    #    lambert_conformal_conic:longitude_of_central_meridian = -100. ;
    #    lambert_conformal_conic:latitude_of_projection_origin = 42.5 ;
    #    lambert_conformal_conic:false_easting = 0. ;
    #    lambert_conformal_conic:false_northing = 0. ;
    #    lambert_conformal_conic:standard_parallel = 25., 60. ;
    #    lambert_conformal_conic:semi_major_axis = 6378137. ;
    #    lambert_conformal_conic:inverse_flattening = 298.257223563 ;
    #Proj4: +proj=lcc +lon_0=-100 +lat_1=25 +lat_2=60 +k=1 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m  +no_defs
    geoxy_proj_str = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
    geoxyProj = CRS.from_proj4(geoxy_proj_str)

    # EPSG: 4326
    # Proj4: +proj=longlat +datum=WGS84 +no_defs
    lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
    
    # only if 2 grids are in different projections, do tansformation
    if (Grid2ifxy and not Grid1ifxy):
        Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)
        Grid2_gxx,Grid2_gyy = Txy2lonlat.transform(Grid2_xx,Grid2_yy)
        
        ij=np.where(Grid2_gxx<0.0)
        if(len(ij[0])>0): Grid2_gxx[ij]=Grid2_gxx[ij]+360.0 # for convenience, longitude from 0~360
        ij=np.where(Grid1_x<0.0)
        if(len(ij[0])>0): Grid1_x[ij]=Grid1_x[ij]+360.0 # for convenience, longitude from 0~360

    elif (not Grid2ifxy and Grid1ifxy):
        Tlonlat2xy = Transformer.from_proj(lonlatProj, geoxyProj, always_xy=True)
        Grid2_gxx,Grid2_gyy = Tlonlat2xy.transform(Grid2_xx,Grid2_yy)
    
    else:
        Grid2_gxx = Grid2_xx
        Grid2_gyy = Grid2_yy

    # DAYMET grids' index (Grid2) included in each ELM land-grid (Grid1)
    Grid2in1_indx = {}
    if (len(Grid1_cells)<=0): 
        Grid1_ij = np.where(~np.isnan(Grid1_xx[:-1,:-1])) # cell-index rather than mesh-line index
    else:
        Grid1_ij = Grid1_cells
        
    for indx in range(len(Grid1_ij[0])): # Grid1 grid-cell no.
        j = Grid1_ij[0][indx]  # ELM output data is in (t,elmy,elmx) dimensional-order
        i = Grid1_ij[1][indx]
        
        iwst = np.min(Grid1_x[i:i+2])
        iest = np.max(Grid1_x[i:i+2])
        jsth = np.min(Grid1_y[j:j+2])
        jnth = np.max(Grid1_y[j:j+2])
        ij = np.where( ((Grid2_gxx<=iest) & (Grid2_gxx>iwst)) & \
                       ((Grid2_gyy<=jnth) & (Grid2_gyy>jsth)) )
        Grid2in1_indx[str(indx)] = deepcopy(ij)
            
        if False: # comment out the following - not correct
        #if(len(ij[0])<1):
             # none of DAYMET cell centroid inside a ELM grid, find the close one instead
            closej  = np.where((Grid2_gyy<=jnth) & (Grid2_gyy>jsth)) # do lat/y first, due to evenly-intervaled along lat/y
            if closej[0].size<=0:
                closei  = np.where((Grid2_gxx<=iest) & (Grid2_gxx>iwst)) # do lon/x first
                if(closei[0].size>0):
                    closeiy = np.argmin(abs(Grid2_gyy[closei]-(jnth+jsth)/2.0))
                    closeij = (np.asarray(closei[0][closeiy]),np.asarray(closei[1][closeiy]))
                else:
                    closeij = deepcopy(closei)
            else:
                closejx  = np.argmin(abs(Grid2_gxx[closej]-(iwst+iest)/2.0))
                closeij = (np.asarray(closej[0][closejx]),np.asarray(closej[1][closejx]))
            if len(closeij[0]>0):
                Grid2in1_indx[str(indx)] = deepcopy(closeij)
    
    # done with all grids
    return Grid2in1_indx, Grid2_gxx, Grid2_gyy
 
#------------------------------------------------------------------------------------------------------------------

#-------------------Parse options-----------------------------------------------

# CONVERT 1-D gridcell output NC file to 2-D geox/geoy output NC file 
# e.g. for plotting or comparison
#
parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
parser.add_option("--workdir2", dest="workdir2", default="", \
                  help = "data work directory2 (default = "", i.e., NO other directory for data-merge)")
parser.add_option("--daymet_elm_mapfile", dest="gridmap", default="daymet_elm_mapping.txt", \
                  help = "DAYMET tile 2D to 1D landmasked grid mapping file ")

(options, args) = parser.parse_args()

#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)
cwdir = options.workdir

if (options.workdir2 != ''):
    print('outputs will be merged from: '+ options.workdir2+' into' +options.workdir)
    workdir2=options.workdir2.strip().split(',') # multiple directories, separated by ','

if (options.gridmap == ''):
    print('MUST input daymet_elm_mapping txt file, including fullpath, by " --daymet_elm_mapfile=???"')
    sys.exit()
# 

if not cwdir.endswith('/'): cwdir = cwdir+'/'

#------------------------------------------------------------------------------

# ELM out data reading
# Note: the output data is in 1D of gridcell or landgridcell
if (options.gridmap != ""):

    # mapping file reading for the first (primary) workdir
    mapfile = options.gridmap.strip()
    [lon, lat, geox, geoy, xidx, yidx, gidx] = Daymet_ELM_mapinfo(mapfile)
    count_gid = np.asarray(len(gidx))
    if (len(workdir2)>0):
        for dir2 in workdir2:
            # mapfile is NOT with path, but in different dir2
            [lon2, lat2, geox2, geoy2, xidx2, yidx2, gidx2] = Daymet_ELM_mapinfo(dir2.strip()+'/'+mapfile)
            
            # x/yidx of the 1st in original/individual tiles
            xidx1st_1 = xidx[0]
            xidx1st_2 = xidx2[0]
            yidx1st_1 = yidx[0]
            yidx1st_2 = yidx2[0]
            geox1st_1 = geox[xidx1st_1]
            geox1st_2 = geox2[xidx1st_2]
            geoy1st_1 = geoy[yidx1st_1]
            geoy1st_2 = geoy2[yidx1st_2]
            
            #appending and removal of duplicates in in 2D grid x/y lines
            geox = np.append(geox, geox2)
            geoy = np.append(geoy, geoy2)
            geox = np.unique(geox)
            geoy = np.unique(geoy)
            
            #
            lon = np.append(lon, lon2)
            lat = np.append(lat, lat2)
            
            # gidx are those in individual nc files, so keep it as original
            gidx = np.append(gidx, gidx2)
            count_gid = np.append(count_gid, count_gid+len(gidx2))
            
            # xidx/yidx are in to-be-merged nc file, so need to redo them by adding x/y offsets
            xidx1st_1new = np.where(geox==geox1st_1)[0][0]  #np.where() output is a tuple, but really here needed is an interger index
            xidx = xidx + (xidx1st_1new - xidx1st_1)
            yidx1st_1new = np.where(geoy==geoy1st_1)[0][0]
            yidx = yidx + (yidx1st_1new - yidx1st_1)

            xidx1st_2new = np.where(geox==geox1st_2)[0][0]
            xidx2 = xidx2 + (xidx1st_2new - xidx1st_2)
            yidx1st_2new = np.where(geoy==geoy1st_2)[0][0]
            yidx2 = yidx2 + (yidx1st_2new - yidx1st_2)

            xidx = np.append(xidx, xidx2)
            yidx = np.append(yidx, yidx2)

    # ascii file write
    file_new = 'merged-'+mapfile
    f = open(file_new, 'w')
    fheader='      lon          lat           geox            geoy        i     j     g '
    f.write('%s \n' %fheader)
    
    # +1 so indexing starting from 1
    xidx = xidx + 1
    yidx = yidx + 1
    gidx = gidx + 1
    with open(file_new, 'a') as f:
        for lno in range(count_gid[-1]):
            f.write("%12.5f " % lon[lno])
            f.write("%12.6f " % lat[lno])
            f.write("%15.1f " % geox[xidx[lno]-1])
            f.write("%15.1f " % geoy[yidx[lno]-1])
            f.write("%5i " % xidx[lno])
            f.write("%5i " % yidx[lno])
            f.write("%5i\n" % gidx[lno])
        f.close()
# end of 'if (options.elmheader != ""):' 
#-------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------