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
from Modules_CLMoutput_nc4 import CLM_NcRead_1simulation
from Modules_CLMoutput_nc4 import CLMvar_1Dtseries

from pyproj import Proj
from pyproj import Transformer
from pyproj import CRS

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
parser.add_option("--elmheader", dest="elmheader", default="", \
                  help = "ELM output Netcdf file header with path but no .nc ")
parser.add_option("--daymet_elm_mapfile", dest="gridmap", default="daymet_elm_mapping.txt", \
                  help = "DAYMET tile 2D to 1D landmasked grid mapping file ")
parser.add_option("--elm_varname", dest="elm_varname", default="ALL", \
                  help = "ELM output Netcdf file's variable name to process, ALL for all")

(options, args) = parser.parse_args()

#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)
cwdir = options.workdir

if (options.elmheader == ''):
    print('MUST input file name header, including fullpath, by " --elmheader=???"')
    sys.exit()

if (options.gridmap == ''):
    print('MUST input daymet_elm_mapping txt file, including fullpath, by " --daymet_elm_mapfile=???"')
    sys.exit()
# 

if not cwdir.endswith('/'): cwdir = cwdir+'/'

elm_varname = options.elm_varname # 'ALL' by default

#------------------------------------------------------------------------------

# ELM out data reading
# Note: the output data is in 1D of gridcell or landgridcell
if (options.elmheader != ""):
    elmpathfileheader = options.elmheader
    
    ftype = 'nc'

    alldirfiles = sorted(glob.glob("%s*.%s" % (elmpathfileheader, ftype)))
    if(len(alldirfiles)<=0):
        sys.exit("No file exists - %s*.%s IN %s" %(elmpathfileheader, ftype, cwdir))
    else:
        print('Total Files of ELM outputs: '+str(len(alldirfiles)))

    ncfileheader = elmpathfileheader.split('/')[-1]
    elm_odir   = elmpathfileheader.replace(ncfileheader,'')
    if(elm_odir.strip()==''):elm_odir='./'
    
    # read-in mapping file
    mapfile = options.gridmap.strip()
    with open(mapfile, 'r') as f:
        # remove header txts
        next(f)
        try:
            data = [x.strip().split() for x in f]
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


    # some tests, should be commented out normally
    #zone=0
    #if zone==1:
    #    gidx=gidx[:10460]
    #    yidx=yidx[:10460]
    #    xidx=xidx[:10460]
    #elif zone==2:
    #    gidx=gidx[10460:]-10460
    #    yidx=yidx[10460:]
    #    xidx=xidx[10460:]

    # read-in datasets one by one
    for ncfile in alldirfiles:

        src = Dataset(ncfile,'r')
        ncformat = src.file_format
        
        ncfile_out = 'daymet-'+ncfile
        dst = Dataset(ncfile_out, mode='w',format=ncformat)

        print ('Processing - ', ncfile, '==> ', ncfile_out)
        
        #
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
    
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dname=name
            if name=='DTIME': dname='time'# rename 'DTIME' to 'time' so that VISIT can regcon. it
            dst.createDimension(dname, (len(dimension) if not dimension.isunlimited() else None))
            # unlimited dimension name to be included dst output
            if dimension.isunlimited:
                unlimited_dim = dname
    
        # add geox/geoy dimensions
        geox_dim = dst.createDimension('geox',  len(geox))
        geoy_dim = dst.createDimension('geoy',  len(geoy))
        
        vgeox = dst.createVariable('geox', np.float32, ('geox',))
        vgeox.units = 'meters'
        vgeox.long_name = 'Easting (Lambert Conformal Conic projection)'
        vgeox.standard_name = "projection_x_coordinate"
        vgeox[:] = geox
            
        vgeoy = dst.createVariable('geoy', np.float32, ('geoy',))
        vgeoy.units = 'meters'
        vgeoy.long_name = 'Northing (Lambert Conformal Conic projection)'
        vgeoy.standard_name = "projection_y_coordinate"
        vgeoy[:] = geoy

        vproj = dst.createVariable('lambert_conformal_conic', np.int32)
        vproj.grid_mapping_name = "lambert_conformal_conic"
        vproj.longitude_of_central_meridian = -100.
        vproj.latitude_of_projection_origin = 42.5
        vproj.false_easting = 0.
        vproj.false_northing = 0.
        vproj.standard_parallel = 25., 60.
        vproj.semi_major_axis = 6378137.
        vproj.inverse_flattening = 298.257223563 

        # copy all data in src, and do 1D-grid --> 2D-geox/geoy copy
        for name, variable in src.variables.items():
            vname = name
            if name=='DTIME': vname='time' # rename 'DTIME' to 'time' so that VISIT can regcon. it

            if ('ALL' not in elm_varname) and (name not in elm_varname.split(',')): 
                # Must include 'unlimited_dim' variable, otherwise output NC is empty
                if (name!=unlimited_dim): continue

            
            #if(ncfile == alldirfiles[0]): print(name)
            
            src_dims = variable.dimensions
            
            #
            if ('gridcell' in src_dims or 'lndgrid' in src_dims or 'n' in src_dims \
                or 'landunit' in src_dims \
                or 'column' in src_dims \
                or 'pft' in src_dims):
                try:
                    FillValue = variable._FillValue
                except Exception as e:
                    print(e)
                    FillValue = -9999
                src_data = np.asarray(src[name])# be ready for original data to re-shape if any below 

                new_dims = []
                i = -1
                SKIPPED = False
                for vdim in src_dims:
                    dim=vdim
                    if vdim=='DTIME': dim='time'# rename 'DTIME' to 'time' so that VISIT can regcon. it
                    i = i + 1
                    if dim in ('gridcell','lndgrid', 'n'): 
                        idim = i
                        new_dims.append('geoy')
                        new_dims.append('geox')
                    elif (dim in ('landunit','column','pft')) and \
                         ('gridcell' in src.dimensions.keys() \
                           or 'lndgrid' in src.dimensions.keys()
                           or 'n' in src.dimensions.keys()):
                        # in ELM output, column/pft dims are gridcells*col/patch
                        idim = i
                        #
                        # reshape of dim column/pft
                        # note: have to exclude NON-vegetated/bare land unit (i.e. lunit=1)
                        lunit=1
                        if (dim == 'landunit'): 
                            idx_lun = np.where(src.variables['land1d_ityplun'][...]==lunit)[0]
                        if (dim == 'column'): 
                            idx_lun = np.where(src.variables['cols1d_ityplun'][...]==lunit)[0]
                        if (dim == 'pft'): 
                            idx_lun = np.where(src.variables['pfts1d_ityplun'][...]==lunit)[0]
                        len_dim = idx_lun.size
                        if idim==0:
                            src_data = src_data[idx_lun,...]
                        elif idim==1:
                            src_data = src_data[:,idx_lun,...]
                        elif idim==2:
                            src_data = src_data[:,:,idx_lun,...]
                        elif idim==3:
                            src_data = src_data[:,:,:,idx_lun,...]
                        else:
                            print('Error - more than 4 dimension variable, not supported yet')
                            sys.exit(-1)
                        
                        
                        if('gridcell' in src.dimensions.keys()):
                            len_grid = src.dimensions['gridcell'].size
                        elif('lndgrid' in src.dimensions.keys()):
                            len_grid = src.dimensions['lndgrid'].size
                            
                        if(math.fmod(len_dim, len_grid)==0):
                            len_dim = int(len_dim/len_grid)
                            name_dim = dim+'_index'
                            
                            # only need to create a new dimension ONCE, if not yet
                            if (name_dim not in dst.dimensions.keys()):
                                dst.createDimension(name_dim, len_dim)
                            new_dims.append(name_dim)
                            
                            #reshape original data to be like from (totalpft=gidx*col/pft) --> (gidx,lunit/col/pft)
                            src_shp = src_data.shape
                            if idim == 0:
                                re_shp = (len(gidx),len_dim,)+src_shp[1:]
                            elif idim >= 1:
                                re_shp = src_shp[0:idim-1]+(len(gidx),len_dim,)+src_shp[idim+1:]
                            src_data = np.reshape(src_data, re_shp)
                            # better move col/pft dim forward, so that in order of (col/pft, gidx), 
                            # then when remapping gidx --> y/x, it's in order of (col/pft, geoy, geox) which VISIT can plot correctly
                            src_data = np.swapaxes(src_data, idim, idim+1)
                            idim = idim + 1 #swaping above actually moved 'gidx' 1 dimension backwardly
                        else:
                            print('size of', dim, len_dim,'is NOT multiple of grids',len_grid)
                            print('variable', name,'excluded in 2-D merging')
                            SKIPPED=True
                        
                        new_dims.append('geoy')
                        new_dims.append('geox')

                    else:
                        new_dims.append(dim)
                
                if SKIPPED: continue
                vdtype = variable.datatype
                if(vdtype!=src_data.dtype):
                    vdtype=src_data.dtype
                x = dst.createVariable(vname, vdtype, \
                                        dimensions=new_dims, \
                                        zlib=True, fill_value=FillValue)
                
                # re-assign data from gidx to (yidx,xidx)
                # gidx should be in order already
                dst_data = np.asarray(dst[vname])
                
                if idim == 0:
                    dst_data[yidx,xidx,] = src_data[gidx,]
                    dst[vname][:,] = dst_data
                elif idim == 1:
                    dst_data[:,yidx,xidx,] = src_data[:,gidx,]
                    dst[vname][:,:,] = dst_data
                elif idim == 2:
                    dst_data[:,:,yidx,xidx,] = src_data[:,:,gidx,]
                    dst[vname][:,:,:,] = dst_data
                elif idim == 3:
                    dst_data[:,:,:,yidx,xidx,] = src_data[:,:,:,gidx,]
                    dst[vname][:,:,:,:,] = dst_data
                elif idim == 4:
                    dst_data[:,:,:,:,yidx,xidx,] = src_data[:,:,:,:,gidx,]
                    dst[vname][:,:,:,:,:,] = dst_data
                else:
                    print('Error - more than at least 4 dimension variable, not supported yet')
                    sys.exit(-1)

                #for g in gidx:  # this is not good way
                #    if i==0:
                #        dst[vname][yidx[g],xidx[g],] = src[name][g,]
                #    elif i==1:
                #        dst[vname][:,yidx[g],xidx[g-1]-1,] = src[name][:,g-1,]
                #    elif i==2:
                #        dst[vname][:,:,yidx[g]-1,xidx[g-1]-1,] = src[name][:,:,g-1,]
                #    elif i==3:
                #        dst[vname][:,:,:,yidx[g-1]-1,xidx[g-1]-1,] = src[name][:,:,g-1,]
                #    else:
                #        print('Error - more than at least 4 dimension variable, not supported yet')
                #        sys.exit(-1)
                     
            else:
                vdims = src_dims
                if 'DTIME' in src_dims: 
                    i=src_dims.index('DTIME')
                    vdims = src_dims[:i]+('time',)+src_dims[i+1:]
                x = dst.createVariable(vname, variable.datatype, \
                                        dimensions = vdims, \
                                        zlib=True)
                dst[vname][:] = src[name][:]
        
            # copy variable attributes all at once via dictionary
            dst[vname].setncatts(src[name].__dict__)
        
        #
        src.close()
        dst.close()
        print ('DONE with ncfile: ', ncfile)
    # end of 'for ncfile in alldirfiles:'
    
# end of 'if (options.elmheader != ""):' 
#-------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------