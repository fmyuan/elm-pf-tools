#!/usr/bin/env python

import os, sys
import glob

from matplotlib.dates import date2num, num2date
from datetime import datetime, date
from copy import deepcopy

import numpy as np
import netCDF4

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

#------------------------------------------------------------------

# Write to geo-referenced CF compliant nc file, if filename given

def Write1GeoNc(vars, vardatas, ptxy=[], ncfname='', newnc=True, FillValue=None):
    # INPUTS: vars      - variable names, separated by ','. if 'all' means every var-key in 'vardatas'
    #         vardatas  - python np list, with dicts/data
    #    (optional) ptxy    - paired lon/lat (x/y), in [x/lon,y/lat] (only 1 point)
    #    (optional) ncfname - if not empty, write processed data into (geo)NC file
    #    (optional) newnc   - if not True, write into existed 'ncfname' nc file
    # OUTPUTS (optional): if NOT write to NC file and 'vars' only has 1 variable
    #                     Output year, doy, lon, lat, data for 'vars' in np.array  
    
    if ncfname !='':
        if ncfname.endswith('.nc4'): ncfname = ncfname.replace('.nc4','.nc')
        if not ncfname.endswith('.nc'): ncfname = ncfname+'.nc'
        if newnc:
            if os.path.isfile(ncfname): os.system('rm -rf '+ncfname)
            ncfile = netCDF4.Dataset(ncfname, mode='w',format='NETCDF4') 
            print('Create and Write NC file: '+ncfname)
        else:
            ncfile = netCDF4.Dataset(ncfname, mode='a',format='NETCDF4') 
            print('Write NC file: '+ncfname)


    # mid of day
    mid_day = vardatas['date']
    daynums = np.asarray([date2num(x) for x in mid_day])
    year    = np.asarray([x.year for x in mid_day])
    doy0    = np.asarray([date2num(date(x,1,1)) for x in year])
    doy     = daynums-doy0+1
    YEARLY  = False
    if(np.min(doy) == np.max(doy)): YEARLY = True
    try:
        nt = len(date2num(mid_day))
    except Exception as e:
        print(e)
        nt = 1

    # Construct the grid in lat/lon.
    if 'lat' in vardatas.keys():
        xlon = vardatas['lon']
        xlat = vardatas['lat']
        PROJECTED = False
    elif 'geox' in vardatas.keys():
        xlon = vardatas['geox']
        xlat = vardatas['geoy']
        PROJECTED = True
        
    # extracting pts, if specified
    if(len(ptxy)>1):
        d = abs(xlon-ptxy[0])
        ix = np.where(d==np.amin(d))
        lon = xlon[ix]
        d = abs(xlat-ptxy[1])
        iy = np.where(d==np.amin(d))
        lat = xlat[iy]
    else:
        lon = xlon
        lat = xlat

    #write to nc file
    DONE_header = False
    if not newnc: DONE_header = True
    DONE_time = False
    
    if vars[0]=='all': 
        vars=vardatas.keys()
        
    for varname in vars:
        #varname = 'Day_CMG_Snow_Cover'

        # header only needs to be done once
        if not DONE_header:
            if ncfname !='':
                # dimensions for nc file
                time_dim= ncfile.createDimension('time', None)

                if PROJECTED:
                    lon_dim = ncfile.createDimension('geox',  lon.size)
                    lat_dim = ncfile.createDimension('geoy',  lat.size)

                    vlat = ncfile.createVariable('geoy', np.float32, ('geoy',))
                    vlat.units = 'meters'
                    vlat.long_name = 'Northing (Lambert Conformal Conic projection)'
                
                    vlon = ncfile.createVariable('geox', np.float32, ('geox',))
                    vlon.units = 'meters'
                    vlon.long_name = 'Easting (Lambert Conformal Conic projection)'

                    vproj = ncfile.createVariable('lambert_conformal_conic', np.int32)
                    vproj.grid_mapping_name = "lambert_conformal_conic"
                    vproj.longitude_of_central_meridian = -100.
                    vproj.latitude_of_projection_origin = 42.5
                    vproj.false_easting = 0.
                    vproj.false_northing = 0.
                    vproj.standard_parallel = 25., 60.
                    vproj.semi_major_axis = 6378137.
                    vproj.inverse_flattening = 298.257223563 

                else:
                    lon_dim = ncfile.createDimension('lon',  lon.size)
                    lat_dim = ncfile.createDimension('lat',  lat.size)
                    vlat = ncfile.createVariable('lat', np.float32, ('lat',))
                    vlat.units = 'degrees_north'
                    vlat.long_name = 'latitude'
                
                    vlon = ncfile.createVariable('lon', np.float32, ('lon',))
                    vlon.units = 'degrees_east'
                    vlon.long_name = 'longitude'

                vlat[:] = lat
                vlon[:] = lon
        
                # time, create only
                vdaysnum = ncfile.createVariable('daysnum', np.float32, ('time',))
                vdaysnum.units = 'days'
                vdaysnum.long_name = 'days since 1980-01-01 UTC + 1'

                vdate = ncfile.createVariable('date', np.int32, ('time',))
                vdate.units = ''
                vdate.long_name = 'date in format of yyyymmdd'

                vdoy = ncfile.createVariable('doy', np.int16, ('time',))
                vdoy.units = 'doy'
                vdoy.long_name = 'day of year'

                vyear = ncfile.createVariable('year', np.int16, ('time',))
                vyear.units = 'year'
                vyear.long_name = 'year in format yyyy'

                # when using VISIT, it shows 'time' as cycle of multiple time-series
                vtime = ncfile.createVariable('time', np.int16, ('time',))
                if YEARLY:
                    vtime.units = 'year'
                    vtime.long_name = 'year in format yyyy'
                    
                else:
                    # DOY is much useful info in this case
                    vtime.units = 'doy'
                    vtime.long_name = 'day of year'

            
                #global attributes
                ncfile.data_source = ('Daymet Software Version 3.0,' +
                                      'Please see http://daymet.ornl.gov/ for current Daymet data citation information')
                if ('ELM' in ncfname):
                    ncfile.data_source2 = ('E3SM v1.1 Land Model offline simulations for No60 and above, ' +
                                           'For NGEE-Arctic Project sponsored by DOE Office of Science')

                if ('ELM' in ncfname):
                    ncfile.history = '2020-06-29: ELM forcing data assinged for each Daymet grid-cell.'
                else:
                    ncfile.history = '2020-06-29: Daymet data aggregating and re-projecting into ELM grid.'
                    
                ncfile.contact = 'F.-M. Yuan, CCSI/ESD-ORNL. yuanf@ornl.gov'
            # done if ncfname !='':
            
            DONE_header = True
        
        # write time
            
        if not DONE_time:
            
            yyyymmdd = [str(x).replace('-','') for x in mid_day ]
            if ncfname !='':
                if newnc:
                    vdaysnum[0:nt] = daynums
                    vdate[0:nt] = np.asarray([np.int32(x) for x in yyyymmdd])
                    vdoy[0:nt]  = doy
                    vyear[0:nt] = year
                    if YEARLY:
                        vtime[0:nt] = year
                    else:
                        vtime[0:nt] = doy
                        
                else:
                    vdaysnum = ncfile.variables['daysnum']
                    prv_nt = len(vdaysnum)
                    vdaysnum[prv_nt:prv_nt+nt] = daynums
                
                    vdate = ncfile.variables['date']
                    vdate[prv_nt:prv_nt+nt] = np.asarray([np.int32(x) for x in yyyymmdd])
                
                    vdoy = ncfile.variables['doy']
                    vdoy[prv_nt:prv_nt+nt] = doy

                    vyear = ncfile.variables['year']
                    vyear[prv_nt:prv_nt+nt] = year

                    vtime = ncfile.variables['time']
                    if YEARLY:
                        vtime[prv_nt:prv_nt+nt] = year
                    else:
                        vtime[prv_nt:prv_nt+nt] = doy

            # done write to nc (if ncfname !='':)
            
            DONE_time = True
        
        # 
        data = vardatas[varname]
        data = np.float32(data) # data type is 'uint8', convert to short (othwise cannot be read by Visit)
        
        if ncfname !='':
            if newnc:
                if PROJECTED:
                    vtemp = ncfile.createVariable(varname, np.float32, \
                                                  dimensions=('time','geoy','geox'), \
                                                  zlib=True, fill_value=FillValue) # note: unlimited dimension is leftmost
                else:
                    vtemp = ncfile.createVariable(varname, np.float32, \
                                                  dimensions=('time','lat','lon'), \
                                                  zlib=True, fill_value=FillValue) # note: unlimited dimension is leftmost
                
                vtemp.units = ''
                vtemp.standard_name = varname.strip() # this is a CF standard name
                
                pref='daily'
                if ('acc_' in varname): pref = 'accumulative'
                if ('prcp' in varname or 'SNOW' in varname or 'RAIN' in varname):
                    vtemp.long_name = pref+" precipitation/SNOW/RAIN"
                    vtemp.units = "mm"
                    vtemp.missing_value = -9999.
                    vtemp.coordinates = "geoy geox"
                    vtemp.grid_mapping = "lambert_conformal_conic"
                    vtemp.cell_methods = "area: mean time: sum"
                    vtemp.Key = "non-negative = daily/accumulative prcp/rain/snow, -22/nan = sea cell, -11 = land cell beyond data available" ;

                elif ('tmax' in varname):
                    vtemp.long_name = pref+ " max. temperature"
                    vtemp.units = "degree C"
                    vtemp.missing_value = -99.
                    vtemp.coordinates = "geoy geox"
                    vtemp.grid_mapping = "lambert_conformal_conic"
                    vtemp.cell_methods = "area: mean time: sum"
                    vtemp.Key = "real value = daily/accumulative max, -999/nan = sea cell, -99 = land cell beyond data available" ;

                elif ('tmin' in varname):
                    vtemp.long_name = pref+" min. temperature"
                    vtemp.units = "degree C"
                    vtemp.missing_value = -99.
                    vtemp.coordinates = "geoy geox"
                    vtemp.grid_mapping = "lambert_conformal_conic"
                    vtemp.cell_methods = "area: mean time: sum"
                    vtemp.Key = "real value = daily/accumulative min., -999/nan = sea cell, -99 = land cell beyond data available" ;

                elif ('tave' in varname):
                    vtemp.long_name = pref+ " average temperature"
                    vtemp.units = "degree C"
                    vtemp.missing_value = -99.
                    vtemp.coordinates = "geoy geox"
                    vtemp.grid_mapping = "lambert_conformal_conic"
                    vtemp.cell_methods = "area: mean time: sum"
                    vtemp.Key = "real value = daily/accumulative average, -999/nan = sea cell, -99 = land cell beyond data available" ;
                    
                vtemp[0:nt,:,:] = data.astype(np.float32)
            else:
                vtemp = ncfile.variables[varname]
                vtemp[prv_nt:prv_nt+nt,:,:] = data.astype(np.float32)
    # done with 'for varname in vars'
    
    # end of for varname in vars:
    if ncfname !='': ncfile.close()
    
def Daymet_TileInfo(inputdir, years=[1980],tiles=[],ncf0='TPHWLHrly/clmforc.Daymet4.1km.TPQWL.1980-01.nc',ncfout='daymet_tiles.nc', FillValue=None):
    # INPUTS: inputdir    - daymet data directory, under which data structured as: yyyy/tileno/
    #    (optional) years - by default, year 1980 under 'inputdir'
    #    (optional) tiles - by default, all tiles under 'inputdir/1980/'
    #    (optional) ncf0  - nc file name (only one needed) containing daymet tile info
    #    (optional) ncfout - 
    # OUTPUTS: by default, return and write to NC file for each pts
    #                     lon, lat, geox, geoy, tileno, xindx, yindx
    
    if not inputdir.endswith('/'): inputdir = inputdir+'/'
        
    # tiles to be checked
    if (len(tiles)<1):
        alldirs = sorted(glob.glob(os.path.join(inputdir+str(years[0])+'/','*')))
    else:
        alldirs = []
        for t in tiles:
            alldirs = np.concatenate(alldirs, 
                sorted(glob.glob(os.path.join((inputdir+str(years[0])+'/'+t,'*')))))
    
        alldirs = sorted(alldirs)

    # outputs  
    tile_no = np.empty(0)
    tile_gindx = np.empty(0)
    tile_xindx = np.empty(0)
    tile_yindx = np.empty(0)
                                  
    tile_lat = np.empty(0)
    tile_lon = np.empty(0)
    tile_geox = np.empty(0)
    tile_geoy = np.empty(0)

    
    #
    for idir in alldirs:
        itile = idir.split('/')[-1]
        ncfile = idir+'/'+ncf0
                                
        if os.path.isfile(ncfile):

            f = netCDF4.Dataset(ncfile,'r')
            geox = f.variables['x'][...]
            geoy = f.variables['y'][...]
            lats = f.variables['lat'][...]
            lons = f.variables['lon'][...]
            if 'FSDS' in f.variables.keys(): v='FSDS'
            if 'TBOT' in f.variables.keys(): v='TBOT'            
            data = f.variables[v][0,...]
            data_fill = f.variables[v]._FillValue
            
            # all points
            tlen = lats.size

            yindx, xindx = np.where((~data.mask))
            gindx = np.where((~data.reshape(tlen).mask))[0]
            tindx = np.empty_like(gindx); tindx[:]=itile
                
            tile_no = np.concatenate((tile_no,tindx))
            tile_gindx = np.concatenate((tile_gindx,gindx))
            tile_xindx = np.concatenate((tile_xindx,xindx))
            tile_yindx = np.concatenate((tile_yindx,yindx))
                                  
            tile_lat  = np.concatenate((tile_lat,lats[yindx,xindx]))
            tile_lon  = np.concatenate((tile_lon,lons[yindx,xindx]))
            tile_geox = np.concatenate((tile_geox,geox[xindx]))
            tile_geoy = np.concatenate((tile_geoy,geoy[yindx]))
            
            # 
            f.close()
              
    # end of for idir in alldirs
        
    # write to nc file

    if ncfout !='' and tile_gindx.size>0:
        
        if ncfout.endswith('.nc4'): ncfout = ncfout.replace('.nc4','.nc')
        if not ncfout.endswith('.nc'): ncfout = ncfout+'.nc'
        
        if os.path.isfile(ncfout): os.system('rm -rf '+ncfout)
        ncfo = netCDF4.Dataset(ncfout, mode='w',format='NETCDF4') 
        print('Create and Write NC file: '+ncfout)
            
        # dimensions for nc file

        grdim = ncfo.createDimension('grid',  tile_gindx.size)
        vgrd = ncfo.createVariable('grid', np.int64, ('grid',))
        vgrd.units = '-'
        vgrd.long_name = 'grid order'
        vgrd[:] = np.asarray(range(tile_gindx.size))

        vproj = ncfo.createVariable('lambert_conformal_conic', np.int32)
        vproj.grid_mapping_name = "lambert_conformal_conic"
        vproj.longitude_of_central_meridian = -100.
        vproj.latitude_of_projection_origin = 42.5
        vproj.false_easting = 0.
        vproj.false_northing = 0.
        vproj.standard_parallel = 25., 60.
        vproj.semi_major_axis = 6378137.
        vproj.inverse_flattening = 298.257223563 

        vgrdx = ncfo.createVariable('gindx', np.int64, ('grid',))
        vgrdx.units = '-'
        vgrdx.long_name = 'grid index in a flatten 2dx2d tile'
        vgrdx[:] = tile_gindx
        vxindx = ncfo.createVariable('xindx', np.int64, ('grid',))
        vxindx.units = '-'
        vxindx.long_name = 'x-coord index in a 2dx2d tile'
        vxindx[:] = tile_xindx
        vyindx = ncfo.createVariable('yindx', np.int64, ('grid',))
        vyindx.units = '-'
        vyindx.long_name = 'y-coord index in a flatten 2dx2d tile'
        vyindx[:] = tile_yindx

        vgeox = ncfo.createVariable('geox', np.float64, ('grid',))
        vgeox.units = 'meters'
        vgeox.long_name = 'Easting (Lambert Conformal Conic projection)'
        vgeox[:] = tile_geox
        vgeoy = ncfo.createVariable('geoy', np.float64, ('grid',))
        vgeoy.units = 'meters'
        vgeoy.long_name = 'Northing (Lambert Conformal Conic projection)'
        vgeoy[:] = tile_geoy                
        vlat = ncfo.createVariable('lat', np.float64, ('grid',))
        vlat.units = 'degrees_north'
        vlat.long_name = 'latitude'
        vlat[:] = tile_lat                
        vlon = ncfo.createVariable('lon', np.float64, ('grid',))
        vlon.units = 'degrees_east'
        vlon.long_name = 'longitude'
        vlon[:] = tile_lon
        
    
    # end of if ncfout != ''

    #
    return tile_no,tile_gindx,tile_xindx,tile_yindx,tile_geox,tile_geoy,tile_lat,tile_lon
    #
#----------------------------------------------------------------------------------------------

'''
print('Testing: ')
Daymet_TileInfo('./')
'''                    
