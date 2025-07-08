
# utility to generate 1D and 2D domain
# use previous dataPartition module. Code need to be clean up

import os
import math
import netCDF4 as nc
import numpy as np
from pyproj import Transformer
from pyproj import CRS

from pytools.commons_utils.Modules_netcdf import geotiff2nc

def domain_unstructured_fromraster_cavm(output_pathfile='./domain.lnd.pan-arctic_CAVM.1km.1d.c241018.nc', \
                                   rasterfile='./raster_cavm_v1.tif', outdata=False):
    
    # 1000mx1000m grids extracted from CAVM image (geotiff) file
    # e.g. raster_cavm_v1.tif
    #  

    file=rasterfile
    bandinfos={'bands':['grid_code']}
    allxc, allyc, crs_res, crs, alldata = geotiff2nc(file, bandinfos=bandinfos, outdata=True)
    alldata = alldata[0] # only 1 banded data

    # truncate non-data points as much as possible
    # cavm: 91 - freshwater, 92 - sea water (ocean), 93 - glacier, 99 - non-arctic
    data = np.ma.masked_equal(alldata,99)  # non-arctic grids
    data = np.ma.masked_equal(data,92)  # sea grids
    nondata = data.mask
    #along x-axis, checking non-data in columns
    ny_nodata = np.sum(nondata,0)
    x0=0; x1=len(allxc)-1
    for i in range(1,len(allxc)-1):
        if ny_nodata[i]<len(allyc) or ny_nodata[i-1]<len(allyc): break 
        x0 = i
    for i in range(len(allxc)-2,1,-1):
        if ny_nodata[i]<len(allyc) or ny_nodata[i+1]<len(allyc): break 
        x1 = i
    # note: column x0/x1 shall be all non-data, but x0+1/x1-1 not
    # since the projection is centered in northern pole, better to make truncated x still centered at pole
    X_cut = min(x0, len(allxc)-1-x1)
    x0 = X_cut; x1 = len(allxc)-1-X_cut
    X_axis = allxc[x0:x1]
    data = data[:,x0:x1]   
    
    #along y-axis, checking non-data in rows
    y0=0; y1=len(allyc)-1
    nx_nodata = np.sum(nondata,1)
    for i in range(1,len(allyc)-1):
        if nx_nodata[i]<len(allxc) or nx_nodata[i-1]<len(allxc): 
            if crs.is_projected:
                break
            elif allyc[i]>=50.0: #truncating southern lat of 50 degN
                break 
        y0 = i
    for i in range(len(allyc)-2,1,-1):
        if nx_nodata[i]<len(allxc) or nx_nodata[i+1]<len(allxc): break 
        y1 = i
    # note: row y0/y1 shall be all non-data, but y0+1/y1-1 not
    # since the projection is centered in northern pole, better to make truncated y still centered at pole
    Y_cut = min(y0, len(allyc)-1-y1)
    y0 = Y_cut; y1 = len(allyc)-1-Y_cut
    Y_axis = allyc[y0:y1]      
    data = data[y0:y1,:]   

    XC, YC = np.meshgrid(X_axis, Y_axis)

    # lon/lat from XC/YC, if in projected coordinates
    if crs.is_projected:    
        #Proj4: +proj=laea +lon_0=-180 +lat_0=90 +x_0=0 +y_0=0 +R=6370997 +f=0 +units=m +no_defs
        geoxy_proj_str = "+proj=laea +lon_0=-180 +lat_0=90 +x_0=0 +y_0=0 +R=6370997 +f=0 +units=m +no_defs"
        geoxyProj = CRS.from_proj4(geoxy_proj_str)
        
        # EPSG: 4326
        # Proj4: +proj=longlat +datum=WGS84 +no_defs
        epsg_code = 4326
        lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
        Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)
        #Tlonlat2xy = Transformer.from_proj(lonlatProj, geoxyProj, always_xy=True)

        lon,lat = Txy2lonlat.transform(XC,YC)
    
    else:
        lon = XC
        lat = YC

    # add the gridcell IDs. and its x/y indices
    total_rows = len(Y_axis)
    total_cols = len(X_axis)
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
    grid_ids = grid_ids.reshape(lat.shape)
    grid_xids = np.indices(grid_ids.shape)[1]
    grid_yids = np.indices(grid_ids.shape)[0]


    #create land gridcell mask, area, and landfrac (otherwise 0)
    masked = np.where(~data.mask)
    landmask = np.where(~data.mask, 1, 0)
    landfrac = landmask.astype(float)*1.0
    
    lat[lat==90.0]=lat[lat==90.0]-0.00001
    lat[lat==-90.0]=lat[lat==-90.0]-0.00001
    kmratio_lon2lat = np.cos(np.radians(lat))
    re_km = 6370.997
    area = np.empty_like(lat)
    area_km2 = np.empty_like(kmratio_lon2lat)
    if crs.is_projected:    
        # area in km2 --> in arcrad2
        area_km2[...] = crs_res[0]*crs_res[1]*1.0e-6
        yscalar = crs_res[1]/1000.0/(math.pi*re_km/180.0)
        xscalar = crs_res[0]/1000.0/(math.pi*re_km/180.0*kmratio_lon2lat)
        area[...] = xscalar*yscalar

    else:
        area[...] = np.radians(crs_res[0])*np.radians(crs_res[1])
        yscalar = crs_res[1]*(math.pi*re_km/180.0)
        xscalar = crs_res[0]*(math.pi*re_km/180.0*kmratio_lon2lat)
        area_km2[...] = xscalar*yscalar

    
    # 2d --> 1d, with masked only
    grid_id_arr = grid_ids[masked]
    grid_xids_arr = grid_xids[masked]
    grid_yids_arr = grid_yids[masked]
    mask_arr = landmask[masked]
    landfrac_arr = landfrac[masked]
    landcov_arr = data[masked]
    area_arr = area[masked]
    area_km2_arr = area_km2[masked]
    lat_arr = lat[masked]
    lon_arr = lon[masked]
    XC_arr = XC[masked]
    YC_arr = YC[masked]

    
    file_name = output_pathfile
    print("The domain file is " + file_name)

    # Open a new NetCDF file to write the domain information. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
    w_nc_fid.title = '1D domain file for pan arctic region, based on CAVM vegetation map'

    # Create new dimensions of new coordinate system
    w_nc_fid.createDimension('x', total_cols)
    w_nc_fid.createDimension('y', total_rows)
    dst_var = w_nc_fid.createVariable('x', np.float32, ('x'))
    if crs.is_projected:
        dst_var.units = "m"
        dst_var.long_name = "x coordinate of projection"
        dst_var.standard_name = "projection_x_coordinate"
    else:
        dst_var.units = 'degrees_east'
        dst_var.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
        dst_var.standard_name = "longitude"
    w_nc_fid['x'][...] = np.copy(XC[0,:])
    dst_var = w_nc_fid.createVariable('y', np.float32, ('y'))
    if crs.is_projected:
        dst_var.units = "m"
        dst_var.long_name = "y coordinate of projection"
        dst_var.standard_name = "projection_y_coordinate"
    else:
        dst_var.units = 'degrees_north'
        dst_var.long_name = 'latitude of land gridcell center (GCS_WGS_84), decreasing from north to south'
        dst_var.standard_name = 'latitude'

    w_nc_fid['y'][...] = np.copy(YC[:,0])
    w_nc_var = w_nc_fid.createVariable('lon', np.float64, ('y','x'))
    w_nc_var.long_name = 'longitude of 2D land gridcell center (GCS_WGS_84), increasing from west(-180) to east(180)'
    w_nc_var.units = "degrees_east"
    w_nc_fid.variables['lon'][...] = lon
    w_nc_var = w_nc_fid.createVariable('lat', np.float64, ('y','x'))
    w_nc_var.long_name = 'latitude of 2D land gridcell center (GCS_WGS_84), increasing from south to north'
    w_nc_var.units = "degrees_north"
    w_nc_fid.variables['lat'][...] = lat
    w_nc_var = w_nc_fid.createVariable('grid_code', np.float64, ('y','x'))
    w_nc_var.long_name = 'code of 2D land gridcells'
    w_nc_var.units = "-"
    w_nc_fid.variables['grid_code'][...] = data

    # create the gridIDs, lon, and lat variable
    w_nc_fid.createDimension('ni', grid_id_arr.size)
    w_nc_fid.createDimension('nj', 1)
    w_nc_fid.createDimension('nv', 4)

    w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId in the pan-Arctic domain of CAVM v1 map'
    w_nc_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells" 
    w_nc_fid.variables['gridID'][...] = grid_id_arr

    w_nc_var = w_nc_fid.createVariable('gridXID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId x in the pan-Arctic domain of CAVM v1 map'
    w_nc_var.decription = "start from #0 at the upper left corner and from west to east of the domain, with gridID=gridXID+gridYID*x_dim" 
    w_nc_fid.variables['gridXID'][...] = grid_xids_arr

    w_nc_var = w_nc_fid.createVariable('gridYID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId y in the pan-Arctic domain of CAVM v1 map'
    w_nc_var.decription = "start from #0 at the upper left corner and from north to south of the domain, with gridID=gridXID+gridYID*y_dim" 
    w_nc_fid.variables['gridYID'][...] = grid_yids_arr

    w_nc_var = w_nc_fid.createVariable('xc', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_var.bounds = "xv"
    w_nc_fid.variables['xc'][...] = lon_arr
        
    w_nc_var = w_nc_fid.createVariable('yc', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell center (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"
    w_nc_fid.variables['yc'][...] = lat_arr
        
    # create the XC, YC variable
    w_nc_var = w_nc_fid.createVariable('xc_LAEA', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'X of land gridcell center (Lambert Azimuthal Equal Area), increasing from west to east'
    w_nc_var.units = "m"
    w_nc_fid.variables['xc_LAEA'][...] = XC_arr
        
    w_nc_var = w_nc_fid.createVariable('yc_LAEA', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'Y of land gridcell center (Lambert Azimuthal Equal Area), decreasing from north to south'
    w_nc_var.units = "m"
    w_nc_fid.variables['yc_LAEA'][...] = YC_arr

    #
    w_nc_var = w_nc_fid.createVariable('xv', np.float64, ('nj','ni','nv'))
    w_nc_var.long_name = 'longitude of land gridcell verticles'
    w_nc_var.units = "degrees_east"
    w_nc_var.comment = 'by GCS_WGS_84, increasing from west (-180) to east (180), vertices ordering anti-clock from left-low corner'
    w_nc_var = w_nc_fid.createVariable('yv', np.float64, ('nj','ni','nv'))
    w_nc_var.long_name = 'latitude of land gridcell verticles'
    w_nc_var.units = "degrees_north"
    w_nc_var.comment = 'by GCS_WGS_84, decreasing from north to south, vertices ordering anti-clock from left-low corner'

    if crs.is_projected:
        xv_arr,yv_arr = Txy2lonlat.transform(XC_arr-crs_res[0]/2.0,YC_arr+crs_res[1]/2.0)
    else:
        xv_arr = XC_arr-crs_res[0]/2.0
        yv_arr = YC_arr+crs_res[1]/2.0
    w_nc_fid.variables['xv'][...,0] = xv_arr
    w_nc_fid.variables['yv'][...,0] = yv_arr
    if crs.is_projected:
        xv_arr,yv_arr = Txy2lonlat.transform(XC_arr+crs_res[0]/2.0,YC_arr+crs_res[1]/2.0)
    else:
        xv_arr = XC_arr+crs_res[0]/2.0
        yv_arr = YC_arr+crs_res[1]/2.0
    w_nc_fid.variables['xv'][...,1] = xv_arr
    w_nc_fid.variables['yv'][...,1] = yv_arr
    if crs.is_projected:
        xv_arr,yv_arr = Txy2lonlat.transform(XC_arr+crs_res[0]/2.0,YC_arr-crs_res[1]/2.0)
    else:
        xv_arr = XC_arr+crs_res[0]/2.0
        yv_arr = YC_arr-crs_res[1]/2.0
    w_nc_fid.variables['xv'][...,2] = xv_arr
    w_nc_fid.variables['yv'][...,2] = yv_arr
    if crs.is_projected:
        xv_arr,yv_arr = Txy2lonlat.transform(XC_arr-crs_res[0]/2.0,YC_arr-crs_res[1]/2.0)
    else:
        xv_arr = XC_arr-crs_res[0]/2.0
        yv_arr = YC_arr-crs_res[1]/2.0
    w_nc_fid.variables['xv'][...,3] = xv_arr
    w_nc_fid.variables['yv'][...,3] = yv_arr

    w_nc_var = w_nc_fid.createVariable('area', np.float64, ('nj','ni'))
    w_nc_var.long_name = "area of grid cell in radians squared" ;
    w_nc_var.coordinates = "xc yc" ;
    w_nc_var.units = "radian2" ;
    w_nc_fid.variables['area'][...] = area_arr

    w_nc_var = w_nc_fid.createVariable('area_LAEA', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells (Lambert Azimuthal Equal Area)'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "km^2"
    w_nc_fid.variables['area_LAEA'][...] = area_km2_arr

    w_nc_var = w_nc_fid.createVariable('landcov_code', np.int32, ('nj','ni'))
    w_nc_var.long_name = "land-unit or veg. code" ;
    w_nc_var.note = "unitless" ;
    w_nc_var.coordinates = "xc yc" ;
    w_nc_var.comment = "CAVM land/veg code, masked-off seawater and non-arctic grids" ;
    w_nc_fid.variables['landcov_code'][...] = landcov_arr

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('nj','ni'))
    w_nc_var.long_name = "domain mask" ;
    w_nc_var.note = "unitless" ;
    w_nc_var.coordinates = "xc yc" ;
    w_nc_var.comment = "0 value indicates cell is not active" ;
    w_nc_fid.variables['mask'][...] = mask_arr

    w_nc_var = w_nc_fid.createVariable('frac', np.float64, ('nj','ni'))
    w_nc_var.long_name = "fraction of grid cell that is active" ;
    w_nc_var.coordinates = "xc yc" ;
    w_nc_var.note = "unitless" ;
    w_nc_var.filter1 = "error if frac> 1.0+eps or frac < 0.0-eps; eps = 0.1000000E-11" ;
    w_nc_var.filter2 = "limit frac to [fminval,fmaxval]; fminval= 0.1000000E-02 fmaxval=  1.000000" ;
    w_nc_fid.variables['frac'][...] = landfrac_arr

    if crs.is_projected:
        w_nc_var = w_nc_fid.createVariable('lambert_azimuthal_equal_area', np.short)
        w_nc_var.grid_mapping_name = "lambert_azimuthal_equal_area"
        w_nc_var.crs_wkt = crs.wkt
        w_nc_var.longitude_of_center = 90.
        w_nc_var.latitude_of_center = -180.
        w_nc_var.false_easting = 0.
        w_nc_var.false_northing = 0.
        w_nc_var.semi_major_axis = 6370997.
        w_nc_var.inverse_flattening = 0.
        w_nc_var.proj4_str = geoxy_proj_str
        w_nc_var.epsg = epsg_code

    w_nc_fid.close()  # close the new file 
    
    if outdata:
        return XC, YC, data 

def domain_unstructured_fromdaymet(output_pathfile='./domain.lnd.Daymet4.1km.1d.c231120.nc', \
                                   tileinfo_ncfile='./daymet_tiles.nc', outdata=False):

    # 
    # Daymet LCC proj4 string
    #Proj4: +proj=lcc +lon_0=-100 +lat_1=25 +lat_2=60 +k=1 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m  +no_defs
    geoxy_proj_str = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
    geoxyProj = CRS.from_proj4(geoxy_proj_str)
    
    # EPSG: 4326
    # Proj4: +proj=longlat +datum=WGS84 +no_defs
    epsg_code = 4326
    lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
    Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)
    #Tlonlat2xy = Transformer.from_proj(lonlatProj, geoxyProj, always_xy=True)

    # centroids from input
    if tileinfo_ncfile!='':
        ncf = nc.Dataset(tileinfo_ncfile)
        lon = ncf.variables['lon'][...]
        lat = ncf.variables['lat'][...]
        geox = ncf.variables['geox'][...]
        geoy = ncf.variables['geoy'][...]
        #gidx = ncf.variables['gindx'][...]
        #xidx = ncf.variables['xindx'][...]
        #yidx = ncf.variables['yindx'][...]
        
        
    else:
        print('No daymet grid information file, e.g.', tileinfo_ncfile)
        return

    #create land gridcell mask and landfrac (assuming all 1 unless having inputs)
    landmask = np.ones_like(lon)
    
    landfrac = landmask.astype(float)*1.0   
    # area in km2 --> in arcrad2
    area_km2 = 1.0 # this is by default from daymet
    side_km = math.sqrt(float(area_km2))
    lat[lat==90.0]=lat[lat==90.0]-0.00001
    lat[lat==-90.0]=lat[lat==-90.0]-0.00001
    kmratio_lon2lat = np.cos(np.radians(lat))
    re_km = 6371.22
    yscalar = side_km/(math.pi*re_km/180.0)
    xscalar = side_km/(math.pi*re_km/180.0*kmratio_lon2lat)
    area = xscalar*yscalar

    # build a coordinate system for ALL (geox, geoy) of centroids
    # this is good to mapping
    [X_axis, i] = np.unique(geox, return_inverse=True)
    [Y_axis, j] = np.unique(geoy, return_inverse=True)
    YY, XX = np.meshgrid(Y_axis, X_axis, indexing='ij')
    LON2D,LAT2D = Txy2lonlat.transform(XX,YY)
    GridXID = np.indices(XX.shape)[1]
    GridYID = np.indices(XX.shape)[0]
    GridID = np.indices(XX.flatten().shape)[0]
    GridID = np.reshape(GridID, GridXID.shape)
    
    # (2d --> 1d) grid indices with pts (lat,lon) only
    grid_id_arr = GridID[(j,i)].flatten()
    grid_xids_arr = GridXID[(j,i)].flatten()
    grid_yids_arr = GridYID[(j,i)].flatten()
    mask_arr = landmask.flatten()
    landfrac_arr = landfrac.flatten()
    area_arr = area.flatten()
    lat_arr = lat.flatten()
    lon_arr = lon.flatten()
    geox_arr = geox.flatten()
    geoy_arr = geoy.flatten()


    # write to domain.nc    
    file_name = output_pathfile
    print("The domain file is " + file_name)

    # Open a new NetCDF file to write the domain information. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
    w_nc_fid.title = '1D (unstructured) domain file for the Daymet grid coordinate system'
    w_nc_fid.Conventions = "CF-1.0"
    w_nc_fid.note = "a 2d-grid coordinates system, ('y','x'), added, with 2d lon/lat variables. " + \
                    " With this, GridXID, GridYID, and GridID added for actual grids(nj,ni) " + \
                    " so that transform between 2d and 1d may be possible"

    # Create new dimensions of new coordinate system
    w_nc_fid.createDimension('x', len(X_axis))
    w_nc_fid.createDimension('y', len(Y_axis))
    dst_var = w_nc_fid.createVariable('x', np.float32, ('x'))
    dst_var.units = "m"
    dst_var.long_name = "x coordinate of projection"
    dst_var.standard_name = "projection_x_coordinate"
    w_nc_fid['x'][...] = np.copy(X_axis)
    dst_var = w_nc_fid.createVariable('y', np.float32, ('y'))
    dst_var.units = "m"
    dst_var.long_name = "y coordinate of projection"
    dst_var.standard_name = "projection_y_coordinate"
    w_nc_fid['y'][...] = np.copy(Y_axis)
    w_nc_var = w_nc_fid.createVariable('lon', np.float64, ('y','x'))
    w_nc_var.long_name = 'longitude of 2D land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_fid.variables['lon'][...] = LON2D
    w_nc_var = w_nc_fid.createVariable('lat', np.float64, ('y','x'))
    w_nc_var.long_name = 'latitude of 2D land gridcell center (GCS_WGS_84), increasing from south to north'
    w_nc_var.units = "degrees_north"
    w_nc_fid.variables['lat'][...] = LAT2D

    # create the gridIDs, lon, and lat variable
    w_nc_fid.createDimension('ni', grid_id_arr.size)
    w_nc_fid.createDimension('nj', 1)
    w_nc_fid.createDimension('nv', 4)

    w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId in daymet coordinate domain'
    w_nc_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells" 
    w_nc_fid.variables['gridID'][...] = grid_id_arr

    w_nc_var = w_nc_fid.createVariable('gridXID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId x in daymeet coordinate domain'
    w_nc_var.decription = "start from #0 at the upper left corner and from west to east of the domain, with gridID=gridXID+gridYID*x_dim" 
    w_nc_fid.variables['gridXID'][...] = grid_xids_arr

    w_nc_var = w_nc_fid.createVariable('gridYID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId y in daymet coordinate domain'
    w_nc_var.decription = "start from #0 at the upper left corner and from north to south of the domain, with gridID=gridXID+gridYID*y_dim" 
    w_nc_fid.variables['gridYID'][...] = grid_yids_arr

    w_nc_var = w_nc_fid.createVariable('xc', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_var.bounds = "xv"
    w_nc_fid.variables['xc'][...] = lon_arr
        
    w_nc_var = w_nc_fid.createVariable('yc', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell center (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"
    w_nc_fid.variables['yc'][...] = lat_arr
        
    # create the XC, YC variable
    w_nc_var = w_nc_fid.createVariable('xc_LCC', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'X of land gridcell center (Lambert Conformal Conic), increasing from west to east'
    w_nc_var.units = "m"
    w_nc_fid.variables['xc_LCC'][...] = geox_arr
        
    w_nc_var = w_nc_fid.createVariable('yc_LCC', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'Y of land gridcell center (Lambert Conformal Conic), decreasing from north to south'
    w_nc_var.units = "m"
    w_nc_fid.variables['yc_LCC'][...] = geoy_arr

    #
    w_nc_var = w_nc_fid.createVariable('xv_LCC', np.float64, ('nj','ni','nv'))
    w_nc_var.long_name = 'X of land gridcell verticles (Lambert Conformal Conic), increasing from west to east'
    w_nc_var.units = "m"
    w_nc_var = w_nc_fid.createVariable('yv_LCC', np.float64, ('nj','ni','nv'))
    w_nc_var.long_name = 'Y of land gridcell verticles (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "m"

    xv_arr,yv_arr = np.asarray([geox_arr-side_km*500.,geoy_arr-side_km*500.]) #left-lower vertice
    w_nc_fid.variables['xv_LCC'][...,0] = xv_arr
    w_nc_fid.variables['yv_LCC'][...,0] = yv_arr
    xv_arr,yv_arr = np.asarray([geox_arr+side_km*500.,geoy_arr-side_km*500.]) #right-lower vertice
    w_nc_fid.variables['xv_LCC'][...,1] = xv_arr
    w_nc_fid.variables['yv_LCC'][...,1] = yv_arr
    xv_arr,yv_arr = np.asarray([geox_arr+side_km*500.,geoy_arr+side_km*500.]) #right-upper vertice
    w_nc_fid.variables['xv_LCC'][...,2] = xv_arr
    w_nc_fid.variables['yv_LCC'][...,2] = yv_arr
    xv_arr,yv_arr = np.asarray([geox_arr-side_km*500.,geoy_arr+side_km*500.]) #left-upper vertice
    w_nc_fid.variables['xv_LCC'][...,3] = xv_arr
    w_nc_fid.variables['yv_LCC'][...,3] = yv_arr

    w_nc_var = w_nc_fid.createVariable('xv', np.float64, ('nj','ni','nv'))
    w_nc_var.long_name = 'longitude of land gridcell verticles (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_var = w_nc_fid.createVariable('yv', np.float64, ('nj','ni','nv'))
    w_nc_var.long_name = 'latitude of land gridcell verticles (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"
    # ideally, [lontv,latv] will be transformed from above [xv_LCC, yv_LCC] 
    # but there is data mismatch issue in original dataset. 
    # the following is a temporay adjustment, which may have issue of grid-bounds
    lonx, laty = Txy2lonlat.transform(geox_arr, geoy_arr)
    xdiff = lon_arr - lonx
    ydiff = lat_arr - laty
    lonx, laty = Txy2lonlat.transform(w_nc_fid.variables['xv_LCC'][...,0], 
                                      w_nc_fid.variables['yv_LCC'][...,0])    
    w_nc_fid.variables['xv'][...,0] = lonx + xdiff # this adjustment likely makes vertices are not shared by neighboring grids
    w_nc_fid.variables['yv'][...,0] = laty + ydiff
    lonx, laty = Txy2lonlat.transform(w_nc_fid.variables['xv_LCC'][...,1], 
                                      w_nc_fid.variables['yv_LCC'][...,1])    
    w_nc_fid.variables['xv'][...,1] = lonx + xdiff
    w_nc_fid.variables['yv'][...,1] = laty + ydiff
    lonx, laty = Txy2lonlat.transform(w_nc_fid.variables['xv_LCC'][...,2], 
                                      w_nc_fid.variables['yv_LCC'][...,2])    
    w_nc_fid.variables['xv'][...,2] = lonx + xdiff
    w_nc_fid.variables['yv'][...,2] = laty + ydiff
    lonx, laty = Txy2lonlat.transform(w_nc_fid.variables['xv_LCC'][...,3], 
                                      w_nc_fid.variables['yv_LCC'][...,3])    
    w_nc_fid.variables['xv'][...,3] = lonx + xdiff
    w_nc_fid.variables['yv'][...,3] = laty + ydiff

    #
    w_nc_var = w_nc_fid.createVariable('area', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "radian^2"
    w_nc_fid.variables['area'][...] = area_arr

    w_nc_var = w_nc_fid.createVariable('area_LCC', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells (Lambert Conformal Conic)'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "km^2"
    w_nc_fid.variables['area_LCC'][...] = area_km2

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'mask of land gridcells (1 means land)'
    w_nc_var.units = "unitless"
    w_nc_fid.variables['mask'][...] = mask_arr

    w_nc_var = w_nc_fid.createVariable('frac', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'fraction of land gridcell that is active'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "unitless"
    w_nc_fid.variables['frac'][...] = landfrac_arr

    w_nc_var = w_nc_fid.createVariable('lambert_conformal_conic', np.short)
    w_nc_var.grid_mapping_name = "lambert_conformal_conic"
    w_nc_var.longitude_of_central_meridian = -100.
    w_nc_var.latitude_of_projection_origin = 42.5
    w_nc_var.false_easting = 0.
    w_nc_var.false_northing = 0.
    w_nc_var.standard_parallel = 25., 60.
    w_nc_var.semi_major_axis = 6378137.
    w_nc_var.inverse_flattening = 298.257223563
    w_nc_var.proj4_str = geoxy_proj_str
    w_nc_var.epsg = epsg_code

    w_nc_fid.close()  # close the new file  


''' 
 standardilized ELM domain.nc write
'''
def domain_ncwrite(elmdomain_data, WRITE2D=True, ncfile='domain.nc', coord_system=True):
        
    # (optional) full 2D spatial extent
    if coord_system and 'X_axis' in elmdomain_data.keys():
        X_axis = elmdomain_data['X_axis']  # either longitude or projected geox, indiced as gridXID
        Y_axis = elmdomain_data['Y_axis']  # either latitude or projected geoy, indiced as gridYID
        XX = elmdomain_data['XX']  # longitude[Y_axis,X_axis], indiced as gridID
        YY = elmdomain_data['YY']  # latitude[Y_axis,X_axis], indiced as gridID
              
        # (optional) the following is for map conversion of 1D <==> 2D   
        grid_id_arr = elmdomain_data['gridID']
        grid_xid_arr = elmdomain_data['gridXID']
        grid_yid_arr = elmdomain_data['gridYID']
        
    # standard
    lon_arr = elmdomain_data['xc']
    lat_arr = elmdomain_data['yc']
    xv_arr  = elmdomain_data['xv']
    yv_arr  = elmdomain_data['yv']
    area_arr = elmdomain_data['area']
    mask_arr = elmdomain_data['mask']
    landfrac_arr = elmdomain_data['frac']
  
    # write to nc file
    file_name = ncfile
    print("The domain file to be write is " + file_name)

    # Open a new NetCDF file to write the domain information. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
    if WRITE2D and len(grid_id_arr.shape)>1:
        w_nc_fid.title = '2D domain file for running ELM'
    else:
        w_nc_fid.title = '1D (unstructured) domain file for running ELM'
    w_nc_fid.Conventions = "CF-1.0"
    if coord_system and 'X_axis' in elmdomain_data.keys():
        w_nc_fid.note = "a 2d-grid coordinates system, ('y','x'), added, with 2d lon/lat variables. " + \
                    " With this, GridXID, GridYID, and GridID added for actual grids(nj,ni) " + \
                    " so that transform between 2d and 1d may be possible"

    # Create new dimensions of new coordinate system
    if coord_system and 'X_axis' in elmdomain_data.keys():

        # create the gridIDs, lon, and lat variable
        w_nc_fid.createDimension('x', len(X_axis))
        w_nc_fid.createDimension('y', len(Y_axis))
        w_nc_var = w_nc_fid.createVariable('x', np.float64, ('x'))
        w_nc_var.long_name = 'longitude of x-axis'
        w_nc_var.units = "degrees_east"
        w_nc_var.comment = "could be in unit of projected x such as meters"
        w_nc_fid.variables['x'][...] = X_axis
    
        w_nc_var = w_nc_fid.createVariable('y', np.float64, ('y'))
        w_nc_var.long_name = 'latitude of y-axis'
        w_nc_var.units = "degrees_north"
        w_nc_var.comment = "could be in unit of projected y such as meters"
        w_nc_fid.variables['y'][...] = Y_axis
    
        w_nc_var = w_nc_fid.createVariable('lon', np.float64, ('y','x'))
        w_nc_var.long_name = 'longitude of 2D land gridcell center (GCS_WGS_84), increasing from west to east'
        w_nc_var.units = "degrees_east"
        w_nc_fid.variables['lon'][...] = XX
    
        w_nc_var = w_nc_fid.createVariable('lat', np.float64, ('y','x'))
        w_nc_var.long_name = 'latitude of 2D land gridcell center (GCS_WGS_84), increasing from south to north'
        w_nc_var.units = "degrees_north"
        w_nc_fid.variables['lat'][...] = YY

    if WRITE2D and (not 1 in lon_arr.shape):
        w_nc_fid.createDimension('ni', lon_arr.shape[1])
        w_nc_fid.createDimension('nj', lon_arr.shape[0])
        w_nc_fid.createDimension('nv', 4)
    else:
        w_nc_fid.createDimension('ni', lon_arr.size)
        w_nc_fid.createDimension('nj', 1)
        w_nc_fid.createDimension('nv', 4)

    if coord_system and 'X_axis' in elmdomain_data.keys():
        # for 2D <--> 1D indices
        w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('nj','ni'))
        w_nc_var.long_name = '1d unique gridId in boxed range of coordinates y, x'
        w_nc_var.decription = "start from #0 at the lower left corner of the domain, covering all land and ocean gridcells" 
        w_nc_fid.variables['gridID'][...] = grid_id_arr
    
        w_nc_var = w_nc_fid.createVariable('gridXID', np.int32, ('nj','ni'))
        w_nc_var.long_name = 'gridId x in boxed range of coordinate x'
        w_nc_var.decription = "start from #0 at the lower left corner and from west to east of the domain, with gridID=gridXID+gridYID*len(y)" 
        w_nc_fid.variables['gridXID'][...] = grid_xid_arr
    
        w_nc_var = w_nc_fid.createVariable('gridYID', np.int32, ('nj','ni'))
        w_nc_var.long_name = 'gridId y in boxed range of coordinate y'
        w_nc_var.decription = "start from #0 at the lower left corner and from south to north of the domain, with gridID=gridXID+gridYID*len(y)" 
        w_nc_fid.variables['gridYID'][...] = grid_yid_arr

    #
    w_nc_var = w_nc_fid.createVariable('xc', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell center'
    w_nc_var.units = "degrees_east"
    w_nc_var.bounds = "xv"
    w_nc_var.comment = 'by GCS_WGS_84, increasing from west to east'
    w_nc_fid.variables['xc'][...] = lon_arr
        
    w_nc_var = w_nc_fid.createVariable('yc', np.float64, ('nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell center'
    w_nc_var.units = "degrees_north"
    w_nc_var.comment = 'by GCS_WGS_84, decreasing from north to south'
    w_nc_fid.variables['yc'][...] = lat_arr
        
    # create the XC, YC variable        
    #
    w_nc_var = w_nc_fid.createVariable('xv', np.float64, ('nj','ni','nv'))
    w_nc_var.long_name = 'longitude of land gridcell verticles'
    w_nc_var.units = "degrees_east"
    w_nc_var.comment = 'by GCS_WGS_84, increasing from west (-180) to east (180), vertices ordering anti-clock from left-low corner'
    w_nc_var = w_nc_fid.createVariable('yv', np.float64, ('nj','ni','nv'))
    w_nc_var.long_name = 'latitude of land gridcell verticles'
    w_nc_var.units = "degrees_north"
    w_nc_var.comment = 'by GCS_WGS_84, decreasing from north to south, vertices ordering anti-clock from left-low corner'

    w_nc_fid.variables['xv'][...]= xv_arr
    w_nc_fid.variables['yv'][...]= yv_arr

    w_nc_var = w_nc_fid.createVariable('area', np.float64, ('nj','ni'))
    w_nc_var.long_name = "area of grid cell in radians squared" ;
    w_nc_var.coordinates = "xc yc" ;
    w_nc_var.units = "radian2" ;
    w_nc_fid.variables['area'][...] = area_arr

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('nj','ni'))
    w_nc_var.long_name = "domain mask" ;
    w_nc_var.note = "unitless" ;
    w_nc_var.coordinates = "xc yc" ;
    w_nc_var.comment = "0 value indicates cell is not active" ;
    w_nc_fid.variables['mask'][...] = mask_arr

    w_nc_var = w_nc_fid.createVariable('frac', np.float64, ('nj','ni'))
    w_nc_var.long_name = "fraction of grid cell that is active" ;
    w_nc_var.coordinates = "xc yc" ;
    w_nc_var.note = "unitless" ;
    w_nc_var.filter1 = "error if frac> 1.0+eps or frac < 0.0-eps; eps = 0.1000000E-11" ;
    w_nc_var.filter2 = "limit frac to [fminval,fmaxval]; fminval= 0.1000000E-02 fmaxval=  1.000000" ;
    w_nc_fid.variables['frac'][...] = landfrac_arr

    w_nc_fid.close()  # close the new file 

def domain_remeshbynewmesh(input_pathfile='./domain.nc', newmesh='', \
                     LONGXY360=False, ncwrite_coords=True):
    '''
    re-meshing ELM domain grids, by providing a new mesh
        input_pathfile - source ELM domain nc file, with required variables of 'xc', 'yc', 'xv', 'yv', 'mask','frac'
                         grid system is lat/lon with dimension names of (nj, ni) of nj for latitude and ni for longitude
        LONGXY360      - true if longitude from 0-360; otherwise -180 ~ 180
        ncwrite_coords - when output type is 'domain_ncwrite', write coordinates and info OR not
    '''
    
    newmesh = 'projection_elm_mappings.txt'
    if newmesh.endswith('.txt'):
        with open(newmesh) as f:
            # remove header txts
            next(f)
            try:
                data = [x.strip().split() for x in f]
            except Exception as e:
                print(e)
                print('Error in reading - '+newmesh)
                os.sys.exit(-1)
        data = np.asarray(data,float)
        lon=data[:,0]
        lat=data[:,1]
        geox=data[:,2]
        geoy=data[:,3]
        xidx=np.asanyarray(data[:,4],int)-1  # in mappings.txt, those index are 1-based
        yidx=np.asanyarray(data[:,5],int)-1
        gidx=np.asanyarray(data[:,6],int)-1

    # re-mesh to make lon/lat axised 
    [jdx, j] = np.unique(yidx, return_index=True)
    jj = np.argsort(jdx)
    jdx = j[jj]
    yy = geoy[jdx]

    [idx, i] = np.unique(xidx, return_index=True)
    ii = np.argsort(idx)
    idx = i[ii]
    xx = geox[idx]
   
    
    # 
    # Open the source domain or surface nc file
    # 
    # output arrays
        
    file_name = input_pathfile+'_remeshed'
    print("new domain file is " + file_name)


def domain_remeshbycentroid(input_pathfile='./share/domains/domain.clm/domain.nc', 
                     LONGXY360=False, edge_wider=1.0, out2d=False, ncwrite_coords=True):
    '''
    re-meshing ELM domain grids, by providing centroids
        input_pathfile - source ELM domain nc file, with required variables of 'xc', 'yc', 'xv', 'yv', 'mask','frac'
                         grid system is lat/lon with dimension names of (nj, ni) of nj for latitude and ni for longitude
        LONGXY360      - true if longitude from 0-360; otherwise -180 ~ 180
        edge_wider     - in case need to have a user-defined wider edges of domain
        out2d          - output data in 2D or flatten (so-called unstructured)
        ncwrite_coords - when output type is 'domain_ncwrite', write coordinates and info OR not
    '''
    
    # 
    # Open the source domain or surface nc file, 
    # from which trunked a box with limits of [mask_yc,mask_xc]
    src = nc.Dataset(input_pathfile, 'r')
    # points = [y,x] coordinates for src's grid
    if 'domain' in input_pathfile:
        # directly from an ELM domain nc file
        src_yc = src.variables['yc'][...]
        src_xc = src.variables['xc'][...]
        src_xv = src.variables['xv'][...]
        src_yv = src.variables['yv'][...]
        #nv dim position
        nvdim = src.variables['xv'].dimensions.index('nv')

        src_mask = src.variables['mask'][...]
        src_area = src.variables['area'][...]
        src_landfrac = src.variables['frac'][...]
    
    elif 'surfdata' in input_pathfile:
        # redo domain from an ELM surface data file
        src_yc = src.variables['LONGXY'][...]
        src_xc = src.variables['LATIXY'][...]
        
        src_mask = src.variables['LANDFRAC_MASK'][...]
        src_area = src.variables['AREA'][...] # needs to convert from km^2 to rad^2
        src_landfrac = src.variables['LANDFRAC_PFT'][...]
    
    else:
        print('Source domain file should be either a "domain*.nc" or a "surfdata*.nc" ')
        os.sys.exit(-1)

    if not LONGXY360:
        src_xc[src_xc>=180.0] = src_xc[src_xc>=180.0]-360.0
        src_xv[src_xv>=180.0] = src_xv[src_xv>=180.0]-360.0
    else:
        src_xc[src_xc<0.0]=360+src_xc[src_xc<0.0]
        src_xv[src_xv<0.0]=360+src_xv[src_xv<0.0]
        
    
    landmask = src_mask
    landfrac = src_landfrac
    area = src_area
    xc = src_xc
    yc = src_yc
    if nvdim==0: # nv dim may be in the first or the last
        yv1 = src_yv[0,...]; xv1 = src_xv[0,...]
        yv2 = src_yv[1,...]; xv2 = src_xv[1,...]
        yv3 = src_yv[2,...]; xv3 = src_xv[2,...]
        yv4 = src_yv[3,...]; xv4 = src_xv[3,...]
    else:
        yv1 = src_yv[...,0]; xv1 = src_xv[...,0]
        yv2 = src_yv[...,1]; xv2 = src_xv[...,1]
        yv3 = src_yv[...,2]; xv3 = src_xv[...,2]
        yv4 = src_yv[...,3]; xv4 = src_xv[...,3]        
    
    # re-meshing, by new mesh-griding with xc/yc as centroids 
    X_axis = xc
    Y_axis = yc
    [X_axis, i] = np.unique(X_axis, return_inverse=True)
    [Y_axis, j] = np.unique(Y_axis, return_inverse=True)
    YY, XX = np.meshgrid(Y_axis, X_axis, indexing='ij')
    xid = np.indices(XX.shape)[1]
    yid = np.indices(XX.shape)[0]
    xyid = np.indices(XX.flatten().shape)[0]
    xyid = np.reshape(xyid, xid.shape)
    
    # indices of original grids in new-XX/YY coordinate system 
    if (1 not in src_xc.shape):
        # structured grids yc[nj,ni], xc[nj,ni]      
        # new indices of [yc,xc] in grid-mesh YY/XX
        # note: here don't override xc,yc, so xc=X_axis[i], yc=Y_axis[j]
        polygonized = np.isin(XX, X_axis[i]) & np.isin(YY, Y_axis[j])
        ji = np.nonzero(polygonized) 
        grid_xid = xid[ji]     
        grid_yid = yid[ji]    
        grid_id = xyid[ji]  

    else:
        # unstructured  yc[1,ni], xc[1,ni]
        # new indices of paired [yc,xc] in grid-mesh YY/XX
        grid_xid = xid[j,i]     
        grid_yid = yid[j,i] 
        grid_id = xyid[j,i]
    #
    
    # new grid vertices, middle of neiboured centroids
    if len(X_axis)>1:
        xdiff = np.diff(X_axis)/2.0
    else:
        xdiff = (max(src_xv)-min(src_xv))/2.0
    xdiff = np.insert(xdiff,0, xdiff[0]*edge_wider)
    xv=np.append(X_axis-xdiff,X_axis[-1]+xdiff[-1]*edge_wider)  # new X bounds for grids
    if len(Y_axis)>1:
        ydiff = np.diff(Y_axis)/2.0
    else:
        ydiff = (max(src_yv)-min(src_yv))/2.0        
    ydiff = np.insert(ydiff,0,ydiff[0]*edge_wider)
    yv=np.append(Y_axis-ydiff,Y_axis[-1]+ydiff[-1]*edge_wider)  # new Y bounds for grids
    # re-fill vertices
    yvv, xvv = np.meshgrid(yv[:-1], xv[:-1], indexing='ij')  # left-lower vertice
    xv1 = xvv[(grid_yid,grid_xid)]; yv1 = yvv[(grid_yid,grid_xid)]
    yvv, xvv = np.meshgrid(yv[:-1], xv[1:], indexing='ij')   # right-lower vertice
    xv2 = xvv[(grid_yid,grid_xid)]; yv2 = yvv[(grid_yid,grid_xid)]
    yvv, xvv = np.meshgrid(yv[1:], xv[1:], indexing='ij')    # right-upper vertice
    xv3 = xvv[(grid_yid,grid_xid)]; yv3 = yvv[(grid_yid,grid_xid)]
    yvv, xvv = np.meshgrid(yv[1:], xv[:-1], indexing='ij')   # left-upper vertice
    xv4 = xvv[(grid_yid,grid_xid)]; yv4 = yvv[(grid_yid,grid_xid)]

    # to 2D but masked    
    if out2d:
        grid_id_arr = np.reshape(grid_id,xyid.shape)
        grid_xid_arr = np.reshape(grid_xid,xyid.shape)
        grid_yid_arr = np.reshape(grid_yid,xyid.shape)
        
        xc_arr = np.reshape(xc,xyid.shape)
        yc_arr = np.reshape(yc,xyid.shape)
        xv_arr = np.stack((np.reshape(xv1,xyid.shape), \
                           np.reshape(xv2,xyid.shape), \
                           np.reshape(xv3,xyid.shape), \
                           np.reshape(xv4,xyid.shape)), axis=2)
        yv_arr = np.stack((np.reshape(yv1,xyid.shape), \
                           np.reshape(yv2,xyid.shape), \
                           np.reshape(yv3,xyid.shape), \
                           np.reshape(yv4,xyid.shape)), axis=2)
        area_arr = np.reshape(area,xyid.shape)
        mask_arr = np.reshape(landmask,xyid.shape)
        landfrac_arr = np.reshape(landfrac,xyid.shape)
    else:
    # to 1D only masked
        masked = np.where((landmask==1))
        grid_id_arr = grid_id[masked]
        grid_xid_arr = grid_xid[masked]
        grid_yid_arr = grid_yid[masked]
        xc_arr = xc[masked]
        yc_arr = yc[masked]
        xv_arr = np.stack((xv1[masked],xv2[masked],xv3[masked],xv4[masked]), axis=1)
        yv_arr = np.stack((yv1[masked],yv2[masked],yv3[masked],yv4[masked]), axis=1)
        area_arr = area[masked]
        mask_arr = landmask[masked]
        landfrac_arr = landfrac[masked]
 
    # 
    # output arrays
    subdomain = {}
    subdomain['X_axis'] = X_axis
    subdomain['Y_axis'] = Y_axis
    subdomain['XX'] = XX
    subdomain['YY'] = YY
            
    subdomain['gridID']  = grid_id_arr
    subdomain['gridXID'] = grid_xid_arr
    subdomain['gridYID'] = grid_yid_arr
    subdomain['xc'] = xc_arr
    subdomain['yc'] = yc_arr
    subdomain['xv'] = xv_arr
    subdomain['yv'] = yv_arr
    subdomain['area'] = area_arr
    subdomain['mask'] = mask_arr
    subdomain['frac'] = landfrac_arr
        
    file_name = input_pathfile+'_remeshed'
    print("new domain file is " + file_name)
    domain_ncwrite(subdomain, WRITE2D=out2d, ncfile=file_name, \
                       coord_system=ncwrite_coords)

def domain_remask(input_pathfile='./share/domains/domain.lnd.r05_RRSwISC6to18E3r5.240328.nc', 
                     masked_pts={}, reorder_src=False, keep_duplicated=False, \
                     unlimit_xmin=False, unlimit_xmax=False, unlimit_ymin=False, unlimit_ymax=False,\
                     out2d=False, out='subdomain', ncwrite_coords=True):
    '''
    re-masking ELM domain grids, by providing centroids which also may be masked
        input_pathfile - source ELM domain nc file, with required variables of 'xc', 'yc', 'xv', 'yv', 'mask','frac'
                         grid system is lat/lon with dimension names of (nj, ni) of nj for latitude and ni for longitude
        masked_pts     - list of np arrays of 'xc', 'yc', 'mask' of user-provided grid-centroids
        reorder_src    - trunked source domain will re-order following that in masked_pts, 
                        otherwise NOT except 'mask' and maybe 'frac'
        out2d          - output data in 2D or flatten (so-called unstructured)
        unlimit_x(y)min(max) - x/y axis bounds are from source domain, when True; otherwise, from masked_pts
        out            - output type: 'subdomain' for outputing a domain-style dataset,
                                      'mask' for outputing a np.where style tuple of indices of grids, with a new mask
                                      'domain_ncwrite' for writing out a domain.nc file. 
        ncwrite_coords - when output type is 'domain_ncwrite', write coordinates and info OR not
    '''

    from pytools.commons_utils.gridlocator import grids_nearest_or_within
    
    srcnc = nc.Dataset(input_pathfile, 'r')
    src_grids = {}
    src_grids['xc'] = srcnc['xc'][...]
    src_grids['yc'] = srcnc['yc'][...]
    if 'xv' in srcnc.variables.keys() \
        and 'yv' in srcnc.variables.keys():
        src_grids['xv'] = srcnc['xv'][...]
        src_grids['yv'] = srcnc['yv'][...]
        remask_approach = 'within'
    else:
        remask_approach = 'nearest'
    
    subdomain, boxed_idx, grids_uid, grids_maskptsid = grids_nearest_or_within( \
                    src_grids=src_grids, masked_pts=masked_pts, \
                    remask_approach = remask_approach, \
                    unlimit_xmin=unlimit_xmin, unlimit_xmax=unlimit_xmax, \
                    unlimit_ymin=unlimit_ymin, unlimit_ymax=unlimit_ymax,\
                    out2d=out2d, reorder_src=reorder_src, keep_duplicated=keep_duplicated)
     
    if reorder_src and keep_duplicated:
        # in this case, sub-domain's xc/yc pair should be using those of masked_pts rather than of src_grids
        if 'xc' in masked_pts.keys(): subdomain['xc'] = masked_pts['xc']
        if 'yc' in masked_pts.keys():subdomain['yc'] = masked_pts['yc']
        if 'xv' in masked_pts.keys():subdomain['xv'] = masked_pts['xv']
        if 'yv' in masked_pts.keys():subdomain['yv'] = masked_pts['yv']
                

    # re-do frac of landed mask, if option ON, i.e. pts in a source grid are km2 of land
    landfrac = srcnc['frac'][...]
    if 'km2perpt' in masked_pts.keys() and not reorder_src:
        km2perpt = masked_pts['km2perpt']
    
        # grid area in km2, assuming a rectangle grid (TODO - refining this later)
        yc = srcnc['yc'][...].flatten()
        kmratio_lon2lat = np.cos(np.radians(yc))
        re_km = 6370.997
        xv = subdomain['xv']; yv = subdomain['yv']
        xv0=xv[...,0];xv1=xv[...,1]; xv2=xv[...,2];xv3=xv[...,3]
        yv0=yv[...,0];yv1=yv[...,1]; yv2=yv[...,2];yv3=yv[...,3]
        xside_km = (abs(xv0-xv1)+abs(xv2-xv3))/2.0*(math.pi*re_km/180.0)   # 
        yside_km = (abs(yv0-yv2)+abs(yv1-yv3))/2.0*(math.pi*re_km/180.0*kmratio_lon2lat)
        area_km2 = xside_km*yside_km     
    
        frac = landfrac.flatten()
        for i in range(len(grids_uid)):
            igrid = grids_uid[i]
            pts_counts = len(grids_maskptsid[igrid])
            frac[igrid] = pts_counts*km2perpt/area_km2[igrid]
        landfrac = frac.reshape(landfrac.shape)
        landfrac[np.where(landfrac>1.0)]=1.0
    
    # truncked landmask, fraction, and area 
    landmask = srcnc['mask'][...][boxed_idx]
    landmask[...] = 0 # remasked
    landmask[np.where(subdomain['remasked'])]= 1

    landfrac = landfrac[...][boxed_idx]
    # in case landmask changed to 0, i.e. NOT a land cell, better reassign its fraction to 0
    # but not do so on 'area' which are for a whole grid
    landfrac[landmask==0] = 0.0
    
    landarea=srcnc['area'][...][boxed_idx]  
    srcnc.close()
    
    # trunked domain with updated land mask and land fraction only
    # and exit
    if out == 'mask':
        return boxed_idx, landmask, landfrac, subdomain['xc'], subdomain['yc']

    # 
    # additional subdomain variables
    subdomain['area'] = landarea
    subdomain['mask'] = landmask
    subdomain['frac'] = landfrac
        
    if out == 'subdomain':
        return subdomain
        # a new domain, may or may not same as boxed-truncked
    elif out == 'domain_ncwrite':
        file_name = input_pathfile+'_remasked'
        print("new domain file is " + file_name)            
                 
        domain_ncwrite(subdomain, WRITE2D=out2d, ncfile=file_name, \
                       coord_system=ncwrite_coords)
    #          
#

''' 
 subset an ELM domain.nc by another user provided domain (but only needs: xc, yc, mask)
'''

def domain_subsetbymaskncf(srcdomain_pathfile='./share/domains/domain.lnd.r05_RRSwISC6to18E3r5.240328.nc', \
                        maskncf='./share/domains/domain.clm/domain.lnd.pan-arctic_CAVM.1km.1d.c241018.nc', \
                        maskncv='mask', masknc_area=-999.99, reorder_src=False, keep_duplicated=False, \
                        unlimit_xmin=False, unlimit_xmax=False, unlimit_ymin=False, unlimit_ymax=False,\
                        out2D=True, out_type='subdomain'):
    
    # user provided mask file, which may or may not in same resolution as source domain,
    #                          and only needed are: xc, yc, mask
    mask_new = {}

    mask_f = nc.Dataset(maskncf,'r')
    mask_v = mask_f[maskncv][...]
    mask_checked = np.where(mask_v==1)
    
    mask_new['xc']= mask_f['xc'][...][mask_checked] # this will flatten xc/yc, if in 2D
    mask_new['yc']= mask_f['yc'][...][mask_checked]
    mask_new['mask']= mask_v[mask_checked]
    if masknc_area != -999.99: 
        # if need to convert new masked grid land area or fraction
        # the unit is km^2 per points included in the source domain
        mask_new['km2perpt'] = masknc_area
    
    #
    if out_type == 'subdomain':
        return domain_remask(input_pathfile=srcdomain_pathfile, \
                      masked_pts=mask_new, reorder_src=reorder_src, \
                      keep_duplicated=keep_duplicated, \
                      unlimit_xmin=unlimit_xmin, unlimit_xmax=unlimit_xmax, \
                      unlimit_ymin=unlimit_ymin, unlimit_ymax=unlimit_ymax, \
                      out2d=out2D, out=out_type)
    elif out_type == 'domain_ncwrite':
        domain_remask(input_pathfile=srcdomain_pathfile, \
                      masked_pts=mask_new, reorder_src=reorder_src, \
                      keep_duplicated=keep_duplicated, \
                      unlimit_xmin=unlimit_xmin, unlimit_xmax=unlimit_xmax, \
                      unlimit_ymin=unlimit_ymin, unlimit_ymax=unlimit_ymax, \
                      out2d=out2D, out=out_type)
    elif out_type == 'mask':
        return domain_remask(input_pathfile=srcdomain_pathfile, \
                      masked_pts=mask_new, reorder_src=reorder_src, \
                      keep_duplicated=keep_duplicated, \
                      unlimit_xmin=unlimit_xmin, unlimit_xmax=unlimit_xmax, \
                      unlimit_ymin=unlimit_ymin, unlimit_ymax=unlimit_ymax, \
                      out2d=out2D, out=out_type)
            
    # 
#--------------------------------------------------------------------------------------------------------

''' 
subset an ELM domain.nc by user provided paired latlons[[lats],[lons]] (but only needs: xc, yc)
'''
def domain_subsetbylatlon(srcdomain_pathfile='./share/domains/domain.lnd.r05_RRSwISC6to18E3r5.240328.nc', \
                          latlons=np.empty((0,0)), reorder_src=False, keep_duplicated=False, \
                          unlimit_xmin=False, unlimit_xmax=False, unlimit_ymin=False, unlimit_ymax=False, \
                          out2D=False, out_type='subdomain'):

    
    # new masked domain pts in np.array [[lats],[lons]]
    if (latlons.shape[0]<2):
        print('latlons must have paired location points: y/x or latitude/longidue',latlons.shape[0])
        return
    
    mask_new = {}
    mask_new['yc'] = latlons[0] # y or latitudes
    mask_new['xc'] = latlons[1] # x or longitudes
    mask_new['mask']= np.ones_like(latlons[0]) # assume all pts are masked land grids    
    #
    if out_type == 'subdomain':
        return domain_remask(input_pathfile=srcdomain_pathfile, \
                      masked_pts=mask_new, reorder_src=reorder_src, keep_duplicated=keep_duplicated, \
                      out2d=out2D, out=out_type)
    elif out_type == 'domain_ncwrite':
        domain_remask(input_pathfile=srcdomain_pathfile, \
                      masked_pts=mask_new, reorder_src=reorder_src, keep_duplicated=keep_duplicated, \
                      unlimit_xmin=unlimit_xmin, unlimit_xmax=unlimit_xmax, \
                      unlimit_ymin=unlimit_ymin, unlimit_ymax=unlimit_ymax,\
                      out2d=out2D, out=out_type)
    elif out_type == 'mask':
        return domain_remask(input_pathfile=srcdomain_pathfile, \
                      masked_pts=mask_new, reorder_src=reorder_src, keep_duplicated=keep_duplicated, \
                      out2d=out2D, out=out_type)
            
    # 
#--------------------------------------------------------------------------------------------------------

def ncdata_subsetbynpwhereindex(npwhere_indices, npwhere_mask=np.empty((0,0)), npwhere_frac=np.empty((0,0)), \
                        newxc_box=np.empty((0,0)), newyc_box=np.empty((0,0)), reordered_box=False, \
                        srcnc_pathfile='./lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c240308_TOP.nc', \
                        indx_dim=['lsmlat','lsmlon'], indx_dim_flatten=''):
    '''
    npwhere_box      - a tuple of np.where outputs from domain or surfdata, on 'indx_dim'
    npwheremask_box  - mask of npwhere_box, if provided
    indx_dim_flatten - if not empty, it's the new dimension name for flatten 'indx_dim' if 2d 
    '''
    
    # check if empty indices: 
    if len(npwhere_indices[0])>0:
        #        
        #
        ncfilein  = srcnc_pathfile
        ncfileout = srcnc_pathfile.split('/')[-1]
        ncfileout = ncfileout.split('.nc')[0]+'_subset.nc'
                    
        #subset and write to ncfileout
        with nc.Dataset(ncfilein,'r') as src, nc.Dataset(ncfileout, "w") as dst:
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # dimensions for dst
            for dname, dimension in src.dimensions.items():
                len_dimension = len(dimension)
                if dname in indx_dim:
                    # indx_dim are multiple-D, needs to flatten, i.e. forcing to be length of 1 
                    # because input-indices are for points rather than block
                    dim_len=1
                    if indx_dim.index(dname)==len(indx_dim)-1: 
                        dim_len = sum(npwhere_mask>0)
                    
                    if (dim_len>0):
                        len_dimension = dim_len
                    elif (dim_len==0):
                        dimension=None
                    
                    if indx_dim_flatten!='':
                        if len(dst.dimensions.keys())<=0:
                            dst.createDimension(indx_dim_flatten, len_dimension if not dimension.isunlimited() else None)
                        elif indx_dim_flatten not in dst.dimensions.keys():
                            #only need once, i.e. indx_dim will be flatten into index_dim_new
                            dst.createDimension(indx_dim_flatten, len_dimension if not dimension.isunlimited() else None)
                    else:
                        dst.createDimension(dname, len_dimension if not dimension.isunlimited() else None)                        
                    
                else:
                    dst.createDimension(dname, len_dimension if not dimension.isunlimited() else None)
            #   
        
            # all variables data
            for vname, variable in src.variables.items():

                if all(d in variable.dimensions for d in indx_dim):
                    vdim = np.array(variable.dimensions)
                    if indx_dim_flatten!='':
                        vdim = [d.replace(indx_dim[0], indx_dim_flatten) for d in vdim]
                        for i in range(1,len(indx_dim)):
                            vdim.remove(indx_dim[i])
                else:
                    vdim = variable.dimensions
    
                # create variables, but will update its values later 
                # NOTE: here the variable value size updated due to newly-created dimensions above
                dst.createVariable(vname, variable.datatype, vdim)
                # copy variable attributes all at once via dictionary after created
                dst[vname].setncatts(src[vname].__dict__)
                  
                # values
                print(vname)

                # dimensional slicing to extract data and fill into dst
                if all(d in variable.dimensions for d in indx_dim):
                    vd = variable.dimensions
                    ix = range(vd.index(indx_dim[0]),
                               vd.index(indx_dim[-1])+1) # position of dims in 'indx_dim'
                    
                    # build a tuple of indices, with masked only in indx_dim
                    ix_mask = 0
                    # for unstructured surfdata, ix_mask will be fixed
                    if 'nj' in vd and 'ni' in vd:
                        if src.dimensions['nj'].size == 1: ix_mask = 1
                    elif len(npwhere_indices)>1:
                        if 'gridcell' in vdim or 'n' in vdim: ix_mask = 1
                    for i in range(len(vd)):
                        if i==0:
                            newidx = (slice(None),)
                            if (i in ix and i==ix_mask+ix[0]) or \
                                (i in ix and 'gridcell' in vdim):
                                newidx = (npwhere_indices[ix_mask][npwhere_mask>0],)
                                # for structured 2D, ix_mask will move 1 dimension 
                                if 'nj' in vd and 'ni' in vd:
                                    if src.dimensions['nj'].size != 1: ix_mask = ix_mask+1
                                else:
                                    if not 'gridcell' in vd or not 'n' in vd:ix_mask=ix_mask+1                          
          
                        else:
                            if (i in ix and i==ix_mask+ix[0]) or \
                                (i in ix and 'gridcell' in vdim):
                                newidx = newidx+(npwhere_indices[ix_mask][npwhere_mask>0],)
                                # for structured 2D, ix_mask will move 1 dimension 
                                if 'nj' in vd and 'ni' in vd:
                                    if src.dimensions['nj'].size != 1: ix_mask = ix_mask+1
                                else:
                                    if not 'gridcell' in vd or not 'n' in vd:ix_mask=ix_mask+1
                            else:
                                newidx = newidx+(slice(None),)
                                               
                    varvals = src[vname][...][newidx]
                    #
                   
                else:
                    varvals = src[vname][...]
                    #
                #
                if reordered_box:
                    if not np.array_equal(newxc_box,np.empty((0,0))) \
                        and (vname=='LONGXY' or vname=='xc'): 
                        varvals[...] = newxc_box
                    if not np.array_equal(newyc_box,np.empty((0,0))) \
                        and (vname=='LATIXY' or vname=='yc'): 
                        varvals[...] = newyc_box
                    # the following is hard-weired now - TODO from inputs 
                    if 'AREA' in vname:
                        varvals[...] = 1.0 # 1km^2 
                
                #                            
                dst[vname][...] = varvals
                    
            # variable-loop        
                
        
        print('subsetting done successfully!')
        
    else:
        print('NO subsetting done due to indices NOT provided!')
    #
#
#--------------------------------------------------------------------------------------------------------
   
    
#--------------------------------------------------------------------------------------------------------
  
def main():

    import glob
    try:
        from mpi4py import MPI
        HAS_MPI4PY=True
    except ImportError:
        HAS_MPI4PY=False

    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """  
    input_path= '../../surfdata_0.5x0.5_simyr1850_c240308_TOP_cavm1d'
    #input_path= './domain.lnd.pan-arctic_CAVM.0.01deg.1D.c250623'

    # create an unstructured domain from a raster image, e.g. CAVM image of land cover type, including veg
    #domain_unstructured_fromraster_cavm()
    #domain_unstructured_fromraster_cavm(output_pathfile='./domain.lnd.pan-arctic_CAVM.0.01deg.1D.c250623.nc', \
    #                               rasterfile='./raster_cavm_v1_01d.tif')
    #domain_unstructured_fromraster_cavm(output_pathfile='./domain.lnd.original.1D.c250624_TFSarcticpfts.nc', \
    #                               rasterfile='./deciduous_shrub_toolik.tif')
    #return # for only do domain.nc writing
    
    # create an unstructured domain from grids of daymet tile
    #domain_unstructured_fromdaymet(tileinfo_ncfile='./daymet_tiles.nc')
    #return # for only do domain.nc writing

    # vertices redo, either incorrect (e.g. from daymet), or unknown (e.g. from surfdata)
    #domain_remeshbycentroid(input_pathfile='./mysurfnc_dir/domain.nc', 
    #                 LONGXY360=False, edge_wider=1.0, out2d=False, ncwrite_coords=False)
    #return # for only do domain.nc writing
    
    
    #-------------
    # within input_path, directory and/or file header, dataset's domain file to extract subset of data
    inputdomain_ncfile = '../../domain.lnd.r05_RRSwISC6to18E3r5.240328_cavm1d.nc'
    #inputdomain_ncfile= './domain.lnd.pan-arctic_CAVM.0.01deg.1D.c250623.nc'
    #inputdomain_ncfile = './TFSarcticpfts/domain.lnd.pan-arctic_CAVM.1km.1D.c241018_TFSarcticpfts.nc'
    output_path = './'
    SUBDOMAIN_ONLY = False
    SUBDOMAIN_REORDER = True  # True: masked file re-ordered by userdomain below, otherwise only mask and trunck
    KEEP_DUPLICATED = True
    NC2D = False # for domain.nc writing only, so doesn't matter if SUBDOMAIN_ONLY is False
    
    
    # -----------
    # user-provided  lat/lon of domain or sites or geotiff  to extract subset data for
    #userdomain = './domain.lnd.pan-arctic_CAVM.0.01deg.1D.c250623.nc'
    userdomain = './domain.lnd.0.01deg.1D.c250708_TFSarcticpfts.nc'
    km2perpt = -999.99 #1.0 # user-grid area in km^2

    '''
    #userdomain = '/Users/f9y/mygithub/E3SM_REPOS/pt-e3sm-inputdata/atm/datm7/'+ \
    #             'atm_forcing.datm7.GSWP3.0.5d.v2.c180716_NGEE-Grid/'+ \
    #             'atm_forcing.datm7.GSWP3.0.5d.v2.c180716_ngee-TFS-Grid/info_TFS_meq2_sites.txt'
    
    userdomain='./info_TFS_meq2_sites.txt'
    #userdomain = './zone_mappings.txt'
    #km2perpt = 1.0 # user-grid area in km^2
    # e.g. 
    #  [f9y@baseline-login3 atm_forcing.datm7.GSWP3.0.5d.v2.c180716_ngee4]$ cat info_TFS_meq2_sites.txt 
    #    site_name lat lon
    #    MEQ2-MAT 68.6611 -149.3705
    #    MEQ2-DAT 68.607947 -149.401596
    #    MEQ2-PF 68.579315 -149.442279
    #    MEQ2-ST 68.606131 -149.505794
    '''
    lats=[]; lons=[]
    '''
    with open(userdomain) as f:
        dtxt=f.readlines()
        
        dtxt=filter(lambda x: x.strip(), dtxt)
        for d in dtxt:
            allgrds=np.array(d.split()[1:3],dtype=float)
            lons.append(allgrds[1])
            lats.append(allgrds[0])
    f.close()
    lons = np.asarray(lons)
    lats = np.asarray(lats)
    '''
    
    '''
    userdomain = './TFSarcticpfts/graminoid_toolik_extent.tif'
    xx, yy, crs_res, crs_wkt, alldata = geotiff2nc(userdomain, outdata=True)
    allyc, allxc = np.meshgrid(yy,xx,indexing='ij')
    lons = allxc[~alldata[0].mask]
    lats = allyc[~alldata[0].mask]
    '''
    
    
    # search for source-domain indices masked by user-domain
    if SUBDOMAIN_ONLY:
        # only write new domain file
        if userdomain.endswith('.nc'):
            domain_subsetbymaskncf( \
            srcdomain_pathfile=inputdomain_ncfile, \
            maskncf=userdomain, masknc_area=km2perpt, reorder_src=SUBDOMAIN_REORDER, \
            unlimit_xmin=True, unlimit_xmax=True, unlimit_ymax=True, \
            out2D=NC2D, out_type='domain_ncwrite')
        
        elif userdomain.endswith('.txt') or userdomain.endswith('.tif'):
            domain_subsetbylatlon( \
            srcdomain_pathfile=inputdomain_ncfile,  \
            latlons=np.asarray([lats, lons]), reorder_src=SUBDOMAIN_REORDER, \
            out2D=NC2D, out_type='domain_ncwrite')
        
        return

    # continue for subsetting, if option is ON
    if userdomain.endswith('.nc'):
        idx_box, newmask_box, newfrac_box, newxc_box, newyc_box = \
            domain_subsetbymaskncf( \
            srcdomain_pathfile=inputdomain_ncfile, \
            maskncf=userdomain, masknc_area=km2perpt, reorder_src=SUBDOMAIN_REORDER, \
            keep_duplicated=KEEP_DUPLICATED, \
            #unlimit_xmin=True, unlimit_xmax=True, unlimit_ymax=True, \
            out2D=NC2D, out_type='mask')
    
    elif userdomain.endswith('.txt') or userdomain.endswith('.tif'):
        idx_box, newmask_box, newfrac_box, newxc_box, newyc_box = \
            domain_subsetbylatlon( \
            srcdomain_pathfile=inputdomain_ncfile,  \
            latlons=np.asarray([lats, lons]), reorder_src=SUBDOMAIN_REORDER, \
            keep_duplicated=KEEP_DUPLICATED, \
            out2D=NC2D, out_type='mask')

    # run with srun

    # all source nc files to be subset
    allncfiles = sorted(glob.glob("%s*.%s" % (input_path,'nc')))
    if HAS_MPI4PY:
        mycomm = MPI.COMM_WORLD
        myrank = mycomm.Get_rank()
        mysize = mycomm.Get_size()
        
        len_total = len(allncfiles)
        len_myrank = int(math.floor(len_total/mysize))
        len_mod = int(math.fmod(len_total,mysize))

        n_myrank = np.full([mysize], int(1)); n_myrank = np.cumsum(n_myrank)*len_myrank
        x_myrank = np.full([mysize], int(0))
        if(len_mod>0): x_myrank[:len_mod] = 1
        n_myrank = n_myrank + np.cumsum(x_myrank) - 1        # ending index, starting 0, for each rank
        n0_myrank = np.hstack((0, n_myrank[0:mysize-1]+1))   # starting index, starting 0, for each rank
    
        #print('on ',myrank, 'indx: ',len_total, n0_myrank[myrank], n_myrank[myrank], \
        #        allncfiles[n0_myrank[myrank]],allncfiles[n_myrank[myrank]])
        allncfiles_byrank = allncfiles[n0_myrank[myrank]:n_myrank[myrank]+1] # [n0:n] IS not inclusive of [n]
    
    else:
        mycomm = 0
        myrank = 0
        mysize = 1
        allncfiles_byrank = allncfiles
    
    # subsetting other datasets
    # (TODO) shall exclude those sources of grid-system NOT consistent with source-domain
    if (myrank>=mysize-3): print(allncfiles_byrank[0], allncfiles_byrank[-1], 'ON: ', myrank)
    for ncfile in allncfiles_byrank:
        ncfile_dims = ['lsmlat', 'lsmlon']
        dims_new =''

        if 'surfdata' in ncfile or 'landuse' in ncfile:
            # the following is for surfdata, standard grid dims are ['lsmlat', 'lsmlon'], 
            # while unstructured dim name is 'gridcell'
            ncfile_dims = ['gridcell']#['lsmlat', 'lsmlon']
            dims_new= 'gridcell' # if want to reduce dimensions to single and given a new dimension
        
        elif 'domain' in ncfile: 
            ncfile_dims = ['nj', 'ni']
            dims_new = '' 
            # if no change of dimensions, but will force dimension len to be 1, 
            # except the last which will be total unmasked grid number
            # (TODO - maybe another option of masked only?)

        elif 'Daymet_ERA5' in ncfile or 'GSWP3_daymet4' in ncfile: 
            ncfile_dims = ['n']
            dims_new = '' 
        
        else:
            print('NO subsetting done for surfdata, domain, or forcing data - file NOT exists')
            return

        ncdata_subsetbynpwhereindex(idx_box, npwhere_mask=newmask_box, npwhere_frac=newfrac_box, \
                                    newxc_box=newxc_box, newyc_box=newyc_box, reordered_box=SUBDOMAIN_REORDER, \
                                    srcnc_pathfile=ncfile, \
                                    indx_dim=ncfile_dims, indx_dim_flatten=dims_new)

        #

    
if __name__ == '__main__':
    main()
    
