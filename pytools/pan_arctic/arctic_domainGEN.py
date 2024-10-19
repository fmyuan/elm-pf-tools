
# utility to generate 1D and 2D domain
# use previous dataPartition module. Code need to be clean up

import os
import math
import netCDF4 as nc
import numpy as np
from pyproj import Proj
from pyproj import Transformer
from pyproj import CRS
from pytools.commons_utils.Modules_netcdf import geotiff2nc
from numpy import meshgrid

def domain_unstructured_fromraster(output_pathfile='./domain.lnd.pan-arctic_CAVM.1km.1d.c241018.nc', rasterfile='./raster_cavm_v1.tif', outdata=False):
    
    # 1000mx1000m grids extracted from CAVM image (geotiff) file
    # e.g. raster_cavm_v1.tif
    #  

    file=rasterfile
    bandinfos={'bands':['grid_code']}
    allxc, allyc, crs_wkt, alldata = geotiff2nc(file, bandinfos, outdata=True)
    alldata = alldata[0] # only 1 banded data

    # truncate non-data points as much as possible
    # cavm: 91 - freshwater, 92 - sea water (ocean), 93 - glacier, 99 - non-arctic
    data = np.ma.masked_equal(alldata,99)  # non-arctic grids
    data = np.ma.masked_equal(data,92)  # sea grids
    nondata = data.mask
    #along x-axis, checking non-data in columns
    ny_nodata = sum(nondata,0)
    for i in range(1,len(allxc)-1):
        if ny_nodata[i]<len(allyc) or ny_nodata[i-1]<len(allyc): break 
        x0 = i
    for i in range(len(allxc)-2,1,-1):
        if ny_nodata[i]<len(allyc) or ny_nodata[i+1]<len(allyc): break 
        x1 = i
    X_axis = allxc[x0:x1]
    data = data[:,x0:x1]   # note: column x0/x1 shall be all non-data, but x0+1/x1-1 not
    
    #along y-axis, checking non-data in rows
    nx_nodata = sum(nondata,1)
    for i in range(1,len(allyc)-1):
        if nx_nodata[i]<len(allxc) or nx_nodata[i-1]<len(allxc): break 
        y0 = i
    for i in range(len(allyc)-2,1,-1):
        if nx_nodata[i]<len(allxc) or nx_nodata[i+1]<len(allxc): break 
        y1 = i
    Y_axis = allyc[y0:y1]      
    data = data[y0:y1,:]   # note: row y0/y1 shall be all non-data, but y0+1/y1-1 not

    XC, YC = np.meshgrid(X_axis, Y_axis)

    # lon/lat from XC/YC
    #Proj4: +proj=laea +lon_0=-180 +lat_0=90 +x_0=0 +y_0=0 +R=6370997 +f=0 +units=m  +no_defs
    geoxy_proj_str = "+proj=laea +lon_0=-180 +lat_0=90  +x_0=0 +y_0=0 +R=6370997 +f=0 +units=m +no_defs"
    geoxyProj = CRS.from_proj4(geoxy_proj_str)
    # EPSG: 4326
    # Proj4: +proj=longlat +datum=WGS84 +no_defs
    lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
    Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)
    #Tlonlat2xy = Transformer.from_proj(lonlatProj, geoxyProj, always_xy=True)


    lon,lat = Txy2lonlat.transform(XC,YC)

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
    
    # area in km2 --> in arcrad2
    area_km2 = 1.0 # this is by default
    side_km = math.sqrt(float(area_km2))
    lat[lat==90.0]=lat[lat==90.0]-0.00001
    lat[lat==-90.0]=lat[lat==-90.0]-0.00001
    kmratio_lon2lat = np.cos(np.radians(lat))
    re_km = 6370.997
    yscalar = side_km/(math.pi*re_km/180.0)
    xscalar = side_km/(math.pi*re_km/180.0*kmratio_lon2lat)
    area = xscalar*yscalar
    
    # 2d --> 1d, with masked only
    grid_id_arr = grid_ids[masked]
    grid_xids_arr = grid_xids[masked]
    grid_yids_arr = grid_yids[masked]
    mask_arr = landmask[masked]
    landfrac_arr = landfrac[masked]
    area_arr = area[masked]
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
    x_dim = w_nc_fid.createDimension('x', total_cols)
    y_dim = w_nc_fid.createDimension('y', total_rows)
    dst_var = w_nc_fid.createVariable('x', np.float32, ('x'))
    dst_var.units = "m"
    dst_var.long_name = "x coordinate of projection"
    dst_var.standard_name = "projection_x_coordinate"
    w_nc_fid['x'][...] = np.copy(XC[0,:])
    dst_var = w_nc_fid.createVariable('y', np.float32, ('y'))
    dst_var.units = "m"
    dst_var.long_name = "y coordinate of projection"
    dst_var.standard_name = "projection_y_coordinate"
    w_nc_fid['y'][...] = np.copy(YC[:,0])
    w_nc_var = w_nc_fid.createVariable('lon', np.float32, ('y','x'))
    w_nc_var.long_name = 'longitude of 2D land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_fid.variables['lon'][...] = lon
    w_nc_var = w_nc_fid.createVariable('lat', np.float32, ('y','x'))
    w_nc_var.long_name = 'latitude of 2D land gridcell center (GCS_WGS_84), increasing from south to north'
    w_nc_var.units = "degrees_north"
    w_nc_fid.variables['lat'][...] = lat

    # create the gridIDs, lon, and lat variable
    ni_dim = w_nc_fid.createDimension('ni', grid_id_arr.size)
    nj_dim = w_nc_fid.createDimension('nj', 1)
    nv_dim = w_nc_fid.createDimension('nv', 4)

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

    w_nc_var = w_nc_fid.createVariable('xc', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_var.bounds = "xv"
    w_nc_fid.variables['xc'][...] = lon_arr
        
    w_nc_var = w_nc_fid.createVariable('yc', np.float32, ('nj','ni'))
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
    w_nc_var = w_nc_fid.createVariable('xv', np.float32, ('nv','nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell verticles (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_var = w_nc_fid.createVariable('yv', np.float32, ('nv','nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell verticles (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"

    xv_arr,yv_arr = Txy2lonlat.transform(XC_arr-side_km*500.,YC_arr+side_km*500.)
    w_nc_fid.variables['xv'][0,...] = xv_arr
    w_nc_fid.variables['yv'][0,...] = yv_arr
    xv_arr,yv_arr = Txy2lonlat.transform(XC_arr+side_km*500.,YC_arr+side_km*500.)
    w_nc_fid.variables['xv'][1,...] = xv_arr
    w_nc_fid.variables['yv'][1,...] = yv_arr
    xv_arr,yv_arr = Txy2lonlat.transform(XC_arr+side_km*500.,YC_arr-side_km*500.)
    w_nc_fid.variables['xv'][2,...] = xv_arr
    w_nc_fid.variables['yv'][2,...] = yv_arr
    xv_arr,yv_arr = Txy2lonlat.transform(XC_arr-side_km*500.,YC_arr-side_km*500.)
    w_nc_fid.variables['xv'][3,...] = xv_arr
    w_nc_fid.variables['yv'][3,...] = yv_arr

    w_nc_var = w_nc_fid.createVariable('area', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "radian^2"
    w_nc_fid.variables['area'][...] = area_arr

    w_nc_var = w_nc_fid.createVariable('area_LAEA', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells (Lambert Azimuthal Equal Area)'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "km^2"
    w_nc_fid.variables['area_LAEA'][...] = area_arr

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'mask of land gridcells (1 means land)'
    w_nc_var.units = "unitless"
    w_nc_fid.variables['mask'][...] = mask_arr

    w_nc_var = w_nc_fid.createVariable('frac', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'fraction of land gridcell that is active'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "unitless"
    w_nc_fid.variables['frac'][...] = landfrac_arr

    w_nc_var = w_nc_fid.createVariable('lambert_azimuthal_equal_area', np.short)
    w_nc_var.grid_mapping_name = "lambert_azimuthal_equal_area"
    w_nc_var.longitude_of_center = 90.
    w_nc_var.latitude_of_center = -180.
    w_nc_var.false_easting = 0.
    w_nc_var.false_northing = 0.
    w_nc_var.semi_major_axis = 6370997.
    w_nc_var.inverse_flattening = 0.

    w_nc_fid.close()  # close the new file 
    
    if outdata:
        return XC, YC, data 
    
def main():
    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """
    
    input_path= './'
    output_path = input_path
    domain_unstructured_fromraster()
    
if __name__ == '__main__':
    main()
    