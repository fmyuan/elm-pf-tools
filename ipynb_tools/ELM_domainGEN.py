
# utility to generate 1D and 2D domain
# use previous dataPartition module. Code need to be clean up

import math
import netCDF4 as nc
import numpy as np
from pyproj import Transformer
from pyproj import CRS


#----------------------------------------------------------------------------             
# convert geotiff to ncfile
def geotiff2nc(file, bandinfos=None, outdata=False):
    # file: file name
    # bandinfos in 2-D strings: tiff file' bands and its NC name, units, long_name, standard_name
    import rasterio

    f = rasterio.open(file, mode='r')
    alldata = f.read()
    nband,ny,nx = np.shape(alldata)
    try:
        filled_value = f.nodata
        alldata = np.ma.masked_equal(alldata,filled_value)
    except:
        filled_value = -9999
        
    
    # geox/y coordinates in 1-D (axis), centroid
    geox = np.arange(nx)*f.transform[0]+f.transform[2]+f.transform[0]/2.0
    geoy = np.arange(ny)*f.transform[4]+f.transform[5]+f.transform[4]/2.0
        
    
    if not bandinfos is None and 'bands' in bandinfos.keys():
        nvars = len(bandinfos['bands'])

        # create NetCDF output file
        ncof = nc.Dataset(file+'.nc','w',clobber=True)
       
        # create dimensions, variables and attributes:
        ncof.createDimension('geox',nx)
        ncof.createDimension('geoy',ny)
    
        geoxo = ncof.createVariable('geox',np.double,('geox'))
        if f.crs.is_projected:
            geoxo.units = 'm'
            geoxo.long_name = 'x coordinate of projection'
            geoxo.standard_name = 'projection_x_coordinate'
        else:
            geoxo.units = 'degrees_east'
            geoxo.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
            geoxo.standard_name = "longitude"
    
        geoyo = ncof.createVariable('geoy',np.double,('geoy'))
        if f.crs.is_projected:
            geoyo.units = 'm'
            geoyo.long_name = 'y coordinate of projection'
            geoyo.standard_name = 'projection_y_coordinate'
        else:
            geoyo.units = 'degrees_north'
            geoyo.long_name = 'latitude of land gridcell center (GCS_WGS_84), decreasing from north to south'
            geoyo.standard_name = 'latitude'
    
        georef = ncof.createVariable('crs_wkt','int')
        georef.CRS = f.crs.wkt
    
        bounds = ncof.createVariable('bounds','int')
        bounds.box=f.bounds
    
        res = ncof.createVariable('resolution','int')
        res.box=f.res

        #write x,y
        geoxo[:]=geox
        geoyo[:]=geoy
    
        # create variables        
        varo = {}
        for iv in range(0,nvars):
            v = bandinfos['bands'][iv]
            varo[v] = ncof.createVariable(v, np.double,  ('geoy', 'geox'), fill_value=filled_value)
            if 'units' in bandinfos.keys():
                varo[v].units = bandinfos['units'][iv]
            if 'long_name' in bandinfos.keys():
                varo[v].long_name = bandinfos['long_name'][iv]
            if 'standard_name' in bandinfos.keys():
                varo[v].standard_name = bandinfos['standard_name'][iv]
    
            # write variable
            varo[v][:,:] = alldata[iv]
        #        
        ncof.close()
    
    if outdata:
        return geox,geoy,f.res,f.crs,alldata
#

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

    #          
#    
#--------------------------------------------------------------------------------------------------------
  
    
if __name__ == '__main__':
    # create an unstructured domain from a raster image, e.g. CAVM image of land cover type, including veg
    
    #domain_unstructured_fromraster_cavm(output_pathfile='./domain.lnd.pan-arctic_CAVM.0.01deg.1D.c250623.nc', \
    #                               rasterfile='./raster_cavm_v1_01d.tif')
    
    domain_unstructured_fromraster_cavm(output_pathfile='./domain.lnd.original.1D.c250624_TFSarcticpfts.nc', \
                                   rasterfile='./deciduous_shrub_toolik.tif')
    
    # create an unstructured domain from grids of daymet tile
    #domain_unstructured_fromdaymet(tileinfo_ncfile='./daymet_tiles.nc')

    return # for only do domain.nc writing
    
