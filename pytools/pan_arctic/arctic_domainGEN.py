
# utility to generate 1D and 2D domain
# use previous dataPartition module. Code need to be clean up

import os
import math
import netCDF4 as nc
import numpy as np
from pyproj import Transformer
from pyproj import CRS
from pytools.commons_utils.Modules_netcdf import geotiff2nc

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
    nx_nodata = sum(nondata,1)
    for i in range(1,len(allyc)-1):
        if nx_nodata[i]<len(allxc) or nx_nodata[i-1]<len(allxc): break 
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
    w_nc_var = w_nc_fid.createVariable('grid_code', np.float32, ('y','x'))
    w_nc_var.long_name = 'code of 2D land gridcells'
    w_nc_var.units = "-"
    w_nc_fid.variables['grid_code'][...] = data

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
    w_nc_var.crs_wkt = crs_wkt
    w_nc_var.longitude_of_center = 90.
    w_nc_var.latitude_of_center = -180.
    w_nc_var.false_easting = 0.
    w_nc_var.false_northing = 0.
    w_nc_var.semi_major_axis = 6370997.
    w_nc_var.inverse_flattening = 0.

    w_nc_fid.close()  # close the new file 
    
    if outdata:
        return XC, YC, data 

def domain_subsetbymask(input_pathfile='./share/domains/domain.clm/domain.lnd.360x720_cruncep.c20190221.nc', mask1dncf='./share/domains/domain.clm/domain.lnd.pan-arctic_CAVM.1km.1d.c241018.nc', maskncv='mask',WRITE2D=True,outdata=False):

    import geopandas as geopd
    from shapely.geometry import Point
    from geopandas.tools import sjoin
    from shapely.geometry import MultiLineString, LineString
    from shapely.ops import polygonize

    
    # targetting mask file in 1D
    mask_f = nc.Dataset(mask1dncf,'r')
    mask_v = mask_f[maskncv][...]
    mask_checked = np.where(mask_v==1)
    mask_xc= mask_f['xc'][...][mask_checked]
    mask_yc= mask_f['yc'][...][mask_checked]
    mask_xc[mask_xc<0.0]=360+mask_xc[mask_xc<0.0] # better for interpolating or search if 0~360
    
    # 
    # Open the source file
    src = nc.Dataset(input_pathfile, 'r')
    # points = [y,x] coordinates for src's grid
    src_yc = src.variables['yc'][...]
    src_xc = src.variables['xc'][...]
    src_xc[src_xc<0.0]=360+src_xc[src_xc<0.0]
    src_xv = src.variables['xv'][...]
    src_xv[src_xv<0.0]=360+src_xv[src_xv<0.0]
    src_yv = src.variables['yv'][...]
    #nv dim position
    nvdim = src.variables['xv'].dimensions.index('nv')

    src_mask = src.variables['mask'][...]
    src_area = src.variables['area'][...]
    src_landfrac = src.variables['frac'][...]

    
    # for pan-arctic: longitutde is all, and latitude is above min. lat in masked_yc
    idx = np.where(src_yc>=min(mask_yc)-0.5)  # 0.5 allows a wider range of latitude
    landmask = src_mask[idx]
    landfrac = src_landfrac[idx]
    area = src_area[idx]
    xc = src_xc[idx]
    yc = src_yc[idx]
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
    yv1=yv1[idx]; xv1=xv1[idx]
    yv2=yv2[idx]; xv2=xv2[idx]
    yv3=yv3[idx]; xv3=xv3[idx]
    yv4=yv4[idx]; xv4=xv4[idx]
    #yv=np.stack((yv1,yv2,yv3,yv4), axis=0)
    #xv=np.stack((xv1,xv2,xv3,xv4), axis=0)
    
    # re-do landmask from targetted area
    points = geopd.GeoDataFrame({"x":mask_xc,"y":mask_yc})
    points['geometry'] = points.apply(lambda p: Point(p.x, p.y), axis=1)
 
    # xv,yv to polygons
    vpts1=[(x,y) for x,y in zip(xv1,yv1)]
    vpts2=[(x,y) for x,y in zip(xv2,yv2)]
    vpts3=[(x,y) for x,y in zip(xv3,yv3)]
    vpts4=[(x,y) for x,y in zip(xv4,yv4)]
    lines = [(vpts1[i], vpts2[i], vpts3[i], vpts4[i], vpts1[i]) for i in range(len(vpts1))]
    lines_str = [LineString(lines[i]) for i in range(len(lines))]
    polygons = list(polygonize(MultiLineString(lines_str)))
    # if regular x/y mesh may be do like following: 
    #x = np.linspace(src_xv[0,0], xv[0,-1], len(src_xc[1]))         
    #y = np.linspace(src_yv[0,0], yv[0,-1], len(src_yc[0]))    
    #hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(x[:-1], x[1:]) for yi in y]
    #vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(y[:-1], y[1:]) for xi in x]
    #lines_str = vlines+hlines
    #polygons = list(polygonize(MultiLineString(lines_str)))

    # ELM land domain in dims of [nj,ni] or [lat, lon]
    # new bound-box and grid indices
    if len(src_mask.shape)>1:
        X_axis = src_xc[0,:]
        X_axis = X_axis[np.where((X_axis>=np.min(xc)) & (X_axis<=np.max(xc)))]
        Y_axis = src_yc[:,0]
        Y_axis = Y_axis[np.where((Y_axis>=np.min(yc)) & (Y_axis<=np.max(yc)))]
    else:
        X_axis = src_xc[0,:]
        X_axis = X_axis[np.where((X_axis>=np.min(xc)) & (X_axis<=np.max(xc)))]
        Y_axis = np.median(yc[:,0])  # unstructured grids may have varied lat/yc, so this is the mid-point
    YY, XX = np.meshgrid(Y_axis,X_axis, indexing='ij')
    xid = np.indices(XX.shape)[1]
    if len(XX.shape)>1: 
        yid = np.indices(XX.shape)[0]
    else:
        yid = [0]
    xyid = np.indices(XX.flatten().shape)[0]
    xyid = np.reshape(xyid, xid.shape)
    
    # new indices of polygons in newly-created bounding-box
    polygonized = np.isin(XX, xc) & np.isin(YY, yc)
    ji = np.nonzero(polygonized) 
    grid_xid = xid[ji]     
    grid_yid = yid[ji]    
    grid_id = xyid[ji]  
    
    grids = geopd.GeoDataFrame({"gid":grid_id,"xid":grid_xid,"yid":grid_yid,"geometry":polygons})
 
    pointsInPolys = sjoin(points, grids[['gid','xid','yid','geometry']], how='inner')
    points_gid = pointsInPolys.gid
    pts_gid, pts_counts  = np.unique(points_gid, return_counts=True)
    
    landmask[...] = 0
    landmask[pts_gid]= 1

    # re-do frac of landed mask
    kmratio_lon2lat = np.cos(np.radians(yc))
    re_km = 6370.997
    xside_km = (abs(xv1-xv2)+abs(xv3-xv4))/2.0*(math.pi*re_km/180.0)   # 
    yside_km = (abs(yv1-yv3)+abs(yv2-yv4))/2.0*(math.pi*re_km/180.0*kmratio_lon2lat)
    area_km2 = xside_km*yside_km     # grid area in km2
    landfrac[...] = 0.0
    landfrac[pts_gid] = pts_counts/area_km2[pts_gid]
    landfrac[np.where(landfrac>1.0)]=1.0


    # output arrays and exit (i.e. NO nc file writing)
    if outdata:
        subdomain = {}
        subdomain['x'] = X_axis
        subdomain['y'] = Y_axis
        subdomain['gridID']  = grid_id
        subdomain['gridXID'] = grid_xid
        subdomain['gridYID'] = grid_yid
        subdomain['xc'] = xc
        subdomain['yc'] = yc
        subdomain['xv'] = np.stack((xv1,xv2,xv3,xv4), axis=0)
        subdomain['yv'] = np.stack((yv1,yv2,yv3,yv4), axis=0)
        subdomain['area'] = area
        subdomain['mask'] = landmask
        subdomain['frac'] = landfrac
        
        return idx, subdomain
        # this will stop here



    # to 2D but masked    
    if WRITE2D:
        grid_id_arr = np.reshape(grid_id,xyid.shape)
        grid_xid_arr = np.reshape(grid_xid,xyid.shape)
        grid_yid_arr = np.reshape(grid_yid,xyid.shape)
        
        lon_arr = np.reshape(xc,xyid.shape)
        lat_arr = np.reshape(yc,xyid.shape)
        xv_arr = np.stack((np.reshape(xv1,xyid.shape), \
                           np.reshape(xv2,xyid.shape), \
                           np.reshape(xv3,xyid.shape), \
                           np.reshape(xv4,xyid.shape)), axis=0)
        yv_arr = np.stack((np.reshape(yv1,xyid.shape), \
                           np.reshape(yv2,xyid.shape), \
                           np.reshape(yv3,xyid.shape), \
                           np.reshape(yv4,xyid.shape)), axis=0)
        area_arr = np.reshape(area,xyid.shape)
        mask_arr = np.reshape(landmask,xyid.shape)
        landfrac_arr = np.reshape(landfrac,xyid.shape)
    else:
    # to 1D only masked
        masked = np.where((landmask==1))
        grid_id_arr = grid_id[masked]
        grid_xid_arr = grid_xid[masked]
        grid_yid_arr = grid_yid[masked]
        lon_arr = xc[masked]
        lat_arr = yc[masked]
        xv_arr = np.stack((xv1[masked],xv2[masked],xv3[masked],xv4[masked]), axis=0)
        yv_arr = np.stack((yv1[masked],yv2[masked],yv3[masked],yv4[masked]), axis=0)
        area_arr = area[masked]
        mask_arr = landmask[masked]
        landfrac_arr = landfrac[masked]    
 
    # write to nc file
    file_name = input_pathfile+'_remasked'
    print("The domain file is " + file_name)

    # Open a new NetCDF file to write the domain information. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF4')
    if WRITE2D:
        w_nc_fid.title = '2D domain file for pan arctic region, based on CAVM vegetation map'
    else:
        w_nc_fid.title = '1D domain file for pan arctic region, based on CAVM vegetation map'

    # Create new dimensions of new coordinate system

    # create the gridIDs, lon, and lat variable
    xdim = w_nc_fid.createDimension('x', len(X_axis))
    ydim = w_nc_fid.createDimension('y', len(Y_axis))
    w_nc_var = w_nc_fid.createVariable('x', np.float64, ('x'))
    w_nc_var.long_name = 'longitude of x-axis'
    w_nc_var.units = "degrees_east"
    w_nc_fid.variables['x'][...] = X_axis

    w_nc_var = w_nc_fid.createVariable('y', np.float64, ('y'))
    w_nc_var.long_name = 'latitude of y-axis'
    w_nc_var.units = "degrees_north"
    w_nc_fid.variables['y'][...] = Y_axis

    w_nc_var = w_nc_fid.createVariable('lon', np.float64, ('y','x'))
    w_nc_var.long_name = 'longitude of 2D land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_fid.variables['lon'][...] = XX

    w_nc_var = w_nc_fid.createVariable('lat', np.float64, ('y','x'))
    w_nc_var.long_name = 'latitude of 2D land gridcell center (GCS_WGS_84), increasing from south to north'
    w_nc_var.units = "degrees_north"
    w_nc_fid.variables['lat'][...] = YY

    if WRITE2D:
        ni_dim = w_nc_fid.createDimension('ni', grid_id_arr.shape[1])
        nj_dim = w_nc_fid.createDimension('nj', grid_id_arr.shape[0])
        nv_dim = w_nc_fid.createDimension('nv', 4)
    else:
        ni_dim = w_nc_fid.createDimension('ni', grid_id_arr.size)
        nj_dim = w_nc_fid.createDimension('nj', 1)
        nv_dim = w_nc_fid.createDimension('nv', 4)

    # for 2D <--> 1D indices
    w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId in the pan-Arctic domain'
    w_nc_var.decription = "start from #0 at the lower left corner of the domain, covering all land and ocean gridcells" 
    w_nc_fid.variables['gridID'][...] = grid_id_arr

    w_nc_var = w_nc_fid.createVariable('gridXID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId x in the pan-Arctic domain'
    w_nc_var.decription = "start from #0 at the lower left corner and from west to east of the domain, with gridID=gridXID+gridYID*len(y)" 
    w_nc_fid.variables['gridXID'][...] = grid_xid_arr

    w_nc_var = w_nc_fid.createVariable('gridYID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId y in the pan-Arctic domain'
    w_nc_var.decription = "start from #0 at the lower left corner and from south to north of the domain, with gridID=gridXID+gridYID*len(y)" 
    w_nc_fid.variables['gridYID'][...] = grid_yid_arr

    #
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
    #
    w_nc_var = w_nc_fid.createVariable('xv', np.float32, ('nv','nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell verticles (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_var = w_nc_fid.createVariable('yv', np.float32, ('nv','nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell verticles (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"

    w_nc_fid.variables['xv']= xv_arr
    w_nc_fid.variables['yv']= yv_arr

    w_nc_var = w_nc_fid.createVariable('area', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "radian^2"
    w_nc_fid.variables['area'][...] = area_arr

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'mask of land gridcells (1 means land)'
    w_nc_var.units = "unitless"
    w_nc_fid.variables['mask'][...] = mask_arr

    w_nc_var = w_nc_fid.createVariable('frac', np.float32, ('nj','ni'))
    w_nc_var.long_name = 'fraction of land gridcell that is active'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "unitless"
    w_nc_fid.variables['frac'][...] = landfrac_arr

    w_nc_fid.close()  # close the new file 
    
def ncdata_subsetbyindx(indx, mask=[], input_pathfile='./lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c211019.nc', indx_dim=['lsmlat','lsmlon'], indx_dim_new='gridcell'):
    
    # check if empty indices: indx is a tuple of np.where outputs
    if len(indx[0])>0:
        #
        if not len(indx)==len(indx_dim):
            print('inconsistent dimensions: ', len(indx), len(indx_dim))
            os.sys.exit(-1)
        
        #
        ncfilein  = input_pathfile
        ncfileout = input_pathfile.split('/')[-1]
        ncfileout = ncfileout.split('.nc')[0]+'_subset.nc'

        # masked only
        indx_masked = {}
        if len(mask)>0:            
            for i in range(len(indx)):
                indx_masked[i]=indx[i][mask>0]
        else:
            indx_masked = indx
         
        #subset and write to ncfileout
        with nc.Dataset(ncfilein,'r') as src, nc.Dataset(ncfileout, "w") as dst:
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # dimensions
            for dname, dimension in src.dimensions.items():
                len_dimension = len(dimension)
                if dname in indx_dim:
                    dim_len = len(indx_masked[indx_dim.index(dname)])
                    if (dim_len>0):
                        len_dimension = dim_len
                    elif (dim_len==0):
                        dimension=None
                    
                    if indx_dim_new not in dst.dimensions.keys():
                        #only need once, i.e. indx_dim will be flatten into index_dim_new
                        dst.createDimension(indx_dim_new, len_dimension if not dimension.isunlimited() else None)
                    
                else:
                    dst.createDimension(dname, len_dimension if not dimension.isunlimited() else None)
            #   
        
            # all variables data
            for vname, variable in src.variables.items():

                if all(d in variable.dimensions for d in indx_dim):
                    vdim = np.array(variable.dimensions)
                    vdim = [d.replace(indx_dim[0], indx_dim_new) for d in vdim]
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
                varvals = src[vname][...]

                # dimensions to extract data
                if all(d in variable.dimensions for d in indx_dim):
                    vd = variable.dimensions
                    i = vd.index(indx_dim[-1]) # position of last dim in 'indx_dim'
                    if i==len(vd)-1:
                        if len(indx_masked)==1: varvals = varvals[...,indx_masked[0]]
                        if len(indx_masked)==2: varvals = varvals[...,indx_masked[0],indx_masked[1]]
                        if len(indx_masked)==3: varvals = varvals[...,indx_masked[0],indx_masked[1],indx_masked[2]]
                    elif i==len(vd)-2:
                        if len(indx_masked)==1: varvals = varvals[...,indx_masked[0],:]
                        if len(indx_masked)==2: varvals = varvals[...,indx_masked[0],indx_masked[1],:]
                        if len(indx_masked)==3: varvals = varvals[...,indx_masked[0],indx_masked[1],:]
                    elif i==len(vd)-3:
                        if len(indx_masked)==1: varvals = varvals[...,indx_masked[0],:,:]
                        if len(indx_masked)==2: varvals = varvals[...,indx_masked[0],indx_masked[1],:,:]
                        if len(indx_masked)==3: varvals = varvals[...,indx_masked[0],indx_masked[1],:,:]
                    elif i==len(vd)-4:
                        if len(indx_masked)==1: varvals = varvals[...,indx_masked[0],:,:,:]
                        if len(indx_masked)==2: varvals = varvals[...,indx_masked[0],indx_masked[1],:,:,:]
                        if len(indx_masked)==3: varvals = varvals[...,indx_masked[0],indx_masked[1],indx_masked[2],:,:,:]
                    else:
                        print('unsupported dimension indices: ', i, indx_dim, vd)
                        os.sys.exit(-1)
                #                            
                dst[vname][...] = varvals
                    
            # variable-loop        
                
        
        print('subsetting done successfully!')
        
    else:
        print('NO subsetting done due to indices NOT provided!')
    
  
def main():
    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """
    
    input_path= './'
    output_path = input_path
    
    #domain_unstructured_fromraster()
    
    idx_original, newdomain = domain_subsetbymask(outdata=True)
    #domain_subsetbymask(input_pathfile='domain.lnd.Daymet_NA.1km.2d.c240327.nc')
    
    ncdata_subsetbyindx(idx_original, mask=newdomain['mask'])
    
if __name__ == '__main__':
    main()
    