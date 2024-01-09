#!/usr/bin/env python

import os, sys, time, math, glob
import re
import numpy as np

from datetime import datetime, date

#from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.dates import date2num, num2date

from pyproj import Proj
from pyproj import Transformer
from pyproj import CRS

from numpy import flipud
import gzip

import netCDF4
from builtins import str
from copy import deepcopy

#--------------------------------------------------

# Pre-process NIC IMS snowcover ASCII or GeoTiff datasets, with binary grids

def Prep_ims_grid(res, minlat=None, lonlat2xy=False):
    # INPUTS: 
    #         res    - '1km', '4km', or '24km'
    #   (optional) minlat    - min. latitude to truncate grids/data 
    #   (optional) lonlat2xy - converting grids in lon/lat to projected x/y in meters
    # OUTPUTS: 
    #         Output lon, lat

    
    # grid file names
    file_imslons = 'imslon_'+str(res)+'.bin'
    file_imslats = 'imslat_'+str(res)+'.bin'

    if res in ['1km','4km','24km']:
        imslons = np.fromfile(file_imslons,dtype=np.float32) # flat binary 4-byte floating-point
        imslats = np.fromfile(file_imslats,dtype=np.float32)
        # the grids are actually in square, with centroid in North-Polar point nearby (-80/90)
        n = np.int(np.sqrt(imslats.size))
        if n!=np.sqrt(imslats.size):
            print('NOT a squared grid-domain: x, y - ', n, np.sqrt(imslats.size))
            sys.exit()
        imslats = np.reshape(imslats, (n,n))
        imslons = np.reshape(imslons, (n,n))
        
        # truncating along latitudal, i.e. cut-off out-ring of squared-domain
        if not minlat==None:
            
            ix = np.where(np.any(imslats>=minlat,axis=0))
            ix = np.squeeze(ix)
            iy = np.where(np.any(imslats>=minlat,axis=1))
            iy = np.squeeze(iy)

            imslats = imslats[ix,:]
            imslons = imslons[ix,:]
            
            imslats = imslats[:,iy]
            imslons = imslons[:,iy]
    else:
        print('NOT supported Resolution, must be of one of 1km, 4km, or 24km')
        sys.exit(-1)

    # if needed, projection coversion from LON/LAT -> X/Y in meters
    if lonlat2xy:
        # projection used in output data layer
        # Polar stereographic ellipsoidal projection with WGS-84 ellipsoid
        #Proj4: +proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs
        #out_proj_str = "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"
        
        # for testing daymet projection
        out_proj_str = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"

        outProj = CRS.from_proj4(out_proj_str)

        # EPSG: 4326
        # Proj4: +proj=longlat +datum=WGS84 +no_defs
        inProj = CRS.from_epsg(4326)
    
        lonlat2xyOP = Transformer.from_proj(inProj, outProj, always_xy=True)
        imsxm,imsym = lonlat2xyOP.transform(imslons,imslats)
        
        # since in X/Y meters, it shall be equally intervalled as 'res', e.g. 4km
        imsxm = imsxm[0,:].astype(np.int) # grid-x nodes
        imsxm = (np.round(imsxm/10))*10 # there is ~5m error from original data
        imsym = imsym[:,0].astype(np.int) # grid-y nodes
        imsym = (np.round(imsym/10))*10 # there is ~5m error from original data
        
        
        if minlat!=None:
            return ix, iy, imsxm, imsym, imslons, imslats
        else:
            return imsxm, imsym, imslons, imslats  # note: imsxm/imsym in 1D-axis, imslons/imslats in 2-D for each cell
    
    else:
        
        if minlat!=None:
            return ix, iy, imslons, imslats
        else:
            return imslons, imslats  # this is in 2-D, and fully paired  lon/lat for each cells, But NOT evenly longitude-interval along latitudal-axis

#------------------------------------------------------------
# processing data file *.asc or *.tiff each a day
def Prep_ims_snowcov(ifile, fileheader, varname, alldata={}):
    # 
    # file named as 'imsYYYYDOY_{res}_v{version}.asc.gz/tif.gz'
    if os.path.isfile(ifile):
        fpath,fname = os.path.split(ifile)
        # date
        yrdoy = fname.split('_')[0]
        yrdoy = yrdoy.replace('ims','')

        year  = np.int(yrdoy[:4])
        doy   = np.int(yrdoy[4:])
        doy0  = date(year,1,1)
        daynums = date2num(doy0)+doy-1
        mid_day = num2date(daynums).date()

        try:
            data = np.genfromtxt(ifile,dtype=np.int8,delimiter=1,skip_header=30)
        except:
            
            # unpacked format for '24km' resolution dataset (some files)
            if '24km' in ifile:
                
                with gzip.open(ifile, 'rt') as f:
                    for line in f:
                        # remove header txts
                        meta = line.split()
                        if len(meta)>0 and meta[0].isnumeric(): 
                            data = line.strip()+f.read() # line already read, must be included
                            break # out of for loop
                if len(data)>0:
                    data = re.split('\W+|\bt',data.strip())
                    data = np.asarray(data,dtype=np.float32)
                    n = np.int(np.sqrt(len(data)))
                    if n==np.sqrt(len(data)):
                        data = np.reshape(data,(n,n))
                    else:
                        print('Error: Data size is not a squared format')
                        return -1
                else:
                    print ('Error: Data size is zero')
                    return -1
                
            
            else:
                print ('Error: data reading from file')
                return -1
            
        # appears vardatas are in S-N ordered, so must be flip upside down to match with grids
        data = np.flipud(data)
        

        # data re-grouping as follows, so that it can be read in Visit using blue-scale :
        # 0 outside the coverage area  ==> none (-0.5)
        # 1 sea                        ==> water-body (-2)
        # 2 land without snow          ==> snow-free land (0)
        # 3 Sea ice                    ==> water-body (-1), 164 for unpacked format data
        # 4 land covered with snow     ==> snow-covered land (1), 165 for unpacked format data
        data[data == 0] = -0.5 # when plot in Visit, negative could be blocked as continuing of non-snow
        data[data == 1] = -2
        data[(data == 3) | (data == 164)] = -1
        data[data == 2] = 0
        data[(data == 4) | (data == 165)] = 1
        
        if len(alldata)<=0:
            
            alldata["time"] = np.empty(1,dtype=type(daynums))
            alldata["time"][0] = deepcopy(daynums)
            
            alldata[varname] = np.empty((1,)+data.shape)
            temp = np.reshape(data, (1,)+data.shape)
            alldata[varname][0,] = deepcopy(temp)
        
        else: # append data into old one, when input
            alldata["time"] = np.vstack((alldata["time"], daynums))
            temp = np.reshape(data, (1,)+data.shape)
            alldata[varname] =  np.vstack((alldata[varname], temp))
    
    del data
    
    return mid_day, year, doy, alldata

#--------------------------------------------------------------------------------

def IMS_ELM_gridmatching(Grid1_Xdim, Grid1_Ydim, Grid2_x, Grid2_y, Grid2_dxy=None, \
                         Grid1prjstr='', Grid2prjstr='', Grid1_cells=()):
    
    # Match Grid2 within each Grid1, so return Grid2-xy index for each Grid1 cell
    # by default, (1) Grid1 (parent) in xlon/ylat mesh, Grid2 in lonx/laty 1-D grids
    #             (2) Grid1 XY is grid-mesh-nodes, Grid2xy is grid-centroids. Then is good for searching Grid2 in Grid1
    #             (3) all cells in Grid1 are assigned Grid2's cell-index - maybe only those indiced in 'Grid1_cells'
    
    # it's supposed that: Grid1 X/Y are in 1-D regularly-intervaled (may not evenly) nodes along axis
    # while, Grid2 might be either like Grid2 or in 2-D mesh, with optional interval (Grid2_dxy in length of 2, e.g. [1000 1000] )
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
    # Polar stereographic ellipsoidal projection with WGS-84 ellipsoid
    #Proj4: +proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs
    # geoxy_proj_str = "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"
    if (Grid2prjstr==''):
        grid2Proj = CRS.from_epsg(4326) # in lon/lat coordinates
    else:
        grid2Proj = CRS.from_proj4(Grid2prjstr)

    if (Grid1prjstr==''):
        # EPSG: 4326
        # Proj4: +proj=longlat +datum=WGS84 +no_defs
        grid1Proj = CRS.from_epsg(4326) # in lon/lat coordinates
    else:
        # e.g. Grid1 is in geox/y as well, e.g. DAYMET geox/y here
        #lambert_conformal_conic:grid_mapping_name = "lambert_conformal_conic" ;
        #lambert_conformal_conic:longitude_of_central_meridian = -100. ;
        #lambert_conformal_conic:latitude_of_projection_origin = 42.5 ;
        #lambert_conformal_conic:false_easting = 0. ;
        #lambert_conformal_conic:false_northing = 0. ;
        #lambert_conformal_conic:standard_parallel = 25., 60. ;
        #lambert_conformal_conic:semi_major_axis = 6378137. ;
        #lambert_conformal_conic:inverse_flattening = 298.257223563 ;
        #Proj4: +proj=lcc +lon_0=-100 +lat_1=25 +lat_2=60 +k=1 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m  +no_defs
        #grid1prjstr = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
        grid1Proj = CRS.from_proj4(Grid1prjstr)
    
    # only if 2 grids are in different projections, do tansformation
    if (Grid2prjstr.strip() != Grid1prjstr.strip()):
        T2to1 = Transformer.from_proj(grid2Proj, grid1Proj, always_xy=True)
        Grid2_gxx,Grid2_gyy = T2to1.transform(Grid2_xx,Grid2_yy)
        
        if(Grid1prjstr==''):  # because grid2 will transformed as grid1's proj
            ij=np.where(Grid2_gxx<0.0)
            if(len(ij[0])>0): Grid2_gxx[ij]=Grid2_gxx[ij]+360.0 # for convenience, longitude from 0~360
            ij=np.where(Grid1_x<0.0)
            if(len(ij[0])>0): Grid1_x[ij]=Grid1_x[ij]+360.0 # for convenience, longitude from 0~360
    
    else:
        Grid2_gxx = Grid2_xx
        Grid2_gyy = Grid2_yy

    # IMS grids' index (Grid2) included in each ELM land-grid (Grid1)
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
        
        # none of IMS cell centroid inside a ELM grid (e.g. IMS cell is larger and centroid beyond any ELM grids), 
        # find the close one instead
        if(len(ij[0])<1):
            # only search within Grid2' range
            if (Grid2_dxy!=None):
                print()
                
            dx=abs(Grid2_gxx-(iwst+iest)/2.0)
            dy=abs(Grid2_gyy-(jnth+jsth)/2.0)
            dxy = dx*dx+dy*dy
            ij = np.unravel_index(np.argmin(dxy, axis=None), dxy.shape)
            
        
        Grid2in1_indx[str(indx)] = deepcopy(ij)

    # done with all grids
    return Grid2in1_indx, Grid2_gxx, Grid2_gyy
   
#-------------------Writing submodule-----------------------------------------------
# Write to geo-referenced CF compliant nc file, if filename given

def Write1GeoNc(vars, vardatas, ptxy=[], ncfname='', newnc=True, FillValue=None, geoprj=''):
    # INPUTS: vars      - variable names, separated by ','. if 'all' means every var-key in 'vardatas'
    #         vardatas  - python np list, with dicts/data
    #    (optional) ptxy    - paired lon/lat (x/y), in [x/lon,y/lat] (only 1 point)
    #    (optional) ncfname - if not empty, write processed data into (geo)NC file
    #    (optional) newnc   - if not True, write into existed 'ncfname' nc file
    #    (optional) FillValue - 
    #    (optional) geoprj   - if not empty, GeoX/Y projection names
    # OUTPUTS (optional): if NOT write to NC file and 'vars' only has 1 variable
    #                     Output year, doy, lon, lat, data for 'vars' in np.array  
    
    if ncfname !='':
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
                    vlat.long_name = 'Northing ('+geoprj+')'
                
                    vlon = ncfile.createVariable('geox', np.float32, ('geox',))
                    vlon.units = 'meters'
                    vlon.long_name = 'Easting ('+geoprj+')'

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
                vdaysnum.long_name = 'days since 0000-01-01 UTC + 1'

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
                if '1km' in ncfname:
                    ncfile.description = 'Northern Hemisphere daily snow coverage @1km resolution '
                elif '24km' in ncfname:
                    ncfile.description = 'Northern Hemisphere daily snow coverage @24km resolution '
                elif '4km' in ncfname:
                    ncfile.description = 'Northern Hemisphere daily snow coverage @4km resolution '
                ncfile.data_source = ('US National Ice Center (USNIC), ' +
                    'Interactive Multsensor Snow and Ice Mapping Service (IMS),' +
                    'Version 1.2 ')
                if ('ELM' in ncfname):
                    ncfile.data_source2 = ('E3SM v1.1 Land Model offline simulations for No60 and above, ' +
                                           'For NGEE-Arctic Project sponsored by DOE Office of Science')

                ncfile.data_citation = ('U.S. National Ice Center. 2008. Updated daily. ' +
                    'IMS Daily Northern Hemisphere Snow and Ice Analysis at 1 km, 4 km, and 24 km Resolutions, Version 1.2 [4km data].' +
                    'Boulder, Colorado USA. '+
                    'NSIDC: National Snow and Ice Data Center. '+
                    'doi: https://doi.org/10.7265/N52R3PMC. '+
                    '2020-03-24.')
                if ('ELM' in ncfname):
                    ncfile.history = ('2020-04-14: IMS snow data conversion from ASCII format. ' +
                                     '2020-05-15: IMS data either aggrated into ELM land-cell, OR, ' +
                                     ' ELM data assinged for each IMS grid-cell, upon which griding-system adopted')
                else:
                    ncfile.history = '2020-04-14: conversion from ASCII format.'
                    
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
                                            dimensios=('time','lat','lon'), \
                                            zlib=True, fill_value=FillValue) # note: unlimited dimension is leftmost
                
                vtemp.units = ''
                vtemp.standard_name = varname.strip() # this is a CF standard name
                if '24km' in ncfname:
                    vtemp.long_name =  varname.strip()+' at 24km resolution'
                elif '4km' in ncfname:
                    vtemp.long_name =  varname.strip()+' at 4km resolution'
                elif '1km' in ncfname:
                    vtemp.long_name =  varname.strip()+' at 1km resolution'
                if ('ELM' in ncfname):
                    if ('Day' in varname or 'Doy' in varname): 
                        vtemp.Key = "0~365 = days or doy for land with snow, -22 = sea ice (if any), -55/nan = sea, -11 = beyond map coverage" ;
                    else:
                        vtemp.Key = "0 = land without snow, 0~1.0 = land with snow fraction, -2/nan = sea, -1 = beyond map coverage" ;
                else:
                    if ('Day' in varname or 'Doy' in varname): 
                        vtemp.Key = "0~365 = days or doy for land with snow, -22 = sea ice (if any), -55/nan = sea, -11 = beyond map coverage" ;
                    else:
                        vtemp.Key = "0 = land without snow, 0~1 = land with snow, -1 = sea with ice, -2/nan = sea, -0.50 = beyond map coverage" ;
                    
                vtemp[0:nt,:,:] = data.astype(np.float32)
            else:
                vtemp = ncfile.variables[varname]
                vtemp[prv_nt:prv_nt+nt,:,:] = data.astype(np.float32)
    # done with 'for varname in vars'
    
    # end of for varname in vars:
    if ncfname !='': ncfile.close()
    
    

#-------------------PLotting submodule-----------------------------------------------
# plotting 1 graph with 1 map

def SingleMap(varname, vardatas, figno=None):
    '''
    fig = plt.figure(figsize=(11.5, 8.5))

    # grid box
    #['WESTBOUNDINGCOORDINATE','EASTBOUNDINGCOORDINATE',
    #'NORTHBOUNDINGCOORDINATE','SOUTHBOUNDINGCOORDINATE']

    lon_w = np.float(vardatas['WESTBOUNDINGCOORDINATE'])
    lon_e = np.float(vardatas['EASTBOUNDINGCOORDINATE'])
    lat_n = np.float(vardatas['NORTHBOUNDINGCOORDINATE'])
    lat_s = np.float(vardatas['SOUTHBOUNDINGCOORDINATE'])
    #if(lon_w<0.): lon_w = lon_w + 360.
    #if(lon_e<0.): lon_e = lon_e + 360.
    
    nlat, nlon = vardatas[varname].shape
    dlon = (lon_w - lon_e) / nlon
    dlat = (lat_n - lat_s) / nlat
        
    # Construct the grid in lat/lon.
    x = np.linspace(lon_e, lon_e + dlon*nlon, nlon)
    y = np.linspace(lat_s, lat_s + dlat*nlat, nlat)
    lon, lat = np.meshgrid(x, y)
    
    # appears vardatas are in S-N/E-W ordered, so shall be flip over
    data = vardatas[varname]
    #data = np.fliplr(data) # not with nic-ims data ??? - be super-cautious here
    #data = np.flipud(data) # not with nic-ims data ??? - be super-cautious here
    
    # mid of day
    mid_day = vardatas['RANGEBEGINNINGDATE']+ (vardatas['RANGEENDINGDATE']-vardatas['RANGEENDINGDATE'])/2.0

    
    m = Basemap(projection='npstere', resolution='l',
                boundinglat=30,
                lon_0=-90)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(30., 90., 10.), labels=[False,False,False,False])
    m.drawmeridians(np.arange(-180, 180., 20.), labels=[True,True,True,True])
    
    # Bin the data as follows:
    # 0% snow
    # 1-99% snow
    # 100% snow
    # lake ice (107)
    # night (111)                ==> missing (255)
    # inland water (237)         ==> water-body (255)
    # ocean (239)                ==> water-body (255)
    # cloud-obscured water (250) ==> water-body (255)
    # data not mapped (253)      ==> missing (255)
    # fill (255)                 ==> missing (255)
    data[abs(data-250)<0.1] = 254 # data type is 'uint8', cannot be over 255
    data[abs(data-237)<0.1] = 254
    data[abs(data-239)<0.1] = 254
    
    data[abs(data-255)<0.1] = 255
    data[abs(data-253)<0.1] = 255
    data[abs(data-111)<0.1] = 255

    lst = ['#00ff00', # 0% snow
           '#888888', # 1-99% snow
           '#ffffff', # 100% snow
           '#63c6ff', # lake ice
           '#0000cc', # water-body
           '#00ffcc'] # missing
    cmap = mpl.colors.ListedColormap(lst)
    bounds = [0, 1, 100, 107, 254, 255, 256]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # since use 'polar' projection, we need to truncate edge, otherwise 'pcolrmesh' (below) complaining non-finite x/y
    lat=lat[1:,:] # remove the -90 latitude row
    lon=lon[1:,:]
    data=data[1:,:]
    
    
    # Render the image in the projected coordinate system.
    m.pcolormesh(lon, lat, data,
                 latlon=True, cmap=cmap, norm=norm)

    basename = 'Daily CMG Snow Cover'
    long_name = mid_day
    plt.title('{0}\n{1}'.format(basename, long_name))
    fig = plt.gcf()

    color_bar = plt.colorbar(orientation='horizontal')
    color_bar.set_ticks([0.5, 50.5, 103.5, 180.5, 254.5, 255.5])
    color_bar.set_ticklabels(['0%\nsnow', '1-99%\nsnow', '100%\nsnow', 'lake\nice',
                              'water\nbody',
                              'data\nnot available'])
    color_bar.draw_all()

    #if saving the figure
    if (not figno is None): 
        ofname = 'Figure_'+figno+'.pdf'
        plt.savefig(ofname)
    
    plt.show()
    plt.close('all')
    
    '''


#####################################################################
# module testing

# mapping
#SingleMap('Day_CMG_Snow_Cover', vardatas)

# write to geo-referenced NC
#Write1GeoNc(ncfname, vardatas)

# read NIC-IMS snow cover data
#dataformat = 'ascii'
#year = '2004'
#res = '4km'
#imslons, imslats = \
#    Prep_ims_grid(res, lonlat2xy=False)

#mid_day, year, doy, data = \
#    Prep_ims_snowcov(filename)

#print('done')
