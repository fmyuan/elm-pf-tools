#!/usr/bin/env python

import os, sys, time, math
import re
import numpy as np

import h5py as h5
from datetime import datetime, date

from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser
from numpy import flipud

import netCDF4
from matplotlib.dates import date2num
from builtins import str
from copy import deepcopy

#--------------------------------------------------

#-------------------EOS hdf4 reading submodule-----------------------------------------------
# 
def h5path_to_dataset(hdf5_file):
    def h5_dataset_iterator(h5grp, prefix=''):
        for key in h5grp.keys():
            ipathvar = h5grp[key]
            path = f'{prefix}/{key}'
            if isinstance(ipathvar, h5.Dataset): # test for dataset
                yield (path, ipathvar)
            elif isinstance(ipathvar, h5.Group): # test for group (go down)
                yield from h5_dataset_iterator(ipathvar, path)

    for path, dset in h5_dataset_iterator(hdf5_file):
        yield path

# read hdf5 or hdf4 file
def Read1hdf(hdfname, vars, hdftools_path='', demo=False):
    # if hdf4 file, convert to h5
    if hdfname.endswith('.hdf'):
        h5name = hdfname.replace('.hdf','.h5')
        h5name = h5name.split('/')[-1]
        if hdftools_path != '':
            if not hdftools_path.endswith('/'):
                hdftools_path = hdftools_path+'/'
            if not hdftools_path.endswith('bin/'):
                hdftools_path = hdftools_path+'bin/'
            h4toh5 = hdftools_path+'/'+'h4toh5'
                
        else:
            h4toh5 = 'h4toh5'

        print(h4toh5+' '+hdfname+ ' '+h5name)
        os.system(h4toh5+' '+hdfname+ ' '+h5name)
    elif hdfname.endswith('.h5'):
        h5name = hdfname
    

    vardatas = {}
    # read data 
    if len(vars)>0:
        f0 = h5.File(h5name, 'r')

        # data ranges of date/time & bounding-rectangle box
        for iatt in f0.attrs.keys():
            str_attr = str(f0.attrs[iatt],'utf-8').splitlines()
            blocs = ['RANGEBEGINNINGTIME', 'RANGEENDINGTIME',
                 'RANGEBEGINNINGDATE','RANGEENDINGDATE',
                 'WESTBOUNDINGCOORDINATE','EASTBOUNDINGCOORDINATE',
                 'NORTHBOUNDINGCOORDINATE','SOUTHBOUNDINGCOORDINATE']
        
            in_bloc = False
            for istr in str_attr:
                for ibloc in blocs:
                    if ibloc in istr:
                        if 'END_' not in istr: 
                            in_bloc = True
                            bloc_key = ibloc
                        else:
                            in_bloc = False
                if in_bloc and 'VALUE' in istr:
                    bloc_val = istr.split('=')[-1].strip().replace('"','')
                    if 'TIME' in bloc_key:
                        bloc_val = datetime.strptime(bloc_val,'%H:%M:%S').time()
                    elif 'DATE' in bloc_key:
                        bloc_val = datetime.strptime(bloc_val,'%Y-%m-%d').date()
                    else:
                        bloc_val = float(bloc_val)
                    vardatas[bloc_key] = bloc_val
                    #print(hdfname, ':', iatt, bloc_key,bloc_val)
        
        # dataset
        for ipath in h5path_to_dataset(f0):
        
            varname = ipath.split("/")[-1]
        
            vdata = np.asarray(f0[ipath])
                        
            # vars individually                
            if(vars[0] == 'ALL' or vars[0] == 'all'  # all variables
               or varname in vars):                  # specific variable(s) 
                if varname not in vardatas.keys(): 
                    vardatas[varname] = vdata
                else:
                    vardatas[varname]=np.vstack([vardatas[varname],vdata])
        
        f0.close()
    # end of if (len(vars)>0)
        
    #visualization demo
    if demo and 'Day_CMG_Snow_Cover' in vardatas.keys():
        data = vardatas['Day_CMG_Snow_Cover']
        plt.imshow(data)
        plt.show()
    

    return vardatas

#--------------------------------------------------
# Pre-process modis snowcover h5 datasets, including smoothing

def Prep_modis_snowcov(vars, vardatas, ptxy=[]):
    # INPUTS: vars      - variable names, separated by ','. if 'all' means every var-key in 'vardatas'
    #         vardatas  - python np list, with dicts/data
    #   (optional) ptxy - paired lon/lat (x/y), in [x/lon,y/lat] (only 1 point)
    # OUTPUTS: if 'vars' only has 1 variable
    #                     Output year, doy, lon, lat, data for 'vars' in np.array  
    
    

    # mid of day
    mid_day = vardatas['RANGEBEGINNINGDATE']+ (vardatas['RANGEENDINGDATE']-vardatas['RANGEENDINGDATE'])/2.0
    try:
        nt = len(date2num(mid_day))
    except:
        nt = 1

    
    # grid box
    #['WESTBOUNDINGCOORDINATE','EASTBOUNDINGCOORDINATE',
    #'NORTHBOUNDINGCOORDINATE','SOUTHBOUNDINGCOORDINATE']

    lon_w = np.float(vardatas['WESTBOUNDINGCOORDINATE'])
    lon_e = np.float(vardatas['EASTBOUNDINGCOORDINATE'])
    lat_n = np.float(vardatas['NORTHBOUNDINGCOORDINATE'])
    lat_s = np.float(vardatas['SOUTHBOUNDINGCOORDINATE'])

    DONE_time = False
    
    if vars[0]=='all': 
        vars=vardatas.keys()
        

    # Construct the grid in lat/lon.
    nlat, nlon = vardatas[vars[0]].shape
    dlon = (lon_w - lon_e) / nlon
    dlat = (lat_n - lat_s) / nlat
    xlon = np.linspace(lon_e, lon_e + dlon*nlon, nlon)
    xlat = np.linspace(lat_s, lat_s + dlat*nlat, nlat)

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

    alldata = {}
    for varname in vars:
        #varname = 'Day_CMG_Snow_Cover'

        # write time
            
        if not DONE_time:
            daynums = date2num(mid_day)
            year    = mid_day.year
            doy0    = date(year,1,1)
            doy     = daynums-date2num(doy0)+1
                        
            DONE_time = True
        
        # 
        # appears vardatas are in S-N/E-W ordered, so must be flip over
        data = vardatas[varname]
        data = np.fliplr(data)
        data = np.flipud(data)
        if(len(ptxy)>1): data = data[(iy,ix)] # data space is in N-S/W-E (i.e. lat/lon)
    

        # data re-grouping as follows, so that it can be read in Visit using blue-scale :
        # 0% snow
        # 1-99% snow
        # 100% snow
        # lake ice (107)             ==> water-body (-20)
        # night (111)                ==> undecided  (-10): it shall be what at previous-time or interpolation btw 
        # inland water (237)         ==> water-body (-20)
        # ocean (239)                ==> water-body (-20)
        # cloud-obscured water (250) ==> water-body (-20)
        # data not mapped (253)      ==> missing (-99)
        # fill (255)                 ==> missing (-99)
        data = np.int16(data) # data type is 'uint8', convert to short (othwise cannot be read by Visit)
        data[abs(data-107)<0.1] = -20 # when plot in Visit, negative could be blocked as continuing of non-snow
        data[abs(data-111)<0.1] = -10
        data[abs(data-237)<0.1] = -20
        data[abs(data-239)<0.1] = -20
        data[abs(data-250)<0.1] = -20
        
        data[abs(data-253)<0.1] = -99
        data[abs(data-255)<0.1] = -99
    
        alldata[varname] = deepcopy(data)
    
    return mid_day, year, doy, lon, lat, alldata
    
#-------------------Writing submodule-----------------------------------------------
# Write to geo-referenced CF compliant nc file, if filename given

def Write1GeoNc(vars, vardatas, ptxy=[], ncfname='', newnc=True):
    # INPUTS: vars      - variable names, separated by ','. if 'all' means every var-key in 'vardatas'
    #         vardatas  - python np list, with dicts/data
    #    (optional) ptxy    - paired lon/lat (x/y), in [x/lon,y/lat] (only 1 point)
    #    (optional) ncfname - if not empty, write processed data into (geo)NC file
    #    (optional) newnc   - if not True, write into existed 'ncfname' nc file
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
    try:
        nt = len(date2num(mid_day))
    except:
        nt = 1

    # Construct the grid in lat/lon.
    xlon = vardatas['lon']
    xlat = vardatas['lat']
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
                lon_dim = ncfile.createDimension('lon',  len(lon))
                lat_dim = ncfile.createDimension('lat',  len(lat))
                time_dim= ncfile.createDimension('time', None)

                vlat = ncfile.createVariable('lat', np.float32, ('lat',))
                vlat.units = 'degrees_north'
                vlat.long_name = 'latitude'
                vlat[:] = lat

                vlon = ncfile.createVariable('lon', np.float32, ('lon',))
                vlon.units = 'degrees_east'
                vlon.long_name = 'longitude'
                vlon[:] = lon
        
                # time, create only
                vdaysnum = ncfile.createVariable('daysnum', np.float32, ('time',))
                vdaysnum.units = 'days'
                vdaysnum.long_name = 'days since 0000-01-01 UTC + 1'

                vdate = ncfile.createVariable('date', np.unicode_, ('time',))
                vdate.units = ''
                vdate.long_name = 'date in standard python-datetime calendar'

                vdoy = ncfile.createVariable('doy', np.int16, ('time',))
                vdoy.units = 'day'
                vdoy.long_name = 'day of year'

                vtime = ncfile.createVariable('time', np.float64, ('time',))
                vtime.units = 'day'
                vtime.long_name = 'doy in a year in format yyyydoy'
            
                #global attributes
                ncfile.description = 'daily snow coverage @0.05 degree resolution'
                ncfile.data_source = ('National Snow and Ice Data Center (NSIDC),' +
                    'MODIS/Terra Snow Cover Daily L3 Global 0.05Deg CMG, Version 6 ')
                ncfile.data_citation = ('Hall, D. K. and G. A. Riggs. 2016.' +
                    'MODIS/Terra Snow Cover Daily L3 Global 0.05Deg CMG, Version 6.' +
                    'Boulder, Colorado USA.'+
                    'NASA National Snow and Ice Data Center Distributed Active Archive Center.'+
                    'doi: https://doi.org/10.5067/MODIS/MOD10C1.006.'+
                    '2020-03-15.')
                ncfile.history = '2020-03-19: conversion from h5 format.'
                ncfile.contact = 'F.-M. Yuan, CCSI/ESD-ORNL. yuanf@ornl.gov'
            # done if ncfname !='':
            
            DONE_header = True
        
        # write time
            
        if not DONE_time:
            daynums = date2num(mid_day)
            year    = mid_day.year
            doy0    = date(year,1,1)
            doy     = daynums-date2num(doy0)+1
            ydoy    = year+doy/1000
            
            if ncfname !='':
                if newnc:
                    vdaysnum[0] = daynums
                    vdate[0] = np.unicode_(mid_day)
                    vdoy[0]  = doy
                    vtime[0] = ydoy
                else:
                    vdaysnum = ncfile.variables['daysnum']
                    prv_nt = len(vdaysnum)
                    vdaysnum[prv_nt] = daynums
                
                    vdate = ncfile.variables['date']
                    vdate[prv_nt] = np.unicode_(mid_day)
                
                    vdoy = ncfile.variables['doy']
                    vdoy[prv_nt] = doy

                    vtime = ncfile.variables['time']
                    vtime[prv_nt] = ydoy
            # done write to nc (if ncfname !='':)
            
            DONE_time = True
        
        # 
        # appears vardatas are in S-N/E-W ordered, so must be flip over
        data = vardatas[varname]
        data = np.int16(data) # data type is 'uint8', convert to short (othwise cannot be read by Visit)
        
        if ncfname !='':
            if newnc:
                vtemp = ncfile.createVariable(varname, np.int16, ('time','lat','lon')) # note: unlimited dimension is leftmost
                vtemp.units = '%'
                vtemp.standard_name = varname.strip() # this is a CF standard name
                vtemp.long_name =  varname.strip()+' at 0.05 degree resolution'
                vtemp.Key = "0-100=percent of snow in cell, -10=undecided (fully-night or too short daytime), -20=water body (ocean, inland water/lake), -99=missing (e.g. not mapped, filled data)" ;
        
                vtemp[0,:,:] = np.int16(data)
            else:
                vtemp = ncfile.variables[varname]
                vtemp[prv_nt,:,:] = np.int16(data)
    
    # end of for varname in vars:
    if ncfname !='': ncfile.close()
    
    

#-------------------PLotting submodule-----------------------------------------------
# plotting 1 graph with 1 map

def SingleMap(varname, vardatas, figno=None):
    
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
    
    # appears vardatas are in S-N/E-W ordered, so must be flip over
    data = vardatas[varname]
    data = np.fliplr(data)
    data = np.flipud(data)
    
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



#####################################################################
# module testing
#hdfname = 'MOD10C1.A2019365.006.2020002233131.h5'
#varnames=['all']
#vardatas = \
#    Read1hdf(hdfname, varnames, '/usr/local/gcc-x/h4h5tools')

# mapping
#SingleMap('Day_CMG_Snow_Cover', vardatas)

# write to geo-referenced NC
#ncfname = hdfname.replace('.h5','.nc')
#Write1GeoNc(ncfname, vardatas)


