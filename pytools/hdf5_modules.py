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

import netCDF4
from matplotlib.dates import date2num
from builtins import str
from copy import deepcopy

#--------------------------------------------------

#-------------------hdf5/hdf4 reading submodule-----------------------------------------------
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
        

    return vardatas

#--------------------------------------------------
# Pre-process modis snowcover h5 datasets, including smoothing

   
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
    except Exception as e:
        print(e)
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


#####################################################################
# module testing
#hdfname = 'CESM-RCP8.5-2006-2100_dm1985-2015.h5'
#varnames=['all']
#vardatas = \
#    Read1hdf(hdfname, varnames)

# write to geo-referenced NC
#ncfname = hdfname.replace('.h5','.nc')
#Write1GeoNc(ncfname, vardatas)


