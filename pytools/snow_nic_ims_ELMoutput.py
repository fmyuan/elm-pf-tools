#!/usr/bin/env python

import os, sys, time, math
import glob
import re
import numpy as np
from datetime import datetime, date
from matplotlib.dates import date2num, num2date

from optparse import OptionParser

import Modules_nic_ims as nicims
from numpy import long, int16
from netCDF4 import Dataset
from copy import deepcopy

from Modules_CLMoutput_nc4 import CLM_NcRead_1simulation
from Modules_CLMoutput_nc4 import CLMvar_1Dtseries

#------------------------------------------------------------------------------------------------------------------
# after sorting data, looking-up for DOYs with 'daily snow coverage' < 25% for a year
def seeking_YRLYsnowfreeseason(daynums, snowcov, no_leap=False, snow_yearly={}):
    crit_snowfree = 0.25 # 0 - 1

    # timely-ordered data ('snowcov')
    year_prv = -9999
    for it in range(len(daynums)):
        if no_leap:
            year = int(daynums[it]/365.0)
            doy  = daynums[it] - year*365.0 + 1.0
        else:
            mid_day = num2date(daynums[it]).date()
            year = mid_day.year
            doy  = daynums[it] - date2num(date(year,1,1)) + 1

        temp_snowcov = snowcov[it,]

        #------------------------------------------------------------------------------------------------
         
        if it==0:
            snow_start = np.full(temp_snowcov.shape,np.int16(-99))
            snow_end   = np.full(temp_snowcov.shape,np.int16(-99))
            year_prv   = year
            snow_free  = np.full(temp_snowcov.shape,np.int16(0))
            totdays_yr = 0
        
        #------------------------------------------------------------------------------------------------
        # save previous-year data when a new year starts, note this will not be for end of 'daynums'
        if year != year_prv: 
            
            # replace '-xx' with those same cells from snowcov (i.e. non-land)
            ij = np.where((snow_end<0) & (temp_snowcov<0))
            if len(snow_end[ij])>0: snow_end[ij] = temp_snowcov.astype(int16)[ij]
            
            ij = np.where((snow_start<0) & (temp_snowcov<0))
            if len(snow_start[ij])>1: snow_start[ij] = temp_snowcov.astype(int16)[ij]

            # have to flag never-snow/snowfree cells (snow_free within 'totdays_yr-30d' for a full year)
            if totdays_yr>0:
                temp_crit = np.int(-30.0*np.float(totdays_yr)/366.0)
                
                # never-snow cells
                temp = snow_free - totdays_yr
                ij = np.where((temp>=temp_crit))
                if len(snow_end[ij])>0: snow_end[ij] = 0
                ij = np.where((temp>=temp_crit))
                if len(snow_start[ij])>0: snow_start[ij] = 0
                ij = np.where((temp>=temp_crit))
                if len(snow_free[ij])>0: snow_free[ij] = 365
                
                 # never-snowfree cells
                ij = np.where((snow_free<-temp_crit))
                if len(snow_end[ij])>0: snow_end[ij] = 365
                ij = np.where((snow_free<-temp_crit))
                if len(snow_start[ij])>0: snow_start[ij] = 1
                ij = np.where((snow_free<-temp_crit)  & (snow_free>0))
                if len(snow_free[ij])>0: snow_free[ij] = 0
                
                # sea/out-of-bound cells
                ij = np.where((temp_snowcov==-1)) # sea-iced
                if len(snow_end[ij])>0: snow_end[ij] = -22
                if len(snow_start[ij])>0: snow_start[ij] = -22
                if len(snow_free[ij])>0: snow_free[ij] = -22
                ij = np.where((temp_snowcov==-2)) # sea-iced
                if len(snow_end[ij])>0: snow_end[ij] = -55
                if len(snow_start[ij])>0: snow_start[ij] = -55
                if len(snow_free[ij])>0: snow_free[ij] = -55
                ij = np.where((temp_snowcov==-99))#out-of-bounds
                if len(snow_end[ij])>0: snow_end[ij] = -99
                if len(snow_start[ij])>0: snow_start[ij] = -99
                if len(snow_free[ij])>0: snow_free[ij] = -99
                 
            # save yearly data
            if (len(snow_yearly)<=0):
                snow_yearly['Year'] = np.empty(1)
                snow_yearly['Snow_ending']   = np.empty((1,)+snow_end.shape)
                snow_yearly['Snow_starting'] = np.empty((1,)+snow_start.shape)
                snow_yearly['Snow_freedays'] = np.empty((1,)+snow_free.shape)
                
                snow_yearly['Year'][0]           = year_prv
                temp = np.reshape(snow_end,(1,)+snow_end.shape)
                snow_yearly['Snow_ending'][0,]   = deepcopy(temp) # to avoid assign-values by reference
                temp = np.reshape(snow_start,(1,)+snow_start.shape)
                snow_yearly['Snow_starting'][0,] = deepcopy(temp)
                temp = np.reshape(snow_free,(1,)+snow_free.shape)
                snow_yearly['Snow_freedays'][0,] = deepcopy(temp)
        
            else:
                
                snow_yearly['Year'] = np.vstack((snow_yearly['Year'], year_prv))
                temp = np.reshape(snow_end,(1,)+snow_end.shape)
                snow_yearly['Snow_ending'] = np.vstack((snow_yearly['Snow_ending'], deepcopy(temp)))
                temp = np.reshape(snow_start,(1,)+snow_start.shape)
                snow_yearly['Snow_starting'] = np.vstack((snow_yearly['Snow_starting'], deepcopy(temp)))
                temp = np.reshape(snow_free,(1,)+snow_free.shape)
                snow_yearly['Snow_freedays'] = np.vstack((snow_yearly['Snow_freedays'], deepcopy(temp)))
            
            snow_free[:,:] = 0 # reset to 0 for accumulation
            totdays_yr = 0
            snow_start[:,:] = -99
            snow_end[:,:] = -99
        #------------------------------------------------------------------------------------------------
        

        # snow-melting ends
        if it==0:
            ij = np.where( (temp_snowcov<=crit_snowfree) & (temp_snowcov>=0) )        #snow-free non-water cells
        else:
            ij = np.where( ((temp_snowcov<=crit_snowfree) & (temp_snowcov>=0)) &      #snow-free non-water cells
                           ((temp_snowcov_prv>crit_snowfree) & (temp_snowcov_prv>=0)) )               #snow-covered previously
        if len(snow_end[ij])>0 and doy<213: # snow-melting ends must be not over end of July
            snow_end[ij] = doy
            snow_free[ij]= 1

        
        # snow-covering starts
        if it==0:
            ij = np.where( (temp_snowcov>=crit_snowfree) & (temp_snowcov>=0) )            #snow-covered cells
        else:
            ij = np.where( ((temp_snowcov>=crit_snowfree) & (temp_snowcov>=0)) &          #snow-covered cells
                           ((temp_snowcov_prv<crit_snowfree) & (temp_snowcov_prv>=0)) &   #snow-free previously
                           ((snow_start-snow_end)<=0) )                                     #snow-melting done
        if len(snow_start[ij])>0: 
            snow_start[ij] = doy

        # snow-free days in a year (need this to flag never-snowed cells)
        ij = np.where( (temp_snowcov<crit_snowfree) & (temp_snowcov>=0) &                  # snow-free cells
                       ((snow_start-snow_end)<=0) )                                        # snow-covering not yet
        if len(snow_free[ij])>0: snow_free[ij] = snow_free[ij] + 1
        totdays_yr = totdays_yr + 1

         #------------------------------------------------------------------------------------------------
        # save data at ending of time-steps, otherwise at a new-year starts (see above)
        if it==len(daynums)-1: 
                
            # replace '-xx' with those same cells from snowcov (indicating water bodies or missings)
            ij = np.where((snow_end<0) & (temp_snowcov<0))
            if len(snow_end[ij])>0: snow_end[ij] = temp_snowcov.astype(int16)[ij]
            
            ij = np.where((snow_start<0) & (temp_snowcov<0))
            if len(snow_start[ij])>1: snow_start[ij] = temp_snowcov.astype(int16)[ij]

            # have to flag never-snow/snowfree cells (snow_free within 'totdays_yr-30d' for a full year)
            if totdays_yr>0:
                temp_crit = np.int(-30.0*np.float(totdays_yr)/366.0)
                    
                # never-snow cells
                temp = snow_free - totdays_yr
                ij = np.where((temp>=temp_crit))
                if len(snow_end[ij])>0: snow_end[ij] = 0
                ij = np.where((temp>=temp_crit))
                if len(snow_start[ij])>0: snow_start[ij] = 0
                ij = np.where((temp>=temp_crit))
                if len(snow_free[ij])>0: snow_free[ij] = 365
                
                 # never-snowfree cells
                ij = np.where((snow_free<-temp_crit))
                if len(snow_end[ij])>0: snow_end[ij] = 365
                ij = np.where((snow_free<-temp_crit))
                if len(snow_start[ij])>0: snow_start[ij] = 1
                ij = np.where((snow_free<-temp_crit)  & (snow_free>0))
                if len(snow_free[ij])>0: snow_free[ij] = 0
                
                # sea/out-of-bound cells
                ij = np.where((temp_snowcov==-1)) # sea-iced
                if len(snow_end[ij])>0: snow_end[ij] = -22
                if len(snow_start[ij])>0: snow_start[ij] = -22
                if len(snow_free[ij])>0: snow_free[ij] = -22
                ij = np.where((temp_snowcov==-2)) # sea
                if len(snow_end[ij])>0: snow_end[ij] = -55
                if len(snow_start[ij])>0: snow_start[ij] = -55
                if len(snow_free[ij])>0: snow_free[ij] = -55
                ij = np.where((temp_snowcov==-99))#out-of-bounds
                if len(snow_end[ij])>0: snow_end[ij] = -99
                if len(snow_start[ij])>0: snow_start[ij] = -99
                if len(snow_free[ij])>0: snow_free[ij] = -99
 
            # save into yearly-data
            if (len(snow_yearly)<=0):
                snow_yearly['Year'] = np.empty(1)
                snow_yearly['Snow_ending']   = np.empty((1,)+snow_end.shape)
                snow_yearly['Snow_starting'] = np.empty((1,)+snow_start.shape)
                snow_yearly['Snow_freedays'] = np.empty((1,)+snow_free.shape)
            
                snow_yearly['Year'][0]           = year
                temp = np.reshape(snow_end,(1,)+snow_end.shape)
                snow_yearly['Snow_ending'][0,]   = deepcopy(temp) # to avoid assign-values by reference
                temp = np.reshape(snow_start,(1,)+snow_start.shape)
                snow_yearly['Snow_starting'][0,] = deepcopy(temp)
                temp = np.reshape(snow_free,(1,)+snow_free.shape)
                snow_yearly['Snow_freedays'][0,] = deepcopy(temp)
        
            else:
                
                snow_yearly['Year'] = np.vstack((snow_yearly['Year'], year))
                temp = np.reshape(snow_end,(1,)+snow_end.shape)
                snow_yearly['Snow_ending'] = np.vstack((snow_yearly['Snow_ending'], deepcopy(temp)))
                temp = np.reshape(snow_start,(1,)+snow_start.shape)
                snow_yearly['Snow_starting'] = np.vstack((snow_yearly['Snow_starting'], deepcopy(temp)))
                temp = np.reshape(snow_free,(1,)+snow_free.shape)
                snow_yearly['Snow_freedays'] = np.vstack((snow_yearly['Snow_freedays'], deepcopy(temp)))
            snow_free[:,:] = 0 # reset to 0 for accumulation
            totdays_yr = 0
         #------------------------------------------------------------------------------------------------
 
        # save previous day for filling night/missing values at next day
        year_prv = year
        temp_snowcov_prv = deepcopy(temp_snowcov)
        # end for 1 day here
        
    # 
    return snow_yearly
    # 

#------------------------------------------------------------------------------------------------------------------

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
parser.add_option("--imsheader", dest="imsheader", default="", \
                  help = "NIC-IMS output Netcdf file header with path but no .nc ")
parser.add_option("--elmheader", dest="elmheader", default="", \
                  help = "ELM output Netcdf file header with path but no .nc ")
parser.add_option("--res", dest="res", default="24km", \
                  help = " 1km, 4km or 24km for IMS dataset, default '24km' ")
parser.add_option("--ptxy", dest="ptxy", default=None, \
                  help = "one point's x/y or lon/lat pair, in [x,y] format, to extract data \
                          default 'None' ")
parser.add_option("--clmout_timestep", dest="clmout_ts", default="daily", \
                  help="clm output variable timestep (default = 'daily', other option monthly)")
parser.add_option("--startyr", dest="startyr", default="1", \
                  help="clm output starting year to plot (default = 1, i.e. first year of simulation" \
                   " and can be user-defined)")
parser.add_option("--endyr", dest="endyr", default="", \
                  help="clm output ending year to plot (default = none, i.e. end of simulation)")
parser.add_option("--lookup_snowfreeseason", dest="lookup_snowfreeseason", default=False, \
                  help = " lookup snow ending/starting doy in a year and write a Netcdf4 file, default: FALSE ", action="store_true")

(options, args) = parser.parse_args()

#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)
cwdir = options.workdir

if (options.elmheader == '' and options.imsheader == ''):
    print('MUST input file name header, including fullpath, by " --elmheader=??? or --imsheader"')
    sys.exit()

# 
if options.ptxy is not None:
    ptxy = re.split(',|:|;| ',options.ptxy)
    ptxy = np.asarray(ptxy,dtype=np.float)
else:
    ptxy = []
##

startdays = (int(options.startyr)-1)*365.0
if(options.clmout_ts=='daily'): 
    startdays=startdays+1.0
elif(options.clmout_ts=='monthly'): 
    startdays=startdays-1.0
    
enddays = -9999
if(options.endyr !=""): enddays = int(options.endyr)*365.0


if not cwdir.endswith('/'): cwdir = cwdir+'/'

res = options.res
ims_varname = 'snowcov'
elm_varname = 'FSNO'

ftype = 'nc'

#------------------------------------------------------------------------------
snow_yearly = {}
# ELM snow coverages reading
if (options.elmheader != ""):
    elmpathfileheader = options.elmheader

    alldirfiles = sorted(glob.glob("%s*.%s" % (elmpathfileheader, ftype)))
    if(len(alldirfiles)<=0):
        sys.exit("No file exists - %s*.*.%s IN %s" %(elmpathfileheader, ftype, cwdir))
    else:
        print('Total Files: '+str(len(alldirfiles)))

    ncfileheader = elmpathfileheader.split('/')[-1]
    elm_odir   = elmpathfileheader.replace(ncfileheader,'')
    if(elm_odir.strip()==''):elm_odir='./'
    elmfincl = 'h0'
    if ('.h1.' in alldirfiles[0]): elmfincl='h1'
    
    # read-in datasets from 1 simulation
    ny, nx, nlgrnd, nldcmp, ncolumn, npft, varsdata, varsdims = \
        CLM_NcRead_1simulation(elm_odir, \
                           ncfileheader, \
                           elmfincl, \
                           False, \
                           [elm_varname], \
                           startdays, enddays, \
                           False)

    if2dgrid = True
    if(len(varsdata['topo'].shape)==1): if2dgrid = False
    # names for variable and its time-axis
    vars_list = list(varsdata.keys())    # var_list is in format of 'h*_varname', in which 'varname' is the real variable names
    for hv in vars_list:
        if re.search(elm_varname, hv): 
            var_h = hv
            hinc  = hv.replace(elm_varname,'')
            var_t = '%stime' %hinc
            break
    vdata = varsdata[var_h]
    vdims = varsdims[var_h]
    tt = varsdata[var_t]   # time unit: days
    # processing original data
    t, gdata, sdata, zdim_indx, pdim_indx = \
        CLMvar_1Dtseries(tt, vdims, vdata, nx, ny, -999, -999, if2dgrid, \
                    False, False)

    # snow relevant ELM output usually is grid-wised
    if(if2dgrid):
        snowcov = np.reshape(gdata,(len(t),ny,nx))
    else:
        snowcov = gdata
    landmask = varsdata['landmask']
    elmsea_ij = np.where(landmask==0)
    
    # time is in unit of days since '1850-01-01 00:00:00', without leap-year (SO CANNOT use date/time functions from python)
    t0 = 1850*365.0
    daynums = [x+t0 for x in t]
    
    # lon/lat of ELM output
    lonlat2xy = False
    elmx = varsdata['lon']
    elmy = varsdata['lat']

    # looking-up snow-free seasons over years
    if (options.lookup_snowfreeseason):
        snow_yearly = \
          seeking_YRLYsnowfreeseason(daynums, snowcov, no_leap=True)
    for it in range(len(snow_yearly['Year'])):
        snow_yearly['Snow_freedays'][it][elmsea_ij]=-55
        snow_yearly['Snow_starting'][it][elmsea_ij]=-55
        snow_yearly['Snow_ending'][it][elmsea_ij]=-55
    
# NIC-IMS snow coverages reading from daily NC snowcov dataset
if (options.imsheader != ""):
    imspathfileheader = options.imsheader

    alldirfiles = sorted(glob.glob("%s*.%s" % (imspathfileheader, ftype)))
    if(len(alldirfiles)<=0):
        sys.exit("No file exists - %s*.*.%s IN %s" %(imspathfileheader, ftype, cwdir))
    else:
        print('Total Files: '+str(len(alldirfiles)))

    ncfileheader = imspathfileheader.split('/')[-1]
    ims_odir   = imspathfileheader.replace(ncfileheader,'')
    if(ims_odir.strip()==''):ims_odir='./'

    for ncfile in alldirfiles:
        f = Dataset(ncfile,'r')
        
        if ncfile == alldirfiles[0]:
            imsx = np.asarray(f.variables['geox'])
            imsy = np.asarray(f.variables['geoy'])
            lonlat2xy = True
        
        # read-in datasets from one or multiple nc file 
        tt = np.asarray(f.variables['daysnum']) #units - days since date/time 0000-00-00 00:00:00
        vdata = np.asarray(f.variables[ims_varname])
        vdata = np.flip(vdata,1) # weired that have to flip data (TODO checking)
        vdata = np.flip(vdata,2)

        # 
        if ncfile == alldirfiles[0]:
            daynums = deepcopy(tt)
            snowcov = deepcopy(vdata)
        else:
            daynums = np.hstack((daynums, tt))
            snowcov = np.vstack((snowcov, vdata))
    # done with all nc files

    # looking-up snow-free seasons over years
    if (options.lookup_snowfreeseason):
        snow_yearly = \
          seeking_YRLYsnowfreeseason(daynums, snowcov)


#------------------------------------------------------------------------------------------------
    

# write yearly snow-free season data into 1 NC file
if len(snow_yearly)>0:

    if lonlat2xy:
        snow_yearly['geox'] = deepcopy(imsx)
        snow_yearly['geoy'] = deepcopy(imsy)
    else:
        snow_yearly['lon']  = deepcopy(elmx)
        snow_yearly['lat']  = deepcopy(elmy)

    
    #-----------------------------------------------------------
    # write all snow_yearly data to NC file
    if(options.imsheader != "" and options.elmheader==""):
        ncfname = 'NSIDC_yearly_snowfree.nc'
    elif(options.elmheader != "" and options.imsheader==""):
        ncfname = 'ELM20181101_N60_yearly_snowfree.nc'

    if os.path.isfile(ncfname): os.system('rm -rf '+ncfname)
    ncfile = Dataset(ncfname, mode='w',format='NETCDF4') 
    print('Create and Write NC file: '+ncfname)

    yr_dim  = ncfile.createDimension('time', None)
    if 'lat' in snow_yearly.keys():
        lon = snow_yearly['lon']
        lat = snow_yearly['lat']
        lon_dim = ncfile.createDimension('lon',  len(lon))
        lat_dim = ncfile.createDimension('lat',  len(lat))
        vlon = ncfile.createVariable('lon', np.float32, ('lon',))
        vlon.units = 'degrees_east'
        vlon.long_name = 'longitude'
        vlat = ncfile.createVariable('lat', np.float32, ('lat',))
        vlat.units = 'degrees_north'
        vlat.long_name = 'latitude'
    elif 'geox' in snow_yearly.keys():
        lon = snow_yearly['geox']
        lat = snow_yearly['geoy']
        lon_dim = ncfile.createDimension('geox',  len(lon))
        lat_dim = ncfile.createDimension('geoy',  len(lat))
        vlat = ncfile.createVariable('geox', np.float32, ('geox',))
        vlat.units = 'meters'
        vlat.long_name = 'Northing (Polar stereographic ellipsoidal projection)'
        vlon = ncfile.createVariable('geoy', np.float32, ('geoy',))
        vlon.units = 'meters'
        vlon.long_name = 'Easting (Polar stereographic ellipsoidal projection)'

    vlat[:] = lat
    vlon[:] = lon

    vtime = ncfile.createVariable('time', np.int16, ('time',))
    vtime.units = 'year'
    vtime.long_name = 'year'
    vtime[:] = snow_yearly['Year']

    if 'lat' in snow_yearly.keys():
        vdoy_end = ncfile.createVariable('Doy_snowmelted', np.int16, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vdoy_end = ncfile.createVariable('Doy_snowmelted', np.int16, ('time','geoy','geox'))
    vdoy_end.units = 'doy'
    vdoy_end.long_name = 'day of year when snow fully melted on ground'
    vdoy_end.Key = "0-365 = land with/without snow, -22 = sea with ice, -55 = sea, -99 = beyond map coverage" ;
    vdata = snow_yearly['Snow_ending']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdoy_end[:,:,:] = vdata

    if 'lat' in snow_yearly.keys():
        vdoy_start = ncfile.createVariable('Doy_snowcovered', np.int16, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vdoy_start = ncfile.createVariable('Doy_snowcovered', np.int16, ('time','geoy','geox'))
    vdoy_start.units = 'doy'
    vdoy_start.long_name = 'day of year when snow starts covered on ground'
    vdoy_start.Key = "0-365 = land with/without snow, -22 = sea with ice, -55 = sea, -99 = beyond map coverage" ;
    vdata = snow_yearly['Snow_starting']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdoy_start[:,:,:] = vdata

    if 'lat' in snow_yearly.keys():
        vdays_free = ncfile.createVariable('Days_snowfree', np.int16, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vdays_free = ncfile.createVariable('Days_snowfree', np.int16, ('time','geoy','geox'))
    vdays_free.units = 'days'
    vdays_free.long_name = 'days in year from snow ends  to starts covered on ground'
    vdays_free.Key = "0-365 = land with/without snow, -22 = sea with ice, -55 = sea, -99 = beyond map coverage" ;
    vdata = snow_yearly['Snow_freedays']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdays_free[:,:,:] = vdata
            
    #
    #global attributes
    if(options.elmheader !=""):
        ncfile.description = 'ELM simulated No60 and above Northern High-latitude Region yearly snow-cover gone/start day and snowfree length @0.5deg resolution '
        ncfile.data_source = ('Offline E3SM Land Model, ' +
        'forced by GSWP3 v2, half-degree, over >=N60, fully CNP coupled' +
        'Master @ 2018-11-01 ')
        ncfile.history = '2020-05-06: calculated from  daily snow coverage FSNO simulations.'
        
    if(options.imsheader!=""):
        if (res=='1km'):
            ncfile.description = 'Northern Hemisphere yearly snow-cover gone/start day and snowfree length @1km resolution '
        elif (res=='24km'):
            ncfile.description = 'Northern Hemisphere daily snow-cover gone/start day and snowfree length  @24km resolution '
        elif (res=='4km'):
            ncfile.description = 'Northern Hemisphere daily snow-cover  gone/start day and snowfree length @4km resolution '
        ncfile.data_source = ('US National Ice Center (USNIC), ' +
        'Interactive Multsensor Snow and Ice Mapping Service (IMS),' +
        'Version 1.2 ')
        ncfile.data_citation = ('U.S. National Ice Center. 2008. Updated daily. ' +
        'IMS Daily Northern Hemisphere Snow and Ice Analysis at 1 km, 4 km, and 24 km Resolutions, Version 1.2 [4km data].' +
        'Boulder, Colorado USA. '+
        'NSIDC: National Snow and Ice Data Center. '+
        'doi: https://doi.org/10.7265/N52R3PMC. '+
        '2020-03-24.')
        ncfile.history = '2020-04-14: conversion from ASCII format.'
    
    ncfile.contact = 'F.-M. Yuan, CCSI/ESD-ORNL. yuanf@ornl.gov'
 
    ncfile.close()
    
    