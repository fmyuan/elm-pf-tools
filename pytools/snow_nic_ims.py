#!/usr/bin/env python

import os, sys, time, math
import glob
import re
import numpy as np
from datetime import datetime, date
from matplotlib.dates import date2num

from optparse import OptionParser

import Modules_nic_ims as nicims
from numpy import long, int16
import netCDF4
from copy import deepcopy

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
parser.add_option("--fileheader", dest="fileheader", default="", \
                  help = "NIC-IMS file header without .asc or .tif ")
parser.add_option("--filetype", dest="ftype", default="ascii", \
                  help = " ascii ('ascii') or geotiff file ('tiff'), default 'ascii' ")
parser.add_option("--res", dest="res", default="4km", \
                  help = " 1km, 4km or 24km, default '4km' ")
parser.add_option("--ptxy", dest="ptxy", default=None, \
                  help = "one point's x/y or lon/lat pair, in [x,y] format, to extract data \
                          default 'None' ")
parser.add_option("--file2nc", dest="file2nc", default=False, \
                  help = " convert an ascii/geotiff file to Netcdf4 file, default: FALSE ", action="store_true")
parser.add_option("--lookup_snowfreeseason", dest="lookup_snowfreeseason", default=False, \
                  help = " lookup snow ending/starting doy in a year and write a Netcdf4 file, default: FALSE ", action="store_true")
parser.add_option("--days_smoothing", dest="smoothingdays", default=7, \
                  help = "smoothing window-size of time, default: 7 days ")

(options, args) = parser.parse_args()

#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)
cwdir = options.workdir

if (options.fileheader == ''):
    print('MUST input file name header by " --fileheader=??? "')
    sys.exit()

# 
if options.ptxy is not None:
    ptxy = re.split(',|:|;| ',options.ptxy)
    ptxy = np.asarray(ptxy,dtype=np.float)
else:
    ptxy = []
##

if not cwdir.endswith('/'): cwdir = cwdir+'/'
fileheader = cwdir+'/'+ options.fileheader
dataformat = options.ftype

if dataformat == 'ascii':
    ftype = 'asc.gz'
elif dataformat == 'geotiff':
    ftype = 'tiff'

res = options.res
varname = 'snowcov'
lonlat2xy = True
minlat = 22.5#None#30.0 # min. latitude to truncate the grids/dataset

alldirfiles = sorted(glob.glob("%s*.%s" % (fileheader, ftype)))
if(len(alldirfiles)<=0):
    sys.exit("No file exists - %s*.*.%s IN %s" %(fileheader, ftype, cwdir))
else:
    print('Total Files: '+str(len(alldirfiles)))

    # convert to geo-referenced NC, with 1 ncfile for each year, without smoothing
if options.file2nc and not options.lookup_snowfreeseason:
    # grids
    if minlat==None:
        if lonlat2xy:
            imsx,imsy, imslons, imslats = \
                nicims.Prep_ims_grid(res, minlat, lonlat2xy=lonlat2xy)
        else:
            imsx,imsy = \
                nicims.Prep_ims_grid(res, minlat, lonlat2xy=lonlat2xy)

    else:
        if lonlat2xy:
            xind, yind, imsx,imsy, imslons, imslats = \
                nicims.Prep_ims_grid(res, minlat, lonlat2xy=lonlat2xy)
        else:
            xind, yind, imsx, imsy = \
                nicims.Prep_ims_grid(res, minlat, lonlat2xy=lonlat2xy)
    
    # daily files
    year_prv = -9999
    for ifile in alldirfiles:

        try: 
            mid_day, year, doy, snowdata = \
                nicims.Prep_ims_snowcov(ifile, fileheader, varname, alldata={})

            # test lon/lat --> geox/y conversion, if Not comment out
            #if lonlat2xy: 
            #    varname = 'ylat'
            #    snowdata[varname] = deepcopy(imslons)
            #    snowdata[varname] = deepcopy(imslats)

        except:
            print (ifile + ' reading Error!')
            continue
        
        
        if minlat!=None:
            snowdata[varname] = snowdata[varname][:,xind,:]
            snowdata[varname] = snowdata[varname][:,:,yind]
            
        snowdata['date'] = [mid_day]
        if lonlat2xy:
            snowdata['geox']  = deepcopy(imsx)
            snowdata['geoy']  = deepcopy(imsy)
        else:
            snowdata['lon']  = deepcopy(imsx)
            snowdata['lat']  = deepcopy(imsy)

        # write NC file(s)
        ncfname = fileheader.split('/')[-1] # remove directory name if any
        ncfname = './'+ncfname+'_'+res+'_yr'+str(year)+'.nc'
        if os.path.isfile(ncfname) and (year==year_prv): 
            nicims.Write1GeoNc([varname], snowdata, ptxy=[], ncfname=ncfname, newnc=False)                    
        else:
            nicims.Write1GeoNc([varname], snowdata, ptxy=[], ncfname=ncfname, newnc=True)

        # save previous year counter
        year_prv = year

        # end of 1 file reading/processing (normally for 1 day here)
        print('DONE with data from: '+ifile)
#done with 'file2nc'

# convert to geo-referenced NC, with 1 ncfile for each year
# processing 1 variable only at daily but for multiple year-doy, 
# e.g. here for snow-free seasons (DOYs for snow melted fully and snow starts)
snow_yearly = {}
if options.lookup_snowfreeseason and (options.ftype=='ascii' or options.ftype=='tif'): 

    # after sorting data, looking-up for DOYs with 'daily snow coverage' < 25% for a year
    crit_snowfree = 0.25 # 0 - 1
    crit_days = int(options.smoothingdays)
    snowcov_days = []
    
    # grids
    if minlat==None:
        imsx,imsy = \
            nicims.Prep_ims_grid(res, minlat, lonlat2xy=lonlat2xy)

    else:
        xind, yind, imsx, imsy = \
            nicims.Prep_ims_grid(res, minlat, lonlat2xy=lonlat2xy)
    
    # daily files
    year_prv = -9999
    for ifile in alldirfiles:
        print('Processing data from: '+ifile+' ......')

        try: 
            mid_day, year, doy, snowcov = \
                nicims.Prep_ims_snowcov(ifile, fileheader, varname, alldata={})

        except:
            print ('FILE: '+ifile + ' reading Error! Skipped')
            continue

        snowcov = np.squeeze(snowcov[varname])
        if minlat!=None:
            snowcov = snowcov[xind,:]
            snowcov = snowcov[:,yind]

        #------------------------------------------------------------------------------------------------
        # if have to do moving-window smoothing, prapare a few days' data
        if ifile==alldirfiles[0]:
            snowcov_days = np.full((crit_days,)+snowcov.shape,np.float(-99))
            snowcov_days[0:crit_days,] = snowcov[:,:]
        snowcov_days = np.delete(snowcov_days,0,0)  # remove [0,], i.e. oldest
        snowcov_days = np.vstack((snowcov_days, np.reshape(snowcov,(1,)+snowcov.shape))) # append the newest

        # smoothing, upon previous day's coverage ('temp_snowcov_prv')
        if ifile==alldirfiles[0]:
            temp_snowcov     = snowcov.astype(float)
            temp_snowcov_prv = snowcov.astype(float)
        else:
            #
            ij = np.where(snowcov_days<0)
            if len(snowcov_days[ij])>1: 
                snowcov_days[ij]=np.nan # excluding non-snow-extent data
                
                # in general, taking max. snow coverage in moving-window
                n_prv = np.int(crit_days/2)
                temp_snowcov = np.nanmax(snowcov_days,0)
                
                # previously non-snowed cells, re-assign mean value in them so that smoothing occasional outliers
                prv1 = np.nanmean(snowcov_days[0:n_prv,],0)
                prv2 = np.nanmean(snowcov_days[n_prv+1:,],0)
                #ij = np.where( ((prv1-prv2)>=0) | 
                #              ((prv1<crit_snowfree) & (prv2<crit_snowfree)) )
                ij = np.where(prv2<crit_snowfree)
                if len(temp_snowcov[ij])>0: 
                    temp_snowcov[ij] = np.nanmin(snowcov_days,0)[ij]

                # since 'np.nanxxx' operations with 'snowcov_days', need to reassign those missing values
                # otherwise '0' would be assigned to those cells, which causes fake values for 'nights'
                ij = np.where(np.isnan(temp_snowcov))
                if len(temp_snowcov[ij])>0: 
                    temp_snowcov[ij] = temp_snowcov_prv[ij]
                    
        #------------------------------------------------------------------------------------------------
        
        # write smoothed 'temp_snowcov' to NC
        if options.file2nc:
            ncdata = {}
            ncdata['date'] = [mid_day]
            if lonlat2xy:
                ncdata['geox']  = deepcopy(imsx)
                ncdata['geoy']  = deepcopy(imsy)
            else:
                ncdata['lon']  = deepcopy(imsx)
                ncdata['lat']  = deepcopy(imsy)
            ncdata[varname] = deepcopy(temp_snowcov)

            # write NC file(s)
            ncfname = fileheader.split('/')[-1] # remove directory name if any
            ncfname = './'+ncfname+'_'+res+'_yr'+str(year)+'.nc'
            if os.path.isfile(ncfname) and (year==year_prv): 
                nicims.Write1GeoNc([varname], ncdata, ptxy=[], ncfname=ncfname, newnc=False)                    
            else:
                nicims.Write1GeoNc([varname], ncdata, ptxy=[], ncfname=ncfname, newnc=True)
            
        #------------------------------------------------------------------------------------------------
         
        if ifile==alldirfiles[0]:
            snow_start = np.full(snowcov.shape,np.int16(-99))
            snow_end   = np.full(snowcov.shape,np.int16(-99))
            year_prv   = year
            snow_free  = np.full(snowcov.shape,np.int16(0))
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
                ij = np.where((temp>=temp_crit) & (snow_end>0))
                if len(snow_end[ij])>0: snow_end[ij] = 0
                ij = np.where((temp>=temp_crit) & (snow_start>0))
                if len(snow_start[ij])>0: snow_start[ij] = 0
                ij = np.where((temp>=temp_crit))
                if len(snow_free[ij])>0: snow_free[ij] = 365
                
                 # never-snowfree cells
                ij = np.where((snow_free<-temp_crit)  & (snow_end>0))
                if len(snow_end[ij])>0: snow_end[ij] = 213
                ij = np.where((snow_free<-temp_crit)  & (snow_start>0))
                if len(snow_start[ij])>0: snow_start[ij] = 1
                ij = np.where((snow_free<-temp_crit)  & (snow_free>0))
                if len(snow_free[ij])>0: snow_free[ij] = 0
                
                # sea/out-of-bound cells
                ij = np.where((temp_snowcov==-1)) # sea-iced
                if len(snow_end[ij])>0: snow_end[ij] = -10
                if len(snow_start[ij])>0: snow_start[ij] = -10
                if len(snow_free[ij])>0: snow_free[ij] = -10
                ij = np.where((temp_snowcov==-2)) # sea-iced
                if len(snow_end[ij])>0: snow_end[ij] = -20
                if len(snow_start[ij])>0: snow_start[ij] = -20
                if len(snow_free[ij])>0: snow_free[ij] = -20
                ij = np.where((temp_snowcov==-99))#out-of-bounds
                if len(snow_end[ij])>0: snow_end[ij] = temp_snowcov[ij]
                if len(snow_start[ij])>0: snow_start[ij] = temp_snowcov[ij]
                if len(snow_free[ij])>0: snow_free[ij] = temp_snowcov[ij]
                 
            # save yearly data
            if (len(snow_yearly)<=0):
                if lonlat2xy:
                    snow_yearly['geox']  = deepcopy(imsx)
                    snow_yearly['geoy']  = deepcopy(imsy)
                else:
                    snow_yearly['lon']  = deepcopy(imsx)
                    snow_yearly['lat']  = deepcopy(imsy)
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
        if ifile==alldirfiles[0]:
            ij = np.where( (temp_snowcov<=crit_snowfree) & (temp_snowcov>=0) )        #snow-free non-water cells
        else:
            ij = np.where( ((temp_snowcov<=crit_snowfree) & (temp_snowcov>=0)) &      #snow-free non-water cells
                           ((temp_snowcov_prv>crit_snowfree) & (temp_snowcov_prv>=0)) )               #snow-covered previously
        if len(snow_end[ij])>0 and doy<213: # snow-melting ends must be not over end of July
            snow_end[ij] = doy
            snow_free[ij]= 1

        
        # snow-covering starts
        if ifile==alldirfiles[0]:
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
        if ifile==alldirfiles[-1]: 
                
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
                ij = np.where((temp>=temp_crit) & (snow_end>0))
                if len(snow_end[ij])>0: snow_end[ij] = 0
                ij = np.where((temp>=temp_crit) & (snow_start>0))
                if len(snow_start[ij])>0: snow_start[ij] = 0
                ij = np.where((temp>=temp_crit))
                if len(snow_free[ij])>0: snow_free[ij] = 365
                
                 # never-snowfree cells
                ij = np.where((snow_free<-temp_crit)  & (snow_end>0))
                if len(snow_end[ij])>0: snow_end[ij] = 213
                ij = np.where((snow_free<-temp_crit)  & (snow_start>0))
                if len(snow_start[ij])>0: snow_start[ij] = 1
                ij = np.where((snow_free<-temp_crit)  & (snow_free>0))
                if len(snow_free[ij])>0: snow_free[ij] = 0
                
                # sea/out-of-bound cells
                ij = np.where((temp_snowcov==-1)) # sea-iced
                if len(snow_end[ij])>0: snow_end[ij] = -10
                if len(snow_start[ij])>0: snow_start[ij] = -10
                if len(snow_free[ij])>0: snow_free[ij] = -10
                ij = np.where((temp_snowcov==-2)) # sea
                if len(snow_end[ij])>0: snow_end[ij] = -20
                if len(snow_start[ij])>0: snow_start[ij] = -20
                if len(snow_free[ij])>0: snow_free[ij] = -20
                ij = np.where((temp_snowcov==-99))#out-of-bounds
                if len(snow_end[ij])>0: snow_end[ij] = temp_snowcov[ij]
                if len(snow_start[ij])>0: snow_start[ij] = temp_snowcov[ij]
                if len(snow_free[ij])>0: snow_free[ij] = temp_snowcov[ij]
 
            # save into yearly-data
            if (len(snow_yearly)<=0):
                if lonlat2xy:
                    snow_yearly['geox'] = deepcopy(imsx)
                    snow_yearly['geoy'] = deepcopy(imsy)
                else:
                    snow_yearly['lon']  = deepcopy(imsx)
                    snow_yearly['lat']  = deepcopy(imsy)
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
        # end of 1 file reading/processing (normally for 1 day here)
        
    # end of 'for ifile in alldirfiles (i.e. all files)
#Done with 'lookup_snowfreeseason' option for daily 'ascii' or 'tiff' data

#------------------------------------------------------------------------------------------------
    

# write yearly snow-free season data into 1 NC file
if len(snow_yearly)>0:
    
    del snow_end
    del snow_start
    del snowcov
    del snowcov_days
    del temp_snowcov
    del temp_snowcov_prv

    #-----------------------------------------------------------
    yr0 = int(snow_yearly['Year'][0])
    yr1 = int(snow_yearly['Year'][-1])
    # write all snow_yearly data to NC file
    ncfname = 'NSIDC_yearly_snowfree.'+str(yr0)+'-'+str(yr1)+'.nc'
    if os.path.isfile(ncfname): os.system('rm -rf '+ncfname)
    ncfile = netCDF4.Dataset(ncfname, mode='w',format='NETCDF4') 
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
    vdoy_end.Key = "0-365 = land with/without snow, -10 = sea with ice, -20 = sea, -99 = beyond map coverage" ;
    vdata = snow_yearly['Snow_ending']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdoy_end[:,:,:] = vdata

    if 'lat' in snow_yearly.keys():
        vdoy_start = ncfile.createVariable('Doy_snowcovered', np.int16, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vdoy_start = ncfile.createVariable('Doy_snowcovered', np.int16, ('time','geoy','geox'))
    vdoy_start.units = 'doy'
    vdoy_start.long_name = 'day of year when snow starts covered on ground'
    vdoy_start.Key = "0-365 = land with/without snow, -10 = sea with ice, -20 = sea, -99 = beyond map coverage" ;
    vdata = snow_yearly['Snow_starting']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdoy_start[:,:,:] = vdata

    if 'lat' in snow_yearly.keys():
        vdays_free = ncfile.createVariable('Days_snowfree', np.int16, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vdays_free = ncfile.createVariable('Days_snowfree', np.int16, ('time','geoy','geox'))
    vdays_free.units = 'days'
    vdays_free.long_name = 'days in year from snow ends  to starts covered on ground'
    vdays_free.Key = "0-365 = land with/without snow, -10 = sea with ice, -20 = sea, -99 = beyond map coverage" ;
    vdata = snow_yearly['Snow_freedays']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdays_free[:,:,:] = vdata
            
    #global attributes
    #global attributes
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
    
    