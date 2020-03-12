#!/usr/bin/env python

import os, sys, time, math
import glob
import re
import numpy as np
from datetime import datetime, date
from matplotlib.dates import date2num

from optparse import OptionParser

import Modules_eos_hdf as eoshdf
from numpy import long, int16
import netCDF4
from copy import deepcopy


#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
parser.add_option("--eosh4file", dest="eosh4file", default="", \
                  help = "NASA EOS hdf4 file name without .hdf ")
parser.add_option("--hdftype", dest="hdftype", default="hdf", \
                  help = " hdf4('hdf') or hdf5 file ('h5'), default 'hdf' ")
parser.add_option("--varname", dest="vars", default="Day_CMG_Snow_Cover", \
                  help = "variable name(s), ':' separated lists, or, with 'ALL' for all vars to be processed \
                          default 'Day_CMG_Snow_cover' ")
parser.add_option("--ptxy", dest="ptxy", default=None, \
                  help = "one point's x/y or lon/lat pair, in [x,y] format, to extract data \
                          default 'None' ")
parser.add_option("--hdf2nc", dest="hdf2nc", default=False, \
                  help = " convert a hdf file to Netcdf4 file, default: FALSE ", action="store_true")
parser.add_option("--lookup_snowfreeseason", dest="lookup_snowfreeseason", default=False, \
                  help = " lookup snow ending/starting doy in a year and write a Netcdf4 file, default: FALSE ", action="store_true")

(options, args) = parser.parse_args()

#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)
cwdir = options.workdir

if (options.eosh4file == ''):
    print('MUST input file name header by " --eosh4file=??? "')
    sys.exit()
else:
    files = options.eosh4file.split(':')

if (options.vars == ''):
    print('No variable name by " --varname=???" or "ALL/all" ')
    print('Will only do conversion!')
    varnames = ''
elif (options.vars == 'ALL'):
    varnames = 'all'
else:
    varnames = options.vars.split(':')  

# 
if options.ptxy is not None:
    ptxy = re.split(',|:|;| ',options.ptxy)
    ptxy = np.asarray(ptxy,dtype=np.float)
else:
    ptxy = []
##
fileheader = cwdir+'/'+ options.eosh4file
hdftype = options.hdftype

# build all input *.hdf files
alldirfiles = glob.glob("%s*.%s" % (fileheader, hdftype))
if(len(alldirfiles)<=0):
    sys.exit("No hdf file exists - %s*.*.%s in: %s" %(fileheader, hdftype, cwdir))


all_filenames = {} # to be sorted by 'DATE' (No guranttee of DATE is in-order
for hdfname in alldirfiles:
    # readin hdf5 data, and/or if in hdf4 format, convert h4 to h5 before
    vtimename = ['RANGEBEGINNINGDATE']
    vardatas = \
        eoshdf.Read1hdf(hdfname, vtimename, '/usr/local/gcc-x/h4h5tools')
    daynums = date2num(vardatas[vtimename[0]])
    if(len(all_filenames)<1):
        all_filenames['daynums'] = daynums
        all_filenames['fnames']  = hdfname
    else:
        all_filenames['daynums'] = np.vstack((all_filenames['daynums'], daynums))
        all_filenames['fnames']  = np.vstack((all_filenames['fnames'], hdfname))
# sorting data by 'time' (i.e. daynums), so that datasets are in-order by time-series
if len(all_filenames)>0:
    tt = all_filenames['daynums']
    it = sorted(range(len(tt)), key=lambda k: tt[k])
    daynums = tt[it]
    fnames  = all_filenames['fnames'][it]
    print('DONE with Sorting files of total: '+str(len(all_filenames['fnames'])))
del all_filenames
del tt
del it

# convert to geo-referenced NC, with 1 ncfile for each year, without smoothing
if options.hdf2nc and not options.lookup_snowfreeseason:
    year_prv = -9999
    for it in range(len(daynums)):
        # readin hdf5 data, and/or if in hdf4 format, convert h4 to h5 before
        hdfname = fnames[it][0]
        vardatas = \
            eoshdf.Read1hdf(hdfname, varnames, '/usr/local/gcc-x/h4h5tools')

        # if NOT .h5 file, will skip the following so that only h4->h5 done 
        if(len(vardatas)>0 and hdfname.endswith('.h5')):

            # some data processing
            if len(ptxy)>1:
                mid_day, year, doy, lon, lat, snowdata = eoshdf.Prep_modis_snowcov(varnames, vardatas, ptxy=ptxy)
            else:
                mid_day, year, doy, lon, lat, snowdata = eoshdf.Prep_modis_snowcov(varnames, vardatas)
            
            # first, if it's night data in a cell, assign value from previous-day (reasonably)
            if it>0:
                
                for iv in snowdata.keys():
                    ij = np.where( ((snowdata[iv]!=-20) & (snowdata[iv]<0)) 
                                   & (~np.isnan(snowdata_prv[iv])) )  # night (111) ==> undecided  (-10)
                    if len(snowdata[iv][ij])>0: 
                        snowdata[iv][ij] = snowdata_prv[iv][ij]
            snowdata['date'] = mid_day
            snowdata['lon']  = deepcopy(lon)
            snowdata['lat']  = deepcopy(lat)

            # write NC file(s)
            ncfname = hdfname.split('/')[-1] # remove directory name if any
            ncfname = './'+ncfname.split('.')[0]+'_yr'+str(year)+'.nc'
            if os.path.isfile(ncfname) and (it>0 and year==year_prv): 
                if len(ptxy)>1:
                    eoshdf.Write1GeoNc(varnames, snowdata, ptxy=ptxy, ncfname=ncfname, newnc=False)
                else:
                    eoshdf.Write1GeoNc(varnames, snowdata, ptxy=[], ncfname=ncfname, newnc=False)                    
            else:
                if len(ptxy)>1:
                    eoshdf.Write1GeoNc(varnames, snowdata, ptxy=ptxy, ncfname=ncfname, newnc=True)
                else:
                    eoshdf.Write1GeoNc(varnames, snowdata, ptxy=[], ncfname=ncfname, newnc=True)

        
            # save previous data
            snowdata_prv = deepcopy(snowdata)
            year_prv = year
#done with 'hdf2nc'

# convert to geo-referenced NC, with 1 ncfile for each year
# processing 1 variable only but for multiple year-doy, 
# e.g. here for snow-free seasons (DOYs for snow melted fully and snow starts)
if options.lookup_snowfreeseason and len(varnames)==1: 

    # after sorting data, looking-up for DOYs with 'daily snow coverage' < 10% for a year
    crit_snowfree = 25.0
    crit_days = 15  # half of this value as moving-window size to average/max/min snow-coverage  
    snowcov_days = []
    
    snowdata = {}
    snow_yearly = {}

    year_prv   = -9999
    #filetxt = open('d.txt', 'w')
    for it in range(len(daynums)):

        
        # readin hdf5 data, and/or if in hdf4 format, convert h4 to h5 before
        hdfname = fnames[it][0]
        vardatas = \
            eoshdf.Read1hdf(hdfname, varnames, '/usr/local/gcc-x/h4h5tools')

        if(len(vardatas)>0 and hdfname.endswith('.h5')):
            # some data processing
            if len(ptxy)>1:
                mid_day, year, doy, lon, lat, snowcov = eoshdf.Prep_modis_snowcov(varnames, vardatas, ptxy=ptxy)
            else:
                mid_day, year, doy, lon, lat, snowcov = eoshdf.Prep_modis_snowcov(varnames, vardatas)
            snowcov = snowcov[varnames[0]]
            
            #print(year, doy, ': ', snowcov[0][0], file=filetxt)

            # first, if it's missing/night data in a cell, assign value from previous-day (reasonably)
            if it>0:
                ij = np.where( ((snowcov!=-20) & (snowcov<0)) 
                               & (~np.isnan(temp_snowcov_prv)))  # -20 - waterbodies
                               #& (~np.isnan(snowcov_days[crit_days-1,])))  # -20 - waterbodies
                if len(snowcov[ij])>0: 
                    #snowcov[ij] = np.int16(snowcov_days[crit_days-1,][ij])
                    snowcov[ij] = np.int16(temp_snowcov_prv[ij])

            # looks like original daily-data NOT consistent in time, have to do moving-window smoothing
            if it==0:
                snowcov_days = np.full((crit_days,)+snowcov.shape,np.float(-999))
                snowcov_days[0:crit_days,] = snowcov[:,:]
            snowcov_days = np.delete(snowcov_days,0,0)  # remove [0,], i.e. oldest
            snowcov_days = np.vstack((snowcov_days, np.reshape(snowcov,(1,)+snowcov.shape))) # append the newest

            #------------------------------------------------------------------------------------------------
            # smoothing, upon previous day's coverage ('temp_snowcov_prv')
            if it==0:
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
                    temp_snowcov[ij] = snowcov[ij]
                    
            #------------------------------------------------------------------------------------------------
            # write smoothed 'temp_snowcov' to NC
            if options.hdf2nc:
                snowdata['date'] = mid_day
                snowdata['lon']  = deepcopy(lon)
                snowdata['lat']  = deepcopy(lat)
                snowdata[varnames[0]] = deepcopy(temp_snowcov)

                # write NC file(s)
                ncfname = hdfname.split('/')[-1] # remove directory name if any
                ncfname = './'+ncfname.split('.')[0]+'_yr'+str(year)+'_smoothed.nc'
                if os.path.isfile(ncfname) and (it>0 and year==year_prv): 
                    if len(ptxy)>1:
                        eoshdf.Write1GeoNc(varnames, snowdata, ptxy=ptxy, ncfname=ncfname, newnc=False)
                    else:
                        eoshdf.Write1GeoNc(varnames, snowdata, ptxy=[],ncfname=ncfname, newnc=False)                    
                else:
                    if len(ptxy)>1:
                        eoshdf.Write1GeoNc(varnames, snowdata, ptxy=ptxy, ncfname=ncfname, newnc=True)
                    else:
                        eoshdf.Write1GeoNc(varnames, snowdata, ptxy=[], ncfname=ncfname, newnc=True)
                    
                
            #------------------------------------------------------------------------------------------------
             
            if it==0:
                snow_start = np.full(snowcov.shape,np.int16(-999))
                snow_end   = np.full(snowcov.shape,np.int16(-999))
                year_prv   = year
                snow_free  = np.full(snowcov.shape,np.int16(0))
                totdays_yr = 0
           
            if year != year_prv: # save previous-year data when a new year starts, note this will not be for end of 'daynums'
                
                # replace '-20' with those same cells from snowcov (indicating water bodies)
                ij = np.where((snow_end<0) & (temp_snowcov==-20))
                if len(snow_end[ij])>0: snow_end[ij] = temp_snowcov.astype(int16)[ij]
                
                ij = np.where((snow_start<0) & (temp_snowcov==-20))
                if len(snow_start[ij])>1: snow_start[ij] = temp_snowcov.astype(int16)[ij]

                # have to flag never-snow/snowfree cells (snow_free within 'totdays_yr-30d' for a full year)
                if totdays_yr>0:
                    temp_crit = np.int(-30.0*np.float(totdays_yr)/366.0)
                    
                    # never-snow cells
                    temp = snow_free - totdays_yr
                    ij = np.where((temp>=temp_crit) & (snow_end>0))
                    if len(snow_end[ij])>0: snow_end[ij] = -999
                
                    ij = np.where((temp>=temp_crit) & (snow_start>0))
                    if len(snow_start[ij])>0: snow_start[ij] = -999
                
                    # never-snowfree cells
                    ij = np.where((snow_free<-temp_crit)  & (snow_end>0))
                    if len(snow_end[ij])>0: snow_end[ij] = -999
                    ij = np.where((snow_free<-temp_crit)  & (snow_start>0))
                    if len(snow_start[ij])>0: snow_start[ij] = -999
                    ij = np.where((snow_free<-temp_crit)  & (snow_free>0))
                    if len(snow_free[ij])>0: snow_free[ij] = 0
                
                if (len(snow_yearly)<=0):
                    snow_yearly['lon'] = deepcopy(lon)
                    snow_yearly['lat'] = deepcopy(lat)
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
                snow_start[:,:] = -999
                snow_end[:,:] = -999
            

            # snow-melting ends
            if it==0:
                ij = np.where( (temp_snowcov<=crit_snowfree) & (temp_snowcov!=-20) )        #snow-free non-water cells
            else:
                ij = np.where( ((temp_snowcov<=crit_snowfree) & (temp_snowcov!=-20)) &      #snow-free non-water cells
                               ((temp_snowcov_prv>crit_snowfree) & (temp_snowcov_prv!=-20)) )               #snow-covered previously
            if len(snow_end[ij])>0 and doy<213: # snow-melting ends must be not over end of July
                snow_end[ij] = doy
                snow_free[ij]= 1

        
            # snow-covering starts
            if it==0:
                ij = np.where( (temp_snowcov>=crit_snowfree) & (temp_snowcov!=-20) )            #snow-covered cells
            else:
                ij = np.where( ((temp_snowcov>=crit_snowfree) & (temp_snowcov!=-20)) &          #snow-covered cells
                               ((temp_snowcov_prv<crit_snowfree) & (temp_snowcov_prv!=-20)) &   #snow-free previously
                               ((snow_start-snow_end)<=0) )                                     #snow-melting done
            if len(snow_start[ij])>0: 
                snow_start[ij] = doy

            # snow-free days in a year (need this to flag never-snowed cells)
            ij = np.where( (temp_snowcov<crit_snowfree) & (temp_snowcov>=0) &                  # snow-free cells
                           ((snow_start-snow_end)<=0) )                                        # snow-covering not yet
            if len(snow_free[ij])>0: snow_free[ij] = snow_free[ij] + 1
            totdays_yr = totdays_yr + 1

            # checking data
            #if it>0:
            #    print('      ---- ', snowcov[0][0], np.int(temp_snowcov[0][0]), 
            #      np.int(prv1[0][0]), np.int(prv2[0][0]),
            #      snow_end[0][0], snow_start[0][0], snow_free[0][0], file=filetxt)

            if it==len(daynums)-1: # save data at ending of time-steps
                
                # replace '-99' with those same cells from snowcov (indicating water bodies or missings)
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
                    if len(snow_end[ij])>0: snow_end[ij] = -999
                
                    ij = np.where((temp>=temp_crit) & (snow_start>0))
                    if len(snow_start[ij])>0: snow_start[ij] = -999
                
                    # never-snowfree cells
                    ij = np.where((snow_free<-temp_crit)  & (snow_end>0))
                    if len(snow_end[ij])>0: snow_end[ij] = -999
                    ij = np.where((snow_free<-temp_crit)  & (snow_start>0))
                    if len(snow_start[ij])>0: snow_start[ij] = -999
                    ij = np.where((snow_free<-temp_crit)  & (snow_free>0))
                    if len(snow_free[ij])>0: snow_free[ij] = 0
 
                if (len(snow_yearly)<=0):
                    snow_yearly['lon'] = deepcopy(lon)
                    snow_yearly['lat'] = deepcopy(lat)
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
 
            # save previous day for filling night/missing values at next day
            year_prv = year
            temp_snowcov_prv = deepcopy(temp_snowcov)
            
        # end of 1 hdf5 file reading/processing (normally for 1 daynums here)
        print('DONE with data from: '+hdfname)
        
    # end of 'for it in range(daynums) (i.e. all files)

    
    del vardatas
    
    del snow_end
    del snow_start
    del snowcov
    del snowcov_days
    del temp_snowcov
    del temp_snowcov_prv

    #-----------------------------------------------------------
    # write all snow_yearly data to NC file
    ncfname = 'NSIDC_yearly_snowfree.nc'
    if os.path.isfile(ncfname): os.system('rm -rf '+ncfname)
    ncfile = netCDF4.Dataset(ncfname, mode='w',format='NETCDF4') 
    print('Create and Write NC file: '+ncfname)

    lon = snow_yearly['lon']
    lat = snow_yearly['lat']
    lon_dim = ncfile.createDimension('lon',  len(lon))
    lat_dim = ncfile.createDimension('lat',  len(lat))
    yr_dim  = ncfile.createDimension('time', None)

    vlat = ncfile.createVariable('lat', np.float32, ('lat',))
    vlat.units = 'degrees_north'
    vlat.long_name = 'latitude'
    vlat[:] = lat

    vlon = ncfile.createVariable('lon', np.float32, ('lon',))
    vlon.units = 'degrees_east'
    vlon.long_name = 'longitude'
    vlon[:] = lon

    vtime = ncfile.createVariable('time', np.int16, ('time',))
    vtime.units = 'year'
    vtime.long_name = 'year'
    vtime[:] = snow_yearly['Year']

    vdoy_end = ncfile.createVariable('Doy_snowmelted', np.int16, ('time','lat','lon'))
    vdoy_end.units = 'doy'
    vdoy_end.long_name = 'day of year when snow fully melted on ground'
    vdata = snow_yearly['Snow_ending']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdoy_end[:,:,:] = vdata

    vdoy_start = ncfile.createVariable('Doy_snowcovered', np.int16, ('time','lat','lon'))
    vdoy_start.units = 'doy'
    vdoy_start.long_name = 'day of year when snow starts covered on ground'
    vdata = snow_yearly['Snow_starting']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdoy_start[:,:,:] = vdata

    vdays_free = ncfile.createVariable('Days_snowfree', np.int16, ('time','lat','lon'))
    vdays_free.units = 'days'
    vdays_free.long_name = 'days in year from snow ends  to starts covered on ground'
    vdata = snow_yearly['Snow_freedays']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdays_free[:,:,:] = vdata
            
    #global attributes
    ncfile.description = 'ground snow-free starting/ending DOYs from daily snow coverage @0.05 degree resolution'
    ncfile.data_source = ('National Snow and Ice Data Center (NSIDC),' +
                          'MODIS/Terra Snow Cover Daily L3 Global 0.05Deg CMG, Version 6 ')
    ncfile.data_citation = ('Hall, D. K. and G. A. Riggs. 2016.' +
                    'MODIS/Terra Snow Cover Daily L3 Global 0.05Deg CMG, Version 6.' +
                    'Boulder, Colorado USA.'+
                    'NASA National Snow and Ice Data Center Distributed Active Archive Center.'+
                    'doi: https://doi.org/10.5067/MODIS/MOD10C1.006.'+
                    '2020-03-15.')
    ncfile.history = '2020-03-24: Estimated from converted h5, with 10% snow-coverage in a cell as criteria.'
    ncfile.contact = 'F.-M. Yuan, CCSI/ESD-ORNL. yuanf@ornl.gov'

    ncfile.close()
    
    