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
    for it in range(len(daynums)):
        
        if no_leap:
            year = int(daynums[it]/365.0)
            doy  = daynums[it] - year*365.0 + 1.0
            
            # need this to know if a year to be ending
            year_next = year
            if it<len(daynums)-1: year_next = int(daynums[it+1]/365.0)

        else:
            mid_day = num2date(daynums[it]).date()
            year = mid_day.year
            doy  = daynums[it] - date2num(date(year,1,1)) + 1
            
            year_next = year
            if it<len(daynums)-1: year_next = num2date(daynums[it+1]).date().year

        temp_snowcov = snowcov[it,]
        # NaN/Inf usually cause warning in 'numpy' operation
        temp_snowcov[np.where(np.isnan(temp_snowcov) | np.isinf(temp_snowcov))] = -2 
        
        #------------------------------------------------------------------------------------------------
        # set/reset a few yearly variables in the begining
        if it==0 or doy == 1:
            snow_start = np.full(temp_snowcov.shape,np.int16(-55))
            snow_end   = np.full(temp_snowcov.shape,np.int16(-55))
            snow_free  = np.full(temp_snowcov.shape,np.int16(-55))
            totdays_yr = 0
            
            snowcov_yrx= np.full((3,)+temp_snowcov.shape, np.float32(-2)) # snow-coverage of max in spring, min, and max in fall/winter in a year
        
        
        #------------------------------------------------------------------------------------------------
        # sum of total snow-covered fraction in landed area for every-day in a year
        ij = np.where((temp_snowcov>=0))
        temp_sum = np.asarray([it, doy,np.nan])
        if len(temp_snowcov[ij])>0:
            temp_sum[2]=np.nansum(temp_snowcov[ij])/len(temp_snowcov[ij])
        if it==0 or doy==1:
            snowcov_sum = deepcopy(temp_sum)
        else:
            snowcov_sum = np.vstack((snowcov_sum, temp_sum))
            

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
        # save data at ending of all-timesteps or a yeear
        if (it==len(daynums)-1) or (year_next!=year): 
                
            # replace '-xx' with those same cells from snowcov (i.e. non-land)
            ij = np.where((snow_end<0) & (temp_snowcov<0))
            if len(snow_end[ij])>0: snow_end[ij] = temp_snowcov.astype(int16)[ij]
            
            ij = np.where((snow_start<0) & (temp_snowcov<0))
            if len(snow_start[ij])>1: snow_start[ij] = temp_snowcov.astype(int16)[ij]

            # have to flag never-snow/snowfree cells (snow_free within 'totdays_yr-30d' for a full year)
            if totdays_yr>0:
                # max. in spring, min., and max. in fall/winter
                if snowcov_sum.size>0:
                    # snowcov_sum: [it, doy, sum]
                    imin = np.argmin(snowcov_sum[:,2])
                    imax1= np.argmax(snowcov_sum[0:imin,2])
                    imax2= np.argmax(snowcov_sum[imin:,2]) # note: this starts from 'imin'
                    
                    snowcov_yrx_doy = snowcov_sum[(imax1,imin,imin+imax2),1:3]
                    snowcov_yrx_doy = np.reshape(snowcov_yrx_doy,(snowcov_yrx_doy.size))
                    temp = snowcov[int(snowcov_sum[imax1,0]),]
                    temp[np.where(np.isnan(temp))] = -2
                    #temp[np.where(temp<0)] = -2
                    snowcov_yrx[0,] = temp # this way remain non-land as same value for all (-2 is for sea)
                    temp = snowcov[int(snowcov_sum[imin,0]),]
                    temp[np.where(np.isnan(temp))] = -2
                    #temp[np.where(temp<0)] = -2
                    snowcov_yrx[1,] = temp
                    temp = snowcov[int(snowcov_sum[imin+imax2,0]),]
                    temp[np.where(np.isnan(temp))] = -2
                    #temp[np.where(temp<0)] = -2
                    snowcov_yrx[2,] = temp
                
                
                temp_crit = np.int(-30.0*np.float(totdays_yr)/366.0)
                
                # never-snow cells
                temp = snow_free - totdays_yr # non-snowfree days in a year
                ij = np.where((temp>=temp_crit))
                if len(snow_end[ij])>0: snow_end[ij] = 0
                if len(snow_start[ij])>0: snow_start[ij] = 0
                if len(snow_free[ij])>0: snow_free[ij] = 365
                
                 # never-snowfree cells
                ij = np.where((snow_free<-temp_crit))
                if len(snow_end[ij])>0: snow_end[ij] = 365
                if len(snow_start[ij])>0: snow_start[ij] = 1
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
                if len(snow_end[ij])>0: snow_end[ij] = -11
                if len(snow_start[ij])>0: snow_start[ij] = -11
                if len(snow_free[ij])>0: snow_free[ij] = -11
                 
            # save yearly data
            if (len(snow_yearly)<=0):
                snow_yearly['Year'] = np.empty(1)
                snow_yearly['Snow_ending']   = np.empty((1,)+snow_end.shape)
                snow_yearly['Snow_starting'] = np.empty((1,)+snow_start.shape)
                snow_yearly['Snow_freedays'] = np.empty((1,)+snow_free.shape)

                snow_yearly['Snowcov_xdoy'] = np.empty((1,snowcov_yrx_doy.size))
                snow_yearly['Snowcov_max1'] = np.empty((1,)+temp_snowcov.shape)
                snow_yearly['Snowcov_min']  = np.empty((1,)+temp_snowcov.shape)
                snow_yearly['Snowcov_max2'] = np.empty((1,)+temp_snowcov.shape)

                snow_yearly['Year'][0]           = year
                temp = np.reshape(snow_end,(1,)+snow_end.shape)
                snow_yearly['Snow_ending'][0,]   = deepcopy(temp) # to avoid assign-values by reference
                temp = np.reshape(snow_start,(1,)+snow_start.shape)
                snow_yearly['Snow_starting'][0,] = deepcopy(temp)
                temp = np.reshape(snow_free,(1,)+snow_free.shape)
                snow_yearly['Snow_freedays'][0,] = deepcopy(temp)
                
                snow_yearly['Snowcov_xdoy'][0,] = snowcov_yrx_doy
                snow_yearly['Snowcov_max1'][0,] = snowcov_yrx[0,]
                snow_yearly['Snowcov_min'][0,]  = snowcov_yrx[1,]
                snow_yearly['Snowcov_max2'][0,] = snowcov_yrx[2,]
            else:
                
                snow_yearly['Year'] = np.vstack((snow_yearly['Year'], year))
                temp = np.reshape(snow_end,(1,)+snow_end.shape)
                snow_yearly['Snow_ending'] = np.vstack((snow_yearly['Snow_ending'], deepcopy(temp)))
                temp = np.reshape(snow_start,(1,)+snow_start.shape)
                snow_yearly['Snow_starting'] = np.vstack((snow_yearly['Snow_starting'], deepcopy(temp)))
                temp = np.reshape(snow_free,(1,)+snow_free.shape)
                snow_yearly['Snow_freedays'] = np.vstack((snow_yearly['Snow_freedays'], deepcopy(temp)))
                
                snow_yearly['Snowcov_xdoy'] = np.vstack((snow_yearly['Snowcov_xdoy'], snowcov_yrx_doy))
                temp = np.reshape(snowcov_yrx[0,],(1,)+snow_free.shape)
                snow_yearly['Snowcov_max1'] = np.vstack((snow_yearly['Snowcov_max1'], deepcopy(temp)))
                temp = np.reshape(snowcov_yrx[1,],(1,)+snow_free.shape)
                snow_yearly['Snowcov_min']  = np.vstack((snow_yearly['Snowcov_min'], deepcopy(temp)))
                temp = np.reshape(snowcov_yrx[2,],(1,)+snow_free.shape)
                snow_yearly['Snowcov_max2'] = np.vstack((snow_yearly['Snowcov_max2'], deepcopy(temp)))
            
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
parser.add_option("--ims_varname", dest="ims_varname", default="snowcov", \
                  help = "NIC-IMS output Netcdf file's variable name  to process")
parser.add_option("--elm_varname", dest="elm_varname", default="FSNO", \
                  help = "ELM output Netcdf file's variable name to process")
parser.add_option("--res", dest="res", default="24km", \
                  help = " 1km, 4km or 24km for IMS dataset, default '24km' ")
parser.add_option("--ptslon", dest="ptx", default=None, \
                  help = "point(s)/range longitude, in point or range(:) format, to extract data \
                          default 'None' ")
parser.add_option("--ptslat", dest="pty", default=None, \
                  help = "point(s)/range latitude, in point or range(:) format, to extract data \
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
if options.ptx is not None:
    ptx = re.split(',|:|;| ',options.ptx)
    ptx = np.asarray(ptx,dtype=np.float)
else:
    ptx = []
if options.pty is not None:
    pty = re.split(',|:|;| ',options.pty)
    pty = np.asarray(pty,dtype=np.float)
else:
    pty = []
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
ims_varname = options.ims_varname #'Days_snowfree'#'snowcov'
elm_varname = options.elm_varname #'Days_snowfree'#'FSNO'
# sea cell constant
FillValue_SEA = -2.0  # for 'snow coverage' ranging 0-1.0, -2 is reserved for sea, -1 for sea-ice if any
FillValue_LND = -0.5   # for out-of-map land cell (in IMS original, -99)

if('Days' in ims_varname or 'Days' in elm_varname or \
   'Doy' in ims_varname or 'Doy' in elm_varname): 
    # non-snow-coverage, here for snow-free Days or DOYs (ranging 0-365)
    # for better numerical range, set sea-cell filled value to -55
    FillValue_SEA = -55.0  # -22.0 is for sea-ice cell
    FillValue_LND = -11.0

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
        print('Total Files of ELM outputs: '+str(len(alldirfiles)))

    ncfileheader = elmpathfileheader.split('/')[-1]
    elm_odir   = elmpathfileheader.replace(ncfileheader,'')
    if(elm_odir.strip()==''):elm_odir='./'
    elmfileheader = elmpathfileheader.replace(elm_odir,'')

    elmfincl = 'h0'
    if ('.h1.' in alldirfiles[0]): 
        elmfincl='h1'
    elif ('.h2.' in alldirfiles[0]): 
        elmfincl='h2'
    elif ('.h3.' in alldirfiles[0]): 
        elmfincl='h3'
    elif ('.h4.' in alldirfiles[0]): 
        elmfincl='h4'
    
    # read-in datasets from 1 simulation
    nx, ny, nlgrnd, nldcmp, ncolumn, npft, varsdata, varsdims, ttunits = \
        CLM_NcRead_1simulation(elm_odir, \
                           ncfileheader, \
                           elmfincl, \
                           False, \
                           [elm_varname], \
                           startdays, enddays, \
                           False)

    if2dgrid = True
    if('topo' in varsdata.keys()):
        if(len(varsdata['topo'].shape)==1): if2dgrid = False
    
    # names for variable and its time-axis
    vars_list = list(varsdata.keys())    # var_list is in format of 'h*_varname', in which 'varname' is the real variable names

    # time-invariant variables
    constdata = {}
    for v in vars_list:
        if 'time' not in varsdims[v]: 
            constdata[v]=varsdata[v]
    
    for hv in vars_list:
        if re.search(elm_varname, hv): 
            var_h = hv
            hinc  = hv.replace(elm_varname,'')
            var_t = '%stime' %hinc
            break
    vdata = varsdata[var_h]
    vdims = varsdims[var_h]
    tt = varsdata[var_t]   # time unit: days (default)
    
    # processing original data
    t, gdata, sdata, zdim_indx, pdim_indx = \
        CLMvar_1Dtseries(tt, vdims, vdata, constdata, nx, ny, -999, -999, if2dgrid, \
                    False, False)
    t = np.asarray(t)

    # snow relevant ELM output usually is grid-wised
    if(if2dgrid):
        elm_snowcov = np.reshape(gdata,(len(t),ny,nx))
    else:
        elm_snowcov = gdata
        
    # lon/lat of ELM output
    if ( ('lon' in varsdata.keys() and 'lat' in varsdata.keys()) and 
         ('lon' in varsdims['lon'] and 'lat' in varsdims['lat']) ):
        elmx = varsdata['lon']
        elmy = varsdata['lat']
        elmx_res = 0.50
        elmy_res = 0.50
        elmxy = False
        elmprjstr = ''
    elif ( ('geox' in varsdata.keys() and 'geoy' in varsdata.keys()) and 
           ('geox' in varsdims['geox'] and 'geoy' in varsdims['geoy']) ):
        elmx = varsdata['geox']
        elmy = varsdata['geoy']
        elmx_res = np.nanmean(np.diff(elmx))
        elmy_res = np.nanmean(np.diff(elmy))
        elmxy = True
        # elm proj is lambert conformal conic, when using daymet scaled-forcing
        elmprjstr = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
    else:
        print('Error: could not find lat/lon or geox/geoy dimension')
        sys.exit()
    
    if ('landmask' in varsdata.keys()):
        landmask = varsdata['landmask'] # ELM 2-D data is in [y,x] order
    else:
        landmask = np.full((ny,nx),0.0)
        landmask[np.where(elm_snowcov[0]>=0.0)] = 1.0

    # extract point/range data, if any
    if (len(ptx)==1): # point
        if ptx<0.0: ptx=ptx+360.0
        pti = np.squeeze(np.argmin(abs(elmx-ptx)))
    elif(len(ptx)==2): # range
        ptx[np.where(ptx<0.0)] = ptx[np.where(ptx<0.0)]+360.0
        pti = np.squeeze(np.where((elmx>=min(ptx)) & (elmx<=max(ptx))))
    else:
        pti = np.asarray([])
    if pti.size>0:
        if pti.size==1:
            elmx = elmx[pti,None]
            elm_snowcov = elm_snowcov[:,:,pti,None]
            landmask = landmask[:,pti,None]
        else:
            elmx = elmx[pti]
            elm_snowcov = elm_snowcov[:,:,pti]
            landmask = landmask[:,pti]
    
    if (len(pty)==1): # point
        ptj = np.squeeze(np.argmin(abs(elmy-pty)))
    elif(len(pty)==2): # range
        ptj = np.squeeze(np.where((elmy>=min(pty)) & (elmy<=max(pty))))
    else:
        ptj = np.asarray([])
    if ptj.size>0:
        if ptj.size==1:
            elmy = elmy[None,ptj]
            elm_snowcov = elm_snowcov[:,None,ptj,:]
            landmask = landmask[None,ptj,:]
        else:
            elmy = elmy[ptj]
            elm_snowcov = elm_snowcov[:,ptj,:]
            landmask = landmask[ptj,:]

    elm_seaij = np.where(landmask==0)
    elm_lndij = np.where(landmask==1)

    
    # by default, time is in unit of days since '1850-01-01 00:00:00', without leap-year (SO CANNOT use date/time functions from python)
    if ('days' in ttunits):
        elm_t0 = 1850*365.0
        elm_daynums = t+elm_t0
    elif('year' in ttunits):
        elm_daynums = t*365.0 # first day of year t (i.e. year)
    else:
        elm_daynums = t
    elm_yr=np.trunc(elm_daynums/365.0).astype(int)
    

    # looking-up snow-free seasons over years
    if (options.lookup_snowfreeseason):
        snow_yearly = \
          seeking_YRLYsnowfreeseason(elm_daynums, elm_snowcov, no_leap=True)
        
        # save snow_yearly for later-use
        if (options.imsheader != ""):
            elm_snow_yearly = deepcopy(snow_yearly)
    
    # elm data for use in matching with IMS observation
    elm_vdata = elm_snowcov
    
#-------------------------------------------------------------------------
# NIC-IMS snow coverages reading from daily NC snowcov dataset
if (options.imsheader != ""):
    imspathfileheader = options.imsheader
    if (options.elmheader == ""): elmxy=False
    
    alldirfiles = sorted(glob.glob("%s*.%s" % (imspathfileheader, ftype)))
    if(len(alldirfiles)<=0):
        sys.exit("No file exists - %s*.%s IN %s" %(imspathfileheader, ftype, cwdir))
    else:
        print('Total Files of NIC-IMS datasets: '+str(len(alldirfiles)))

    ncfileheader = imspathfileheader.split('/')[-1]
    ims_odir   = imspathfileheader.replace(ncfileheader,'')
    if(ims_odir.strip()==''):ims_odir='./'
    imsfileheader = imspathfileheader.replace(ims_odir,'')

    snow_yearly = {}
    for ncfile in alldirfiles:
        print ('Processing - ', ncfile)
        f = Dataset(ncfile,'r')
        
        if ncfile == alldirfiles[0]:
            imsx = np.asarray(f.variables['geox'])
            imsy = np.asarray(f.variables['geoy'])
            # Polar stereographic ellipsoidal projection with WGS-84 ellipsoid for IMS data
            #Proj4: +proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs
            imsprjstr = "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"
        
        # read-in datasets from one or multiple nc file 
        if ('daysnum' in f.variables.keys()):
            tt = np.asarray(f.variables['daysnum']) #units - days since date/time 0000-00-00 00:00:00
        elif('time' in f.variables.keys()):
            tt = np.asarray(f.variables['time'])
            if ('year' in f.variables['time'].units):
                tt = [date2num(date(x,1,1)) for x in tt]  # year-01-01
                tt = np.asarray(tt)
        vdata = np.asarray(f.variables[ims_varname])

        # ---------------------------------------------------------------------
        # need to merge ELM snowcov datasets for comparison
        if (options.elmheader !=""):
            
            # ELM centroid lon/lat <==> IMS left/bottom corner's geox/geoy
            if(ncfile==alldirfiles[0]): #only need to do once

                # Construct the IMS grid centroids in lon/lat pairs
                dx = np.diff(imsx)
                xres = np.mean(dx)
                dx = np.hstack((dx,xres)) 
                dy = np.diff(imsy)
                yres = np.mean(dy)
                dy = np.hstack((dy,yres)) 
                # imsx/imsy IS grid centroids, shifting half of intervals for mesh 
                x = imsx - 0.5*dx
                y = imsy - 0.5*dy
                ims_gx, ims_gy = np.meshgrid(x, y)

                # ELM grid-mesh nodes: elmx/y are grid-centroids
                if elmx.size>1:
                    halfx = np.mean(np.diff(elmx))/2.0
                    elmnodex = elmx + 0.50*np.hstack((np.diff(elmx),2.0*halfx))
                    elmnodex = np.hstack((elmx[0]-halfx, elmnodex))
                elif elmx.size == 1:
                    halfx = elmx_res/2.0
                    elmnodex = np.asarray([elmx[0]-halfx, elmx[0]+halfx])
                if elmy.size>1:
                    halfy = np.mean(np.diff(elmy))/2.0
                    elmnodey = elmy + 0.50*np.hstack((np.diff(elmy),2.0*halfy))
                    elmnodey = np.hstack((elmy[0]-halfy, elmnodey))
                elif elmy.size == 1:
                    halfy = elmy_res/2.0
                    elmnodey = np.asarray([elmy[0]-halfy, elmy[0]+halfy])
                
                # matching ims grid-centroids within ELM grid-mesh
                #  note: upon 'elmxy' (True/False), ims_lon/ims_lat may be in geox/y too as elm's
                imsinelm_indx, ims_lon, ims_lat = \
                    nicims.IMS_ELM_gridmatching(elmnodex, elmnodey, \
                                                ims_gx, ims_gy, \
                                                Grid1prjstr=elmprjstr, Grid2prjstr=imsprjstr, \
                                                Grid1_cells=elm_lndij)

                    
                # need ALL time-values from ELM output, but only calculated once
                yr_elm   = np.int32(elm_daynums/365.0)
                doy_elm  = elm_daynums-np.int32(elm_daynums/365.0)*365.0+1.0
            # done with if (ncfile is the first of alldirfiles)
            
            daynums_from_elm = np.empty((0),np.float32)
            elm_it_all = np.empty((0),np.int)
            ims_it_all = np.empty((0),np.int)
            
            # (1) ELM and IMS grid-matching for each time-interval of IMS data
            for it in range(tt[0:].size):
                # match YEAR/DOY btw 'vdata' and 'elm_vdata' (NOT date/time due to no_leap in ELM time)
                day_tt = num2date(tt[it]).date()
                yr_tt  = day_tt.year
                doy_tt = tt[it] - date2num(date(yr_tt,1,1)) + 1
                elm_it = np.squeeze(np.where((yr_elm==yr_tt) & (doy_elm==doy_tt)))
                elm_it_all = np.hstack((elm_it_all,elm_it))
                
                vdata_it = np.float32(vdata[it,])
                vdata_fromelm=np.full(np.shape(vdata_it),np.float32(FillValue_SEA)) # to hold from one-time ELM data (-2/-55 is for sea) in IMS cells
                vdata_fromelm[np.where(vdata_it>=0)] = np.float32(FillValue_LND) # pre-filling IMS land cell with '0'(aka beyond data coverage)

                vdata_1delmcell = np.full((len(elm_lndij[0])),np.float32(FillValue_SEA)) # to hold 1-D aggregated IMS 'vdata' by ELM lnd-grid (-2/-55 is for sea)
                vdata_std_1delmcell = np.full((len(elm_lndij[0])),np.float32(FillValue_SEA)) # to hold 1-D aggregated IMS 'vdata' by ELM lnd-grid (-2/-55 is for sea)
                for idx in range(len(elm_lndij[0])):
                    ij = imsinelm_indx[str(idx)] #  paired-tuple index
                    if (len(ij)<=0):
                        vdata_1delmcell[idx] = FillValue_LND
                        vdata_std_1delmcell[idx] = FillValue_LND
                        continue
                    
                    if(ij[0].size>0):
                        # assign IMS averaged to elm cell
                        i = imsinelm_indx[str(idx)][0]
                        j = imsinelm_indx[str(idx)][1]
                        if (vdata_it[ij].size>0):
                            vdata_1delmcell[idx] = np.nanmean(vdata_it[ij]) # IMS ==> ELM
                            vdata_std_1delmcell[idx] = np.nanstd(vdata_it[ij]) # IMS ==> ELM
                            #vdata_1delmcell[idx] = np.mean(ims_lon[ij]) # IMS ==> ELM, for testing if imsx/y correctly mapped into ELM grid
                            #vdata_std_1delmcell[idx] = np.nanstd(ims_lon[ij])
                            #vdata_1delmcell[idx] = np.mean(imsy[j]) # IMS ==> ELM. for test if imsx/y correctly mapped into ELM grid (note the difference from above line)
                            #vdata_std_1delmcell[idx] = np.nanstd(imsy[j])

                    # assign ELM output to IMS' grids, if time matches
                    if(elm_it.size>0):
                        temp = elm_vdata[elm_it]
                        #if(elm_it.size)>1: temp = np.mean(temp,axis=0)
                        j = elm_lndij[0][idx]  # ELM output data is in (t,elmy,elmx) dimensional-order
                        i = elm_lndij[1][idx]
                        imslnd = np.where(vdata_it[ij]>=0.0)[0] # extract sub-set's IMS land-cell
                        if imslnd.size>0:
                            if ij[0].size<=1:
                                vdata_fromelm[(ij[0],ij[1])] = temp[(j,i)] # ELM ==> IMS (Note: here ELM ji flipped to match with IMS)
                            else:
                                vdata_fromelm[(ij[0][imslnd],ij[1][imslnd])] = temp[(j,i)] # ELM ==> IMS (Note: here ELM ji flipped to match with IMS)
                            #vdata_fromelm[ij] = elmy[j] # ELM ==> IMS, for testing if latitude from ELM mapping to IMS grids
                            #vdata_fromelm[ij] = elmx[i] # ELM ==> IMS, for testing if longitude from ELM mapping to IMS grids
                
                # when all grids are done 
                # reshape vdata_1delmcell according to ELM grids
                temp = np.full((elmy.size,elmx.size),np.float32(FillValue_SEA)) # -2/-55 is for sea
                temp[elm_lndij] = vdata_1delmcell
                temp = np.reshape(temp,(1,)+temp.shape)
                temp2 = np.full((elmy.size,elmx.size),np.float32(FillValue_SEA)) # -2/-55 is for sea
                temp2[elm_lndij] = vdata_std_1delmcell
                temp2 = np.reshape(temp2,(1,)+temp2.shape)
                if(it==0):
                    vdata_2delmcell = np.full((1,elmy.size,elmx.size),np.float32(FillValue_SEA))
                    vdata_2delmcell = deepcopy(temp)
                    vdata_std_2delmcell = np.full((1,elmy.size,elmx.size),np.float32(FillValue_SEA))
                    vdata_std_2delmcell = deepcopy(temp2)
                else:
                    vdata_2delmcell = np.vstack((vdata_2delmcell, temp))
                    vdata_std_2delmcell = np.vstack((vdata_std_2delmcell, temp2))
                    
                # if ELM timing matches with IMS
                if(elm_it.size>0):
                    vdata_fromelm = np.reshape(vdata_fromelm,(1,)+vdata_fromelm.shape)
                    if(len(daynums_from_elm)<=0):
                        snowcov_from_elm = deepcopy(vdata_fromelm)
                    else:
                        snowcov_from_elm = np.vstack((snowcov_from_elm, vdata_fromelm))
                    daynums_from_elm = np.hstack((daynums_from_elm,tt[it]))
                    ims_it_all = np.hstack((ims_it_all,it))
            
            # (2) IMS into ELM grid
            if (vdata_2delmcell.shape[0]>0):  # IMS data aggregated into ELM cells, if sucessfully calculated
                varname = elm_varname+'_ims'
                ncdata = {}

                if (elmxy):
                    ncdata['geox']  = deepcopy(elmx)
                    ncdata['geoy']  = deepcopy(elmy)
                else:
                    ncdata['lon']  = deepcopy(elmx)
                    ncdata['lat']  = deepcopy(elmy)
                
                # the following seems not in right place (TODO)
                if(elm_it_all[0]==elm_it and elm_it>0): # matching-time is NOT the first one
                    # may need to write non-matched ELM data, assuming ELM time is in-order
                    for iy in yr_elm[:elm_it]:
                        for idx in range(len(elm_lndij[0])):
                            ij = imsinelm_indx[str(idx)] #  paired-tuple index

                            # assign ELM output to IMS' grids
                            temp = elm_vdata[elm_it]
                            j = elm_lndij[0][idx]  # ELM output data is in (t,elmy,elmx) dimensional-order
                            i = elm_lndij[1][idx]
                            imslnd = np.where(vdata_it[ij]>=0.0)[0] # extract sub-set's IMS land-cell
                            if imslnd.size>0:
                                if ij[0].size<=1:
                                    vdata_fromelm[(0,ij[0],ij[1])] = temp[(j,i)] # ELM ==> IMS (Note: here ELM ji flipped to match with IMS)
                                else:
                                    vdata_fromelm[(0,ij[0][imslnd],ij[1][imslnd])] = temp[(j,i)] # ELM ==> IMS (Note: here ELM ji flipped to match with IMS)
                        
                
                ncdata['date'] = [num2date(it).date() for it in tt[0:]]
                ncdata[elm_varname] = np.full(vdata_2delmcell.shape,np.float32(FillValue_SEA)) # shape as IMS data, but filled with data from ELM next-line
                ncdata[elm_varname][ims_it_all,] = deepcopy(elm_vdata[elm_it_all])
                ncdata[varname] = deepcopy(vdata_2delmcell)
                ncdata[varname+'_std'] = deepcopy(vdata_std_2delmcell)

                # difference btw ELM and IMS
                temp = np.full(vdata_2delmcell.shape,0.0) # shape as IMS data, but will be filled with data from ELM next-line
                
                temp = ncdata[elm_varname] - ncdata[varname]
                # since 'diff' could be negative, better to mark non-land in either dataset as ZERO
                ij = np.where( (ncdata[elm_varname]<0)  | \
                               (np.isnan(ncdata[elm_varname])) | (np.isinf(ncdata[elm_varname])) )
                temp[ij] = 0.0 
                ij = np.where( (ncdata[varname]<0)  | \
                               (np.isnan(ncdata[varname])) | (np.isinf(ncdata[varname])) )
                temp[ij] = 0.0
                # for best visual effect, mask sea-cell  of either dataset
                ij = np.where( (ncdata[elm_varname]==FillValue_SEA)  | \
                               (ncdata[varname]==FillValue_SEA) )
                temp[ij] = FillValue_SEA
                ij = np.where( (ncdata[elm_varname]==FillValue_LND)  | \
                               (ncdata[varname]==FillValue_LND) )
                temp[ij] = FillValue_LND
                ncdata[elm_varname+'_diff'] = temp

                ncfname = ncfile.split('/')[-1] # remove directory name if any
                ncfname = './ELM_obs_from_'+ncfname
                if (elmxy):
                    prjname = 'lambert conformal conic projection: '+ \
                    '"+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"'
                else:
                    prjname = ''
                    
                nicims.Write1GeoNc([elm_varname, elm_varname+'_diff', varname,varname+'_std'], \
                                   ncdata, ptxy=[], ncfname=ncfname, newnc=True, \
                                   FillValue=FillValue_SEA, geoprj=prjname)
                
                
            # (3) ELM into IMS
            if (len(daynums_from_elm)>0): # ELM simulation assigned into IMS grids
                varname = ims_varname+'_elm'
                ncdata = {}
                ncdata['date'] = [num2date(it).date() for it in daynums_from_elm]
                # truncating too-much beyond-ELM grids (i.e. lower-latitudes)
                if (elmxy):
                    ij=np.where( (ims_lat>=np.min(elmnodey)) &
                                 (ims_lat<=np.max(elmnodey)) &
                                 (ims_lon>=np.min(elmnodex)) &
                                 (ims_lon<=np.max(elmnodex)) )
                else:
                    ij=np.where(ims_lat>=(np.min(elmnodey)-5.0))  # since it's polar stereo, no max limit 
                imax=np.max(ij[1])
                imin=np.min(ij[1])
                jmax=np.max(ij[0])
                jmin=np.min(ij[0])
                ncdata['geox']  = deepcopy(imsx[imin:imax+1])
                ncdata['geoy']  = deepcopy(imsy[jmin:jmax+1])
                vdata_ij = vdata[ims_it_all,]
                # the following 2 lines ONLY good for snowcov from NIC-IMS dataset
                vdata_ij[np.where(vdata_ij==-1)] = -2.0 # originally -1 is sea-ice, which not available from ELM
                vdata_ij[np.where(vdata_ij==-99)] = -0.5 # originally -99 is land out-of-map coverage, which not available from ELM
                
                #
                imslndmask = deepcopy(vdata_ij[0,]) # for IMS land-cells, only need once
                ijout = np.where( ((ims_lon<np.min(elmnodex)) | (ims_lon>np.max(elmnodex)) |
                                  (ims_lat<np.min(elmnodey)) | (ims_lat>np.max(elmnodey))) &
                                  (imslndmask>=0.0) ) # land-cells, only need once
                iout = np.asarray(ijout[0]) #ijout is paired indices in tuple
                jout = np.asarray(ijout[1])
                vdata_ij[:,iout,jout] = np.float32(FillValue_LND)
                
                vdata_ij = deepcopy(vdata_ij[:,jmin:jmax+1, imin:imax+1]) # IMS-data in [t,i,j] order of dimension
                ncdata[ims_varname] = deepcopy(vdata_ij) 
                vdata_elm_ij = snowcov_from_elm[:,jmin:jmax+1, imin:imax+1] # ELM-data in [t,j,i] order of dimension, but flip already
                ncdata[varname] = deepcopy(vdata_elm_ij) 

                # difference btw ELM and IMS
                temp = np.full(ncdata[ims_varname].shape,np.float32(FillValue_SEA)) # shape as IMS data, but filled with data from ELM next-line
                temp = ncdata[varname] - ncdata[ims_varname] # ELM - IMS
                # since 'diff' could be negative, better to mark non-land in either dataset as ZERO
                ij = np.where( (ncdata[ims_varname]<0) | \
                              (np.isnan(ncdata[ims_varname])) | (np.isinf(ncdata[ims_varname])) )
                temp[ij] = 0.0 
                ij = np.where( (ncdata[varname]<0) | \
                               (np.isnan(ncdata[varname])) | (np.isinf(ncdata[varname])) )
                temp[ij] = 0.0 
                # for best visual effect, mask sea-cell  of either dataset
                ij = np.where( (ncdata[ims_varname]==FillValue_SEA)  | \
                               (ncdata[varname]==FillValue_SEA) )
                temp[ij] = FillValue_SEA
                ij = np.where( (ncdata[ims_varname]==FillValue_LND)  | \
                               (ncdata[varname]==FillValue_LND) )
                temp[ij] = FillValue_LND
                ncdata[ims_varname+'_diff'] = temp

                # write NC file(s)
                ncfname = ncfile.split('/')[-1] # remove directory name if any
                ncfname = './ELM_sim_for_'+ncfname
                prjname = 'Polar stereographic ellipsoidal projection: '+ \
                    '"+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs" '
                nicims.Write1GeoNc([ims_varname, ims_varname+'_diff', varname], \
                                   ncdata, ptxy=[], ncfname=ncfname, newnc=True, \
                                   FillValue=FillValue_SEA,geoprj=prjname)
        
        # end of if (option.elmheader != '') # i.e. merging two datasets of ELM and IMS
        # but only with 1 IMS nc file.
        # -----------------------------------------------------------------------------
        
        # at end of 1 nc file, need to concat it into 1 all-time datasets for yearly lookup snow statistics
        if ncfile == alldirfiles[0]:
            daynums = deepcopy(tt)
            snowcov = deepcopy(vdata)
        else:
            daynums = np.hstack((daynums, tt))
            snowcov = np.vstack((snowcov, vdata))
    
        # ONLY needs at most 2 year data for 'lookup snowfree season' really - this will save a lot of memory if more than 2-year data
        # (NOTE: here assuming data time-series is in order, i.e. continuously)
        if (options.lookup_snowfreeseason):
            daynums_yr = np.asarray([num2date(x).year for x in daynums])
            if(daynums_yr[-1]>daynums_yr[0] or ncfile==alldirfiles[-1]): # year number changed or last file
                if (ncfile==alldirfiles[-1]):
                    ij_1yr = np.where(daynums_yr!=np.nan)[0] # all T-series to end, likely 1 or 2-years
                else:
                    ij_1yr = np.where(daynums_yr==daynums_yr[0])[0] # this is 1-D tuple index, better to get its [0]
                print('Snow Statistics on Data of T-length: '+str(len(ij_1yr))+ \
                      ' From: '+ str(num2date(daynums[ij_1yr[0]]))+ \
                      ' To: ' + str(num2date(daynums[ij_1yr[-1]])) )
                if (len(snow_yearly)<=0):
                    snow_yearly = \
                        seeking_YRLYsnowfreeseason(daynums[ij_1yr], snowcov[ij_1yr], no_leap=False)
                else:
                    snow_yearly = \
                        seeking_YRLYsnowfreeseason(daynums[ij_1yr], snowcov[ij_1yr], no_leap=False,
                                                   snow_yearly=snow_yearly)
                
                
                # only need the rest of data for next year
                daynums = deepcopy(daynums[ij_1yr[-1]+1:])
                snowcov = deepcopy(snowcov[ij_1yr[-1]+1:,])
    
    # done with all nc files

    # looking-up snow-free seasons over years (THIS requires huge amount memory)
    #if (options.lookup_snowfreeseason):
    #    snow_yearly = \
    #      seeking_YRLYsnowfreeseason(daynums, snowcov)


#------------------------------------------------------------------------------------------------
    

# write yearly snow-free season data into 1 NC file
if len(snow_yearly)>0:
    yr0 = int(snow_yearly['Year'][0])
    yr1 = int(snow_yearly['Year'][-1])
    
    #-----------------------------------------------------------
    # write all snow_yearly data to NC file
    if(options.imsheader != "" and options.elmheader==""):
        ncfname = 'Yly_snowstats-'+imsfileheader+str(yr0)+'-'+str(yr1)+'.nc'

        snow_yearly['geox'] = deepcopy(imsx)
        snow_yearly['geoy'] = deepcopy(imsy)
        elmxy = False
    elif(options.elmheader != "" and options.imsheader==""):
        ncfname = 'Yly_snowstats-'+elmfileheader+'.'+elmfincl+'.'+str(yr0)+'-'+str(yr1)+'.nc'
        if(elmxy):
            snow_yearly['geox']  = deepcopy(elmx)
            snow_yearly['geoy']  = deepcopy(elmy)
        
        else:
            snow_yearly['lon']  = deepcopy(elmx)
            snow_yearly['lat']  = deepcopy(elmy)

    if os.path.isfile(ncfname): os.system('rm -rf '+ncfname)
    ncfile = Dataset(ncfname, mode='w',format='NETCDF4') 
    print('Create and Write NC file: '+ncfname)

    yr_dim  = ncfile.createDimension('time', None)
    if 'lon' in snow_yearly.keys():
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
        vlat = ncfile.createVariable('geoy', np.float32, ('geoy',))
        vlat.units = 'meters'
        if (elmxy):
            vlat.long_name = 'Northing (lambert conformal conic projection)'
            
            vproj = ncfile.createVariable('lambert_conformal_conic', np.int32)
            vproj.grid_mapping_name = "lambert_conformal_conic"
            vproj.longitude_of_central_meridian = -100.
            vproj.latitude_of_projection_origin = 42.5
            vproj.false_easting = 0.
            vproj.false_northing = 0.
            vproj.standard_parallel = 25., 60.
            vproj.semi_major_axis = 6378137.
            vproj.inverse_flattening = 298.257223563 

        else:
            vlat.long_name = 'Northing (Polar stereographic ellipsoidal projection)'
        vlon = ncfile.createVariable('geox', np.float32, ('geox',))
        vlon.units = 'meters'
        if (elmxy):
            vlon.long_name = 'Easting (lambert conformal conic projection)'
        else:
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
    vdoy_end.Key = "0-365 = land with/without snow, -22 = sea with ice (if any), -55 = sea, -11 = beyond map coverage (if any)" ;
    vdata = snow_yearly['Snow_ending']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdoy_end[:,:,:] = vdata

    if 'lon' in snow_yearly.keys():
        vdoy_start = ncfile.createVariable('Doy_snowcovered', np.int16, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vdoy_start = ncfile.createVariable('Doy_snowcovered', np.int16, ('time','geoy','geox'))
    vdoy_start.units = 'doy'
    vdoy_start.long_name = 'day of year when snow starts covered on ground'
    vdoy_start.Key = "0-365 = land with/without snow, -22 = sea with ice (if any), -55 = sea, -11 = beyond map coverage (if any)" ;
    vdata = snow_yearly['Snow_starting']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdoy_start[:,:,:] = vdata

    if 'lon' in snow_yearly.keys():
        vdays_free = ncfile.createVariable('Days_snowfree', np.int16, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vdays_free = ncfile.createVariable('Days_snowfree', np.int16, ('time','geoy','geox'))
    vdays_free.units = 'days'
    vdays_free.long_name = 'days in year from snow ends  to starts covered on ground'
    vdays_free.Key = "0-365 = land with/without snow, -22 = sea with ice (if any), -55 = sea, -11 = beyond map coverage (if any)" ;
    vdata = snow_yearly['Snow_freedays']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vdays_free[:,:,:] = vdata

    #--------------------------------------------
    if 'lon' in snow_yearly.keys():
        vmax1 = ncfile.createVariable('Snowcov_Yrlymax1', np.float32, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vmax1 = ncfile.createVariable('Snowcov_Yrlymax1', np.float32, ('time','geoy','geox'))
    vmax1.units = ' '
    vmax1.long_name = 'Snow coverage when regionally-summed at max. in spring '
    vmax1.Key = "0 = land without snow, 0~1.0 = land with snow, -1 = sea with ice, -2 = sea, -0.5 = beyond map coverage (if remarked)" ;
    vdata = snow_yearly['Snowcov_max1']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vmax1[:,:,:] = vdata

    if 'lon' in snow_yearly.keys():
        vmin = ncfile.createVariable('Snowcov_Yrlymin', np.float32, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vmin = ncfile.createVariable('Snowcov_Yrlymin', np.float32, ('time','geoy','geox'))
    vmin.units = ' '
    vmin.long_name = 'Snow coverage when regionally-summed at min. in a year '
    vmin.Key = "0 = land without snow, 0~1.0 = land with snow, -1 = sea with ice, -2 = sea, -0.5 = beyond map land coverage (if remarked)" ;
    vdata = snow_yearly['Snowcov_min']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vmin[:,:,:] = vdata

    if 'lon' in snow_yearly.keys():
        vmax2 = ncfile.createVariable('Snowcov_Yrlymax2', np.float32, ('time','lat','lon'))
    if 'geox' in snow_yearly.keys():
        vmax2 = ncfile.createVariable('Snowcov_Yrlymax2', np.float32, ('time','geoy','geox'))
    vmax2.units = ' '
    vmax2.long_name = 'Snow coverage when regionally-summed at max. in Fall/Winter '
    vmax2.Key = "0 = land without snow, 0~1.0 = land with snow, -1 = sea with ice, -2 = sea, -0.5 = beyond map land-coverage (if remarked)" ;
    vdata = snow_yearly['Snowcov_max2']
    if len(vdata.shape)<3: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vmax2[:,:,:] = vdata

    vdata = snow_yearly['Snowcov_xdoy']
    nxval_dim = ncfile.createDimension('nxvals', vdata.shape[1])
    vxdoy = ncfile.createVariable('Snowcov_xdoy', np.float32, ('time','nxvals'))
    vxdoy.units = 'doy/fraction'
    vxdoy.long_name = 'DOY/Fraction when Snow coverage regionally-summed at max. in Spring, min., and in late-Fall/Winter '
    if len(vdata.shape)<2: vdata = np.reshape(vdata,(1,)+vdata.shape) # in case that data only 1 time-step
    vxdoy[:,:] = vdata

            
    #
    #global attributes
    if(options.elmheader !=""):
        ncfile.description = 'ELM simulated No60 and above Northern High-latitude Region yearly snow-cover gone/start day and snowfree length @0.5deg resolution '
        if (elmxy):
            ncfile.data_source = ('Offline E3SM Land Model, ' +
                'forced by GSWP3 v2/DaymetNA, ~1km, over Northern America, fully CNP coupled' +
                'Master @ 2018-11-01 ')
            ncfile.history = '2020-12-21: calculated from  daily snow coverage FSNO simulations.'
        else:
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
    
    