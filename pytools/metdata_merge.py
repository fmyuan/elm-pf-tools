#!/usr/bin/env python

import sys
import glob
import re
import math
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages
from copy import deepcopy

import Modules_metdata
#from Modules_metdata import clm_metdata_cplbypass_read    # CPL_BYPASS
#from Modules_metdata import clm_metdata_read              # GSWP3
#from Modules_metdata import singleNCDCReadCsvfile         # NCDC 
#from Modules_metdata import subsetDaymetRead1NCfile       # DAYMET in nc4 format
#from Modules_metdata import singleDaymetReadCsvfile       # DAYMET in csv format

from hdf5_modules import Read1hdf                         # for ATS metdata in h5 format

from Modules_plots import SinglePlotting, One2OnePlotting,TimeGridedVarPlotting



# ---------------------------------------------------------------
#  Data sub-timing (down) precipitation (rate-type data), from up-integration (known)
def RateRedistribution(sdata_src, sdata_up, rh=[]):
    # 'sdata_src' is providing sub-timely variation, 
    # while 'sdata_up' is rate integrated-up, which should be conserved after re-distribution
    # For 'precipitation', for example, if no sub-timely pattern (variation, e.g. no rain in a day), 
    #  it's still possible to re-distribute if relative humidity known
    
    if (len(rh) >0):
        if (len(rh) != len(sdata_src)):
            print('Error: redistributing data size different: ', len(sdata_src), len(rh))
            os.exit(-1)
    
    # it's required that len(sdata_src) is muliple of len(sdata_up)
    if(np.mod(len(sdata_src),len(sdata_up))!=0):
        print('Error: length of "sdata_src" not multiple of "sdata_up" - ', len(sdata_src), len(sdata_up))
        os.exit(-2)
    else:
        nsub = int(len(sdata_src)/len(sdata_up))
    
    #
    sdata_adj = []
    for i in range(len(sdata_up)):
        isub = i*nsub
        frac = sdata_src[isub:isub+nsub]
        tot = np.sum(frac)
        if tot!=0: 
            frac = frac/tot
        else:
            frac[:] = 1.0/nsub
            # this is temporarily set, and would be using 'rh' variation (TODO)
        
        sdata_adj = np.hstack((sdata_adj, frac*sdata_up[i]*nsub))  # '*nsub' will allow rate unit same
        
    #
    return sdata_adj

#-------------------------------------------------------------------------------
#  Data sub-timing (down)
def DataTimeDown(sdata_jointting, src, target, data_method='offset'):
    
    if (len(src) != len(target)):
        print('Error: src/target data size different: ', len(src), len(target))
        os.exit(-1)
    
    if (data_method == 'offset'):
        offset = target - src# 
        sdata_jointting = sdata_jointting + offset
    elif (data_method == 'ratio'):
        if (not np.any(src==0.0)):
            ratio = target/src# 
            sdata_jointting = sdata_jointting*ratio
        else:
            print('Error: ratio cannot calculated duo to 0.0 value in data ')
            os.exit(-2)
    else:
        print('Error: not supported method - ',  data_method)
    
    
    #
    return sdata_jointting

# ---------------------------------------------------------------
#  Dataset reduction for temporal patterns
def DataTimePatterns(t, data, SEASONALLY, ANNUALLY, t_unit='Days', ts_yrly=365, detrending=''):
    
    # 'ts_yrly' - timesteps per year in t
    # "t_unit" is for 't', default 'Days' 
    # here 't' unit is 'Days' (i.e. default)
    if t_unit.startswith("H"): 
        t2 = t/24.0
    elif t_unit.startswith("Y"): 
        t2 = t*365.0
    else:
        t2 = deepcopy(t)

    #reshape data in (year,season) shape
    t2=np.asarray(t2)/365
    dim_yr=int(math.ceil(max(t2))-math.floor(min(t2)))
    dim_season = t2.reshape(dim_yr,-1).shape[1]

    t3=np.floor(t2)          # in 'years'
    t2=(t2-np.floor(t2))*365 # still in original time-unit, but in DOY

    
    # if 'season' length in 'data' is greater than pre-defined yrly timesteps
    # need to aggregate first
    if(dim_season>int(ts_yrly)):
        data2d = data.reshape([dim_yr,int(ts_yrly),-1])
        data2d = np.nanmean(data2d,axis=2)
        data = data2d.reshape(dim_yr*int(ts_yrly))    # this will be output

        data2d = t2.reshape([dim_yr,int(ts_yrly),-1])
        data2d = np.nanmean(data2d,axis=2)
        t2 = data2d.reshape([dim_yr*int(ts_yrly)])  # in DOY

        data2d = t3.reshape([dim_yr,int(ts_yrly),-1])
        data2d = np.nanmean(data2d,axis=2)
        t3 = data2d.reshape([dim_yr*int(ts_yrly)])  # in YEARS
        #
        t = t2+t3*365.0                             # this will be output
                
        #reset seasonally length
        dim_season = ts_yrly


    # de-trending, if any
    if (detrending!='yearly' and detrending!=''):
        fit_model = np.polyfit(t, data, 1)
        f_fitted = np.poly1d(fit_model)
        # Note: (1) adjusted 'data_adj' will be used below
        #       (2) 'f_fitted(t[0])' allows adjusted and original 'data' matching at the first data point
        #            i.e. de-trending is starting from 1st point. 
        data_adj = (data - f_fitted(t))+f_fitted(t[0])
    else:
        data_adj = deepcopy(data)
        
    #---------------------------------------------------------
    # ANNUALITY
    t_yrly=[]
    data_yrly = []
    if(ANNUALLY):
        t3=t3.reshape(dim_yr,-1)
        t_yrly=np.mean(t3,axis=1)
        
        if(len(data_adj)>0):
            shp=np.hstack(([-1, dim_season],data_adj.shape[1:]))
            shp=np.asarray(shp,dtype=int)
            data2d=data_adj.reshape(shp)
            data_yrly=np.nanmean(data2d,axis=1)

        # de-trending, if any and NOT YET done above
        if (detrending=='yearly'):
            t_ave = t.reshape(dim_yr,-1)
            t_ave = np.mean(t_ave,axis=1)
            fit_model = np.polyfit(t_ave, data_yrly, 1)
            f_fitted = np.poly1d(fit_model)
            
            # Note: (1) adjusted 'data_adj' will be used below
            #       (2) 'f_fitted(t[0])' allows adjusted and original 'data' matching at the first data point
            #            i.e. de-trending is starting from 1st point. 
            data_adj = (data_adj - f_fitted(t))+f_fitted(t[0])
            
            #
            data2d=data_adj.reshape(shp)
            data_yrly=np.nanmean(data2d,axis=1)

    #---------------------------------------------------------
    # SEASONALITY
    t_seasonlly=[]
    data_seasonally = []
    if(SEASONALLY):
        t2=t2.reshape(dim_yr,-1)
        t_seasonally=np.mean(t2,axis=0)
        if(t_seasonally[0]!=0.0): t_seasonally=np.where(t_seasonally==0,365.0,t_seasonally)  # day 0 shall be day 365 if not first element, otherwise plotting X axis not good
        
        if(len(data)>0):
            shp=np.hstack(([dim_yr,-1],data_adj.shape[1:]))
            shp=np.asarray(shp,dtype=int)
            data2d=data_adj.reshape(shp)
            data_seasonally=np.nanmean(data2d,axis=0)


    #---------------------------------------------------------
    #
    if(detrending!=''):
        return t, data, t_yrly, data_yrly, t_seasonally, data_seasonally, f_fitted
    else:
        return t, data, t_yrly, data_yrly, t_seasonally, data_seasonally

#---------------------------------
#  Dataset 1 & 2 similarity indicators

def DataSimilarity(data1, data2, t=None, t_unit='', v_name='', v_unit='', plt=None):


    
    if len(data1.shape)>1:
        data1=np.reshape(data1,(data1.size))
        data2=np.reshape(data2,(data2.size))
    l11 = [np.nanmin(data1),np.nanmin(data2),np.nanmax(data1),np.nanmax(data2)]
    
    
    # ME/RMSE
    me=np.nanmean(data1-data2)
    rmse=np.nanmean([x*x for x in (data1-data2)])
    rmse=np.sqrt(rmse)

    fit_model = np.polyfit(data1, data2, 1, cov=True)
    slope = fit_model[0]
    intcp = fit_model[1]

    # plotting, if option on
    if plt!=None:
        datalabel = v_name+v_unit
        
        nrow = 2   # sub-plot vertically arranged number (row no.)
        ncol = 1   # sub-plot horizontally arranged number (column no.)
 
        # 1:1 plotting
        isub = 1
        plotlabel = '1:1 PLOT'
        One2OnePlotting(plt, nrow, ncol, isub, data1, data2, datalabel, plotlabel)

        
        # T-series or NONE
        isub = isub + 1
        plotlabel = 'T-series'
        if (t==None):
            t_unit = 'none'
            t = range(0,len(data1)-1)
        
        sdata = np.swapaxes(np.vstack((data1,data2)),0,1)
        varnames = [v_name+'_1',v_name+'_2']
        TimeGridedVarPlotting(plt, nrow, ncol, isub, t, t_unit, sdata, \
                    varnames, plotlabel)


    # return indicators, if any
    return rmse, me

##*****************************************************************************


#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--mettype", dest="met_type", default="CRU", \
                  help="e3sminput met data type (default = 'CRU', options: CRU, cplbypass, GSWP3, GSWP3v1, Site, ATS_h5, NCDC, daymet_nc4)")
# E3SM met data 
parser.add_option("--e3sm_metdir", dest="met_idir", default="./", \
                  help="e3sminput met directory (default = ./, i.e., under current directory)")
parser.add_option("--e3sm_metdomain", dest="met_domain", default="", \
                  help="e3sminput met domain file (default = a file in current met_idir")
parser.add_option("--e3sm_metheader", dest="met_header", default="", \
                  help="e3sminput directory (default = '', i.e., standard file name)")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name to be merged/jointed: ONLY one of 'TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND' ")
parser.add_option("--lon", dest="lon", default=-999, \
                  help = " longitude to be reading/plotting, default -999 for first one")
parser.add_option("--lat", dest="lat", default=-999, \
                  help = " latitude to be reading/plotting, default -999 for first one")
# other sources of data
parser.add_option("--user_metdir", dest="user_metdir", default="./", \
                  help="user-defined met directory (default = ./, i.e., under current directory)")
parser.add_option("--user_metfile", dest="user_metfile", default="", \
                  help="user-defined met file(s) under user_metdir (default = '', i.e., None)")
parser.add_option("--source_adjusting", dest="src_adj", default=False, \
                  help = "adjusting source data by detrended mean", action="store_true")
parser.add_option("--plotting", dest="plotting", default=False, \
                  help = "plotting data for checking", action="store_true")
#
(options, args) = parser.parse_args()


#---------------------------------------------------------------------------

startdays = -9999
enddays = -9999

#-----------------------------------------------------------------------------------
if (options.vars == ''):
    print('No variable name by " --varname=???", one of following ')
    print('TBOT, PRECT, QBOT, RH, FSDS, FLDS, PSRF, WIND')
    sys.exit(-1)
elif(not options.vars in ['TBOT', 'PRECT', 'QBOT', 'RH', 'FSDS', 'FLDS', 'estFLDS', 'PSRF', 'WIND']):
    print('NOT supported variable name, should be one of : ')
    print('TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND')
    sys.exit(-1)


#--------------------------------------------------------------------------------------
# standard met variables for ELM
if(options.vars=='TBOT'):
    vname_elm = 'TBOT'
    varunit = 'K'
elif(options.vars=='PRECT'):
    vname_elm = 'PRECTmms'
    varunit = 'mm/s'
elif(options.vars=='QBOT'):
    vname_elm = 'QBOT'
    varunit = 'kg/kg'
elif(options.vars=='RH'):
    vname_elm = 'RH'
    varunit = '-'
elif(options.vars=='FSDS'):
    vname_elm = 'FSDS'
    varunit = 'W/m2'
elif(options.vars=='FLDS'):
    vname_elm = 'FLDS'
    varunit = 'W/m2'
elif(options.vars=='estFLDS'):
    vname_elm = 'estFLDS'
    varunit = 'W/m2'
elif(options.vars=='PSRF'):
    vname_elm = 'PSRF'
    varunit = 'Pa'
elif(options.vars=='WIND'):
    vname_elm = 'WIND'
    varunit = 'm/s'
else:
    print('NOT supported variable name, should be one of : ')
    print('TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND')
    sys.exit(-1)

t_elm = []
vardatas = []
# read-in metdata from CPL_BYPASS_FULL
if ('cplbypass' in options.met_type or 'CPL' in options.met_type):
    cplbypass_dir=options.met_idir+'/cpl_bypass_full'
    cplbypass_fileheader=''
    cplbypass_mettype='GSWP3_daymet'
    lon = float(options.lon)
    lat = float(options.lat)

    # read in
    ix,iy, varsdims, vardatas = \
        Modules_metdata.clm_metdata_cplbypass_read( \
                            cplbypass_dir,cplbypass_fileheader, cplbypass_mettype, lon, lat, varnames)

    # assign data to plotting variables
    nx=len(ix)
    ny=len(iy)
    if2dgrid = False
    nxy=nx*ny

    tvarname = 'DTIME'    # variable name for time/timing
    #cplvars =  ['DTIME','tunit','t0_datetime','LONGXY','LATIXY',
    #            'FLDS','FSDS','PRECTmms','PSRF','QBOT/RH','TBOT','WIND']

#--------------------------------------------------------------------------------------
# read-in metdata from full met directory, except for CPL_BYPASS
# E3SM met data
elif (('CRU' in options.met_type or \
    'GSWP3' in options.met_type or \
    'Site' in options.met_type) and \
    ('cplbypass' not in options.met_type and 'CPL' not in options.met_type)):
    if (options.met_idir == './'):
        print('e3sm met. directory is the current')
    else:
        print('e3sm met.  directory: '+ options.met_idir)

    metdir=options.met_idir
    metfileheader=options.met_header
    met_type=options.met_type
    met_domain=options.met_domain
    lon = float(options.lon)
    lat = float(options.lat)

    # read in
    ix,iy, varsdims, vardatas = \
        Modules_metdata.clm_metdata_read(metdir,metfileheader, met_type, met_domain, lon, lat,'')
    # assign data to plotting variables
    nx=len(ix)
    ny=len(iy)
    if2dgrid = False
    nxy=nx*ny

    tvarname = 'time'    # variable name for time/timing
    #cplvars =  ['time','tunit','LONGXY','LATIXY',
    #            'FLDS','FSDS','PRECTmms','PSRF','QBOT/RH','TBOT','WIND']

#------------------------------
# if read-in data successfully
if len(vardatas)>0:
    vars_list = list(vardatas.keys())
    
    t = vardatas[tvarname] # days since 1901-01-01-00000
    time_unit = 'days-since-0000-01-01-00:00:00'
    t0 = 1900*365 # converted to days since 0000-01-01-00000

    for varname in vars_list:
        if vname_elm not in varname: continue
        
        if 'PRECTmms' in varname: 
            sdata_elm = deepcopy(vardatas[varname])*86400 # mm/s -> mm/day
        else:
            sdata_elm = deepcopy(vardatas[varname])
        sdata_elm = np.squeeze(sdata_elm)
        #SinglePlotting(t, time_unit, varname, varunit, sdata_elm) # for checking
        
        vname_elm = varname
        break
    #
    t_elm = np.asarray(t)+t0
    del t
    #       
    # usually  either 'RH' or 'QBOT' is available, but if the other is required
    if (vname_elm=='RH'):
        if 'QBOT' in vars_list and 'RH' not in vars_list:
            qbot = np.squeeze(vardatas['QBOT'])
            if 'TBOT' in vars_list:
                tk = np.squeeze(vardatas['TBOT'])
                if 'PSRF' in vars_list:
                    pres_pa = np.squeeze(vardatas['PSRF'])
                else:
                    pres_pa = 101325.0
                sdata_elm = Modules_metdata.converHumidity(tk, pres_pa, q_kgkg=qbot)
            else:
                print('ERROR: for RH coverting from QBOT, air temperature is required')
                sys.exit(-1)
    elif (vname_elm=='QBOT'):
        if 'RH' in vars_list and 'QBOT' not in vars_list:
            rh = np.squeeze(vardatas['RH'])
            if 'TBOT' in vars_list:
                tk = np.squeeze(vardatas['TBOT'])
                if 'PSRF' in vars_list:
                    pres_pa = np.squeeze(vardatas['PSRF'])
                else:
                    pres_pa = 101325.0
            else:
                print('ERROR: for RH coverting from QBOT, air temperature is required')
                sys.exit(-1)
            sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, rh_100=rh)
    # if FLDS estimated from humidity and temperature
    elif (vname_elm=='estFLDS'):
        #Longwave radiation (calculated from air temperature, humidity)
        if 'TBOT' in vars_list:
            tk = np.squeeze(vardatas['TBOT'])
        else:
            print('ERROR: for calculating FLDS, air temperature is required')
            sys.exit(-1)
        if 'PSRF' in vars_list:
            pres_pa = np.squeeze(vardatas['PSRF'])
        else:
            pres_pa = 101325.0
        if 'QBOT' in vars_list:
            qbot = np.squeeze(vardatas['QBOT'])
            rh = []
        elif 'RH' in vars_list:
            rh = np.squeeze(vardatas['RH'])
            qbot = []
        else:
            print('ERROR: for calculating FLDS, either RH or QBOT is required')
            sys.exit(-1)
        sdata_elm = \
            Modules_metdata.calcFLDS(tk, pres_pa, q_kgkg=qbot, rh_100=rh)
    
    # clean-up
    del vardatas, vars_list


#--------------------------------------------------------------------------------------
# read-in metdata from NCDC daily Tmax/Tmin, Precipitation data

if ('NCDC' in options.met_type):
    site,odata_header,odata = \
        Modules_metdata.singleNCDCReadCsvfile('2059560_Alert_initproc.csv','NCDC_metric','-999.99') # 'metric' refers to NCDC data converted to metric system already
        #singleNCDCReadCsvfile('2022211_MysVanKarem_initproc.csv','NCDC_metric','-9999.99') # 'metric' refers to NCDC data converted to metric system already
        #singleNCDCReadCsvfile('GHCND_Alert_initproc.csv','GHCND','-999.99')
    #site_header = ("LATITUDE","LONGITUDE","ELEVATION")
    #data_header = ("YEAR","MONTH","DOY","PRCP","SNWD","TAVG","TMAX","TMIN")
    
    # need to recalculate time
    tvarname = 'time'    # variable name for time/timing
    
    vardatas = {}
    vardatas['time']=(odata['YEAR']-1900.0)*365.0+(odata['DOY']-1.0)  # days since 1900-01-01 00:00:00
    lpyrindex=[i for i, x in enumerate(odata['DOY']) if x == 366.0 ]  # NCDC data is in leap-year calender. Here simply remove last day of the year
    vardatas['time']=np.delete(vardatas['time'],lpyrindex)

    varsdims = {}
    varsdims['time']  = tvarname
    varsdims['TBOTd'] = tvarname
    varsdims['TMAXd'] = tvarname
    varsdims['TMINd'] = tvarname
    varsdims['PRCPd'] = tvarname
    varsdims['SNWDd'] = tvarname

    vardatas['TBOTd']=np.delete(odata['TAVG'],lpyrindex)+273.15
    vardatas['TMAXd']=np.delete(odata['TMAX'],lpyrindex)+273.15
    vardatas['TMINd']=np.delete(odata['TMIN'],lpyrindex)+273.15
    vardatas['PRCPd']=np.delete(odata['PRCP'],lpyrindex)
    vardatas['SNWDd']=np.delete(odata['SNWD'],lpyrindex)

    
    nx=1#len(site['LONGITUDE'])
    ny=1#len(site['LATITUDE'])
    if2dgrid = False
    nxy=nx*ny
    
    ix = 0
    iy = 0
    vars_list = list(varsdata.keys())

    t = vardatas[tvarname] # days since 1900-01-01-00000
    time_unit = 'days-since-0000-01-01-00:00:00'
    t0 = 1900*365 # converted to days since 0000-01-01-00000

    if(options.vars=='TBOT'):
        vname_elm = 'TBOTd'
        varunit = 'K'
    elif(options.vars=='PRECT'):
        vname_elm = 'PRCPd'
        varunit = 'mm/d'
    else:
        print('NOT supported variable name, should be one of : ')
        print('TBOT, PRECT, QBOT, RH, FSDS, FLDS, PSRF, WIND')
        sys.exit(-1)

    for varname in vars_list:
        if vname_ncdc not in varname: continue
        
        vardata = vardatas[varname]
        sdata_ncdc = deepcopy(vardata)
        #SinglePlotting(t, time_unit, varname, varunit, sdata_ncdc) # for checking
        
        vname_ncdc = varname
        break
        
    #
    t_ncdc = np.asarray(t)+t0
    del vardata, vardatas, t
    
#--------------------------------------------------------------------------------------
# read-in metdata from ATS daily forcing

if ('ATS_h5' in options.met_type):
    #hdfname = 'CESM-RCP8.5-2006-2100_dm1985-2015.h5'
    hdfname = options.user_metdir+'/'+options.user_metfile
    varnames=['all']
    vardatas = \
        Read1hdf(hdfname, varnames)

    # 
    t_name='time [s]'
    time_unit='[s]'
    t = vardatas[t_name]
    varnames = [x for x in vardatas.keys() if x!=t_name]
    nvars = len(varnames)
    nt = len(t)
    sdata_ats = np.zeros(nt)

    if(options.vars=='TBOT'):
        vname_ats = 'air temperature'
        varunit = 'K'
    elif(options.vars=='PRECT'):
        vname_ats = 'precipitation [m s^-1]'
        varunit = 'm/s'
    elif(options.vars=='QBOT'):
        vname_ats = 'specific humidity' # note that must do convertion below
        varunit = '-'
    elif(options.vars=='RH'):
        vname_ats = 'relative humidity'
        varunit = '-'
    elif(options.vars=='FSDS'):
        vname_ats = 'incoming shortwave radiation'
        varunit = 'W/m2'
    elif(options.vars=='FLDS' or options.vars=='estFLDS'):
        vname_ats = 'incoming longwave radiation'
        varunit = 'W/m2'
    elif(options.vars=='PSRF'):
        vname_ats = 'NONE'  # a specific situation
        varunit = 'Pa'
    elif(options.vars=='WIND'):
        vname_ats = 'wind speed'
        varunit = 'm/s'
    else:
        print('NOT supported variable name, should be one of : ')
        print('TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND')
        sys.exit(-1)

    t0 = 2015*365 # ==>days since 0000-01-01-000000
    
    for ivar in range(0,nvars):
        if vname_ats not in varnames[ivar]: 
            continue #skip
        else:
            print ('ATS varname: ', varnames[ivar])
        sdata_ats = vardatas[varnames[ivar]]
        varunit = varnames[ivar].strip().split(" ")[-1]
        varname = varnames[ivar].strip().replace(varunit,'')

        if 'rain' in varname or 'snow' in varname or 'precipitation' in varname: 
            sdata_ats = deepcopy(sdata_ats)*86400000 # m/s -> mm/day
        #SinglePlotting(t, time_unit, varname, varunit, sdata_ats) # for checking
        
        if 'relative humidity' in varname:
            sdata_ats = sdata_ats*100.
        
        vname_ats = varnames[ivar]
        break # exit for loop
    
    #
    if(options.vars=='QBOT' or options.vars=='estFLDS'):
        for ivar in range(0,nvars):
            if 'air temperature' in varnames[ivar]: 
                tk = vardatas[varnames[ivar]]
            elif 'relative humidity' in varnames[ivar]:
                rh = vardatas[varnames[ivar]]*100
            else:
                continue
        
        #
        pres_pa = 101325.0
        if(options.vars=='QBOT'):
            sdata_ats = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=[], rh_100=rh)
        #
        if(options.vars=='estFLDS'):
            sdata_ats = Modules_metdata.calcFLDS(tk, pres_pa, q_kgkg=[], rh_100=rh)
    
    #
    t_ats = np.asarray(t)/86400.0+t0
    del vardatas, t
    #


#---------------------------------------------------------------------------------------------------------
# 2 data sets for similiarity and temporal up-scaling/down-scaling

if ('ATS_h5' in options.met_type):
    data1 = sdata_ats
    t1 = t_ats
    vname1=vname_ats
elif ('NCDC' in options.met_type):
    data1 = sdata_ncdc
    t1 = t_ncdc
    vname1=vname_ncdc
if('CRU' in options.met_type or \
   'Site' in options.met_type or \
   'GSWP3' in options.met_type or \
   'CPL' in options.met_type or \
   'cplbypass' in options.met_type): # ELM offline forcing
    data2 = np.squeeze(sdata_elm)
    t2 = t_elm
    vname2=vname_elm

#--------------------------------------------------------------------------------------
# SEASONALITY
SEASONALLY=True
# ANNUALITY
ANNUALLY=True
#
detrending='yearly'
#detrending=''
if (vname1=='NONE' or vname2=='NONE'): detrending=''

if (detrending!=''):
    t1, data1, t1_yrly, data1_yrly, t1_seasonally, data1_seasonally, fitfunc_data1 = \
        DataTimePatterns(t1, data1, SEASONALLY, ANNUALLY, ts_yrly=365, detrending=detrending) 
    # (optional) detrending: '', 'yearly','all', with output 'fitfunc_data1' (np.poly1d, i.e. fitfunc_data1(t1))
    #SinglePlotting(t1, '', ['original','de-trended'], ['-','-'], [sdata_ats, data1]) # for checking 'detrending'
    t2, data2, t2_yrly, data2_yrly, t2_seasonally, data2_seasonally, fitfunc_data2 = \
        DataTimePatterns(t2, data2, SEASONALLY, ANNUALLY, ts_yrly=365, detrending=detrending)
    noted = 'detrended'
else:
    
    t1, data1, t1_yrly, data1_yrly, t1_seasonally, data1_seasonally = \
        DataTimePatterns(t1, data1, SEASONALLY, ANNUALLY, ts_yrly=365)
    # t1/data1 may be integrated if its 'ts_yrly' over 365 (i.e. sub-daily TS)
    #SinglePlotting(t1, '', ['original','de-trended'], ['-','-'], [sdata_ats, data1]) # for checking
    
    t2, data2, t2_yrly, data2_yrly, t2_seasonally, data2_seasonally = \
        DataTimePatterns(t2, data2, SEASONALLY, ANNUALLY, ts_yrly=365)
    noted = 'orig'
    
    

#--------------------------------------------------------------------------------------
# Year Matching

# by comparing statistics/regression fitting
#rmse, se, r2, offset, slope = DataSimilarity(data1, data2, t=t, t_unit=time_unit, v_name='PRECP', v_unit='mm/day', plt=plt)
#plt.show()

#--------------------------------------------------------------------------------------
# Merging

# by jointing 2 datasets directly
BY_JOINTING=True
if (BY_JOINTING):
    noted = noted+'_jointed'

    # Source time/data for up-/-down temporally-scaling
    t_src = deepcopy(t_elm)
    sdata_src = deepcopy(sdata_elm)
    
    # joint-smooth-transition
    if detrending!='':
        # prior to adjust data, save if need to backwardly adjust 'data2' from 'data1'
        if (options.src_adj):

            # precipitation (rate): re-distruting but conserved daily rate
            if (options.vars=='PRECT'):
                for i in range(len(t2_yrly)):
                    iyr2 = t2_yrly[i]
                    iyr1 = t1_yrly[i]  # this is temporary, should be from 'similarity' analysis above
                    idx2 = np.where((t_src>=iyr2*365.0) & (t_src<(iyr2+1.0)*365.0))
                    idx1 = np.where((t1>=iyr1*365.0) & (t1<(iyr1+1.0)*365.0))
                    sdata_iyr = sdata_src[idx2]
                    sdata_iyr = RateRedistribution(sdata_iyr, data1[idx1])
                    sdata_src[idx2] = sdata_iyr
            
            elif (options.vars=='TBOT' or options.vars=='PSRF'):
                # T/P variables
                offset12_yrly = np.nanmean(data1_yrly)-np.nanmean(data2_yrly)
                sdata_src = sdata_src-offset12_yrly
            else:
                # humidity, radiations, wind speedy
                try:
                    ratio12_yrly = np.nanmean(data1_yrly)/np.nanmean(data2_yrly)
                except:
                    # in case np.nanmen(data2_yrly) = 0
                    ratio12_yrly = 1.0
                sdata_src = sdata_src*ratio12_yrly
    
            # data scaled, requiring re-do time-pattern analysis
            t2, data2, t2_yrly, data2_yrly, t2_seasonally, data2_seasonally, fitfunc_data2 = \
                DataTimePatterns(t_src, sdata_src, SEASONALLY, ANNUALLY, ts_yrly=365, detrending=detrending)
            
            # no-need to do scaling below
            jointgap_data = 0.0
            jointgap_yrly = 0.0
            jointgap_seasonally = 0.0
        
        else:
            
            # time-interval sync'ed data2/data1
            tj=(t2[-1]+t1[0])/2.0
            jointgap_data = fitfunc_data2(tj)-fitfunc_data1(tj)
            data2 = data2 - jointgap_data
    
            # integrated data2/data1
            jointgap_yrly = np.nanmean(data2_yrly)-np.nanmean(data1_yrly)
            data2_yrly = data2_yrly - jointgap_yrly
            jointgap_seasonally = np.nanmean(data2_seasonally)-np.nanmean(data1_seasonally)
            data2_seasonally = data2_seasonally - jointgap_seasonally
    
    else:
        jointgap_data = 0.0
        jointgap_yrly = 0.0
        jointgap_seasonally = 0.0
    
    
    # to be jointed/temporally-scaling
    t_jointed = deepcopy(t_src)
    sdata_jointed = deepcopy(sdata_src)
    
    # loop through 't1' by year
    iyr_t_src = 0
    for iyr in np.floor(t1_yrly):
        
        for iday in np.floor(t1_seasonally):
            td0_t1 = iyr*365.0+iday                  # starting point of a day in t1
            td0_t2 = t2_yrly[iyr_t_src]*365.0+iday   # starting point of a day in t2
        
            idx_src = np.where((t_src>=td0_t2) & (t_src<(td0_t2+1.0)))
            
            t_jointting = td0_t1 + (t_src[idx_src]-td0_t2) # the second portion is sub-daily (fraction) time
            t_jointed = np.hstack((t_jointed, t_jointting))
            
            # temporal down-scaling:
            sdata_jointting = sdata_src[idx_src]
            # 'sdata_jointing' - in/out (for sub-time variation), 
            # 'src/target' are time-sync'ed datasets to provide scalors (a moving window)
            twidth = 1.0  # 1-d moving-window (for 'data_method' of 'smooth', t spans 2 or more timesteps)
            idx = np.where((t2>=td0_t2) & (t2<(td0_t2+twidth)))
            src = data2[idx]
            
            if (vname1!='NONE'):
                idx = np.where((t1>=td0_t1) & (t1<(td0_t1+twidth)))
                target = data1[idx]
                
                # 'data_method': 'offset' (default), 'ratio','smooth'
                sdata_jointting = DataTimeDown(sdata_jointting, src, target, data_method='offset')
                
            #else:
                # do nothing if 'data1' not exists, then just simply copy data
            
            sdata_jointed = np.hstack((sdata_jointed, sdata_jointting))

        #DONE 'for iday in np.floor(t2_seasonally)'
        #cycling 'sdata_elm' throughout by year-order
        iyr_t_src = iyr_t_src + 1
        if(iyr_t_src>=len(t2_yrly)): iyr_t_src = 0

    #DONE 'for iyr in np.floor(t1_yrly)'
    if (options.plotting):
        SinglePlotting(t_jointed, 'Jointed', ['-'], ['-'], sdata_jointed)

        txx, dataxx, txx_yrly, dataxx_yrly, txx_seasonally, dataxx_seasonally = \
            DataTimePatterns(t_jointed, sdata_jointed, SEASONALLY, ANNUALLY, ts_yrly=365) 
        SinglePlotting(txx, 'Jointed', ['-'], ['-'], dataxx, figno='XX-1')
        SinglePlotting(txx_yrly, 'Jointed', ['-'], ['-'], dataxx_yrly, figno='XX-2')
        SinglePlotting(txx_seasonally, 'Jointed', ['-'], ['-'], dataxx_seasonally, figno='XX-3')

#--------------------------------------------------------------------------------------
# printing plot in PDF
# for checking data-jointing
if (options.plotting):
    tt=np.hstack((t_elm,t_ats))
    dd=np.hstack((np.squeeze(sdata_elm),sdata_ats))
    SinglePlotting(tt, 'Original', ['-'], ['-'], dd, figno='1-orig')
    
    tt=np.hstack((t2,t1))
    dd=np.hstack((data2,data1))
    SinglePlotting(tt, 'Original-aggregated', ['-'], ['-'], dd, figno='1A-'+noted)
    
    tt=np.hstack((t2_yrly,t1_yrly))
    dd=np.hstack((data2_yrly,data1_yrly))
    SinglePlotting(tt, 'Yearly', ['-'], ['-'], dd, figno='2-'+noted)
    
    tt=np.hstack((t2_seasonally,t1_seasonally))
    dd=np.hstack((data2_seasonally,data1_seasonally))
    SinglePlotting(tt, 'Seasonally-DOY', ['-'], ['-'], dd, figno='3-'+noted)

