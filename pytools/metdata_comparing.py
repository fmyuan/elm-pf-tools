#!/usr/bin/env python

import os, sys
from datetime import datetime
import glob
import re
import math
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages
from copy import deepcopy

# customized modules
import Modules_metdata
#from Modules_metdata import clm_metdata_cplbypass_read    # CPL_BYPASS
#from Modules_metdata import clm_metdata_read              # GSWP3
#from Modules_metdata import singleNCDCReadCsvfile         # NCDC 
#from Modules_metdata import subsetDaymetRead1NCfile       # DAYMET in nc4 format
#from Modules_metdata import singleDaymetReadCsvfile       # DAYMET in csv format
from hdf5_modules import Read1hdf                         # for ATS metdata in h5 format
from Modules_plots import SubPlotting, SinglePlotting, One2OnePlotting,TimeGridedVarPlotting
import netcdf_modules as nfmod

# ---------------------------------------------------------------
# read in standard ELM met data
def ELMmet_1vardatas(vars, met_dir, met_type, met_fileheader='', met_domain='',lon=-999,lat=-999):
    # standard met variables for ELM
    if(vars=='TBOT'):
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
    if ('cplbypass' in met_type or 'cplbypass_site' in met_type):
        cplbypass_dir=met_dir
        cplbypass_mettype='GSWP3'
        if ('daymet' in met_type): cplbypass_mettype='GSWP3_daymet'
        if ('v1' in met_type): cplbypass_mettype='GSWP3v1'
        
        vnames=[vname_elm]
        if('cplbypass_Site' in met_type): 
            cplbypass_mettype='Site'
            if (vname_elm=='QBOT' or vname_elm=='estFLDS'): 
                vnames=['RH','TBOT','PSRF']
        else:
            if (vname_elm=='RH' or vname_elm=='estFLDS'): 
                vnames=['QBOT','TBOT','PSRF']
    

        # read in
        ix,iy, varsdims, vardatas = \
            Modules_metdata.clm_metdata_cplbypass_read( \
                            cplbypass_dir,cplbypass_mettype, lon, lat, vnames)

        # assign data to plotting variables
        nx=len(ix)
        ny=len(iy)
        if2dgrid = False
        nxy=nx*ny

        tvarname = 'DTIME'    # variable name for time/timing
        #cplvars =  ['DTIME','tunit','LONGXY','LATIXY',
        #            'FLDS','FSDS','PRECTmms','PSRF','QBOT/RH','TBOT','WIND']
        LATIXY = vardatas['LATIXY']
        LONGXY = vardatas['LONGXY']

    #--------------------------------------------------------------------------------------
    # read-in metdata from full met directory, except for CPL_BYPASS
    # E3SM met data
    elif (('CRU' in met_type or \
           'GSWP3' in met_type or \
           'Site' in met_type) and \
           ('cplbypass' not in met_type and 'cplbypass_Site' not in met_type)):
        if (met_dir == './'):
            print('e3sm met. directory is the current')
        else:
            print('e3sm met.  directory: '+ met_dir)

        # read in
        ix,iy, varsdims, vardatas = \
            Modules_metdata.clm_metdata_read(met_dir, met_fileheader, met_type, met_domain, lon, lat,'')
        # assign data to plotting variables
        nx=len(ix)
        ny=len(iy)
        if2dgrid = False
        nxy=nx*ny

        tvarname = 'time'    # variable name for time/timing
        #vars =  ['time','tunit','LONGXY','LATIXY',
        #            'FLDS','FSDS','PRECTmms','PSRF','QBOT/RH','TBOT','WIND']
        LATIXY = vardatas['LATIXY']
        LONGXY = vardatas['LONGXY']

    #------------------------------
    # if read-in ELM forcing data successfully
    if len(vardatas)>0:
        vars_list = list(vardatas.keys())
        
        t = vardatas[tvarname] # days since 1901-01-01 00:00
        tunit = vardatas['tunit']
        tunit = tunit.replace('1901-','0001-')
        t0 = 1901*365 # converted from 'days-since-1901-01-01-00:00' to days since 0001-01-01 00:00

        for varname in vars_list:
            if vname_elm not in varname: continue
        
            if 'PRECTmms' in varname: 
                sdata_elm = deepcopy(vardatas[varname])*86400 # mm/s -> mm/day
            else:
                sdata_elm = deepcopy(vardatas[varname])
            sdata_elm = np.squeeze(sdata_elm)
            #SubPlotting(t, time_unit, varname, varunit, sdata_elm) # for checking
        
            vname_elm = varname
            break
        #
        t_elm = np.asarray(t)+t0
        tunit_elm = tunit
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
                    sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=qbot)
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
            vname_elm='FLDS'
        
        # clean-up
        del vardatas, vars_list
        
        # in case when needs cut-off (TO comment out, change the 'if True' to 'if False'
        if False:
            idx=np.where(t_elm>=1985*365.0)
            t_elm = t_elm[idx]
            sdata_elm = sdata_elm[idx]
        # aggregate half-hourly to 3-hourly (e.g. to GSWP3 datasets)
        if False:
            t_elm = np.reshape(t_elm, [-1,6])
            t_elm = np.squeeze(t_elm[:,0])
            sdata_elm = np.reshape(sdata_elm, [-1,6])
            sdata_elm = np.mean(sdata_elm,axis=1)
            
    #
    #
    return vname_elm, t_elm, tunit_elm, sdata_elm, LONGXY, LATIXY

# ---------------------------------------------------------------
# ---------------------------------------------------------------
#  Dataset time aggregation
def DataTimeAggregation(t, data, t_unit='Days',DWMY='daily'):
        
    # 'DWMY' - aggregation option: daily, weekly, monthly, or, yearly
    # "t_unit" is for 't', default 'Days' 
    # here 't' unit is 'Days' (i.e. default)
    if t_unit.startswith("H"): 
        t2 = t/24.0
    elif t_unit.startswith("Y"): 
        t2 = t*365.0
    else:
        t2 = deepcopy(t)

    #
    t2=np.asarray(t2)/365
    t3=np.floor(t2)               # in 'years'
    t2=(t2-np.floor(t2))*365      # still in original time-unit, but in DOY only (zero-based here)
    # doy 0, if from above, implies it be 365 and years be previous year, when it's not the first.
    # when it's the starting, it's OK, otherwise it may cause issue for regrouping and plotting
    if (t2[0]!=0):  
        idx = np.where(t2==0)
        t3[idx] = t3[idx]-1
        t2[idx] = 365-1.0e-8
    if(abs(t2[1]-t2[0])<1.0):
        t2 = t2 + 1                  # convert to 1-based DOY, if sub-daily data

    
    if (DWMY=='daily' or DWMY =='d' or DWMY=='D'):
        # doy range in format of as [1,32) for Jan. (1-based)
        mm=deepcopy(t2)
        for m in range(366):
            mm[np.where((t2>=m) & (t2<(m+1)))]=m
        #
        yrs = t3+mm/365.0
    elif (DWMY=='weekly' or DWMY =='w' or DWMY=='W'):
        # doy range in format of as [1,32) for Jan. (1-based)
        mm=deepcopy(t2)
        for m in range(51):
            mm[np.where((t2>=7*m) & (t2<7*(m+1)))]=m
        mm[np.where((t2>=7*(m+1)) & (t2<366))]=m+1  # day 365 as in week 52
        #
        yrs = t3+mm/52.0
    elif(DWMY=='monthly' or DWMY =='m' or DWMY=='M'):
        # doy range in format of as [1,32) for Jan. (1-based)
        doym = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
        mm=deepcopy(t2)
        for m in range(12):
            mm[np.where((t2>=doym[m]) & (t2<doym[m+1]))]=m
        #
        yrs = t3+mm/12.0
    
    #
    yrs = yrs*365.0  # yrs --> days
    tt = np.unique(yrs)
    sdata = np.full(np.shape(tt), np.nan, dtype=data.dtype)
    for i in range(len(tt)):
        idx = np.where(yrs == tt[i])
        sdata[i]=np.nanmean(data[idx])

    #
    return tt, sdata

# ---------------------------------------------------------------
#  Dataset reduction for temporal patterns
def DataTimePatterns(t, data, SEASONALLY, ANNUALLY, t_unit='Days', ts_yrly=365, yrmin=-9999, yrmax=-9999):
    
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
    # doy 0, if from above, implies it be 365 and years be previous year, when it's not the first.
    # when it's the starting, it's OK, otherwise it may cause issue for regrouping and plotting
    if (t2[0]!=0):  
        idx = np.where(t2==0)
        t3[idx] = t3[idx]-1
        t2[idx] = 365-1.0e-8
    if(abs(t2[1]-t2[0])<1.0):
        t2 = t2 + 1                  # convert to 1-based DOY, if sub-daily data

    
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
    
    #---------------------------------------------------------
    # ANNUALITY
    t_yrly=[]
    data_yrly = []
    if(ANNUALLY):
        t3=t3.reshape(dim_yr,-1)
        t_yrly=np.mean(t3,axis=1)
        
        if(len(data)>0):
            shp=np.hstack(([-1, dim_season],data.shape[1:]))
            shp=np.asarray(shp,dtype=int)
            data2d=data.reshape(shp)
            data_yrly=np.nanmean(data2d,axis=1)

    #---------------------------------------------------------
    # SEASONALITY
    t_seasonlly=[]
    data_seasonally = []
    if(SEASONALLY):
        t2=t2.reshape(dim_yr,-1)
        # year range, if any (NOTE: didn't do so for annually above)
        if (yrmin!=-9999 or yrmax!=-9999):
            if (not ANNUALLY): t3=t3.reshape(dim_yr,-1)
            if (yrmin!=-9999 and yrmax!=-9999): 
                t2=t2[np.where((t3[:,0]>=yrmin) & (t3[:,0]<=yrmax))]
            elif (yrmin!=-9999): 
                t2=t2[np.where(t3[:,0]>=yrmin)]
            elif (yrmax!=-9999): 
                t2=t2[np.where(t3[:,0]<=yrmax)]
            
        t_seasonally=np.mean(t2,axis=0)
        
        if(len(data)>0):
            shp=np.hstack(([dim_yr,-1],data.shape[1:]))
            shp=np.asarray(shp,dtype=int)
            data2d=data.reshape(shp)
            # year range, if any (NOTE: didn't do so for annually above)
            if (yrmin!=-9999 or yrmax!=-9999):
                if (not ANNUALLY): t3=t3.reshape(dim_yr,-1)
                if (yrmin!=-9999 and yrmax!=-9999): 
                    data2d=data2d[np.where((t3[:,0]>=yrmin) & (t3[:,0]<=yrmax))]
                elif (yrmin!=-9999): 
                    data2d=data2d[np.where(t3[:,0]>=yrmin)]
                elif (yrmax!=-9999): 
                    data2d=data2d[np.where(t3[:,0]<=yrmax)]

            data_seasonally=np.nanmean(data2d,axis=0)


    #---------------------------------------------------------
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

# E3SM met data 
parser.add_option("--mettype", dest="met_type", default="CRU", \
                  help="e3sminput met data type (default = 'CRU', options: CRU, GSWP3, GSWP3v1, Site, cplbypass, cplbypass_Site")
parser.add_option("--metdir", dest="met_dir", default="./", \
                  help="e3sminput met directory (default = ./, i.e., under current directory)")
parser.add_option("--metdomain", dest="met_domain", default="", \
                  help="e3sminput met domain file (default = a file in current met_idir")
parser.add_option("--metheader", dest="met_header", default="", \
                  help="e3sminput directory (default = '', i.e., standard file name)")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name to be merged/jointed: ONLY one of 'TBOT, PRECT, QBOT, RH, FSDS, FLDS, PSRF, WIND' ")
#2nd ELM met data
parser.add_option("--mettype2", dest="met_type2", default="", \
                  help="e3sminput met data type 2 (default = '', options: CRU, GSWP3, GSWP3v1, Site, cplbypass, cplbypass_Site")
parser.add_option("--metdir2", dest="met_dir2", default="./", \
                  help="e3sminput met directory 2 (default = ./, i.e., under current directory)")
parser.add_option("--metdomain2", dest="met_domain2", default="", \
                  help="e3sminput met domain file (default = a file in current met_idir")
parser.add_option("--metheader2", dest="met_header2", default="", \
                  help="e3sminput directory (default = '', i.e., standard file name)")
parser.add_option("--varname2", dest="vars2", default="", \
                  help = "variable name to be merged/jointed: ONLY one of 'TBOT, PRECT, QBOT, RH, FSDS, FLDS, PSRF, WIND' ")
# other sources of data
parser.add_option("--user_mettype", dest="user_mettype", default="", \
                  help="user met data type (default = '', options: CRU, cplbypass, GSWP3, GSWP3v1, Site, ATS_h5, NCDC, daymet_nc4)")
parser.add_option("--user_metdir", dest="user_metdir", default="./", \
                  help="user-defined met directory (default = ./, i.e., under current directory)")
parser.add_option("--user_metfile", dest="user_metfile", default="", \
                  help="user-defined met file(s) under user_metdir (default = '', i.e., None)")
parser.add_option("--user_metdomain", dest="user_metdomain", default="", \
                  help="user-defined met domain file (default = a file in current met_dir")
parser.add_option("--user_varname", dest="user_vars", default="", \
                  help = "variable name to be merged/jointed: ONLY one of 'TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND' ")
# output options
parser.add_option("--lon", dest="lon", default=-999, \
                  help = " longitude to be reading/plotting, default -999 for first one")
parser.add_option("--lat", dest="lat", default=-999, \
                  help = " latitude to be reading/plotting, default -999 for first one")
parser.add_option("--yrmin", dest="yrmin", default=-9999, \
                  help = " longitude to be reading/plotting, default -999 for first one")
parser.add_option("--yrmax", dest="yrmax", default=-9999, \
                  help = " latitude to be reading/plotting, default -999 for first one")
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# First set of metdata, in ELM forcing format
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ts_yrly = 365 # daily (365), weekly(52), monthly(12)
DWM = 'daily'   # 'daily', 'weekly', 'monthly'
if (options.met_dir != '' and options.met_type != ''):
    vname_elm, t_elm, tunit_elm, sdata_elm, LONGXY, LATIXY = \
        ELMmet_1vardatas(options.vars, options.met_dir, options.met_type, \
                         options.met_header, options.met_domain, \
                         options.lon, options.lat)

    if ts_yrly > 0:
        tt, sdata = \
            DataTimeAggregation(t_elm, sdata_elm, DWMY=DWM)
        t_elm = deepcopy(tt)
        sdata_elm = deepcopy(sdata)
    
#--------------------------------------------------------------------------------------
# read-in 2nd metdata in ELM standard forms 
if (options.met_dir2 != '' and options.met_type2 != ''):
    if (options.vars2==''): options.vars2 = options.vars
    vname_elm2, t_elm2, tunit_elm2, sdata_elm2, LONGXY, LATIXY = \
        ELMmet_1vardatas(options.vars2, options.met_dir2, options.met_type2, \
                         options.met_header2, options.met_domain2, \
                         options.lon, options.lat)

    if ts_yrly > 0:
        tt, sdata = \
            DataTimeAggregation(t_elm2, sdata_elm2, DWMY=DWM)
        t_elm2 = deepcopy(tt)
        sdata_elm2 = deepcopy(sdata)

#--------------------------------------------------------------------------------------
# read-in 3rd metdata in ELM standard forms 
if (options.user_metdir != '' and options.user_mettype != ''):
    if (options.user_vars==''): options.user_vars = options.vars
    vname_elm3, t_elm3, tunit_elm3, sdata_elm3, LONGXY, LATIXY = \
        ELMmet_1vardatas(options.user_vars, options.user_metdir, options.user_mettype, \
                         options.user_metfile, options.user_metdomain, \
                         options.lon, options.lat)

    if ts_yrly >0:
        tt, sdata = \
            DataTimeAggregation(t_elm3, sdata_elm3, DWMY=DWM)
        t_elm3 = deepcopy(tt)
        sdata_elm3 = deepcopy(sdata)

#--------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------
# 3 data sets for similiarity and temporal up-scaling/down-scaling
if (options.met_dir != '' and options.met_type != ''):
    data1  = sdata_elm
    t1     = t_elm
    vname1 = vname_elm
    tunit1 = tunit_elm
else:
    vname1 = ''

if (options.met_dir2 != '' and options.met_type2 != ''):
    data2  = sdata_elm2
    t2     = t_elm2
    vname2 = vname_elm2
    tunit2 = tunit_elm2
else:
    vname2 = ''

if (options.user_metdir != '' and options.user_mettype != ''):
    data3  = sdata_elm3
    t3     = t_elm3
    vname3 = vname_elm3
    tunit3 = tunit_elm3
else:
    vname3 = ''

#--------------------------------------------------------------------------------------
# SEASONALITY
SEASONALLY=True
# ANNUALITY
ANNUALLY=True
#
if (vname1!=''):
    t1, data1, t1_yrly, data1_yrly, t1_seasonally, data1_seasonally = \
            DataTimePatterns(t1, data1, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly, \
                             yrmin=int(options.yrmin), yrmax=int(options.yrmax))
if (vname2!=''):
    t2, data2, t2_yrly, data2_yrly, t2_seasonally, data2_seasonally = \
            DataTimePatterns(t2, data2, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly, \
                             yrmin=int(options.yrmin), yrmax=int(options.yrmax))
if (vname3!=''):
    t3, data3, t3_yrly, data3_yrly, t3_seasonally, data3_seasonally = \
            DataTimePatterns(t3, data3, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly, \
                             yrmin=int(options.yrmin), yrmax=int(options.yrmax))

    
    

#--------------------------------------------------------------------------------------
# Year Matching

# by comparing statistics/regression fitting
#rmse, se, r2, offset, slope = DataSimilarity(data1, data2, t=t, t_unit=time_unit, v_name='PRECP', v_unit='mm/day', plt=plt)
#plt.show()

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# printing plot in PDF
if (options.plotting):
    snames = ['NOAA Station', 'GSWP3']
    if(vname3!=''): snames = ['NOAA Station', 'GSWP3', 'GSWP3-Daymet4']
    
    tt={}
    tt[snames[0]] = deepcopy(t_elm)
    tt[snames[1]] = deepcopy(t_elm2)
    if(vname3!=''): tt[snames[2]] = deepcopy(t_elm3)
    sdatas={}
    sdatas[snames[0]] = deepcopy(sdata_elm)
    sdatas[snames[1]] = deepcopy(sdata_elm2)
    if(vname3!=''): sdatas[snames[2]] = deepcopy(sdata_elm3)
    SinglePlotting(tt, 'days-since-0001-01-01', snames, vname_elm, sdatas, figno='1-orig')
    
    tt={}
    tt[snames[0]] = deepcopy(t1_yrly)
    tt[snames[1]] = deepcopy(t2_yrly)
    if(vname3!=''): tt[snames[2]] = deepcopy(t3_yrly)
    sdatas={}
    sdatas[snames[0]] = deepcopy(data1_yrly)
    sdatas[snames[1]] = deepcopy(data2_yrly)
    if(vname3!=''): sdatas[snames[2]] = deepcopy(data3_yrly)
    SinglePlotting(tt, 'Yearly', snames, vname_elm, sdatas, figno='2-yearly')
    
    tt={}
    tt[snames[0]] = deepcopy(t1_seasonally)
    tt[snames[1]] = deepcopy(t2_seasonally)
    if(vname3!=''): tt[snames[2]] = deepcopy(t3_seasonally)
    sdatas={}
    sdatas[snames[0]] = deepcopy(data1_seasonally)
    sdatas[snames[1]] = deepcopy(data2_seasonally)
    if(vname3!=''): sdatas[snames[2]] = deepcopy(data3_seasonally)
    SinglePlotting(tt, 'Seasonally-DOY', snames, vname_elm, sdatas, figno='3-seasonally')

#--------------------------------------------------------------------------------------


