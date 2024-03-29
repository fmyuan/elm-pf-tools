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
from Modules_CLMoutput_nc4 import CLM_NcRead_1simulation
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
    vnames=[vname_elm]
    # read-in metdata from CPL_BYPASS_FULL
    if ('cplbypass_' in met_type):
        cplbypass_dir=met_dir
        
        #cplbypass subtype
        cplbypass_mettype='GSWP3'
        if ('v1' in met_type): cplbypass_mettype='GSWP3v1'
        if ('daymet' in met_type): cplbypass_mettype='GSWP3_daymet'
        if ('v1' in met_type and 'daymet' in met_type): cplbypass_mettype='GSWP3v1_daymet'
        
        if(met_type=='cplbypass_Site'): 
            cplbypass_mettype='Site'
            if (vname_elm=='QBOT' or vname_elm=='estFLDS'): 
                vnames=['RH','TBOT','PSRF']
        else:
            if (vname_elm=='RH' or vname_elm=='estFLDS'): 
                vnames=['QBOT','TBOT','PSRF']
    
        if(met_type=='cplbypass_GSWP3' or met_type=='cplbypass_GSWP3v1'):
            # GSWP3 data QBOT needs redo, because its RH is NOT freezing adjusted
            # (When calculating RH from QBOT, RH often over 100 in freezing winter)
            if (vname_elm=='QBOT'):
                vnames=['QBOT','TBOT','PSRF']
        
        # read in
        zones,zlines, varsdims, vardatas = \
            Modules_metdata.clm_metdata_cplbypass_read( \
                            cplbypass_dir,cplbypass_mettype, vnames, lons=[lon], lats=[lat])

        # assign data to plotting variables
        nx=1
        ny=0
        for iz in zlines.keys():
            if np.isscalar(zlines[iz]):
                ny=ny+1
            else:
                ny=ny+len(zlines[iz])
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
           ('cplbypass_' not in met_type)):
        if (met_dir == './'):
            print('e3sm met. directory is the current')
        else:
            print('e3sm met.  directory: '+ met_dir)

        if(met_type=='Site'): 
            if (vname_elm=='QBOT' or vname_elm=='estFLDS'): 
                vnames=['RH','TBOT','PSRF']
        elif (vname_elm=='RH' or vname_elm=='estFLDS'): 
                vnames=['QBOT','TBOT','PSRF']
        if(met_type=='GSWP3' or met_type=='GSWP3v1'):
            # GSWP3 data QBOT needs redo, because its RH is NOT freezing adjusted
            # (When calculating RH from QBOT, RH often over 100 in freezing winter)
            if (vname_elm=='QBOT'):
                vnames=['QBOT','TBOT','PSRF']


        # read in
        ix,iy, varsdims, vardatas = \
            Modules_metdata.clm_metdata_read(met_dir, met_fileheader, met_type, met_domain, [lon], [lat],'')
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
    
    
    # read-in actual forcings from ELM simulation outputs
    elif (met_fileheader != '' and met_type == 'ELM'):
        tvarname = 'h0_time'    # variable name for time/timing
        nx, ny, nlgrnd, nldcmp, ncolumn, npft, vardatas, vardims, var_tunits = \
            CLM_NcRead_1simulation(met_dir, \
                               met_fileheader, \
                               'h0', \
                               False, \
                               vnames, \
                               -9999, -9999, \
                               False)
        
        vardatas['tunit'] = var_tunits
        
        # dimension max.
        if2dgrid = False
        if('topo' in vardatas.keys()):
            if(len(vardatas['topo'].shape)==2): if2dgrid = True
        
        nxy = nx*ny
        LATIXY = lat
        LONGXY = lon
        if('lat' in vardatas.keys()):LATIXY = vardatas['lat']
        if('lon' in vardatas.keys()):LONGXY = vardatas['lon']
    
    #------------------------------
    # if read-in ELM forcing data successfully
    if len(vardatas)>0:
        vars_list = list(vardatas.keys())
        
        t = vardatas[tvarname] # days since 1901-01-01 00:00
        tunit = vardatas['tunit']
        t0 = 0.0
        if ('1901-01-01' in tunit):
            tunit = tunit.replace('1901-','0001-')
            t0 = 1901*365 # converted from 'days-since-1901-01-01-00:00' to days since 0001-01-01 00:00

        for varname in vars_list:
            if vname_elm not in varname: continue
        
            if 'PRECTmms' in varname: 
                sdata_elm = deepcopy(vardatas[varname])*86400 # mm/s -> mm/d for convenient plotting
                varunit = 'mm/d'
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
                    vpsat_frzing=True
                    if 'GSWP3' in met_type and 'daymet' not in met_type: vpsat_frzing=False
                    #GSWP3 data of QBOT is NOT calculated from RH with vpsat adjusted for freezing air
                    sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=qbot,vpsat_frz=vpsat_frzing)
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
            #
            if 'QBOT' in vars_list and ('GSWP3' in met_type and 'daymet' not in met_type):
                #correction for GSWP3 QBOT data
                qbot = np.squeeze(vardatas['QBOT'])
                if 'TBOT' in vars_list:
                    tk = np.squeeze(vardatas['TBOT'])
                    if 'PSRF' in vars_list:
                        pres_pa = np.squeeze(vardatas['PSRF'])
                    else:
                        pres_pa = 101325.0
                else:
                    print('ERROR: for RH coverting from QBOT, air temperature is required')
                    sys.exit(-1)
                rh = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=qbot,vpsat_frz=False) # restore RH to orginal
                sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, rh_100=rh)  # re-cal. QBOT for GSWP3 dataset
            
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
    
    # there is difference for float32 or double type
    t=np.double(t)
    data=np.double(data)
    if t_unit.startswith("H"): 
        t2 = t/24.0
    elif t_unit.startswith("Y"): 
        t2 = t*365.0
    else:
        t2 = deepcopy(t)

    #
    t2=np.asarray(t2)/365.0
    t3=np.floor(t2)                 # in 'years'
    t2=(t2-np.floor(t2))*365.0      # still in original time-unit, but in DOY only (zero-based here)
    # doy 0, if from above, implies it be 365 and years be previous year, when it's not the first.
    # when it's the starting, it's OK, otherwise it may cause issue for regrouping and plotting
    t2 = np.round(t2,8)  # for unkown reason, t2 has some rounding error after conversions above
    #if(abs(t2[1]-t2[0])<1.0):
    #    t2 = t2 + 1                  # convert to 1-based DOY, if sub-daily data
    if (t2[0]!=0):  
        idx = np.where(t2==0)
        t3[idx] = t3[idx]-1
        t2[idx] = 365-1.0e-8
    else:
        t2 = t2 + 1

    
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

    idx1=np.where(np.isnan(data1))
    if len(idx1[0]>0):
        dd1=np.delete(data1,idx1)
        dd2=np.delete(data2,idx1)
        
        idx2=np.where(np.isnan(dd2))
        if len(idx2[0]>0):
            dd1=np.delete(dd1,idx2)
            dd2=np.delete(dd2,idx2)
    else:
        dd1=data1
        dd2=data2

    fit_model = np.polyfit(dd1, dd2, 1, cov=True)
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
        One2OnePlotting(plt, nrow, ncol, isub, dd1, dd2, datalabel, plotlabel)

        
        # T-series or NONE
        isub = isub + 1
        plotlabel = 'T-series'
        if (len(t)<=0):
            t_unit = 'none'
            t = range(0,len(data1)-1)
        
        sdata = np.swapaxes(np.vstack((data1,data2)),0,1)
        varnames = [v_name+'_1',v_name+'_2']
        TimeGridedVarPlotting(plt, nrow, ncol, isub, t, t_unit, sdata, \
                    varnames, plotlabel)


    # return indicators, if any
    return rmse, me, intcp, slope

##*****************************************************************************


#-------------------Parse options-----------------------------------------------


# ---------------------------------------------------------------

parser = OptionParser()

# E3SM met data 
parser.add_option("--mettype", dest="met_type", default="CRU", \
                  help="e3sminput met data type (default = 'CRU', options: CRU, GSWP3, GSWP3v1, Site, GSWP3_daymet, cplbypass_*")
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
                  help="e3sminput met data type 2 (default = '', options: CRU, GSWP3, GSWP3v1, Site, GSWP3_daymet, cplbypass_*")
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
                  help="user met data type (default = '', options: CRU, GSWP3, GSWP3v1, Site, GSWP3_daymet, cplbypass_*, ATS_h5, NCDC, daymet_nc4)")
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
parser.add_option("--similarity", dest="datasimilarity", default=False, \
                  help = "comparing 2 datasets", action="store_true")
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
    while (len(sdata_elm.shape)>1): #2-D data, requring aggregating over spatial dimension
        idim = np.where(np.asarray(sdata_elm.shape)!=len(t_elm))[0]
        sdata_elm = np.nanmean(sdata_elm, axis=idim[0])

    if ts_yrly > 0:
        tt, sdata = \
            DataTimeAggregation(t_elm, sdata_elm, DWMY=DWM)
        t_elm = deepcopy(tt)
        sdata_elm = deepcopy(sdata)
        del tt, sdata
    
#--------------------------------------------------------------------------------------
# read-in 2nd metdata in ELM standard forms 
if (options.met_dir2 != '' and options.met_type2 != ''):
    if (options.vars2==''): options.vars2 = options.vars
    if (options.met_header2==''): options.met_header2 = options.met_header
    if (options.met_domain2==''): options.met_domain2 = options.met_domain
    
    vname_elm2, t_elm2, tunit_elm2, sdata_elm2, LONGXY, LATIXY = \
        ELMmet_1vardatas(options.vars2, options.met_dir2, options.met_type2, \
                         options.met_header2, options.met_domain2, \
                         options.lon, options.lat)
    
    while (len(sdata_elm2.shape)>1): #2-D data, requring aggregating over spatial dimension
        idim = np.where(np.asarray(sdata_elm2.shape)!=len(t_elm2))[0]
        sdata_elm2 = np.nanmean(sdata_elm2, axis=idim[0])

    if ts_yrly > 0:
        tt, sdata = \
            DataTimeAggregation(t_elm2, sdata_elm2, DWMY=DWM)
        t_elm2 = deepcopy(tt)
        sdata_elm2 = deepcopy(sdata)
        del tt, sdata

#--------------------------------------------------------------------------------------
# read-in 3rd metdata in ELM standard forms 
if (options.user_metdir != '' and options.user_mettype != ''):
    if (options.user_vars==''): options.user_vars = options.vars
    if (options.user_metfile==''): options.user_metfile = options.met_header
    if (options.user_metdomain==''): options.user_metdomain = options.met_domain
    
    vname_elm3, t_elm3, tunit_elm3, sdata_elm3, LONGXY, LATIXY = \
        ELMmet_1vardatas(options.user_vars, options.user_metdir, options.user_mettype, \
                         options.user_metfile, options.user_metdomain, \
                         options.lon, options.lat)
    while (len(sdata_elm3.shape)>1): #2-D data, requring aggregating over spatial dimension
        idim = np.where(np.asarray(sdata_elm3.shape)!=len(t_elm3))[0]
        sdata_elm3 = np.nanmean(sdata_elm3, axis=idim[0])
        

    if ts_yrly >0:
        tt, sdata = \
            DataTimeAggregation(t_elm3, sdata_elm3, DWMY=DWM)
        t_elm3 = deepcopy(tt)
        sdata_elm3 = deepcopy(sdata)
        del tt, sdata

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

# 
# comparing statistics/regression fitting
if options.datasimilarity and (vname1!='' and vname2!='' and vname3==''):
    tt, idx1, idx2 = np.intersect1d(t1,t2, return_indices=True)
    if (tt.size>5):
        dd1 = data1[idx1]
        dd2 = data2[idx2]
        
        rmse, se, offset, slope = DataSimilarity(dd1, dd2, t=tt, t_unit=tunit1, v_name=vname1, v_unit='', plt=plt)
        plt.show()


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# printing plot in PDF
if (options.plotting):
    snames = ['1-'+options.met_type, '2-'+options.met_type2]
    if(vname3!=''): 
        snames = ['1-'+options.met_type, '2-'+options.met_type2, '3-'+options.user_mettype]
    
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


