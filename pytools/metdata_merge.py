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
        sys.exit(-1)
    
    if (data_method == 'offset'):
        offset = target - src# 
        sdata_jointting = sdata_jointting + offset
    elif (data_method == 'nanmean'):# 
        if (not np.any(src==0.0)):
            ratio = np.nanmean([src,target], axis=0)/src
            sdata_jointting = sdata_jointting*ratio
        else:
            offset = target - src# 
            sdata_jointting = sdata_jointting + offset/2.0
    elif (data_method == 'ratio'):
        if (not np.any(src==0.0)):
            ratio = target/src# 
            sdata_jointting = sdata_jointting*ratio
        else:
            offset = target - src# 
            sdata_jointting = sdata_jointting + offset
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
    # doy 0, if from above, implies it be 365 and years be previous year, when it's not the first.
    # when it's the starting, it's OK, otherwise it may cause issue for regrouping and plotting
    if (t2[0]!=0):  
        idx = np.where(t2==0)
        t3[idx] = t3[idx]-1
        t2[idx] = 365-1.0e-8
        
    
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

parser.add_option("--e3sm_metdomain", dest="met_domain", default="", \
                  help="e3sminput met domain file (default = a file in current met_dir")
# E3SM met data 
parser.add_option("--mettype", dest="met_type", default="CRU", \
                  help="e3sminput met data type (default = 'CRU', options: CRU, GSWP3, GSWP3v1, GSWP3, Site, and *_daymet/cplbypass_*)")
parser.add_option("--e3sm_metdir", dest="met_dir", default="./", \
                  help="e3sminput met directory (default = ./, i.e., under current directory)")
parser.add_option("--e3sm_metheader", dest="met_header", default="", \
                  help="e3sminput directory (default = '', i.e., standard file name)")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name to be merged/jointed: ONLY one of 'TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND' ")
parser.add_option("--var_depended", dest="var_depended", default="", \
                  help = "variable name to be merged/jointed: ONLY one of 'TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND' ")
parser.add_option("--lon", dest="lon", default=-999, \
                  help = " longitude to be reading/plotting, default -999 for first one")
parser.add_option("--lat", dest="lat", default=-999, \
                  help = " latitude to be reading/plotting, default -999 for first one")

# other sources of data to be pull into merged E3SM data
parser.add_option("--user_mettype", dest="user_mettype", default="", \
                  help="e3sminput met data type (default = 'CRU', options: CRU, GSWP3, GSWP3v1, GSWP3, Site, and *_daymet/cplbypass_*, None)")
parser.add_option("--user_metdir", dest="user_metdir", default="./", \
                  help="user-defined met directory (default = ./, i.e., under current directory)")
parser.add_option("--user_metfile", dest="user_metfile", default="", \
                  help="user-defined met file(s) under user_metdir (default = '', i.e., None)")
#
parser.add_option("--extending", dest="extending", default=False, \
                  help = "extending merged-data into maximal period of time", action="store_true")
parser.add_option("--plotting", dest="plotting", default=False, \
                  help = "plotting data for checking", action="store_true")
parser.add_option("--nc_create", dest="nc_create", default=False, \
                  help = "output new nc file", action="store_true")
parser.add_option("--nc_write", dest="nc_write", default=False, \
                  help = "output to a existed nc file", action="store_true")
parser.add_option("--ncout_mettype", dest="nc_write_mettype", default="", \
                  help = "output to nc files in defined format (default = '', i.e., as original)")
parser.add_option("--ncout_stdmetdir", dest="nc_write_stdmetdir", default="", \
                  help = "output nc met files templates directory (default = '', i.e., as in original)")
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

elmmettype = ['CRU', 'GSWP3', 'GSWP3v1', 'Site', \
              'GSWP3v1_Daymet','GSWP3_Daymet', 'GSWP3_daymet4']  # v1 is from 1980-2010, _Daymet is in half-degree, and _daymet4 is in 1km res, and NA only
  # + 'cplbypass_*' of above all elmmettype
lon = float(options.lon)
lat = float(options.lat)
met_domain=options.met_domain

metdir=[options.met_dir]
metfileheader=[options.met_header]
mettype=[options.met_type]

# if for merging mets are both ELM met-type
if options.user_mettype in elmmettype or 'cplbypass_' in options.user_mettype:
    metdir=[options.met_dir, options.user_metdir]
    if(options.met_header=='' and options.user_metfile==''):
        metfileheader=[' ', ' ']
    else:
        metfileheader=[options.met_header, options.user_metfile]
    mettype=[options.met_type, options.user_mettype]
    
    
for imet in range(len(mettype)):
    imet_type = mettype[imet]
    imet_dir = metdir[imet]
    imet_header = metfileheader[imet]
    t = []
    vardatas = []
    # read-in metdata from CPL_BYPASS_FULL
    if ('cplbypass_' in imet_type):
        cplbypass_dir=imet_dir
        cplbypass_mettype='GSWP3'
        if ('daymet' in imet_type.lower()): cplbypass_mettype='GSWP3_daymet'
        if ('v1' in imet_type): cplbypass_mettype='GSWP3v1'
        if ('v1' in imet_type and 'daymet' in imet_type.lower()): cplbypass_mettype='GSWP3v1_daymet'
        
        vnames=[vname_elm]
        if('cplbypass_Site' in imet_type): 
            cplbypass_mettype='Site'
            if (vname_elm=='QBOT' or vname_elm=='estFLDS'): 
                vnames=['RH','TBOT','PSRF']
        else:
            if (vname_elm=='RH' or vname_elm=='estFLDS'): 
                vnames=['QBOT','TBOT','PSRF']
        if(imet_type=='cplbypass_GSWP3' or imet_type=='cplbypass_GSWP3v1'):
            # GSWP3 data QBOT needs redo, because its RH is NOT freezing adjusted
            # (When calculating RH from QBOT, RH often over 100 in freezing winter)
            if (vname_elm=='QBOT'):
                vnames=['QBOT','TBOT','PSRF']
        
        if(options.var_depended!='' and options.var_depended not in vnames): 
            vnames=vnames+[options.var_depended]
        
        #if (vname_elm=='PRECT'): # for excluding data in winter
        #    vnames=['PRECT','TBOT']
    
        # read in
        zones, zlines, varsdims, vardatas = \
            Modules_metdata.clm_metdata_cplbypass_read( \
                                cplbypass_dir, cplbypass_mettype, vnames, lons=[lon], lats=[lat])
        
        # in case PRECT is actually rainfall only
        if 'TBOT' in vardatas.keys() and 'PRECT' in vardatas.keys():
            idx=np.where(vardatas['TBOT']<=273.15)
            vardatas['PRECT'][idx] = np.nan
        
    
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
    elif (('CRU' in imet_type or \
        'GSWP3' in imet_type or \
        'Site' in imet_type) and \
        ('cplbypass_' not in imet_type)):
        if (imet_dir == './'):
            print('e3sm met. directory is the current')
        else:
            print('e3sm met.  directory: '+ imet_dir)

        if('Site' in imet_type): 
            if (vname_elm=='QBOT' or vname_elm=='estFLDS'): 
                vnames=['RH','TBOT','PSRF']
        elif (vname_elm=='RH' or vname_elm=='estFLDS'): 
                vnames=['QBOT','TBOT','PSRF']
        if(imet_type=='GSWP3' or imet_type=='GSWP3v1'):
            # GSWP3 data QBOT needs redo, because its RH is NOT freezing adjusted
            # (When calculating RH from QBOT, RH often over 100 in freezing winter)
            if (vname_elm=='QBOT'):
                vnames=['QBOT','TBOT','PSRF']
    
        # read in
        ix,iy, varsdims, vardatas = \
            Modules_metdata.clm_metdata_read(imet_dir,imet_header, imet_type, met_domain, [lon], [lat],'')
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
    # if read-in data successfully
    if len(vardatas)>0:
        vars_list = list(vardatas.keys())
        
        t_elm = vardatas[tvarname]
        tunit_elm = vardatas['tunit']
        if ('1901-01-01' in tunit_elm):
            tunit_elm = tunit_elm.replace('1901-','0001-')
            t0 = 1901*365 # converted from 'days-since-1901-01-01-00:00' to days since 0001-01-01 00:00
            t_elm = np.asarray(t_elm)+t0
        tvarname_elm = tvarname
        for varname in vars_list:
            if vname_elm not in varname: continue
            
            sdata_elm = deepcopy(vardatas[varname])
            sdata_elm = np.squeeze(sdata_elm)
            
            vname_elm = varname
            break
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
                    if ('GSWP3' in imet_type and 'daymet' not in imet_type):
                        sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=qbot, vpsat_frz=False)
                    else:
                        sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=qbot)
                    del qbot#, tk, pres_pa
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
                del rh #, tk, pres_pa
            #
            if 'QBOT' in vars_list and ('GSWP3' in imet_type and 'daymet' not in imet_type):
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
                del qbot, rh #, tk, pres_pa

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
            del pres_pa, qbot, rh
        
        #
        if(options.var_depended!=''): 
            sdata_depended = vardatas[options.var_depended]
        
        # clean up large data in memory
        del vardatas
        
        # aggregate half-hourly to 3-hourly (e.g. to GSWP3 datasets)
        if imet_type == 'GSWP3' or imet_type=='GSWP3v1':
            t_elm = np.reshape(t_elm, [-1,6])
            t_elm = np.squeeze(t_elm[:,0])
            sdata_elm = np.reshape(sdata_elm, [-1,6])
            sdata_elm = np.mean(sdata_elm,axis=1)
            if(options.var_depended!=''): 
                sdata_depended = np.reshape(sdata_depended, [-1,6])
                sdata_depended = np.mean(sdata_depended,axis=1)

        
        if len(imet_type)>1:
            if imet == 0:
                t_elm1 = deepcopy(t_elm)
                tunit_elm1 = deepcopy(tunit_elm)
                tvarname_elm1 = tvarname_elm
                vname_elm1 = vname_elm
                sdata_elm1 = deepcopy(sdata_elm)
                if(options.var_depended!=''): 
                    sdata_depended1 = deepcopy(sdata_depended)
                    pres_elm = deepcopy(pres_pa)
                    tk_elm = deepcopy(tk)
            
            elif imet == 1: # if this is the case, the following user_mettype won't read
                t_user = deepcopy(t_elm)
                tunit_user = deepcopy(tunit_elm)
                tvarname_user = tvarname_elm
                vname_user = vname_elm
                sdata_user = deepcopy(sdata_elm)
                if(options.var_depended!=''): sdata_depended_user = deepcopy(sdata_depended)
                 
                #copy back elm1
                t_elm = deepcopy(t_elm1)
                tunit_elm = deepcopy(tunit_elm1)
                tvarname_elm = tvarname_elm1
                vname_elm = vname_elm1
                sdata_elm = deepcopy(sdata_elm1)
                if(options.var_depended!=''): 
                    sdata_depended_elm = deepcopy(sdata_depended1)
                    del sdata_depended1
                del t_elm1, tunit_elm1, tvarname_elm1, vname_elm1, sdata_elm1
            
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# read-in metdata from totally user-defined data
if ('NONE' in options.user_mettype):
    odata_header,odata = \
        singleReadCsvfile('5_minute_data_v21.csv')
    
    
        #singleNCDCReadCsvfile('2022211_MysVanKarem_initproc.csv','NCDC_metric','-9999.99') # 'metric' refers to NCDC data converted to metric system already
        #singleNCDCReadCsvfile('GHCND_Alert_initproc.csv','GHCND','-999.99')
    #site_header = ("LATITUDE","LONGITUDE","ELEVATION")
    #data_header = ("YEAR","MONTH","DOY","PRCP","SNWD","TAVG","TMAX","TMIN")
    
    # need to recalculate time
    tvarname = 'time'    # variable name for time/timing
    
    vardatas = {}
    vardatas['time']=(odata['YEAR']-1901.0)*365.0+(odata['DOY']-1.0)  # days since 1901-01-01 00:00:00
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

    t = vardatas[tvarname] # days since 1901-01-01-00000
    tunit = 'days since 0001-01-01 00:00'
    t0 = 1901*365 # converted to days since 0000-01-01-00000

    if(options.vars=='TBOT'):
        vname_elm = 'TBOTd'
        varunit = 'K'
    elif(options.vars=='PRECT'):
        vname_elm = 'PRCPd'
        varunit = 'mm/d'
    else:
        print('NOT supported variable name, should be one of : ')
        print('TBOT, PRECT')
        sys.exit(-1)

    for varname in vars_list:
        if vname_user not in varname: continue
        
        vardata = vardatas[varname]
        sdata_user = deepcopy(vardata)
        #SubPlotting(t, time_unit, varname, varunit, sdata_user) # for checking
        
        vname_user = varname
        break
        
    #
    t_user = np.asarray(t)+t0
    tunit_user = tunit
    del vardata, vardatas, t
    
#--------------------------------------------------------------------------------------
# read-in metdata from ATS daily forcing

elif ('ATS_h5' in options.user_mettype):
    #hdfname = 'CESM-RCP8.5-2006-2100_dm1985-2015.h5'
    hdfname = options.user_metdir+'/'+options.user_metfile
    varnames=['all']
    vardatas = \
        Read1hdf(hdfname, varnames)

    # 
    t_name='time [s]'
    tunit='[s]'
    t = vardatas[t_name]
    varnames = [x for x in vardatas.keys() if x!=t_name]
    nvars = len(varnames)
    nt = len(t)
    sdata_user = np.zeros(nt)

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

    t0 = 2006*365 # ==>days since 0001-01-01-000000, USER-defined
    
    for ivar in range(0,nvars):
        if vname_ats not in varnames[ivar]: 
            continue #skip
        else:
            print ('ATS varname: ', varnames[ivar])
        sdata_user = vardatas[varnames[ivar]]
        varunit = varnames[ivar].strip().split(" ")[-1]
        varname = varnames[ivar].strip().replace(varunit,'')

        if 'rain' in varname or 'snow' in varname or 'precipitation' in varname: 
            sdata_user = deepcopy(sdata_user)*86400000 # m/s -> mm/day
        #SinglePlotting(t, time_unit, varname, varunit, sdata_user) # for checking
        
        if 'relative humidity' in varname:
            sdata_user = sdata_user*100.
        
        vname_user = varnames[ivar]
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
            sdata_user = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=[], rh_100=rh)
        #
        if(options.vars=='estFLDS'):
            sdata_user = Modules_metdata.calcFLDS(tk, pres_pa, q_kgkg=[], rh_100=rh)
    
    #
    t_user = np.asarray(t)/86400.0+t0
    tunit_user = 'days since 0001-01-01 00:00'

    del vardatas, t
    #

#---------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
# 2 data sets for similiarity and temporal up-scaling/down-scaling
if (options.user_mettype !=''):
    data1 = sdata_user
    if (options.var_depended!=''): 
        data1_depended=np.squeeze(sdata_depended_user)
        #del sdata_depended_user
    t1 = t_user
    vname1=vname_user
    tunit1 = tunit_user
else:
    t_user=[]
    sdata_user=[]
    data1=[]
    t1=[]
    data1_yrly=[]
    t1_yrly=[]
    data1_seasonally=[]
    t1_seasonally=[]
    tunit1=''
    vname1=''

if options.met_type in elmmettype or 'cplbypass_' in options.met_type:
    data2 = np.squeeze(sdata_elm)
    if (options.var_depended!=''): 
        data2_depended=np.squeeze(sdata_depended_elm)
        #del sdata_depended_elm
    t2 = t_elm
    vname2=vname_elm
    tunit2=tunit_elm
else:
    t_elm=[]
    sdata_elm=[]
    data2=[]
    t2=[]
    data2_yrly=[]
    t2_yrly=[]
    data2_seasonally=[]
    t2_seasonally=[]
    tunit2=''
    vname2=''

#--------------------------------------------------------------------------------------
# SEASONALITY
SEASONALLY=True
# ANNUALITY
ANNUALLY=True
#
dts1 = np.unique(np.diff(np.round(t1,10)))  # must be in 'days'
dts2 = np.unique(np.diff(np.round(t2,10)))  # must be in 'days'
ts_yrly = 365
#
detrending=''
if (vname1=='NONE' or vname2=='NONE'): detrending=''

if (detrending !=''):
    if (vname1!=''):
        t1, data1, t1_yrly, data1_yrly, t1_seasonally, data1_seasonally, fitfunc_data1 = \
            DataTimePatterns(t1, data1, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly, detrending=detrending) 
        # (optional) detrending: '', 'yearly','all', with output 'fitfunc_data1' (np.poly1d, i.e. fitfunc_data1(t1))
    if (vname2!=''):
        t2, data2, t2_yrly, data2_yrly, t2_seasonally, data2_seasonally, fitfunc_data2 = \
            DataTimePatterns(t2, data2, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly, detrending=detrending)
    noted = 'detrended'
else:
    
    if (vname1!=''):
        t1, data1, t1_yrly, data1_yrly, t1_seasonally, data1_seasonally = \
            DataTimePatterns(t1, data1, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly)
        if(options.var_depended!=''): 
            t1, data1_depended, t1_yrly, data1_depended_yrly, t1_seasonally, data1_depended_seasonally = \
                DataTimePatterns(t1, data1_depended, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly)
    if (vname2!=''):
        t2, data2, t2_yrly, data2_yrly, t2_seasonally, data2_seasonally = \
            DataTimePatterns(t2, data2, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly)
        if(options.var_depended!=''): 
            t2, data2_depended, t2_yrly, data2_depended_yrly, t2_seasonally, data2_depended_seasonally = \
                DataTimePatterns(t2, data2_depended, SEASONALLY, ANNUALLY, ts_yrly=ts_yrly)

    noted = 'orig'
    
    

#--------------------------------------------------------------------------------------
# Merging

# by jointing 2 datasets directly, from data1-->data2
if(vname1!='' and vname2!=''):
    BY_JOINTING=True
else:
    BY_JOINTING=False
if (BY_JOINTING):
    noted = noted+'_jointed'

    # Source time/data for up-/-down temporally-scaling
    t_src2 = deepcopy(t_elm)
    sdata_src2 = deepcopy(sdata_elm)
    
    # joint-smooth-transition
    if detrending!='':
        # time-interval sync'ed data2/data1
        tj=(t2[-1]+t1[0])/2.0
        jointgap_data = fitfunc_data2(tj)-fitfunc_data1(tj)
        data2 = data2 - jointgap_data
        sdata_src2 = sdata_src2 - jointgap_data
        
        # integrated data2/data1
        jointgap_yrly = np.nanmean(data2_yrly)-np.nanmean(data1_yrly)
        data2_yrly = data2_yrly - jointgap_yrly
        jointgap_seasonally = np.nanmean(data2_seasonally)-np.nanmean(data1_seasonally)
        data2_seasonally = data2_seasonally - jointgap_seasonally
        
    else:
        jointgap_data = 0.0
        jointgap_yrly = 0.0
        jointgap_seasonally = 0.0
    
    
    # to be jointed/temporally-scaling (pulling data1 into data2)
    t_jointed = deepcopy(t_src2[t_src2<t_user[0]-dts1/2.0])
    sdata_jointed = deepcopy(sdata_src2[t_src2<t_user[0]-dts1/2.0])
    
    # loop through 't1' by year
    iyr_t_src = np.where(t2_yrly==t1_yrly[0])[0][0]
    for iyr in np.floor(t1_yrly):
        
        if iyr<t2_yrly[0]: continue  # skip data prior to min of t2
        
        dt = np.round(np.max(np.diff(t1_seasonally)),8)
        for it in range(len(t1_seasonally)):
            iday = t1_seasonally[it]
            td0_t1 = iyr*365.0+it*dt                  # starting point of a day in t1
            td0_t2 = t2_yrly[iyr_t_src]*365.0+it*dt   # starting point of a day in t2
            
            idx_target = np.where((t_user>=td0_t1) & (t_user<(td0_t1+dt)))
            idx_src = np.where((t_src2>=td0_t2) & (t_src2<(td0_t2+dt)))
            
            t_jointting = td0_t1 + (t_src2[idx_src]-td0_t2) # the second portion is sub-daily (fraction) time
            t_jointed = np.hstack((t_jointed, t_jointting))
            
            # temporal down-scaling:
            sdata_jointting = sdata_src2[idx_src]
            
            # timestep interval is not same, need to use aggregated data (data1, data2)
            if len(idx_src[0]) != len(idx_target[0]):
                # 'sdata_jointing' - in/out (for sub-time variation), 
                # 'src/target' are time-sync'ed datasets to provide scalors (a moving window)
                twidth = 1.0*dt  # 1-d moving-window (for 'data_method' of 'smooth', t spans 2 or more timesteps -TODO)
                idx = np.where((t2>td0_t2) & (t2<=(td0_t2+twidth)))
                src = data2[idx]
                if(options.var_depended!=''): src_depended = data2_depended[idx]
                
                if (vname1!='NONE' and src.size>0):
                    idx = np.where((t1>td0_t1) & (t1<=(td0_t1+twidth)))
                    target = data1[idx]
                    if(options.var_depended!=''): target_depended = data1_depended[idx]
                    
                    if(len(src)!=len(target)):
                        print('Error: src/target data NOT in same length!')
                        sys.exit(-1)
                    
                    # 'data_method': 'offset' (default),'nanmean', 'ratio', etc
                    if sum(np.isnan(target))<len(target): # if aimed data available
                        sdata_jointting = DataTimeDown(sdata_jointting, src, target, data_method='offset')
                
                #else:
                    # do nothing if 'data1' not exists, then just simply copy data
            
            else: # didn't do de-trending
                src = deepcopy(sdata_jointting)
                target = sdata_user[idx_target]
                
                if sum(np.isnan(target))<len(target): # if aimed data available
                    #sdata_jointting = DataTimeDown(sdata_jointting, src, target, data_method='nanmean')
                    sdata_jointting = np.nanmean([src,target], axis=0)
                # very specific to correct daymet-derived QBOT or RH datasets
                if((options.vars=='QBOT' or options.vars=='RH') \
                       and options.var_depended=='TBOT' and 'daymet' in options.met_type.lower()):
                        idx=np.where(target_depended>273.15) # i.e. if air NOT freezing, pick the merging data, otherwise averaging as above
                        if len(idx[0])>0: sdata_jointting[idx]=target[idx]
                    
                
            
            sdata_jointed = np.hstack((sdata_jointed, sdata_jointting))
            

        #DONE 'for iday in np.floor(t2_seasonally)'
        #cycling 'sdata_elm' throughout by year-order
        iyr_t_src = iyr_t_src + 1
        if(iyr_t_src>=len(t2_yrly)): 
            if not options.extending: break
            iyr_t_src = np.where(t2_yrly==t1_yrly[0])[0][0]
    
    #DONE 'for iyr in np.floor(t1_yrly)'
    # some simple checkings
    if (options.vars=='RH' and 'Site' not in options.nc_write_mettype):
        # var name has modified for non-Site met-type
        options.vars='QBOT'
        sdata_jointed = Modules_metdata.convertHumidity(tk_elm, pres_elm, rh_100=sdata_jointed)

    if(options.vars.strip()=='RH'):
        idx=np.where(sdata_jointed>100.0)
        sdata_jointed[idx]=100.0
        idx=np.where(sdata_jointed<5.0)
        sdata_jointed[idx]=5.0
    if('FLDS' in options.vars.strip() or 'FSDS' in options.vars.strip() \
       or 'QBOT' in options.vars.strip() or 'PRECT' in options.vars.strip() \
       or 'WIND' in options.vars.strip()):
        idx=np.where(sdata_jointed<0.0)
        sdata_jointed[idx]=0.0
    if(options.vars.strip()=='TBOT'):
        idx=np.where(sdata_jointed<203.15)
        sdata_jointed[idx]=203.15
    if(options.vars.strip()=='PSRF'):
        idx=np.where(sdata_jointed<50000.0)
        sdata_jointed[idx]=101325.0
        
    
    if (options.plotting):
        snames = ['Jointed','Original-1','Original-2']
        tt={}
        tt[snames[0]] = deepcopy(t_jointed)
        tt[snames[1]] = deepcopy(t_elm)
        tt[snames[2]] = deepcopy(t_user)
        sdatas={}
        sdatas[snames[0]] = deepcopy(sdata_jointed)
        sdatas[snames[1]] = deepcopy(sdata_elm)
        sdatas[snames[2]] = deepcopy(sdata_user)
        
        SinglePlotting(tt, 'days-since-0001-01-01', snames, vname_elm, sdatas)
        
        txx, dataxx, txx_yrly, dataxx_yrly, txx_seasonally, dataxx_seasonally = \
            DataTimePatterns(t_jointed, sdata_jointed, SEASONALLY, ANNUALLY, ts_yrly=365) 
        tt={}
        tt[snames[0]] = deepcopy(txx)
        tt[snames[1]] = deepcopy(t2)
        tt[snames[2]] = deepcopy(t1)
        sdatas={}
        sdatas[snames[0]] = deepcopy(dataxx)
        sdatas[snames[1]] = deepcopy(data2)
        sdatas[snames[2]] = deepcopy(data1)
        SinglePlotting(tt, 'days-since-0001-01-01', snames, vname_elm, sdatas, figno='XX-1')

        tx2, datax2, tx2_yrly, datax2_yrly, tx2_seasonally, datax2_seasonally = \
            DataTimePatterns(t_src2, sdata_elm, SEASONALLY, ANNUALLY, ts_yrly=365) 
        tx1, datax1, tx1_yrly, datax1_yrly, tx1_seasonally, datax1_seasonally = \
            DataTimePatterns(t_user, sdata_user, SEASONALLY, ANNUALLY, ts_yrly=365) 

        tt={}
        tt[snames[0]] = deepcopy(txx_yrly)
        tt[snames[1]] = deepcopy(tx2_yrly)
        tt[snames[2]] = deepcopy(tx1_yrly)
        sdatas={}
        sdatas[snames[0]] = deepcopy(dataxx_yrly)
        sdatas[snames[1]] = deepcopy(datax2_yrly)
        sdatas[snames[2]] = deepcopy(datax1_yrly)
        SinglePlotting(tt, 'YEAR', snames, vname_elm, sdatas, figno='XX-2')

        tt={}
        tt[snames[0]] = deepcopy(txx_seasonally)
        tt[snames[1]] = deepcopy(tx2_seasonally)
        tt[snames[2]] = deepcopy(tx1_seasonally)
        sdatas={}
        sdatas[snames[0]] = deepcopy(dataxx_seasonally)
        sdatas[snames[1]] = deepcopy(datax2_seasonally)
        sdatas[snames[2]] = deepcopy(datax1_seasonally)
        SinglePlotting(tt, 'DOY', snames, vname_elm, sdatas, figno='XX-3')

else: # no jointing of datasets
    t_jointed = deepcopy(t_elm)
    sdata_jointed = deepcopy(sdata_elm)

    
#--------------------------------------------------------------------------------------
# printing plot in PDF
# for checking data-jointing
if (options.plotting):
    snames = ['elm','user']
    tt={}
    tt[snames[0]] = deepcopy(t_elm)
    tt[snames[1]] = deepcopy(t_user)
    sdatas={}
    sdatas[snames[0]] = deepcopy(sdata_elm)
    sdatas[snames[1]] = deepcopy(sdata_user)
    SinglePlotting(tt, 'days-since-0001-01-01', snames, vname_elm, sdatas, figno='1-orig')
    
    tt={}
    tt[snames[0]] = deepcopy(t2)
    tt[snames[1]] = deepcopy(t1)
    sdatas={}
    sdatas[snames[0]] = deepcopy(data2)
    sdatas[snames[1]] = deepcopy(data1)
    SinglePlotting(tt, 'days-since-0001-01-01', snames, vname_elm, sdatas, figno='1A-orig')
    
    tt={}
    tt[snames[0]] = deepcopy(t2_yrly)
    tt[snames[1]] = deepcopy(t1_yrly)
    sdatas={}
    sdatas[snames[0]] = deepcopy(data2_yrly)
    sdatas[snames[1]] = deepcopy(data1_yrly)
    SinglePlotting(tt, 'Yearly', snames, vname_elm, sdatas, figno='2-'+noted)
    
    tt={}
    tt[snames[0]] = deepcopy(t2_seasonally)
    tt[snames[1]] = deepcopy(t1_seasonally)
    sdatas={}
    sdatas[snames[0]] = deepcopy(data2_seasonally)
    sdatas[snames[1]] = deepcopy(data1_seasonally)
    SinglePlotting(tt, 'Seasonally-DOY', snames, vname_elm, sdatas, figno='3-'+noted)

#--------------------------------------------------------------------------------------
# save in ELM forcing data format
if (options.nc_create or options.nc_write):

    # if cutoff data for whatever reason to output
    NCOUT_CUTOFF=False
    if NCOUT_CUTOFF:
        yr=np.floor(t_jointed/365.0)
        idx=np.where(yr<=2015)
        t_jointed = t_jointed[idx]
        sdata_jointed = sdata_jointed[idx]
    
    
    out_mettype = options.nc_write_mettype
    if out_mettype =='': out_mettype = options.met_type
    if out_mettype==options.met_type and 'cplbypass_' not in options.met_type:
        NCOUT_ORIGINAL=True
    else:
        NCOUT_ORIGINAL=False
    if 'cplbypass_' in options.nc_write_mettype:
        NCOUT_CPLBYPASS=True
    else:
        NCOUT_CPLBYPASS=False

    if (options.nc_create and options.nc_write):
        print('Error: cannot have both "--nc_create" and "--nc_write" ')
        sys.exit(-1)
    elif(NCOUT_ORIGINAL or NCOUT_CPLBYPASS):
        if (options.nc_create):
            print('Create new ELM forcing data ? ', options.nc_create)
        elif (options.nc_write):
            print('Write to existed ELM forcing data ? ', options.nc_write)
    
    # get a template ELM forcing data nc file
    # 
    varname = vname_elm
    if('Site' not in out_mettype and varname=='RH'): varname='QBOT'
    ncfilein_cplbypass = ''
    out_metdir = options.nc_write_stdmetdir
    if out_metdir=='': out_metdir='./'
    if 'GSWP3' in out_mettype:
        if (varname == 'FSDS'):
            fdir = options.met_dir+'/Solar3Hrly/' # where to access template nc file
        elif (varname == 'PRECTmms'):
            fdir = options.met_dir+'/Precip3Hrly/'
        else:
            fdir = options.met_dir+'/TPHWL3Hrly/'
        
        fdirheader = options.met_dir+'/GSWP3_'+varname+'_'
        if ('daymet4' in out_mettype.lower()):# from 1980-2014, and 1km resolution
            fdirheader = options.met_dir+'/GSWP3_daymet4_'+varname+'_'
        elif('daymet' in options.met_dir.lower() and 'v1' in out_mettype.lower()):# from 1980-2010, half-degree
            fdirheader = options.met_dir+'/GSWP3v1_Daymet_'+varname+'_'
        elif('daymet' in options.met_dir.lower()):# from 1980-2014, half-degree
            fdirheader = options.met_dir+'/GSWP3_Daymet_'+varname+'_'
        ncfilein_cplbypass = sorted(glob.glob("%s*.*" % fdirheader))
        ncfilein_cplbypass = ncfilein_cplbypass[0]
        
    elif 'site' in out_mettype.lower():
        fdir = options.met_dir+'/'
        # So 'metdir' must be full path, e.g. ../atm/datm7/CLM1PT_data/1x1pt_US-Brw
        ncfilein_cplbypass=fdir+'all_hourly.nc'
    
    # new met nc files to create or write
    
    # (1) if in original output format
    if NCOUT_ORIGINAL:
        if (options.met_header==''):
            dirfiles = sorted(os.listdir(options.met_dir))
        else:
            fdirheader=options.met_dir+'/'+options.met_header.strip()
            dirfiles = sorted(glob.glob("%s*.*" % fdirheader))  # maybe file pattern in 'fileheader'
        if len(dirfiles)<=0: 
            print('Error: No template met nc file FOUND to create new one as following - ')
            print(fdirheader)
            sys.exit(-1)
        elif (os.path.isfile(dirfiles[0]) and str(dirfiles[0]).endswith('.nc')):
            ncfilein = dirfiles[0]
        else:
            print('Error: No template met nc file FOUND to create new one as following - ')
            print(dirfiles[0])
            sys.exit(-2)
            
        
        yyyymm=ncfilein.split('/')[-1]
        yyyymm=yyyymm.split('.')[-2]
        # elm met file usually ending like '*.1980-01.nc', while in CPL_BYPASS format, ending like '*_1901-2014_z??.nc'
        if ('GSWP3' in options.met_type): yyyymm = yyyymm.split('_')[-2]
        yyyy=yyyymm.split('-')[-2]
        mm=yyyymm.split('-')[-1]  # for 'GSWP3', 'mm' is the end-year while 'yyyy' is the starting-year
        
        mdoy=[0,31,59,90,120,151,181,212,243,273,304,334,365]#monthly starting DOY
        tyr = np.asarray(np.floor(t_jointed/365.0))
        tdoy = t_jointed-tyr*365.0
        
        tyrly = np.asarray(np.sort(np.unique(tyr)))
        
        for iyr in tyrly:
            print('YEAR: ', int(iyr), varname)
            for imm in range(len(mdoy)-1):
                #
                idx = np.where((tdoy>=mdoy[imm]) & (tdoy<mdoy[imm+1]) & (tyr==iyr))
                
                t=t_jointed[idx]
                varvals = sdata_jointed[idx]
                #
                ncfileout = ncfilein.split('/')[-1]
                ncfileout = ncfileout.replace(str(yyyy)+'-',str(int(iyr)).zfill(4)+'-') #tip: prefix '-' to prevent yyyy-mm messing-up
                ncfileout = ncfileout.replace('-'+str(mm),'-'+str(int(imm)+1).zfill(2))
                #
                if(options.nc_create):
                    # no-expanding for any dimension but re-size
                    nfmod.dupexpand(ncfilein, ncfileout, dim_name=tvarname, dim_len=len(t))
                    
                    
                    # time
                    try:
                        tunit = Dataset(ncfilein).variables[tvarname].getncattr('units')
                        t0=str(tunit.lower()).strip('days since')
                        t0=datetime.strptime(t0,'%Y-%m-%d %X')
                        t=t-(iyr*365+mdoy[imm])
                        tunit = tunit.replace(str(t0.year).zfill(4)+'-', str(int(iyr)).zfill(4)+'-')
                        tunit = tunit.replace('-'+str(t0.month).zfill(2)+'-', '-'+str(int(imm)+1).zfill(2)+'-')
                        if(tunit.endswith(' 00') and not tunit.endswith(' 00:00:00')):
                            tunit=tunit+':00:00'
    
                    except Exception as e:
                        print(e)
                        tunit = tunit_elm
                    error=nfmod.putvar(ncfileout, [tvarname], t, varatts=tvarname+'::units='+tunit)
                    # varatts must in format: 'varname::att=att_val; varname::att=att_val; ...'
                    if error!=0: sys.exit('nfmod.putvar WRONG')
                
                elif(options.nc_write):
                    if (not os.path.isfile(ncfileout)): 
                        print('Error: File '+ncfileout+ ' NOT exists, so CANNOT write data into.')
                        sys.exit(-3)
                #
                if (varname=='PRECTmms'): 
                    varvals=varvals/86400.0 # mm/day -> mm/s
                
                error=nfmod.putvar(ncfileout, [varname], varvals)
                if error!=0: sys.exit('nfmod.putvar WRONG')
            # end of for imm
        #end fo for iyr
    # end of if (NCOUT_ORIGINAL)
    
    #(2) in cplbypass format (NOTE: can output both original and cplbypass)
    if (NCOUT_CPLBYPASS): print('cpl_bypass template file: '+ ncfilein_cplbypass)
    if (NCOUT_CPLBYPASS and ncfilein_cplbypass!=''):
        if 'site' in out_mettype.lower()  or 'GSWP3' in out_mettype: 
            if 'site' in out_mettype.lower():
                ncfileout_cplbypass='all_hourly.nc'
            elif 'GSWP3' in out_mettype:
                ncfileout_cplbypass=ncfilein_cplbypass.split('/')[-1]
            else:
                print('currently only support 2 types of cpl_bypass: Site or GSWP3*')
                sys.exit(-1)
                
            tvarname = 'DTIME'
            
            if (options.nc_create):
                nfmod.dupexpand(ncfilein_cplbypass, ncfileout_cplbypass, dim_name=tvarname, dim_len=len(t_jointed))
                # ONLY create new nc ONCE
                
                
                # time
                try:
                    tunit = Dataset(ncfilein_cplbypass).variables[tvarname].getncattr('units')
                    t0=str(tunit.lower()).strip('days since')
                    if (t0.endswith(' 00')):
                        t0=t0+':00:00'
                    elif(t0.endswith(' 00:00')):
                        t0=t0+':00'
                    t0=datetime.strptime(t0,'%Y-%m-%d %X') # time must be in format '00:00:00'
                    yr0=np.floor(t_jointed[0]/365.0)
                    t=t_jointed - yr0*365.0
                    tunit = tunit.replace(str(t0.year).zfill(4)+'-', str(int(yr0)).zfill(4)+'-')
                    tunit = tunit.replace('-'+str(t0.month).zfill(2)+'-', '-01-')
                    if(tunit.endswith(' 00') and not tunit.endswith(' 00:00:00')):
                        tunit=tunit+':00:00'
                except Exception as e:
                    print(e)
                    tunit = tunit_elm
                error=nfmod.putvar(ncfileout_cplbypass, [tvarname], t, varatts=tvarname+'::units='+tunit)
                # varatts must in format: 'varname::att=att_val; varname::att=att_val; ...'
                if error!=0: sys.exit('nfmod.putvar WRONG')
                
                error=nfmod.putvar(ncfileout_cplbypass,['LONGXY'], LONGXY)
                error=nfmod.putvar(ncfileout_cplbypass,['LATIXY'], LATIXY)
                if 'site' in out_mettype.lower():
                    error=nfmod.putvar(ncfileout_cplbypass,['start_year'], np.floor(t_jointed[0]/365.0))
                    error=nfmod.putvar(ncfileout_cplbypass,['end_year'], np.floor(t_jointed[-1]/365.0))
            #
            elif(options.nc_write):
                if (not os.path.isfile(ncfileout_cplbypass)): 
                    print('Error: File '+ncfileout_cplbypass+ ' NOT exists, so CANNOT write data into.')
                    sys.exit(-3)
            
            if (varname=='PRECTmms' and varunit=='mm/d'): sdata_jointed=sdata_jointed/86400.0 # mm/d -> mm/s
            # scaling data as initeger
            if (varname=='PRECTmms'):
                data_ranges = [-0.04, 0.04]
            elif (varname=='FSDS'):
                data_ranges = [-20.0, 2000.0]
            elif (varname=='TBOT'):
                data_ranges = [175.0, 350.0]
            elif (varname=='RH'):
                data_ranges = [0.0, 100.0]
            elif (varname=='QBOT'):
                data_ranges = [0.0, 0.10]
            elif (varname=='FLDS'):
                data_ranges = [0.0, 1000.0]
            elif (varname=='PSRF'):
                data_ranges = [20000.0, 120000.0]
            elif (varname=='WIND'):
                data_ranges = [-1.0, 100.0]
            
            add_offset = (data_ranges[1]+data_ranges[0])/2.0
            scale_factor = (data_ranges[1]-data_ranges[0])*1.1/(2**15)
            
            # TIP: when data written, the input valule is in unpacked and the python nc4 program will do packing
            # the following line IS WRONG
            # varvals = (sdata_jointed-add_offset)/scale_factor # this IS WRONG when written, i.e. NOT NEEDED absolutely
            varvals = np.reshape(sdata_jointed, (1,-1)) # in CPL_BYPASS forcing-data format, dimension in (gridcell, DTIME) 
            
            # varatts must in format: 'varname::att=att_val; varname::att=att_val; ...'
            varatts = varname+'::add_offset='+str(add_offset)+ \
                ';'+varname+'::scale_factor='+str(scale_factor)
            error=nfmod.putvar(ncfileout_cplbypass, [varname], varvals, varatts=varatts)
            #error=nfmod.putvar(ncfileout_cplbypass, [varname], varvals)
            if error!=0: sys.exit('nfmod.putvar WRONG-'+varname+'-'+ncfileout_cplbypass)
            
        
        else:
            print('TODO: CPL_BYPASS format other than "Site/GSWP3" not yet supported')
            sys.exit(-1)
    
    elif(NCOUT_CPLBYPASS):
        print('CPL_BYPASS format output required but cannot find a template netcdf file, such as: all_hourly.nc or GSWP3_TBOT_1901-2014_z14.nc')
        sys.exit(-1)
    


