#!/usr/bin/env python

import sys
import math
from optparse import OptionParser
import numpy as np
from copy import deepcopy
from math import nan
from types import SimpleNamespace

# customized modules
from src.pytools.metdata_processing import elm_metdata_read
from src.pytools.metdata_processing import met_utils
from src.pytools.metdata_processing.met_utils import vpsat_pa
from src.pytools.metdata_processing.elm_metdata_write import elm_metdata_write

# ---------------------------------------------------------------
#  Data sub-timing (down) precipitation (rate-type data), from up-integration (known)
def Site_RateRedistribution(sdata_src, sdata_up, rh=[]):
    # 'sdata_src' is providing sub-timely variation, 
    # while 'sdata_up' is rate integrated-up, which should be conserved after re-distribution
    # For 'precipitation', for example, if no sub-timely pattern (variation, e.g. no rain in a day), 
    #  it's still possible to re-distribute if relative humidity known
    
    if (len(rh) >0):
        if (len(rh) != len(sdata_src)):
            print('Error: redistributing data size different: ', len(sdata_src), len(rh))
            exit(-1)
    
    # it's required that len(sdata_src) is muliple of len(sdata_up)
    if(np.mod(len(sdata_src),len(sdata_up))!=0):
        print('Error: length of "sdata_src" not multiple of "sdata_up" - ', len(sdata_src), len(sdata_up))
        exit(-2)
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
        sdata_jointting = np.nanmean([src,target], axis=0)
    elif (data_method == 'ratio'):
        if (not np.any(src==0.0)):
            ratio = target/src# 
            sdata_jointting = sdata_jointting*ratio
        else:
            print('Error: ratio cannot calculated duo to 0.0 value in data ')
            sys.exit(-2)
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
    t2 = np.round(t2,8)  # for unkown reason, t2 has some rounding error after conversions above
    #if(abs(t2[1]-t2[0])<1.0):
    #    t2 = t2 + 1                  # convert to 1-based DOY, if sub-daily data
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
    t_seasonally=[]
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
    
    
def DataSimilarity(data1, data2, t=None, t_unit='', v_name='', v_unit='', plt=None):
    
    if len(data1.shape)>1:
        data1=np.reshape(data1,(data1.size))
        data2=np.reshape(data2,(data2.size))    
    
    # ME/RMSE
    me=np.nanmean(data1-data2)
    rmse=np.nanmean([x*x for x in (data1-data2)])
    rmse=np.sqrt(rmse)
    
    """
    fit_model = np.polyfit(data1, data2, 1, cov=True)
    slope = fit_model[0]
    intcp = fit_model[1]

    l11 = [np.nanmin(data1),np.nanmin(data2),np.nanmax(data1),np.nanmax(data2)]

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
    """

    # return indicators, if any
    return rmse, me

##*****************************************************************************


#------------------- MAIN Program -----------------------------------------------

parser = OptionParser()

# E3SM met data 
parser.add_option("--e3sm_mettype", dest="met_type", default="CRUJRA", \
                  help="e3sminput met data type (default = 'GSWP3', options: CRU, CRUJRA, GSWP3, GSWP3v1, GSWP3, Site, and *_daymet/cplbypass_*)")
parser.add_option("--e3sm_metdir", dest="met_dir", default="./", \
                  help="e3sminput met directory (default = ./, i.e., under current directory)")
parser.add_option("--e3sm_metheader", dest="met_header", default="", \
                  help="e3sminput directory (default = '', i.e., standard file name)")
parser.add_option("--e3sm_metvars", dest="met_vars", default="ALL", \
                  help = "variable name to be merged/jointed: ALL or one of 'TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND' ")
#
# other sources of data to be pull into merged E3SM data
parser.add_option("--user_mettype", dest="user_mettype", default="daymet", \
                  help="user met data type (default = 'daymet', options: CRU, GSWP3, GSWP3v1, GSWP3, Site, and *_daymet/cplbypass_*, ATS_h5, NCDC, daymet, None)")
parser.add_option("--user_metdir", dest="user_metdir", default="./", \
                  help="user-defined met directory (default = ./, i.e., under current directory)")
parser.add_option("--user_metfile", dest="user_metfile", default="", \
                  help="user-defined met file(s) under user_metdir (default = '', i.e., None)")
#
parser.add_option("--joint_adjusting", dest="joint_adj", default=False, \
                  help = "joint data by smoothing transition ", action="store_true")
parser.add_option("--nc_create", dest="nc_create", default=True, \
                  help = "output new nc file", action="store_true")
parser.add_option("--nc_write", dest="nc_write", default=False, \
                  help = "output to a existed nc file", action="store_true")
#
(options, args) = parser.parse_args()


#---------------------------------------------------------------------------

startdays = -9999
enddays = -9999

#-----------------------------------------------------------------------------------
if (options.met_vars == '' or options.met_vars.lower() == 'all'):
    print('No variable name specified by " --varname=???", all of following ')
    print('TBOT, PRECTmms, QBOT/RH, FSDS, FLDS, PSRF, WIND')
elif(not options.met_vars in ['TBOT', 'PRECT', 'QBOT', 'RH', 'FSDS', 'FLDS', 'estFLDS', 'PSRF', 'WIND']):
    print('NOT supported variable name, should be one of : ')
    print('TBOT, PRECTmms, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND')
    sys.exit(-1)


#--------------------------------------------------------------------------------------
# standard met variables for ELM
if(options.met_vars=='TBOT'):
    vname_elm = ['TBOT']
    vunit_elm = ['K']
elif(options.met_vars=='PRECTmms'):
    vname_elm = ['PRECTmms']
    vunit_elm = ['mm/s']
elif(options.met_vars=='QBOT'):
    vname_elm = ['QBOT']
    vunit_elm = ['kg/kg']
elif(options.met_vars=='RH'):
    vname_elm = ['RH']
    vunit_elm = ['-']
elif(options.met_vars=='FSDS'):
    vname_elm = ['FSDS']
    vunit_elm = ['W/m2']
elif(options.met_vars=='FLDS'):
    vname_elm = ['FLDS']
    vunit_elm = ['W/m2']
elif(options.met_vars=='estFLDS'):
    vname_elm = ['estFLDS']
    vunit_elm = ['W/m2']
elif(options.met_vars=='PSRF'):
    vname_elm = ['PSRF']
    vunit_elm = ['Pa']
elif(options.met_vars=='WIND'):
    vname_elm = ['WIND']
    vunit = ['m/s']
elif(options.met_vars == '' or options.met_vars.lower() == 'all'):
    vname_elm =['FSDS', 'PRECTmms', 'TBOT','QBOT', 'RH', 'FLDS', 'PSRF', 'WIND']    
    vunit_elm =['W/m2', 'mm/s', 'K','kg/kg', '-', 'W/m2', 'Pa', 'm/s']    
else:
    print('NOT supported variable name, should be one of : ')
    print('TBOT, PRECTmms, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND')
    sys.exit(-1)

if (('RH' in options.met_type and 'QBOT' not in vname_elm) or \
    ('QBOT' in options.met_type and 'RH' not in vname_elm) or \
    'estFLDS' in options.met_vars):
    if ('TBOT' not in vname_elm): 
        vname_elm = vname_elm.append('RH')  # needed for conversion
    if ('PSRF' not in vname_elm): 
        vname_elm = vname_elm.append('PSRF')  # needed for conversion


#--- (1) ELM known-format met data reading


t_elm = []
vardatas = []
# read-in metdata from CPL_BYPASS_FULL
if ('cplbypass' in options.met_type or 'cplbypass_site' in options.met_type):
    cplbypass_dir=options.met_dir
    cplbypass_fileheader=options.met_header
    if 'crujra' in options.met_type.lower():
        cplbypass_mettype='CRUJRA'
    
    vnames=vname_elm
    if 'cplbypass_site' in options.met_type and 'QBOT' in vname_elm:
        vnames=np.delete(np.array(vname_elm), np.argwhere(np.array(vname_elm)=='QBOT'))
    elif 'RH' in vname_elm: # for other, RH must be excluded
        vnames=np.delete(np.array(vname_elm), np.argwhere(np.array(vname_elm)=='RH'))
        
    # read in
    ix,iy, varsdims, vardatas = \
        elm_metdata_read.clm_metdata_cplbypass_read( \
                            cplbypass_dir, cplbypass_mettype, vnames, [-999], [-999])

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
    
    del vnames

#--------------------------------------------------------------------------------------
# read-in metdata from full met directory, except for CPL_BYPASS
# E3SM met data
elif (('CRUJRA' in options.met_type) and \
    ('cplbypass' not in options.met_type and 'cplbypass_site' not in options.met_type)):
    if (options.met_dir == './'):
        print('e3sm met. directory is the current')
    else:
        print('e3sm met.  directory: '+ options.met_dir)

    # read in (TODO)
"""
    ix,iy, varsdims, vardatas = \
        Modules_metdata.clm_metdata_read(options.met_dir, options.met_header, options.met_type, 
                                         met_domain='NONE', [-999], [-999],'')
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
"""

#------------------------------
# if read-in ELM met-data successfully, need to do some conversion of dataset
t_elm     = []
tunit_elm = []
sdata_elm = {}
if len(vardatas)>0:
    vars_list = list(vardatas.keys())
    
    t_elm = vardatas[tvarname] # days since 1901-01-01 00:00
    tunit_elm = vardatas['tunit']
    if 'days since' in tunit_elm:
        t0 = str(tunit_elm.lower()).strip('days since')
        t0 = t0.split(' ')[0]
        yr0 = t0.split('-')[0]
        mm0 = t0.split('-')[1]
        dd0 = t0.split('-')[2]

        tunit_elm = tunit_elm.replace(yr0,'0001')
        t_elm = t_elm + int(yr0)*365.0 # converted from 'days-since-1901-01-01-00:00' to days since 0001-01-01 00:00

    #
    #       
    # usually  either 'RH' or 'QBOT' is available, but if the other is required
    if ('RH' in vname_elm):
        if 'QBOT' in vars_list and 'RH' not in vars_list:
            qbot = vardatas['QBOT']
            idx_zeros = np.where(qbot<=0)
            if len(idx_zeros[0])>0:
                for i in range(len(idx_zeros[0])):
                    qbot[idx_zeros[0][i],idx_zeros[1][i]] = (qbot[idx_zeros[0][i],idx_zeros[1][i]-1]+ \
                         qbot[idx_zeros[0][i],idx_zeros[1][i]+1])/2.0 # need something better to gap-fill
                
            if 'TBOT' in vars_list:
                tk = vardatas['TBOT']
                if 'PSRF' in vars_list:
                    pres_pa = vardatas['PSRF']
                else:
                    pres_pa = 101325.0
                sdata_elm['RH'] = met_utils.convertHumidity(tk, pres_pa, q_kgkg=qbot)
            else:
                print('ERROR: for RH coverting from QBOT, air temperature is required')
                sys.exit(-1)
                
    if ('QBOT' in vname_elm):
        if 'RH' in vars_list and 'QBOT' not in vars_list:
            rh = vardatas['RH']
            if 'TBOT' in vars_list:
                tk = np.squeeze(vardatas['TBOT'])
                if 'PSRF' in vars_list:
                    pres_pa = vardatas['PSRF']
                else:
                    pres_pa = 101325.0
            else:
                print('ERROR: for RH coverting from QBOT, air temperature is required')
                sys.exit(-1)
            sdata_elm['QBOT'] = met_utils.convertHumidity(tk, pres_pa, rh_100=rh)
    
    # if FLDS to be estimated from humidity and temperature
    if ('FLDS' in vname_elm and 'estFLDS' in options.met_vars):
        #Longwave radiation (calculated from air temperature, humidity)
        if 'TBOT' in vars_list:
            tk = vardatas['TBOT']
            print('ERROR: for calculating FLDS, air temperature is required')
            sys.exit(-1)
        if 'PSRF' in vars_list:
            pres_pa = vardatas['PSRF']
        else:
            pres_pa = 101325.0
        if 'QBOT' in vars_list:
            qbot = vardatas['QBOT']
            rh = []
        elif 'RH' in vars_list:
            rh = vardatas['RH']
            qbot = []
        else:
            print('ERROR: for calculating FLDS, either RH or QBOT is required')
            sys.exit(-1)
        sdata_elm['FLDS'] = \
            met_utils.calcFLDS(tk, pres_pa, q_kgkg=qbot, rh_100=rh)
    
    for v in vname_elm:
        if len(sdata_elm)>0:
            if v in sdata_elm.keys(): continue   # if already in as above
        
        if v in vardatas.keys():
            sdata_elm[v] = vardatas[v]
        else:
            print('ERROR: data ', v, 'NOT  in diretory: ', options.met_dir)
            sys.exit(-1)
            
    
    # clean-up
    del vardatas, vars_list
    
#--------------------------------------------------------------------------------------
#--- (2) read-in metdata from daymet (daily)
if (options.user_mettype=='daymet'):
    user_metfilename = options.user_metdir+'/'+options.user_metfile
    
    if ('.csv' in options.user_metfile):    
        siteinfo,vardatas = \
            elm_metdata_read.singleDaymetReadCsvfile(user_metfilename)


    # 
    t = (vardatas['year']-1)*365.0 + vardatas['yday']   # this is the end of day, and a full-year=year-1
    t_user = np.asarray(t)
    tunit_user = 'days since 0001-01-01 00:00'

    if(options.met_vars=='TBOT'):
        vname_user = ['tmax (deg c)', 'tmin (deg c)']
        vunit_user = ['deg c', 'deg c']
    elif(options.met_vars=='PRECT'):
        vname_user = ['prcp (mm/day)']
        vunit_user = ['mm/day']
    elif(options.met_vars=='QBOT'):
        vname_user = ['vp (Pa)','tmax (deg c)', 'tmin (deg c)']  # will do conversion
        vunit_user = ['Pa','deg c', 'deg c']
    elif(options.met_vars=='RH'):
        vname_user = ['vp (Pa)','tmax (deg c)', 'tmin (deg c)']  # will do conversion
        vunit_user = ['Pa','deg c', 'deg c']
    elif(options.met_vars=='FSDS'):
        vname_user = ['srad (W/m^2)']
        vunit_user = ['W/m^2']
    elif(options.met_vars=='FLDS' or options.met_vars=='estFLDS'):
        vname_user = ['estFLDS']
        vunit_user = ['W/m^2']
    elif(options.met_vars=='PSRF'):
        vname_user = ['NONE']  # a specific situation
        vunit_user = ['Pa']
    elif(options.met_vars=='WIND'):
        vname_user = ['NONE']
        vunit_user = ['m/s']
    elif(options.met_vars.lower()=='all' or options.met_vars==''):
        vname_user = ['tmax (deg c)', 'tmin (deg c)', 'prcp (mm/day)','vp (Pa)','srad (W/m^2)']
        vunit_user = ['deg c', 'deg c', 'mm/day', 'Pa', 'W/m^2']
    else:
        print('NOT supported variable name, should be one of : ')
        print('TBOT, PRECT, QBOT, RH, FSDS, FLDS, estFLDS, PSRF, WIND')
        sys.exit(-1)
            
    sdata_user = {}
    for v in vname_user:
        if v in vardatas.keys():
            sdata_user[v] = vardatas[v]
        else:
            print('ERROR: data ', v, 'NOT  in diretory: ', user_metfilename)
            sys.exit(-1)
            
    #

    del vardatas, t
    #


#--------------------------------------------------------------------------------------
# Year/day Matching and cutoff
dt_elm = np.mean(np.diff(t_elm))
dt_user= np.mean(np.diff(t_user))
t_min = max(min(t_elm)-dt_elm, min(t_user)-dt_user) # assuming t_elm/t_user is end of timestep
t_max = min(max(t_elm), max(t_user)) 
tidx_elm  = np.where((t_elm>=t_min) & (t_elm<t_max))
tidx_user = np.where((t_user>=t_min) & (t_user<=t_max)) # test shows both ends must be inclusive

t_elm = t_elm[tidx_elm]
t_user= t_user[tidx_user]

ts_dly = 4
t2d_elm = t_elm.reshape(-1,ts_dly)   # (days, sub-days), for daily data processing
sdata2d_elm = {}
for v in sdata_elm.keys():
    sdata_elm[v] = sdata_elm[v][...,tidx_elm]
    sdata2d_elm[v] = sdata_elm[v].reshape(nxy,-1,ts_dly) #[grid, day, time]
for v in sdata_user.keys():
    sdata_user[v] = sdata_user[v][...,tidx_user]


#---------------------------------------------------------------------------------------------------------
#---- Sub-daily Downscaling
#
sdata_elm_adj = {}

#---------------------------------------------------------------------------------------------------------
#---- PSRF, WIND: no downscaling

sdata_elm_adj['PSRF'] = sdata_elm['PSRF']

sdata_elm_adj['WIND'] = sdata_elm['WIND']



#---------------------------------------------------------------------------------------------------------
#---- do TBOT first
varadj = 'TBOT'
if varadj in vname_elm:
    # daily max/min of elm data
    print("Adjusting sub-daily TBOT to match with Daymet daily Tmax/Tmin")
        
    td_max = np.max(sdata2d_elm[varadj],axis=2)
    idx_max= np.argmax(sdata2d_elm[varadj],axis=2)
    idx_day= np.asarray(range(sdata2d_elm[varadj].shape[1]))
    if td_max.shape == sdata_user['tmax (deg c)'].shape:
        vadj = (sdata_user['tmax (deg c)']+273.15)/td_max
        tadj = t_elm[idx_day*ts_dly+idx_max[0]]
    else:
        print('day numubers are NOT same: ', td_max.shape, sdata_user['tmax (deg c)'].shape)
        sys.exit(-1)

    td_min = np.min(sdata2d_elm[varadj],axis=2)
    idx_min= np.argmin(sdata2d_elm[varadj],axis=2)
    if td_min.shape == sdata_user['tmin (deg c)'].shape:
        vadj = np.hstack((vadj, (sdata_user['tmin (deg c)']+273.15)/td_min)) # will sort later
        tadj = np.hstack((tadj, t_elm[idx_day*ts_dly+idx_min[0]]))
    else:
        print('day numubers are NOT same: ', td_max.shape, sdata_user['tmin (deg c)'].shape)
        sys.exit(-1)
    
    idx_sorted = np.argsort(tadj)
    tadj = tadj[idx_sorted]
    vadj = vadj[...,idx_sorted]
    #linear interpolating vadj(tadj)
    sdata_elm_adj[varadj] = sdata_elm[varadj]*np.interp(t_elm, tadj, vadj[0])

    
#---------------------------------------------------------------------------------------------------------
#---- do QBOT/RH, with dependcy on TBOT, PSRF
varadj = 'QBOT'
if varadj not in vname_elm: varadj = 'RH'
if varadj in vname_elm:
    print("Adjusting sub-daily QBOT/RH to match with Daymet daily vp")
        
    # interpolation is done by QBOT, specific humidity, because RH interpolation may have math issue
    vd_mean = np.mean(sdata2d_elm[varadj],axis=2)
    if varadj=='RH':
        tk = np.mean(sdata2d_elm['TBOT'],axis=2)
        pres = np.mean(sdata2d_elm['PSRF'], axis=2) 
        vd_mean = met_utils.convertHumidity(tk, pres, rh_100=vd_mean) # rh -> qbot in kg/kg
    idx_mean= np.argmax(sdata2d_elm[varadj],axis=2)+np.argmax(sdata2d_elm[varadj],axis=2)
    idx_mean= np.int32(idx_mean/2)
    idx_day= np.asarray(range(sdata2d_elm[varadj].shape[1]))
    if vd_mean.shape == sdata_user['vp (Pa)'].shape:
        tk = (sdata_user['tmax (deg c)']+sdata_user['tmin (deg c)'])/2.0+273.15         # daily Tair in K
        pres = np.mean(sdata2d_elm['PSRF'], axis=2)                           # no daymet air pressure, using elm one but at daily
        vadj = sdata_user['vp (Pa)']/vpsat_pa(tk)*100.0   # rh = vp/vp_sat
        vadj = met_utils.convertHumidity(tk, pres, rh_100=vadj) # rh -> qbot in kg/kg
        vadj = vadj/vd_mean
        tadj = t_elm[idx_day*ts_dly+idx_mean[0]]
    else:
        print('day numubers are NOT same: ', vd_mean.shape, sdata_user['vp (Pa)'].shape)
        sys.exit(-1)

    
    idx_sorted = np.argsort(tadj)
    tadj = tadj[idx_sorted]
    vadj = vadj[...,idx_sorted]
    #linear interpolating vadj(tadj)
    sdata_elm_adj[varadj] = sdata_elm[varadj]*np.interp(t_elm, tadj, vadj[0]) #QBOT even if varadj=='RH'
    if varadj=='RH' and 'TBOT' in sdata_elm.keys():
        tk = sdata_elm['TBOT']
        pres_pa = sdata_elm['PSRF']
        rh = sdata_elm_adj['RH']   # from above, but which actually is QBOT
        sdata_elm_adj[varadj] = \
            met_utils.convertHumidity(tk, pres_pa, rh_100=rh)
    

#---------------------------------------------------------------------------------------------------------
#---- do PRECTmms, with dependcy on RH
varadj = 'PRECTmms'

if varadj in vname_elm:
    print("Adjusting sub-daily PRECTmms to match with Daymet daily PRCT")
        
    #
    if 'RH' not in sdata_elm_adj.keys():
        if 'TBOT' in sdata_elm_adj.keys() and 'QBOT' in sdata_elm_adj.keys() \
            and 'PSRF' in sdata_elm_adj.keys():
            
            tk = sdata_elm_adj['TBOT']
            pres_pa = sdata_elm_adj['PSRF']
            qbot = sdata_elm_adj['QBOT']
            sdata_elm_adj['RH'] = \
                met_utils.convertHumidity(tk, pres_pa, q_kgkg=qbot)        
        else:
            print("requiring RH for redistributing rate")
            sys.exit(-1)
    
    prct = Site_RateRedistribution(sdata_elm[varadj][0,0,:], \
                                    sdata_user['prcp (mm/day)'][0,:], \
                                    rh=sdata_elm_adj['RH'][0,0,:])

    sdata_elm_adj[varadj] = prct.reshape(1,1,-1)/86400.0  #mm/day to mm/s


#---------------------------------------------------------------------------------------------------------
#---- do FSDS
varadj = 'FSDS'

if varadj in vname_elm:
    print("Adjusting sub-daily FSDS to match with Daymet daily srad")
        
    # only adjusting NON-zeros (zeros implies night-time)
    idx_zeros = np.where(sdata2d_elm[varadj]<=0.0)
    v_nonzeros = sdata2d_elm[varadj]
    v_nonzeros[idx_zeros] = nan
    vd_mean = np.nanmean(v_nonzeros,axis=2)

    if vd_mean.shape == sdata_user['srad (W/m^2)'].shape:
        vadj = sdata_user['srad (W/m^2)']/vd_mean
    else:
        print('day numubers are NOT same: ', td_max.shape, sdata_user['tmax (deg c)'].shape)
        sys.exit(-1)
    
    #simple daily factor adjusting 
    for it in range(ts_dly): 
        sdata2d_elm[varadj][...,it] = sdata2d_elm[varadj][...,it]*vadj
    sdata_elm_adj[varadj] = sdata2d_elm[varadj].reshape(sdata_elm[varadj].shape)


#---------------------------------------------------------------------------------------------------------
#---- do FLDS
varadj = 'FLDS'

if varadj in vname_elm:
    print("Estimating sub-daily FLDS")
        
    #
    if 'TBOT' in sdata_elm_adj.keys() and 'QBOT' in sdata_elm_adj.keys() \
        and 'PSRF' in sdata_elm_adj.keys():

        tk = sdata_elm_adj['TBOT']
        pres = sdata_elm_adj['PSRF']
        qbot = sdata_elm_adj['QBOT']
        rh = sdata_elm_adj['RH']

        sdata_elm_adj['FLDS'] = \
            met_utils.calcFLDS(tk, pres_pa, q_kgkg=qbot, rh_100=rh)
        
        del tk, pres, qbot, rh

    else:
        print("requiring TBOT, QBOT, PSRF for estimating FLDS")
        sys.exit(-1)

 
#--------------------------------------------------------------------------------------
#
#------------------- (3) Save into ELM offline forcing data ---------------------------
#                    
"""
options.nc_create=False
options.nc_write=False
    
"""
write_options = SimpleNamespace( \
            met_idir = options.met_dir, \
            nc_create = options.nc_create, \
            nc_write = options.nc_write, \
            nc_write_mettype = options.met_type+'_daymet' )

sdata_elm_adj['LONGXY'] = LONGXY
sdata_elm_adj['LATIXY'] = LATIXY
sdata_elm_adj['time']   = t_elm
sdata_elm_adj['tunit']  = tunit_elm

elm_metdata_write(write_options, sdata_elm_adj)


