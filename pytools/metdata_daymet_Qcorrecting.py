#!/usr/bin/env python

import os, sys
from datetime import datetime
import glob
import re
import math
from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset
from copy import deepcopy
from scipy import interpolate

# customized modules
import Modules_metdata
#from Modules_metdata import clm_metdata_cplbypass_read    # CPL_BYPASS
#from Modules_metdata import clm_metdata_read              # GSWP3
#from Modules_metdata import singleNCDCReadCsvfile         # NCDC 
#from Modules_metdata import subsetDaymetRead1NCfile       # DAYMET in nc4 format
#from Modules_metdata import singleDaymetReadCsvfile       # DAYMET in csv format
from hdf5_modules import Read1hdf                         # for ATS metdata in h5 format
import netcdf_modules as nfmod

# ---------------------------------------------------------------
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

#---------------------------------

parser = OptionParser()

parser.add_option("--e3sm_metdomain", dest="met_domain", default="", \
                  help="e3sminput met domain file (default = a file in current met_dir")
# E3SM GSWP3-daymet data
parser.add_option("--mettype", dest="met_type", default="cplbypass_GSWP3_daymet4", \
                  help="e3sminput met data type (default = 'cplbypass_GSWP3_daymet4', options:  cplbypass_GSWP3v1_Daymet, cplbypass_GSWP3_Daymet)")
parser.add_option("--e3sm_metdir", dest="met_dir", default="./", \
                  help="e3sminput met directory (default = ./, i.e., under current directory)")
parser.add_option("--e3sm_metheader", dest="met_header", default="", \
                  help="e3sminput directory (default = '', i.e., standard file name)")
parser.add_option("--varname", dest="vars", default="QBOT", \
                  help = "variable name to be merged/jointed: ONLY one of 'TBOT, PRECT, QBOT, RH, FSDS, FLDS, PSRF, WIND' ")
# other sources of data to be pull into merged E3SM data
parser.add_option("--user_mettype", dest="user_mettype", default="cplbypass_GSWP3", \
                  help="e3sminput met data type (default = 'cplbypass_GSWP3', options: cplbypass_GSWP3, GSWP3, GSWP3v1)")
parser.add_option("--user_metdir", dest="user_metdir", default="./", \
                  help="user-defined met directory (default = ./, i.e., under current directory)")
parser.add_option("--user_metfile", dest="user_metfile", default="", \
                  help="user-defined met file(s) under user_metdir (default = '', i.e., standard file name)")
#
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
elif(not options.vars in ['TBOT', 'PRECT', 'QBOT', 'RH', 'FSDS', 'FLDS', 'PSRF', 'WIND']):
    print('NOT supported variable name, should be one of : ')
    print('TBOT, PRECT, QBOT, RH, FSDS, FLDS, PSRF, WIND')
    sys.exit(-1)


#--------------------------------------------------------------------------------------
# 
# v1 is from 1980-2010, otherwise from 1980-2014
#_Daymet is in half-degree, and _daymet4 is in 1km res, and NA only
elmmettype = ['cplbypass_GSWP3', 'cplbypass_GSWP3v1', \
            'cplbypass_GSWP3_daymet4', \
            'cplbypass_GSWP3v1_Daymet','cplbypass_GSWP3_Daymet', \
            'GSWP3', 'GSWP3v1']

met_domain=options.met_domain

metdir=options.met_dir
metfileheader=options.met_header
mettype=options.met_type

vname_elm = options.vars
vnames=[vname_elm]
if (vname_elm=='QBOT' or vname_elm=='RH'):
    vnames=['QBOT','TBOT','PSRF']

# if for merging mets are both ELM met-type
if options.user_mettype in elmmettype:
    metdir=[options.met_dir, options.user_metdir]
    if(options.met_header=='' and options.user_metfile==''):
        metfileheader=['./', './']
    else:
        metfileheader=[options.met_header, options.user_metfile]
    mettype=[options.met_type, options.user_mettype]
    
    
for imet in range(len(mettype)):
    
    imet_type = mettype[imet]
    if imet_type not in elmmettype:
        print(imet_type, 'NOT yet supported! Must be one of followint -')
        print(elmmettype)
        sys.exit(-1)
    imet_dir = metdir[imet]
    imet_header = metfileheader[imet]
    t = []
    vardatas = []
    
    if 'GSWP3' in imet_type and 'daymet' not in imet_type.lower():
        # this must be the second met-type options
        dx=0.5
        dy=0.5
        lon=[-999]
        lat=[-999]
        # as long as Daymet-1km-cells' centers known
        if isinstance(LONGXY, np.ndarray): lon=LONGXY
        if isinstance(LATIXY, np.ndarray): lat=LATIXY
    else:
        dx=-999.9
        dy=-999.9
        lon=[-999]
        lat=[-999]
        
    
    # read-in metdata from CPL_BYPASS_FULL
    if ('cplbypass_' in imet_type):
        cplbypass_dir=imet_dir
        cplbypass_mettype='GSWP3'
        if ('daymet' in imet_type.lower()): cplbypass_mettype='GSWP3_daymet'
        if ('v1' in imet_type): cplbypass_mettype='GSWP3v1'
        if ('v1' in imet_type and 'daymet' in imet_type.lower()): cplbypass_mettype='GSWP3v1_daymet'
        
        
        # read in
        zones,zlines, varsdims, vardatas = \
            Modules_metdata.clm_metdata_cplbypass_read( \
                                cplbypass_dir, cplbypass_mettype, vnames, lons=lon, lats=lat)
        
    
        # 
        nx=1
        ny=0
        for iz in zlines.keys():
            if np.isscalar(zlines[iz]):
                ny=ny+1
            else:
                ny=ny+len(zlines[iz])
        nxy=nx*ny
    
        tvarname = 'DTIME'
        LATIXY = vardatas['LATIXY']
        LONGXY = vardatas['LONGXY']
        #--------------------------------------------------------------------------------------
    # read-in metdata from full met directory, except for CPL_BYPASS
    # E3SM met data
    elif ('GSWP3' in imet_type and \
           ('cplbypass_' not in imet_type)):
        
        # GSWP3 data QBOT needs redo, because its RH is NOT freezing adjusted
        # (When calculating RH from QBOT, RH often over 100 in freezing winter)
        if (vname_elm=='QBOT'):
            vnames=['QBOT','TBOT','PSRF']
        
        # read in
        #
        ix,iy, varsdims, vardatas = \
            Modules_metdata.clm_metdata_read(imet_dir, imet_header, imet_type, met_domain, lon, lat,vnames)
        # 
        nx=1
        ny=len(iy)
        nxy=nx*ny

        tvarname = 'time'    # variable name for time/timing
        #vars =  ['time','tunit','LONGXY','LATIXY',
        #            'FLDS','FSDS','PRECTmms','PSRF','QBOT','TBOT','WIND']
        LATIXY = vardatas['LATIXY']
        LONGXY = vardatas['LONGXY']

    #--------------------------------------------------------------------------------------
    
    #------------------------------
    # if read-in data successfully
    if len(vardatas)>0:
        vars_list = list(vardatas.keys())
        
        t_elm = vardatas[tvarname]
        tunit_elm = vardatas['tunit']
        t0 = 0.0
        if ('1901-01-01' in tunit_elm):
            tunit_elm = tunit_elm.replace('1901-','0001-')
            t0 = 1901*365 # converted from 'days-since-1901-01-01-00:00' to days since 0001-01-01 00:00
            t_elm = np.asarray(t_elm)+t0
        tvarname_elm = tvarname
       # 
        if (vname_elm=='QBOT' or vname_elm=='RH'):
            #
            if 'TBOT' in vars_list:
                tk = np.squeeze(vardatas['TBOT'])
                if 'PSRF' in vars_list:
                    pres_pa = np.squeeze(vardatas['PSRF'])
                else:
                    print('ERROR: for RH coverting from QBOT, air presssure is required')
                    sys.exit(-1)
            else:
                print('ERROR: for RH coverting from QBOT, air temperature is required')
                sys.exit(-1)

            sdata_elm = np.squeeze(vardatas['QBOT'])
            
            if ('GSWP3' in imet_type) and ('daymet' not in imet_type):
                #correction for GSWP3 QBOT data, 
                # especially during winter freezing period, RH from original QBOT often over 100%
                # a possible reason: RH-->QBOT using 1 single QSAT function for freezing air
                sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=sdata_elm,vpsat_frz=False) # original RH
                if (vname_elm=='QBOT'):
                    sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, rh_100=sdata_elm)  # re-cal. QBOT for GSWP3 dataset
            elif(vname_elm=='RH'):
                sdata_elm = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=sdata_elm,vpsat_frz=True)
        
        
        else:
            sdata_elm = deepcopy(vardatas[varname])
            sdata_elm = np.squeeze(sdata_elm)

        
        # clean up large data in memory
        del vardatas
        
        #
        idx=np.where(LONGXY<=0)
        if (len(idx[0])>0): LONGXY[idx]=360.0+LONGXY[idx]
        #
        
        if len(sdata_elm.shape)==1:
            sdata_elm = np.reshape(sdata_elm, (1,-1)) # 1D --> 2D (n,DTIME/time)
            if(vname_elm=='QBOT' or vname_elm=='RH'): 
                tk = np.reshape(tk, (1,-1))
                pres_pa = np.reshape(pres_pa, (1,-1))
        elif imet_type=='GSWP3' or imet_type=='GSWP3v1':
            sdata_elm = np.swapaxes(sdata_elm, 0, 1) # data in (time,n) --> (n,time)
            
        if len(imet_type)>1:
            if imet == 0:
                lon_elm=LONGXY
                lat_elm=LATIXY
                nxy_elm=nxy
                zones_elm=deepcopy(zones)
                zlines_elm=deepcopy(zlines)
                is2d_elm = True
                if(nx==1 or ny==1):is2d_elm=False
                
                t_elm1 = deepcopy(t_elm)
                tunit_elm1 = deepcopy(tunit_elm)
                tvarname_elm1 = tvarname_elm
                vname_elm1 = vname_elm
                sdata_elm1 = deepcopy(sdata_elm)
                if(vname_elm=='QBOT' or vname_elm=='RH'):
                    pres_elm = deepcopy(pres_pa)
                    tk_elm = deepcopy(tk)
            
            elif imet == 1: # if this is the case, the following user_mettype won't read
                if (len(t_elm1))>0:
                    dt=max(np.diff(t_elm1))
                    idx=np.where((t_elm>=t_elm1[0]-dt/2.0) & (t_elm<=t_elm1[-1]+dt/2.0))
                    if (len(idx[0])>0):
                        sdata_elm=sdata_elm[:,idx[0]]
                        t_elm=t_elm[idx[0]]
                        
                if(len(t_elm)>0):
                    dt=max(np.diff(t_elm))
                    idx=np.where((t_elm1>=t_elm[0]-dt/2.0) & (t_elm1<=t_elm[-1]+dt/2.0))
                    if (len(idx[0])>0):
                        t_elm1=t_elm1[idx[0]]
                        sdata_elm1=sdata_elm1[:,idx[0]]
                        if(vname_elm=='QBOT' or vname_elm=='RH'): 
                            pres_elm = pres_elm[:,idx[0]]
                            tk_elm = tk_elm[:,idx[0]]
                
                
                lon_user=LONGXY
                lat_user=LATIXY
                if('cplbypass' in options.user_mettype):
                    zones_user=deepcopy(zones)
                    zlines_user=deepcopy(zlines)
                nxy_user=nxy
                is2d_user = True
                if(nx==1 or ny==1):is2d_user=False
                t_user = deepcopy(t_elm)
                tunit_user = deepcopy(tunit_elm)
                tvarname_user = tvarname_elm
                vname_user = vname_elm
                sdata_user = deepcopy(sdata_elm)
                
                #copy back elm1
                t_elm = deepcopy(t_elm1)
                tunit_elm = deepcopy(tunit_elm1)
                tvarname_elm = tvarname_elm1
                vname_elm = vname_elm1
                sdata_elm = deepcopy(sdata_elm1)
                
                
                del t_elm1, tunit_elm1, tvarname_elm1, vname_elm1, sdata_elm1, tk, pres_pa
                del lon, lat, zones, zlines
            
#--------------------------------------------------------------------------------------
#look-up half-degree grid indx for daymet4 1km cells

if (dx>0.0 and dy>0.0 and (not is2d_user)): 
    gidx = []
    cidx_elm = {}
    g=-1
    sum_cidx = 0
    for ig in range(nxy_user):
        idx=np.where((lon_elm>lon_user[ig]-dx/2.0) & (lon_elm<=lon_user[ig]+dx/2.0) & \
                     (lat_elm>lat_user[ig]-dy/2.0) & (lat_elm<=lat_user[ig]+dy/2.0))
        
        if len(idx[0])>0:
            g=g+1
            gidx.append(int(ig))  # half-degree grid index
            cidx_elm[g]=idx       # 1km cell index, which within the half-degree grid idx 'ig'
        
        sum_cidx = sum_cidx + len(idx[0])
        if (ig==nxy_user-1) and sum_cidx!=len(lon_elm):
            print('Inconsistent numbers in dividing daymet cells into GSWP3 grid: ')
            print('summed cells: ', sum_cidx, 'MUST be equal: ', len(lon_elm))
else:
    print('NOT yet supported user met-type:', options.user_mettype)
    sys.exit(-1)

#---------------------------------------------------------------------------------------------------------
ratio_elm = np.ones_like(sdata_elm)  # adjusting ratio (multiplier) for sdata_elm
for i in range(len(gidx)):
    #
    ig = gidx[i]
    if (options.user_mettype !=''): # assuming this data covering half-degree range
        t1 = t_user
        sdata1 = sdata_user[ig,:]
    else:
        print('Must have user met-data')
        sys.exit(-1)
    
    ic = cidx_elm[i][0]
    if options.met_type in elmmettype:
        sdata2 = sdata_elm[ic,:]
        t2 = t_elm
        sdata2_tk = tk_elm[ic,:]
        sdata2 = np.nanmean(sdata2, axis=0)  # averaging over data within half-degree range
        sdata2_tk = np.nanmean(sdata2_tk, axis=0)
    else:
        print('Must have met-data')
        sys.exit(-1)
        
    
    #--------------------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------------------
    # Merging
    dts1 = np.max(np.diff(t1))
    dts2 = np.max(np.diff(t2))
    
    # by jointing 2 datasets directly, from data1-->data2
    T_SMOOTHING = False
    if (T_SMOOTHING):
        
        # to be jointed/temporally-scaling (pulling data1 into data2)
        t_jointed = deepcopy(t2[t2<t1[0]])  # probably empty initially
        sdata_jointed = deepcopy(sdata2[:,[t2<t1[0]][0]])
        
        # loop through 't2'
        for it in range(len(t2)):
            
            td0_t1 = t1[it]
            td0_t2 = t2[it]
            t_jointting = t2[it]
            tidx = np.where((t1>=td0_t2-dts2/2.0) & (t1<td0_t2+dts2/2.0))
            if len(tidx[0])>0:
                it1=tidx[0][0]  #usually only 1, but just in case
                sdata=sdata1[:,it:it+1]
                sdata=np.asarray(sdata)
                sdata_jointting = np.nanmean([sdata2[:,it:it+1],sdata], axis=0)
                sdata_jointting[sdata<=0.0] = sdata[sdata<=0.0]
                # very specific to correct daymet-derived QBOT
                if((options.vars=='QBOT' or options.vars=='RH') \
                    and 'daymet' in options.met_type.lower()):
                        idx=np.where(sdata2_tk[:,it:it+1]>273.15) # i.e. if air NOT freezing, pick the merging data, otherwise averaging as above
                        if len(idx[0])>0: 
                            sdata_jointting[idx]=sdata[idx]
            else:
                sdata_jointting = sdata2[:,it:it+1]
            
            t_jointed = np.hstack((t_jointed, t_jointting))
            sdata_jointed = np.hstack((sdata_jointed, sdata_jointting))
        
        #DONE 'for it in range(t2)'
        
    else: # no jointing of datasets
        t_jointed = deepcopy(t2)
        idx1=np.where(t1<=t2[-1]+dt/2.0)[0]
        sdata=sdata1[...,idx1]
        sdata_jointed = np.nanmean([sdata2,sdata], axis=0)
        sdata_jointed[sdata<=0.0] = sdata[sdata<=0.0]  # just in case some negative values as filling
        # very specific to correct daymet-derived QBOT
        if((options.vars=='QBOT' or options.vars=='RH') \
            and 'daymet' in options.met_type.lower()):
                idx=np.where(sdata2_tk>273.15) # i.e. if air NOT freezing, pick the merging data, otherwise averaging as above
                if len(idx[0])>0: sdata_jointed[idx]=sdata1[idx]
        
    # end of if
    #
    
    r=sdata_jointed/sdata2
    idx=np.where((sdata2<=0.0) | (sdata_jointed<=0.0))
    if (len(idx[0])>0): r[idx]=1.0
    ratio_elm[ic,...] = r
    print('checking adjusting ratio:',i, ig, np.nanmean(ratio_elm[ic,...]), np.nanmean(ratio_elm), r.shape)
    
# end of for i in range(len(gidx))


# by this point, ratio_elm shall be filled values fully
sdata_jointed = ratio_elm*sdata_elm

#-------------------------------------------------------------------------------------    
# some simple checkings
if (options.vars=='RH' and 'Site' not in options.nc_write_mettype):
    # var name has modified for non-Site met-type
    options.vars='QBOT'
    vname_elm = 'QBOT'
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
        
    
    

    
#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
# save in ELM forcing data format
if (options.nc_create or options.nc_write):
    met_type = options.nc_write_mettype
    if met_type =='': met_type = options.met_type
    if 'cplbypass_' in options.nc_write_mettype:
        NCOUT_CPLBYPASS=True
    else:
        NCOUT_CPLBYPASS=False
    
    if (options.nc_create and options.nc_write):
        print('Error: cannot have both "--nc_create" and "--nc_write" ')
        sys.exit(-1)
    elif (options.nc_create):
        print('Create new ELM forcing data ? ', options.nc_create)
    elif (options.nc_write):
        print('Write to existed ELM forcing data ? ', options.nc_write)
    
    # get a template ELM forcing data nc file
    # 
    varname = options.vars
    ncfilein_cplbypass = ''
    metdir = options.nc_write_stdmetdir
    if metdir=='': metdir=options.met_idir
    if 'GSWP3' in met_type and 'cplbypass' in met_type:
        
        for iz in range(len(zones_elm)):
            z=zones_elm[iz]
            if('v1' in met_type):
                ncfilein_cplbypass=metdir+'/GSWP3_'+varname+'_1901-2010_z'+str(int(z)).zfill(2)+'.nc'
            elif('daymet' in met_type):
                ncfilein_cplbypass=metdir+'/GSWP3_daymet4_'+varname+'_1980-2014_z'+str(int(z)).zfill(2)+'.nc'
            else:
                ncfilein_cplbypass=metdir+'/GSWP3_'+varname+'_1901-2014_z'+str(int(z)).zfill(2)+'.nc'
            
            # new met nc files to create or write
            
            if (NCOUT_CPLBYPASS and ncfilein_cplbypass!=''):
                print('cpl_bypass template file: '+ ncfilein_cplbypass)
                ncfileout_cplbypass=options.nc_write_stdmetdir+'/'+ncfilein_cplbypass.split('/')[-1]
                
                tvarname = 'DTIME'
            
                if options.nc_write:
                    os.system('cp '+ncfilein_cplbypass+' '+ncfilein_cplbypass+'-orig')
                    if not os.path.isfile(ncfileout_cplbypass): ncfileout_cplbypass = ncfilein_cplbypass
                elif (options.nc_create):
                    if(os.path.abspath(ncfilein_cplbypass)==os.path.abspath(ncfileout_cplbypass)):
                        os.system('cp '+ncfilein_cplbypass+' '+ncfilein_cplbypass+'-orig')
                    #
                    nfmod.dupexpand(ncfilein_cplbypass, ncfileout_cplbypass, dim_name=tvarname, dim_len=len(t_jointed))
                    # ONLY create new nc ONCE
                    options.nc_create = False
                    
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
                    if error!=0: sys.exit('nfmod.putvar-varatts WRONG')
                    
                    error=nfmod.putvar(ncfileout_cplbypass,['LONGXY'], LONGXY)
                    error=nfmod.putvar(ncfileout_cplbypass,['LATIXY'], LATIXY)
                
                # data writing
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
                ## in CPL_BYPASS forcing-data format, dimension in (gridcell, DTIME) 
                
                # for multiple zones', read-in data are stacked altogether along grid idx
                if iz==0:
                    zg0=0
                else:
                    zg0=zg0+len(zlines_elm[z-1])
                varvals = deepcopy(sdata_jointed[(zlines_elm[z]-1)+zg0,:])
                
                # varatts must in format: 'varname::att=att_val; varname::att=att_val; ...'
                varatts = varname+'::add_offset='+str(add_offset)+ \
                    ';'+varname+'::scale_factor='+str(scale_factor)
                error=nfmod.putvar(ncfileout_cplbypass, [varname], varvals, varatts=varatts)
                #error=nfmod.putvar(ncfileout_cplbypass, [varname], varvals)
                if error!=0: sys.exit('nfmod.putvar WRONG-'+varname+'-'+ncfileout_cplbypass)
            # end if NCOUT_CPLBYPASS
        #end of for loop of iz in cplbypass zones
    #
    else:
        print('TODO: CPL_BYPASS format other than "GSWP3" not yet supported')
        sys.exit(-1)
    
#end if nc_write or nc_create


