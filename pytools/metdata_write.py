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

parser = OptionParser()

# E3SM met data directory for template
parser.add_option("--e3sm_metdir", dest="met_idir", default="./", \
                  help="e3sminput met directory (default = ./, i.e., under current directory)")
# other sources of data
parser.add_option("--user_metdir", dest="user_metdir", default="./", \
                  help="user-defined met directory (default = ./, i.e., under current directory)")
parser.add_option("--user_metfile", dest="user_metfile", default="", \
                  help="user-defined met file(s) under user_metdir (default = '', i.e., None)")
parser.add_option("--estFLDS", dest="estFLDS", default=False, \
                  help = "Estimating FLDS", action="store_true")
parser.add_option("--nc_create", dest="nc_create", default=False, \
                  help = "output new nc file", action="store_true")
parser.add_option("--nc_write", dest="nc_write", default=False, \
                  help = "output to a existed nc file", action="store_true")
parser.add_option("--ncout_mettype", dest="nc_write_mettype", default="", \
                  help = "output to nc files in defined format (default = '', i.e., as original)")
#
(options, args) = parser.parse_args()

#--------------------------------------------------------------------------------------

if('Site' in options.nc_write_mettype or 'cplbypass_Site' in options.nc_write_mettype): 
    vnames=['TBOT', 'PRECTmms', 'RH', 'FSDS', 'FLDS', 'PSRF', 'WIND']
else:
    vnames=['TBOT', 'PRECTmms', 'QBOT', 'FSDS', 'FLDS', 'PSRF', 'WIND']


#--------------------------------------------------------------------------------------
# read-in metdata from totally user-defined data

if (options.user_metfile!=''):
    metfile = options.user_metdir+'/'+options.user_metfile
    odata_header,odata = \
        Modules_metdata.singleReadCsvfile(metfile)
    
    #LONGXY = [210.4033]
    #LATIXY = [68.6333]
    LONGXY = [196.4215]
    LATIXY = [65.44037]
    # need to recalculate time
    tvarname = 'time'    # variable name for time/timing
    
    vardatas = {}
    vardatas['time']=(odata['YEAR']-1.0)*365.0+(odata['DOY']-1.0) + \
                     (odata['HOUR']/24.0)  # days since 0001-01-01 00:00:00
    lpyrindex=[i for i, x in enumerate(odata['DOY']) if x == 366.0 ]  # NCDC data is in leap-year calender. Here simply remove last day of the year
    vardatas['time']=np.delete(vardatas['time'],lpyrindex)
    tunit = 'days since 0001-01-01 00:00'

    varsdims = {}
    for ivar in odata_header:
        varsdims[ivar]  = tvarname
        vardatas[ivar]=np.delete(odata[ivar],lpyrindex)
    del odata_header, odata

    # usually  either 'RH' or 'QBOT' is available, but if the other is required
    if 'QBOT' in vardatas.keys() and 'RH' not in vardatas.keys():
        qbot = np.squeeze(vardatas['QBOT'])
        if 'TBOT' in vardatas.keys():
            tk = np.squeeze(vardatas['TBOT'])
            if 'PSRF' in vardatas.keys():
                pres_pa = np.squeeze(vardatas['PSRF'])
            else:
                pres_pa = 101325.0
            vardatas['RH'] = Modules_metdata.convertHumidity(tk, pres_pa, q_kgkg=qbot)
    elif 'RH' in vardatas.keys() and 'QBOT' not in vardatas.keys():
        rh = np.squeeze(vardatas['RH'])
        if 'TBOT' in vardatas.keys():
            tk = np.squeeze(vardatas['TBOT'])
            if 'PSRF' in vardatas.keys():
                pres_pa = np.squeeze(vardatas['PSRF'])
            else:
                pres_pa = 101325.0
            vardatas['QBOT'] = Modules_metdata.convertHumidity(tk, pres_pa, rh_100=rh)
    # if FLDS estimated from humidity and temperature
    if (options.estFLDS):
            #Longwave radiation (calculated from air temperature, humidity)
            if 'TBOT' in vardatas.keys():
                tk = np.squeeze(vardatas['TBOT'])
            else:
                print('ERROR: for calculating FLDS, air temperature is required')
                sys.exit(-1)
            if 'PSRF' in vardatas.keys():
                pres_pa = np.squeeze(vardatas['PSRF'])
            else:
                pres_pa = 101325.0
            if 'QBOT' in vardatas.keys():
                qbot = np.squeeze(vardatas['QBOT'])
                rh = []
            elif 'RH' in vardatas.keys():
                rh = np.squeeze(vardatas['RH'])
                qbot = []
            else:
                print('ERROR: for calculating FLDS, either RH or QBOT is required')
                sys.exit(-1)
            
            vardatas['estFLDS'] = \
                Modules_metdata.calcFLDS(tk, pres_pa, q_kgkg=qbot, rh_100=rh)

#--------------------------------------------------------------------------------------
# save in ELM forcing data format
if (options.nc_create or options.nc_write):

    
    met_type = options.nc_write_mettype
    if 'cplbypass' in met_type:
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
    ncfilein = ''
    ncfilein_cplbypass = ''
    metdir = options.met_idir
    if metdir=='': metdir='./'
    mdoy=[0,31,59,90,120,151,181,212,243,273,304,334,365]#monthly starting DOY
    
    for varname in vnames:
        
        if 'GSWP3' in met_type:
            if (varname == 'FSDS'):
                fdir = metdir+'/Solar3Hrly/'
            elif (varname == 'PRECTmms'):
                fdir = metdir+'/Precip3Hrly/'
            else:
                fdir = metdir+'/TPHWL3Hrly/'
            ncfilein = sorted(glob.glob("%s*.nc" % fdir))
            ncfilein = ncfilein[0]
            
            fdirheader = metdir+'/GSWP3_'+varname+'_'
            ncfilein_cplbypass = sorted(glob.glob("%s*.nc" % fdirheader))
            ncfilein_cplbypass = ncfilein_cplbypass[0]
            
        elif 'site' in met_type.lower():
            fdir = metdir+'/'
            # So 'metdir' must be full path, e.g. ../atm/datm7/CLM1PT_data/1x1pt_US-Brw
            ncfilein_cplbypass=fdir+'all_hourly.nc'
            
            ncfilein = sorted(glob.glob("%s/????-??.nc" % fdir))
            ncfilein = ncfilein[0]
        
        # new met nc files to create or write
        
        #(1) in non-cplbypass format
        if (ncfilein!=''):
            if 'site' in met_type.lower():
                tname = 'time'
                for iyr in range(int(np.min(vardatas['YEAR'])), int(np.max(vardatas['YEAR']))+1):
                    for imon in range(1,13):
                        ncfileout = str(iyr)+'-'+str(imon).zfill(2)+'.nc'
                        
                        tidx = np.argwhere((vardatas['YEAR']==iyr) & 
                                           (vardatas['DOY']>mdoy[imon-1]) & (vardatas['DOY']<=mdoy[imon])) #DOY starting from 1
                        t_jointed=vardatas[tvarname][tidx]
                        sdata = vardatas[varname][tidx]
                        sdata = sdata[..., np.newaxis] # in 'Site' forcing-data format, dimension in (time, lat, lon), while 'sdata' is in (time, gridcell) 
                        
                        # create nc file
                        if options.nc_create:
                            nfmod.dupexpand(ncfilein, ncfileout, dim_name=tname, dim_len=len(t_jointed))

                            # time and locations
                            try:
                                tunit = Dataset(ncfilein).variables[tname].getncattr('units')
                                t0=str(tunit.lower()).strip('days since')
                                if (t0.endswith(' 00')):
                                    t0=t0+':00:00'
                                elif(t0.endswith(' 00:00')):
                                    t0=t0+':00'
                                t0=datetime.strptime(t0,'%Y-%m-%d %X') # time must be in format '00:00:00'
                                yr0=np.floor(t_jointed[0]/365.0)+1
                                t=t_jointed - (yr0-1)*365.0 - mdoy[imon-1]
                                tunit = tunit.replace(str(t0.year).zfill(4)+'-', str(int(iyr)).zfill(4)+'-')
                                tunit = tunit.replace('-'+str(t0.month).zfill(2)+'-', '-'+str(imon).zfill(2)+'-')
                                if(tunit.endswith(' 00') and not tunit.endswith(' 00:00:00')):
                                    tunit=tunit+':00:00'
                            except:
                                tunit = ''
                            error=nfmod.putvar(ncfileout, [tname], t, varatts=tname+'::units='+tunit)
                            error=nfmod.putvar(ncfileout,['LONGXY'], LONGXY)
                            error=nfmod.putvar(ncfileout,['LATIXY'], LATIXY)
                        #end of if create nc file
                        #
                        error=nfmod.putvar(ncfileout, [varname], sdata)
                        if error!=0: sys.exit('nfmod.putvar WRONG-'+varname+'-'+ncfileout)
                    #end of month loop
                #end of year loop
                if (not NCOUT_CPLBYPASS and options.nc_create): # If not write cpl_bypass format file
                    options.nc_write = True
                    options.nc_create= False
        
        #(2) in cplbypass format (NOTE: can output both original and cplbypass)
        if (NCOUT_CPLBYPASS): print('cpl_bypass template file: '+ ncfilein_cplbypass)
        if (NCOUT_CPLBYPASS and ncfilein_cplbypass!=''):
            if 'site' in met_type.lower()  or 'GSWP3' in met_type: 
                if 'site' in met_type.lower():
                    ncfileout_cplbypass='all_hourly.nc'
                elif 'GSWP3' in met_type:
                    ncfileout_cplbypass=ncfilein_cplbypass.split('/')[-1]
                else:
                    print('currently only support 2 types of cpl_bypass: Site or GSWP3*')
                    sys.exit(-1)
                    
                tname = 'DTIME'
                t_jointed=vardatas[tvarname]
                sdata = vardatas[varname]
                
                if (options.nc_create):
                    nfmod.dupexpand(ncfilein_cplbypass, ncfileout_cplbypass, dim_name=tname, dim_len=len(t_jointed))
                    # ONLY create new nc ONCE
                    options.nc_create = False
                    options.nc_write = True
                    
                    # time
                    try:
                        tunit = Dataset(ncfilein_cplbypass).variables[tname].getncattr('units')
                        t0=str(tunit.lower()).strip('days since')
                        if (t0.endswith(' 00')):
                            t0=t0+':00:00'
                        elif(t0.endswith(' 00:00')):
                            t0=t0+':00'
                        t0=datetime.strptime(t0,'%Y-%m-%d %X') # time must be in format '00:00:00'
                        yr0=np.floor(t_jointed[0]/365.0)+1
                        t=t_jointed - (yr0-1)*365.0
                        tunit = tunit.replace(str(t0.year).zfill(4)+'-', str(int(yr0)).zfill(4)+'-')
                        tunit = tunit.replace('-'+str(t0.month).zfill(2)+'-', '-01-')
                        if(tunit.endswith(' 00') and not tunit.endswith(' 00:00:00')):
                            tunit=tunit+':00:00'
                    except:
                        tunit = ''
                    error=nfmod.putvar(ncfileout_cplbypass, [tname], t, varatts=tname+'::units='+tunit)
                    # varatts must in format: 'varname::att=att_val; varname::att=att_val; ...'
                    if error!=0: sys.exit('nfmod.putvar WRONG')
                    
                    error=nfmod.putvar(ncfileout_cplbypass,['LONGXY'], LONGXY)
                    error=nfmod.putvar(ncfileout_cplbypass,['LATIXY'], LATIXY)
                    if 'site' in met_type.lower():
                        error=nfmod.putvar(ncfileout_cplbypass,['start_year'], np.floor(t[0]/365.0))
                        error=nfmod.putvar(ncfileout_cplbypass,['end_year'], np.floor(t[-1]/365.0))
                    
                
                #elif(options.nc_write):
                if (varname=='PRECTmms'): sdata=sdata/86400.0 # mm/day -> mm/s
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
                varvals = np.reshape(sdata, (1,-1)) # in CPL_BYPASS forcing-data format, dimension in (gridcell, DTIME) 
                
                # varatts must in format: 'varname::att=att_val; varname::att=att_val; ...'
                varatts = varname+'::add_offset='+str(add_offset)+ \
                    ';'+varname+'::scale_factor='+str(scale_factor)
                error=nfmod.putvar(ncfileout_cplbypass, [varname], varvals, varatts=varatts)
                #error=nfmod.putvar(ncfileout_cplbypass, [varname], varvals)
                if error!=0: sys.exit('nfmod.putvar WRONG-'+varname+'-'+ncfileout_cplbypass)
                
            
            else:
                print('TODO: CPL_BYPASS format other than "Site/GSWP3" not yet supported')
                sys.exit(-1)
        #
        elif(NCOUT_CPLBYPASS):
            print('CPL_BYPASS format output required but cannot find a template netcdf file, such as: all_hourly.nc or GSWP3_TBOT_1901-2014_z14.nc')
            sys.exit(-1)
    # for varname in vnames

print('DONE!')

#

