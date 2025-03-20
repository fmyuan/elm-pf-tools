#!/usr/bin/env python

import sys
from optparse import OptionParser
import numpy as np
from types import SimpleNamespace

# customized modules
from pytools.metdata_processing import elm_metdata_read
from pytools.metdata_processing import met_utils
from pytools.metdata_processing import elm_metdata_write as wrt



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
parser.add_option("--user_metvars", dest="user_metvars", default="", \
                  help="user-defined met file(s) var names by exact order of 'LONGXY,LATIXY,time,TBOT,PRECTmms,QBOT,FSDS,FLDS,PSRF,WIND' ")
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

if('site' in options.nc_write_mettype.lower() \
   or 'cplbypass_site' in options.nc_write_mettype.lower()): 
    vnames=['LONGXY','LATIXY','time', \
            'TBOT', 'PRECTmms', 'RH', 'FSDS', 'FLDS', 'PSRF', 'WIND']
else:
    vnames=['LONGXY','LATIXY','time', \
            'TBOT', 'PRECTmms', 'QBOT', 'FSDS', 'FLDS', 'PSRF', 'WIND']


#--------------------------------------------------------------------------------------
# read-in metdata from totally user-defined data

if (options.user_metfile!=''):

    metfile = options.user_metdir+'/'+options.user_metfile
    if (options.user_metvars ==''):
        user_metvars = vnames
    else:
        user_metvars = options.user_metvars.split(',')

    vardatas = {}
    if ('csv' in metfile):
        odata_header,odata = \
            elm_metdata_read.singleReadCsvfile(metfile)
        vardatas['time']=(odata['YEAR']-1.0)*365.0+(odata['DOY']-1.0) + \
                     (odata['HOUR']/24.0)  # days since 0001-01-01 00:00:00 and date already in no-leap
        vardatas['tunit'] = 'days since 0001-01-01 00:00'
        
        vardatas['prect_unit'] = 'mm/s'
        
        # may need to manually set lat/lon
        vardatas['LONGXY'] = [-70.899]
        vardatas['LATIXY'] = [42.757]
        #vardatas['LONGXY'] = [196.4215]
        #vardatas['LATIXY'] = [65.44037]
    elif ('nc' in metfile):
        odata_header,odata = \
            elm_metdata_read.singleReadNcfile(metfile, \
                                             uservars=user_metvars)
        vardatas['time']= odata['time'] # odata is in days since 1901-01-01 00:00:00
                        
        vardatas['tunit'] = 'days since 1901-01-01 00:00'
        vardatas['prect_unit'] = 'mm/s'

        vardatas['LONGXY']= odata['LONGXY']
        vardatas['LATIXY']= odata['LATIXY']

        for iv in vnames:
            if iv not in odata_header :
                sys.exit(iv + ' or its equivalent NOT in: '+metfile)
            if iv!='time': vardatas[iv] = odata[iv]
    
    #     
    # removal of leap-year doy366
    if ('DOY' in odata.keys()):
        lpyrindex=[i for i, x in enumerate(odata['DOY']) if x == 366.0 ]  
        # some datasets, eg. NCDC data,  is in leap-year calender. Here simply remove last day of the year
        vardatas['time']=np.delete(vardatas['time'],lpyrindex)
        for ivar in odata_header:
            if not ivar in ['time','tunit','prect_unit','LONGXY','LATIXY']:
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
            vardatas['RH'] = met_utils.convertHumidity(tk, pres_pa, q_kgkg=qbot)
    elif 'RH' in vardatas.keys() and 'QBOT' not in vardatas.keys():
        rh = np.squeeze(vardatas['RH'])
        if 'TBOT' in vardatas.keys():
            tk = np.squeeze(vardatas['TBOT'])
            if 'PSRF' in vardatas.keys():
                pres_pa = np.squeeze(vardatas['PSRF'])
            else:
                pres_pa = 101325.0
            vardatas['QBOT'] = met_utils.convertHumidity(tk, pres_pa, rh_100=rh)
    
    # if FLDS estimated from humidity and temperature
    if (options.estFLDS or 'FLDS' not in vardatas.keys()):
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
            
            vardatas['FLDS'] = \
                met_utils.calcFLDS(tk, pres_pa, q_kgkg=qbot, rh_100=rh)

#--------------------------------------------------------------------------------------
# save in ELM forcing data format
if (options.nc_create or options.nc_write):
    
    options_wrt = SimpleNamespace( \
            met_idir = options.user_metdir, \
            nc_create = options.nc_create, \
            nc_write = options.nc_write, \
            nc_write_mettype = options.nc_write_mettype )
    wrt.elm_metdata_write(options_wrt, vardatas)

#

