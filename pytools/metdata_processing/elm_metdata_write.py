#!/usr/bin/env python

import sys
from datetime import datetime
import glob
#from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset

# customized modules
from Modules_netcdf import putvar

def elm_metdata_write(options, metdata, time_dim=0):
    
    """
    options.*
            "met_idir, default="./", \
                E3SM met data directory for template files to look up
            "nc_create", default=False, \
                help = "output new nc file", action="store_true")
            "nc_write", default=False, \
                help = "output to a existed nc file", action="store_true")
            "nc_write_mettype", default="", \
                help = "output to nc files in defined format (default = '', i.e., as original)")

    metdata[*]
            LATIXY
            LONGXY
            YEAR
            DOY
            time
            FSDS
            FLDS
            PRECTmms
            PRSF
            RH or QBOT
            TBOT
            WIND
     
    time_dim=0
            time dimension indice in met data, by default 0, i.e. [time, ...]
            Note: if not, may need to swap for either standard ELM format or cpl_bypass format     
    
    """
   

    if('Site' in options.nc_write_mettype or 'cplbypass_Site' in options.nc_write_mettype): 
        vnames=['LONGXY','LATIXY','time', \
                'TBOT', 'PRECTmms', 'RH', 'FSDS', 'FLDS', 'PSRF', 'WIND']
    else:
        vnames=['LONGXY','LATIXY','time', \
                'TBOT', 'PRECTmms', 'QBOT', 'FSDS', 'FLDS', 'PSRF', 'WIND']
    
    
    #--------------------------------------------------------------------------------------
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
        ncfilein = ''; ncfilein_prv = ''
        ncfilein_cplbypass = ''; ncfilein_cplbypass_prv = ''
        metidir = options.met_idir
        if metidir=='': metidir='./'
        mdoy=[0,31,59,90,120,151,181,212,243,273,304,334,365]#monthly starting DOY
        
        for varname in vnames:
            if varname in ['LONGXY','LATIXY','time']: continue  #skip and always write if new output nc file
            if varname not in metdata.keys(): continue

            
            if 'GSWP3' in met_type:
                if (varname == 'FSDS'):
                    fdir = metidir+'/Solar3Hrly/'
                elif (varname == 'PRECTmms'):
                    fdir = metidir+'/Precip3Hrly/'
                else:
                    fdir = metidir+'/TPHWL3Hrly/'
                ncfilein = sorted(glob.glob("%s*.nc" % fdir))
                if len(ncfilein)>0: ncfilein = ncfilein[0]
                
                if '/cpl_bypass' not in metidir:                
                    fdirheader = metidir+'/cpl_bypass_template/GSWP3_daymet4_'+varname+'_'
                else:
                    fdirheader = metidir+'/GSWP3_daymet4_'+varname+'_'
                ncfilein_cplbypass = sorted(glob.glob("%s*.nc" % fdirheader))
                if len(ncfilein)<=0 and len(ncfilein_cplbypass)<=0:
                    sys.exit('there is NO file as -'+fdirheader)
                ncfilein_cplbypass = ncfilein_cplbypass[0]
                
            elif 'site' in met_type.lower():
                fdir = metidir+'/'
                # So 'metdir' must be full path, e.g. ../atm/datm7/CLM1PT_data/1x1pt_US-Brw
                ncfilein_cplbypass=fdir+'all_hourly.nc'
                
                ncfilein = sorted(glob.glob("%s/????-??.nc" % fdir))
                ncfilein = ncfilein[0]
            elif 'crujra' in met_type or 'CRUJRA' in met_type:
                if (varname == 'FSDS'):
                    fdir = metidir+'/Solar6Hrly/'
                elif (varname == 'PRECTmms'):
                    fdir = metidir+'/Precip6Hrly/'
                else:
                    fdir = metidir+'/TPHWL6Hrly/'
                ncfilein = sorted(glob.glob("%s*.nc" % fdir))
                if len(ncfilein)>0: ncfilein = ncfilein[0]

                if '/cpl_bypass' not in metidir:                
                    fdirheader = metidir+'/cpl_bypass_full/CRUJRAV2.3.c2023.0.5x0.5_'+varname+'_'
                else:
                    fdirheader = metidir+'/CRUJRAV2.3.c2023.0.5x0.5_'+varname+'_'
                ncfilein_cplbypass = sorted(glob.glob("%s*.nc" % fdirheader))
                if len(ncfilein)<=0 and len(ncfilein_cplbypass)<=0:
                    sys.exit('there is NO file as -'+fdirheader)
                ncfilein_cplbypass = ncfilein_cplbypass[0]

            elif 'ERA5' in met_type or 'era5' in met_type:
                if (varname == 'FSDS'):
                    fdir = metidir+'/SolarHrly/'
                elif (varname == 'PRECTmms'):
                    fdir = metidir+'/PrecipHrly/'
                else:
                    fdir = metidir+'/TPHWLHrly/'
                ncfilein = sorted(glob.glob("%s*.nc" % fdir))
                if len(ncfilein)>0: ncfilein = ncfilein[0]
                
                if '/cpl_bypass' not in metidir:                
                    fdirheader = metidir+'/cpl_bypass_full/Daymet_ERA5.1km_'+varname+'_'
                else:
                    fdirheader = metidir+'/Daymet_ERA5.1km_'+varname+'_'
                ncfilein_cplbypass = sorted(glob.glob("%s*.nc" % fdirheader))
                if len(ncfilein)<=0 and len(ncfilein_cplbypass)<=0:
                    sys.exit('there is NO file as -'+fdirheader)
                ncfilein_cplbypass = ncfilein_cplbypass[0]
            
            # new template nc file, and has a prv file 
            if options.nc_write and ncfilein_prv != ncfilein and \
                                    ncfilein_prv != '':
                options.nc_create = True
                options.nc_write = False
            if options.nc_write and ncfilein_cplbypass_prv != ncfilein_cplbypass and \
                                    ncfilein_cplbypass_prv != '':
                options.nc_create = True
                options.nc_write = False
                        
            # new met nc files to create or write
            
            #(1) in non-cplbypass format
            if (ncfilein!='' and len(ncfilein))>0:
                if 'site' in met_type.lower():
                    tname = 'time'
                    for iyr in range(int(np.min(metdata['YEAR'])), int(np.max(metdata['YEAR']))+1):
                        for imon in range(1,13):
                            ncfileout = str(iyr)+'-'+str(imon).zfill(2)+'.nc'
                            
                            tidx = np.argwhere((metdata['YEAR']==iyr) & 
                                               (metdata['DOY']>mdoy[imon-1]) & (metdata['DOY']<=mdoy[imon])) #DOY starting from 1
                            t_jointed=metdata['time'][tidx]
                            sdata = metdata[varname][tidx]
                            sdata = sdata[..., np.newaxis] # in 'Site' forcing-data format, dimension in (time, lat, lon), while 'sdata' is in (time, gridcell) 
                            
                            # create nc file
                            if options.nc_create:
                                with Dataset(ncfilein,'r') as src, Dataset(ncfileout, "w") as dst:
                                    
                                    #------- new nc dimensions
                                    for dname, dimension in src.dimensions.items():
                                        len_dimension = len(dimension)
                                        if dname == 'LATIXY': len_dimension = len(metdata['LATIXY'])
                                        if dname == 'LONGXY': len_dimension = len(metdata['LONGXY'])                                    
                                        if dname == tname: len_dimension = len(t_jointed)
                                        
                                        dst.createDimension(dname, len_dimension if not dimension.isunlimited() else None)
    
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
                                except Exception as e:
                                    print(e)
                                    tunit = ''
                                error=putvar(ncfileout, [tname], t, varatts=tname+'::units='+tunit)
                                error=putvar(ncfileout,['LONGXY'], metdata['LONGXY'])
                                error=putvar(ncfileout,['LATIXY'], metdata['LATIXY'])
                            #end of if create nc file
                            #
                            error=putvar(ncfileout, [varname], sdata)
                            if error!=0: sys.exit('nfmod.putvar WRONG-'+varname+'-'+ncfileout)
                        #end of month loop
                    #end of year loop
                    if (not NCOUT_CPLBYPASS and options.nc_create): # If not create cpl_bypass format file
                        options.nc_write = True
                        options.nc_create= False
            
            #(2) in cplbypass format (NOTE: can output both original and cplbypass)
            if (NCOUT_CPLBYPASS): print('cpl_bypass template file: '+ ncfilein_cplbypass)
            if (NCOUT_CPLBYPASS and ncfilein_cplbypass!=''):
                if 'site' in met_type.lower() or \
                   'GSWP3' in met_type or \
                   'crujra' in met_type.lower() or \
                   'ERA5' in met_type: 
                    if 'site' in met_type.lower():
                        ncfileout_cplbypass='all_hourly.nc'
                    elif 'gswp3' in met_type.lower() or \
                         'crujra' in met_type.lower() or \
                         'era5' in met_type.lower():
                        ncfileout_cplbypass=ncfilein_cplbypass.split('/')[-1]
                        if 'daymet' in met_type.lower():
                            ncfileout_cplbypass=ncfileout_cplbypass.replace('1901', '1980')
                    else:
                        print('currently only support 3 types of cpl_bypass: Site or GSWP3* or CRUJRA or ERA5*')
                        sys.exit(-1)
                        
                    tname = 'DTIME'
                    t_jointed=metdata['time']
                        
                    sdata = metdata[varname]
                    
                    if (options.nc_create):
                        with Dataset(ncfilein_cplbypass,'r') as src, Dataset(ncfileout_cplbypass, "w") as dst:
                            
                            #------- new nc dimensions
                            for dname, dimension in src.dimensions.items():
                                len_dimension = len(dimension)
                                if dname == 'n': len_dimension = len(metdata['LONGXY'])
                                if dname == tname: len_dimension = len(t_jointed)
                                
                                dst.createDimension(dname, len_dimension if not dimension.isunlimited() else None)
    
                            # variables
                            for vname in src.variables.keys():
                                vtype = src.variables[vname].datatype
                                vdim = src.variables[vname].dimensions
                                dst.createVariable(vname, vtype, vdim)
                                dst[vname].setncatts(src[vname].__dict__)
    
    
                        # ONLY create new nc ONCE, and save file I/O name
                        options.nc_create = False
                        options.nc_write = True
                        ncfilein_cplbypass_prv = ncfilein_cplbypass
                        #ncfileout_cplbypass_prv= ncfileout_cplbypass
                        
                        # time
                        try:
                            # tunit from template file, which need to be replaced by that from metdata
                            tunit = Dataset(ncfilein_cplbypass).variables[tname].getncattr('units')
                            t0=str(tunit.lower()).strip('days since')
                            if(t0.endswith(' 00:00:00')):
                                t0=t0
                            elif(t0.endswith(' 00:00')):
                                t0=t0+':00'
                            elif (t0.endswith(' 00')):
                                t0=t0+':00:00'
                            else:
                                t0=t0+' 00:00:00'
                            t0=datetime.strptime(t0,'%Y-%m-%d %X') # time must be in format '00:00:00'
                            
                            
                            if 'tunit' in metdata.keys():
                                if 'days since' in metdata['tunit']:
                                    #t_jointed alread in unit of days since ...
                                    tunit = metdata['tunit']
                                    t = t_jointed
                            
                            else:
                                #default, assuming days since year 0
                                yr0=np.floor(t_jointed[0]/365.0)+1
                                t=t_jointed - (yr0-1)*365.0
                                tunit = tunit.replace(str(t0.year).zfill(4)+'-', str(int(yr0)).zfill(4)+'-')
                                tunit = tunit.replace('-'+str(t0.month).zfill(2)+'-', '-01-')
                            
                            if(tunit.endswith(' 00') and not tunit.endswith(' 00:00:00')):
                                tunit=tunit+':00:00'
                        except Exception as e:
                            print(e)
                            tunit = ''
                        error=putvar(ncfileout_cplbypass, [tname], np.asarray(t), varatts=tname+'::units='+tunit)
                        # varatts must in format: 'varname::att=att_val; varname::att=att_val; ...'
                        if error!=0: sys.exit('nfmod.putvar WRONG')
                        
                        error=putvar(ncfileout_cplbypass,['LONGXY'], metdata['LONGXY'])
                        error=putvar(ncfileout_cplbypass,['LATIXY'], metdata['LATIXY'])
                        if 'site' in met_type.lower():
                            error=putvar(ncfileout_cplbypass,['start_year'], np.floor(t[0]/365.0))
                            error=putvar(ncfileout_cplbypass,['end_year'], np.floor(t[-1]/365.0))

                        else:
                            # except for 'site', other type of cpl_bypass requires zone_mapping.txt file
                            lon = metdata['LONGXY']
                            lat = metdata['LATIXY']
                            g_zno = 1
                            f = open('./zone_mappings.txt', 'w')
                            for ig in range(len(lon)):
                                f.write('%12.5f ' % lon[ig] )
                                f.write('%12.6f ' % lat[ig] )
                                f.write('%5d ' % g_zno )
                                f.write('%5d ' % (ig+1) )
                                f.write('\n')
                            f.close()


                                            
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
                    
                    # note if in CPL_BYPASS forcing-data format, dimension is in (gridcell, DTIME)
                    # but source met data may not.
                    varvals = sdata
                    if time_dim!=len(sdata.shape)-1:
                        varvals = np.swapaxes(sdata, time_dim, -1)
                    
                    # varatts must in format: 'varname::att=att_val; varname::att=att_val; ...'
                    varatts = varname+'::add_offset='+str(add_offset)+ \
                        ';'+varname+'::scale_factor='+str(scale_factor)
                    error=putvar(ncfileout_cplbypass, [varname], varvals, varatts=varatts)
                    #error=nfmod.putvar(ncfileout_cplbypass, [varname], varvals)
                    if error!=0: sys.exit('nfmod.putvar WRONG - '+varname+', - '+ncfileout_cplbypass)
                    
                
                else:
                    print('TODO: CPL_BYPASS format other than "Site/GSWP3/crujra/ERA5" not yet supported')
                    sys.exit(-1)
            #
            elif(NCOUT_CPLBYPASS):
                print('CPL_BYPASS format output required but cannot find a template netcdf file, such as: all_hourly.nc or GSWP3_TBOT_1901-2014_z14.nc')
                sys.exit(-1)
        # for varname in vnames
    
    print('DONE!')

#

# ---------------------------------------------------------------
"""
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

elm_metdata_write(options)

#--------------------------------------------------------------------------------------
"""

