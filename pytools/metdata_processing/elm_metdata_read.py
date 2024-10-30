#!/usr/bin/env python
import os, sys, math
from datetime import datetime, timedelta
#import calendar
import numpy as np
from netCDF4 import Dataset
from calendar import isleap
from numpy import intersect1d
import glob
import dask.dataframe as dd

from types import SimpleNamespace

try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False

####################################################################################################
#
# -------subset reading DAYMET data in nc format  ----------------
#

def subsetDaymetRead1NCfile(ncfile, lat_range=[], lon_range=[],SUBSETNC=False, ncopath=''):
#dimensions:
#    x = 7814 ;
#    y = 8075 ;
#    time = UNLIMITED ; // (365 currently)
#    nv = 2 ;
#variables:
#    float x(x) ;
#        x:units = "m" ;
#        x:long_name = "x coordinate of projection" ;
#        x:standard_name = "projection_x_coordinate" ;
#    float y(y) ;
#        y:units = "m" ;
#        y:long_name = "y coordinate of projection" ;
#        y:standard_name = "projection_y_coordinate" ;
#    float lat(y, x) ;
#        lat:units = "degrees_north" ;
#        lat:long_name = "latitude coordinate" ;
#        lat:standard_name = "latitude" ;
#    float lon(y, x) ;
#        lon:units = "degrees_east" ;
#        lon:long_name = "longitude coordinate" ;
#        lon:standard_name = "longitude" ;
#    float time(time) ;
#        time:long_name = "time" ;
#        time:calendar = "standard" ;
#        time:units = "days since 1980-01-01 00:00:00 UTC" ;
#        time:bounds = "time_bnds" ;
#    short yearday(time) ;
#        yearday:long_name = "yearday" ;
#    float time_bnds(time, nv) ;
#    short lambert_conformal_conic ;
#        lambert_conformal_conic:grid_mapping_name = "lambert_conformal_conic" ;
#        lambert_conformal_conic:longitude_of_central_meridian = -100. ;
#        lambert_conformal_conic:latitude_of_projection_origin = 42.5 ;
#        lambert_conformal_conic:false_easting = 0. ;
#        lambert_conformal_conic:false_northing = 0. ;
#        lambert_conformal_conic:standard_parallel = 25., 60. ;
#        lambert_conformal_conic:semi_major_axis = 6378137. ;
#        lambert_conformal_conic:inverse_flattening = 298.257223563 ;
#    float prcp(time, y, x) ;
#        prcp:_FillValue = -9999.f ;
#        prcp:long_name = "daily total precipitation" ;
#        prcp:units = "mm/day" ;
#        prcp:missing_value = -9999.f ;
#        prcp:coordinates = "lat lon" ;
#        prcp:grid_mapping = "lambert_conformal_conic" ;
#        prcp:cell_methods = "area: mean time: sum" ;
#
#// global attributes:
#        :start_year = 1980s ;
#        :source = "Daymet Software Version 3.0" ;
#        :Version_software = "Daymet Software Version 3.0" ;
#        :Version_data = "Daymet Data Version 3.0" ;
#        :Conventions = "CF-1.6" ;
#        :citation = "Please see http://daymet.ornl.gov/ for current Daymet data citation information" ;
#        :references = "Please see http://daymet.ornl.gov/ for current information on Daymet references" ;
#}
    
    daymet_vars = ['tmax', 'tmin', 'prcp', 'srad', 'vp','swe','dayl']
    #
    try:
        f0 = Dataset(ncfile,'r')
        print('\n read nc FILE: '+ncfile+' ------- ')
        
    except Exception as e:
        print(e)

    for ivar in f0.variables.keys():
        if ivar in daymet_vars:
            data_varname = ivar
            break

    geox = np.asarray(f0.variables['x'])
    geoy = np.asarray(f0.variables['y'])
    
    alllons = np.asarray(f0.variables['lon'])
    i=np.where(alllons<0.0)
    if(i[0].size>0): alllons[i]=alllons[i]+360.0
    alllats = np.asarray(f0.variables['lat'])
    
    if not SUBSETNC:
        alldata = np.asarray(f0.variables[data_varname])

    data = {}
    # timing
    data['time'] = np.asarray(f0.variables['time'])
    data['yearday'] = np.asarray(f0.variables['yearday'])
    t_indx = np.asarray(range(0,data['time'].size))

    #truncating spatially, if any
    lon_range = np.asarray(lon_range)
    lat_range = np.asarray(lat_range)
    if (len(lon_range)>0 and len(lat_range)<=0):
        lon_range[np.where(lon_range<0.0)]=lon_range[np.where(lon_range<0.0)]+360.0
        ij_indx = np.where( (alllons<=np.max(lon_range)) &
                            (alllons>=np.min(lon_range)) )
   
    elif (len(lat_range)>0 and len(lon_range)<=0):
        ij_indx = np.where( (alllats<=np.max(lat_range)) &
                            (alllats>=np.min(lat_range)) )

    elif (len(lon_range)>0 and len(lat_range)>0):
        lon_range[np.where(lon_range<0.0)]=lon_range[np.where(lon_range<0.0)]+360.0
        ij_indx = np.where( (alllons<=np.max(lon_range)) &
                            (alllons<=np.max(lon_range)) &
                            (alllats>=np.min(lat_range)) &
                            (alllats>=np.min(lat_range)) )
    else:
        ij_indx = []
   
    if len(ij_indx)>0: 
        #ij_indx is paired tuple for y/x, if not empty
        if len(ij_indx[0])>0:
            # but for x/y, they must be dimensionally continuous
            i_indx = np.asarray(range(np.min(ij_indx[1]), np.max(ij_indx[1])))
            j_indx = np.asarray(range(np.min(ij_indx[0]), np.max(ij_indx[0])))
        
        
            # subsetting NC file, if option on
            if SUBSETNC:
                x1=np.int(np.min(alllons[ij_indx]))
                x2=np.int(np.max(alllons[ij_indx]))
                y1=np.int(np.min(alllats[ij_indx]))
                y2=np.int(np.max(alllats[ij_indx]))
                
                ncfilename = ncfile.split('/')[-1]
                if ncfilename.startswith('subset'): # re-subsetting already subsetted data
                    oldhead = ncfilename.split('_')[0]+'_'+ncfilename.split('_')[1]
                    newhead = 'subsetX'+str(x1)+'-'+str(x2)+'_Y'+str(y1)+'-'+str(y2)
                    subncfile = ncfilename.replace(oldhead, newhead)
                else:
                    subncfile = 'subsetX'+str(x1)+'-'+str(x2)+'_Y'+str(y1)+'-'+str(y2)+'_' \
                        +ncfile.split('/')[-1]
                if (os.path.isfile(subncfile)):
                    print('Warning:  Overwriting existing subset ncfile '+subncfile)
                    os.system('rm -rf '+subncfile)
                else:
                    print('Warning:  Creating subset ncfile '+subncfile)
                    
                os.system(ncopath+'ncks -d x,'+str(np.min(i_indx))+','+str(np.max(i_indx))+ \
                                      ' -d y,'+str(np.min(j_indx))+','+str(np.max(j_indx))+ \
                        ' '+ncfile+' '+subncfile)

            
            else:
                data['geox'] = geox[i_indx]
                data['geoy'] = geoy[j_indx]
        
                data['lon']  = alllons[j_indx,][:,i_indx]
                data['lat']  = alllats[j_indx,][:,i_indx]

                data[data_varname]  = alldata[t_indx,][:,j_indx,][:,:,i_indx]


    elif not SUBSETNC:
        
        # all data
        data['geox'] = geox
        data['geoy'] = geoy
        
        data['lon']  = alllons
        data['lat']  = alllats
        data[data_varname]  = alldata
    
    return data
    
####################################################################################################
#
# -------Read multiple location DAYMET data, *.csv one by one ----------------
#

def singleDaymetReadCsvfile(filename):
    
    siteinfo = []

    #
    try:
        f0=open(filename,'r')
        for line in f0:
            if('Latitude' in line):
                linedata=line.rsplit(' ')
                indx=[i for i, s in enumerate(linedata) if 'Latitude' in s]
                lat = float(linedata[indx[0]+1])
                
            if('Longitude' in line):
                linedata=line.rsplit(' ')
                indx=[i for i, s in enumerate(linedata) if 'Longitude' in s]
                lon = float(linedata[indx[0]+1])

            if('Elevation' in line):
                linedata=line.rsplit(' ')
                indx=[i for i, s in enumerate(linedata) if 'Elevation' in s]
                elev = float(linedata[indx[0]+1])
                
            if('year,yday' in line):
                data_header=line.strip().rsplit(',')
                
                dataline=next(f0).strip()
                data_value=np.array([float(i) for i in dataline.rsplit(',')])
                
                dataline=next(f0).strip()
                while True:
                    data=np.array([float(i) for i in dataline.rsplit(',')])
                    temp=data_value
                    data_value=np.vstack((temp,data))
                    try: 
                        dataline=next(f0).strip()
                    except Exception as e:
                        print(e)
                        break
        
        #
        vardatas = {}
        for iv in range(len(data_header)):
            vardatas[data_header[iv]] = data_value[...,iv]
        del data_header, data_value     
        
    except Exception as e:
        print(e)

    
    #

    siteinfo.append(lon)
    siteinfo.append(lat)
    siteinfo.append(elev)
    
    return siteinfo, vardatas
    

####################################################################################################
#
# -------Read 1 site NCDC or GHCND met data, *.csv ----------------
#

def singleNCDCReadCsvfile(filename, ncdc_type, missing):
    

    # HEADER
    std_header = ("YEAR","MONTH","DOY","PRCP","SNWD","TAVG","TMAX","TMIN","PRES","RH") # standard data header
    if ('GHCND' in ncdc_type):
        site_header = ("")
        data_header = ("YEAR","MONTH","DOY","PRCP","SNDP","TEMP","MAX","MIN","","") # must be in same order as 'std_header'
    elif('NCDC' in ncdc_type):
        site_header = ("LATITUDE","LONGITUDE","ELEVATION")
        data_header = ("YEAR","MONTH","DOY","PRCP","SNWD","TAVG","TMAX","TMIN","","")
    
    #
    try:
        f0=open(filename,'r')

        indx_site =[]
        key_site=[]
        indx_data =[]
        key_data=[]

        siteinfo = {} # in case don't have site data

        for line in f0:
            linedata=line.rsplit(',')

            # HEADER line
            if any( ((s in line and len(s)>0) for s in site_header) or
                    ((s in line and len(s)>0) for s in data_header) ):
                
                for j in site_header:
                    indx=[i for i, s in enumerate(linedata) if j in s]
                    if (indx[0]>=0): 
                        indx_site.append(indx[0])
                        key_site.append(j)
                siteinfo = {}
                
                for j in data_header:
                    if (len(j)>0):
                        jidx = data_header.index(j)
                        indx=[i for i, s in enumerate(linedata) if j.strip() == s.strip()]
                        if (indx[0]>=0): 
                            indx_data.append(indx[0])
                            j_std = std_header[jidx] # the standard header corresponding to 'data_header'
                            key_data.append(j_std)
                data_value = {}
                
                continue
        
            # data line
            else:
                if(len(indx_site)>0 and len(siteinfo)<=0): # only read once
                    for i,indx in enumerate(indx_site):
                        key = key_site[i]
                        val = linedata[indx].strip()
                        siteinfo[key] = np.float(val)
                
                for i,indx in enumerate(indx_data):
                    key = key_data[i]
                    val = linedata[indx].strip()
                    if(val ==''):
                        val = np.float(missing)
                    else:
                        val = np.float(val)
                    if(val<=float(missing)): val=np.nan
                    if (not key in data_value.keys()):
                        data_value[key] = val
                    else:
                        data_value[key] = np.append(data_value[key],val)
                    
                
                

    except ValueError:
        print("Error in reading data")
    #

    # missing day(s) filling, if any
    yrmin=int(np.min(data_value['YEAR']))
    yrmax=int(np.max(data_value['YEAR']))
    tdata={}
    for ivar in data_value.keys(): tdata[ivar] = []
    for yr in range(yrmin,yrmax):
        days = 365
        if(isleap(yr)): days=366
        for doy in range(1,days+1):
            tdata['YEAR'] = np.append(tdata['YEAR'],yr)
            tdata['DOY'] = np.append(tdata['DOY'],doy)
            curdate=datetime(yr,1,1)+timedelta(doy-1)
            tdata['MONTH'] = np.append(tdata['MONTH'],curdate.date)

            yridx = np.argwhere(data_value['YEAR']==yr)
            dayidx = np.argwhere(data_value['DOY']==doy)
            idx=intersect1d(yridx, dayidx)
            for ivar in data_value.keys():
                if ivar not in ['YEAR','MONTH','DOY']:
                    if (len(idx)>0): 
                        if(len(idx)>1): # in one NCDC file, multiple stations may be available
                            ivar_val = np.nanmean(data_value[ivar][idx])
                        else:
                            ivar_val = data_value[ivar][idx[0]]
                        
                        if(ivar=='TAVG' and (not np.isfinite(ivar_val))): # daily average missing but max/min not
                            if( np.isfinite(any(data_value['TMAX'][idx])) and 
                                np.isfinite(any(data_value['TMIN'][idx])) ):
                                ivar_val = np.nanmean([data_value['TMAX'][idx],data_value['TMIN'][idx]])
                        
                        tdata[ivar] = np.append(tdata[ivar],ivar_val)
                    else:
                        tdata[ivar] = np.append(tdata[ivar],np.nan)
    
    # data source type dependent unit-relevant conversion
    if('metric' in ncdc_type and 'NCDC' in ncdc_type):
        if('SNWD' in tdata.keys()): tdata['SNWD'] = tdata['SNWD']*0.1 # from mm --> cm
        
    if('metric' not in ncdc_type and 'GHCND' in ncdc_type):
        if('TAVG' in tdata.keys()): tdata['TAVG'] = (tdata['TAVG']-32.0)*5.0/9.0 # from F --> degree C
        if('TMAX' in tdata.keys()): tdata['TMAX'] = (tdata['TMAX']-32.0)*5.0/9.0 # from F --> degree C
        if('TMIN' in tdata.keys()): tdata['TMIN'] = (tdata['TMIN']-32.0)*5.0/9.0 # from F --> degree C

        if('PRCP' in tdata.keys()): tdata['PRCP'] = tdata['PRCP']*2.54 # from inches per day --> mm/day
        if('SNWD' in tdata.keys()): tdata['SNWD'] = tdata['SNWD']*0.254 # from inches --> cm
    
    if('metric' in ncdc_type and 'GHCND' in ncdc_type):
        if('SNWD' in tdata.keys()): tdata['SNWD'] = tdata['SNWD']*0.1 # from mm --> cm
    

    return siteinfo, data_header, tdata
####################################################################################################
#
#------ ------ Read a met variables from CPL_BYPASS directory 
#
def clm_metdata_cplbypass_read(filedir,met_type, varnames, lons=[-999], lats=[-999]):
    #
    if (len(lons)!=len(lats)):
        print('paired "lons=,lats=" must be in same length', len(lons), len(lats) )
        sys.exit(-1)
    #
    if('GSWP3' in met_type or 'ESM' in met_type or 'CRUJRA' in met_type):
        # zone_mapping.txt
        f_zoning = filedir+'/zone_mappings.txt'
        all_lons=[]
        all_lats=[]
        all_zones=[]
        all_lines=[]
        with open(f_zoning) as f:
            dtxt=f.readlines()
            dtxt=filter(lambda x: x.strip(), dtxt)
            for d in dtxt:
                allgrds=np.array(d.split(),dtype=float)
                if(allgrds[0]<0.0): allgrds[0]=360.0+allgrds[0] # convert longitude in format of 0 - 360
                all_lons.append(allgrds[0])
                all_lats.append(allgrds[1])
                all_zones.append(int(allgrds[2]))
                all_lines.append(int(allgrds[3]))
                
        all_lons = np.asarray(all_lons)
        all_lats = np.asarray(all_lats)
        all_zones= np.asarray(all_zones)
        all_lines= np.asarray(all_lines)
        
        # look-up zone/line to extract
        zone = []; zline = {}
        i={}; j={}
        if(lons[0]<0):
            i[0] = np.asarray(range(len(all_lons)))
        elif(lons[0]==0):
            i[0] = [0]
        else:
            idx=np.where(lons<0.0)
            if len(idx[0])>0: lons[idx]=360.0+lons[idx] # convert longitude in format of 0 - 360
            for ix in range(len(lons)):
                lon_i=all_lons[np.argmin(abs(all_lons-lons[ix]))]
                i[ix]=np.where(all_lons == lon_i)[0]# in general, there are many numbers of 'lon' in 'all_lons' 
        
        if(lats[0]<0):
            j[0] = np.asarray(range(len(all_lats)))
        elif(lats[0]==0):
            j[0] = [0]
        else:
            for iy in range(len(lats)):
                lat_j=all_lats[np.argmin(abs(all_lats-lats[iy]))]
                j[iy]=np.where(all_lats == lat_j)[0]# in general, there are many numbers of 'lat' in 'all_lats' 
        
        if(len(i[0])>0 and len(j[0])>0): 
            for ixy in range(len(i)):
                if ixy==0:
                    ij = np.intersect1d(i[ixy], j[ixy])
                else:
                    ij = np.hstack((ij,np.intersect1d(i[ixy], j[ixy])))
        else:
            print('NOT found point: lon - '+str(lons)+ ' lat - '+str(lats))
            sys.exit()
        ij=np.unique(ij)
        zone = np.unique(np.array(all_zones)[ij])
        for iz in zone:
            idx = np.where(np.array(all_zones)[ij]==iz)[0]
            zline[iz] = np.unique(np.array(all_lines)[ij][idx])
    
    elif ('Site' in met_type):
        zone=[1]
        zline={}; zline[1]=1
    else:
        print('NOT supported CPLBYPASS format - ', met_type)
        sys.exit()
        

    #file
    met ={}
    vdims={}
        
    
    mdoy=[0,31,59,90,120,151,181,212,243,273,304,334] #monthly starting DOY
    for iv in range(len(varnames)):
        #if('GSWP3' in met_type or 'Site' in met_type or 'ESM' in met_type):
            v = varnames[iv]
            varlist=['FLDS','FSDS','PRECTmms','PSRF','QBOT','TBOT','WIND']
            if 'Site' in met_type:
                varlist=['FLDS','FSDS','PRECTmms','PSRF','RH','TBOT','WIND']
            if v not in varlist:
                print('NOT found variable:  '+ v)
                sys.exit()
                
            #
            for iz in zone:
                
                if ('GSWP3' in met_type):
                    if('v1' in met_type):
                        file=filedir+'/GSWP3_'+v+'_1901-2010_z'+str(int(iz)).zfill(2)+'.nc'
                    elif('daymet' in met_type):
                        file=filedir+'/GSWP3_daymet4_'+v+'_1980-2014_z'+str(int(iz)).zfill(2)+'.nc'
                    else:
                        file=filedir+'/GSWP3_'+v+'_1901-2014_z'+str(int(iz)).zfill(2)+'.nc'
                elif ('ESM' in met_type):
                    if('daymet' in met_type):
                        file=filedir+'/ESM_daymet4_'+v+'_1980-1980_z'+str(int(iz)).zfill(2)+'.nc'
                    else:
                        file=filedir+'/ESM_'+v+'_2021-2099_z'+str(int(iz)).zfill(2)+'.nc'
                
                elif('CRUJRA' in met_type):
                    file=filedir+'/CRUJRAV2.3.c2023.0.5x0.5_'+v+'_1901-2021_z'+str(int(iz)).zfill(2)+'.nc'

                elif('Site' in met_type):
                    file=filedir+'/all_hourly.nc'
                
                #
                fnc = Dataset(file,'r')
    
                if ('LATIXY' not in met.keys() and iv==0):
                    met['LATIXY']=np.asarray(fnc.variables['LATIXY'])[zline[iz]-1]
                elif(iv==0):
                    met['LATIXY']=np.hstack((met['LATIXY'],np.asarray(fnc.variables['LATIXY'])[zline[iz]-1]))
                if ('LONGXY' not in met.keys() and iv==0):
                    met['LONGXY']=np.asarray(fnc.variables['LONGXY'])[zline[iz]-1]
                elif(iv==0):
                    met['LONGXY']=np.hstack((met['LONGXY'],np.asarray(fnc.variables['LONGXY'])[zline[iz]-1]))
            
                if ('DTIME' not in met.keys()):
                    t=np.asarray(fnc.variables['DTIME'])
                    vdims['DTIME'] = fnc.variables['DTIME'].dimensions
                    if('units' in fnc.variables['DTIME'].ncattrs()):
                        tunit= fnc.variables['DTIME'].getncattr('units')
                        tunit=tunit.lower()
                        if ('days since' in tunit):
                            t0=str(tunit).strip('days since')
                            if (t0.endswith(' 00:00')): t0=t0+':00' # in case time written NOT exactly in 00:00:00
                            if (t0.endswith(' 00')): t0=t0+':00:00'
                            t0=datetime.strptime(t0,'%Y-%m-%d %X')
                            y0=t0.year
                            m0=t0.month
                            d0=t0.day-1.0
                            if ('daymet' in met_type):
                                days0 = (y0-1980)*365+mdoy[m0-1]+d0+t0.second/86400.0 # force days since 1901-01-01 00:00:00
                            else:
                                days0 = (y0-1901)*365+mdoy[m0-1]+d0+t0.second/86400.0 # force days since 1901-01-01 00:00:00
                            t = t + days0
                            met['DTIME']=t
                            met['tunit']='days since 1901-01-01 00:00:00'
                    
                    #correct 'DTIME' (when do this modify fnc to be 'r+')
                    if ('GSWP3' in met_type):
                        try:
                            dt=np.diff(met['DTIME'])
                            idx=np.where(dt!=0.125)
                            if len(idx)>1:
                                dt[idx]=0.125
                                dt=np.insert(dt, 0, met['DTIME'][0])
                                dt=np.add.accumulate(dt)
                                met['DTIME']=dt
                                daysnum = fnc.variables['DTIME']
                                daysnum[:] = dt
                        except Exception as e:
                            print(e)
    
                if(v in fnc.variables.keys()): 
                    d = fnc.variables[v][zline[iz]-1,:]
                    if (d.dtype==float or d.dtype==np.float32  or d.dtype==np.float64): 
                        d[np.where(d.mask)] = np.nan
                    else:
                        d[np.where(d.mask)] = -9999
                    d=np.asarray(d)
                    if(v not in met.keys()):
                        met[v] = d # 'scale_factor' and 'add_offset' have already used when reading nc data.
                    else:
                        met[v] = np.vstack((met[v],d))
                    if(v not in vdims.keys()):
                        vdims[v] = fnc.variables[v].dimensions
    
                # when done ncfile, close it.
                fnc.close()
                #
            # all vars done
        # end if of met-type
    # end of for loop of zone(s)
    #    
    return zone, zline, vdims,met

#################################################################################
#    
#------ ------ read a met variables from a standard CLM full met directory 
#
def clm_metdata_read(metdir,fileheader, met_type, met_domain, lon=[-999], lat=[-999], varnames=[]):
    #
    # datm domain file  
    all_lons=[]
    all_lats=[]
    UNSTRUCTURED = False
    try:
        f = Dataset(met_domain,'r')
    
        all_lons=np.asarray(f.variables['xc'])
        ni = np.where(all_lons<0.0)
        if(len(ni)>1): all_lons[ni] = all_lons[ni]+360.0 # convert longitude in format of 0 - 360

        all_lats=np.asarray(f.variables['yc'])
             
        all_lons = np.asarray(all_lons)
        all_lats = np.asarray(all_lats)
        if all_lats.shape[0]==1 and all_lats.shape[1]>1: UNSTRUCTURED=True
    except ValueError:
        print("Error in reading datm domain file:" + met_domain)
        
    i={}; j={}
    if(lon[0]==-999):
        i[0] = np.asarray(range(all_lons.size))
    else:
        idx=np.where(lon<0.0)
        if len(idx[0])>0: lon[idx]=360.0+lon[idx] # convert longitude in format of 0 - 360
        for ix in range(len(lon)):
            idx=np.unravel_index(np.argmin(abs(all_lons-lon[ix])), all_lons.shape)
            i[ix]=np.ravel_multi_index(np.where(all_lons == all_lons[idx]), all_lons.shape)
            # in general, there are many numbers of 'lon' in 'all_lons' 
        
    if(lat[0]==-999):
        j[0] = np.asarray(range(all_lats.size))
    else:
        for iy in range(len(lat)):
            jdx=np.unravel_index(np.argmin(abs(all_lats-lat[iy])), all_lats.shape)
            j[iy]=np.ravel_multi_index(np.where(all_lats == all_lats[jdx]), all_lats.shape)
            # in general, there are many numbers of 'lat' in 'all_lats' 
        
    if(len(i[0])>0 and len(j[0])>0): 
        for ixy in range(len(i)):
            if ixy==0:
                ij = np.intersect1d(i[ixy], j[ixy])
            else:
                ij = np.hstack((ij,np.intersect1d(i[ixy], j[ixy])))
    else:
        print('NOT found point: lon - '+str(lon)+ ' lat - '+str(lat))
        sys.exit()
    ij=np.unique(ij)
    ij=np.unravel_index(ij, all_lats.shape)
    
    if UNSTRUCTURED:
        iy = [0]
        ix = ij[1]
    else:
        iy=ij[0]
        if len(all_lats.shape)>1: 
            ix=ij[1]
        else:
            ix=[0]

    print('len ix,iy: ', len(ix), len(iy))

    #files data-reading
    met ={}
    vdims={}
    mdoy=[0,31,59,90,120,151,181,212,243,273,304,334]#monthly starting DOY
    if('GSWP3' in met_type or 'ESM' in met_type or \
       'Site' in met_type):
        
        if 'GSWP3' in met_type or 'ESM' in met_type:
            varlist=['FSDS','PRECTmms','FLDS','PSRF','QBOT','TBOT','WIND']
        if 'Site' in met_type:
            varlist=['FSDS','PRECTmms','FLDS','PSRF','RH','TBOT','WIND']

        if (len(varnames)<=0): varnames=varlist
        for v in varnames:
            if v not in varlist:
                print('NOT found variable:  '+ v)
                print('IN: '+ varlist)
                sys.exit()
                 
            
            # file directories
            if 'GSWP3' in met_type:
                if (v == 'FSDS'):
                    fdir = metdir+'/Solar3Hrly/'
                elif (v == 'PRECTmms'):
                    fdir = metdir+'/Precip3Hrly/'
                else:
                    fdir = metdir+'/TPHWL3Hrly/'
            elif 'ESM' in met_type:
                if (v == 'FSDS'):
                    fdir = metdir+'/Solar1Hrly/'
                elif (v == 'PRECTmms'):
                    fdir = metdir+'/Precip1Hrly/'
                else:
                    fdir = metdir+'/TPHWL1Hrly/'
            elif 'Site' in met_type:
                fdir = metdir+'/'
                # So 'metdir' must be full path, e.g. ../atm/datm7/CLM1PT_data/1x1pt_US-Brw

            #
            tvar = np.asarray([])
            vdata= np.asarray([])
            if (fileheader==''):
                dirfiles = sorted(os.listdir(fdir))
            else:
                fdirheader=fdir+fileheader.strip()
                dirfiles = sorted(glob.glob("%s*.*" % fdirheader))  # maybe file pattern in 'fileheader'
            if len(dirfiles)<=0:
                print('NO metdata file Found: ',v, ' IN ', fdir+fileheader)
                continue
            
            for dirfile in dirfiles:
                # in metdata directory
                if(os.path.isfile(dirfile)):
                    fnc = Dataset(dirfile,'r')

                    # only need to read once for all files
                    if ('LATIXY' not in met.keys()):
                        if 'LATIXY' in fnc.variables.keys():
                            met['LATIXY']=np.asarray(fnc.variables['LATIXY'])
                        elif 'lat' in fnc.variables.keys():
                            met['LATIXY']=np.asarray(fnc.variables['lat'])
                        if UNSTRUCTURED:
                            met['LATIXY'] = met['LATIXY'][ix]
                        else:
                            met['LATIXY'] = met['LATIXY'][iy,ix]

                    if ('LONGXY' not in met.keys()):
                        if 'LONGXY' in fnc.variables.keys():
                            met['LONGXY']=np.asarray(fnc.variables['LONGXY'])
                        elif 'lon' in fnc.variables.keys():
                            met['LONGXY']=np.asarray(fnc.variables['lon'])       
                        if UNSTRUCTURED:
                            met['LONGXY'] = met['LONGXY'][ix]
                        else:
                            met['LONGXY'] = met['LONGXY'][iy,ix]
                        
                    
                    # non-time variables from multiple files
                    if(v in fnc.variables.keys()): 
                        # for 'time', need to convert tunit AND must do for each var, due to from different files
                        if ('time' in fnc.variables.keys()):
                            t=np.asarray(fnc.variables['time'])
                            tdims = fnc.variables['time'].dimensions
                            if('units' in fnc.variables['time'].ncattrs()):
                                tunit= fnc.variables['time'].getncattr('units')
                                if ('days since' in tunit):
                                    t0=str(tunit).strip('days since')
                                    if(t0.endswith(' 00') and not t0.endswith(' 00:00:00')):
                                        t0=t0+':00:00'
                                    t0=datetime.strptime(t0,'%Y-%m-%d %X')
                                    #t0=t0-datetime.strptime('1901-01-01 00:00:00','%Y-%m-%d %X') 
                                    # the above is not right, because  'datetime' is a leap-year system
                                    y0=t0.year
                                    m0=t0.month
                                    d0=t0.day-1.0
                                    days0 = (y0-1901)*365+mdoy[m0-1]+d0+t0.second/86400.0
                                    t = t + days0

                                    tvar = np.hstack((tvar, t))
                            
                        # values
                        if(len(vdata)<=0):
                            if UNSTRUCTURED:
                                vdata=np.asarray(fnc.variables[v])[:,ix]
                            else:
                                vdata=np.asarray(fnc.variables[v])[:,iy,ix]
                        else:                                              
                            if UNSTRUCTURED:
                                d=np.asarray(fnc.variables[v])[:,ix]
                            else:
                                d=np.asarray(fnc.variables[v])[:,iy,ix]
                            vdata=np.vstack((vdata,d))
        
                        if(v not in vdims.keys()):
                            vdim = fnc.variables[v].dimensions

            # done with 1 variable
            # must sort both time and data due to 'dirfiles' are not time-sorted (and maybe varied for each variable)
            if(len(tvar)>0): 
                tt  = sorted(tvar)
                it = sorted(range(len(tvar)), key=lambda k: tvar[k])
                
                
                if ('time' not in met.keys()): # only need once, assuming all variables exactly same timing
                    met['time'] = tt
                    vdims['time'] = tdims
                    met['tunit']='days since 1901-01-01 00:00:00' # Force all time unit to be this one

                if(v not in met.keys()):
                    met[v]=vdata[it]
        
                if(v not in vdims.keys()):
                    vdims[v] = vdim
                
        # all vars done
    
    else:
        print('NOT supported Met_type: '+ met_type)
        sys.exit(-1)
        
    return ix,iy,vdims,met
    
####################################################################################################

#------ ------ Read a single ELM met NC file user-defined ---------------------------------

def singleReadNcfile(metdirfile, uservars=None, \
                    lons=[-999], lats=[-999]):
    #
    if (len(lons)!=len(lats)):
        print('paired "lons=,lats=" must be in same length', len(lons), len(lats) )
        sys.exit(-1)
    
    # standard ELM met vars
    elmvars=['LONGXY','LATIXY','time', \
            'TBOT', 'PRECTmms', 'QBOT', 'FSDS', 'FLDS', 'PSRF', 'WIND']
    
    #
    metdata = {}
    
    #file        
    fnc = Dataset(metdirfile,'r')

    # user-input met vars
    if uservars is None:
        # in this case, var names are exactly as 'elmvars'
        uservars = list(fnc.variables.keys())
        for v in elmvars:
            if v not in uservars:
                print('Error: '+v+' not in: '+metdirfile )
                exit(-1)
    else:
        print('user-defined var name: ', uservars)
        print('make sure they are in order of: ', elmvars)       
    
    # grid x/y variable name in nc file
    var_x = ['LONGXY','lon','longitude','x','geox']
    var_y = ['LATIXY','lat','latitude','y','geoy']
    gridxname = uservars[1]
    gridyname = uservars[0]
    if (gridyname in fnc.variables.keys()):        
        metdata['LATIXY']=np.asarray(fnc.variables[gridyname])
    else:
        print('Warning: grid lat/y varname is NOT one of ',var_y)
    if (gridxname in fnc.variables.keys()):        
        metdata['LONGXY']=np.asarray(fnc.variables[gridxname])
    else:
        print('Warning: grid lon/x varname in NOT one of ', var_x)
    
    # time variable
    vtime = uservars[2] #'time'
    mdoy=[0,31,59,90,120,151,181,212,243,273,304,334]#monthly starting DOY

    if (vtime in fnc.variables.keys()):
        t=np.asarray(fnc.variables[vtime])
        
        if('units' in fnc.variables[vtime].ncattrs()):
            tunit= fnc.variables[vtime].getncattr('units')
            tunit=tunit.lower()
            if ('days since' in tunit or 'hours since' in tunit):
                if 'days since' in tunit: t0=str(tunit).strip('days since')
                if 'hours since' in tunit: t0=str(tunit).strip('hours since')
                
                # t starting point
                if (t0.endswith(' 00:00')): 
                    t0=t0+':00' # in case time written NOT exactly in 00:00:00
                elif (t0.endswith(' 00')): 
                    t0=t0+':00:00'
                else:
                    t0=t0+' 00:00:00'
                t0=datetime.strptime(t0,'%Y-%m-%d %X')
                y0=t0.year
                m0=t0.month
                d0=t0.day-1.0
                days0 = (y0-1901)*365+mdoy[m0-1]+d0+t0.second/86400.0 # force days since 1901-01-01 00:00:00
                
                if 'days since' in tunit: t = t + days0
                if 'hours since' in tunit: t = t/24.0 + days0
                
                metdata['time']=t
                metdata['tunit']='days since 1901-01-01 00:00:00'

    # met variables
    for iv in range(len(uservars)):
        v = uservars[iv].strip()   # user-provided var name
        elmv = elmvars[iv] # standard var name
        if v in [vtime, gridxname, gridyname]: continue  # skip 
                
        #
        if(v in fnc.variables.keys()): 
            d = fnc.variables[v][...]
            if(v not in metdata.keys()):
                metdata[elmv] = d 
            else:
                metdata[elmv] = np.vstack((metdata[elmv],d))
    
    # when done ncfile, close it.
    fnc.close()
    #
    #  
    metdata_header = metdata.keys()
   
    return metdata_header, metdata
    
####################################################################################################
#
# -------Read 1 site met data, *.csv, TOTALLY user-defined ----------------
#
def singleReadCsvfile(filename, hrly=99, tc2tk=False, PRCP2s=False):
    # filename: csv data with one line header
    # hrly: optional, hours if aggregating data to more longer interval than orginal data 
    # tc2tk: optional, convert Tair from degreeC to degreeK
    # PRECP2s: optional, covert Preciptation unit to per second (from per interval or period) 
    
    # std elm var and its headers to be extracted from 'data_header', 
    # for ELM offline forcing 
    # 'ELM std name:csv column header;csv column header2'
    
    # Toolik '5_minute_data_v21.csv'
    #std_header = ('data:date', 'hour:hour', \
    #              'TBOT:air_temp_3m_5min;air_temp_5m_5min', \
    #              'FLDS:lw_dn_avg_5min', \
    #              'FSDS:sw_dn_avg_5min', \
    #              'RH:rh_3m_5min;rh_5m_5min', \
    #              'WIND:wind_sp_5m_5min')
    
    # Toolik '1_hour_data.csv'
    #std_header = ('date:date', 'hour:hour', \
    #              'TBOT:air_temp_1m;air_temp_3m;air_temp_5m', \
    #              'FLDS:lw_dn_avg', \
    #              'FSDS:sw_dn_avg', \
    #              'PRECTmms:rain;rain2', \
    #              'RH:rh_1m;rh_3m;rh_5m', \
    #              'WIND:wind_sp_1m;wind_sp_5m')
    
    # WERC-ATLAS data in SP,AK
    #std_header = ('date:date', 'hour:hour', \
    #'TBOTc:Air Temperature @ 1 m (C); Air Temperature @ 3 m (C)', \
    #'RH:Relative Humidity @ 1 m (%);Relative Humidity @ 3 m (%)', \
    #'WIND:Wind Speed (m/s)', \
    #'PRECTmms:Rainfall(mm)')

    # PIE met towers 'Met_Final_for ELM.csv'
    # 'ELM std name:csv column header;csv column header2'
    std_header = ('date:date', 'hour:hour', \
                  'TBOT:TBOT', \
                  'RH:RH', \
                  'WIND:Wind', \
                  'PSRF:PSRF', \
                  'FSDS:FSDS', \
                  'FLDS:', \
                  'PRECTmms:PRECTmms')

    
    #missing ='nan'
    missing =''
    #
    try:
        f0 = dd.read_csv(filename, parse_dates=['date'])
        #f0 = dd.read_csv(filename, dtype={'hour': 'float'}) # if hour in format hh:mm:ss
        
        data_name  = []
        data_value = {}

        cols_header = f0.columns.values
        for j in std_header:
            j_std = j.strip().split(':')[0]
            j_cols= j.strip().split(':')[1]
            j_cols= j_cols.strip().split(';')
            if (len(j_cols)>0 and j_cols[0]!=''):
                data_name.append(j_std)# the standard header corresponding to 'data_header'
                key_data =[]
                for cols in j_cols:
                    print(cols)
                    indx = np.argwhere(cols_header==cols)[0]
                    if len(indx)>0:key_data.append(cols)
        
                for key in key_data:
                    col_val = np.asarray(f0[key])
                    if (not j_std in data_value.keys()):
                        data_value[j_std] = col_val
                    else:
                        if (np.ndim(data_value[j_std])==1):
                            data_value[j_std] = np.append(data_value[j_std][:,None],col_val[:,None],axis=1) # mutiple data source for variable j_std
                        else:
                            data_value[j_std] = np.append(data_value[j_std],col_val[:,None],axis=1) # mutiple data source for variable j_std
                            
                #
                if np.ndim(data_value[j_std])>1:
                    data_value[j_std] = np.nanmean(data_value[j_std], axis=1)
                

    except ValueError:
        print("Error in reading data")
    #

    #
    del f0
    
    # datetime conversion
    datestr=np.datetime_as_string(data_value['date'], unit='D')
    y=[int(s.split('-')[0]) for s in datestr]
    data_value['YEAR']=np.asarray(y)
    m=[int(s.split('-')[1]) for s in datestr]
    data_value['MONTH']=np.asarray(m)
    d=[int(s.split('-')[2]) for s in datestr]
    data_value['DAY']=np.asarray(d)
    data_value['HOUR']=np.asarray([math.floor(t) for t in data_value['hour']])
    data_value['MM']=np.asarray([(t-math.floor(t))*60.0 for t in data_value['hour']])
    ix=np.argwhere(data_value['HOUR']==24)
    if (len(ix)>0): # data @24:00 is dated as day of @00:00, so have to do the following hack
        data_value['HOUR'][ix]=23
        data_value['MM'][ix]=59
    

    # 
    yrmin=int(np.min(data_value['YEAR']))
    yrmax=int(np.max(data_value['YEAR']))
    tdata={}
    for ivar in data_value.keys(): 
        if (ivar not in ['date', 'hour','MM']): tdata[ivar] = []
    if 'DOY' not in data_value.keys(): tdata['DOY']=[]
    
    # for hourly data aggregation, if needed
    if hrly==99: # no temporal data aggregation
        hrly=data_value['hour'][1]-data_value['hour'][0]
    # redo timing in sequence SO THAT any missing data will be caught OR integrating if needed
    for yr in range(yrmin,yrmax+1):
        days = 365
        if(isleap(yr)): days=366
        for doy in range(1,days+1):
            for hour in np.arange(0,24,hrly):
                tdata['YEAR'] = np.append(tdata['YEAR'],yr)
                tdata['DOY'] = np.append(tdata['DOY'],doy)
                curdate=datetime(yr,1,1)+timedelta(doy-1)
                tdata['MONTH'] = np.append(tdata['MONTH'],curdate.month)
                tdata['DAY'] = np.append(tdata['DAY'],curdate.day)
                tdata['HOUR'] = np.append(tdata['HOUR'],hour)
    
                yridx = np.argwhere(data_value['YEAR']==yr)
                dayidx = np.argwhere((data_value['MONTH']==curdate.month) & \
                                     (data_value['DAY']==curdate.day))
                tidx = np.argwhere((data_value['HOUR']>=hour) & \
                                   (data_value['HOUR']<(hour+hrly)))   # So 3-hrly is the starting time
                idx=intersect1d(yridx, dayidx)
                idx=intersect1d(idx, tidx)
                for ivar in tdata.keys():
                    if ivar not in ['YEAR','MONTH','DAY','DOY','HOUR']:
                        if (len(idx)>0): 
                            if(len(idx)>1): 
                                ivar_val = np.nanmean(data_value[ivar][idx])
                            else:
                                ivar_val = data_value[ivar][idx[0]]
                                                    
                            tdata[ivar] = np.append(tdata[ivar],ivar_val)
                        else:
                            tdata[ivar] = np.append(tdata[ivar],np.nan)
    
    #
    # pass to clm data format
    del data_value

    clm_header = tdata.keys()
    # some misc. data conversion
    if tc2tk:     
        tdata['TBOT'] = tdata['TBOT']+273.15
    if PRCP2s:
        tdata['PRECTmms'] = tdata['PRECTmms']/hrly
    
    return clm_header, tdata


#
# ------- CRU-JRA met forcing data convertion-extraction for CLM/ELM ------------------------------------
#

#
def clm_metdata_CRUJRA(crujra_dirfilehead, template_clm_metfile=''):
    # CRU-JRA data formats
    #  file naming (a full year, 6-hourly): 
        # crujra.v2.3.5d.$VAR.yyyy.365d.noc.nc
        #  where, $VAR - tmp, spfh, pre, pres, dswrf, dlwrf, ugrid, vgrd
    # CLM/ELM data format from 'template', but generally following CRU-NCEP format
    #    e.g.  under datm7/atm_forcing.datm7.CRUJAR.0.5d.v2.c2023/
    #          clmforc.CRUJRAV2.3.c2023.0.5x0.5.TPQWL.yyyy-mm.nc
    #          clmforc.CRUJRAV2.3.c2023.0.5x0.5.Solr.yyyy-mm.nc
    #          clmforc.CRUJRAV2.3.c2023.0.5x0.5.Prec.yyyy-mm.nc
    #   where,
    #     Solr  - from 'dswrf' (FSDS)
    #     Prec  - from 'pre' (PRECTmms) 
    #     TPQWL - from 'tmp' (TBOT), 'pres' (PSRF), 'spfh' (QBOT), 'ugrd/vgrd' (WIND), 'dlwrf' (FLDS)
    crujra_vars1 = ['time','lat','lon','dswrf']
    clm_vars1    = ['time','LATIXY','LONGXY','FSDS']

    crujra_vars2 = ['time','lat','lon','pre']
    clm_vars2    = ['time','LATIXY','LONGXY','PRECTmms']

    crujra_vars3 = ['time','lat','lon','tmp', 'pres', 'spfh', 'dlwrf', 'ugrd', 'vgrd']   # don't change the order here
    clm_vars3    = ['time','LATIXY','LONGXY','TBOT','PSRF', 'QBOT', 'FLDS', 'WIND']
    
    #checking where is the original data
    filehead = crujra_dirfilehead.split('/')[-1]
    crujradir  = crujra_dirfilehead.replace(filehead, '')
    dirfiles = os.listdir(crujradir)

    if template_clm_metfile=='':
        metfile_tpqwl0 = 'clmforc.CRUJRAV2.3.c2023.0.5x0.5.TPQWL.yyyy-mm.nc'
        metfile_solr0 = 'clmforc.CRUJRAV2.3.c2023.0.5x0.5.Solr.yyyy-mm.nc'
        metfile_prec0 = 'clmforc.CRUJRAV2.3.c2023.0.5x0.5.Prec.yyyy-mm.nc'
    else:
        metfile_new0 = template_clm_metfile
    
    for dirfile in dirfiles:
        # in metdata directory
        if(os.path.isfile(crujradir+'/'+dirfile)): 
            if(filehead in dirfile):
            
                print('\n file: '+dirfile)
    
                metfile = crujradir+'/'+dirfile
                metvar  = dirfile.split('.')[-5] # *.$VAR.yyyy.365d.noc.nc
                metyyyy = int(dirfile.split('.')[-4])
                
                fmetdata = Dataset(metfile,'r')
                metdata_t= np.asarray(fmetdata.variables['time']) # days since 1901-01-01 00:00:00
                if 'ugrd' in metfile:
                    metfile2 = metfile.replace('ugrd','vgrd')
                    fmetdata2 = Dataset(metfile2,'r')
                    
                # clm offline met filenames
                if template_clm_metfile =='':
                    if metvar in crujra_vars1: 
                        metfile_new0 = metfile_solr0
                        clm_vars    = clm_vars1
                        crujra_vars = crujra_vars1
                    if metvar in crujra_vars2: 
                        metfile_new0 = metfile_prec0
                        clm_vars    = clm_vars2
                        crujra_vars = crujra_vars2
                    if metvar in crujra_vars3: 
                        metfile_new0 = metfile_tpqwl0
                        clm_vars    = clm_vars3
                        crujra_vars = crujra_vars3
                    fmetdata_tmpl = Dataset(metfile_new0, 'r')
    
                # monthly files for clm offline met data
                mdays_range = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
                metdata_t_yr = metdata_t - (metyyyy - 1901)*365.0
                for im in range(1,13):
                    
                    v_ym = metfile_new0.split('.')[-3:]
                    metfile_new = metfile_new0.replace(v_ym[0]+'.'+v_ym[1]+'.'+v_ym[2], '')
                    metfile_new = metfile_new+v_ym[0]+'.' \
                        +str(metyyyy)+'-'+str(im).zfill(2)+'.nc'

                    #
                    print('dirfile: '+dirfile + '  =======>  '+metfile_new)

                    tindx = np.where((metdata_t_yr>=mdays_range[im-1]) & (metdata_t_yr<mdays_range[im]))

                    if not os.path.isfile(metfile_new):
                        dst = Dataset(metfile_new, mode='w',format=fmetdata_tmpl.file_format)
                        # create new nc file from template
                        # dimensions
                        #unlimited_dim=''
                        #unlimited_size=0
                        for dname, dimension in fmetdata_tmpl.dimensions.items():
                            if dname in fmetdata.dimensions:
                                if dname=='time':
                                    ldim = tindx[0].size
                                else:
                                    ldim = fmetdata.dimensions[dname].size
                                
                                dst.createDimension(dname, (ldim if not dimension.isunlimited() else None))
                                # unlimited dimension name to be included dst output
                                #if dimension.isunlimited():
                                #    unlimited_dim = dname
                                #    unlimited_size = ldim

                        #global attributes
                        if 'Solr' in metfile_new:
                            dst.case_title = "offline ELM format of CRUJRA V2.3 6-Hourly Atmospheric Forcing: Downward Shortwave Radiation"
                        elif 'Prec' in metfile_new:
                            dst.case_title = "offline ELM format of CRUJRA V2.3 6-Hourly Atmospheric Forcing: Total Precipitation"
                        elif 'TPQWL' in metfile_new:
                            dst.case_title = "offline ELM format of CRUJRA V2.3 6-Hourly Atmospheric Forcing: temperature, pressure, specific humidity, wind, downward Longwave Radiation"
                        dst.data_source = "Data is provided from the Japanese 55-year Reanalysis (JRA-55) project carried out by the Japan Meteorological Agency (JMA), and from CRU, UEA, Norwich UK"
                        dst.data_source_ftp = "ftp://ftp.ceda.ac.uk/badc/cru/data/cru_jra/cru_jra_2.3/"
                        dst.data_source_citation = ("University of East Anglia Climatic Research Unit; "+
                                             " Harris, I.C. (2022): CRU JRA v2.3: A forcings dataset of gridded land surface blend of Climatic Research Unit (CRU) and Japanese reanalysis (JRA) data;" +
                                             " Jan.1901 - Dec.2021.. NERC EDS Centre for Environmental Data Analysis,"+
                                             " date of citation. https://catalogue.ceda.ac.uk/uuid/38715b12b22043118a208acd61771917")
                        dst.references = "data source contact: i.harris@uea.ac.uk"
                        dst.history = '2023-04-21: CRU-JRA v2.3 data format covertion to offline ELM style'
                        dst.contact = 'F.-M. Yuan, CCSI/ESD-ORNL. yuanf@ornl.gov'


                        
                        # variables
                        for vname in fmetdata_tmpl.variables.keys():
                            
                            vdtype = fmetdata_tmpl.variables[vname].datatype
                            vdims  = np.asarray(fmetdata_tmpl.variables[vname].dimensions)
                            # new Filling value
                            try:
                                vFillvalue = fmetdata_tmpl.variables[vname]._FillValue
                            except:
                                vFillvalue = np.finfo(vdtype).max
                            
                            dst.createVariable(vname, vdtype, \
                                        dimensions=vdims, \
                                        zlib=True, fill_value=vFillvalue)
                            dst[vname].setncatts(fmetdata_tmpl[vname].__dict__)
                            if vname == 'time':
                                dst[vname].units = "days since "+str(metyyyy)+'-'+str(im).zfill(2)+"-01 00:00:00"

                            # data name/type in crujra
                            v_idx = clm_vars.index(vname)
                            if crujra_vars[v_idx] not in fmetdata.variables.keys(): continue
                            if vname!=crujra_vars[v_idx]: # checking variable name in crujra data
                                print(vname,' is in CRU-JRA data file: ', dirfile, 'as: ', crujra_vars[v_idx])
                            # missing/Filling in crujra dataset
                            try:
                                vFillMissing = fmetdata.variables[crujra_vars[v_idx]]._FillValue
                            except:
                                try:
                                    vFillMissing = fmetdata.variables[crujra_vars[v_idx]].missing_value
                                except:
                                    vFillMissing = np.finfo(fmetdata.variables[crujra_vars[v_idx]].datatype).max

                            # write data
                            if vname == 'time':
                                dst[vname][...] = metdata_t_yr[tindx] - mdays_range[im-1]
                            elif vname == 'LATIXY':
                                dst[vname][...] = np.transpose(np.tile(np.asarray(fmetdata[crujra_vars[v_idx]]),(fmetdata.dimensions['lon'].size,1)))
                            elif vname == 'LONGXY':
                                dst[vname][...] = np.tile(np.asarray(fmetdata[crujra_vars[v_idx]]),(fmetdata.dimensions['lat'].size,1))
                            elif 'time' in vdims and crujra_vars[v_idx] in fmetdata.variables.keys():
                                tmp_data = np.asarray(fmetdata[crujra_vars[v_idx]])[tindx,]
                                idx_err = np.where((tmp_data==vFillMissing) | (tmp_data>1.e10)) #due to data storage in nc, might has Fillvalue issue 
                                if vname=='WIND':
                                    tmp_data2= np.asarray(fmetdata2[crujra_vars[v_idx+1]])[tindx,]
                                    idx_err = np.where((tmp_data==vFillMissing) | (tmp_data>1.e10) \
                                                      | (tmp_data2==vFillMissing) | (tmp_data2>1.e10)) #due to data storage in nc, might has Fillvalue issue 
                                    tmp_data = tmp_data/2.0 + tmp_data2/2.0
                                elif vname=='FSDS':
                                    tmp_data = tmp_data/6.0/3600.0  # unit: J/m2/6hr --> w/m2
                                elif vname=='PRECTmms':
                                    tmp_data = tmp_data/6.0/3600.0  # unit: mm/6hr --> mm/s
                                if len(idx_err[0])>0: tmp_data[idx_err]=vFillvalue
                                dst[vname][...] = tmp_data
           
                        
                    else:
                        print('writing data into exiting file: ', metfile_new)
                        # TPQWL file actually combines 5 vars
                        dst = Dataset(metfile_new, mode='a')
                    
                        # write met data only
                        for vname in fmetdata_tmpl.variables.keys():
                            if vname not in ['time','LATIXY','LONGXY']:
                            
                                # data name/type in crujra
                                v_idx = clm_vars.index(vname)
                                if crujra_vars[v_idx] not in fmetdata.variables.keys(): continue
                                # missing/Filling in crujra dataset
                                try:
                                    vFillMissing = fmetdata.variables[crujra_vars[v_idx]]._FillValue
                                except:
                                    try:
                                        vFillMissing = fmetdata.variables[crujra_vars[v_idx]].missing_value
                                    except:
                                        vFillMissing = np.finfo(fmetdata.variables[crujra_vars[v_idx]].datatype).max

                                
                                vdtype = fmetdata_tmpl.variables[vname].datatype
                                vdims  = np.asarray(fmetdata_tmpl.variables[vname].dimensions)
                                # new missing/Filling
                                try:
                                    vFillvalue = fmetdata_tmpl.variables[vname]._FillValue
                                except:
                                    vFillvalue = np.finfo(vdtype).max
                                if 'time' in vdims and crujra_vars[v_idx] in fmetdata.variables.keys():
                                    if vname!=crujra_vars[v_idx]: # checking variable name in crujra data
                                        print(vname,' is in CRU-JRA data file: ', dirfile, 'as: ', crujra_vars[v_idx])
                                    tmp_data = np.asarray(fmetdata[crujra_vars[v_idx]])[tindx,]
                                    idx_err = np.where((tmp_data==vFillMissing) | (tmp_data>1.e10)) #due to data storage in nc, might has Fillvalue issue 
                                    if vname=='WIND':
                                        tmp_data2= np.asarray(fmetdata2[crujra_vars[v_idx+1]])[tindx,]
                                        idx_err = np.where((tmp_data==vFillMissing) | (tmp_data>1.e10) \
                                                          | (tmp_data2==vFillMissing) | (tmp_data2>1.e10)) #due to data storage in nc, might has Fillvalue issue 
                                        tmp_data = tmp_data/2.0 + tmp_data2/2.0
                                    elif vname=='FSDS':
                                        tmp_data = tmp_data/6.0/3600.0  # unit: J/m2/6hr --> w/m2
                                    elif vname=='PRECTmms':
                                        tmp_data = tmp_data/6.0/3600.0  # unit: mm/6hr --> mm/s
                                    if len(idx_err[0])>0: tmp_data[idx_err]=vFillvalue
                                    dst[vname][...] = tmp_data
                        #
                    #monthly ending
                    dst.close()
                # yearly ending
                fmetdata_tmpl.close()
                fmetdata.close()
            #file list checking (multiple years, with 1 file for each year)
        #file list ending
    #
    print('DONE!')
    #
#

def clm_metdata_Daymet_downscaled_read(daymetera5_dir, ts_hr=1, fileheader='', ptxyind=[], varnames=[]):
    #
    # Daymet_ERA5 or GSWP3 data formats
    #  file naming (a full year, hourly) for a 2degx2deg tile: 
        # year/tile_no/SolarHrly/clmforc.Daymet4.1km.Solr.yyyy-mm.nc
        # year/tile_no/PrecipHrly/clmforc.Daymet4.1km.Prec.yyyy-mm.nc
        # year/tile_no/TPHWLHrly/clmforc.Daymet4.1km.TPQWL.yyyy-mm.nc
        
        # (1) data in 2D, with x/y coordinates of daymet projection
        # (2) there is a file 'daymet_tiles.nc' for easy searching tile information,
        #        which includes: tile_no,tile_gindx,tile_xindx,tile_yindx,tile_geox,tile_geoy,tile_lat,tile_lon
    #

    #
    fileheader0 = fileheader  # save the input for multiple varnames

    #
    if len(ptxyind)>0:
        ix = ptxyind[:,0]
        iy = ptxyind[:,1]

    #files data-reading
    met ={}
    vdims={}
    mdoy=[0,31,59,90,120,151,181,212,243,273,304,334,365]#monthly starting DOY
    
    varlist=['FSDS','PRECTmms','FLDS','PSRF','QBOT','TBOT','WIND']

    ts_str=''
    if ts_hr!=1: ts_str=str(ts_hr)

    if (len(varnames)<1): varnames=varlist
    for v in varnames:
        if v not in varlist:
            print('NOT found variable:  '+ v)
            print('IN: '+ varlist)
            sys.exit()
                 
            
        # file directories
        if (v == 'FSDS'):
            fdir = daymetera5_dir+'/Solar'+ts_str+'Hrly/'
            if fileheader0=='': fileheader = 'clmforc.Daymet4.1km.Solr'
        elif (v == 'PRECTmms'):
            fdir = daymetera5_dir+'/Precip'+ts_str+'Hrly/'
            if fileheader0=='': fileheader = 'clmforc.Daymet4.1km.Prec'
        else:
            fdir = daymetera5_dir+'/TPHWL'+ts_str+'Hrly/'
            if fileheader0=='': fileheader = 'clmforc.Daymet4.1km.TPQWL'

        #
        tvar = np.asarray([])
        vdata= np.asarray([])
        if (fileheader==''):
            dirfiles = sorted(os.listdir(fdir))
        else:
            fdirheader=fdir+fileheader.strip()
            dirfiles = sorted(glob.glob("%s*.*" % fdirheader))  # maybe file pattern in 'fileheader'
        
        # find the mask
        if len(dirfiles)<1: return 
        
        fnc = Dataset(dirfiles[0],'r')
        vdata0 = fnc.variables[v][0,...]
        yindx, xindx = np.where((~vdata0.mask))
        # only need to read once for all files
        if ('LATIXY' not in met.keys()):
            if 'LATIXY' in fnc.variables.keys():
                met['LATIXY']=np.asarray(fnc.variables['LATIXY'])[yindx,xindx]
            elif 'lat' in fnc.variables.keys():
                met['LATIXY']=np.asarray(fnc.variables['lat'])[yindx,xindx]
        if ('LONGXY' not in met.keys()):
            if 'LONGXY' in fnc.variables.keys():
                met['LONGXY']=np.asarray(fnc.variables['LONGXY'])[yindx,xindx]
            elif 'lon' in fnc.variables.keys():
                met['LONGXY']=np.asarray(fnc.variables['lon'])[yindx,xindx]
        fnc.close()

        for dirfile in dirfiles:
            # in metdata directory
            if (not os.path.isfile(dirfile)): continue
            fnc = Dataset(dirfile,'r')
            # non-time variables from multiple files
            if(v in fnc.variables.keys()): 
                # for 'time', need to convert tunit AND must do for each var, due to from different files
                if ('time' in fnc.variables.keys()):
                    t=fnc.variables['time'][...]
                    if len(np.where(t.mask)[0])>0:
                        tindx = np.where(~t.mask)[0]
                        t=t[tindx]
                    else:
                        tindx = []
                    #
                    tdims = fnc.variables['time'].dimensions
                    if('units' in fnc.variables['time'].ncattrs()):
                        tunit= fnc.variables['time'].getncattr('units')
                        if ('days since' in tunit):
                            t0=str(tunit).strip('days since')
                            if(t0.endswith(' 00') and not t0.endswith(' 00:00:00')):
                                t0=t0+':00:00'
                            t0=datetime.strptime(t0,'%Y-%m-%d %X')
                            #t0=t0-datetime.strptime('1950-01-01 00:00:00','%Y-%m-%d %X') 
                            # the above is not right, because  'datetime' is a leap-year system
                            y0=t0.year
                            m0=t0.month
                            d0=t0.day-1.0
                            days0 = (y0-1950)*365+mdoy[m0-1]+d0+t0.second/86400.0
                            t = t + days0
                        #
                    #
                    tvar = np.hstack((tvar, t))
                            
                # values
                if 'time' in fnc.variables[v].dimensions:
                    d = np.asarray(fnc.variables[v])[:,yindx,xindx]
                    if len(tindx)>0: d = d[tindx,...]
                
                    if(len(vdata)<=0):
                        vdata=d
                    else:                                              
                        vdata=np.vstack((vdata,d))
                    
                    if(v not in vdims.keys()):
                        vdim = fnc.variables[v].dimensions
                #
            #
            fnc.close()
            
            # done with 1 dirfile

        # done with 1 variable for all dirfiles
        # must sort both time and data due to 'dirfiles' are not time-sorted (and maybe varied for each variable)
        if(len(tvar)>0): 
            tt  = sorted(tvar)
                
            it = sorted(range(len(tvar)), key=lambda k: tvar[k])
                                
            if ('time' not in met.keys()): # only need once, assuming all variables exactly same timing
                met['time'] = np.asarray(tvar[it])
                vdims['time'] = tdims
                print('checking timing length: ', len(met['time']),len(tt),v)
                met['tunit']='days since 1950-01-01 00:00:00' # Force all time unit to be this one
            elif (len(tvar) != len(met['time'])):
                met['time_'+v] = np.asarray(tvar[it])
                vdims['time_'+v] = tdims
                met['tunit_'+v]='days since 1950-01-01 00:00:00' # Force all time unit to be this one

            if(v not in met.keys()):
                met[v]=vdata[it]
        
            if(v not in vdims.keys()):
                vdims[v] = vdim
                
    # all vars done
                
    return vdims,met
    
####################################################################################################


####################################################################################################
#
#
# test modules
#
##################################################################################

#odata = \
#    subsetDaymetReadNCfile('daymet_v3_prcp_1980_na.nc', lon_range=[-170.0,-141.0], lat_range=[60.0, 90.0], SUBSETNC=True)


#site,odata_header,odata = \
#    singleDaymetReadCsvfile('daymet_kougarok-NGEE00.csv')

# read-in metdata from CPL_BYPASS_FULL
#cplbypass_dir='./cpl_bypass_full'
#cplbypass_fileheader=''
#cplbypass_mettype='GSWP3'
#lon = 195.25
#lat = 65.25  
#met_cplbypass = clm_metdata_cplbypass_read(cplbypass_dir,cplbypass_fileheader, cplbypass_mettype, 'TBOT', lon, lat)

#site,odata_header,odata = \
#    singleNCDCReadCsvfile('2016100_initproc.csv')


#odata_header,odata = \
#    singleReadCsvfile('5_minute_data_v21.csv')
    

#clm_metdata_CRUJRA('./rawdata/crujra.v2.3.5d')


"""
if HAS_MPI4PY:
    mycomm = MPI.COMM_WORLD
    myrank = mycomm.Get_rank()
    mysize = mycomm.Get_size()
else:
    mycomm = 0
    myrank = 0
    mysize = 1

print('Testing ', myrank)

# variable names
#vars=['QBOT','PSRF','FSDS','FLDS']
vars_name=['TBOT','WIND','PRECTmms']
vars_rank = []
for i in range(len(vars_name)):
    if myrank == i%mysize:
        vars_rank=np.append(vars_rank, int(i))


from pytools.metdata_processing.elm_metdata_write import elm_metdata_write


#--------------


#--------------
#tno=['11749','11750','11751','11929','11930','11931','11932','11933','11934','11935','11936','11937','11938']
tno=['14244']
for i in tno:
    for iv in range(len(vars_rank)):
        v = vars_name[int(vars_rank[iv])]
        print('TILE: ', i, int(vars_rank[iv]), v, myrank)
        if not os.path.isdir('./TILE'+i): continue
        metvd, metdata = clm_metdata_Daymet_downscaled_read('./TILE'+i+'/', ts_hr=1, fileheader='', ptxyind=[], varnames=[v])

        write_options = SimpleNamespace( \
            met_idir = './', \
            nc_create = True, \
            nc_write = False, \
            nc_write_mettype = 'cplbypass_ERA5_daymet' )

        elm_metdata_write(write_options, metdata) 

    if HAS_MPI4PY: mycomm.Barrier()
    if myrank==0:
        os.system('mv *.nc ./TILE'+i+'/cpl_bypass_full/')
        os.system('mv zone_mappings.txt ./TILE'+i+'/cpl_bypass_full/')

        print('Done Testing ./TILE',i)
"""
#--------------
"""
for iv in range(len(vars_rank)):
    v = vars_name[int(vars_rank[iv])]
    print('variable ', int(vars_rank[iv]), v, ' ON: ', myrank)
    ix, iy, metvd, metdata = clm_metdata_read('./', 'clmforc.GSWP3.c2011.0.5x0.5.', 'GSWP3', './domain.lnd.360x720_cruncep.c20190221_cavm1d.nc', lon=[-999], lat=[-999], varnames=[v])

    if len(metdata)<=0:
        print('data not found for variable ', v)
        continue
        
    write_options = SimpleNamespace( \
            met_idir = './', \
            nc_create = True, \
            nc_write = False, \
            nc_write_mettype = 'cplbypass_GSWP3' )

    elm_metdata_write(write_options, metdata)

if HAS_MPI4PY: mycomm.Barrier()
if myrank==0:
    os.system('mv GSWP3*.nc ./cpl_bypass_full/')
    os.system('mv zone_mappings.txt ./cpl_bypass_full/')

    print('Done Testing ./')

"""
