#!/usr/bin/env python
import os, sys, csv, time, math
from datetime import datetime, timedelta
import calendar
import numpy as np
from netCDF4 import Dataset
from pip._vendor.distlib.util import CSVReader
from calendar import isleap
from numpy import intersect1d
from matplotlib.pyplot import axis
import glob

#ncopath = '/usr/local/nco/bin/'
#####################################################################################################

# ---------------------------------------------------------------
# vapor pressure (pa) at saturated from air temperature (K)
def vpsat_pa(tk):
    a = (6.107799961, 4.436518521e-01, \
         1.428945805e-02,2.650648471e-04, \
         3.031240396e-06, 2.034080948e-08, \
         6.136820929e-11)

    b = (6.109177956, 5.034698970e-01, \
         1.886013408e-02, 4.176223716e-04, \
         5.824720280e-06, 4.838803174e-08, \
         1.838826904e-10)

    vpsat = np.zeros_like(tk)
    idx=np.where(tk>273.15)
    if (len(idx[0])>0):
        tc = tk[idx]-273.15
        vpsat[idx] = 100.*(a[0]+tc*(a[1]+tc*(a[2]+tc*(a[3]+tc*(a[4]+tc*(a[5]+tc*a[6]))))))
    
    idx=np.where(tk<=273.15)
    if (len(idx[0])>0):
        tc = tk[idx]-273.15
        vpsat[idx] = 100.*(b[0]+tc*(b[1]+tc*(b[2]+tc*(b[3]+tc*(b[4]+tc*(b[5]+tc*b[6]))))))

    return vpsat

# ---------------------------------------------------------------
# conversion btw specific (kg/kg) and relative humidity (percentage), 
# known air temperature (K) and pressure (pa)
def convertHumidity(tk, pres_pa, q_kgkg=[], rh_100=[]):
    
    vpsat = vpsat_pa(tk)
    d_pres = pres_pa - 0.378*vpsat
    qsat = 0.622*vpsat / d_pres
    if rh_100.size>0:
        q_kgkg = qsat*rh_100/100.0
        
        return q_kgkg
    
    elif q_kgkg.size>0:
        rh_100 = q_kgkg/qsat*100.0
        
        return rh_100
    else:
        print('ERROR: must provide either "q_kgkg" or "rh_100" ')

# ---------------------------------------------------------------
#Longwave radiation (calculated from air temperature, in K, specific humidity in kg/kg or RH in percentage)
def calcFLDS(tk, pres_pa, q_kgkg=[], rh_100=[]):
    
    CONST_STEBOL  = 5.67e-8      # Stefan-Boltzmann constant ~ W/m^2/K^4 
    
    if rh_100.size>0:
        q_kgkg = convertHumidity(tk, pres_pa, rh_100=rh_100)
    elif q_kgkg.size<=0:
        print('ERROR: must provide either "q_kgkg" or "rh_100" ')
    
    es =  pres_pa * q_kgkg /(0.622 + 0.378 * q_kgkg )
    ea = 0.70 + 5.95e-7 * es * np.exp(1500.0/tk)
    FLDS = ea * CONST_STEBOL * (tk**4)
    
    #
    return FLDS



#
# ------- GSWP3 met forcing data extraction for site or sites ------------------------------------
#
def clm_metdata_extraction(metdomainfile, metfiles, sites, ncopath=''):
    # ---- sites
    if (sites.len<2):
        print('sites must have paired location points: x/y or longitude/latidue')
        return
    else:
        sitex = sites[0] # x or longitudes
        sitey = sites[1] # y or latitudes
        if(sites.len>2): 
            sitez = sites[2]
        else:
            sitez=[]
    
    # ---- meteorological data domain
    domain = metdomainfile.rsplit('/')
    filelength=domain.len
    if(filelength>1):
        domain_dir=domainfile[0]
    else:
        domain_dir='./'   
    domain_file = domain[filelength-1]
    if(domain_file.endwith('.nc')):
       domain_file=domain[filelength-1].replace('.nc','') # removal of suffix of .nc
    
    #
    ncfile = metdomainfile
    try:
        f = Dataset(ncfile,'r')
        print('\n FILE: '+ncfile+' ------- ')
    except:
        print('\n Error in READING File: '+ncfile)
            
    allx=[]
    if('LONGXY' in f.variables.keys()):
        xkey = 'LONGXY'       
    elif('xc' in f.variables.keys()):
        xkey = 'xc'
    else:
        print('cannot find longitude coordinates: '+ncfile)
        exit()
    allx=np.asarray(f.variables[xkey])
    dimx=f.variables[xkey].dimensions
    if('lon' in dimx[0] or 'ni' in dimx[0]):
        allx = allx[:,0]                        # only need 1-D, although original data is in 2-D
    elif('lon' in dimx[1] or 'ni' in dimx[1]):
        allx = allx[0,:]
    # longitude in negative for 'W', but site in positive around, or vice versa    
    if(np.min(allx)<0.0):
        ni = np.where(sitex>180.0)
        if(len(ni)>1): sitex[ni] = sitex[ni]-360.0
    elif(np.max(allx)>180.0):
        ni = np.where(sitex<0.0)
        if(len(ni)>1): sitex[ni] = sitex[ni]+360.0 
        
    ally=[]
    if('LATIXY' in f.variables.keys()):
        ykey = 'LATIXY'       
    elif('yc' in f.variables.keys()):
        ykey = 'yc'
    else:
        print('cannot find latitude coordinates: '+ncfile)
        exit()
    ally=np.asarray(f.variables[ykey])
    dimy=f.variables[ykey].dimensions
    if('lat' in dimx[0] or 'nj' in dimy[0]):
        ally = ally[:,0]
    elif('lat' in dimx[1] or 'nj' in dimy[1]):
        ally = ally[0,:]
    
    
    ni = 0
    numxpts=0
    if(len(sitex)==1):
        numxpts = 1
        x = abs(allx-sitex);
        ni = np.argmin(x)    
    else:
        ni = np.where((allx>=sitex[0] and allx>=sitex[1]))
        numxpts = len(ni)
            
    nj = 0
    numypts=0
    if(len(sitey)==1):
        numypts = 1
        y = abs(ally-sitey);
        nj = np.argmin(y)    
    else:
        nj = np.where((ally>=sitey[0] and ally>=sitey[1]))
        numypts = len(nj)
    
    
    #---------------------------------------------------------------------------------------------------------
    print('Extracting domain data for: Site - ', sitex, sitey)
    print('Extracted grid: ni,nj - ', ni, nj, 'lon,lat -', allx[ni], ally[nj])

    domainfile_new = domain_dir+ \
             +domain_file+'_ptindxy_'+str(ni)+'_'+str(nj)+'.nc'
    
    
    if (os.path.isfile(domainfile_new)):
        print('Warning:  Removing existing domain file')
        os.system('rm -rf '+domainfile_new)
    
    os.system(ncopath+'ncks -d ni,'+str(ni)+','+str(ni+numxpts-1)+ \
                  ' -d nj,'+str(nj)+','+str(nj+numypts-1)+ \
              ' '+domainfile_orig+' '+domainfile_new)
    
    
    #--------------------------------------------------------------------------------------------------------
    #
    print('Extracting met-forcing data for: Site - ', site, sitex, sitey)
    
    pt_name = str(numxpts)+'x'+str(numypts)+'pt_'+site
    met_input_new = casedatadir+'/atm/datm7/CLM1PT_data/'+ pt_name
    os.system('mkdir -p ' + met_input_new)
    
    #checking where is the original data
    dirfiles = os.listdir(met_input)
    for dirfile in dirfiles:        
        filehead = metfilehead
        # in metdata directory
        if(os.path.isfile(met_input+'/'+dirfile)): 
            if(filehead in dirfile):
            
                print('\n file: '+dirfile)
    
                metfile_old = met_input+'/'+dirfile
                dirfile_new = dirfile.replace(filehead,pt_name)
                metfile_new = met_input_new+'/'+ dirfile_new
    
                #extracting data
                print('dirfile: '+dirfile + '  =======>  '+dirfile_new)
    
                os.system(ncopath+'ncks -a -O -d lon,'+str(ni)+','+str(ni+numxpts-1)+ \
                          ' -d lat,'+str(nj)+','+str(nj+numypts-1)+ \
                          ' '+metfile_old+' '+metfile_new)
        # in subdirectory of metdata directory
        elif(os.path.isdir(met_input+'/'+dirfile)):
            subfiles = os.listdir(met_input+'/'+dirfile)
            for subfile in subfiles:
    
                if(os.path.isfile(met_input+'/'+dirfile+'/'+subfile)):
                    if(filehead in subfile):
     
                        metfile_old = met_input+'/'+dirfile+'/'+subfile
                        
                        subfile_new = subfile.replace(filehead,pt_name)
                        metfile_new = met_input_new+'/'+subfile_new
     
                        #extracting data
                        print('Subfile: '+subfile + '  =======>  '+subfile_new)
                    
                        os.system(ncopath+'ncks -a -O -d lon,'+str(ni)+','+str(ni+numxpts-1)+ \
                                  ' -d lat,'+str(nj)+','+str(nj+numypts-1)+ \
                                  ' '+metfile_old+' '+metfile_new)
        
    print('DONE!')

####################################################################################################
#
# -------Read a subset of DAYMET data, in NC4 original format ----------------
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
        
    except:
        print('\n Error in READING File: '+ncfile)

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
                data_header=line.rsplit(',')
                
                dataline=next(f0).strip()
                data_value=np.array([float(i) for i in dataline.rsplit(',')])
                
                dataline=next(f0).strip()
                while True:
                    data=np.array([float(i) for i in dataline.rsplit(',')])
                    temp=data_value
                    data_value=np.vstack((temp,data))
                    try: 
                        dataline=next(f0).strip()
                    except:
                        break

    except ValueError:
        print("Error in reading data")

    
    #

    siteinfo.append(lon)
    siteinfo.append(lat)
    siteinfo.append(elev)
    
    return siteinfo, data_header, data_value
    

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

# ---------------------------------------------------------------
# read a met variables from CPL_BYPASS directory 
def clm_metdata_cplbypass_read(filedir,fileheader, met_type, lon, lat, vars):
    #
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
            all=np.array(d.split(),dtype=float)
            if(all[0]<0.0): all[0]=360.0+all[0] # convert longitude in format of 0 - 360
            all_lons.append(all[0])
            all_lats.append(all[1])
            all_zones.append(int(all[2]))
            all_lines.append(int(all[3]))
            
    all_lons = np.asarray(all_lons)
    all_lats = np.asarray(all_lats)
    all_zones= np.asarray(all_zones)
    all_lines= np.asarray(all_lines)
    
    
    if(lon==-999):
        i = 0
    else:
        if(lon<0.0): lon=360.0+lon # convert longitude in format of 0 - 360
        i=np.argmin(abs(all_lons-lon))
    lon_i=all_lons[i]
    i=np.where(all_lons == lon_i)[0]# in general, there are many numbers of 'lon' in 'all_lons' 

    if(lon==-999):
        j = 0
    else:
        j=np.argmin(abs(all_lats-lat))
    lat_j=all_lats[j]
    j=np.where(all_lats == lat_j)[0]# in general, there are many numbers of 'lat' in 'all_lats' 
        
    if(len(i)>0 and len(j)>0): 
        ij = np.intersect1d(i, j)
    else:
        print('NOT found point: lon - '+str(lon)+ ' lat - '+str(lat))
        sys.exit()
    zone = np.array(all_zones)[ij]
    line = np.array(all_lines)[ij]
    #because extraction going to do for zone/line, the grid is orderly cleared
    ix=np.asarray(range(len(ij))) 
    iy=np.asarray(range(len(ij))) 

    #file
    met ={}
    vdims={}
    if('GSWP3' in met_type):
        varlist=['FLDS','FSDS','PRECTmms','PSRF','QBOT','TBOT','WIND']
        for v in vars:
            if('v1' in met_type):
                file=filedir+'/GSWP3_'+v+'_1901-2010_z'+str(int(zone[0])).zfill(2)+'.nc'
            if('daymet' in met_type):
                file=filedir+'/GSWP3_Daymet3_'+v+'_1980-2010_z'+str(int(zone[0])).zfill(2)+'.nc'
            else:
                file=filedir+'/GSWP3_'+v+'_1901-2014_z'+str(int(zone[0])).zfill(2)+'.nc'

            if v not in varlist:
                print('NOT found variable:  '+ v)
                print('IN: '+ varlist)
                sys.exit()
                
            #
            fnc = Dataset(file,'r')

            if ('LATIXY' not in met.keys()):
                met['LATIXY']=np.asarray(fnc.variables['LATIXY'])[line[0]-1]
            if ('LONGXY' not in met.keys()):
                met['LONGXY']=np.asarray(fnc.variables['LONGXY'])[line[0]-1]
        
            if ('DTIME' not in met.keys()):
                met['DTIME']=np.asarray(fnc.variables['DTIME'])
                vdims['DTIME'] = fnc.variables['DTIME'].dimensions
                if('units' in fnc.variables['DTIME'].ncattrs()):
                    tunit= fnc.variables['DTIME'].getncattr('units')
                    if ('days since' in tunit):
                        met['tunit']='days'
                        t0=str(tunit).strip('days since')
                        t0=datetime.strptime(t0,'%Y-%m-%d %X')
                        met['t0_datetime']=t0
                
            #correct 'DTIME' (when do this modify fnc to be 'r+')
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
            except:
                print('NC file not writable')

            if(v in fnc.variables.keys()): 
                d=np.asarray(fnc.variables[v])[line[0]-1,:]
                met[v] = d # 'scale_factor' and 'add_offset' have already used when reading nc data.
                if(v not in vdims.keys()):
                    vdims[v] = fnc.variables[v].dimensions

            # when done ncfile, close it.
            fnc.close()
            
        # all vars done
        
    #because extraction going to do for zone/line, the grid is orderly cleared
    ix=np.asarray(range(len(ix))) 
    iy=np.asarray(range(len(iy))) 

    return ix,iy, vdims,met
    
# ---------------------------------------------------------------
# read a met variables from CLM full met directory 
def clm_metdata_read(metdir,fileheader, met_type, met_domain, lon, lat, vars):
    #
    # datm domain file  
    all_lons=[]
    all_lats=[]
    try:
        f = Dataset(met_domain,'r')
    
        all_lons=np.asarray(f.variables['xc'])
        ni = np.where(all_lons<0.0)
        if(len(ni)>1): all_lons[ni] = all_lons[ni]+360.0 # convert longitude in format of 0 - 360

        all_lats=np.asarray(f.variables['yc'])
             
        all_lons = np.asarray(all_lons)
        all_lats = np.asarray(all_lats)
    except ValueError:
        print("Error in reading datm domain file:" + met_domain)
        
    if(lon==-999):
        i = 0
    else:
        if(lon<0.0): lon=360.0+lon # convert longitude in format of 0 - 360
        i=np.argmin(abs(all_lons-lon))
    lon_i=all_lons[i]
    i=np.where(all_lons == lon_i)[0]# in general, there are many numbers of 'lon' in 'all_lons' 

    if(lon==-999):
        j = 0
    else:
        j=np.argmin(abs(all_lats-lat))
    lat_j=all_lats[j]
    j=np.where(all_lats == lat_j)[0]# in general, there are many numbers of 'lat' in 'all_lats' 
        
    if(len(i)>0 and len(j)>0): 
        ij,ix,iy = np.intersect1d(i, j, return_indices=True)
    else:
        print('NOT found point: lon - '+str(lon)+ ' lat - '+str(lat))
        sys.exit()

    #files data-reading
    met ={}
    vdims={}
    mdoy=[0,31,59,90,120,151,181,212,243,273,304,334]#monthly starting DOY
    if('GSWP3' in met_type or \
       'Site' in met_type):
        
        if 'GSWP3' in met_type:
            varlist=['FSDS','PRECTmms','FLDS','PSRF','QBOT','TBOT','WIND']
        if 'Site' in met_type:
            varlist=['FSDS','PRECTmms','FLDS','PSRF','RH','TBOT','WIND']

        if (len(vars)<=0): vars=varlist
        for v in vars:
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
            if 'Site' in met_type:
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
            for dirfile in dirfiles:
                # in metdata directory
                if(os.path.isfile(dirfile)):
                    fnc = Dataset(dirfile,'r')

                    # only need to read once for all files
                    if ('LATIXY' not in met.keys()):
                        met['LATIXY']=np.asarray(fnc.variables['LATIXY'])
                    if ('LONGXY' not in met.keys()):
                        met['LONGXY']=np.asarray(fnc.variables['LONGXY'])
        
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
                            vdata=np.asarray(fnc.variables[v])[:,ix,iy]
                        else:                                              
                            d=np.asarray(fnc.variables[v])[:,ix,iy]
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
        
    return ix,iy, vdims,met
    
    
##################################################################################
# test modules

#clm_metdata_extraction()

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
#met_cplbypass = clm_metdata_cplbypass_read(cplbypass_dir,cplbypass_fileheader, cplbypass_mettype, lon, lat)

#site,odata_header,odata = \
#    singleNCDCReadCsvfile('2016100_initproc.csv')


    
    
    



