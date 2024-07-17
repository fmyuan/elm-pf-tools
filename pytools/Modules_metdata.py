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
#import dask.dataframe as dd
from src.pytools.commons_utils import Modules_netcdf as ncmod
#from src.pytools.metdata_checking import tunit

#ncopath = '/usr/local/nco/bin/'
#####################################################################################################

# ---------------------------------------------------------------
# vapor pressure (pa) at saturated from air temperature (K)
def vpsat_pa(tk, freezing_eff=True):
    a = (6.107799961, 4.436518521e-01, \
         1.428945805e-02,2.650648471e-04, \
         3.031240396e-06, 2.034080948e-08, \
         6.136820929e-11)

    b = (6.109177956, 5.034698970e-01, \
         1.886013408e-02, 4.176223716e-04, \
         5.824720280e-06, 4.838803174e-08, \
         1.838826904e-10)

    vpsat = np.zeros_like(tk)
    
    tc = tk-273.15
    vpsat = 100.*(a[0]+tc*(a[1]+tc*(a[2]+tc*(a[3]+tc*(a[4]+tc*(a[5]+tc*a[6]))))))
    
    idx=np.where(tk<=273.15)
    if (len(idx[0])>0 and freezing_eff):
        tc = tk[idx]-273.15
        vpsat[idx] = 100.*(b[0]+tc*(b[1]+tc*(b[2]+tc*(b[3]+tc*(b[4]+tc*(b[5]+tc*b[6]))))))

    return vpsat

# ---------------------------------------------------------------
# conversion btw specific (kg/kg) and relative humidity (percentage), 
# known air temperature (K) and pressure (pa)
def convertHumidity(tk, pres_pa, q_kgkg=[], rh_100=[], vpsat_frz=True):
    
    vpsat = vpsat_pa(tk, vpsat_frz)
    d_pres = pres_pa - 0.378*vpsat
    qsat = 0.622*vpsat / d_pres
    if len(rh_100)>0:
        q_kgkg = qsat*rh_100/100.0
        idx=np.where(rh_100>100.0)
        if (len(idx[0])>0): q_kgkg[idx]=qsat[idx]
        
        return q_kgkg
    
    elif len(q_kgkg)>0:
        rh_100 = q_kgkg/qsat*100.0
        idx=np.where(rh_100>100.0)
        if (len(idx[0])>0): rh_100[idx]=100.0
        
        return rh_100
    else:
        print('ERROR: must provide either "q_kgkg" or "rh_100" ')

# ---------------------------------------------------------------
#Longwave radiation (calculated from air temperature, in K, specific humidity in kg/kg or RH in percentage)
def calcFLDS(tk, pres_pa, q_kgkg=[], rh_100=[]):
    
    CONST_STEBOL  = 5.67e-8      # Stefan-Boltzmann constant ~ W/m^2/K^4 
    
    if len(rh_100)>0:
        q_kgkg = convertHumidity(tk, pres_pa, rh_100=rh_100)
    elif len(q_kgkg)<=0:
        print('ERROR: must provide either "q_kgkg" or "rh_100" ')
    
    es =  pres_pa * q_kgkg /(0.622 + 0.378 * q_kgkg )
    ea = 0.70 + 5.95e-7 * es * np.exp(1500.0/tk)
    FLDS = ea * CONST_STEBOL * (tk**4)
    
    #
    return FLDS

#
#------- ------- extract a set of met variables from CPL_BYPASS directory 
#
def clm_metdata_cplbypass_extraction(filedir,met_type, lon, lat, ncopath='', z=0, l=0):
    #
    if ('Site' in met_type):
        zone=1
        line=1
        ni = 0

    else:
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
                if (not ('domain' in met_type or 'surfdata' in met_type)):
                    all_zones.append(int(all[2]))
                    all_lines.append(int(all[3]))
                else:
                    all_zones.append(1)
                    all_lines.append(1)                    
        f.close()
        all_lons = np.asarray(all_lons)
        all_lats = np.asarray(all_lats)
        all_zones= np.asarray(all_zones)
        all_lines= np.asarray(all_lines)
        
        if(lon<0.0): lon=360.0+lon # convert longitude in format of 0 - 360
        
        if len(all_lons)>1:
            dist2 = (all_lats-lat)*(all_lats-lat) + \
                    (all_lons-lon)*(all_lons-lon)
            ni=np.argmin(dist2)
            print('Nearest grid: ', ni, 
                  'dist:', math.sqrt(np.min(dist2)),
                  'xdist:', all_lons[ni]-lon, 
                  'xdist max:', max(np.fabs(np.diff(np.sort(all_lons)))),
                  'ydist:', all_lats[ni]-lat,
                  'ydist max:', max(np.fabs(np.diff(np.sort(all_lats)))))
            
            #out of bounds
            d=all_lons[ni]-lon
            dmax=np.fabs(np.diff(np.sort(all_lons)))
            if lon>np.max(all_lons) or lon<np.min(all_lons):
                if math.fabs(d)>max(dmax): return
            d=all_lats[ni]-lat
            dmax=np.fabs(np.diff(np.sort(all_lats)))
            if lat>np.max(all_lats) or lat<np.min(all_lats):
                if math.fabs(d)>max(dmax): return

            if (not ('domain' in met_type or 'surfdata' in met_type)):
                zone = np.array(all_zones)[ni]
                line = np.array(all_lines)[ni]
                print(' In Zone: ',zone, ' @line: ',line)
            
        else:
            # only 1 data
            ni = 0
            zone=np.array(all_zones)[ni]
            line=np.array(all_lines)[ni]
    ## 
    
    #files
    filedir_new = os.path.abspath('./subset')
    if (os.path.isdir(filedir_new)):
        os.system('rm -rf '+filedir_new)
    os.makedirs(filedir_new)
    
    # new zone_mappings
    zone_new = all_zones[ni]
    line_new = 1
    # using input zone/line
    if z>0: zone_new=z
    if l>0: line_new=l

    f_zoning = filedir_new+'/zone_mappings.txt'
    f = open(f_zoning, 'w')
    f.write("%12.5f " % all_lons[ni])
    f.write("%12.6f " % all_lats[ni])
    f.write("%5i " % zone_new)
    f.write("%5i\n" % line_new)
    #f.write("%5i\n" % line)
    f.close()
    
    # domain.nc/surfdata.nc/surfdata.nc, if for GSWP3_daymet4 or only for surfdata/domain nc(s)
    if (('GSWP3' in met_type and 'daymet4' in met_type) or 'domain' in met_type or 'surfdata' in met_type):
        os.system(ncopath+'ncks --no_abc -O -d ni,'+str(ni)+','+str(ni)+ \
                        ' '+filedir+'/domain.nc '+filedir_new+'/domain.nc')
        os.system(ncopath+'ncks --no_abc -O -d gridcell,'+str(ni)+','+str(ni)+ \
                        ' '+filedir+'/surfdata.nc '+filedir_new+'/surfdata.nc')
        os.system(ncopath+'ncks --no_abc -O -d gridcell,'+str(ni)+','+str(ni)+ \
                        ' '+filedir+'/surfdata.pftdyn.nc '+filedir_new+'/surfdata.pftdyn.nc')

      
    if('GSWP3' in met_type or 'Site' in met_type or 'CRUJRA' in met_type):
        varlist=['FLDS','FSDS','PRECTmms','PSRF','QBOT','TBOT','WIND']
        if 'Site' in met_type:
            varlist=['FLDS','FSDS','PRECTmms','PSRF','RH','TBOT','WIND']
        
        for v in varlist:
            if ('GSWP3' in met_type):
                if('v1' in met_type):
                    file='./GSWP3_'+v+'_1901-2010_z'+str(int(zone)).zfill(2)+'.nc'
                    file_new='./GSWP3_'+v+'_1901-2010_z'+str(int(zone_new)).zfill(2)+'.nc'
                elif('daymet' in met_type):
                    file='./GSWP3_daymet4_'+v+'_1980-2014_z'+str(int(zone)).zfill(2)+'.nc'
                    file_new='./GSWP3_daymet4_'+v+'_1980-2014_z'+str(int(zone_new)).zfill(2)+'.nc'
                else:
                    file='./GSWP3_'+v+'_1901-2014_z'+str(int(zone)).zfill(2)+'.nc'
                    file_new='./GSWP3_'+v+'_1901-2014_z'+str(int(zone_new)).zfill(2)+'.nc'
            
            elif('CRUJRA' in met_type):
                if('V2.3' in met_type):
                    file=met_type.strip()+'_'+v+'_1901-2021_z'+str(int(zone)).zfill(2)+'.nc'
                    file_new=met_type.strip()+'_'+v+'_1901-2021_z'+str(int(zone_new)).zfill(2)+'.nc'
                elif('V2.4' in met_type):
                    file=met_type.strip()+'_'+v+'_1901-2022_z'+str(int(zone)).zfill(2)+'.nc'
                    file_new=met_type.strip()+'_'+v+'_1901-2022_z'+str(int(zone_new)).zfill(2)+'.nc'
            
            elif('Site' in met_type):
                file='./all_hourly.nc'
                file_new='./all_hourly.nc'
            
            file_new = filedir_new+'/'+file_new
            file=filedir+'/'+file
            #
            #
            #extracting data
            print('extracting file: '+file + '  =======>  '+file_new)
            
            if ('Site' in met_type):
                os.system(ncopath+'ncks --no_abc -O -d n,'+str(ni)+','+str(ni)+ \
                                  ' '+file+' '+file_new)
            else:
                os.system(ncopath+'ncks --no_abc -O -d n,'+str(line-1)+','+str(line-1)+ \
                                  ' '+file+' '+file_new)
        # all vars done
        print('DONE!')
    
     
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
        domain_dir=metdomainfile[0]
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
    except Exception as e:
        print(e)
            
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
                    except Exception as e:
                        print(e)
                        break

    except Exception as e:
        print(e)

    
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

#------ ------ Read a met variables from CPL_BYPASS directory 
def clm_metdata_cplbypass_read(filedir,met_type, vars, lons=[-999], lats=[-999]):
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
    
    # the following is for testing
    #met['DTIME']= []
    #met['tunit']='days since 1901-01-01 00:00:00'
    #met['LATIXY']=[]
    #met['LONGXY']=[]
    #for iv in range(len(vars)):
    #    v = vars[iv]
    #    met[v] = []
    #    for iz in zone:
    #        if iv==0: 
    #            idx=np.where((all_zones==iz) & (np.isin(all_lines,zline[iz])))
    #            met['LONGXY']=np.hstack((met['LONGXY'], all_lons[idx]))
    #            met['LATIXY']=np.hstack((met['LATIXY'], all_lats[idx]))
    #return zone, zline, vdims, met
    
    
    mdoy=[0,31,59,90,120,151,181,212,243,273,304,334]#monthly starting DOY
    for iv in range(len(vars)):
        #if('GSWP3' in met_type or 'Site' in met_type or 'ESM' in met_type):
            v = vars[iv]
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
    
#------ ------ read a met variables from CLM full met directory 
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
    iy=ij[0]
    if len(all_lats.shape)>1: 
        ix=ij[1]
    else:
        ix=[0]

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
            for dirfile in dirfiles:
                # in metdata directory
                if(os.path.isfile(dirfile)):
                    fnc = Dataset(dirfile,'r')

                    # only need to read once for all files
                    if ('LATIXY' not in met.keys()):
                        if 'LATIXY' in fnc.variables.keys():
                            met['LATIXY']=np.asarray(fnc.variables['LATIXY'])[iy,ix]
                        elif 'lat' in fnc.variables.keys():
                            met['LATIXY']=np.asarray(fnc.variables['lat'])[iy,ix]
                    if ('LONGXY' not in met.keys()):
                        if 'LONGXY' in fnc.variables.keys():
                            met['LONGXY']=np.asarray(fnc.variables['LONGXY'])[iy,ix]
                        elif 'lon' in fnc.variables.keys():
                            met['LONGXY']=np.asarray(fnc.variables['lon'])[iy,ix]
        
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
                            vdata=np.asarray(fnc.variables[v])[:,iy,ix]
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
        
    return ix,iy, vdims,met
    
####################################################################################################

#------ ------ Read a single ELM met NC file user-defined 

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
    
#
# -------Read 1 site met data, *.csv, TOTALLY user-defined ----------------
#

def singleReadCsvfile(filename):
    
    
    # headers to be extracted from 'data_header', for ELM offline forcing 
    
    # Toolik '5_minute_data_v21.csv'
    #std_header = ('date', 'hour', \
    #              'air_temp_3m_5min', 'air_temp_5m_5min', \
    #              'lw_dn_avg_5min', \
    #              'sw_dn_avg_5min', \
    #              'rh_3m_5min', 'rh_5m_5min', \
    #              'wind_sp_5m_5min')
    
    # Toolik '1_hour_data.csv'
    std_header = ('date', 'hour', \
                  'air_temp_1m', 'air_temp_3m', 'air_temp_5m', \
                  'lw_dn_avg', \
                  'sw_dn_avg', \
                  'rain', 'rain2', \
                  'rh_1m','rh_3m', 'rh_5m', \
                  'wind_sp_1m','wind_sp_5m')
    
    # WERC-ATLAS data in SP,AK
    #std_header = ('date', 'hour', \
    #'Air Temperature @ 1 m (C)', 'Air Temperature @ 3 m (C)', \
    #'Relative Humidity @ 1 m (%)', 'Relative Humidity @ 3 m (%)', \
    #'Wind Speed (m/s)', \
    #'Wind Direction (degrees from true North)', 'standard deviation of Wind Direction', \
    #'Rainfall(mm)')

    
    missing ='nan'
    #
    try:
        #f0 = dd.read_csv(filename, dtype={'hour': 'float'})
        
        indx_data =[]
        key_data=[]

        data_header = f0.columns.values
        for j in data_header:
            if (len(j.strip())>0)  and (j.strip() in std_header):
                jidx = std_header.index(j.strip())
                indx = np.argwhere(data_header==j)[0]
                if (indx[0]>=0): 
                    indx_data.append(indx[0])
                    j_std = std_header[jidx] # the standard header corresponding to 'data_header'
                    key_data.append(j_std)
        
        data_value = {}
        for i,indx in enumerate(indx_data):
            key = key_data[i]
            val = np.asarray(f0[key])
            if (not key in data_value.keys()):
                data_value[key] = val
            else:
                data_value[key] = np.append(data_value[key],val)
                

    except ValueError:
        print("Error in reading data")
    #

    #
    del f0
    
    # timing coversion
    ymd=[int(s.split('/')[0]) for s in data_value['date']]
    data_value['MONTH']=np.asarray(ymd)
    ymd=[int(s.split('/')[1]) for s in data_value['date']]
    data_value['DAY']=np.asarray(ymd)
    ymd=[int(s.split('/')[2]) for s in data_value['date']]
    data_value['YEAR']=np.asarray(ymd)
    data_value['HOUR']=np.asarray([math.floor(t/100.0) for t in data_value['hour']])
    data_value['MM']=np.asarray([math.fmod(t,100.0) for t in data_value['hour']])
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
    hrly=3  # 3-hourly data aggregation
    for yr in range(yrmin,yrmax):
        days = 365
        if(isleap(yr)): days=366
        for doy in range(1,days+1):
            for hour in range(0,24,hrly):
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

    # standard varlists from ELM offline forcing
    # ['DTIME'], ['tunit']='days since 1901-01-01 00:00:00'
    #clm_header = ('DTIME','tunit', \
    #              'FLDS','FSDS','PRECTmms','PSRF','QBOT','TBOT','WIND')

    clm_data = {}
    clm_data['YEAR'] = tdata['YEAR']
    clm_data['MONTH'] = tdata['MONTH']
    clm_data['DAY'] = tdata['DAY']
    clm_data['HOUR'] = tdata['HOUR']
    clm_data['DOY'] = tdata['DOY']
    #clm_data['TBOT'] = np.nanmean([tdata['Air Temperature @ 1 m (C)'], \
                       #            tdata['Air Temperature @ 3 m (C)']], axis=0)+273.15
    clm_data['TBOT'] = np.nanmean([tdata['air_temp_1m'], tdata['air_temp_3m'], \
                                    tdata['air_temp_5m']], axis=0)+273.15

    clm_data['FLDS'] = np.empty_like(clm_data['TBOT']); clm_data['FLDS']= tdata['lw_dn_avg']
    clm_data['FSDS'] = np.empty_like(clm_data['TBOT']); clm_data['FSDS']=tdata['sw_dn_avg']
    
    #clm_data['PRECTmms'] = tdata['Rainfall(mm)']
    clm_data['PRECTmms'] = np.nanmean([tdata['rain'], tdata['rain2']], axis=0)
    
    clm_data['PSRF'] = np.empty_like(clm_data['TBOT']); clm_data['PSRF'][:]=101325.0
    
    #clm_data['RH'] = np.nanmean([tdata['Relative Humidity @ 1 m (%)'], \
    #                             tdata['Relative Humidity @ 3 m (%)']], axis=0)
    clm_data['RH'] = np.nanmean([tdata['rh_1m'], tdata['rh_3m'], \
                                 tdata['rh_5m']], axis=0)

    #clm_data['WIND'] = tdata['Wind Speed (m/s)']
    clm_data['WIND'] = np.nanmean([tdata['wind_sp_1m'],tdata['wind_sp_5m']], axis=0)

    clm_header = clm_data.keys()

    return clm_header, clm_data


#
# ------- CRU-JRA met forcing data convertion-extraction for CLM/ELM ------------------------------------
#
def clm_metdata_CRUJRA(crujra_dirfilehead, template_clm_metfile=''):
    # CRU-JRA data formats
    #  file naming (a full year, 6-hourly): 
          #crujra.v2.3.5d.$VAR.yyyy.365d.noc.nc
          #   where, $VAR - tmp, spfh, pre, pres, dswrf, dlwrf, ugrid, vgrd
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
                        unlimited_dim=''
                        unlimited_size=0
                        for dname, dimension in fmetdata_tmpl.dimensions.items():
                            if dname in fmetdata.dimensions:
                                if dname=='time':
                                    ldim = tindx[0].size
                                else:
                                    ldim = fmetdata.dimensions[dname].size
                                
                                dst.createDimension(dname, (ldim if not dimension.isunlimited() else None))
                                # unlimited dimension name to be included dst output
                                if dimension.isunlimited():
                                    unlimited_dim = dname
                                    unlimited_size = ldim

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
                            
                            x = dst.createVariable(vname, vdtype, \
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
#------- extract multiple point cpl_bypass dataset, using a list of points: point name, point lon, point lat
#
def multiple_cplbypass_extraction(fsites):
    #fsites = 'README' (txt file includes: site_name lat lon)
    lats=[]; lons=[]
    with open(fsites) as f:
        dtxt=f.readlines()
        
        dtxt=filter(lambda x: x.strip(), dtxt)
        for d in dtxt:
            all=np.array(d.split()[1:3],dtype=float)
            if(all[1]<0.0): all[1]=360.0+all[1] # convert longitude in format of 0 - 360
            lons.append(all[1])
            lats.append(all[0])
    f.close()
    lons = np.asarray(lons)
    lats = np.asarray(lats)
    for i in range(len(lats)):
        #clm_metdata_cplbypass_extraction('/Users/f9y/mygithub/pt-e3sm-inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v2.c180716/cpl_bypass_full/', \
        #clm_metdata_cplbypass_extraction('/lustre/or-scratch/cades-ccsi/proj-shared/project_acme/e3sm_inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v2.c180716/cpl_bypass_full/', \
        #                                  'GSWP3', \
        #                                  ncopath='/software/user_tools/current/cades-ccsi/nco/nco-5.1/bin/', z=1, l=i+1)
        clm_metdata_cplbypass_extraction('/home/fmyuan/mydata/unstructured_permafrost/', \
                                          'domain', \
                                          lons[i], lats[i], \
                                          ncopath='/usr/bin/', z=1, l=i+1)
        subfnc = glob.glob("%s*.%s" % ('./subset/', 'nc'))
        for ifile in subfnc:
            ncoutfile = ifile.split('/')[-1]
            tmpnc = 'tmp_'+ncoutfile
            if i==0:
                os.system('mv '+ifile+' ./'+ncoutfile) 
            else:
                if 'domain' in ifile:
                    ncmod.mergefilesby1dim(tmpnc, ifile, ncoutfile, 'ni')
                elif 'surfdata' in ifile:
                    ncmod.mergefilesby1dim(tmpnc, ifile, ncoutfile, 'gridcell')
                else:
                    ncmod.mergefilesby1dim(tmpnc, ifile, ncoutfile, 'n')
                
            os.system('cp '+ncoutfile+' '+tmpnc)
        
        if i==0:
            os.system('mv ./subset/zone_mappings.txt ./')
        else:
            with open("zone_mappings.txt", "a") as myfile:
                with open("./subset/zone_mappings.txt") as f: apptxt=f.readlines()[0]
                myfile.write(apptxt)
    
    os.system('rm -r ./subset')
    os.system('rm ./tmp_*.nc')


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


#clm_metdata_cplbypass_extraction('./', 'GSWP3_daymet4', 203.1241, 70.5725,ncopath='/usr/local/nco/bin/') #BEO
#clm_metdata_cplbypass_extraction('./', 'GSWP3_daymet4', -157.4089, 70.4696,ncopath='/usr/local/nco/bin/')  #ATQ
##clm_metdata_cplbypass_extraction('./', 'CRUJRAV2.3.c2023.0.5x0.5', -97.0287, 27.9798, ncopath='/software/user_tools/current/cades-ccsi/nco/nco-5.1/bin/')  #test
#multiple_cplbypass_extraction('info_14sites.txt')

#clm_metdata_cplbypass_extraction('./', 'GSWP3_daymet4', -80.25611, 34.97139, ncopath='/usr/local/gcc-x/nco_pacakge/nco-5.2.x/bin/')


