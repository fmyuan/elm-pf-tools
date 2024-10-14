#!/usr/bin/env python
import os, math
import numpy as np
from netCDF4 import Dataset
import glob
from pytools.commons_utils import Modules_netcdf

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
                allgrds=np.array(d.split(),dtype=float)
                if(allgrds[0]<0.0): allgrds[0]=360.0+allgrds[0] # convert longitude in format of 0 - 360
                all_lons.append(allgrds[0])
                all_lats.append(allgrds[1])
                if (not ('domain' in met_type or 'surfdata' in met_type)):
                    all_zones.append(int(allgrds[2]))
                    all_lines.append(int(allgrds[3]))
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
# ------- meteorological forcing data extraction for site or sites ------------------------------------
#
def clm_metdata_extraction(metdomainfile, metfiles, sites, ncopath=''):
    # ---- sites
    if (len(sites)<2):
        print('sites must have paired location points: x/y or longitude/latidue')
        return
    else:
        sitex = np.asarray(sites[0]) # x or longitudes
        sitey = np.asarray(sites[1]) # y or latitudes
        if(len(sites)>2): 
            site = sites[2]   # site name if any
        else:
            site='NEW'
    
    # ---- meteorological data domain
    domain = metdomainfile.rsplit('/')
    filelength=len(domain)
    if(filelength>1):
        domain_dir=domain[0]+'/'
    else:
        domain_dir='./'   
    domain_file = domain[-1]
    
    #
    ncfile = metdomainfile
    try:
        f = Dataset(ncfile,'r')
        print('\n FILE: '+ncfile+' ------- ')
    except Exception as e:
        print(e)
            
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
    
    
    pt_name = str(numxpts)+'x'+str(numypts)+'pt_'+site
    domaindir_new = './'+ pt_name+'/'
    os.system('mkdir -p ' + domaindir_new)

    #---------------------------------------------------------------------------------------------------------
    print('Extracting domain data for: Site - ', sitex, sitey)
    print('Extracted grid: ni,nj - ', ni, nj, 'lon,lat -', allx[ni], ally[nj])

    domainfile_new = domaindir_new+domain_file
    
    
    if (os.path.isfile(domainfile_new)):
        print('Warning:  Removing existing domain file')
        os.system('rm -rf '+domainfile_new)
    
    os.system(ncopath+'ncks -d ni,'+str(ni)+','+str(ni+numxpts-1)+ \
                  ' -d nj,'+str(nj)+','+str(nj+numypts-1)+ \
              ' '+metdomainfile+' '+domainfile_new)
    
    
    #--------------------------------------------------------------------------------------------------------
    #
    print('Extracting met-forcing data for: Site - ', site, sitex, sitey)
    
    met = metfiles.rsplit('/')
    filelength=len(met)
    if(filelength>1):
        met_dir=met[0]+'/'
    else:
        met_dir='./'   
    metfilehead = met[-1]
    if(metfilehead.endswith('.nc')):
        metfilehead=metfilehead.replace('.nc','') # removal of suffix of .nc
 
    pt_name = str(numxpts)+'x'+str(numypts)+'pt_'+site
    metdir_new = './'+ pt_name
    os.system('mkdir -p ' + metdir_new)
     
    #checking where is the original data
    dirfiles = os.listdir(met_dir)
    for dirfile in dirfiles:        
        filehead = metfilehead
        # in metdata directory
        if(os.path.isfile(met_dir+'/'+dirfile)): 
            if(filehead in dirfile):
            
                print('\n file: '+dirfile)
    
                metfile_old = met_dir+'/'+dirfile
                dirfile_new = dirfile
                metfile_new = metdir_new+'/'+ dirfile_new
    
                #extracting data
                print('dirfile: '+metfile_old + '  =======>  '+metfile_new)
    
                os.system(ncopath+'ncks --no_abc -O -d lon,'+str(ni)+','+str(ni+numxpts-1)+ \
                          ' -d lat,'+str(nj)+','+str(nj+numypts-1)+ \
                          ' '+metfile_old+' '+metfile_new)
        # in subdirectory of metdata directory
        elif(os.path.isdir(met_dir+'/'+dirfile)):
            subfiles = os.listdir(met_dir+'/'+dirfile)
            for subfile in subfiles:
    
                if(os.path.isfile(met_dir+'/'+dirfile+'/'+subfile)):
                    if(filehead in subfile):
     
                        metfile_old = met_dir+'/'+dirfile+'/'+subfile
                        
                        subfile_new = subfile
                        if(not os.path.isdir(metdir_new+'/'+dirfile)):
                            os.system('mkdir -p ' + metdir_new+'/'+dirfile)
                        metfile_new = metdir_new+'/'+dirfile+'/'+subfile_new
     
                        #extracting data
                        print('Subfile: '+metfile_old + '  =======>  '+metfile_new)
                    
                        os.system(ncopath+'ncks --no_abc -O -d lon,'+str(ni)+','+str(ni+numxpts-1)+ \
                                  ' -d lat,'+str(nj)+','+str(nj+numypts-1)+ \
                                  ' '+metfile_old+' '+metfile_new)
        
    print('DONE!')

####################################################################################################
#
# -------multiple sites cplbypass data extraction ----------------
#

def multiple_cplbypass_extraction(fsites):
    #fsites = 'README' (txt file includes: site_name lat lon)
    lats=[]; lons=[]
    with open(fsites) as f:
        dtxt=f.readlines()
        
        dtxt=filter(lambda x: x.strip(), dtxt)
        for d in dtxt:
            allgrds=np.array(d.split()[1:3],dtype=float)
            if(allgrds[1]<0.0): allgrds[1]=360.0+allgrds[1] # convert longitude in format of 0 - 360
            lons.append(allgrds[1])
            lats.append(allgrds[0])
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

#clm_metdata_cplbypass_extraction('./', 'GSWP3_daymet4', 203.1241, 70.5725,ncopath='/usr/local/nco/bin/') #BEO
#clm_metdata_cplbypass_extraction('./', 'GSWP3_daymet4', -157.4089, 70.4696,ncopath='/usr/local/nco/bin/')  #ATQ
##clm_metdata_cplbypass_extraction('./', 'CRUJRAV2.3.c2023.0.5x0.5', -97.0287, 27.9798, ncopath='/software/user_tools/current/cades-ccsi/nco/nco-5.1/bin/')  #test
##multiple_cplbypass_extraction('info_14sites.txt')

#clm_metdata_extraction('../domain.T62.050609.nc', './Solar6Hrly/clmforc.Qian.c2006.T62.Solr', [[267.0228], [40.6878]], ncopath='/usr/local/gcc-x/nco_pacakge/nco-5.2.x/bin/')

clm_metdata_extraction('./domain_42_FLUXNETSITES_simyr1850_c170912.nc', './42_FLUXNETSITES/20', [[5.9981], [50.3051]], ncopath='/usr/local/gcc-x/nco_pacakge/nco-5.2.x/bin/')
