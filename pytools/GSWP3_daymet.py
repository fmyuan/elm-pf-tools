#!/usr/bin/env python
import os, sys, csv, time, math
from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset
from pip._vendor.distlib.util import CSVReader
from numbers import Real
from token import ISEOF

# options
parser = OptionParser()
parser.add_option("--casedataroot", dest="casedataroot", default="./")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='./', \
                  help = "input data directory for CESM (required)")
parser.add_option("--metdomainfile", dest="metdomainfile", default="domain.360x720_ORCHIDEE0to360.100409.nc", \
		          help = 'met forcing data domain file, usually under atm/datm7, either under met_directory or domain.clm - so must be indicated with file name')
parser.add_option("--metdir", dest="metdir", default="atm_forcing.datm7.cruncep.0.5d.V4.c111104", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--metfilehead", dest="metfilehead", default="clmforc.cruncep.V4.c2011.0.5d.", \
          help = 'met forcing data file name header, i.e. that without date-time portion')
parser.add_option("--lons", dest="long", default="", \
                      help = "point longitude(s)")
parser.add_option("--lats", dest="lati", default="", \
                      help = "point latitude(s)")
parser.add_option("--yrstart", dest="yrstart", default=1901, \
                      help = "climate data starting year")
parser.add_option("--yrend", dest="yrend", default=2010, \
                      help = "climate data ending year")
parser.add_option("--sitename", dest="site", default="", \
                      help = "point site-name")
parser.add_option("--ncobinpath", dest="ncobinpath", default="", \
                      help = "NCO bin path if not in $PATH")

# ---------------------------------------------------------
(options, args) = parser.parse_args()

casedatadir = os.path.abspath(options.casedataroot)
ccsm_input  = os.path.abspath(options.ccsm_input)
if(options.ccsm_input=='./'):
    met_input   = os.path.abspath(options.ccsm_input+options.metdir)
    metdomain_input   = os.path.abspath(options.ccsm_input)
else:
    met_input   = os.path.abspath(options.ccsm_input+"/atm/datm7/"+options.metdir)
    metdomain_input = os.path.abspath(options.ccsm_input+"/atm/datm7/")

#--------------- nco bin path ------
if(options.ncobinpath==""):
    ncopath = options.ncobinpath
else:
    ncopath = options.ncobinpath+"/"


#------------------- get site information ----------------------------------

scriptdir = os.getcwd()


#####################################################################################################
#
# ------- GSWP3v1 met forcing data extraction -------------------------------------------------
#
def clm_metdata_extraction():
    # met domain
    domainfile_orig = metdomain_input+"/"+options.metdomainfile
    #
    ncfile = domainfile_orig
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
    print('Extracting domain data for: Site - ', options.site, sitex, sitey)
    print('Extracted grid: ni,nj - ', ni, nj, 'lon,lat -', allx[ni], ally[nj])
    if(options.ccsm_input=='./'):
        os.system('mkdir -p '+casedatadir)
        domainfile_new = casedatadir+'/' \
             +'domain.lnd.'+str(numxpts)+'x'+str(numypts)+'pt_'+options.site+'_navy.nc'
    else:
        os.system('mkdir -p '+casedatadir+'/atm/datm7/domain.clm')
        domainfile_new = casedatadir+'/atm/datm7/domain.clm/' \
             +'domain.lnd.'+str(numxpts)+'x'+str(numypts)+'pt_'+options.site+'_navy.nc'
    
    
    if (os.path.isfile(domainfile_new)):
        print('Warning:  Removing existing domain file')
        os.system('rm -rf '+domainfile_new)
    
    os.system(ncopath+'ncks -d ni,'+str(ni)+','+str(ni+numxpts-1)+ \
                  ' -d nj,'+str(nj)+','+str(nj+numypts-1)+ \
              ' '+domainfile_orig+' '+domainfile_new)
    
    
    #--------------------------------------------------------------------------------------------------------
    #
    print('Extracting met-forcing data for: Site - ', options.site, sitex, sitey)
    
    pt_name = str(numxpts)+'x'+str(numypts)+'pt_'+options.site
    met_input_new = casedatadir+'/atm/datm7/CLM1PT_data/'+ pt_name
    os.system('mkdir -p ' + met_input_new)
    
    #checking where is the original data
    dirfiles = os.listdir(met_input)
    for dirfile in dirfiles:        
        filehead = options.metfilehead
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
# -------------------Read multiple location DAYMET data, *.csv one by one ----------------
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
    
##################################################################################

site,odata_header,odata = \
    singleDaymetReadCsvfile('daymet_kougarok-NGEE00.csv')


    
    
    



