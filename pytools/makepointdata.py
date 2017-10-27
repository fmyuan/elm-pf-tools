#!/usr/bin/env python
import os, sys, csv, time, math
from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset

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
parser.add_option("--metfile", dest="metfile", default="clmforc.cruncep.V4.c2011.0.5d.", \
          help = 'met forcing data file name header, i.e. that without date-time portion')
parser.add_option("--lons", dest="long", default="", \
                      help = "number of x points (regional only)")
parser.add_option("--lats", dest="lati", default="", \
                      help = "number of y points (regional only)")
parser.add_option("--yrstart", dest="yrstart", default=1901, \
                      help = "number of y points (regional only)")
parser.add_option("--yrend", dest="yrend", default=2010, \
                      help = "number of y points (regional only)")
parser.add_option("--sitename", dest="site", default="", \
                      help = "number of y points (regional only)")
parser.add_option("--ncobinpath", dest="ncobinpath", default="", \
                      help = "number of y points (regional only)")

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

x=options.long.split(":")
sitex=np.array(x,dtype='float')
y=options.lati.split(":")
sitey=np.array(y,dtype='float')

# ------- met forcing data extraction ------

# met domain
domainfile_orig = metdomain_input+"/"+options.metdomainfile
#
ncfile = domainfile_orig
try:
    f = Dataset(ncfile,'r')
    print('FILE: '+ncfile+' ------- ')
except:
    print('Error in READING File: '+ncfile)
        
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
    if(isfile(met_input+'/'+dirfile)): 
        if(options.metfile in dirfile):
            metfile_old = met_input+'/'+dirfile
            metfile_new = metfile_old.replace(metfile,pt_name)

            os.system(ncopath+'ncks -d ni,'+str(ni)+','+str(ni+numxpts-1)+ \
                      ' -d nj,'+str(nj)+','+str(nj+numypts-1)+ \
                      ' '+metfile_old+' '+metfile_new)
        
    elif(isdir(met_input+'/'+dirfile)):
        subfiles = os.listdir(met_input+'/'+dirfile)
        for subfile in subfiles:
            if(isfile(met_input+'/'+dirfile+'/'+subfile)):
                if(options.metfile in dirfile):
                    metfile_old = met_input+'/'+dirfile
                    metfile_new = metfile_old.replace(metfile,pt_name)
                 
                    os.system(ncopath+'ncks -d ni,'+str(ni)+','+str(ni+numxpts-1)+ \
                              ' -d nj,'+str(nj)+','+str(nj+numypts-1)+ \
                              ' '+metfile_old+' '+metfile_new)
    


