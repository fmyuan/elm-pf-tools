#!/usr/bin/env python
import os, sys, csv, time, math
from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset
from pip._vendor.distlib.util import CSVReader
from numbers import Real
from token import ISEOF

from Modules_metdata import singleDaymetReadCsvfile
from Modules_metdata import clm_metdata_extraction
from Modules_metdata import subsetDaymetRead1NCfile

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
parser.add_option("--daymet_dir", dest="daymetdir", default="", \
                  help = 'directory for daymet data')
parser.add_option("--daymet_filehead", dest="daymetfilehead", default="daymet_v3", \
          help = 'daymet data file name header, i.e. that without year-region name, so can process batch files')
parser.add_option("--daymet_splitprcp", dest="split_prcp", default=False, \
          help = 'split daymet prcp into rain/snow, NOTE that must have both tmax and tmin files', \
          action="store_true")
parser.add_option("--lons", dest="long", default="", \
                      help = "point/range longitude(s)")
parser.add_option("--lats", dest="lati", default="", \
                      help = "point/range latitude(s)")
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



###################################################################################
if(options.daymetdir!=""):# daymet single variable directory
    daymet_input   = os.path.abspath(options.daymetdir)# top daymet directory
    if any(x in daymet_input for x in ['prcp','tmax','tmin']):
        daymet_var=daymet_input.split('/')[-1]
        daymet_input = daymet_input.replace(daymet_var,'')
    else:
        daymet_var = ''

    vardir = daymet_input+'/'+daymet_var
    dirfiles = sorted(os.listdir(vardir))
    filecounter = 0
    for dirfile in dirfiles:        
        filehead = options.daymetfilehead
        # in daymet directory
        if(os.path.isfile(vardir+'/'+dirfile)): 
            if(dirfile.startswith(filehead)):
                daymetfile = vardir+'/'+dirfile
                print('\n file: '+daymetfile)

                # single-point .csv daymet data
                if daymetfile.endswith('.csv'):
                    site,odata_header,odata = \
                        singleDaymetReadCsvfile(daymetfile)
                        #singleDaymetReadCsvfile('daymet_kougarok-NGEE00.csv')

                # multiple gridded daymet nc4 dataset
                if daymetfile.endswith('.nc4') or daymetfile.endswith('.nc'):
                    
                    #subsetting NC4 daymet data for most AK (>=60oN only)
                    #"name": "Alaska",
                    #"min_lat": 52.5964,
                    #"max_lat": 71.5232,
                    #"min_lng": -169.9146,
                    #"max_lng": -129.993
                    box_lon = []#[-170.0, -129.0] or [-170.0, -141.0]
                    if options.long != "":
                        box_lon = np.asarray(options.long.split(':'), dtype=np.float)
                        print ('subsetting Longitude Range: ', box_lon)
                    box_lat = []#[52.0, 90.0] or [60.0, 90.0]
                    if options.lati != "":
                        box_lat = np.asarray(options.lati.split(':'), dtype=np.float)
                        print ('subsetting Latitude Range: ', box_lat)
                    if (len(box_lon)>0 or len(box_lat)>0):
                        odata = \
                            subsetDaymetRead1NCfile(daymetfile, lon_range=box_lon, lat_range=box_lat, SUBSETNC=True, ncopath=ncopath)
                
                    # split precipitation into rain or snow
                    if options.split_prcp:
                        if 'prcp' in dirfile:
                            if (os.path.isdir(daymet_input+'/tmax/')):
                                tmaxfile = daymet_input+'/tmax/'+dirfile.replace('prcp','tmax')
                            else:
                                tmaxfile = daymet_input+'/'+dirfile.replace('prcp','tmax')
                            if(not os.path.isfile(tmaxfile)):
                                print ('NO tmax data file found: '+ tmaxfile)
                                print ('Usually in: '+daymet_input+'/tmax/')
                                print ('cannot split prcp into rain or snow')
                                sys.exit()
                            
                            if (os.path.isdir(daymet_input+'/tmin/')):
                                tminfile = daymet_input+'/tmin/'+dirfile.replace('prcp','tmin')
                            else:
                                tminfile = daymet_input+'/'+dirfile.replace('prcp','tmin')
                            if(not os.path.isfile(tminfile)):
                                print ('NO tmin data file found: '+ tminfile)
                                print ('Usually in: '+daymet_input+'/tmin/')
                                print ('cannot split prcp into rain or snow')
                                sys.exit()
                            
                            tempfile='temp.'+dirfile.split('.')[-1]
                            os.system('cp -rf '+daymetfile+' ./'+tempfile)
                            os.system(ncopath+'ncks -h -A '+tmaxfile+' -o '+tempfile)
                            os.system(ncopath+'ncks -h -A '+tminfile+' -o '+tempfile)
                            # in the following, '--cnk_dmn time,30' is to limit nco netcdf max. chunk size of 4GB, may need customized
                            os.system(ncopath+"ncap2 --cnk_dmn time,30 -O -s 'tave=(tmax+tmin)/2.0' " +tempfile+ " "+tempfile)
                            os.system(ncopath+'ncks -h -O -x -v tmax,tmin ' +tempfile+' temp2.nc')
                            os.system('mv temp2.nc'+' '+tempfile)
                            
                            # rain
                            outfile = dirfile.replace('prcp', 'rain')
                            os.system(ncopath+"ncap2 -O -s 'where(tave<0.0) {prcp=0; }' " +tempfile+ " temp2.nc")
                                #cannot put '<' and '==' together or '<=' in line above
                            os.system(ncopath+"ncap2 -O -s 'where(tave==0.0) {prcp=0; }' temp2.nc " +outfile)
                            os.system(ncopath+'ncks -h -O -x -v tave '+outfile + ' -o'+outfile)
                            os.system(ncopath+'ncrename -O -v prcp,rain '+outfile)
                            os.system('rm -rf temp2.nc')
                            
                            #snow
                            outfile = dirfile.replace('prcp', 'snow')
                            os.system(ncopath+"ncap2 -O -s 'where(tave>0.0) {prcp=0; }' " \
                                      +tempfile+ " " +outfile)
                            os.system(ncopath+'ncks -h -O -x -v tave '+outfile + ' -o'+outfile)
                            os.system(ncopath+'ncrename -O -v prcp,snow '+outfile)
                            os.system('rm -rf '+tempfile)
                            
                            
                            
                        else:
                            continue
                
                filecounter = filecounter + 1
            
    if (filecounter<=0):
        print('Warning: NO daymet nc file found in '+ daymet_input)

