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


###################################################################################

site,odata_header,odata = \
    singleDaymetReadCsvfile('daymet_kougarok-NGEE00.csv')


    
    
    



