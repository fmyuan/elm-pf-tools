#!/usr/bin/env python

#----------------------
# Concatentating ELM monthly/yearly output nc files, for a list of variables
# with each in a single nc file.
# The dafault is to do so for all variables in the first file.
#

import os, sys, time, math
import numpy as np
from optparse import OptionParser
import re
from netCDF4 import Dataset

# the nc file variable(s) (not dimension variable) to be excluded for concatentation
excl_vars = ['mcdate','mcsec','mdcur','mscur','nstep','date_written','time_written','time','time_bounds']

# the nc file variable(s) always to be included together with specific variable for concatentation, if existed
incl_vars = ['lat','lon','mcdate']

#-------------------Parse options-----------------------------------------------
parser = OptionParser()

parser.add_option("--outdir", dest="outdir", default="", \
                  help = 'output directory')
parser.add_option("--ofheader", dest="header", default="none", \
                  help = 'header of output files either in fullpath or filename in default path')
parser.add_option("--varlist", dest="varlist", default="", \
                      help = "variable list, separated by ':,;' default is all variables")
parser.add_option("--varlisttxt", dest="varlisttxt", default="", \
                      help = "variable list txt file")
parser.add_option("--yrstart", dest="yrstart", default="", \
                      help = "user-defined start year")
parser.add_option("--yrend", dest="yrend", default="", \
                      help = "user-defined end year")
parser.add_option("--timer_type", dest="t_type", default="time", \
                      help = "either 'time' or 'date'")
parser.add_option("--ncobinpath", dest="ncobinpath", default="", \
                      help = "NCO bin path if not in $PATH")


(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------
currdir = os.getcwd()

#--------------- nco bin path -----------
if(not options.ncobinpath.strip()==''):
    ncopath = options.ncobinpath
else:
    ncopath = '/usr/local/nco/bin'
if(not ncopath.endswith('/')): ncopath = ncopath.strip()+'/'

###########################################################################
if(options.outdir.startswith('/')):
    outdir = os.path.abspath(options.outdir)
    if(not outdir.endswith('/')): outdir = outdir+'/'
else:
    outdir = './'

#------------get all files if it's what needed -----------------------------------
filehead = options.header
if (os.path.isdir(outdir.strip())):
    #checking where is the original data
    dirfiles = os.listdir(outdir)
    ncfiles = []
    for dirfile in dirfiles:                
        if(dirfile.startswith(filehead)): 
            ncfiles.append(dirfile)

    #
    ncfiles.sort(key = lambda x: x.split('.')[-2])
                
yrstart = ncfiles[0].split('.')[-2]
yrend = ncfiles[-1].split('.')[-2]

varlist = [] # this actually will concatentate all variables, except for 'excl_vars'
if(options.varlisttxt!=''):  # varlist txt file is checked first, so over-ride 'options.varlist'
    ftxt = open(options.varlisttxt,'r')
    varlist = ftxt.read().splitlines()
    ftxt.close()    
elif(options.varlist!=''):
    varlist = re.split('[:,;]',options.varlist)
    
#------------extract all vars one by one as needed -----------------------------------
vars_name={}
vars_dims={}
for nf in ncfiles:            
    print('\n file: '+nf)
    
    outfile_old = outdir+'/'+nf

    #-------------------------------------------------------------
    # (1) get time-varying variable names and their dimensions ONLY once
    if (len(vars_name)==0):
        f = Dataset(outfile_old,'r')
        allvars = f.variables.keys()
        for key in allvars:                       
            if(('time' in f.variables[key].dimensions) and (key not in excl_vars)):
                if(key in varlist or len(varlist)==0):
                    vars_dims[key]=f.variables[key].dimensions
                    vars_name[key]=f.variables[key].name
        f.close()
                      
    #-------------------------------------------------------------
    # (2) 'time' not always a record dimension
    os.system(ncopath+'ncks -O -h --mk_rec_dim time '+ outfile_old + ' -o alltemp.nc') 
           
    #-------------------------------------------------------------
    # (3)
    for key in vars_name.keys():                       
        dirfile_new = 'CONCAT_'+filehead+'_'+vars_name[key]+'.'+yrstart+'_'+yrend+'.nc'        
        outfile_new = './'+dirfile_new
        if(os.path.isfile(outfile_new) and nf==ncfiles[0]): 
            os.system('rm -f '+outfile_new)
                            
        #nco command ' -v' options
        vstring = ' -v time'
        for v in incl_vars:
            if (v not in vstring and v in allvars): 
                vstring += (','+str(v).strip())
                    
        if(str(vars_name[key]) not in vstring): 
                vstring += (','+str(vars_name[key]).strip())
        # needs dimenson variable values, if in nc files
        for dimv in vars_dims[key]:
                if (dimv not in vstring and dimv in allvars): 
                    vstring += (','+str(dimv.strip()))

        #
        os.system(ncopath+'ncks -h -O '+ vstring+ \
                              ' alltemp.nc -o temp.nc')
                    
        #
        if(os.path.isfile(outfile_new)):
            #concatentating
            os.system(ncopath+'ncrcat -O -h '+outfile_new + ' temp.nc -o' +outfile_new)

            # after the final file done, switch 'time' with 'mcdate' if required
            # NOTE: 'time' is the day numbers since starting yyyy-mm-dd-time
            #       which are hard to know what real time is. Then better to use 'mcdate'
            if(nf==ncfiles[-1]):
                if(options.t_type != 'time'):
                    os.system(ncopath+'ncrename -O -h -v time,dtime '+outfile_new+ \
                              ' -o '+outfile_new)
                    os.system(ncopath+'ncrename -O -h -v mcdate,time '+outfile_new+ \
                              ' -o '+outfile_new)
        else:
            #copy the first one
            os.system('mv temp.nc ' +outfile_new)
                                         
    # (END) loop for all vars possibly in the nc file
    os.system('rm -f alltemp.nc')
    os.system('rm -f temp.nc')
                    
#(END) for-loop of all files in outdir directory

           
#------------END of CLM_outvars_ncocat -----------------------------------------------------------------------
    

