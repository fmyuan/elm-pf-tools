#!/usr/bin/env python

import sys
import glob
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Modules_CLM_nc4 import *

# ---------------------------------------------------------------
# Plot soil data with layers
def SoilLayeredVarPlotting(plt, nrow, ncol, isubplot, t, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    if(varname == 'SOILPSI'):
        sdata = -sdata
        ax.set_yscale("log", nonposy='clip')
    
    plt.plot(t, sdata[:,0], 'b-.', t, sdata[:,1], 'g-.', \
                  t, sdata[:,2], 'm-.', t, sdata[:,3], 'c-.', \
                  t, sdata[:,4], 'r-.', \
                  t, sdata[:,5], 'b-', t, sdata[:,6], 'g-', \
                  t, sdata[:,7], 'm-', t, sdata[:,8], 'c-', \
                  t, sdata[:,9], 'r-', \
                  t, sdata[:,10], 'b--', \
                  t, sdata[:,11], 'g--',  t, sdata[:,12], 'm--', \
                  t, sdata[:,13], 'c--',  t, sdata[:,14], 'r--')
    plt.legend(('Layer 15','Layer 14','Layer 13','Layer 12','Layer 11', \
            'Layer 10','Layer 9','Layer 8','Layer 7','Layer 6', \
            'Layer 5','Layer 4','Layer 3','Layer 2','Layer 1'), \
            loc=0, fontsize=12)
    plt.xlabel('Years')
    plt.ylabel(varname)

    lx = 0.40
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

# Plot Grided data
def GridedVarPlotting(plt, nrow, ncol, isubplot, t, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    
    plt.plot(t, sdata)
    
    #plt.legend(loc=0, fontsize=12)
    plt.xlabel('Years')
    plt.ylabel(varname)

    lx = 0.40
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--clmout_dir", dest="clm_odir", default="./", \
                  help="clm output directory (default = ./, i.e., under current directory)")
parser.add_option("--clmfile_head", dest="ncfileheader", default="", \
                  help = "clm output file name header, usually the portion before *.clm2.h[0-5].*.nc")
parser.add_option("--varname_help", dest="varnames_print", default=False, \
                  help = "print out VARIABLES available ")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4) to be reading/plotting, separated by comma ")
parser.add_option("--adspinup", dest="adspinup", action="store_true", default=False, \
                  help="whether results of an ad_spinup run (default = False)")
parser.add_option("--startyr", dest="startyr", default="1", \
                  help="clm run starting year (default = 1, this is for spinup; for transient it should be 1850; " \
                   " and can be user-defined)")
parser.add_option("--endyr", dest="endyr", default="", \
                  help="clm run ending year (default = none, i.e. end of simulation)")
parser.add_option("--Xindex", dest="xindex", default=0, \
                  help = " X direction grid index to be reading/plotting, default 0 ")
parser.add_option("--Yindex", dest="yindex", default=0, \
                  help = " Y direction grid index to be reading/plotting, default 0 ")

(options, args) = parser.parse_args()

#
if (options.clm_odir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.clm_odir)

if (options.ncfileheader == ''):
    print('MUST input clm file header by " --ncfileheader=??? "' //
     ' , which usually is the portion before "*.clm2.[h0-h5].*.nc"! ')
    sys.exit()
else:
    print('clm nc file header: '+ options.ncfileheader)

if (options.vars == ''):
    print('No variable name by " --varname=??? "; So print out ALL variable names')
    varnames =[]
    options.varnames_print = True
else:
    varnames = options.vars.split(':')  
    nvars = len(varnames)
    if(nvars>4 or options.vars=="ALL"):
        print('Currently ONLY support at most 4 variables to be plotted')
        sys.exit()


starttime = int(options.startyr)
if(options.endyr !=""): endtime = int(options.endyr)

ix=int(options.xindex);
iy=int(options.yindex);

# read-in datasets from 1 simulation
tt, nx, ny, nlgrnd, nldcmp, npft, varsdata, varsdims = \
    CLM_NcRead_1simulation(options.clm_odir, \
                           options.ncfileheader, \
                           options.varnames_print, \
                           varnames, \
                           starttime, endtime, \
                           options.adspinup)


#--------------------------------------------------------------------------------------
 
    

# data-sets will be sorted by time-series
t  = sorted(tt)
it = sorted(range(len(tt)), key=lambda k: tt[k])
nt = len(tt)

# dimension max.
if2dgrid = True
if(len(varsdata['topo'].shape)==1): if2dgrid = False
nxy = nx*ny
nl = max(nlgrnd, nldcmp)

# plotting

nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
if(nvars>=2):
    nrow = 2
if(nvars>=3):
    ncol = 2

sdata = np.zeros((nt,nl*nx*ny)) #temporary data holder in 2-D (tt, layers*grids)
gdata = np.zeros((nt,nx*ny))    #temporary data holder in 2-D (tt, grids)

ivar = 0
for var in varnames:
    # plot 1
    vdata = varsdata[var]
    vdims = varsdims[var]
    
    #time dimension
    dim_tt = -999
    if('time' in vdims): 
        dim_tt = vdims.index('time')
        ivar = ivar + 1 # counting vars, with time-series only
    else:
        continue
    
    # Does have vertical dimension?  
    zdim_indx = -999
    if('levgrnd' in vdims): 
        zdim_indx = vdims.index('levgrnd')
        nl = vdims[zdim_indx]
    elif('levdcmp' in vdims): 
        zdim_indx = vdims.index('levdcmp')
        nl = vdims[zdim_indx]

    # pft ? (need further work here)
    if(npft>0):
        pdim_indx = vdims.index('pft')  
    
    # data series
    for i in range(len(tt)):
                
        if(zdim_indx<0):# 2-D grid data
            gdata[i,:] = vdata[it[i]].reshape(nx*ny)
        else:
                       
            if(zdim_indx == 1): # 3-D soil data, z_dim in 1 
                if(if2dgrid): 
                    sdata[i,:] = vdata[it[i]][:,iy,ix]
                else:
                    sdata[i,:] = vdata[it[i]][:,max(iy,ix)]
                
            elif(zdim_indx >= 2): # 3-D soil data, z_dim in 2 or 3 (likely x/y dims before z) 
                if(if2dgrid): 
                    sdata[i,:] = vdata[it[i]][:,iy,ix,]
                else:
                    sdata[i,:] = vdata[it[i]][:,max(iy,ix),]
        
    #plotting
    vname = varnames[varnames.index(var)]
    if(zdim_indx<0):
        GridedVarPlotting(plt, nrow, ncol, ivar, t, gdata, vname, '(a) All Grids')
    else:
        SoilLayeredVarPlotting(plt, nrow, ncol, ivar, t, sdata, vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')
    


# printing plot in PDF
ofname = 'Figure_CLM.pdf'
fig = plt.gcf()
fig.set_size_inches(12, 15)
plt.savefig(ofname)
plt.show()

plt.close('all')


