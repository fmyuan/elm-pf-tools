#!/usr/bin/env python

import sys
import glob
import re
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Modules_CLM_nc4 import CLM_NcRead_1simulation

# ---------------------------------------------------------------
# Plot soil data with layers for ONE specific grid
def SoilLayeredVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, layer_index, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    if(varname == 'SOILPSI'):
        sdata = -sdata
        ax.set_yscale("log", nonposy='clip')
    
    layertext = []
    for il in layer_index:
        layertext.append(("Layer "+str(il)))
        plt.plot(t, sdata[:,il])
 
    plt.legend((layertext), loc=0, fontsize=12)
    plt.xlabel(t_unit)
    plt.ylabel(varname)

    lx = 0.40
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

# Plot PFT fractioned data with layers for ONE specific grid
def PFTVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, pftwt, pft_index, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    
    active_pfts = []
    for ix in pft_index:
        for ip, wt in enumerate(pftwt):
            if(wt>0.0 and (ix<0 or ix==ip)): 
                active_pfts.append(("PFT "+str(ip)))
                plt.plot(t, sdata[:,ip])
    
    plt.legend((active_pfts), loc=0, fontsize=12)
    plt.xlabel(t_unit)
    plt.ylabel(varname)

    lx = 0.40
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)


# Plot Grided data
def GridedVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    
    plt.plot(t, sdata)
    
    #plt.legend(loc=0, fontsize=12)
    plt.xlabel(t_unit)
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
                  help = "print out VARIABLES available ", action="store_true")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4) to be reading/plotting, separated by comma ")
parser.add_option("--adspinup", dest="adspinup", action="store_true", default=False, \
                  help="whether results of an ad_spinup run (default = False)")
parser.add_option("--startyr", dest="startyr", default="1", \
                  help="clm run starting year (default = 1, this is for spinup; for transient it should be 1850; " \
                   " and can be user-defined)")
parser.add_option("--endyr", dest="endyr", default="", \
                  help="clm run ending year (default = none, i.e. end of simulation)")
parser.add_option("--plot_tunit", dest="t_unit", default="Days", \
                  help="X-axis time unit (default = Days, i.e. Days since start-time of simulation)")
parser.add_option("--Xindex", dest="xindex", default=0, \
                  help = " X direction grid index to be reading/plotting, default 0 ")
parser.add_option("--Yindex", dest="yindex", default=0, \
                  help = " Y direction grid index to be reading/plotting, default 0 ")
parser.add_option("--LAYERindex", dest="zindex", default=-999, \
                  help = " SOIL layer index to be reading/plotting, default -999 for all, with indexing from 0 ")
parser.add_option("--PFTindex", dest="pindex", default=-999, \
                  help = " PFT index to be reading/plotting, default -999 for all, with indexing from 0 ")

(options, args) = parser.parse_args()

if(options.pindex==-999):
    pft_index = options.pindex
else:
    pft_index=[]
    oppfts = options.pindex.split(':')
    for ip in oppfts:    
        pft_index.append(int(ip))

if(options.zindex==-999):
    layer_index = options.zindex
else:
    layer_index=[]
    oplayers = options.zindex.split(':')
    for il in oplayers:
        layer_index.append(int(il))

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

if (options.varnames_print):
    print('checking Varname lists')
    varnames = []

startdays = -9999
startdays = (int(options.startyr)-1)*365.0
enddays = -9999
if(options.endyr !=""): enddays = int(options.endyr)*365.0

tunit = str.capitalize(options.t_unit)
if tunit.startswith("H"): 
    day_scaling = 24.0
elif tunit.startswith("Y"): 
    day_scaling = 1.0/365.0
else:
    day_scaling = 1.0



ix=int(options.xindex);
iy=int(options.yindex);

# read-in datasets from 1 simulation
nx, ny, nlgrnd, nldcmp, ncol, npft, varsdata, varsdims = \
    CLM_NcRead_1simulation(options.clm_odir, \
                           options.ncfileheader, \
                           options.varnames_print, \
                           varnames, \
                           startdays, enddays, \
                           options.adspinup)

if (options.varnames_print): sys.exit('Variable Names Printing DONE!')
#--------------------------------------------------------------------------------------
# 
vars_list = list(varsdata.keys())    # var_list is in format of 'h*_varname', in which 'varname' is the real variable names

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

ivar = 0
for var in varnames:

    # names for variable and its time-axis
    for hv in vars_list:
        if re.search(var, hv): 
            var_h = hv
            hinc  = hv.replace(var,'')
            var_t = '%stime' %hinc
            break
        
    vdata = varsdata[var_h]
    vdims = varsdims[var_h]
    
    tt = varsdata[var_t]
    t  = sorted(tt)
    it = sorted(range(len(tt)), key=lambda k: tt[k])
    nt = len(tt)

    
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
        nl = vdata.shape[zdim_indx]
    elif('levdcmp' in vdims): 
        zdim_indx = vdims.index('levdcmp')
        nl = vdata.shape[zdim_indx]

    # pft ? (need further work here)
    pdim_indx = -999
    if('pft' in vdims):
        pdim_indx = vdims.index('pft')  
        npft = vdata.shape[pdim_indx]
        pwt1cell = varsdata['pfts1d_wtgcell']
        
    # data series
    if(zdim_indx<0 and pdim_indx<0):# 2-D grid data
        gdata = np.zeros((nt,nx*ny))    #temporary data holder in 2-D (tt, grids)
    elif(zdim_indx>0):        
        sdata = np.zeros((nt,nl*nx*ny)) #temporary data holder in 2-D (tt, layers*grids)
    elif(pdim_indx>0):        
        sdata = np.zeros((nt,npft*nx*ny)) #temporary data holder in 2-D (tt, pfts*grids)
    else:
        exit("Variable to be plotted has 4-D dimension - NOT YET supported!")
    
    for i in range(len(tt)):
                
        if(zdim_indx<0 and pdim_indx<0):# 2-D grid data
            gdata[i,:] = vdata[it[i]].reshape(nx*ny)
        else:
                       
            if(zdim_indx == 1 or pdim_indx==1): # 3-D soil/pft data, z_dim/p_dim in 1 
                if(if2dgrid): 
                    sdata[i,:] = vdata[it[i]][:,iy,ix]
                else:
                    if(len(vdata[it[i]].shape)>1):
                        sdata[i,:] = vdata[it[i]][:,max(iy,ix)]
                    else:
                        sdata[i,:] = vdata[it[i]][:]
                        
                
            elif(zdim_indx >= 2 or pdim_indx==2): # 3-D soil data, z_dim in 2 or 3 (likely x/y dims before z) 
                if(if2dgrid): 
                    sdata[i,:] = vdata[it[i]][:,iy,ix,]
                else:
                    sdata[i,:] = vdata[it[i]][:,max(iy,ix),]
        
    #plotting
    vname = varnames[varnames.index(var)]
    if(float(day_scaling)!=1.0): t = np.asarray(t)*day_scaling
    if(zdim_indx<0 and pdim_indx<0):
        GridedVarPlotting(plt, nrow, ncol, ivar, t, tunit, gdata, \
                        vname, '(a) All Grids')
    elif(zdim_indx>0):
        if(len(layer_index)==1 and layer_index[0]<0): layer_index = range(0,nl)
        SoilLayeredVarPlotting(plt, nrow, ncol, ivar, t, tunit, layer_index, sdata, \
                        vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')
    
    elif(pdim_indx>0):
        if(len(pft_index)==1 and pft_index[0]<0): pft_index = range(0,npft)
        PFTVarPlotting(plt, nrow, ncol, ivar, t, tunit, pwt1cell, pft_index, sdata, \
                        vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')


# printing plot in PDF
ofname = 'Figure_CLM.pdf'
fig = plt.gcf()
fig.set_size_inches(12, 15)
plt.savefig(ofname)
plt.show()

plt.close('all')


