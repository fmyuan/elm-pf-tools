#!/usr/bin/env python

import sys
import glob
import re
import math
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Modules_CLMoutput_nc4 import CLM_NcRead_1simulation
from Modules_CLMoutput_nc4 import CLMvar_1Dtseries


# ---------------------------------------------------------------
# Plot soil data with layers for ONE specific grid
def SoilLayeredVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, layer_index, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    if(varname == 'SOILPSI'):
        sdata = -sdata
        ax.set_yscale("log", nonposy='clip')
    
    layertext = []
    
    # add a zero-degree-C line for 'TSOI'
    if(varname == 'TSOI'):
        plt.plot([min(t),max(t)],[273.15,273.15],'k--')
        layertext.append("0oC ")

    for il in layer_index:
        layertext.append(("Layer "+str(il)))
        plt.plot(t, sdata[:,il])
 
    plt.legend((layertext), loc=0, fontsize=12)
    plt.xlabel(t_unit, fontsize=12, fontweight='bold')
    plt.ylabel(varname, fontsize=12, fontweight='bold')

    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

# Plot PFT fractioned data with layers for ONE specific grid
def PFTVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, pftwt, pftvidx, pft_index, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    
    active_pfts = []

    # add a zero line for a few variables
    if(varname in ['GPP','NPP','NEE','MR','GR']):
        plt.plot([min(t),max(t)],[0.0,0.0],'k--')
        active_pfts.append("0 line ")
    
    for ix in pft_index:
        for ip, wt in enumerate(pftwt):
            if(wt>0.0 and (pftvidx[ix]>0 and ix==ip)): 
                active_pfts.append(("PFT "+str(ip)))
                plt.plot(t, sdata[:,ip])
    
    plt.legend((active_pfts), loc=0, fontsize=12)
    plt.xlabel(t_unit, fontsize=12, fontweight='bold')
    plt.ylabel(varname, fontsize=12, fontweight='bold')

    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)


# Plot Grided data
def GridedVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)

    plt.subplots_adjust(left=0.065, bottom=None, right=0.99, top=0.98,
                wspace=0.33, hspace=None)
    
    if('TOTSOMC' in varname):         
        varname=varname+' (kgC/m2)'
        sdata = sdata/1000.0

    if(varname == 'SNOW' or varname == 'RAIN'):         
        varname=varname+' (mm/d)'
        sdata = sdata*86400.0
    if(varname in ['SNOW_DEPTH','SNOWDP']):         
        varname=varname+' (m)'

    
    gridtext = []

    # add a zero line for a few variables
    if(varname in ['GPP','NPP','NEE','MR','GR']):
        plt.plot([min(t),max(t)],[0.0,0.0],'k--')
        gridtext.append("0 line ")

    # add a zero-degree-C line for 'TBOT'
    if(varname == 'TBOT'):
        plt.plot([min(t),max(t)],[273.15,273.15],'k--')
        gridtext.append("0oC")

    if(len(sdata.shape)>1):
        for igrd in range(sdata.shape[1]):
            #gridtext.append(("GRID "+str(igrd)))
            plt.plot(t, sdata[:,igrd])
        #gridtext.append(["NAMC","DSLT","AS","WBT","TT-WBT","TT"])
        gridtext.append("Grid_Alert-CANADA")

    else:
        #gridtext.append(("GRID "+str(0)))
        plt.plot(t, sdata)
        
    
    plt.legend((gridtext), loc=0, fontsize=12)
    
    plt.xlabel(t_unit, fontsize=12, fontweight='bold')    
    plt.ylabel(varname, fontsize=12, fontweight='bold')
    
    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    plot_title = ''
    #plot_title = 'ICB-Highlat_pt406x22' #plotlabel
    plt.text(lx, ly, plot_title, transform=ax.transAxes,fontsize=14, fontweight='bold')

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--clmout_dir", dest="clm_odir", default="./", \
                  help="clm output directory (default = ./, i.e., under current directory)")
parser.add_option("--clmout_timestep", dest="clmout_ts", default="monthly", \
                  help="clm output variable timestep (default = 'monthly', other option daily)")
parser.add_option("--clmfile_head", dest="ncfileheader", default="", \
                  help = "clm output file name header, usually the portion before *.clm2.h[0-5].*.nc")
parser.add_option("--clmfile_fincl", dest="ncfincl", default="h0", \
                  help = "clm output file numbering, h[0-5]")
parser.add_option("--varname_help", dest="varnames_print", default=False, \
                  help = "print out VARIABLES available ", action="store_true")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4) to be reading/plotting, separated by comma ")
parser.add_option("--adspinup", dest="adspinup", action="store_true", default=False, \
                  help="whether results of an ad_spinup run (default = False)")
parser.add_option("--year0", dest="yr0", default="1", \
                  help="clm simulation starting year (default = 1, this is for spinup; for transient it should be 1850; " \
                   " and can be user-defined)")
parser.add_option("--startyr", dest="startyr", default="1", \
                  help="clm output starting year to plot (default = 1, i.e. first year of simulation" \
                   " and can be user-defined)")
parser.add_option("--endyr", dest="endyr", default="", \
                  help="clm output ending year to plot (default = none, i.e. end of simulation)")
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
parser.add_option("--COLindex", dest="cindex", default=-999, \
                  help = " COLUMN index to be reading/plotting, default -999 for all, with indexing from 0 ")
parser.add_option("--seasonally", dest="seasonally", default=False, \
                  help = "averaged over years to get seasonal", action="store_true")
parser.add_option("--annually", dest="annually", default=False, \
                  help = "averaged over seasons to get annual", action="store_true")

(options, args) = parser.parse_args()

if(options.pindex==-999):
    pft_index = [options.pindex]
else:
    pft_index=[]
    oppfts = options.pindex.split(':')
    for ip in oppfts:    
        pft_index.append(int(ip))

if(options.cindex==-999):
    col_index = [options.cindex]
else:
    col_index=[]
    opcols = options.cindex.split(':')
    for icol in opcols:    
        col_index.append(int(icol))

if(options.zindex==-999):
    layer_index = [options.zindex]
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
if(options.yr0 != 1 and options.startyr != 1):
    startdays = (int(options.startyr)-int(options.yr0))*365.0
else:
    startdays = (int(options.startyr)-1)*365.0
if(options.clmout_ts=='daily'): 
    startdays=startdays+1.0
elif(options.clmout_ts=='monthly'): 
    startdays=startdays-1.0
    
enddays = -9999
if(options.endyr !=""): enddays = int(options.endyr)*365.0

tunit = str.capitalize(options.t_unit)
if(options.annually): tunit = 'YEAR'
if tunit.startswith("H"): 
    day_scaling = 24.0
elif tunit.startswith("Y"): 
    day_scaling = 1.0/365.0
else:
    day_scaling = 1.0
tunit0 = float(options.yr0)*365.0*day_scaling #the simulation year 0 timing (in day)


ix=int(options.xindex);
iy=int(options.yindex);

# read-in datasets from 1 simulation
nx, ny, nlgrnd, nldcmp, ncolumn, npft, varsdata, varsdims, ttunits = \
    CLM_NcRead_1simulation(options.clm_odir, \
                           options.ncfileheader, \
                           options.ncfincl, \
                           options.varnames_print, \
                           varnames, \
                           startdays, enddays, \
                           options.adspinup)


if (options.varnames_print): sys.exit('Variable Names Printing DONE!')
#--------------------------------------------------------------------------------------
# 
vars_list = list(varsdata.keys())    # var_list is in format of 'h*_varname', in which 'varname' is the real variable names

# plotting
nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
if(nvars==2): ncol = 2
if(nvars==3): ncol = 3
if(nvars==4): nrow=2; ncol=2

# dimension max.
if2dgrid = True
if(len(varsdata['topo'].shape)==1): if2dgrid = False

ivar = 0
for var in varnames:
    if(sys.version_info[0]==3): 
        print (ivar,var)
    #elif(sys.version_info[0]==2):
    #    print ivar,var
        
    # names for variable and its time-axis
    for hv in vars_list:
        if re.search(var, hv): 
            var_h = hv
            hinc  = hv.replace(var,'')
            var_t = '%stime' %hinc
            break
        
    vdata = varsdata[var_h]
    vdims = varsdims[var_h]
    
    tt = varsdata[var_t]   # time unit: days
    
    #time dimension
    dim_tt = -999
    if('time' in vdims): 
        dim_tt = vdims.index('time')
        ivar = ivar + 1 # counting vars, with time-series only
    else:
        continue
    
    # processing original data, if needed
    t, gdata, sdata, zdim_indx, pdim_indx = \
        CLMvar_1Dtseries(tt, vdims, vdata, nx, ny, ix, iy, if2dgrid, \
                    annually=options.annually, \
                    seasonally=options.seasonally)

    #plotting
    vname = varnames[varnames.index(var)]
    t = np.asarray(t)*day_scaling + tunit0
    if(zdim_indx<0 and pdim_indx<0):
        GridedVarPlotting(plt, nrow, ncol, ivar, t, tunit, gdata, \
                        vname, '(a) All Grids')
    elif(zdim_indx>0):
        if(len(layer_index)==1 and layer_index[0]<0): layer_index = range(0,nl)
        SoilLayeredVarPlotting(plt, nrow, ncol, ivar, t, tunit, layer_index, sdata, \
                        vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')
    
    elif(pdim_indx>0):
        if(len(pft_index)==1 and pft_index[0]<0): pft_index = range(0,npft)
        PFTVarPlotting(plt, nrow, ncol, ivar, t, tunit, pwt1cell, pft1vidx, pft_index, sdata, \
                        vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')


# printing plot in PDF
ofname = 'Figure_CLM.pdf'
fig = plt.gcf()
fig.set_size_inches(11.5, 8.5)
plt.savefig(ofname)
plt.show()

plt.close('all')


