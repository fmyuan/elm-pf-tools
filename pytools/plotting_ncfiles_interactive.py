#!/usr/bin/env python

import sys
import os
import glob
import re
import math
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset
import mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D

# ---------------------------------------------------------------
# Plot grid surface data all grids
def One2OnePlotting(plt, nrow, ncol, isubplot, data1, data2, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)

    plt.plot(data1, data2)

    plt.legend('1:1', loc=0, fontsize=12)
    plt.xlabel('data 1', fontsize=12, fontweight='bold')
    plt.ylabel('data 2', fontsize=12, fontweight='bold')

    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

# Plot grid surface data all grids
def Grid2dSurfVarPlotting(plt, nrow, ncol, isubplot, lat1d, lon1d, gdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot, projection='3d')

    X, Y = np.meshgrid(lon1d, lat1d)
    ax.plot_surface(X, Y, gdata, rstride=1, cstride=1,
                       linewidth=0, antialiased=False)

    #plt.legend((layertext), loc=0, fontsize=12)
    #plt.xlabel('LONGITUDE', fontsize=12, fontweight='bold')
    #plt.ylabel('LATITUDE', fontsize=12, fontweight='bold')

    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    #plt.text(lx, ly, plotlabel, transform=ax.transAxes)

# ---------------------------------------------------------------
# Plot time-series soil data with layers for ONE specific grid
def TimeSoilLayeredVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, layer_index, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    if(varname == 'SOILPSI'):
        sdata = -sdata
        ax.set_yscale("log", nonposy='clip')
    
    layertext = []
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
def TimePFTVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, pftwt, pftvidx, pft_index, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    
    active_pfts = []
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
def TimeGridedVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)

    plt.subplots_adjust(left=0.065, bottom=None, right=0.99, top=0.98,
                wspace=0.33, hspace=None)
    
    if('TOTSOMC' in varname):         
        varname=varname+' (kgC/m2)'
        sdata = sdata/1000.0

    if(varname == 'SNOW'):         
        varname=varname+' (mm/d)'
        sdata = sdata*86400.0
    if(varname == 'SNOW_DEPTH'):         
        varname=varname+' (m)'
    
    gridtext = []
    if(len(sdata.shape)>1):
        for igrd in range(sdata.shape[1]):
            gridtext.append(("GRID "+str(igrd)))
            plt.plot(t, sdata[:,igrd])
    else:
        gridtext.append(("GRID "+str(0)))
        plt.plot(t, sdata)
        
    gridtext = ["Grid_Uelen-RUSSIA"]
    #gridtext = ["NAMC","DSLT","AS","WBT","TT-WBT","TT"]
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

parser.add_option("--work_dir", dest="workdir", default="./", \
                  help="work directory (default = ./, i.e., under current directory)")
parser.add_option("--ncfile_header", dest="ncfileheader", default="", \
                  help = "nc file name header, usually the portion before .nc")

(options, args) = parser.parse_args()


#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)


options.ncfileheader = 'MODIS_ELM_HighLat'
if (options.ncfileheader == ''):
    print('MUST  have file header by " --ncfile_header=??? "' //
     ' , which usually is the portion before "*.nc"! ')
    sys.exit()
else:
    print('nc file header: '+ options.ncfileheader)
    ncfile = options.workdir+'/'+options.ncfileheader+'.nc'


# read-in a file
try:
    f = Dataset(ncfile,'r')
except:
   print('nc file NOT exists: ' + options.ncfile)
    
#--------------------------------------------------------------------------------------

vname='timeint_of_lai'
nvars = 1

vars = f.groups['MeanState'].variables
varnames = vars.keys()
gdata = np.asarray(vars[vname]
lat1d = np.asarray(vars['lat']
lon1d = np.asarray(vars['lon']
# 
# plotting
nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
if(nvars==2): ncol = 2
if(nvars==3): ncol = 3
if(nvars==4): nrow=2; ncol=2

#sub-plot no.
ivar = 1

#Grid2dSurfVarPlotting(plt, nrow, ncol, ivar, lat1d, lon1d, gdata, \
#                    vname, '')

ncfile0 = 'domain.lnd.360x720_cruncep_vji.170516_N60.nc'
f0 = Dataset(ncfile0,'r')
vars0 = f0.variables
varnames0 = vars0.keys()
gdata0 = np.asarray(vars0['frac'])
lon1d0 = np.asarray(vars0['xc'][0,:])
lat1d0 = np.asarray(vars0['yc'][:,0])


One2OnePlotting(plt, nrow, ncol, ivar, gdata0, gdata, plotlabel)
    
    
    #GridedVarPlotting(plt, nrow, ncol, ivar, t, tunit, gdata, \
    #                    vname, '(a) All Grids')

    #SoilLayeredVarPlotting(plt, nrow, ncol, ivar, t, tunit, layer_index, sdata, \
    #                    vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')
    
    #PFTVarPlotting(plt, nrow, ncol, ivar, t, tunit, pwt1cell, pft1vidx, pft_index, sdata, \
    #                    vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')


# printing plot in PDF
ofname = 'Figure_CLM.pdf'
fig = plt.gcf()
fig.set_size_inches(11.5, 8.5)
plt.savefig(ofname)
plt.show()

plt.close('all')


