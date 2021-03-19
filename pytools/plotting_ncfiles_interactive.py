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
def One2OnePlotting(plt, nrow, ncol, isubplot, data1, data2, datalabel, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    plt.subplots_adjust(left=0.065, bottom=None, right=0.98, top=0.97,
                wspace=0.33, hspace=None)

    if len(data1.shape)>1:
        data1=np.reshape(data1,(data1.size))
        data2=np.reshape(data2,(data2.size))
    l11 = [np.nanmin(data1)-1,np.nanmin(data2)-1,np.nanmax(data1)+1,np.nanmax(data2)+1]
    
    plt.plot(data1, data2,'x')
    plt.plot(l11,l11,'--')
    
    # ME/RMSE
    me=np.nanmean(data1-data2)
    rmse=np.nanmean([x*x for x in (data1-data2)])
    rmse=np.sqrt(rmse)
    plt.text(0.70, 0.075, 'ME='+str(round(me,1)), fontweight='bold',transform=ax.transAxes)
    plt.text(0.70, 0.030, 'RSME='+str(round(rmse,1)), fontweight='bold',transform=ax.transAxes)

    plt.legend([datalabel,'1:1'], loc=2, fontsize=12)
    plt.xlabel('ELM Simulation', fontsize=12, fontweight='bold')
    plt.ylabel('NIC-IMS Observation', fontsize=12, fontweight='bold')

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
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

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

    plt.subplots_adjust(left=0.065, bottom=None, right=0.98, top=0.97,
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
            plt.plot(t, sdata[:,igrd], 'o-')

        # difference, if 2 time-series dataset only
        if(sdata.shape[1]==2):
            diff=sdata[:,0]-sdata[:,1]
            plt.text(0.70, 0.075, 'Diff_mean='+str(round(np.nanmean(diff),1)), fontweight='bold',transform=ax.transAxes)
            plt.text(0.70, 0.030, 'Diff_Stdev='+str(round(np.nanstd(diff),1)), fontweight='bold',transform=ax.transAxes)

    else:
        gridtext.append(("GRID "+str(0)))
        plt.plot(t, sdata, 'o-')
        
    gridtext = ["ELM Simulation","NIC-IMS Observation"]
    plt.legend((gridtext), loc=2, fontsize=12)
    

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
parser.add_option("--ncfile2_header", dest="ncfile2header", default="", \
                  help = "nc file 2 name header, usually the portion before .nc")
parser.add_option("--varname", dest="varname", default="", \
                  help = "variable name in NC files to be plotting")
parser.add_option("--datalabel", dest="datalabel", default="", \
                  help = "when plot 'varname', labelling it")

(options, args) = parser.parse_args()


#
if (options.workdir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.workdir)


if (options.ncfileheader == ''):
    print('MUST  have file header by " --ncfile_header=??? "'
     ' , which usually is the portion before "*.nc"! ')
    sys.exit()
else:
    print('nc file header: '+ options.ncfileheader)
    ncfile = options.workdir+'/'+options.ncfileheader+'.nc'


#--------------------------------------------------------------------------------------
# read-in a file
try:
    f = Dataset(ncfile,'r')
except:
   print('nc file NOT exists: ' + options.ncfile)
    

vars = f.variables
varnames = vars.keys()

nvars = 2

vname= options.varname #'Days_snowfree'
#vname='Doy_snowmelted'
#vname='Doy_snowcovered'
if options.datalabel == "":
    datalabel = vname
else:
    datalabel = options.datalabel #'Yearly Snow-free Days'
    #datalabel = 'DOY Snow-melted'
    #datalabel = 'DOY Snow-covered'

if ('lat' in varnames or 'lon' in varnames):
    lat1d = np.asarray(vars['lat'])
    lon1d = np.asarray(vars['lon'])

elif('geoy' in varnames or 'geox' in varnames):
    lat1d = np.asarray(vars['geoy'])
    lon1d = np.asarray(vars['geox'])

if (vname+'_elm' in varnames):
    gdata = np.asarray(vars[vname+'_elm'])
    gdata2= np.asarray(vars[vname])
elif (vname+'_ims' in varnames):
    gdata = np.asarray(vars[vname])
    gdata2= np.asarray(vars[vname+'_ims'])

t     = np.asarray(vars['time'])
tunit = vars['time'].units
# 
# observation or dataset 2 in separated nc file, but with same 'varname'
if(options.ncfile2header != ''): # from another nc file(s)
    ncfile2 = options.ncfile2header
    f2 = Dataset(ncfile2,'r')
    vars2 = f2.variables
    varnames2 = vars2.keys()
    gdata2 = np.asarray(vars2[vname])

# IMS data in 1998/2020 not complete
tij=np.where((t==1998) | (t==2020))
gdata2[tij,] = np.nan

# remove 'outlier'
if('snowfree' in vname):
    ij=np.where((gdata<=0) | (gdata>=365))
    gdata[ij] = np.nan
    ij=np.where((gdata2<=0) | (gdata2>=365))
    gdata2[ij]= np.nan
if('snowmelted' in vname):
    ij=np.where((gdata<=0) | (gdata>=365))
    gdata[ij] = np.nan
    ij=np.where((gdata2<=0) | (gdata2>=365))
    gdata2[ij]= np.nan
if('snowcovered' in vname):
    ij=np.where((gdata<=1) | (gdata>=365))
    gdata[ij] = np.nan
    ij=np.where((gdata2<=1) | (gdata2>=365))
    gdata2[ij]= np.nan
#--------------------------------------------------------------------------------------

# plotting
nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
if(nvars==2): ncol = 2
if(nvars==3): ncol = 3
if(nvars==4): nrow=2; ncol=2

#++++++++++++++++++++++++++++++

# 1:1 plotting
#sub-plot no.
isub = 1

ONE2ONE_PLOTTING = True
if(ONE2ONE_PLOTTING):
    plotlabel = ''
    data1 = np.nanmean(gdata,axis=0) # timely-averaged
    data2 = np.nanmean(gdata2,axis=0) # timely-averaged
    One2OnePlotting(plt, nrow, ncol, isub, data1, data2, datalabel, plotlabel)
    
    
#++++++++++++++++++++++++++++++

# time-series plotting
T_PLOTTING = True
isub = 2
if (T_PLOTTING):
    data1 = np.reshape(gdata, (gdata.shape[0],gdata.shape[1]*gdata.shape[2])) # 2-D grid to 1-D
    data1 = np.nanmean(data1, axis=1) # averaged over grids
    data2 = np.reshape(gdata2, (gdata2.shape[0],gdata2.shape[1]*gdata2.shape[2])) # 2-D grid to 1-D
    data2 = np.nanmean(data2, axis=1) # averaged over grids
    sdata = np.swapaxes(np.vstack((data1,data2)),0,1)
    TimeGridedVarPlotting(plt, nrow, ncol, isub, t, tunit, sdata, \
                    vname, 'Mean of all-grids')


# printing plot in PDF
ofname = 'Figure_CLM_obs.pdf'
fig = plt.gcf()
fig.set_size_inches(11.5, 8.5)
plt.savefig(ofname)
plt.show()

plt.close('all')


