#!/usr/bin/env python

import sys
import glob
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ---------------------------------------------------------------
# plotting 1 graph with at most 4 sub-plots 

def SubPlotting(tt, time_unit, varnames, varunits, vardatas, figno=None):
    nvars = len(varnames)
    
    nrow = 1   # sub-plot vertically arranged number (row no.)
    ncol = 1   # sub-plot horizontally arranged number (column no.)
    if(nvars>=2):
        nrow = 2
    if(nvars>=3):
        ncol = 2


    fig = plt.figure(figsize=(11.5,8.5))
     
    # plot 1, or subplot 1 if more than 1 subplot
    if (nvars >= 2):
        sdata = vardatas[0][:]
    else:
        sdata = vardatas
       
    ax0=plt.subplot(nrow, ncol, 1)
    if (not figno is None): plt.suptitle('FIGURE '+figno)
    plt.plot(tt, sdata)
    plt.xlabel('Time ('+time_unit+')')
    plt.ylabel(varnames[0]+varunits[0])

    lx = 0.05
    ly = 0.95
    if(nvars>=2):
        plt.text(lx, ly, '(a) ', transform=ax0.transAxes)
    else:
        plt.text(lx, ly, '', transform=ax0.transAxes)
        
    # subplot 2
    if (nvars >= 2):
        sdata = vardatas[1][:]
        
        ax1=plt.subplot(nrow, ncol, 2)
        plt.plot(tt, sdata)
        plt.xlabel('Time ('+time_unit+')')
        plt.ylabel(varnames[1]+varunits[1])
        lx = 0.05
        ly = 0.95
        plt.text(lx, ly, '(b) ', transform=ax1.transAxes)

    # subplot 3
    if (nvars >= 3):
        sdata = vardatas[2][:]

        ax2=plt.subplot(nrow, ncol, 3)
        plt.plot(tt, sdata)
        plt.xlabel('Time ('+time_unit+')')
        plt.ylabel(varnames[2]+varunits[2])
        lx = 0.05
        ly = 0.95
        plt.text(lx, ly, '(c) ', transform=ax2.transAxes)

    # sub-plot 4
    if (nvars >= 4):
        sdata = vardatas[3][:]

        ax3=plt.subplot(nrow, ncol, 4)
        plt.plot(tt, sdata)
        plt.xlabel('Time ('+time_unit+')')
        plt.ylabel(varnames[3]+varunits[3])
        lx = 0.05
        ly = 0.95
        plt.text(lx, ly, '(d) ', transform=ax3.transAxes)

    #
    ofname = 'Figure_XX.pdf'
    if (not figno is None): ofname = 'Figure_-'+figno+'.pdf'
    plt.savefig(ofname)
    plt.show()

    plt.close('all')

# ---------------------------------------------------------------
# plotting 1 plot with at most 4 varied-length vars 

def SinglePlotting(tt, time_unit, varnames, varunits, vardatas, figno=None):
    nvars = len(varnames)
    
    nrow = 1   # sub-plot vertically arranged number (row no.)
    ncol = 1   # sub-plot horizontally arranged number (column no.)

    fig = plt.figure(figsize=(10.5,6.5))
    ax0=plt.subplot(nrow, ncol, 1)
    #
    # data must be in dictionary, due to length may not be same
    for var in varnames:
        t = tt[var]
        sdata = vardatas[var]
        plt.plot(t, sdata)
    
    plt.xlabel('Time ('+time_unit+')')
    plt.ylabel(varunits)
    plt.legend((varnames), loc=2, fontsize=12)

    lx = 0.05
    ly = 0.95
    plt.text(lx, ly, '', transform=ax0.transAxes)
        

    #
    ofname = 'Figure_XX.pdf'
    if (not figno is None): ofname = 'Figure_-'+figno+'.pdf'
    plt.savefig(ofname)
    plt.show()

    plt.close('all')


# ---------------------------------------------------------------
# Plot grid surface data all grids
def One2OnePlotting(plt, nrow, ncol, isubplot, data1, data2, datalabel, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    plt.subplots_adjust(left=0.065, bottom=None, right=0.98, top=0.97,
                wspace=0.33, hspace=None)

    if len(data1.shape)>1:
        data1=np.reshape(data1,(data1.size))
        data2=np.reshape(data2,(data2.size))
    l11 = [np.nanmin(data1)-0.1*abs(np.nanmin(data1)), \
           np.nanmin(data2)-0.1*abs(np.nanmin(data2)), \
           np.nanmax(data1)+0.1*abs(np.nanmax(data1)), \
           np.nanmax(data2)+0.1*abs(np.nanmax(data2))]
    
    plt.plot(data1, data2,'x')
    plt.plot(l11,l11,'--')
    
    # ME/RMSE
    me=np.nanmean(data1-data2)
    rmse=np.nanmean([x*x for x in (data1-data2)])
    rmse=np.sqrt(rmse)
    plt.text(0.70, 0.075, 'ME='+str(round(me,1)), fontweight='bold',transform=ax.transAxes)
    plt.text(0.70, 0.030, 'RSME='+str(round(rmse,1)), fontweight='bold',transform=ax.transAxes)

    plt.legend([datalabel,'1:1'], loc=2, fontsize=12)
    plt.xlabel(datalabel+'_1', fontsize=12, fontweight='bold')
    plt.ylabel(datalabel+'_2', fontsize=12, fontweight='bold')

    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    plt.text(lx, ly, plotlabel, transform=ax.transAxes)

# ---------------------------------------------------------------
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

# ---------------------------------------------------------------
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

# ---------------------------------------------------------------

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
        
    gridtext = varname#["ELM Simulation","NIC-IMS Observation"]
    plt.legend((gridtext), loc=2, fontsize=12)
    

    plt.xlabel(t_unit, fontsize=12, fontweight='bold')    
    plt.ylabel(varname, fontsize=12, fontweight='bold')
    
    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    #plot_title = ''
    plot_title = plotlabel#'ICB-Highlat_pt406x22' #plotlabel
    plt.text(lx, ly, plot_title, transform=ax.transAxes,fontsize=14, fontweight='bold')


