#!/usr/bin/env python

from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

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
    plt.ylabel('NIC-IMS Derived', fontsize=12, fontweight='bold')

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
        
    gridtext = ['']
    plt.legend((gridtext), loc=2, fontsize=12)
    

    plt.xlabel(t_unit, fontsize=12, fontweight='bold')    
    plt.ylabel(varname, fontsize=12, fontweight='bold')
    #plt.xlim(2003,2022)
    #plt.ylim(100,165)
    
    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    plot_title = ''
    #plot_title = 'ICB-Highlat_pt406x22' #plotlabel
    plt.text(lx, ly, plot_title, transform=ax.transAxes,fontsize=14, fontweight='bold')

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help="work directory (default = ./, i.e., under current directory)")

(options, args) = parser.parse_args()


vname = 'FLDS'

#ncfile = options.workdir+'./CRUJRAV2.3.c2023_daymet4_FSDS_1980-2021_z03.nc'
#ncfile = options.workdir+'./Daymet_ERA5.1km_'+vname+'_1980-2023_z01.nc'
#ncfile = options.workdir+'./ERA5_TBOT_1950-2024_z01.nc'
ncfile = options.workdir+'./ATS-subdaily_'+vname+'_z01.nc'


#--------------------------------------------------------------------------------------
# read-in a file
try:
    f = Dataset(ncfile,'r')
except:
    print('nc file NOT exists: ' + ncfile)

gdata = f.variables[vname][0,...]
gdata = gdata.T  # if need to transpose or swap grid/time axis
gdata2 = []
datalabel = ''

t = f.variables['DTIME'][...]
tunit = ''
# 
#--------------------------------------------------------------------------------------

# plotting
nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
isub = 1
#++++++++++++++++++++++++++++++

# 1:1 plotting

ONE2ONE_PLOTTING = False
if(ONE2ONE_PLOTTING):
    plotlabel = ''
    data1 = np.nanmean(gdata,axis=0) # timely-averaged
    data2 = np.nanmean(gdata2,axis=0) # timely-averaged
    One2OnePlotting(plt, nrow, ncol, isub, data1, data2, datalabel, plotlabel)
    
    
#++++++++++++++++++++++++++++++

# time-series plotting
T_PLOTTING = True
if ONE2ONE_PLOTTING: isub = 2
if (T_PLOTTING):
    TimeGridedVarPlotting(plt, nrow, ncol, isub, t, tunit, gdata, \
                    vname, 'Mean of all-grids')


# printing plot in PDF
ofname = 'Figure_CLM_obs.pdf'
fig = plt.gcf()
fig.set_size_inches(11.5, 8.5)
plt.savefig(ofname)
plt.show()

plt.close('all')


