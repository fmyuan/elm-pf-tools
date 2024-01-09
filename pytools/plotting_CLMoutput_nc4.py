#!/usr/bin/env python

import sys
import glob
import re
import math
from cmath import nan
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import xlim, ylim
from matplotlib.ticker import FormatStrFormatter

from Modules_CLMoutput_nc4 import CLM_NcRead_1simulation
from Modules_CLMoutput_nc4 import CLMvar_1Dtseries


# ---------------------------------------------------------------
# Plot soil data with layers for ONE specific grid
def SoilLayeredVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, layer_index, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)
    if(varname == 'SOILPSI' or varname=='SOILPSI2' or \
       varname == 'SMP' or varname=='SMP2'):
        sdata = -sdata
        ax.set_yscale("log")
        # the following is for plot's purpose - can be commented out
        sdata[sdata<100] = nan
        plt.plot([min(t),max(t)],[100,100],'k--')
        plt.plot([min(t),max(t)],[1.0e8,1.0e8],'k--')
    
    layertext = []
    
    # add a zero-degree-C line for 'TSOI'
    #if(varname == 'TSOI'):
    #    plt.plot([min(t),max(t)],[273.15,273.15],'k--')
    #    layertext.append("0oC ")

    for il in layer_index:
        if il in [0,2,4,6,8,9]:
            layertext.append(("L"+str(il+1)))
            plt.plot(t, sdata[:,il])
        if il>=9: break
        #if (varname!='TSOI') and il>=9: break  # non-soil temperature only need upper 10 real soil layers 
        
    # limits of axis
    #if 'TSOI' in varname: plt.ylim([253,293])
    #if 'SOILICE' in varname: plt.ylim([0,450])

    if 'TSOI' in varname: plt.legend((layertext), loc=0, ncol=2, fontsize=10)
    
    # appending unit
    if varname=='TSOI': varname='Soil-layer Temperature (K)'
    if varname=='SOILICE': varname='Soil-layer Ice (kg/m2)'
 
 
    plt.xlabel(t_unit, fontsize=10, fontweight='bold')
    plt.ylabel(varname, fontsize=10, fontweight='bold')

    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 10)

    lx = 0.10
    ly = 0.90
    plot_title = '' #plotlabel
    plt.text(lx, ly, plot_title, transform=ax.transAxes, fontsize=18, fontweight='bold')

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
    
    if isubplot==1: plt.legend((active_pfts), loc=0, fontsize=10)
    plt.xlabel(t_unit,  fontsize=10, fontweight='bold')
    plt.ylabel(varname, fontsize=10, fontweight='bold')

    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 10)

    lx = 0.10
    ly = 0.90
    plot_title = '' #plotlabel
    plt.text(lx, ly, plot_title, transform=ax.transAxes, fontsize=18, fontweight='bold')


# Plot Grided data
def GridedVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, sdata, sdata_std=None, varname='', plotlabel='', plottype='default'):

    ax=plt.subplot(nrow, ncol, isubplot)

    plt.subplots_adjust(left=0.065, bottom=None, right=0.99, top=0.98,
                wspace=0.33, hspace=None)

    if(varname in ['GPP','NEE','HR']):         
        sdata = sdata*1000000.0 # unit change to ug/m2/s
        if sdata_std is not None: 
            sdata_std = sdata_std*1000000.0
    
    if('TOTSOMC' in varname):         
        varname=varname+' (kgC/m2)'
        sdata = sdata/1000.0

    if(varname == 'SNOW' or varname == 'RAIN'):
        varname=varname+' (mm/d)'
        sdata = sdata*86400.0
    if(varname in ['SNOW_DEPTH','SNOWDP']):
        varname='Snow Depth'+' (m)'

    if(varname in ['ALT']):
        varname='Active Layer Depth (m)'
        #ax.set_yscale("log")
        #sdata[np.where(sdata<=0.0)] = nan
    
    gridtext = []

    
    # add a zero line for a few variables
    if(varname in ['NEE']):
        plt.plot([min(t),max(t)],[0.0,0.0],'k--')
        gridtext.append("0 line ")

    # add a zero-degree-C line for 'TBOT'
    #if(varname == 'TBOT' or varname == 'TSOI'):
        #plt.plot([min(t),max(t)],[273.15,273.15],'k--')
        #plt.plot([1,12],[273.15,273.15],'k--')
        #gridtext.append("0 oC")

    if(len(sdata.shape)>1):
        if sdata_std is not None:
            for i in range(sdata.shape[1]):
                #igrd = i
                igrd = sdata.shape[1]-i-1 # this will reversely plot data point
                gridtext.append(("Site "+str(igrd+1)))
                #gridtext.append(("Area for Site "+str(igrd+4)))
                if (t_unit=='MONTH'):
                    #slightly shift t-axis so that monthly data point gapped
                    plt.errorbar((t+igrd)*12.0/365.0, sdata[:,igrd], yerr=sdata_std[:,igrd], linestyle='-', marker='.', capsize=2)
                    plt.xticks(np.arange(1, 13, 1))
                    #plt.ylim([0.0,3.0])
                else:
                    plt.errorbar(t, sdata[:,igrd], yerr=sdata_std[:,igrd], linestyle='-', marker='.', capsize=2)
        elif plottype=='errorbar':
            plt.errorbar(t, np.nanmean(sdata,1), yerr=np.nanstd(sdata,1), linestyle='-', marker='.', capsize=2)
        else:
            for i in range(sdata.shape[1]):
                #igrd = i
                igrd = sdata.shape[1]-i-1 # this will reversely plot data point
                gridtext.append(("Site "+str(igrd+1)))
                #gridtext.append(("Area for Site "+str(igrd+4)))
                if(varname == 'SNOW' or varname == 'RAIN'):
                    plt.bar(t, sdata[:,igrd])
                else:
                    plt.plot(t, sdata[:,igrd], label=("Site "+str(igrd+1)))
                
                #trending line added
                if (varname=='FSNO' or 'Snow Depth' in varname):
                    tindx=np.where((t>=1980) & (t<=2015))[0]
                    t1=t[tindx]
                    y1=sdata[tindx,igrd]
                    y1_trend = np.polyfit(t1, y1, 1)
                    funcy1 = np.poly1d(y1_trend)
                    plt.plot(t1,funcy1(t1), '--', label=None)

                    #tindx=np.where((t>=2002) & (t<=2014))[0]
                    #t2=t[tindx]
                    #y2=sdata[tindx,igrd]
                    #y2_trend = np.polyfit(t2, y2, 1)
                    #funcy2 = np.poly1d(y2_trend)
                    #plt.plot(t2,funcy2(t2), '--', label=None)

                
            #gridtext.append(["NAMC","DSLT","AS","WBT","TT-WBT","TT"])

    else:
        #gridtext.append(("GRID "+str(0)))
        if sdata_std is not None:
            plt.errorbar(t, sdata, yerr=sdata_std, linestyle='-', marker='.', capsize=2)
        elif(varname == 'SNOW' or varname == 'RAIN'):
            plt.bar(t, sdata)
        else:
            plt.plot(t, sdata)
        
    # user-defined x/y limits
    #if ('YEAR' in t_unit): plt.xlim(1980,2019)
    
    #if ('LAI' in varname): plt.ylim([0.1,0.9])
    #if ('Active Layer Depth' in varname): plt.ylim([-0.1,1.5])
    #if ('SNOW' in varname): plt.ylim([0.0,20.0])
    #if ('Snow Depth' in varname): plt.ylim([0.0,1.0])
    #if ('FSNO' in varname): plt.ylim([0.3,0.8])

     # mannually edit x/y axis labels
    if varname=='TBOT': varname='Air Temperature (K)'
    if varname=='TSOI': varname='Near-surface Soil Temperature (K)'
    if varname=='TLAI': varname='Total LAI (m2/m2)'
    if varname=='GPP': varname='GPP (10^-6 gC/m2/s)'
    if varname=='NEE': varname='NEE (10^-6 gC/m2/s)'
    if varname=='HR': varname='Hetero Respiration (10^-6 gC/m2/s)'
    if varname=='total_SOILICE': varname='total Soil Ice (kg/m2)'
    if varname=='FSNO': varname='Snow cover fraction (-)'
    #varname = ' surface water head (mm)'
    #varname = 'Snow water equivalent (mm)'

    # x/y ticklabel properties
    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'both', labelsize = 16)
    ax_user.set_xticklabels(ax_user.get_xticks(), weight='bold')
    ax_user.xaxis.set_major_formatter(FormatStrFormatter('%i'))
    ax_user.set_yticklabels(ax_user.get_yticks(), weight='bold')
    ax_user.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if ('Active Layer ' in varname): ax_user.invert_yaxis()

    lx = 0.30
    ly = 0.95
    plot_title = ''
    #if plotlabel!='' and isubplot==2: plot_title = 'Site 7'#'Seward Peninsula, AK' #plotlabel, if any
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.text(lx, ly, plot_title, transform=ax.transAxes, fontsize=18, fontweight='bold')
    if isubplot==1: plt.legend(fontsize=14)
    plt.xlabel(t_unit, fontsize=16, fontweight='bold')
    plt.ylabel(varname, fontsize=16, fontweight='bold')

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
    startdays=startdays
elif(options.clmout_ts=='monthly'): 
    startdays=startdays+1.0
elif(options.clmout_ts=='hhourly'): 
    startdays=startdays+1.0/24.0
    
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
iz=int(options.zindex);
ip=int(options.pindex);

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
#if(nvars==2): ncol = 2
if(nvars==2): ncol = 2
if(nvars==3): ncol = 3
if(nvars==4): nrow=2; ncol=2

# dimension max.
if2dgrid = True
if(nx==1 or ny==1): if2dgrid = False

# time-invariant variables
constdata = {}
for v in vars_list:
    if 'time' not in varsdims[v]: 
        constdata[v]=varsdata[v]

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
            hinc  = options.ncfincl+'_' #hv.replace(var,'')
            var_t = '%stime' %hinc
            break
        
    vdata = varsdata[var_h]
    if (var=='ALT'): vdata[np.where(vdata>=3.801)]=3.801 # the 10-soil-layer bottom depth is about 3.80m, not the whole 15 layer bottom (~42m)
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
    if (options.annually or options.seasonally):
        t, gdata, gdata_std, sdata, sdata_std, zdim_indx, pdim_indx, npft, ncolumn = \
            CLMvar_1Dtseries(tt, vdims, vdata, constdata, nx, ny, ix, iy, if2dgrid, \
                    izp=max(iz,ip), \
                    annually=options.annually, \
                    seasonally=options.seasonally)
    else:
        t, gdata, sdata, zdim_indx, pdim_indx, npft, ncolumn  = \
            CLMvar_1Dtseries(tt, vdims, vdata, constdata, nx, ny, ix, iy, if2dgrid, \
                    izp=max(iz,ip), \
                    annually=options.annually, \
                    seasonally=options.seasonally)

    #plotting
    vname = varnames[varnames.index(var)]
    t = np.asarray(t)*day_scaling + tunit0

    #if manually offset time
    t0 = 0#(2014-0)*365 #0
    #if(tunit=='DOY'):t0=(2014-0)*365
    #if(tunit=='YEAR'):t0=0
    t = (t - t0)
    tunit = tunit.upper()

    if((zdim_indx<0 and pdim_indx<0) or max(iz,ip)>=0):
        if (options.annually or options.seasonally):
            #GridedVarPlotting(plt, nrow, ncol, ivar, t, tunit, gdata, sdata_std=gdata_std, \
            #            varname=vname, plotlabel='1kmx1km Resolution', plottype='errorbar')
            GridedVarPlotting(plt, nrow, ncol, ivar, t, tunit, gdata, \
                        #varname=vname, plotlabel='0.5degx0.5deg Resolution', plottype='errorbar')
                        varname=vname, plotlabel='1kmx1km Resolution')
        else:
            GridedVarPlotting(plt, nrow, ncol, ivar, t, tunit, gdata, \
                        #varname=vname, plotlabel='Site '+str(iy+1))
                        varname=vname, plotlabel='1kmx1km Resolution')
                        #varname=vname, plotlabel='0.5degx0.5deg Resolution')
            #GridedVarPlotting(plt, nrow, ncol, ivar, t, tunit, gdata, \
            #                varname=vname, plotlable-'0.5x0.5 degree Resolution')#, \
            #                #'All Grids', plottype='errorbar')
        
        
        
    elif(zdim_indx>0):
        if(len(layer_index)==1 and layer_index[0]<0): layer_index = range(0,nlgrnd)
        SoilLayeredVarPlotting(plt, nrow, ncol, ivar, t, tunit, layer_index, sdata, \
                        vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')
    
    elif(pdim_indx>0):
        pwt1cell = constdata['pfts1d_wtgcell']
        pft1vidx = constdata['pfts1d_itype_veg']
        if(len(pft_index)==1 and pft_index[0]<0): pft_index = range(0,npft)
        PFTVarPlotting(plt, nrow, ncol, ivar, t, tunit, pwt1cell, pft1vidx, pft_index, sdata, \
                        vname, '(a) Grid ( '+str(ix)+', '+str(iy)+')')


# printing plot in PDF
ofname = 'Figure_CLM.pdf'
fig = plt.gcf()
fig.set_size_inches(11.5, 6.0)
plt.savefig(ofname)
plt.show()

plt.close('all')


