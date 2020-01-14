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
        
    gridtext = ["Climate-Grid"] #["NAMC","DSLT","AS","WBT","TT-WBT","TT"]
    plt.legend((gridtext), loc=0, fontsize=12)
    
    plt.xlabel(t_unit, fontsize=12, fontweight='bold')    
    plt.ylabel(varname, fontsize=12, fontweight='bold')
    
    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    plot_title = ''
    if(varname=='TLAI'): plot_title = 'ICB-Highlat_pt406x22' #plotlabel
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
nx, ny, nlgrnd, nldcmp, ncolumn, npft, varsdata, varsdims = \
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

# dimension max.
if2dgrid = True
if(len(varsdata['topo'].shape)==1): if2dgrid = False
nxy = nx*ny
nl = max(nlgrnd, nldcmp)

# plotting
nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
if(nvars==2): ncol = 2
if(nvars==3): ncol = 3
if(nvars==4): nrow=2; ncol=2

ivar = 0
for var in varnames:
    print ivar,var
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

    # pft dim, if existed
    pdim_indx = -999
    if('pft' in vdims):
        pdim_indx = vdims.index('pft')
        npft = vdata.shape[pdim_indx] # this actually is nxy*npft
        pwt1cell = varsdata['pfts1d_wtgcell']
        pft1vidx = varsdata['pfts1d_itype_veg']
        pft1active = varsdata['pfts1d_active']

        if(nxy>1):
            npft = npft/nxy
            if(if2dgrid): 
                vdata = vdata.reshape(-1,ny,nx,npft)
                pwt1cell = pwt1cell.reshape(ny,nx,npft)
                pft1vidx = pft1vidx.reshape(ny,nx,npft)
                pft1active = pft1active.reshape(ny,nx,npft)
                pdim_indx = pdim_indx + 2 # 
            else:
                vdata = vdata.reshape(-1, nxy,npft)
                pwt1cell = pwt1cell.reshape(nxy,npft)
                pft1vidx = pft1vidx.reshape(nxy,npft)
                pft1active = pft1active.reshape(nxy,npft)
                pdim_indx = pdim_indx + 1 # 
            
            if(len(pft_index)==1 and pft_index[0]<0): 
                # when NOT output specific PFT(s) for all grid, sum all pfts
                if(ix<0 and iy<0):
                    vdata = np.sum(vdata*pwt1cell,axis=pdim_indx)
                    pdim_indx = -999 # because weighted-sum, pft-dimension is removed (no more 3-D data)

            if(ix>=0 and iy>=0):
                # when output for a specific grid, need to extract that grid's pft info
                # (BUT don't do so for 'vdata', which will do so when passing to 'sdata' 
                if(if2dgrid): 
                    pwt1cell = pwt1cell[iy,ix,:]
                    pft1vidx = pft1vidx[iy,ix,:]
                    pft1active = pft1active[iy,ix,:]
                else:
                    pwt1cell = pwt1cell[max(iy,ix),:]
                    pft1vidx = pft1vidx[max(iy,ix),:]
                    pft1active = pft1active[max(iy,ix),:]


    # column dim, if existed
    cdim_indx = -999
    if('column' in vdims):
        cdim_indx = vdims.index('column')  
        ncolumn = vdata.shape[cdim_indx] # this actually is nxy*ncolumn
        colwt1cell = varsdata['cols1d_wtgcell']
        col1active = varsdata['cols1d_active']

        if(nxy>1):
            ncolumn = ncolumn/nxy
            if(if2dgrid): 
                vdata = vdata.reshape(-1,ny,nx,ncolumn)
                colwt1cell = colwt1cell.reshape(ny,nx,ncolumn)
                col1active = col1active.reshape(ny,nx,ncolumn)
                cdim_indx = cdim_indx + 2 # 
            else:
                vdata = vdata.reshape(-1, nxy,ncolumn)
                colwt1cell = colwt1cell.reshape(nxy,ncolumn)
                col1active = col1active.reshape(nxy,ncolumn)
                cdim_indx = cdim_indx + 1 # 
        
        if(len(col_index)==1):
            if(col_index[0]<0): # when NOT output specific COLUMN(s), sum all cols
                vdata = np.sum(vdata*colwt1cell,axis=cdim_indx)
                cdim_indx = -999 # because weighted-sum, column-dimension is removed (no more 3-D data)
        
    # data series
    gdata=[]
    sdata=[]
    if(zdim_indx<0 and pdim_indx<0):# 2-D grid data
        gdata = np.zeros((nt,nx*ny))    #temporary data holder in 2-D (tt, grids)
    elif(zdim_indx>0):        
        sdata = np.zeros((nt,nl*nx*ny)) #temporary data holder in 2-D (tt, layers*grids)
    elif(pdim_indx>0):        
        if(iy>=0 and ix>=0): #[ix,iy] is location of 2-D grids, starting from 0
            sdata = np.zeros((nt,npft))      #temporary data holder in 2-D (tt, pfts*1 grid)
        if(iy>=0):
            sdata = np.zeros((nt,npft*nx))   #temporary data holder in 2-D (tt, pfts*nx grid)
        if(ix>=0):
            sdata = np.zeros((nt,npft*ny))   #temporary data holder in 2-D (tt, pfts*ny grid)
        else:
            sdata = np.zeros((nt,npft*nx*ny)) #temporary data holder in 2-D (tt, pfts*grids)
    
    else:
        exit("Variable to be plotted has 4-D dimension - NOT YET supported!")
    
    for i in range(len(tt)):
                
        if(zdim_indx<0 and pdim_indx<0):# 2-D grid data
            if(ix>=0 or iy>=0):
                if(if2dgrid): 
                    if(ix<0):
                        gdata[i,:] = vdata[it[i]][iy,:]
                    elif(iy<0):
                        gdata[i,:] = vdata[it[i]][:,ix]
                    else:
                        gdata[i,:] = vdata[it[i]][iy,ix]
                else:
                    gdata[i,:] = vdata[it[i]][max(iy,ix)]

            else: # all grids
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
                        
                
            elif(zdim_indx >= 2 or pdim_indx==2): # 3-D soil data, z_dim/p_dim in 2 or 3 (likely x/y dims before z) 
                if(if2dgrid): 
                    sdata[i,:] = vdata[it[i]][:,iy,ix,]
                else:
                    sdata[i,:] = vdata[it[i]][max(iy,ix),]

    if(options.seasonally):
        t=np.asarray(t)/365.0
        dim_yr=int(math.ceil(max(t))-math.floor(min(t)))
        t=(t-np.floor(t))*365.0 # still in days
        t=t.reshape(dim_yr,-1)
        t=np.mean(t,axis=0)
        t=np.where(t==0,365.0,t)  # day 0 shall be day 365, otherwise plotting X axis not good
        
        if(len(gdata)>0):
            shp=np.hstack(([dim_yr,-1],gdata.shape[1:]))
            gdata=gdata.reshape(shp)
            gdata=np.mean(gdata,axis=0)
        elif(len(sdata)>0):
            shp=np.hstack(([dim_yr,-1],sdata.shape[1:]))
            sdata=sdata.reshape(shp)
            sdata=np.mean(sdata,axis=0)
    elif(options.annually):
        t=np.asarray(t)/365.0
        dim_yr=int(math.ceil(max(t))-math.floor(min(t)))
        t=np.floor(t)*365.0 # still in days
        t=t.reshape(dim_yr,-1)
        dim_season = t.shape[1]
        t=np.mean(t,axis=1)
        
        if(len(gdata)>0):
            shp=np.hstack(([-1, dim_season],gdata.shape[1:]))
            gdata=gdata.reshape(shp)
            gdata=np.mean(gdata,axis=1)
        elif(len(sdata)>0):
            shp=np.hstack(([-1, dim_season],sdata.shape[1:]))
            sdata=sdata.reshape(shp)
            sdata=np.mean(sdata,axis=1)

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


