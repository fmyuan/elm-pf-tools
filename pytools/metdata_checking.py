#!/usr/bin/env python

import sys
import glob
import re
import math
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages

from Modules_CLMoutput_nc4 import CLM_NcRead_1simulation
from Modules_metdata import clm_metdata_cplbypass_read
from Modules_metdata import clm_metdata_read
from Modules_metdata import singleNCDCReadCsvfile

# ---------------------------------------------------------------
# Plot Grided data
def GridedVarPlotting(plt, nrow, ncol, isubplot, t, t_unit, sdata, varname, plotlabel):

    ax=plt.subplot(nrow, ncol, isubplot)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.99, top=0.99,
                wspace=0.33, hspace=None)

    if(varname in ['SNOW','RAIN']):         
        varname=varname+' (mm/d)'
        sdata = sdata*86400.0
    if(varname in ['SNOW_DEPTH','SNOWDP']):         
        sdata = sdata*100.0
        varname=varname+' (cm)'

    if(varname in ['TBOTd']):         
        varname=varname+' (K)'
    if(varname in ['PRCPd']):         
        varname=varname+' (mm/d)'
    if(varname in ['SNWDd']):         
        varname=varname+' (cm)'

    
    gridtext = []
    # add a zero-degree-C line for 'TSOI'
    if('TBOT' in varname):
        plt.plot([min(t),max(t)],[273.15,273.15],'k--')
        gridtext.append("0oC ")
    
    if(len(sdata.shape)>1):
        for igrd in range(sdata.shape[1]):
            gridtext.append(("GRID "+str(igrd)))
            #gridtext.append("Station")
            plt.plot(t, sdata[:,igrd])
    else:
        gridtext.append(("GRID "+str(0)))
        plt.plot(t, sdata)
        
    #gridtext.append('GRID')#["Climate-Grid"]
    #gridtext.append('['Station']
    plt.legend((gridtext), loc=0, fontsize=12)
    
    plt.xlabel(t_unit, fontsize=12, fontweight='bold')    
    plt.ylabel(varname, fontsize=12, fontweight='bold')
    
    ax_user=plt.gca()
    ax_user.tick_params(axis = 'both', which = 'major', labelsize = 16)

    lx = 0.10
    ly = 0.90
    plot_title = ''
    #plot_title = 'BEO, AK' #plotlabel
    #plot_title = 'VanKarem, RUSSIA' #plotlabel
    plt.text(lx, ly, plot_title, transform=ax.transAxes,fontsize=14, fontweight='bold')

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--e3sm_metdir", dest="met_idir", default="./", \
                  help="e3sminput met directory (default = ./, i.e., under current directory)")
parser.add_option("--e3sm_metdomain", dest="met_domain", default="", \
                  help="e3sminput met domain file (default = a file in current met_idir")
parser.add_option("--e3sm_mettype", dest="met_type", default="CRU", \
                  help="e3sminput met data type (default = 'CRU', options-CRU, cplbypass, GSWP3, GSWP3v1, Site, ELM)")
parser.add_option("--e3sm_metheader", dest="met_header", default="", \
                  help="e3sminput directory (default = '', i.e., standard file name)")
parser.add_option("--clmout_dir", dest="clm_odir", default="./", \
                  help="clm output directory (default = ./, i.e., under current directory)")
parser.add_option("--clmout_header", dest="clmncheader", default="", \
                  help = "clm output file name header, usually the portion before *.clm2.h[0-5].*.nc")
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4) to be reading/plotting, separated by comma ")
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
parser.add_option("--lon", dest="lon", default=-999, \
                  help = " longitude to be reading/plotting, default -999 for all pts or 0 for first pt")
parser.add_option("--lat", dest="lat", default=-999, \
                  help = " latitude to be reading/plotting, default -999 for all pts or 0 for first pt")
parser.add_option("--seasonally", dest="seasonally", default=False, \
                  help = "averaged over years to get seasonal", action="store_true")
parser.add_option("--annually", dest="annually", default=False, \
                  help = "averaged over seasons to get annual", action="store_true")

(options, args) = parser.parse_args()

#
if (options.met_idir == './'):
    print('e3sm met. directory is the current')
else:
    print('e3sm met.  directory: '+ options.met_idir)

if(options.met_type == 'ELM'):
    if (options.clm_odir == './'):
        print('clmdata directory is the current')
    else:
        print('clmdata directory: '+ options.clm_odir)

    if (options.clmncheader == ''):
        print(' NO output clm file header by " --clmncheader=??? "' 
              ' , which usually is the portion before "*.clm2.[h0-h5].*.nc" ')
    else:
        print('clm nc file header: '+ options.clmncheader)


startdays = -9999
if(int(options.yr0) != 1 and int(options.startyr) != 1):
    startdays = (int(options.startyr)-int(options.yr0))*365.0+1.0
else:
    startdays = (int(options.startyr)-1)*365.0+1.0
    
enddays = -9999
if(options.endyr !=""): enddays = (int(options.endyr)+1)*365.0

tunit = str.capitalize(options.t_unit)
if(options.annually): tunit = 'YEAR'
if tunit.startswith("H"): 
    day_scaling = 24.0
elif tunit.startswith("Y"): 
    day_scaling = 1.0/365.0
else:
    day_scaling = 1.0
tunit0 = float(options.yr0)*365.0*day_scaling #the simulation year 0 timing (in day)

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

# read-in actual forcings from ELM simulation outputs
if (options.clmncheader != '' and options.met_type == 'ELM'):
    nvars = len(varnames)
    tvarname = 'time'    # variable name for time/timing
    varnames_print = False
    adspinup = False
    nx, ny, nlgrnd, nldcmp, ncolumn, npft, varsdata, varsdims, vars_tunits = \
        CLM_NcRead_1simulation(options.clm_odir, \
                           options.clmncheader, \
                           'h0', \
                           varnames_print, \
                           varnames, \
                           startdays, enddays, \
                           adspinup)

    vars_list = list(varsdata.keys())    # var_list is in format of 'h*_varname', in which 'varname' is the real variable names

    # dimension max.
    if2dgrid = False
    if('topo' in varsdata.keys()):
        if(len(varsdata['topo'].shape)==2): if2dgrid = True
    nxy = nx*ny

    # if given pts
    if options.lon==0:
        ix = 0 # for the first points
    elif options.lon<0:
        ix = -1 # for all points
    else:
        ix=float(options.lon)
    if options.lat==0:
        iy=0  # for the first points
    elif options.lat<0:
        iy = -1 # for all points
    else:
        iy=float(options.lat)

#--------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------
# read-in metdata from CPL_BYPASS_FULL
if ('cplbypass' in options.met_type):
    cplbypass_dir=options.met_idir

    if ('GSWP3' in options.met_type):
        cplbypass_mettype='GSWP3'
        if ('v1' in options.met_type and 'daymet' in options.met_type): 
            cplbypass_mettype='GSWP3v1_daymet'
        elif ('v1' in options.met_type): 
            cplbypass_mettype='GSWP3v1'
        elif ('daymet' in options.met_type): 
            cplbypass_mettype='GSWP3_daymet'
    elif ('Site' in options.met_type): 
        cplbypass_mettype = 'Site'
    elif ('ESM' in options.met_type): 
        cplbypass_mettype='ESM'
        if ('daymet' in options.met_type): cplbypass_mettype='ESM_daymet'
    elif (options.met_header != ''):
        cplbypass_mettype = options.met_header
    lon = float(options.lon)
    lat = float(options.lat)

    # read in
    zones,zlines, varsdims, varsdata = \
        clm_metdata_cplbypass_read(cplbypass_dir,cplbypass_mettype, varnames, lons=[lon], lats=[lat])

    # assign data to plotting variables
    nx=1
    ny=0
    for iz in zlines.keys():
        if np.isscalar(zlines[iz]):
            ny=ny+1
        else:
            ny=ny+len(zlines[iz])
    if2dgrid = False
    nxy=nx*ny

    ix = 0 # for 1 point
    iy = 0    
    if options.lon<0:
        ix = -1 # for all points
    if options.lat<0:
        iy = -1 # for all points

    tvarname = 'DTIME'    # variable name for time/timing
    #cplvars =  ['DTIME','tunit','t0_datetime','LONGXY','LATIXY',
    #            'FLDS','FSDS','PRECTmms','PSRF','QBOT','TBOT','WIND']
    vars_list = list(varsdata.keys())
    t0 = 1901
    if 'daymet' in options.met_type: t0 = 1980
    tunit0 = t0*365*day_scaling # to convert from 'days-since-1901-01-01-00:00' to days since 0001-01-01 00:00
    if options.seasonally: tunit0 = 0.0 # this is better for plotting
    
    if startdays!=-9999:
        idx=np.where(varsdata[tvarname]>=startdays-t0*365.0-1.0)  # both in 'days', but one from 1901, one from 1
        varsdata[tvarname] = varsdata[tvarname][idx]
        for iv in varnames:
            if (len(varsdata[iv])>1):
                varsdata[iv] =np.moveaxis(varsdata[iv],-1,0)  # move time-axis from last to the first
            varsdata[iv] = varsdata[iv][idx]
    if enddays!=-9999:
        idx=np.where(varsdata[tvarname]<=enddays-t0*365.0)
        varsdata[tvarname] = varsdata[tvarname][idx]
        for iv in varnames:
            varsdata[iv] = varsdata[iv][idx]
    


#--------------------------------------------------------------------------------------
# read-in metdata from full met directory
if (options.met_type == 'GSWP3' and 'cplbypass' not in options.met_type):
    metdir=options.met_idir
    metfileheader=options.met_header
    met_type=options.met_type
    met_domain=options.met_domain
    lon = float(options.lon)
    lat = float(options.lat)

    # read in
    ix,iy, varsdims, varsdata = \
        clm_metdata_read(metdir,metfileheader, met_type, met_domain, lon, lat,varnames)
    # assign data to plotting variables
    nx=len(ix)
    ny=len(iy)
    if2dgrid = False
    nxy=nx*ny

    tvarname = 'time'    # variable name for time/timing
    #cplvars =  ['time','tunit','LONGXY','LATIXY',
    #            'FLDS','FSDS','PRECTmms','PSRF','QBOT','TBOT','WIND']
    vars_list = list(varsdata.keys())
    tunit0 = 1901*365*day_scaling # converted from 'days-since-1901-01-01-00:00' to days since 0001-01-01 00:00
    if options.seasonally: tunit0 = 0.0
    


#--------------------------------------------------------------------------------------
# read-in metdata from NCDC daily Tmax/Tmin, Precipitation data

if (options.met_type == 'NCDC'):
    site,odata_header,odata = \
        singleNCDCReadCsvfile('2059560_Alert_initproc.csv','NCDC_metric','-999.99') # 'metric' refers to NCDC data converted to metric system already
        #singleNCDCReadCsvfile('2022211_MysVanKarem_initproc.csv','NCDC_metric','-9999.99') # 'metric' refers to NCDC data converted to metric system already
        #singleNCDCReadCsvfile('GHCND_Alert_initproc.csv','GHCND','-999.99')
    #site_header = ("LATITUDE","LONGITUDE","ELEVATION")
    #data_header = ("YEAR","MONTH","DOY","PRCP","SNWD","TAVG","TMAX","TMIN")
    
    # need to recalculate time
    tvarname = 'time'    # variable name for time/timing
    
    varsdata = {}
    varsdata['time']=(odata['YEAR']-1900.0)*365.0+(odata['DOY']-1.0)  # days since 1900-01-01 00:00:00
    lpyrindex=[i for i, x in enumerate(odata['DOY']) if x == 366.0 ]  # NCDC data is in leap-year calender. Here simply remove last day of the year
    varsdata['time']=np.delete(varsdata['time'],lpyrindex)

    varsdims = {}
    varsdims['time']  = tvarname
    varsdims['TBOTd'] = tvarname
    varsdims['TMAXd'] = tvarname
    varsdims['TMINd'] = tvarname
    varsdims['PRCPd'] = tvarname
    varsdims['SNWDd'] = tvarname

    varsdata['TBOTd']=np.delete(odata['TAVG'],lpyrindex)+273.15
    varsdata['TMAXd']=np.delete(odata['TMAX'],lpyrindex)+273.15
    varsdata['TMINd']=np.delete(odata['TMIN'],lpyrindex)+273.15
    varsdata['PRCPd']=np.delete(odata['PRCP'],lpyrindex)
    varsdata['SNWDd']=np.delete(odata['SNWD'],lpyrindex)

    
    nx=1#len(site['LONGITUDE'])
    ny=1#len(site['LATITUDE'])
    if2dgrid = False
    nxy=nx*ny
    
    
    ix = 0
    iy = 0
    vars_list = list(varsdata.keys())
    

#--------------------------------------------------------------------------------------
# plotting
nrow = 1   # sub-plot vertically arranged number (row no.)
ncol = 1   # sub-plot horizontally arranged number (column no.)
if(nvars==2): ncol = 2
if(nvars==3): ncol = 3
if(nvars==4): nrow=2; ncol=2

ivar = 0
for var in varnames:
    print (ivar,var)
    # names for variable and its time-axis
    vars_list = list(varsdata.keys())    # var_list is in format of 'h*_varname', in which 'varname' is the real variable names
    for hv in vars_list:
        if re.search(var, hv): 
            var_h = hv
            hinc  = hv.replace(var,'')
            var_t = '%s' %hinc+tvarname
            break
        
    vdata = varsdata[var_h]
    vdims = varsdims[var_h]
    
    tt = varsdata[var_t]   # time unit: days
    t  = sorted(tt)
    it = sorted(range(len(tt)), key=lambda k: tt[k])
    nt = len(tt)

    
    #time dimension
    dim_tt = -999
    if(tvarname in vdims): 
        dim_tt = vdims.index(tvarname)
        ivar = ivar + 1 # counting vars, with time-series only
    else:
        continue
    

       
    # data series
    gdata=[]
    gdata = np.zeros((nt,nx*ny))    #temporary data holder in 2-D (tt, grids)
    
    for i in range(len(tt)):
        if((ix>=0 or iy>=0) and len(vdata.shape)>=2): # multiple grids but only 1 extracted
            if(if2dgrid): 
                if(ix<0):
                    gdata[i,:] = vdata[it[i]][iy,:]
                elif(iy<0):
                    gdata[i,:] = vdata[it[i]][:,ix]
                else:
                    gdata[i,:] = vdata[it[i]][iy,ix]
            else:
                gdata[i,:] = vdata[it[i]][max(iy,ix)]

        else: # all grids or 1 grid
            gdata[i,:] = vdata[it[i]].reshape(nx*ny)
        
        
 
    if(options.seasonally):
        t=np.asarray(t)/365.0
        dim_yr=int(math.ceil(max(t))-math.floor(min(t)))
        t=(t-np.floor(t))*365.0 # still in days
        t=t.reshape(dim_yr,-1)
        t=np.mean(t,axis=0)
        if(t[0]!=0.0): t=np.where(t==0,365.0,t)  # day 0 shall be day 365 if not first element, otherwise plotting X axis not good
        
        if(len(gdata)>0):
            shp=np.hstack(([dim_yr,-1],gdata.shape[1:]))
            gdata=gdata.reshape(shp)
            gdata=np.nanmean(gdata,axis=0)
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
            gdata=np.nanmean(gdata,axis=1)

    #plotting
    vname = varnames[varnames.index(var)]
    t = np.asarray(t)*day_scaling + tunit0
    GridedVarPlotting(plt, nrow, ncol, ivar, t, tunit, gdata, \
                        vname, '(a) All Grids')

# printing plot in PDF
ofname = 'Figure_Met.pdf'
fig = plt.gcf()
fig.set_size_inches(10.0, 6.5)
plt.savefig(ofname)
plt.show()

plt.close('all')


