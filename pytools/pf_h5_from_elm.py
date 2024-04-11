#!/usr/bin/env python

import os, sys, time, math
import re
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from optparse import OptionParser
from Modules_CLMoutput_nc4 import CLM_NcRead_1simulation
from numpy import float64
from cmath import nan


parser = OptionParser()

parser.add_option("--workdir", dest="workdir", default="./", \
                  help = "data work directory (default = ./, i.e., under current dir)")
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
parser.add_option("--varname", dest="vars", default="", \
                  help = "variable name(s) (max. 4) to be read, separated by comma ")
parser.add_option("--startyr", dest="startyr", default="1", \
                  help="clm output starting year to plot (default = 1, i.e. first year of simulation" \
                   " and can be user-defined)")
parser.add_option("--endyr", dest="endyr", default="", \
                  help="clm output ending year to plot (default = none, i.e. end of simulation)")
parser.add_option("--tunit", dest="t_unit", default="Days", \
                  help="output time unit (default = Days, i.e. Days since start-time of simulation)")
parser.add_option("--pfh5file", dest="pfh5file", default="pfinputs_fromelm", \
                  help = "pflotran h5 file name without .h5 ")

(options, args) = parser.parse_args()

#
if (options.clm_odir == './'):
    print('data directory is the current')
else:
    print('data directory: '+ options.clm_odir)

if (options.ncfileheader == ''):
    print('MUST input elm file header by " --ncfileheader=??? "' //
     ' , which usually is the portion before "*.elm.[h0-h5].*.nc"! ')
    sys.exit()
else:
    print('elm nc file header: '+ options.ncfileheader)

if (options.vars == ''):
    print('No variable name by " --varname=??? "; So print out ALL variable names')
    varnames ='TSOI:H2OSOI:SOILLIQ:QFLX_LIQ_vr'.split(':')
else:
    varnames = options.vars.split(':')
    nvars = len(varnames)

startdays = -9999
if(options.startyr !=""): startdays = (int(options.startyr)-1)*365.0
if(options.clmout_ts=='daily'): 
    startdays=startdays
elif(options.clmout_ts=='monthly'): 
    startdays=startdays+1.0
elif(options.clmout_ts=='hhourly'): 
    startdays=startdays+1.0/24.0
    
enddays = -9999
if(options.endyr !=""): enddays = int(options.endyr)*365.0

tunit = str.capitalize(options.t_unit)
if tunit.startswith("H"): 
    day_scaling = 24.0
elif tunit.startswith("Y"): 
    day_scaling = 1.0/365.0
else:
    day_scaling = 1.0

# read-in datasets from 1 simulation
nx, ny, nlgrnd, nldcmp, ncolumn, npft, varsdata, varsdims, ttunits = \
    CLM_NcRead_1simulation(options.clm_odir, \
                           options.ncfileheader, \
                           options.ncfincl, \
                           False, \
                           varnames, \
                           startdays, enddays, \
                           False)

#-------------------------------------------------------------------------
vars_list = list(varsdata.keys())    # var_list is in format of 'h*_varname', in which 'varname' is the real variable names

# time-invariant variables
constdata = {}
for v in vars_list:
    if 'time' not in varsdims[v]: 
        constdata[v]=varsdata[v]

varsdata_sorted = {}
for var in varnames:
    # names for variable and its time-axis
    for hv in vars_list:
        if re.search(var, hv): 
            var_h = hv
            hinc  = options.ncfincl+'_' #hv.replace(var,'')
            var_t = '%stime' %hinc
            break
        
    vdata = varsdata[var_h]
    vdims = varsdims[var_h]
    tt = varsdata[var_t]

    t  = sorted(tt)
    it = sorted(range(len(tt)), key=lambda k: tt[k])
    nt = len(tt)
    varsdata_sorted[var] = {}
    for i in range(nt):
        vdata_dimshift=vdata[it[i]][...]
        if 'levgrnd' in vdims: 
            vdata_dimshift=np.moveaxis(vdata_dimshift,0,-1)  # ,move 'levgrnd' to last
            if var=='SOILLIQ':
                # kg/m2 --> kg/m3
                dz=np.ones_like(vdata_dimshift)
                dz[...,:]=np.moveaxis(constdata['DZSOI'],0,-1)
                vdata_dimshift = vdata_dimshift/dz
        if 'lndgrid' in vdims:
            # normally, lndgrid is arranged along x-axis in ELM output (ni=n, nj=1)
            vdata_dimshift=np.expand_dims(vdata_dimshift, axis=0)  # (y*x,z) -->(y=1,x,z)
        vdata_dimshift=np.swapaxes(vdata_dimshift, 0, 1)  # swap (y,x,z) into (x,y,z)
        varsdata_sorted[var][i] = vdata_dimshift

#-------------------------------------------------------------------------
f0 = h5.File(options.pfh5file+'.h5', mode='w')

pf_dataset_name = ['Internal_Velocity_X(m/s)', 'Internal_Velocity_Y(m/s)', 'Internal_Velocity_Z(m/s)','tsoil(degC)','H2Osoil(kg/m3)','vwc(m3/m3)']
elm_name = ['NONE', 'NONE', 'QFLX_LIQ_vr','TSOI','SOILLIQ','H2OSOI']
unit_scalor=[1.e-3,1.e-3,1.e-3,1,1,1]
unit_offset=[0,0,0,-273.15,0,0]

for iv in range(len(pf_dataset_name)):
    if elm_name[iv] != 'NONE':
        var = elm_name[iv]
        rarray = np.zeros(tt.shape+varsdata_sorted[var][0].shape, dtype=float64)
        for i in it:
            # from ELM column vertical flow (note: no lateral flow from ELM)
            rarray[i,...] = varsdata_sorted[var][i]
        rarray[np.where(np.isnan(rarray))] = np.nan

        # PFLOTRAN h5 data grouped by dataset var name
        group = pf_dataset_name[iv]
        #'Time: '+ str("{:E}".format(tt[i])) +' d'
        f0.create_dataset(group+'/Times', tt.shape, float64, tt[it]*86400.0) # Times in unit of seconds
        
        f0.create_dataset(group+'/Data', rarray.shape, float64, 
                          rarray*unit_scalor[iv]+unit_offset[iv])

        

f0.close()

print('Done!')
