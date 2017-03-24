#!/usr/bin/env python

## Modules/Functions to Process CLM Netcdf Files
## Author: Fengming YUAN, CCSI/ESD-ORNL, Oak Ridge, TN
## Date: 2017-March 

import sys
import glob
import math
from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset

# ---------------------------------------------------------------
# commonly used for all

# CLM pft names
pfts=["not_vegetated", 
      "arctic_lichen", 
      "arctic_bryophyte",
      "needleleaf_evergreen_temperate_tree",
      "needleleaf_evergreen_boreal_tree",
      "needleleaf_deciduous_boreal_tree",
      "broadleaf_evergreen_tropical_tree",
      "broadleaf_evergreen_temperate_tree",
      "broadleaf_deciduous_tropical_tree",
      "broadleaf_deciduous_temperate_tree",
      "broadleaf_deciduous_boreal_tree",
      "broadleaf_evergreen_shrub",
      "broadleaf_deciduous_temperate_shrub",
      "broadleaf_deciduous_boreal_shrub",
      "evergreen_arctic_shrub",
      "deciduous_arctic_shrub",
      "c3_arctic_sedge",
      "c3_arctic_forb",
      "c3_arctic_grass",
      "c3_non-arctic_grass",
      "c4_grass",
      "c3_crop",
      "c3_irrigated",
      "corn",
      "irrigated_corn",
      "spring_temperate_cereal",
      "irrigated_spring_temperate_cereal",
      "winter_temperate_cereal",
      "irrigated_winter_temperate_cereal",
      "soybean",
      "irrigated_soybean"]

# some needed variables NOT output from CLM outputs, but can sum-up from others
totsomc_vr = ['SOIL1C_vr', 'SOIL2C_vr', 'SOIL3C_vr', 'SOIL4C_vr']
totsomn_vr = ['SOIL1N_vr', 'SOIL2N_vr', 'SOIL3N_vr', 'SOIL4N_vr']
totlitc_vr = ['LITR1C_vr', 'LITR2C_vr', 'LITR3C_vr']
totlitn_vr = ['LITR1N_vr', 'LITR2N_vr', 'LITR3N_vr']

# factors to adjust actual values of SOMC, if ad-spinup output
ad_factor =[1,1,10,100]

# -------------------------------------------------------------------
#
# Read variable(s) from 1 CLM nc file into chunk
#
def CLM_NcRead_1file(ncfile, varnames_print, keep_vars, chunk_keys, \
                     startyr, endyr, adspinup):
    nx = 0
    ny = 0
    nldcmp = 0
    nlgrnd = 0
    npft   = 0
    odata  = {}
    odata_dims = {}

    try:
        f = Dataset(ncfile,'r')
        if varnames_print: 
            print('FILE: '+ncfile+' ------- ')

    except:
        return nx,ny,nldcmp,nlgrnd,npft,odata,odata_dims
        
        
    # If key is on the keep list and in nc file, uniquely add to a dictionary    
    for key in f.variables.keys():
            
        # print out all variables contained in nc files
        if varnames_print:
            print (key)
            
        if key not in keep_vars: continue  #cycle the loop, if not exists in nc files
            
        if (len(chunk_keys)<=0) or (key not in chunk_keys): # only needs to read data once, if more than one available
            odata[key]      = np.asarray(f.variables[key])                
            odata_dims[key] = f.variables[key].dimensions

            if key == 'topo' and key not in odata:
                nx = odata[key].shape[0]
                if (len(odata[key].shape)==1):
                    ny = 1
                elif(len(odata[key].shape)==2):
                    ny = odata[key].shape[1]
                
            if key == 'levgrnd' and key not in odata:
                nlgrnd = len(odata[key])
            if key == 'levdcmp' and key not in odata:
                nldcmp = len(odata[key])                
            if key == 'pft' and key not in odata:
                npft = len(odata[key])    
            

        # timing option - we do checking here,
        # because we want to have those CONSTANTs read-out, some of which ONLY available in the first CLM nc file.
        tt   = np.asarray(f.variables['time'])# days since model simulation starting time
        tyr1 = math.floor(min(tt)/365)
        tyr2 = math.ceil(max(tt)/365)
        if(startyr !=''):
            if(startyr>tyr2): return nx,ny,nldcmp,nlgrnd,npft,odata,odata_dims 
        if(endyr !=''):
            if(endyr<tyr1): return nx,ny,nldcmp,nlgrnd,npft,odata,odata_dims 
        
        
        #ad-spinup run, recal. each component of totsomc_vr, if needed
        if adspinup:
            if key in totsomc_vr:
                sub_indx = totsomc_vr.index(key)
                odata[key] = chunk[key]*ad_factor[sub_indx]

            if key in totsomn_vr:
                sub_indx = totsomn_vr.index(key)
                odata[key] = odata[key]*ad_factor[sub_indx]
            
    
    # summing up of total liter/som C/N, if needed and available
    if 'TOTLITC_vr' in keep_vars:
        indx = keep_vars.index('TOTLITC_vr')
        subs = totlitc_vr
        for isub in subs:
            if isub in odata.keys():
                if isub==subs[0]:
                    odata[keep_vars[indx]]      = odata[isub]
                    odata_dims[keep_vars[indx]] = odata_dims[isub] 
                else:
                    odata[keep_vars[indx]]      = odata[keep[indx]] + odata[isub]
            else:
                continue

    if 'TOTLITN_vr' in keep_vars:
        indx = keep_vars.index('TOTLITN_vr')
        subs = totlitn_vr
        for isub in subs:
            if isub in odata.keys():
                if isub==subs[0]:
                    odata[keep_vars[indx]]      = odata[isub]
                    odata_dims[keep_vars[indx]] = odata_dims[isub] 
                else:
                    odata[keep_vars[indx]]      = odata[keep[indx]] + odata[isub]
            else:
                continue
        
    if 'TOTSOMC_vr' in keep_vars:
        indx = keep_vars.index('TOTSOMC_vr')
        subs = totsomc_vr
        for isub in subs:
            if isub in odata.keys():
                if isub==subs[0]:
                    odata[keep_vars[indx]]      = odata[isub]
                    odata_dims[keep_vars[indx]] = odata_dims[isub] 
                else:
                    odata[keep_vars[indx]]      = odata[keep[indx]] + odata[isub]
            else:
                continue

    if 'TOTSOMN_vr' in keep_vars:
        indx = keep_vars.index('TOTSOMN_vr')
        subs = totsomn_vr
        for isub in subs:
            if isub in odata.keys():
                if isub==subs[0]:
                    odata[keep_vars[indx]]      = odata[isub]
                    odata_dims[keep_vars[indx]] = odata_dims[isub] 
                else:
                    odata[keep_vars[indx]]      = odata[keep[indx]] + odata[isub]
            else:
                continue

    # out datasets
    return nx, ny, nldcmp, nlgrnd, npft, odata, odata_dims



# -------------------------------------------------------------------
#
# Read variable(s) from multiple CLM nc files in 1 simulation
#
def CLM_NcRead_1simulation(clm_odir, ncfileheader, varnames_print, \
                           vars, \
                           startyr, endyr, \
                           adspinup):

#--------------------------------------------------------------------------------------
# INPUTS
    #
    clmhead = clm_odir+'/'+ ncfileheader +'.clm2'
    print('nc file set : '+ clmhead + '*.nc')

    VARS_ALL = False
    nvars   = 0    
    if (vars == ''):
        print('No variable name by " --varname=??? "; So print out ALL variable names ONLY')
        varnames = []
        varnames_print = True
    elif(varnames_print):
        print('Print out ALL variable names ONLY')
        varnames = []
    elif(vars == 'ALL'):
        VARS_ALL = True
        print ('Extract ALL variables from simulation')
    else:
        varnames = vars
        nvars = len(varnames)

#--------------------------------------------------------------------------------------
    # variables name dictionary
    keep = []

    # a few constants common to all for plotting
    keep_const = []
    keep_const.append("lat")     # grid center latitude
    keep_const.append("lon")     # grid center longitude
    keep_const.append("topo")    # need this for shape of surface cells, and elevation (if data available)
    keep_const.append("levgrnd") # need this for the number of layers  for PHY
    keep_const.append("levdcmp") # need this for the number of layers for BGC
    keep_const.append("ZSOI")
    keep_const.append("DZSOI")
    keep_const.append("pft")
    keep_const.append("pfts1d_wtgcell")

    keep.append(keep_const)
    keep.append("nstep")         # need this for time step numbers
    keep.append("time")          # need this for clm timing
    keep.extend(varnames)

    # for some variables, they are requiring a group of variables in CLM output to sum up
    if 'TOTSOMC_vr' in varnames: 
        indx_totsomn_vr = keep.index('TOTSOMC_vr')
        keep.extend(totsomc_vr)

    if 'TOTSOMN_vr' in varnames:
        indx_totsomn_vr = keep.index('TOTSOMN_vr')
        keep.extend(totsomn_vr)

    if 'TOTLITC_vr' in varnames:
        indx_totlitc_vr = keep.index('TOTLITC_vr')
        keep.extend(totlitc_vr)

    if 'TOTLITN_vr' in varnames:
        indx_totlitn_vr = keep.index('TOTLITN_vr')
        keep.extend(totlitn_vr)

    # build all nc files into one or two arrays
    fincludes = ['h0','h1','h2','h3','h4','h5']
    allfile = glob.glob("%s.h*.*.nc" % clmhead)
    fchunks  = []
    for filename in allfile:
        filename = filename.split(".")
        if fchunks.count(filename[-2]) == 0: fchunks.append(filename[-2])
    

#--------------------------------------------------------------------------------------
#
    # final datasets initialization
    nstep = []
    tt    = []
    varsdata = {}
    varsdims = {}

    nx     = 0
    ny     = 0
    nldcmp = 0
    nlgrnd = 0
    npft   = 0          

    # process CLM files
    for chunk in fchunks:  # actually one time-series in a file

        # for each chunk of time slices, append variables to a single dictionary
        chunkdata      = {}
        chunkdata_dims = {}
        for finc in fincludes:        
            # try reading all files one by one
            filename = "%s.%s.%s.nc" % (clmhead,finc,chunk)
        
            nxi,nyi,nldcmpi,nlgrndi,npfti,chunkdatai,chunkdatai_dims = \
                CLM_NcRead_1file(filename, varnames_print, keep, chunkdata.keys(), \
                                 startyr, endyr, adspinup)
            
            # the following are constant and shall be same for all files (i.e. only need once)
            if nxi>0 and nx<=0: nx = nxi
            if nyi>0 and ny<=0: ny = nyi
            if nldcmpi>0 and nldcmp<=0: nldcmp = nldcmpi
            if nlgrndi>0 and nlgrnd<=0: nlgrnd = nlgrndi
            if npfti>0 and npft<=0: npft = npfti
        
        
            if len(chunkdatai)>0: 
                if (len(chunkdata)<=0):
                    chunkdata = chunkdatai
                    chunkdata_dims = chunkdatai_dims
                else:
                    for ikey in chunkdatai:
                        if (ikey in chunkdata) and ("time" in chunkdata_dims):
                            chunkdata[ikey].extend(chunkdatai[ikey])
                        else:
                            chunkdata[ikey] = chunkdatai[ikey]
                            chunkdata_dims[ikey] = chunkdatai_dims[ikey]

        #append chunks (multiple files)
        if varnames_print: 
            sys.exit("printing out variable names is DONE")
        else:
            # constants
            for iconst in keep_const:
                if (iconst not in varsdata) and (iconst in chunkdata):
                    varsdata[iconst] = chunkdata[iconst]
                    varsdims[iconst] = chunkdata_dims[iconst]
            
            
            #time-series
            if('nstep' in chunkdata): nstep.extend(chunkdata['nstep'])
            if('time' in chunkdata): tt.extend(chunkdata['time'])
                    
            for ivar in chunkdata: 
                if ((ivar in varnames) or VARS_ALL) and "time" in chunkdata_dims[ivar]:
                    
                    if(ivar in varsdata):
                        tmpt = np.vstack((varsdata[ivar],chunkdata[ivar]))
                        varsdata[ivar] = tmpt
                    else:
                        varsdims[ivar] = chunkdata_dims[ivar]
                        varsdata[ivar] = chunkdata[ivar]
        
            
    

#--------------------------------------------------------------------------------------
#   
    # data-sets output
    return tt, nx, ny, nlgrnd, nldcmp, npft, varsdata, varsdims


