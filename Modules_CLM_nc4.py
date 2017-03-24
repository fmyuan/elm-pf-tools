#!/usr/bin/env python

## Modules/Functions to Process CLM Netcdf Files
## Author: Fengming YUAN, CCSI/ESD-ORNL, Oak Ridge, TN
## Date: 2017-March 

import numpy as np
from netCDF4 import Dataset
import glob
import sys
from optparse import OptionParser

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
def CLM_NcRead_1file(ncfile, varnames_print, keep_vars, chunk_keys, adspinup):
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

            if key == 'topo':
                nx = odata[key].shape[0]
                if (len(odata[key].shape)==1):
                    ny = 1
                elif(len(odata[key].shape)==2):
                    ny = odata[key].shape[1]
                
            if key == 'levgrnd':
                nlgrnd = len(odata[key])
            if key == 'levdcmp':
                nldcmp = len(odata[key])                
            if key == 'pft':
                npft = len(odata[key])

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
def CLM_NcRead_1simulation(clm_odir, ncfileheader, varnames_print, vars, adspinup):

#--------------------------------------------------------------------------------------
# INPUTS
    #
    clmhead = clm_odir+'/'+ ncfileheader +'.clm2'
    print('nc file set : '+ clmhead + '*.nc')

    if (vars == ''):
        print('No variable name by " --varname=??? "; So print out ALL variable names ONLY')
        varnames = []
        varnames_print = True
    elif(varnames_print):
        print('Print out ALL variable names ONLY')
        varnames = []
    else:
        varnames = vars
        
    
    nvars = len(varnames)

#--------------------------------------------------------------------------------------
    # variables name dictionary
    keep = []

    # a few variables common to all for plotting
    keep.append("nstep")   # need this for time step numbers
    keep.append("time")    # need this for clm timing
    keep.append("lat")     # grid center latitude
    keep.append("lon")     # grid center longitude
    keep.append("topo")    # need this for shape of surface cells, and elevation (if data available)
    keep.append("levgrnd") # need this for the number of layers  for PHY
    keep.append("levdcmp") # need this for the number of layers for BGC
    keep.append("pft")
    keep.append("pfts1d_wtgcell")

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
    tt = []
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
                CLM_NcRead_1file(filename, varnames_print, keep, chunkdata.keys(), adspinup)
            
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
            if ('lat' not in varsdata) and ('lat' in chunkdata):
                varsdata['lat'] = chunkdata['lat']
                varsdims['lat'] = chunkdata_dims['lat']
            
            if ('lon' not in varsdata) and ('lon' in chunkdata):
                varsdata['lon'] = chunkdata['lon']
                varsdims['lon'] = chunkdata_dims['lon']

            if ('topo' not in varsdata) and ('topo' in chunkdata):
                varsdata['topo'] = chunkdata['topo']
                varsdims['topo'] = chunkdata_dims['topo']
            
            #time-series
            tt.extend(chunkdata['time'])
                    
            for ivar in varnames:
                if ivar in varsdata: 
                    tmpt = np.vstack((varsdata[ivar],chunkdata[ivar]))
                    varsdata[ivar] = tmpt
                else:
                    varsdims[ivar] = chunkdata_dims[ivar]
                    varsdata[ivar] = chunkdata[ivar]
        
            
    

#--------------------------------------------------------------------------------------
#   
    # data-sets output
    return tt, nx, ny, nlgrnd, nldcmp, npft, varsdata, varsdims


