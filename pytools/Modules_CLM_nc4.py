#!/usr/bin/env python

## Modules/Functions to Process CLM Netcdf Files
## Author: Fengming YUAN, CCSI/ESD-ORNL, Oak Ridge, TN
## Date: 2017-March 

import os, sys
import glob
import math
from optparse import OptionParser
import numpy as np
import bisect
from netCDF4 import Dataset
from _bisect import bisect_right, bisect_left

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
                     startdays, enddays, adspinup, vars_all):
    odata  = {}
    odata_dims = {}

    try:
        startdays = float(startdays)
    except ValueError:
        startdays = -9999

    try:
        enddays = float(enddays)
    except ValueError:
        enddays = -9999

    try:
        f = Dataset(ncfile,'r')
        if varnames_print: 
            print('FILE: '+ncfile+' ------- ')

    except:
        return odata,odata_dims
        
        
    # If key is on the keep list and in nc file, uniquely add to a dictionary    
    for key in f.variables.keys():
            
        # print out all variables contained in nc files
        if varnames_print:
            print (key)
            
        if (key not in keep_vars) and (not vars_all): 
            continue  #cycle the loop, if not required by users and NOT all_vars option
        if (key not in f.variables or f.variables[key].size<=0): 
            continue                 #cycle the loop, if no data 
            
        if (len(chunk_keys)<=0) or (key not in chunk_keys): # only needs to read data once, if more than one available
            odata[key]      = np.asarray(f.variables[key])                
            odata_dims[key] = f.variables[key].dimensions
        else:
            continue
           

        if('time' not in odata_dims[key]): continue  #cycle the loop, if not time-series dataset       
        # timing option - we do checking thereafter
        # because we want to have those CONSTANTs read-out, some of which ONLY available in the first CLM nc file.
        tt   = np.asarray(f.variables['time'])# days since model simulation starting time
        tdays1 = min(tt)
        tdays2 = max(tt)
        if(startdays>=0):
            if(startdays>tdays2): 
                odata = {}
                odata_dims = {}
                return odata, odata_dims
            else:
                s_index = int(bisect_left(tt, startdays))
                odata[key] = odata[key][s_index:,]
             
        if(enddays>0):
            if(enddays<tdays1): 
                odata = {}
                odata_dims = {}
                return odata, odata_dims
            else:
                e_index = int(bisect_right(tt, enddays))
                odata[key] = odata[key][:e_index,]
                
        
        
        #ad-spinup run, recal. each component of totsomc_vr, if needed
        if adspinup:
            if key in totsomc_vr:
                sub_indx = totsomc_vr.index(key)
                odata[key] = odata[key]*ad_factor[sub_indx]

            if key in totsomn_vr:
                sub_indx = totsomn_vr.index(key)
                odata[key] = odata[key]*ad_factor[sub_indx]
            
    
    # summing up of total liter/som C/N, if needed and available
    if 'TOTLITC_vr' in keep_vars:
        indx = keep_vars.index('TOTLITC_vr')
        subs = totlitc_vr
        for isub in subs:
            if isub in odata.keys():
                if keep_vars[indx] not in odata:
                    odata[keep_vars[indx]]      = odata[isub]
                    odata_dims[keep_vars[indx]] = odata_dims[isub] 
                else:
                    odata[keep_vars[indx]]      = odata[keep_vars[indx]] + odata[isub]

    if 'TOTLITN_vr' in keep_vars:
        indx = keep_vars.index('TOTLITN_vr')
        subs = totlitn_vr
        for isub in subs:
            if isub in odata.keys():
                if keep_vars[indx] not in odata:
                    odata[keep_vars[indx]]      = odata[isub]
                    odata_dims[keep_vars[indx]] = odata_dims[isub] 
                else:
                    odata[keep_vars[indx]]      = odata[keep_vars[indx]] + odata[isub]
        
    if 'TOTSOMC_vr' in keep_vars:
        indx = keep_vars.index('TOTSOMC_vr')
        subs = totsomc_vr
        for isub in subs:
            if isub in odata.keys():
                if keep_vars[indx] not in odata:
                    odata[keep_vars[indx]]      = odata[isub]
                    odata_dims[keep_vars[indx]] = odata_dims[isub] 
                else:
                    odata[keep_vars[indx]]      = odata[keep_vars[indx]] + odata[isub]

    if 'TOTSOMN_vr' in keep_vars:
        indx = keep_vars.index('TOTSOMN_vr')
        subs = totsomn_vr
        for isub in subs:
            if isub in odata.keys():
                if keep_vars[indx] not in odata:
                    odata[keep_vars[indx]]      = odata[isub]
                    odata_dims[keep_vars[indx]] = odata_dims[isub] 
                else:
                    odata[keep_vars[indx]]      = odata[keep_vars[indx]] + odata[isub]

    
    # convert 'soilliq' to 'saturation_lia', if available
    # sat = h2osoi_liq(c,j) / (watsat(c,j)*dz(c,j)*denh2o)
    denh2o = 1000.0 #kg/m3
    if 'SOILSAT_LIQ' in keep_vars:
        if 'SOILLIQ' in odata.keys():
            dary = np.array(odata['SOILLIQ'])
            dvec = np.array(porosity*dz*denh2o)
            odata['SOILSAT_LIQ'] = dary/dvec[None,:,:,:]
            odata_dims['SOILSAT_LIQ'] = odata_dims['SOILLIQ'] 
        
    denice = 917.0 #kg/m3
    if 'SOILSAT_ICE' in keep_vars:
        if 'SOILICE' in odata.keys():
            dary = np.array(odata['SOILICE'])
            dvec = np.array(porosity*dz*denice)
            odata['SOILSAT_ICE'] = dary/dvec[None,:,:,:]
            odata_dims['SOILSAT_ICE'] = odata_dims['SOILICE'] 
        
        
    # out datasets
    return odata, odata_dims



# -------------------------------------------------------------------
#
# Read variable(s) from multiple CLM nc files in 1 simulation
#
def CLM_NcRead_1simulation(clm_odir, ncfileheader, varnames_print, \
                           vars, \
                           startdays, enddays, \
                           adspinup):

#--------------------------------------------------------------------------------------
    cwdir = os.getcwd()
    
# INPUTS
    #
    clmhead = clm_odir+'/'+ ncfileheader +'.clm2'
    print('nc file set : '+ clmhead + '.*.nc')

    VARS_ALL = False
    nvars   = 0    
    if (len(vars) <= 0):
        print('No variable name by " --varname=??? "; So print out ALL variable names ONLY')
        varnames = []
        varnames_print = True
    elif(varnames_print):
        print('Print out ALL variable names ONLY')
        varnames = []
    elif(vars[0] == 'ALL'):
        VARS_ALL = True
        varnames = []
        print ('Extract ALL variables from simulation')
    else:
        varnames = vars
        nvars = len(varnames)

    try:
        startdays = float(startdays)
    except ValueError:
        startdays = -9999

    try:
        enddays = float(enddays)
    except ValueError:
        enddays = -9999

#--------------------------------------------------------------------------------------

    # a few constants common to all for plotting
    keep_const = []
    keep_const.append("lat")     # grid center latitude
    keep_const.append("lon")     # grid center longitude
    keep_const.append("topo")    # need this for shape of surface cells, and elevation (if data available)
    keep_const.append("levgrnd") # need this for the number of layers  for PHY
    keep_const.append("levdcmp") # need this for the number of layers for BGC
    keep_const.append("ZSOI")
    keep_const.append("DZSOI")
    keep_const.append("column")
    keep_const.append("cols1d_wtgcell")
    keep_const.append("cols1d_ixy")
    keep_const.append("cols1d_jxy")
    keep_const.append("pft")
    keep_const.append("pfts1d_wtgcell")
    keep_const.append("pfts1d_ixy")
    keep_const.append("pfts1d_jxy")
    keep_const.append("WATSAT")
    keep_const.append("SUCSAT")
    keep_const.append("BSW")
    keep_const.append("HKSAT")

    # variables name dictionary
    keep = keep_const
    keep.append("nstep")         # need this for time step numbers
    keep.append("time")          # need this for clm timing
    keep.extend(varnames)

    # for some variables, they are requiring a group of variables in CLM output to sum up
    if 'TOTSOMC_vr' in varnames or VARS_ALL:
        if VARS_ALL and 'TOTSOMC_vr' not in keep: keep.append('TOTSOMC_vr') 
        indx_totsomc_vr = keep.index('TOTSOMC_vr')
        keep.extend(totsomc_vr)

    if 'TOTSOMN_vr' in varnames or VARS_ALL:
        if VARS_ALL and 'TOTSOMN_vr' not in keep: keep.append('TOTSOMN_vr') 
        indx_totsomn_vr = keep.index('TOTSOMN_vr')
        keep.extend(totsomn_vr)

    if 'TOTLITC_vr' in varnames or VARS_ALL:
        if VARS_ALL and 'TOTLITC_vr' not in keep: keep.append('TOTLITC_vr') 
        indx_totlitc_vr = keep.index('TOTLITC_vr')
        keep.extend(totlitc_vr)

    if 'TOTLITN_vr' in varnames or VARS_ALL:
        if VARS_ALL and 'TOTLITN_vr' not in keep: keep.append('TOTLITN_vr') 
        indx_totlitn_vr = keep.index('TOTLITN_vr')
        keep.extend(totlitn_vr)
        
    #
    if 'SOILLIQ' in varnames or VARS_ALL:
        if 'SOILSAT_LIQ' not in keep: keep.append('SOILSAT_LIQ') 
        
    if 'SOILICE' in varnames or VARS_ALL:
        if 'SOILSAT_ICE' not in keep: keep.append('SOILSAT_ICE') 

    # build all nc files into one or two arrays
    fincludes = ['h0','h1','h2','h3','h4','h5']
    allfile = glob.glob("%s.h*.*.nc" % clmhead)
    if(len(allfile)<=0):
        sys.exit("No nc file exists - %s.h*.*.nc in: %s" %(clmhead, cwdir))
    
    
    fchunks  = []
    styr = -999
    endyr= -999
    if(startdays>=0): styr = math.floor(startdays/365.0)
    if(enddays>=0): endyr= math.ceil(enddays/365.0)
    for filename in allfile:
        filename = filename.split(".")
        filetime = filename[-2]                 # clm nc filename format for 'timing', separated by '-'
        #fileyyyy = filetime.split("-")[0]       # the first part is yyyy
        if fchunks.count(filename[-2]) == 0:
            #if (   (int(fileyyyy)>= styr or styr <0) \
            #   and (int(fileyyyy)<=endyr or endyr<0)  ) :   # unfortunately, 'fileyyyy' IS not always same as simulation year.
                fchunks.append(filename[-2])
    
#--------------------------------------------------------------------------------------
#
    # final datasets initialization
    varsdata = {}
    varsdims = {}

    nx     = 0
    ny     = 0
    nldcmp = 0
    nlgrnd = 0
    npft   = 0
    ncol   = 0          

    hinc_done_print = []
    
    # process CLM files
    for chunk in fchunks:  # actually one time-series in a file

        chunkdata      = {}
        chunkdata_dims = {}

        # START LOOP: try reading files (h0-h5) in ONE time-chunk for each of *.h0 - h5 files
        #   (note: the time frequency may NOT be same for h0~h5)
        for finc in fincludes:        

            # file name in full
            filename = "%s.%s.%s.nc" % (clmhead,finc,chunk)
            if (not os.path.isfile(filename)): continue

            # need only to print ONCE for each of h0-h5
            if (varnames_print):
                if finc in hinc_done_print:
                    continue
                else:
                    hinc_done_print.append(finc)
                    print("-----------------------------------------------------")


            # for each chunk of time slices, append variables to a single dictionary
            if(finc not in chunkdata):
                chunkdata[finc]  = {}
                chunkdata_dims[finc] = {}
        
            chunkdatai,chunkdatai_dims = \
                CLM_NcRead_1file(filename, varnames_print, keep, chunkdata[finc].keys(), \
                                 startdays, enddays, adspinup, VARS_ALL)
            
            # the following are constant and shall be same for all files (i.e. only need once)
            if (nx<=0 and 'topo' in chunkdatai):
                nx = chunkdatai['topo'].shape[0]
                if (len(chunkdatai['topo'].shape)==1):
                    ny = 1
                elif(len(chunkdatai['topo'].shape)==2):
                    ny = chunkdatai['topo'].shape[1]
                
            if (nlgrnd<=0 and 'levgrnd' in chunkdatai):
                nlgrnd = len(chunkdatai['levgrnd'])
            if (nldcmp<=0 and 'levdcmp' in chunkdatai):
                nldcmp = len(chunkdatai['levdcmp'])                

            if (npft<=0):
                if('pft' in chunkdatai ):
                    npft = len(chunkdatai['pft'])
                if('pfts1d_wtgcell' in chunkdatai ):
                    npft = len(chunkdatai['pfts1d_wtgcell'])
                    global pfts1d_wtgcell
                    pfts1d_wtgcell = chunkdatai['pfts1d_wtgcell']
                    global pfts1d_ixy
                    pfts1d_ixy = chunkdatai['pfts1d_ixy']
                    global pfts1d_jxy
                    pfts1d_jxy = chunkdatai['pfts1d_jxy']

            if (ncol<=0):
                if('column' in chunkdatai ):
                    ncol = len(chunkdatai['pft'])
                if('cols1d_wtgcell' in chunkdatai ):
                    ncol = len(chunkdatai['cols1d_wtgcell'])
                    global cols1d_wtgcell
                    cols1d_wtgcell = chunkdatai['cols1d_wtgcell']
                    global cols1d_ixy
                    cols1d_ixy = chunkdatai['cols1d_ixy']
                    global cols1d_jxy
                    cols1d_jxy = chunkdatai['cols1d_jxy']
                     
 
            if ('WATSAT' in chunkdatai): 
                global porosity
                porosity = chunkdatai['WATSAT']
            if ('DZSOI' in chunkdatai):
                global dz 
                dz = chunkdatai['DZSOI']

                
            if len(chunkdatai)>0: 
                for ikey in chunkdatai:
                    if ("time" in chunkdatai_dims[ikey]):
                        if (ikey in chunkdata[finc]):
                            chunkdata[finc][ikey].extend(chunkdatai[ikey])
                        else:
                            chunkdata[finc][ikey] = chunkdatai[ikey]
                            chunkdata_dims[finc][ikey] = chunkdatai_dims[ikey]
                                
                    else:
                        chunkdata[ikey] = chunkdatai[ikey]         # for CONSTANTS, only needs once for all h0-h5
                        chunkdata_dims[ikey] = chunkdatai_dims[ikey]
        # END LOOP: reading all files (h0-h5) in ONE time-chunk

        #append chunks of TIMES
        if not varnames_print: 

            for fvar in chunkdata: 
               
                #
                if len(chunkdata[fvar])>0:  
                    if (fvar not in varsdata) and (fvar not in fincludes):   # constants, excluding str of 'h0~h5'
                        varsdata[fvar] = chunkdata[fvar]
                        varsdims[fvar] = chunkdata_dims[fvar]
                    
                    elif fvar in fincludes:                              # time-series in 'h0~h5'               
                        
                        for ivar in chunkdata[fvar]:
                            
                            # chunkdata is a nested dictionary, flatten it
                            finc_var = "%s_%s"%(fvar, ivar)
                            
                            if "time" in list(chunkdata_dims[fvar][ivar]):
                                
                                if finc_var not in varsdata:
                                    varsdata[finc_var] = chunkdata[fvar][ivar]
                                    varsdims[finc_var] = chunkdata_dims[fvar][ivar]
                                
                                else:
                                                                                                        
                                    if(len(chunkdata_dims[fvar][ivar])>1):
                                        tmpt = np.vstack((varsdata[finc_var],chunkdata[fvar][ivar]))
                                    else:
                                        tmpt = np.concatenate((varsdata[finc_var],chunkdata[fvar][ivar]))
                        
                                    varsdata[finc_var] = tmpt
            

#--------------------------------------------------------------------------------------
#   
    # data-sets output
    return nx, ny, nlgrnd, nldcmp, ncol, npft, varsdata, varsdims


