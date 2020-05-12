#!/usr/bin/env python

## Modules/Functions to Process CLM Netcdf output Files
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
pfts_user=["not_vegetated", 
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
    odata_tunits = ''
    
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
            odata_dims[key] = f.variables[key].dimensions
            
            key_val    = np.asarray(f.variables[key])
            if(hasattr(f.variables[key], '_FillValue')):
                v_missing  = f.variables[key]._FillValue
                if isinstance(v_missing, (np.float, np.float16, np.float32)):
                    key_val[key_val==v_missing] = np.nan
            odata[key] = key_val

        else:
            continue
           

        if('time' not in odata_dims[key]): continue  #cycle the loop, if not time-series dataset       
        # timing option - we do checking thereafter
        # because we want to have those CONSTANTs read-out, some of which ONLY available in the first CLM nc file.
        tt   = np.asarray(f.variables['time'])# days since model simulation starting time
        if(odata_tunits==''): odata_tunits = f.variables['time'].units
        tdays1 = min(tt)
        tdays2 = max(tt)
        if(startdays>=0):
            if(startdays>tdays2): 
                odata = {}
                odata_dims = {}
                return odata, odata_dims, odata_tunits
            else:
                s_index = int(bisect_left(tt, startdays))
                odata[key] = odata[key][s_index:,]
             
        if(enddays>0):
            if(enddays<tdays1): 
                odata = {}
                odata_dims = {}
                return odata, odata_dims, odata_tunits
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
            if (len(dvec.shape)==2):
                odata['SOILSAT_LIQ'] = dary/dvec[None,:,:]
            elif (len(dvec.shape)==3):
                odata['SOILSAT_LIQ'] = dary/dvec[None,:,:,:]
            else:
                odata['SOILSAT_LIQ'] = dary/dvec

            odata_dims['SOILSAT_LIQ'] = odata_dims['SOILLIQ'] 
        
    denice = 917.0 #kg/m3
    if 'SOILSAT_ICE' in keep_vars:
        if 'SOILICE' in odata.keys():
            dary = np.array(odata['SOILICE'])
            dvec = np.array(porosity*dz*denice)
            if (len(dvec.shape)==2):
                odata['SOILSAT_ICE'] = dary/dvec[None,:,:]
            elif (len(dvec.shape)==3):
                odata['SOILSAT_ICE'] = dary/dvec[None,:,:,:]
            else:
                odata['SOILSAT_ICE'] = dary/dvec

            odata_dims['SOILSAT_ICE'] = odata_dims['SOILICE'] 
        
        
    # out datasets
    return odata, odata_dims, odata_tunits



# -------------------------------------------------------------------
#
# Read variable(s) from multiple CLM nc files in 1 simulation
#
def CLM_NcRead_1simulation(clm_odir, ncfileheader, ncfincl, varnames_print, \
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
    keep_const.append("landmask")
    keep_const.append("pftmask")
    keep_const.append("ZSOI")
    keep_const.append("DZSOI")
    keep_const.append("column")
    keep_const.append("cols1d_wtgcell")
    keep_const.append("cols1d_active")
    keep_const.append("cols1d_itype_lunit")
    keep_const.append("cols1d_ixy")
    keep_const.append("cols1d_jxy")
    keep_const.append("pft")
    keep_const.append("pfts1d_wtgcell")
    keep_const.append("pfts1d_active")
    keep_const.append("pfts1d_itype_veg")
    keep_const.append("pfts1d_itype_lunit")
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
    if 'SOILSAT_LIQ' in varnames and 'SOILLIQ' not in keep: keep.append('SOILLIQ') 
    if 'SOILLIQ' in varnames or VARS_ALL:
        if 'SOILSAT_LIQ' not in keep: keep.append('SOILSAT_LIQ') 
        
    if 'SOILSAT_ICE' in varnames and 'SOILICE' not in keep: keep.append('SOILICE') 
    if 'SOILICE' in varnames or VARS_ALL:
        if 'SOILSAT_ICE' not in keep: keep.append('SOILSAT_ICE') 

    # build all nc files into one or two arrays
    allfile = glob.glob("%s.h*.*.nc" % clmhead)
    fincludes = [ncfincl]#['h0','h1','h2','h3','h4','h5']
    if(len(allfile)<=0):
        sys.exit("No nc file exists - %s.h*.*.nc in: %s" %(clmhead, cwdir))
    else:
        allfile = sorted(allfile)
        print('Total Files: '+str(len(allfile)))    
    
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
                
        
            print ('Processing - ', filename)
            chunkdatai,chunkdatai_dims, vars_tunits = \
                CLM_NcRead_1file(filename, varnames_print, keep, chunkdata[finc].keys(), \
                                 startdays, enddays, adspinup, VARS_ALL)
            
            # the following are constant and shall be same for all files (i.e. only need once)
            if (nx<=0):
                if('topo' in chunkdatai):
                    ny = chunkdatai['topo'].shape[0]
                    if (len(chunkdatai['topo'].shape)==1):
                        nx = 1
                    elif(len(chunkdatai['topo'].shape)==2):
                        nx = chunkdatai['topo'].shape[1]
                elif('lat' in chunkdatai and 'lon' in chunkdatai):
                    nx = chunkdatai['lon'].size
                    ny = chunkdatai['lat'].size
                
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
                    ncol = len(chunkdatai['column'])
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
    return nx, ny, nlgrnd, nldcmp, ncol, npft, varsdata, varsdims, vars_tunits

#----------------------------------------------------------------
# processing originally-readin CLM data (1,2,3D) into t-series 1D collpased dataset
def CLMvar_1Dtseries(tt, vdims, vdata, nx, ny, ix, iy, if2dgrid, \
                    annually=False, seasonally=False):

    t  = sorted(tt)
    it = sorted(range(len(tt)), key=lambda k: tt[k])
    nt = len(tt)
    nxy = nx*ny
    
    #-------------------------------------------
    #figure out dimensions at first
    
    # Does have vertical dimension?  
    zdim_indx = -999
    if('levgrnd' in vdims): 
        zdim_indx = vdims.index('levgrnd')
        nl = vdata.shape[zdim_indx]
    elif('levdcmp' in vdims): 
        zdim_indx = vdims.index('levdcmp')
        nl = vdata.shape[zdim_indx]
    elif('levsno' in vdims): 
        zdim_indx = vdims.index('levsno')
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
            elif(ix<0 and iy<0):
                if (pft_index[0]<0 and npft>1):
                    # for all grid, must specify a PFT (that is to say: not yet support 4-D plotting)
                    print ('must specify a PFT index for grid-wised plotting')
                    system.exit()
                else:
                    if(if2dgrid): 
                        vdata = vdata[:,:,:,pft_index[0]]
                        pwt1cell = pwt1cell[:,:,pft_index[0]]
                        pft1vidx = pft1vidx[:,:,pft_index[0]]
                        pft1active = pft1active[:,:,pft_index[0]]
                    else:
                        vdata = vdata[:,:,pft_index[0]]
                        pwt1cell = pwt1cell[:,pft_index[0]]
                        pft1vidx = pft1vidx[:,pft_index[0]]
                        pft1active = pft1active[:,pft_index[0]]
                    pdim_indx = -999 # 


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
        
    #-------------------------------------------
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
        exit("Variable to be processed has 4-D dimension - NOT YET supported!")
    
    # sorting data time-series by time, tt
    for i in range(len(tt)):
                
        if(zdim_indx<0 and pdim_indx<0):# 2-D grid data
            if((ix>=0 or iy>=0) and nxy>1):
                if(if2dgrid): 
                    if(ix<0):
                        gdata[i,:] = vdata[it[i]][iy,:]
                    elif(iy<0):
                        gdata[i,:] = vdata[it[i]][:,ix]
                    else:
                        gdata[i,:] = vdata[it[i]][iy,ix]
                else:
                    gdata[i,:] = vdata[it[i]][max(iy,ix)]

            else: # all grids or only 1 grid
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
                        
                
            elif(zdim_indx >= 2 or pdim_indx>=2): # 3-D soil/pft data, z_dim/p_dim in 2 or 3 (likely x/y dims before z) 
                if(if2dgrid): 
                    sdata[i,:] = vdata[it[i]][:,iy,ix,]
                else:
                    sdata[i,:] = vdata[it[i]][max(iy,ix),]

    if(seasonally):
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
    elif(annually):
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

    #
    return t, gdata, sdata, zdim_indx, pdim_indx

#----------------------------------------------------------------

