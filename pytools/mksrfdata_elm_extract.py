#!/usr/bin/env python
import os
from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset

import pytools.commons_utils.Modules_netcdf as ncmod

#---
#---Trucating/masking domain_surface data for ELM
def masking_domain(fmksrfnc_domain, lmask_dname =['nj','ni'], lmask_vname = 'mask', \
                   domainfile_orig=''):
    print('#--------------------------------------------------#')
    print("re-masking domain for ELM point or regional simulations")
    
    #--------------------------------------
    # converting from partial global to global
    #lmask_dname =['nj','ni']
    #lmask_vname = 'mask'
    if (os.path.isfile(fmksrfnc_domain)==True):
        src0=Dataset(fmksrfnc_domain,'r')
        lmask = np.asarray(src0.variables[lmask_vname]).astype(np.int16)
        if lmask_dname[1] in src0.variables.dimension():
            lx = np.asarray(src0.variables[lmask_dname[1]]).astype(np.float32)
        if lmask_dname[0] in src0.variables.dimension():
            ly = np.asarray(src0.variables[lmask_dname[0]]).astype(np.float32)
        try:
            src0_fillvalue = src0.variables[lmask_vname]._FillValue
        except:
            src0_fillvalue = np.iinfo(np.int16).min
        lmask_idx = np.where((lmask!=src0_fillvalue) & (lmask>0) & (lmask<np.iinfo(np.int16).max))
        src0.close()
        
        #
        lmask[:,:] = 0
        lmask[lmask_idx] = 1
        
        idx_tmp = np.where(lmask==1)
        idx_y   = np.asarray(idx_tmp[0])
        idx_x   = np.asarray(idx_tmp[1])
        idx_g   = np.asarray(np.where(lmask.flatten()==1))[0]
          
        # else if NOT exactly-same gridding for both 
        if (domainfile_orig==''):
            return idx_x, idx_y, idx_g, lmask
        
        elif (os.path.isfile(domainfile_orig)==True):
            # original domain.nc from which to extract must be provided
            src1=Dataset(domainfile_orig,'r')
            
            # xc, yc are lon/lat of ELM domain file
            xc = np.asarray(src1.variables['xc']).astype(np.float32)
            yc = np.asarray(src1.variables['yc']).astype(np.float32)
        
            # TODO: nearest point search and mask
            
            return idx_x, idx_y, lmask
        else:
            return


parser = OptionParser()

parser.add_option("--ccsm_input", dest="e3sm_input", \
                  default='../../../', \
                  help = "input data directory for E3SM (required), default when cwd is surfdata_map/")
parser.add_option("--fmksrf_domain", dest="fmksrf_domain", \
                  default='../../../share/domains/domain.clm/domain.nc', \
                  help = "path/filename of domain.nc for which to make surface data (required)")
parser.add_option("--nco_path", dest="nco_path", default="", \
                     help = 'NCO bin PATH, default "" ')
(options, args) = parser.parse_args()


ccsm_input = os.path.abspath(options.ccsm_input)

#---
#------------------- get cru surface data as templatein ----------------------------------

# the following are what to be modified (as tempalate or original datasets)

domainfile_orig = ccsm_input+'/share/domains/domain.clm/domain.lnd.360x720_cruncep.c20190221.nc'


#surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr1850_c180216.nc'
#surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c211019.nc'
surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c211019.nc'

surfdynfile_orig = ccsm_input+'/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_HIST_simyr1850-2015_c211019.nc'

#get grid cells
longx_orig = np.asarray(Dataset(surffile_orig).variables['LONGXY'])[0,:]
latiy_orig = np.asarray(Dataset(surffile_orig).variables['LATIXY'])[:,0]


#---
# masking
fmksrfnc_domain = options.fmksrf_domain#'/Users/f9y/Documents/Works/BenRECCAP/RECCAP2_permafrost_regions_isimip3.nc'
idx_x, idx_y, idx_g, lmask = masking_domain(fmksrfnc_domain)

unstructured=False
if unstructured:
#---
#---get new domain (xc,yc,xv,yv, mask, area), surfdata, landuse.timeseries
    print('good')
else:
    
    x0 = min(idx_x)
    xl = max(idx_x)-min(idx_x)+1
    y0 = min(idx_y)
    yl = max(idx_y)-min(idx_y)+1
    
    # domain.nc
    dim_names = ['ni','nj']
    if (options.nco_path==""):
        ncmod.nco_extract(domainfile_orig, 'domain.nc', dim_names, 
            [x0,y0], [xl,yl])
    else:
        ncmod.nco_extract(domainfile_orig, 'domain.nc', dim_names, 
            [x0,y0], [xl,yl],ncksdir=options.nco_path)
    
    ncmod.putvar('domain.nc', 'mask', lmask[y0:y0+yl-1,x0:x0+xl-1])
    
    # surfdata.nc
    dim_names = ['lsmlon', 'lsmlat']
    if (options.nco_path==""):
        ncmod.nco_extract(surffile_orig, 'surfdata.nc', dim_names, 
            [x0,y0], [xl,yl])
    else:
        ncmod.nco_extract(surffile_orig, 'surfdata.nc', dim_names, 
            [x0,y0], [xl,yl],ncksdir=options.nco_path)
    
    ncmod.putvar('surfdata.nc', 'PFTDATA_MASK', lmask[y0:y0+yl-1,x0:x0+xl-1])

    # surfdata.pftdyn.nc
    dim_names = ['lsmlon', 'lsmlat']
    if (options.nco_path==""):
        ncmod.nco_extract(surfdynfile_orig, 'landuse_timeseries.nc', dim_names, 
            [x0,y0], [xl,yl])
    else:
        ncmod.nco_extract(surfdynfile_orig, 'landuse_timeseries.nc', dim_names, 
            [x0,y0], [xl,yl],ncksdir=options.nco_path)
    
    ncmod.putvar('landuse_timeseries.nc', 'PFTDATA_MASK', lmask[y0:y0+yl-1,x0:x0+xl-1])

    print('good !')

