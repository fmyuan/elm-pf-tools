#!/usr/bin/env python
import os, sys, time, math
from optparse import OptionParser
import numpy as np
from netCDF4 import Dataset
#from cmath import nan, inf
from builtins import int
from numpy import double, float32

import netcdf_modules as ncmod

#---
#---Trucating/masking domain_surface data for ELM
def masking_domain(fmksrfnc_domain, domainfile_orig='', unstructured=True):
    print('#--------------------------------------------------#')
    print("re-masking domain for ELM point or regional simulations")
    
    #--------------------------------------
    # converting from partial global to global
    lmask_dname =['latitude','longitude']
    lmask_vname = 'permafrost_region_mask'
    if (os.path.isfile(fmksrfnc_domain)==True):
        src0=Dataset(fmksrfnc_domain,'r')
        lmask = np.asarray(src0.variables[lmask_vname]).astype(np.int16)
        lx = np.asarray(src0.variables[lmask_dname[1]]).astype(np.float32)
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
        
        # 
        f = open('geoxy_gridcell_mappings.txt', 'w')
        fheader='   lon          lat            geox            geoy        i     j     g '
        f.write(fheader+'\n')
        for g in range(len(idx_g)): 
            
            # re-write xy_grid_mappings.txt    
            #'(f12.5,1x,f12.6,1x, 2(f15.1, 1x),3(I5,1x))'
            f.write('%12.5f ' % lx[idx_x[g]] )
            f.write('%12.6f ' % ly[idx_y[g]] )
            f.write('%15.5f ' % lx[idx_x[g]] )
            f.write('%15.6f ' % ly[idx_y[g]] )
            f.write('%5d ' % (idx_x[g]+1) )  #x/yidx were 0-based, but need to 1-based for the mapping file
            f.write('%5d ' % (idx_y[g]+1) )
            f.write('%5d ' % (g+1) )
            f.write('\n')
        f.close()
  
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


# nearest_neibour for earth surface using kdtree
def nearest_using_kdtree(data, latname='Latitude', lonname='Longitude',kpt=2):
    import scipy.spatial as spatial
    from math import radians, cos, sin, asin, sqrt
   #"Based on https://stackoverflow.com/q/43020919/190597"
    R = 6367000.0 # meters
    def dist_to_arclength(chord_length):
        """
        https://en.wikipedia.org/wiki/Great-circle_distance
        Convert Euclidean chord length to great circle arc length
        """
        central_angle = 2*np.arcsin(chord_length/(2.0*R)) 
        arclength = R*central_angle
        return arclength

    phi = np.deg2rad(data[latname])
    theta = np.deg2rad(data[lonname])
    data['x'] = R * np.cos(phi) * np.cos(theta)
    data['y'] = R * np.cos(phi) * np.sin(theta)
    data['z'] = R * np.sin(phi)
    points = list(zip(data['x'],data['y'],data['z']))
    tree = spatial.KDTree(points)
    distance, index = tree.query(points, k=kpt)
    
    #return nearest points other than itself
    return dist_to_arclength(distance[:, 1:]), index[:,1:]

#--------------------------------------------------------------------


parser = OptionParser()

parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../../../../ccsm_inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--nco_path", dest="nco_path", default="", \
                     help = 'NCO bin PATH, default "" ')
(options, args) = parser.parse_args()


ccsm_input = os.path.abspath(options.ccsm_input)

#---
#------------------- get cru surface data as templatein ----------------------------------

# the following are what to be modified (as tempalate or original datasets)

domainfile_orig = ccsm_input+'/share/domains/domain.clm/domain.lnd.360x720_cruncep.c20190221.nc'


#surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr1850_c180216.nc'
surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c211019.nc'

surfdynfile_orig = ccsm_input+'/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_HIST_simyr1850-2015_c211019.nc'

#get grid cells
longx_orig = np.asarray(Dataset(surffile_orig).variables['LONGXY'])[0,:]
latiy_orig = np.asarray(Dataset(surffile_orig).variables['LATIXY'])[:,0]


#---
# masking
fmksrfnc_domain = '/Users/f9y/Documents/Works/BenRECCAP/RECCAP2_permafrost_regions_isimip3.nc'
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

