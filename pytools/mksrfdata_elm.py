#!/usr/bin/env python
import os, sys, math, glob
from optparse import OptionParser
import numpy as np
from scipy import interpolate
from netCDF4 import Dataset
from cmath import nan, inf
from pyproj import Transformer
from pyproj import CRS
from builtins import int
from numpy import double
from rasterio import band
try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False


import pytools.commons_utils.Modules_netcdf as ncmod

#---
#---Trucating or unstructured domain_surface data for ELM
def mksrfdata_domain(domainfile_orig, surffile_orig, surfdynfile_orig, fmksrfnc_domain, unstructured=True):
    print('#--------------------------------------------------#')
    print("Trucating domain/surface/land.timeseries data ")
    
    #--------------------------------------
    # converting from partial global to global
    dnames_elm=['lsmlat', 'lsmlon']
    dnames_=['latitude','longitude']
    lmask_name = 'permafrost_region_mask'
    if (os.path.isfile(fmksrfnc_domain)==True):
        src0=Dataset(fmksrfnc_domain,'r')
        lmask = np.asarray(src0.variables[lmask_name]).astype(np.int16)
        try:
            src0_fillvalue = src0.variables[lmask_name]._FillValue
        except:
            src0_fillvalue = np.iinfo(np.int16).min
        lmask_idx = np.where((lmask!=src0_fillvalue) & (lmask>0) & (lmask<np.iinfo(np.int16).max))
        src0.close()

        lmask[:,:] = 0
        lmask[lmask_idx] = 1
        #for trimming row/colum 
        idx_x = np.where(np.sum(lmask,0)>0)
        idx_y = np.where(np.sum(lmask,1)>0)
        
        # domain.nc
        if (os.path.isfile(domainfile_orig)==True):
            src1=Dataset(domainfile_orig,'r')
            
        


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

#
#-------------------- convert 'mksurf urban' to standard ELM surfdata -----------------------
def mksrfdata_urban(fsurfnc_all, fmksrfnc_urban_raw, redo_grid=False): 
    print('#--------------------------------------------------#')
    print('Creating surface data  - urban ...')
    fsurf_urban ='./surfdata_urban.nc'
    
    #--------------------------------------
    # region-ID to grid lat/lon
    # In raw mksrf_urban data, except for PCT_URBAN in grid cell, all other in 'region', 
    # with REGION_ID mapping to grid cell.
    fmksrfnc_new='./grided_'+fmksrfnc_urban_raw.split('/')[-1]
    if (redo_grid or os.path.isfile(fmksrfnc_new)==False):
        with Dataset(fmksrfnc_urban) as src, Dataset(fmksrfnc_new, "w") as dst:
            dst.setncatts(src.__dict__)
            
            for dname, dimension in src.dimensions.items():
                if(dname!='region'): # remove 'region' dimension
                    len_dimension = len(dimension)
                    dst.createDimension(dname, len_dimension if not dimension.isunlimited() else None)
    
            for vname in src.variables.keys():
                var = src.variables[vname]
                vdim = []
                for i in range(len(var.dimensions)):
                    vd = var.dimensions[i]
                    if vd=='region': 
                        region_dim_indx = i
                    else:
                        vdim.append(vd)
                if ('region' in var.dimensions):
                    if('lat' not in var.dimensions): vdim.append('lat')
                    if('lon' not in var.dimensions): vdim.append('lon')
                dst.createVariable(vname, var.datatype, vdim)
                dst[vname].setncatts(src[vname].__dict__)
                
                #
                varvals_dst = np.copy(dst[vname][...])   # this is just data structure with initially zero or nan or fillValue
                varvals = np.copy(src[vname][...])
                if 'region' in var.dimensions:
                    region_id = np.copy(src['REGION_ID'][...])
                    
                    # for easy work, shifting dimension 'region', 'lat', 'lon' to axis 0/1
                    if region_dim_indx > 0:
                        varvals = np.moveaxis(varvals, region_dim_indx, 0)  # note: 'moveaxis' will not change the order of other axes
                    vdim_len = len(vdim)
                    if vdim_len > 2:
                        varvals_dst = np.moveaxis(varvals_dst, -1, 0)  # after this, last'lon' shall be the first one but 'lat' moving back the last
                        varvals_dst = np.moveaxis(varvals_dst, -1, 0)  # moving last 'lat' to the first
                    
                    for rid in np.unique(region_id):
                        region1_idx = np.where(region_id==rid)
                        if rid>0:
                            varvals_dst[region1_idx[0],region1_idx[1],...] = varvals[rid-1,...]
                            
                            
                    # swap dimension of lat/lon back
                    if vdim_len > 2:
                        varvals_dst = np.moveaxis(varvals_dst, 0, -1)  # after this, 'lat' shall be the last one but 'lon' moving in the first
                        varvals_dst = np.moveaxis(varvals_dst, 0, -1)  # moving 'lon' to the last
                    
                    
                else:
                    varvals_dst = np.copy(varvals)
                dst[vname][...] = np.copy(varvals_dst)
                
            # end of for vname,var
        #
    #
    #--------------------------------------------
    # grid-cell urban inputs to standard surface data for ELM
    
    dnames_elm=['lsmlat', 'lsmlon','numurbl','nlevurb','numrad']
    dnames_urb=['lat','lon','density_class','nlevurb','numrad']
    vnames = []
    with Dataset(fmksrfnc_new) as src1, Dataset(fsurfnc_all) as src2, Dataset(fsurf_urban, "w") as dst:
        dst.setncatts(src1.__dict__)
        # new urban vars' dimensions
        for dname2, dimension2 in src2.dimensions.items():
            if dname2 in dnames_elm:                        # dim name from standard ELM surface data
                i=np.where(np.asarray(dnames_elm)==dname2)[0][0]
                dimension1 = src1.dimensions[dnames_urb[i]]
                len_dimension2 = len(dimension1)             # dim length from new data
                dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)

        # variables names to be extracted and modified
        
        for vname in src2.variables.keys():
            vdim =src2.variables[vname].dimensions
            # dim 'numurbl' is common for all urban data
            if 'numurbl' in vdim or \
               vname in ['lsmlat', 'lsmlon','LONGXY','LATIXY', 'URBAN_REGION_ID']:
                
                print(vname)
                vnames.append(vname)
                variable = src2.variables[vname]
                dst.createVariable(vname, variable.datatype, variable.dimensions)
                dst[vname].setncatts(src2[vname].__dict__)
                
                # new data values
                vname1 = vname
                sufix = ''
                if 'ALB_' in vname:
                    sufix = '_'+vname.split('_')[-1]  # standard ELM var name ending with _DIR/_DIF
                    vname1 = vname.replace(sufix,'')
                elif vname in ['lsmlat','lsmlon']:
                    vname1 = vname.replace('lsm','')
                elif vname =='URBAN_REGION_ID':
                    vname1 = 'REGION_ID'
                # checking if urban var
                variable1=src1.variables[vname1]
                if 'density_class' not in variable1.dimensions and \
                  vname1 not in ['lsmlat', 'lsmlon','LONGXY','LATIXY','REGION_ID']:
                    print('ERROR: ',vname1,'IS NOT a urban variable')
                    sys.exit(-1)
                varvals = np.copy(src1[vname1][...])
                if sufix=='_DIR': 
                    varvals_dst = varvals[0,...]
                elif sufix=='_DIF': 
                    varvals_dst = varvals[1,...]
                else:
                    varvals_dst = varvals
                dst[vname][...] = np.copy(varvals_dst)
            #end if '
        #end for
    
#--------------------------------------------------------------------
    
def mksrfdata_SoilGrid(fsurfnc_all, fmksrfnc_soilgrid_dir, var='ORGANIC', redo_grid=False):
    # var: one of 'ORGANIC', 'PCT_SAND', or 'PCT_CLAY'
    
    print('#--------------------------------------------------#')
    print("Creating surface data  - soil om & sand/clay percentage from soilGrid ...")
    fsurf_soils ='./surfdata_'+var.strip()+'.nc'
    
    #--------------------------------------
    #
    dnames_elm=['nlevsoi','lsmlat', 'lsmlon']
    dnames_soilgrid=['layer','lat','lon']
    if (redo_grid or os.path.isfile(fsurf_soils)==False):
        
        #fmksrfnc_soilgrid_dir='/Users/f9y/Documents/Works/NGEE/data/NGEEwatersheds/SoilGrids/SOMgdm3_sp_'
        layers = ['0-5cm', '5-15cm', '15-30cm', '30-60cm', '60-100cm', '100-200cm']
        zi = [0, 0.05, 0.15, 0.30, 0.60, 1.00, 2.00]  #meters
        z = zi[:-1]+np.diff(zi)/2.0     # mid-point
        z = np.insert(z,0, zi[0])           # top
        z = np.append(z,zi[-1])             # bottom. Note inclusive of 'top' and 'bottom' would be good for interpolating
        
        for i in range(len(layers)):
            ifile=fmksrfnc_soilgrid_dir+layers[i]+'.nc'
            src1=Dataset(ifile,'r')
            soil1=np.asarray(src1.variables['Band1']).astype('f4')
            src1_fillvalue = 0.0#src1.variables['Band1']._FillValue
            void1 = np.where(soil1==src1_fillvalue)
            soil1 = np.expand_dims(soil1, axis=0)
            if i==0:
                src1_data = soil1                                      # top
                src1_data = np.concatenate((src1_data, soil1), axis=0) # mid-point of 1st layer
                x=np.asarray(src1.variables['lon']).astype('f4')
                y=np.asarray(src1.variables['lat']).astype('f4')
            else:
                src1_data = np.concatenate((src1_data, soil1), axis=0)
                if i==len(layers)-1: src1_data = np.concatenate((src1_data, soil1), axis=0) # bottom
            src1.close()
        
        fn_interp = interpolate.RegularGridInterpolator((z,y,x), src1_data, bounds_error=False, fill_value=0.0)
        src1_dims={}
        src1_dims[dnames_soilgrid[0]] = z
        src1_dims[dnames_soilgrid[1]] = y
        src1_dims[dnames_soilgrid[2]] = x
        
        # ELM soil layer middle point
        ncells = 10
        jidx = np.array(range(ncells))+1 
        zsoi = 0.025*(np.exp(0.5*(jidx-0.5))-1.0)       #ELM soil layer node depths - somewhere inside a layer but not centroid
        dzsoi= np.zeros_like(zsoi)
        dzsoi[0] = 0.5*(zsoi[0]+zsoi[1])                #thickness b/n two vertical interfaces (vertices)
        for j in range(1,ncells-1):
            dzsoi[j]= 0.5*(zsoi[j+1]-zsoi[j-1])
        dzsoi[ncells-1] = zsoi[ncells-1]-zsoi[ncells-2]
        
        #soilv = np.empty((ncells,len(y),len(x)), dtype=float)
        soilgrid=np.meshgrid(zsoi,y,x, indexing='ij')
        soilgrid_list=np.reshape(soilgrid,(3,-1),order='C').T
        src1_data = fn_interp(soilgrid_list)
        if var=='ORGANIC':
            scaling=0.1  # som from SoilGrid, unit seems to be hg/m3 but shown as g/dm3
            src1_data = src1_data/0.58*scaling   # 0.58gC/gOM
            src1_data[src1_data<=0.0]=0.0
        if var=='PCT_SAND' or var=='PCT_CLAY':
            scaling=0.1 # g/kg to percentage
            src1_data = src1_data*scaling 
        
        src1_data=np.reshape(src1_data,(len(zsoi),len(y),len(x)))
        del fn_interp
        
        with Dataset(fsurfnc_all,'r') as src2, Dataset(fsurf_soils, "w") as dst:
            
            #
            for dname2, dimension2 in src2.dimensions.items():
                if dname2 in dnames_elm:                        # dim name from standard ELM surface data
                    i=np.where(np.asarray(dnames_elm)==dname2)[0][0]
                    dimension1 = src1_dims[dnames_soilgrid[i]]
                    if dname2!='nlevsoi': 
                        len_dimension2 = len(dimension1)             # dim length from new data
                    else:
                        len_dimension2 = len(dimension2)
                    dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
            #
            vname = 'lsmlat'
            vdim = ('lsmlat')
            vtype = src1_dims['lat'].dtype
            laty=dst.createVariable(vname, vtype, vdim)
            laty.units = 'degrees_north'
            laty.standard_name = 'latitude'
            dst[vname][...] = np.copy(src1_dims['lat'][...])
            
            vname = 'lsmlon'
            vdim = ('lsmlon')
            vtype = src1_dims['lon'].dtype
            lonx=dst.createVariable(vname, vtype, vdim)
            lonx.units = 'degrees_east'
            lonx.standard_name = 'longitude'
            dst[vname][...] = np.copy(src1_dims['lon'][...])
            
            vname = 'LATIXY'
            vdim = ('lsmlat','lsmlon')
            vtype = src1_dims['lat'].dtype
            lat2d=dst.createVariable(vname, vtype, vdim, fill_value=-999.)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            vals = np.moveaxis(vals,0,1)
            vals[:,...] = np.copy(src1_dims['lat'][...])
            dst[vname][...] = np.moveaxis(vals,1,0)
            del vals
            
            vname = 'LONGXY'
            vdim = ('lsmlat','lsmlon')
            vtype = src1_dims['lon'].dtype
            lon2d=dst.createVariable(vname, vtype, vdim, fill_value=-999.)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            vals[:,...] = np.copy(src1_dims['lon'][...])
            dst[vname][...] = np.copy(vals)
            del vals
            
            vname = var #'ORGANIC', 'PCT_CLAY','PCT_SAND'
            if vname in src2.variables:
                vdim = src2[vname].dimensions
                vtype = src2[vname].dtype
                var3d=dst.createVariable(vname, vtype, vdim, fill_value=-999.)
                var3d.units = src2[vname].units
                var3d.long_name = src2[vname].long_name
                # 
                src1_data[np.where(src1_data==src1_fillvalue)]=dst.variables[vname]._FillValue
                src1_data[np.where(src1_data==np.NaN)]=dst.variables[vname]._FillValue
                src1_data[np.where(src1_data==255)]=dst.variables[vname]._FillValue
                dst[vname][...] = src1_data

#--------------------------------------------------------------------
    
def mksrfdata_soildtb(fsurfnc_all, fmksrfnc_soildtb, fmksrfnc_soildtb2='', fmksrfnc_soildtb_wt='', redo_grid=False):
    print('#--------------------------------------------------#')
    print("Creating surface data  - soil thickness: 'depth to bedrock' ...")
    fsurf_soildtb ='./surfdata_soildtb.nc'
    
    #--------------------------------------
    # converting from partial global to global
    dnames_elm=['lsmlat', 'lsmlon']
    dnames_dtb=['lat','lon']
    if (redo_grid or os.path.isfile(fsurf_soildtb)==False):
        src1=Dataset(fmksrfnc_soildtb,'r')
        dtb = np.asarray(src1.variables['Band1']).astype('f4')
        src1_fillvalue = src1.variables['Band1']._FillValue
        void1 = np.where(dtb==src1_fillvalue)
        src1.close()
        if fmksrfnc_soildtb2!='' and fmksrfnc_soildtb_wt!='':
            # if provided the second DTB datasets and weight for first-sets of data
            # NOTE: this is for original data sets, which includes 2 sets, one for upland_hillslope and another for Upland_valley&lowland
            f2=Dataset(fmksrfnc_soildtb2,'r')
            f3=Dataset(fmksrfnc_soildtb_wt,'r')
            dtb2 = np.asarray(f2.variables['Band1']).astype('f4')
            dtb_wt = np.asarray(f3.variables['Band1']).astype('f4')
            void2=np.where((dtb2==f2.variables['Band1']._FillValue) | \
                           (dtb_wt==f3.variables['Band1']._FillValue) )
            f2.close()
            f3.close()
            
            dtb = dtb*dtb_wt + dtb2*(1.0-dtb_wt)
            dtb[void2] = src1_fillvalue
            # fill in void cell with 'dtb2'
            dtb[void1]=dtb2[void1]
            del dtb_wt, dtb2, void2
        
        with Dataset(fmksrfnc_soildtb,'r') as src1, Dataset(fsurfnc_all,'r') as src2, Dataset(fsurf_soildtb, "w") as dst:
            
            # new urban vars' dimensions
            for dname2, dimension2 in src2.dimensions.items():
                if dname2 in dnames_elm:                        # dim name from standard ELM surface data
                    i=np.where(np.asarray(dnames_elm)==dname2)[0][0]
                    dimension1 = src1.dimensions[dnames_dtb[i]]
                    len_dimension2 = len(dimension1)             # dim length from new data
                    dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
            #
            vname = 'lsmlat'
            vdim = ('lsmlat')
            vtype = src1.variables['lat'].datatype
            laty=dst.createVariable(vname, vtype, vdim)
            laty.units = 'degrees_north'
            laty.standard_name = 'latitude'
            dst[vname][...] = np.copy(src1.variables['lat'][...])
            
            vname = 'lsmlon'
            vdim = ('lsmlon')
            vtype = src1.variables['lon'].datatype
            lonx=dst.createVariable(vname, vtype, vdim)
            lonx.units = 'degrees_east'
            lonx.standard_name = 'longitude'
            dst[vname][...] = np.copy(src1.variables['lon'][...])
            
            vname = 'LATIXY'
            vdim = ('lsmlat','lsmlon')
            vtype = src1.variables['lat'].datatype
            lat2d=dst.createVariable(vname, vtype, vdim, fill_value=-999.)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            vals = np.moveaxis(vals,0,1)
            vals[:,...] = np.copy(src1.variables['lat'][...])
            dst[vname][...] = np.moveaxis(vals,1,0)
            del vals
            
            vname = 'LONGXY'
            vdim = ('lsmlat','lsmlon')
            vtype = src1.variables['lon'].datatype
            lon2d=dst.createVariable(vname, vtype, vdim, fill_value=-999.)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            vals[:,...] = np.copy(src1.variables['lon'][...])
            dst[vname][...] = np.copy(vals)
            del vals
            
            vname = 'aveDTB'
            vdim = ('lsmlat','lsmlon')
            var2d=dst.createVariable(vname, vtype, vdim, fill_value=-999.)
            var2d.units = 'meters below surface'
            var2d.standard_name = 'aveDTB'
            var2d.long_name = 'mean soil depth to bedrock'
            # still void cells in 'dtb'?
            dtb[np.where(dtb==src1_fillvalue)]=dst.variables[vname]._FillValue
            dtb[np.where(dtb==np.NaN)]=dst.variables[vname]._FillValue
            dtb[np.where(dtb==255)]=dst.variables[vname]._FillValue
            dst[vname][...] = dtb

# create PCT_LAKE, PCT_GLACIER datasets from soil thickness accessory data
def mksrfdata_lndunit(fsurfnc_all, fmksrfnc_lndmask, redo_grid=False):
    print('#--------------------------------------------------#')
    print("Creating surface data  - land units: PCT_LAKE, PCT_GLACIER")
    fsurf_lunit ='./surfdata_lake_icedlnd.nc'
    
    #--------------------------------------
    # 
    dnames_elm=['lsmlat', 'lsmlon']
    dnames_lunit=['lat','lon']
    if (redo_grid or os.path.isfile(fsurf_lunit)==False):
        src1=Dataset(fmksrfnc_lndmask,'r')
        lunittype = np.asarray(src1.variables['Band1']).astype('f4')
        lat = np.asarray(src1.variables['lat'])
        lon = np.asarray(src1.variables['lon'])
        src1.close()
       
        # soil thickness land mask 
        # lunittype: 0 -ocean, 1- upland, 2 -lowland, 3 -lake, 4 - perennial ice
        # Here derive 2 vars for ELM: PCT_LAKE, PCT_GLACIER, i.e. type 3 and 4 respectively and over whole land (as 100%)
        #   STEP 1: aggregating 5x5 grids, and calculating fractions for each type
        #   STEP 2: interpolating back onto the orginal grids
        # land (non-ocean) fraction
        idx = np.where(lunittype!=0.0)
        lnd = np.copy(lunittype)
        lnd[:,:] = 0.0; lnd[idx] = 1.0
        xlen = lnd.shape[0]; xlen_new = 6
        ylen = lnd.shape[1]; ylen_new = 6
        new_shp = (int(xlen/xlen_new), xlen_new, int(ylen/ylen_new), ylen_new)
        lnd_new = lnd.reshape(new_shp).sum(axis=(1,3))  # note: both xlen/ylen are multiple of 5
        lnd_new[:,:] = lnd_new[:,:]/xlen_new/ylen_new
        lnd_false = np.where(lnd_new<=0.0)
        lnd_true = np.where(lnd_new>0.0)
        
        # lake fraction over land
        idx = np.where(lunittype==3)
        lnd[:,:] = 0.0; lnd[idx] = 1.0
        lake_new = lnd.reshape(new_shp).sum(axis=(1,3))  # note: both xlen/ylen are multiple of 5
        lake_new[:,:] = lake_new[:,:]/xlen_new/ylen_new
        lake_new[lnd_true] = lake_new[lnd_true]/lnd_new[lnd_true]*100.0
        lake_new[lnd_false]= 0.0

        # perennial ice, i.e. glacier, fraction over land
        idx = np.where(lunittype==4)
        lnd[:,:] = 0.0; lnd[idx] = 1.0
        glacier_new = lnd.reshape(new_shp).sum(axis=(1,3))  # note: both xlen/ylen are multiple of 5
        glacier_new[:,:] = glacier_new[:,:]/xlen_new/ylen_new
        glacier_new[lnd_true] = glacier_new[lnd_true]/lnd_new[lnd_true]*100.0
        glacier_new[lnd_false]= 0.0
        
        # interpolating
        lat_new = lat.reshape((int(xlen/xlen_new),xlen_new)).mean(axis=1)
        dlat_new = np.mean(np.diff(lat_new))
        lon_new = lon.reshape((int(ylen/ylen_new),ylen_new)).mean(axis=1)
        dlon_new = np.mean(np.diff(lon_new))
        finterp_lake =interpolate.interp2d(lon_new-dlon_new/2, lat_new-dlat_new/2, lake_new, kind='cubic')
        lake = finterp_lake(lon, lat)
        lake[np.where(lake<0.0)]=0.0; lake[np.where(lake>100.0)]=100.0
        del lake_new
        finterp_glacier =interpolate.interp2d(lon_new-dlon_new/2, lat_new-dlat_new/2, glacier_new, kind='cubic')
        glacier = finterp_glacier(lon, lat)
        glacier[np.where(glacier<0.0)]=0.0; glacier[np.where(glacier>100.0)]=100.0
        del glacier_new
        del lnd, lnd_new, lnd_false, lnd_true
        
       # write into nc file
        with Dataset(fmksrfnc_lndmask,'r') as src1, Dataset(fsurfnc_all,'r') as src2, Dataset(fsurf_lunit, "w") as dst:
            
            # new surfdata dimensions
            for dname2, dimension2 in src2.dimensions.items():
                if dname2 in dnames_elm:                        # dim name from standard ELM surface data
                    i=np.where(np.asarray(dnames_elm)==dname2)[0][0]
                    dimension1 = src1.dimensions[dnames_lunit[i]]
                    len_dimension2 = len(dimension1)            # dim length from new data
                    dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
            #
            vname = 'lsmlat'
            vdim = ('lsmlat')
            vtype = src1.variables['lat'].datatype
            laty=dst.createVariable(vname, vtype, vdim)
            laty.units = 'degrees_north'
            laty.standard_name = 'latitude'
            dst[vname][...] = np.copy(src1.variables['lat'][...])
            
            vname = 'lsmlon'
            vdim = ('lsmlon')
            vtype = src1.variables['lon'].datatype
            lonx=dst.createVariable(vname, vtype, vdim)
            lonx.units = 'degrees_east'
            lonx.standard_name = 'longitude'
            dst[vname][...] = np.copy(src1.variables['lon'][...])
            
            vname = 'LATIXY'
            vdim = ('lsmlat','lsmlon')
            vtype = src1.variables['lat'].datatype
            lat2d=dst.createVariable(vname, vtype, vdim, fill_value=nan)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            vals = np.moveaxis(vals,0,1)
            vals[:,...] = np.copy(src1.variables['lat'][...])
            dst[vname][...] = np.moveaxis(vals,1,0)
            del vals
            
            vname = 'LONGXY'
            vdim = ('lsmlat','lsmlon')
            vtype = src1.variables['lon'].datatype
            lon2d=dst.createVariable(vname, vtype, vdim, fill_value=nan)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            vals[:,...] = np.copy(src1.variables['lon'][...])
            dst[vname][...] = np.copy(vals)
            del vals
            
            vname = 'PCT_LAKE'
            vdim =src2.variables[vname].dimensions
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, variable.dimensions, fill_value=inf)
            dst[vname].setncatts(src2[vname].__dict__)
            lake[np.where(lunittype==0)] = dst.variables[vname]._FillValue
            dst[vname][...] = np.copy(lake)
            
            vname = 'PCT_GLACIER'
            vdim =src2.variables[vname].dimensions
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, variable.dimensions, fill_value=inf)
            dst[vname].setncatts(src2[vname].__dict__)
            glacier[np.where(lunittype==0)] = dst.variables[vname]._FillValue
            dst[vname][...] = np.copy(glacier)
            
# create LATIXY/LONGXY/TOPO/STDELEV/SLOPE datasets from DEM data, usually in geotiff format
def mksrfdata_topo(fsurfnc_all, fmksrfnc_lndgeo, bands=['aspect','esl','slope'], redo_grid=False):
    print('#--------------------------------------------------#')
    print("Creating surface data  - LATIXY, LONGXY, TOPO, STDELEV, SLOPE, AREA")
    fsurf_lndtopo ='./surfdata_topo.nc'
    
    #--------------------------------------
    # 
    #dnames_elm=['lsmlat', 'lsmlon']
    dnames_elm=['gridcell']
    dnames_lndtopo=['geox','geoy'] #should be projected [x,y]
    if (redo_grid or os.path.isfile(fsurf_lndtopo)==False):
        src1=Dataset(fmksrfnc_lndgeo,'r')
        band = np.where(np.asarray(bands)=='esl')[0][0]
        if 'esl' in src1.variables.keys():
            esl = np.asarray(src1.variables['esl']).astype('f4')
        else:
            esl = np.asarray(src1.variables['Band'+str(band+1)]).astype('f4')

        band = np.where(np.asarray(bands)=='aspect')[0][0]
        if 'aspect' in src1.variables.keys():
            aspect = np.asarray(src1.variables['aspect']).astype('f4')
        else:
            aspect = np.asarray(src1.variables['Band'+str(band+1)]).astype('f4')

        band = np.where(np.asarray(bands)=='slope')[0][0]
        if 'slope' in src1.variables.keys():
            slope = np.asarray(src1.variables['slope']).astype('f4')
        else:
            slope = np.asarray(src1.variables['Band'+str(band+1)]).astype('f4')

        area = np.copy(esl)
        area[:,:] = 25.e-6 # 5mx5m --> km^2
        x = np.asarray(src1.variables[dnames_lndtopo[0]])
        y = np.asarray(src1.variables[dnames_lndtopo[1]])
        src1.close()
       
       # averaging over a box of 20x20, i.e. 5x5m --> 100x100m in this example
        xlen = esl.shape[1]; xlen_new = 200
        ylen = esl.shape[0]; ylen_new = 200
        xlen1 = int(xlen/xlen_new)*xlen_new   #needed to truncate 
        ylen1 = int(ylen/ylen_new)*ylen_new
        new_shp = (int(ylen1/ylen_new), ylen_new, int(xlen1/xlen_new), xlen_new)
        esl_new = esl[:ylen1,:xlen1].reshape(new_shp).mean(axis=(1,3))  # note: both xlen/ylen are multiple of xlen_new
        esl_std = esl[:ylen1,:xlen1].reshape(new_shp).std(axis=(1,3))   # note: both xlen/ylen are multiple of xlen_new
        slope_new = slope[:ylen1,:xlen1].reshape(new_shp).mean(axis=(1,3))  # note: both xlen/ylen are multiple of xlen_new
        area_new = area[:ylen1,:xlen1].reshape(new_shp).mean(axis=(1,3))  # note: both xlen/ylen are multiple of xlen_new
        
        # x/y projection to lat/lon 
        # NAD83/Alaska Albers, for AK Seward Peninsula, ifsar 5-m DEM data
        #Proj4 = +proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
        geoxy_proj_str = "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
        geoxyProj = CRS.from_proj4(geoxy_proj_str)
        # EPSG: 4326
        # Proj4: +proj=longlat +datum=WGS84 +no_defs
        lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
        Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)

        x_new = x[:xlen1].reshape((int(xlen1/xlen_new),xlen_new)).mean(axis=1)
        y_new = y[:ylen1].reshape((int(ylen1/ylen_new),ylen_new)).mean(axis=1)
        xx_new, yy_new = np.meshgrid(x_new, y_new)
        lon_new,lat_new = Txy2lonlat.transform(xx_new,yy_new)
        ij=np.where(lon_new<0.0)
        if(len(ij[0])>0): lon_new[ij]=lon_new[ij]+360.0 # for convenience, longitude from 0~360
        
        # write into nc file
        with Dataset(fmksrfnc_topo,'r') as src1, Dataset(fsurfnc_all,'r') as src2, Dataset(fsurf_lndtopo, "w") as dst:
            
            # new surfdata dimensions
            for dname2, dimension2 in src2.dimensions.items():
                if dname2 in dnames_elm:
                    i=np.where(np.asarray(dnames_elm)==dname2)[0][0]
                    dimension1 = src1.dimensions[dnames_lndtopo[i]]
                    len_dimension2 = len(dimension1)            # dim length from new data
                    dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
            #
            # 2-D structured grids (lat/lon)
            if 'lsmlat' in dnames_elm or 'lsmlon' in dnames_elm:
                vname = 'lsmlat'
                vdim = ('lsmlat')
                vtype = src1.variables['lat'].datatype
                laty=dst.createVariable(vname, vtype, vdim)
                laty.units = 'degrees_north'
                laty.standard_name = 'latitude'
                dst[vname][...] = np.copy(src1.variables['lat'][...])
                
                vname = 'lsmlon'
                vdim = ('lsmlon')
                vtype = src1.variables['lon'].datatype
                lonx=dst.createVariable(vname, vtype, vdim)
                lonx.units = 'degrees_east'
                lonx.standard_name = 'longitude'
                dst[vname][...] = np.copy(src1.variables['lon'][...])

                gdim = ('lsmlat','lsmlon')
            
            else:
                #1d unstructured domain
                gdim = ('gridcell')
                len_dimension2 = lat_new.size
                dst.createDimension('gridcell', len_dimension2)
                
                xdim = ('geox')
                dst.createDimension('geox', x_new.size)
                vname = 'geox'
                vtype = x_new.dtype
                x1d=dst.createVariable(vname, vtype, xdim, fill_value=nan)
                x1d.units = 'm'
                x1d.standard_name = 'geo-projected coordinate x'
                dst[vname][...] = x_new
                
                ydim = ('geoy')
                dst.createDimension('geoy', y_new.size)
                vname = 'geoy'
                vtype = y_new.dtype
                y1d=dst.createVariable(vname, vtype, ydim, fill_value=nan)
                y1d.units = 'm'
                y1d.standard_name = 'geo-projected coordinate y'
                dst[vname][...] = y_new
                
                vname = 'gridcell_jy'
                vtype = np.int32
                yy=dst.createVariable(vname, vtype, gdim, fill_value=-9999)
                yy.units = '-'
                yy.standard_name = 'geo-projected coordinate y indices, 0-based, for gridcell'
                dst[vname][...] = np.unravel_index(range(yy_new.size), yy_new.shape)[0]
                
                vname = 'gridcell_ix'
                vtype = np.int32
                xx=dst.createVariable(vname, vtype, gdim, fill_value=-9999)
                xx.units = '-'
                xx.standard_name = 'geo-projected coordinate x indices, 0-based, for gridcell'
                dst[vname][...] = np.unravel_index(range(xx_new.size), xx_new.shape)[1]
            
            vname = 'LATIXY'
            vdim = gdim
            vtype = src2.variables['LATIXY'].datatype
            lat2d=dst.createVariable(vname, vtype, vdim, fill_value=nan)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = lat_new.flatten()
                dst[vname][...] = np.copy(vals)
            else:
                vals = np.moveaxis(vals,0,1)
                vals[:,...] = np.copy(src1.variables['lat'][...])
                dst[vname][...] = np.moveaxis(vals,1,0)
            del vals
            
            vname = 'LONGXY'
            vdim = gdim
            vtype = src2.variables['LONGXY'].datatype
            lon2d=dst.createVariable(vname, vtype, vdim, fill_value=nan)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = lon_new.flatten()
            else:
                vals[:,...] = lon_new
            dst[vname][...] = np.copy(vals)
            del vals
            
            vname = 'TOPO'
            vdim = gdim
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, vdim, fill_value=nan)
            dst[vname].setncatts(src2[vname].__dict__)
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = esl_new.flatten()
            else:
                vals[:,...] = esl_new
            dst[vname][...] = np.copy(vals)
            del vals
            
            vname = 'STD_ELEV'
            vdim =gdim
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, vdim, fill_value=nan)
            dst[vname].setncatts(src2[vname].__dict__)
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = esl_std.flatten()
            else:
                vals[:,...] = esl_std
            dst[vname][...] = np.copy(vals)
            del vals
            
            vname = 'SLOPE'
            vdim = gdim
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, vdim, fill_value=nan)
            dst[vname].setncatts(src2[vname].__dict__)
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = slope_new.flatten()
            else:
                vals[:,...] = slope_new
            dst[vname][...] = np.copy(vals)
            del vals
            
            vname = 'AREA'
            vdim = gdim
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, vdim, fill_value=nan)
            dst[vname].setncatts(src2[vname].__dict__)
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = area_new.flatten()
            else:
                vals[:,...] = area_new
            dst[vname][...] = np.copy(vals)
            del vals
            
        # write indices for mapping 
        f = open('projection_elm_mappings.txt', 'w')
        fheader='   lon          lat            geox            geoy        i     j     g '
        f.write(fheader+'\n')
        xidx = np.unravel_index(range(xx_new.size), xx_new.shape)[1]
        yidx = np.unravel_index(range(yy_new.size), yy_new.shape)[0]
        for ig in range(xx_new.size):
             #'(f12.5,1x,f12.6,1x, 2(f15.1, 1x),3(I5,1x))'
            f.write('%12.5f ' % lon_new.flatten()[ig] )
            f.write('%12.6f ' % lat_new.flatten()[ig] )
            f.write('%15.1f ' % xx_new.flatten()[ig] )
            f.write('%15.1f ' % yy_new.flatten()[ig] )
            f.write('%5d ' % (xidx[ig]+1) )  #x/yidx were 0-based, but need to 1-based for the mapping file
            f.write('%5d ' % (yidx[ig]+1) )
            f.write('%5d ' % (ig+1) )
            f.write('\n')
        f.close()

# create arcticpft datasets from Jitu's dataset
def mksrfdata_arcticpft(fsurfnc_all, fmksrfnc_arcticpft, redo=False, OriginType=False):
    
    print('#--------------------------------------------------#')
    print("Creating surface data  - LATIXY, LONGXY, PCT_NAT_PFT")
    fsurf_arcticpft ='./surfdata_arcticpft.nc'
    
    
    #--------------------------------------
    #
    # datasets in original arctic pft tiff or nc file
    bandinfos={'bands':["arctic_lichen",
                    "arctic_bryophyte",
                    "arctic_forb",
                    "arctic_graminoid",
                    "arctic_wet_graminoid",
                    "arctic_evergreen_shrub",
                    "arctic_evergreen_tall_shrub",
                    "arctic_deciduous_dwarf_shrub",
                    "arctic_deciduous_low_shrub",
                    "arctic_low_to_tall_willowbirch_shrub",
                    "arctic_low_to_tall_alder_shrub",
                    "arctic_needleleaf_tree",
                    "arctic_broadleaf_tree",
                    "non_vegetated",
                    "water"],
                'pftnum': [1,2,11,12,13,5,6,7,8,9,10,3,4,0,-1]
               };
    
    if OriginType:
        # lichen as not_vegetated (0), moss/forb/graminoids as c3 arctic grass (12),
        # evergreen shrub(9), deci. boreal_shrub(11),
        # evergreen boreal tree(2), deci boreal tree (3)
        bandinfos['pftnum'] = [0,12,12,12,12,9,9,11,11,11,11,2,3,0,-1]
        natpft = np.asarray(range(17))
    else:
        natpft = np.asarray(range(max(bandinfos['pftnum'])+1)) # this is the real arcticpft order number 
    
    #--------------------------------------
    # 
    dnames_elm=['gridcell']
    dnames_lndtopo=['x','y'] #should be projected [x,y]
    if (redo or os.path.isfile(fsurf_arcticpft)==False):
        src1=Dataset(fmksrfnc_arcticpft,'r')

        
        pct_pft_orig = {}
        for v in src1.variables['pftname'][0:]:
            try:
                ib = bandinfos['bands'].index(v)
                iv = np.where(src1.variables['pftname'][0:]==v)
                pct_pft_orig[ib] = src1.variables['pftfrac'][iv][0]
            except ValueError:
                iv = -9999
            
            if v=='water':
                pct_water = src1.variables['pftfrac'][iv][0]

        x = np.asarray(src1.variables[dnames_lndtopo[0]])
        y = np.asarray(src1.variables[dnames_lndtopo[1]])
        src1.close()
       
        # assign data to elm-ordered array 
        # AND have to make sure all PCTs summed to 100%
        pct_nat_pft = np.zeros((len(natpft),len(y),len(x)),dtype=double)
        for ip in natpft: 
            iv = np.where(bandinfos['pftnum']==ip)[0]
            for i in iv: # in case having multiple classes
                if i>=0 and i in pct_pft_orig.keys(): 
                    pct_nat_pft[ip,...] = pct_nat_pft[ip,...]+pct_pft_orig[i]
        sum_pct = np.sum(pct_nat_pft,0)
        nonlnd_idx = np.where(sum_pct<=0.0)  # useful to mask later, either 'water' or 'non-land'
        lnd_idx = np.where(sum_pct>0.0)
        for ip in natpft: 
            pct_nat_pft[ip][lnd_idx] = pct_nat_pft[ip][lnd_idx]/sum_pct[lnd_idx]*100.0
        
        # water fraction redoing (TODO)
        
        
        
        # x/y projection to lat/lon 
        # NAD83/Alaska Albers, for AK Seward Peninsula, ifsar 5-m DEM data
        #Proj4 = +proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
        geoxy_proj_str = "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
        geoxyProj = CRS.from_proj4(geoxy_proj_str)
        # EPSG: 4326
        # Proj4: +proj=longlat +datum=WGS84 +no_defs
        lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
        Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)

        xx_new, yy_new = np.meshgrid(x, y)
        lon_new,lat_new = Txy2lonlat.transform(xx_new,yy_new)
        ij=np.where(lon_new<0.0)
        if(len(ij[0])>0): lon_new[ij]=lon_new[ij]+360.0 # for convenience, longitude from 0~360
        
        # write into nc file
        with Dataset(fmksrfnc_pft,'r') as src1, Dataset(fsurfnc_all,'r') as src2, Dataset(fsurf_arcticpft, "w") as dst:
            
            # new surfdata dimensions
            for dname2, dimension2 in src2.dimensions.items():
                if dname2 in dnames_elm:
                    i=np.where(np.asarray(dnames_elm)==dname2)[0][0]
                    dimension1 = src1.dimensions[dnames_lndtopo[i]]
                    len_dimension2 = len(dimension1)            # dim length from new data
                    dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
                if dname2 == 'natpft':
                    len_dimension2 = len(natpft)            # dim length from new data
                    dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
            #
            
            # pft dim
            pdim = ('natpft')
            
            # 2-D structured grids (lat/lon)
            if 'lsmlat' in dnames_elm or 'lsmlon' in dnames_elm:
                vname = 'lsmlat'
                vdim = ('lsmlat')
                vtype = src1.variables['lat'].datatype
                laty=dst.createVariable(vname, vtype, vdim)
                laty.units = 'degrees_north'
                laty.standard_name = 'latitude'
                dst[vname][...] = np.copy(src1.variables['lat'][...])
                
                vname = 'lsmlon'
                vdim = ('lsmlon')
                vtype = src1.variables['lon'].datatype
                lonx=dst.createVariable(vname, vtype, vdim)
                lonx.units = 'degrees_east'
                lonx.standard_name = 'longitude'
                dst[vname][...] = np.copy(src1.variables['lon'][...])

                gdim = ('lsmlat','lsmlon')
                pgdim = ('natpft','lsmlat','lsmlon')
            
            else:
                #1d unstructured domain
                gdim = ('gridcell')
                len_dimension2 = lat_new.size
                dst.createDimension('gridcell', len_dimension2)

                pgdim = ('natpft','gridcell')

                
                xdim = ('geox')
                dst.createDimension('geox', x.size)
                vname = 'geox'
                vtype = x.dtype
                x1d=dst.createVariable(vname, vtype, xdim, fill_value=nan)
                x1d.units = 'm'
                x1d.standard_name = 'geo-projected coordinate x'
                dst[vname][...] = x
                
                ydim = ('geoy')
                dst.createDimension('geoy', y.size)
                vname = 'geoy'
                vtype = y.dtype
                y1d=dst.createVariable(vname, vtype, ydim, fill_value=nan)
                y1d.units = 'm'
                y1d.standard_name = 'geo-projected coordinate y'
                dst[vname][...] = y
                
                vname = 'gridcell_jy'
                vtype = np.int32
                yy=dst.createVariable(vname, vtype, gdim, fill_value=-9999)
                yy.units = '-'
                yy.standard_name = 'geo-projected coordinate y indices, 0-based, for gridcell'
                dst[vname][...] = np.unravel_index(range(yy_new.size), yy_new.shape)[0]
                
                vname = 'gridcell_ix'
                vtype = np.int32
                xx=dst.createVariable(vname, vtype, gdim, fill_value=-9999)
                xx.units = '-'
                xx.standard_name = 'geo-projected coordinate x indices, 0-based, for gridcell'
                dst[vname][...] = np.unravel_index(range(xx_new.size), xx_new.shape)[1]
            
            vname = 'LATIXY'
            vdim = gdim
            vtype = src2.variables['LATIXY'].datatype
            lat2d=dst.createVariable(vname, vtype, vdim, fill_value=nan)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = lat_new.flatten()
                dst[vname][...] = np.copy(vals)
            else:
                vals = np.moveaxis(vals,0,1)
                vals[:,...] = np.copy(src1.variables['lat'][...])
                dst[vname][...] = np.moveaxis(vals,1,0)
            del vals
            
            vname = 'LONGXY'
            vdim = gdim
            vtype = src2.variables['LONGXY'].datatype
            lon2d=dst.createVariable(vname, vtype, vdim, fill_value=nan)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            if 'gridcell' in gdim:
                vals[:,...] = lon_new.flatten()
            else:
                vals[:,...] = lon_new
            dst[vname][...] = np.copy(vals)
            del vals
                        
            vname = 'natpft'
            vdim = pdim
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, vdim, fill_value=-99)
            dst[vname].setncatts(src2[vname].__dict__)
            vals = dst.variables[vname][...]
            vals[:] = natpft
            dst[vname][...] = np.copy(vals)
            del vals

            #double PCT_NAT_PFT(natpft, lsmlat, lsmlon) ;
            vname = 'PCT_NAT_PFT'
            vdim = pgdim
            variable = src2.variables[vname]
            dst.createVariable(vname, variable.datatype, vdim, fill_value=nan)
            dst[vname].setncatts(src2[vname].__dict__)
            vals = dst.variables[vname][...]
            if 'gridcell' in vdim:
                for ip in natpft:
                    vals[ip,...] = pct_nat_pft[ip].flatten()
            else:
                vals[:,...] = pct_nat_pft
            dst[vname][...] = np.copy(vals)
            del vals

            #need to re-calculate 'water' fraction from above
            vals_sum = np.sum(pct_nat_pft, axis=0)
            
            # if 'sum' is less than 100, the residue is actually 'water' in original data
            # then need to adjust PCT_NATVEG in a gridcell
            sumpft = np.ones_like(vals_sum)*100.0
            idx=np.where(vals_sum<100.0) # excluding 100 and above
            if len(idx[0])>0:
                sumpft[idx]=vals_sum[idx]
                vname = 'PCT_NATVEG'
                vdim = gdim
                vtype = src2.variables[vname].datatype
                dst.createVariable(vname, vtype, vdim)
                dst[vname].setncatts(src2[vname].__dict__)
                vals = dst.variables[vname][...]
                if 'gridcell' in vdim:
                    vals[:,...] = sumpft.flatten()
                else:
                    vals[:,...] = sumpft
                dst[vname][...] = np.copy(vals)
                del vals
                
                # 'water' in original data shall be called 'lake' in ELM land units
                vname = 'PCT_LAKE'
                vtype = src2.variables[vname].datatype
                dst.createVariable(vname, vtype, vdim)
                dst[vname].setncatts(src2[vname].__dict__)
                vals = dst.variables[vname][...]
                if 'gridcell' in vdim:
                    vals[:,...] = 100.0-sumpft.flatten()
                else:
                    vals[:,...] = 100.0-sumpft
                dst[vname][...] = np.copy(vals)
                del vals
        # 

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
#---get new domain (xc,yc,xv,yv, mask, area), surfdata, landuse.timeseries
# truncate or unstructed (if not yet)
if False:
    fmksrfnc_domain = '/Users/f9y/Documents/Works/BenRECCAP/RECCAP2_permafrost_regions_isimip3.nc'
    mksrfdata_domain(domainfile_orig, surffile_orig, surfdynfile_orig, fmksrfnc_domain, unstructured=False)

#---
#---get new TOPO (topo, std-elev, slope, aspect)
# convert geotiff to nc (if not yet)
fsrfnc_topo = ccsm_input+'/lnd/clm2/surfdata_map/high_res/ngee_SPAK/surfdata_topo_1000mx1000m_simyr2020.c231010.nc'
if False:
    fmksrfnc_topo = './ext100x100m_topo_seward.tif.nc'
    mksrfdata_topo(surffile_orig, fmksrfnc_topo, bands=['aspect','esl','slope'], redo_grid=True)
    os.system('mv surfdata_topo.nc '+fsrfnc_topo)

#---
#---get arctic pft data
#fsrfnc_arcticpft = ccsm_input+'/lnd/clm2/surfdata_map/high_res/ngee_SPAK/surfdata_arcticpft_seward_1000mx1000m_simyr2010.c231010.nc'
fsrfnc_arcticpft = ccsm_input+'/lnd/clm2/surfdata_map/high_res/ngee_SPAK/surfdata_originpft_seward_1000mx1000m_simyr2010.c231010.nc'
if False:
    fmksrfnc_pft = '5modelAvg_2010_l2_pftfraction_1000m_all.nc'
    #mksrfdata_arcticpft(surffile_orig, fmksrfnc_pft, redo=True)
    mksrfdata_arcticpft(surffile_orig, fmksrfnc_pft, redo=True, OriginType=True)
    os.system('mv surfdata_arcticpft.nc '+fsrfnc_arcticpft)


#---
#--- urban data --> surfdata ELM standard format
fsrfnc_urban = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_urban_0.05x0.05_simyr2000.c220127.nc'
if False: # edit as 'True' when needed to redo data
    fmksrfnc_urban = ccsm_input+'/lnd/clm2/surfdata_map/high_res/mksrf_urban_0.05x0.05_simyr2000.c220127.nc'
    mksrfdata_urban(surffile_orig, fmksrfnc_urban, redo_grid=True)                   # from raw data --> grided --> surfdata

#
#----
#---- soil thickness data (30 secs resolution) --> surfdata ELM standard format
fsrfnc_soildtb = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_soildtb_30x30sec_nwh.c220613.nc'
if False: # edit as 'True' when needed to redo data
    #fmksrfnc_soildtb = ccsm_input+'/lnd/clm2/surfdata_map/high_res/average_soil_and_sedimentary-deposit_thickness.nc'  # this datasets NOT really averaged one
    fmksrfnc_soildtb = ccsm_input+'/lnd/clm2/surfdata_map/high_res/upland_hill-slope_soil_thickness.nc'
    fmksrfnc_soildtb2 = ccsm_input+'/lnd/clm2/surfdata_map/high_res/upland_valley-bottom_and_lowland_sedimentary_deposit_thickness.nc'
    fmksrfnc_soildtb_wt = ccsm_input+'/lnd/clm2/surfdata_map/high_res/hill-slope_valley-bottom.nc'
    mksrfdata_soildtb(surffile_orig, fmksrfnc_soildtb, \
                      fmksrfnc_soildtb2=fmksrfnc_soildtb2, fmksrfnc_soildtb_wt=fmksrfnc_soildtb_wt,redo_grid=True)              # from raw data --> grided --> surfdata

#
#
fsrfnc_lake = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_lake_icedlnd_30x30sec.c220617.nc'
fsrfnc_glacier = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_lake_icedlnd_30x30sec.c220617.nc'
#--- land cover type included in soil thickness data (30 secs resolution) --> surfdata ELM standard format
# including: ocean(0), upland(1), lowland(2), lake(3), perennial ice (4), assuming single coverage per grid
# SO, here we aggregate them into 0.05-deg resolution, 
if False: # edit as 'True' when needed to redo data
    fmksrfnc_soildtb_landmask = ccsm_input+'/lnd/clm2/surfdata_map/high_res/mksrf_soilthk_land_cover_mask.nc'
    mksrfdata_lndunit(surffile_orig, fmksrfnc_soildtb_landmask, redo_grid=True)              # from raw data --> grided --> surfdata
#
#
fsrfnc_soilorg = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_ORGANIC_spak.nc'
#----- Soil OM density from SoilGrid
if False: # edit as 'True' when needed to redo data
    fmksrfnc_soilgrid_dirheader = ccsm_input+'/lnd/clm2/surfdata_map/high_res/SoilGrids_SPAK/SOMgdm3_sp_'
    mksrfdata_SoilGrid(surffile_orig, fmksrfnc_soilgrid_dirheader, var='ORGANIC',redo_grid=True)
    os.system('mv ./surfdata_ORGANIC.nc '+fsrfnc_soilorg)

fsrfnc_soiltexture = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_SAND_CLAY_spak.nc'
#----- Soil PCT_SAND from SoilGrid
if False: # edit as 'True' when needed to redo data
    fmksrfnc_soilgrid_dirheader = ccsm_input+'/lnd/clm2/surfdata_map/high_res/SoilGrids_SPAK/SOIL_CLAYgkg_sp_'
    mksrfdata_SoilGrid(surffile_orig, fmksrfnc_soilgrid_dirheader, var='PCT_CLAY',redo_grid=True)
    os.system('mv ./surfdata_PCT_CLAY.nc '+fsrfnc_soiltexture)

#----- Soil PCT_CLAY from SoilGrid
    fmksrfnc_soilgrid_dirheader = ccsm_input+'/lnd/clm2/surfdata_map/high_res/SoilGrids_SPAK/SOIL_SANDgkg_sp_'
    mksrfdata_SoilGrid(surffile_orig, fmksrfnc_soilgrid_dirheader, var='PCT_SAND',redo_grid=True)
    os.system('/usr/local/gcc-clang-darwin/nco/bin/ncks -A -v PCT_SAND ./surfdata_PCT_SAND.nc -o '+fsrfnc_soiltexture)

#
#



#-------####################
#----- Sync data for spatial content and resolution

lat_min = 0.0; lat_max = 90.0
lon_min = -180.0; lon_max = -15.0

UNSTRUCTURED_DOMAIN = False

#----- high-res data files (various scales and extents)
#fsrfnc_topo = ccsm_input+'/lnd/clm2/surfdata_map/data_NGEE-Council/surfdata_topo_100mx100m_simyr2020.c220721.nc'
#fsrfnc_topo = ccsm_input+'/lnd/clm2/surfdata_map/data_NGEE-Kougarok/surfdata_topo_100mx100m_simyr2020.c220721.nc'
fsrfnc_topo = ccsm_input+'/lnd/clm2/surfdata_map/high_res/ngee_SPAK/surfdata_topo_spak_ext1000m.nc'
if (os.path.exists(fsrfnc_topo)):
    # if all other high-res surfdata merged already, like following; otherwise comment the following out
    fsrfnc_urban = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_urb_lake_glacier_avedtb_natpft_0.05x0.05_nwh.c20220725.nc'
 
    fsrfnc_glacier = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_urb_lake_glacier_avedtb_natpft_0.05x0.05_nwh.c20220725.nc'
    fsrfnc_lake  = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_urb_lake_glacier_avedtb_natpft_0.05x0.05_nwh.c20220725.nc'
    
     #fsrfnc_natveg_pft = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_urb_lake_glacier_avedtb_natpft_0.05x0.05_nwh.c20220725.nc'
    
    #fsrfnc_natveg_pft = ccsm_input+'/lnd/clm2/surfdata_map/high_res/NONE' #
    #fsrfnc_natveg_pft = ccsm_input+'/lnd/clm2/surfdata_map/high_res/ngee_SPAK/surfdata_originpft_seward_1000mx1000m_simyr2010.c231010.nc'  # arctic-pft re-classified into original ELM PFTs
    fsrfnc_natveg_pft = ccsm_input+'/lnd/clm2/surfdata_map/high_res/ngee_SPAK/surfdata_arcticpft_seward_1000mx1000m_simyr2010.c231010.nc'  # totally new classes of arctic pft
    # data of arctic pfts, as in file below, include 'water' which is for lake in ELM land units 

    fsrfnc_soildtb = ccsm_input+'/lnd/clm2/surfdata_map/high_res/surfdata_soildtb_30x30sec_nwh.c220613.nc'

    fsrfnc_soilorg = ccsm_input+'/lnd/clm2/surfdata_map/high_res/ngee_SPAK/surfdata_ORGANIC_spak.nc' 
    fsrfnc_soiltexture = ccsm_input+'/lnd/clm2/surfdata_map/high_res/ngee_SPAK/surfdata_SAND_CLAY_spak.nc'
    
    #fsurf_new = 'surfdata_defaultpft_1000mx1000m.nc'   # default pft fraction
    #fsurf_new = 'surfdata_newpft_1000mx1000m.nc'   # new pft fraction but in original ELM classes
    fsurf_new = 'surfdata_arcticpft_1000mx1000m.nc' # new pft fractions AND arctic classes
    
    interp_urb=True; interp_pft=False               # for sync spatial resolutions - original are 0.05deg or 3arcmin
    interp_soiltdb=True; interp_lakeglacier=True   # for sync spatial resolutions - original are 30arcsec
    interp_soilorg=False; interp_soiltexture=False
    # 
    lat_new = np.asarray(Dataset(fsrfnc_topo).variables['LATIXY'][...])
    lon_new = np.asarray(Dataset(fsrfnc_topo).variables['LONGXY'][...])
    
    UNSTRUCTURED_DOMAIN = True

else:
    # The following is NOT always needed to do (ONCE is enough)
    # spatial content: North-Western Hemisphere, lat 0-90deg, lon -180~-15deg
    # resolution: 3minx3min, or, interpolated 30secx30sec
    #fsrfnc_natveg_pft = 'surfdata_urb_lake_glacier_avedtb_natpft_0.05x0.05_nwh.c20220725.nc'
    #fsrfnc_lake_glacier = 'surfdata_urb_lake_glacier_avedtb_natpft_0.05x0.05_nwh.c20220725.nc'
    #fsrfnc_soildtb = 'surfdata_urb_lake_glacier_avedtb_natpft_0.05x0.05_nwh.c20220725.nc'

    if False:
        # 3minx3min
        fsurf_new = 'surfdata_urb_lake_glacier_avedtb_natpft_0.05x0.05_nwh.c20220725_new.nc'
        interp_urb=False; interp_pft=False            # for sync spatial resolutions - original are 0.05deg or 3arcmin
        interp_soiltdb=True; interp_lakeglacier=True  # for sync spatial resolutions - original are 30 arcsec
        #5x5min
        lat_new = np.asarray(Dataset(fsrfnc_natveg_pft).variables['lsmlat'][...])
        lon_new = np.asarray(Dataset(fsrfnc_natveg_pft).variables['lsmlon'][...])
    
    if False:
        #30secx30sec
        fsurf_new = 'surfdata_urb_lake_glacier_avedtb_natpft_30secx30sec_nwh.c20220725.nc'
        interp_urb=True              # for sync spatial resolutions - original are 0.05deg or 3arcmin
        interp_pft=False
        interp_soiltdb=False;        # for sync spatial resolutions - original are 30arcsec
        interp_lakeglacier=False
        # 30 arcsec
        lat_new = np.asarray(Dataset(fsrfnc_soildtb).variables['lsmlat'][...])
        lon_new = np.asarray(Dataset(fsrfnc_soildtb).variables['lsmlon'][...])

# 
lon_new[np.where(lon_new>180.0)]=lon_new[np.where(lon_new>180.0)]-360.0
lat_new = lat_new[np.where((lat_new>=lat_min) & (lat_new<=lat_max))]
lon_new = lon_new[np.where((lon_new>=lon_min) & (lon_new<=lon_max))]

if (os.path.exists(fsrfnc_urban)):
    numurbl_new = Dataset(fsrfnc_urban).dimensions['numurbl']
    nlevurb_new = Dataset(fsrfnc_urban).dimensions['nlevurb']
    numrad_new  = Dataset(fsrfnc_urban).dimensions['numrad']

if (os.path.exists(fsrfnc_natveg_pft)):
    natpft_new = Dataset(fsrfnc_natveg_pft).dimensions['natpft']


'''
 Merge high-res data from above single nc files into ONE 
'''

# coarse original surface data
#surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr1850_c180216.nc'

# interpolated original surface data for high-res
surffile_orig = './default/surfdata.nc'
landuse_timeseries_orig = './default/surfdata.pftdyn.nc'

#----
len_total = len(lon_new)
if HAS_MPI4PY:
    mycomm = MPI.COMM_WORLD
    myrank = mycomm.Get_rank()
    mysize = mycomm.Get_size()
    
    len_myrank = int(math.floor(len_total/mysize))
    len_mod = int(math.fmod(len_total,mysize))
    n_myrank = np.full([mysize], np.int(1));n_myrank = np.cumsum(n_myrank)*len_myrank
    x_myrank = np.full([mysize], np.int(0));x_myrank[:len_mod] = 1
    n_myrank = n_myrank + np.cumsum(x_myrank) - 1        # ending index, starting 0, for each rank
    n0_myrank = np.hstack((0, n_myrank[0:mysize-1]+1))   # starting index, starting 0, for each rank

    lon_new = len_total[n0_myrank:n_myrank]
    fsurf_mynew = str(myrank)+fsurf_new

else:
    mycomm = 0
    myrank = 0
    mysize = 1
    fsurf_mynew = fsurf_new

if True:
#--------- write into nc file
    with Dataset(surffile_orig,'r') as src2, Dataset(fsurf_mynew, "w") as dst:
        
#------- new surfdata dimensions
        for dname2, dimension2 in src2.dimensions.items():
            len_dimension2 = len(dimension2)
            if dname2 == 'lsmlat': len_dimension2 = len(lat_new)
            if dname2 == 'lsmlon':
                if UNSTRUCTURED_DOMAIN:
                    len_dimension2 = 1
                else:
                    len_dimension2 = len(lon_new)
            if(os.path.exists(fsrfnc_urban)):
                if dname2 == 'numurbl':
                    len_dimension2 = len(numurbl_new)
                if dname2 == 'nlevurb':
                    len_dimension2 = len(nlevurb_new)
                if dname2 == 'numbrad':
                    len_dimension2 = len(numrad_new)
            if(os.path.exists(fsrfnc_natveg_pft)):
                if dname2 == 'natpft' or dname2=='lsmpft': len_dimension2 = len(natpft_new)
            
            dst.createDimension(dname2, len_dimension2 if not dimension2.isunlimited() else None)
        #
        if UNSTRUCTURED_DOMAIN:
            gdim = ('gridcell')
            
            if 'gridcell' not in src2.dimensions: 
                dst.createDimension(gdim, len(lat_new))
            vname = 'LATIXY'
            vtype = lat_new.dtype
            lat2d=dst.createVariable(vname, vtype, gdim, fill_value=nan)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            if UNSTRUCTURED_DOMAIN:
                vals[:,...] = np.copy(lat_new)
                dst[vname][...] = np.copy(vals)
            else:
                vals = np.moveaxis(vals,0,1)
                vals[:,...] = np.copy(lat_new)
                dst[vname][...] = np.asarray(np.moveaxis(vals,1,0))
            del vals
                
            vname = 'LONGXY'
            vtype = lon_new.dtype
            lon2d=dst.createVariable(vname, vtype, gdim, fill_value=nan)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            vals[:,...] = np.copy(lon_new)
            dst[vname][...] = np.copy(vals)
            del vals
            
            vnames = ['LATIXY', 'LONGXY']
            
        else:
            gdim = ('lsmlat','lsmlon')

            vname = 'lsmlat'
            vdim = ('lsmlat')
            vtype = lat_new.dtype
            laty=dst.createVariable(vname, vtype, vdim)
            laty.units = 'degrees_north'
            laty.standard_name = 'latitude'
            dst[vname][...] = np.copy(lat_new)
            
            vname = 'lsmlon'
            vdim = ('lsmlon')
            vtype = lon_new.dtype
            lonx=dst.createVariable(vname, vtype, vdim)
            lonx.units = 'degrees_east'
            lonx.standard_name = 'longitude'
            dst[vname][...] = np.copy(lon_new)
            
            vname = 'LATIXY'
            vtype = lat_new.dtype
            lat2d=dst.createVariable(vname, vtype, gdim, fill_value=nan)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            vals = np.moveaxis(vals,0,1)
            vals[:,...] = np.copy(lat_new)
            dst[vname][...] = np.asarray(np.moveaxis(vals,1,0))
            del vals
                
            vname = 'LONGXY'
            vtype = lon_new.dtype
            lon2d=dst.createVariable(vname, vtype, gdim, fill_value=nan)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            vals[:,...] = np.copy(lon_new)
            dst[vname][...] = np.copy(vals)
            del vals
            
            vnames = ['lsmlat', 'lsmlon', 'LATIXY', 'LONGXY']
        
#------ copy TOPO data directly
        if (os.path.exists(fsrfnc_topo)):
            fdata_src1 = Dataset(fsrfnc_topo)
            for vname in fdata_src1.variables.keys():
                if vname in src2.variables.keys() and vname not in vnames and \
                  ('TOPO' in vname or 'STD_ELEV' in vname or 'SLOPE' in vname):
                    print ('variable: ', vname)
                    vdim = fdata_src1.variables[vname].dimensions
                    if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                        vdim_new=[]
                        for i in range(len(vdim)):
                            if vdim[i]=='lsmlon':
                                vdim_new.append('gridcell')
                            elif vdim[i]!='lsmlat':
                                vdim_new.append(vdim[i])
                        vdim = vdim_new
                    
                    vtype = src2.variables[vname].datatype
                    dst.createVariable(vname, vtype, vdim)
                    dst[vname].setncatts(src2[vname].__dict__)
                    vnames.append(vname)
                    vdata = np.asarray(fdata_src1.variables[vname])
                    dst[vname][...] = np.copy(vdata)
                    del vdata
            fdata_src1.close()

#------ urban data
        if (os.path.isfile(fsrfnc_urban)):
            fdata_src1 = Dataset(fsrfnc_urban)
            vlat = np.asarray(fdata_src1.variables['LATIXY'][:,0])
            vlon = np.asarray(fdata_src1.variables['LONGXY'][0,:])

            # URBAN_REGION_ID first, which have no 'numurbl' dim
            vname = 'URBAN_REGION_ID'
            if vname in src2.variables.keys() and vname not in vnames:
                vdim = fdata_src1.variables[vname].dimensions
                if UNSTRUCTURED_DOMAIN: vdim = ('gridcell')
                vtype = src2.variables[vname].datatype
                dst.createVariable(vname, vtype, vdim)
                dst[vname].setncatts(src2[vname].__dict__)
                vnames.append(vname)
                
                vdata = np.asarray(fdata_src1.variables[vname])
                finterp_src1 =interpolate.interp2d(vlon, vlat, vdata, kind='linear')
                if UNSTRUCTURED_DOMAIN:
                    vdata_regionid = [finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                    vdata_regionid = np.squeeze(np.asarray(vdata_regionid))
                else:
                    vdata_regionid = finterp_src1(lon_new, lat_new)
                dst[vname][...] = np.copy(vdata_regionid)
                
            # rest urban vars, with dim of 'numurbl'
            for vname in fdata_src1.variables.keys():
                if vname in src2.variables.keys() and vname not in vnames:
                    vdim = fdata_src1.variables[vname].dimensions
                    if 'numurbl' not in vdim: continue        # not urban variable, exit for loop
                    
                    print ('variable: ', vname)
                    vdim = fdata_src1.variables[vname].dimensions
                    if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                        vdim_new=[]
                        for i in range(len(vdim)):
                            if vdim[i]=='lsmlon':
                                vdim_new.append('gridcell')
                            elif vdim[i]!='lsmlat':
                                vdim_new.append(vdim[i])
                        vdim = vdim_new
                        
                    vtype = src2.variables[vname].datatype
                    dst.createVariable(vname, vtype, vdim)
                    dst[vname].setncatts(src2[vname].__dict__)
                    vnames.append(vname)
                    
                    vdata = np.asarray(fdata_src1.variables[vname])
                    idx_void = np.where(vdata==-999.)
                    vdata_min = np.min(vdata[np.where(vdata!=-999.)])
                    vdata_max = np.max(vdata[np.where(vdata!=-999.)])
                    vdata[idx_void] = vdata_min  # -999f is bad for interpolating
                    vdata_new = np.asarray(dst.variables[vname]) # blanket at this point
                    vshp = vdata.shape
                    for i in range(vshp[-3]):
                        if 'numrad' in vdim or 'nlevurb' in vdim:
                            for j in range(vshp[-4]):
                                if interp_urb and (vname=='PCT_URBAN' or vname =='WTLUNIT_ROOF' or vname == 'WTROAD_PERV'):
                                    finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[j,i,...], kind='cubic')
                                else:
                                    # nearest
                                    finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[j,i,...], kind='linear')
                                if UNSTRUCTURED_DOMAIN: 
                                    temp = [finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                                    temp = np.squeeze(np.asarray(temp))
                                else:
                                    temp = finterp_src1(lon_new, lat_new)
                                temp[np.where(vdata_regionid<=0)]=0.0   # don't go beyond non-urban region, which causes data issue
                                vdata_new[j,i,] = np.copy(temp)
                                del temp
                            
                        else:
                                if interp_urb and (vname=='PCT_URBAN' or vname =='WTLUNIT_ROOF' or vname == 'WTROAD_PERV'):
                                    finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='cubic')
                                else:
                                    # nearest
                                    finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='linear')
                                if UNSTRUCTURED_DOMAIN: 
                                    temp = [finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                                    temp = np.squeeze(np.asarray(temp))
                                else:
                                    temp = finterp_src1(lon_new, lat_new)
                                temp[np.where(vdata_regionid<=0)]=0.0   # don't go beyond non-urban region, which causes data issue
                                vdata_new[i,] = np.copy(temp)
                                del temp
                        #
                    #
                    vdata_new[np.where(vdata_new<vdata_min)]=vdata_min
                    vdata_new[np.where(vdata_new>vdata_max)]=vdata_max
                    # limits of two types of variables (percentage or fraction)
                    if vname == 'PCT_URBAN':
                        vdata_new[np.where(vdata_new<0.0)]=0.0; vdata_new[np.where(vdata_new>100.0)]=100.0
                    elif vname == 'WTLUNIT_ROOF' or vname == 'WTROAD_PERV':
                        vdata_new[np.where(vdata_new<0.0)]=0.0; vdata_new[np.where(vdata_new>100.0)]=100.0
                        
                    #
                    dst[vname][...] = np.copy(vdata_new)
                    del vdata, vdata_new
                    #
                    #
            fdata_src1.close()
        
#------ PCT_LAKE
        if os.path.isfile(fsrfnc_lake):
            fdata_src1 = Dataset(fsrfnc_lake)
            for vname in fdata_src1.variables.keys():
                if vname in src2.variables.keys() and vname not in vnames and \
                  ('PCT_LAKE' in vname):
                    print ('variable: ', vname)
                    vdim = fdata_src1.variables[vname].dimensions
                    if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                        vdim_new=[]
                        for i in range(len(vdim)):
                            if vdim[i]=='lsmlon':
                                vdim_new.append('gridcell')
                            elif vdim[i]!='lsmlat':
                                vdim_new.append(vdim[i])
                        vdim = vdim_new
                    
                    vtype = src2.variables[vname].datatype
                    dst.createVariable(vname, vtype, vdim)
                    dst[vname].setncatts(src2[vname].__dict__)
                    vnames.append(vname)
                    
                    vdata = np.asarray(fdata_src1.variables[vname])
                    try:
                        idx_void = np.where(vdata==fdata_src1.variables[vname]._FillValue)
                        vdata_min = np.min(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata_max = np.max(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata[idx_void] = vdata_min  # nan is bad for interpolating
                    except:
                        vdata_min = 0.0
                        vdata_max = 100.0
                        
                    vlat = np.asarray(fdata_src1.variables['lsmlat'])
                    vlon = np.asarray(fdata_src1.variables['lsmlon'])
                    if interp_lakeglacier:
                        finterp_src1 =interpolate.interp2d(vlon, vlat, vdata, kind='cubic')
                    else:
                        finterp_src1 =interpolate.interp2d(vlon, vlat, vdata, kind='linear')
                    
                    if UNSTRUCTURED_DOMAIN: 
                        vdata_new = [finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                        vdata_new = np.squeeze(np.asarray(vdata_new))
                    else:
                        vdata_new = finterp_src1(lon_new, lat_new)
                        
                    vdata_new[np.where(vdata_new<0.0)]=0.0; vdata_new[np.where(vdata_new>100.0)]=100.0
                    dst[vname][...] = np.copy(vdata_new)
                    del vdata, vdata_new, vlon, vlat
            fdata_src1.close()
        
#------ PCT_GLACIER
        if os.path.isfile(fsrfnc_glacier):
            fdata_src1 = Dataset(fsrfnc_glacier)
            for vname in fdata_src1.variables.keys():
                if vname in src2.variables.keys() and vname not in vnames and \
                  ('PCT_GLACIER' in vname):
                    print ('variable: ', vname)
                    vdim = fdata_src1.variables[vname].dimensions
                    if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                        vdim_new=[]
                        for i in range(len(vdim)):
                            if vdim[i]=='lsmlon':
                                vdim_new.append('gridcell')
                            elif vdim[i]!='lsmlat':
                                vdim_new.append(vdim[i])
                        vdim = vdim_new
                    
                    vtype = src2.variables[vname].datatype
                    dst.createVariable(vname, vtype, vdim)
                    dst[vname].setncatts(src2[vname].__dict__)
                    vnames.append(vname)
                    
                    vdata = np.asarray(fdata_src1.variables[vname])
                    try:
                        idx_void = np.where(vdata==fdata_src1.variables[vname]._FillValue)
                        vdata_min = np.min(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata_max = np.max(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata[idx_void] = vdata_min  # nan is bad for interpolating
                    except:
                        vdata_min = 0.0
                        vdata_max = 100.0
                        
                    vlat = np.asarray(fdata_src1.variables['lsmlat'])
                    vlon = np.asarray(fdata_src1.variables['lsmlon'])
                    if interp_lakeglacier:
                        finterp_src1 =interpolate.interp2d(vlon, vlat, vdata, kind='cubic')
                    else:
                        finterp_src1 =interpolate.interp2d(vlon, vlat, vdata, kind='linear')
                    
                    if UNSTRUCTURED_DOMAIN: 
                        vdata_new = [finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                        vdata_new = np.squeeze(np.asarray(vdata_new))
                    else:
                        vdata_new = finterp_src1(lon_new, lat_new)
                        
                    vdata_new[np.where(vdata_new<0.0)]=0.0; vdata_new[np.where(vdata_new>100.0)]=100.0
                    dst[vname][...] = np.copy(vdata_new)
                    del vdata, vdata_new, vlon, vlat
            fdata_src1.close()

#------ soil 'aveDTB'
        if os.path.isfile(fsrfnc_soildtb):
            fdata_src1 = Dataset(fsrfnc_soildtb)
            vname = 'aveDTB'
            vdim = fdata_src1.variables[vname].dimensions
            if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                vdim_new=[]
                for i in range(len(vdim)):
                    if vdim[i]=='lsmlon':
                        vdim_new.append('gridcell')
                    elif vdim[i]!='lsmlat':
                        vdim_new.append(vdim[i])
                vdim = vdim_new
            
            if vname in src2.variables.keys():
                vtype = src2.variables[vname].datatype
            else:
                vtype = fdata_src1.variables[vname].datatype
            try:
                fillvalue = fdata_src1.variables[vname]._FillValue
            except:
                fillvalue = nan
            dst.createVariable(vname, vtype, vdim, fill_value=fillvalue)
            dst[vname].setncatts(src2[vname].__dict__)
            vnames.append(vname)
            
            print ('variable: ', vname)
            
            vdata = np.asarray(fdata_src1.variables[vname])
            try:
                idx_void = np.where(vdata==fdata_src1.variables[vname]._FillValue)
                vdata_min = np.min(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                vdata_max = np.max(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                vdata[idx_void] = vdata_min  # nan is bad for interpolating
            except:
                idx_void = np.where((np.isnan(vdata) | ~np.isfinite(vdata)))
                vdata_min = np.min(vdata[np.where(~(np.isnan(vdata) | ~np.isfinite(vdata)))])
                vdata_max = np.max(vdata[np.where(~(np.isnan(vdata) | ~np.isfinite(vdata)))])
                vdata[idx_void] = vdata_min  # nan/inf is bad for interpolating
            vlat = np.asarray(fdata_src1.variables['lsmlat'])
            vlon = np.asarray(fdata_src1.variables['lsmlon'])
            if interp_soiltdb:
                finterp_src1 =interpolate.interp2d(vlon, vlat, vdata, kind='cubic')
            else:
                finterp_src1 =interpolate.interp2d(vlon, vlat, vdata, kind='linear')
            if UNSTRUCTURED_DOMAIN: 
                vdata_new = [finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                vdata_new = np.squeeze(np.asarray(vdata_new))
            else:
                vdata_new = finterp_src1(lon_new, lat_new)
            vdata_new[np.where(vdata_new<vdata_min)] = vdata_min
            vdata_new[np.where(vdata_new>vdata_max)] = vdata_max  # 0-50.0 m are the ranges of original datasets
            dst[vname][...] = np.copy(vdata_new)
            del vdata, vdata_new, vlat, vlon
            
            fdata_src1.close()
            
#------- soil organic, clay, sand 
        if os.path.isfile(fsrfnc_soilorg):
            fdata_src1 = Dataset(fsrfnc_soilorg)
            for vname in fdata_src1.variables.keys():
                if vname in src2.variables.keys() and vname not in vnames and \
                  ('ORGANIC' in vname):
                    print ('variable: ', vname)
                    vdim = fdata_src1.variables[vname].dimensions
                    if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                        vdim_new=[]
                        for i in range(len(vdim)):
                            if vdim[i]=='lsmlon':
                                vdim_new.append('gridcell')
                            elif vdim[i]!='lsmlat':
                                vdim_new.append(vdim[i])
                        vdim = vdim_new

                    vtype = src2.variables[vname].datatype
                    dst.createVariable(vname, vtype, vdim)
                    dst[vname].setncatts(src2[vname].__dict__)
                    vnames.append(vname)
                    
                    vdata = np.asarray(fdata_src1.variables[vname])
                    try:
                        idx_void = np.where(vdata==fdata_src1.variables[vname]._FillValue)
                        vdata_min = np.min(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata_max = np.max(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata[idx_void] = vdata_min  # nan is bad for interpolating
                    except:
                        idx_void = np.where((np.isnan(vdata) | ~np.isfinite(vdata)))
                        if (len(idx_void[0])>0):
                            vdata_min = np.min(vdata[np.where(~(np.isnan(vdata) | ~np.isfinite(vdata)))])
                            vdata_max = np.max(vdata[np.where(~(np.isnan(vdata) | ~np.isfinite(vdata)))])
                            vdata[idx_void] = vdata_min  # nan/inf is bad for interpolating
                        else:
                            vdata_min = 0.0
                            vdata_max = 129.0    # 130 is the max. in ELM model
                    vdata[vdata<vdata_min] = vdata_min
                    vdata[vdata>vdata_max] = vdata_max
                    vlat = np.asarray(fdata_src1.variables['lsmlat'])
                    vlon = np.asarray(fdata_src1.variables['lsmlon'])
                    vdata_new = np.asarray(dst.variables[vname][...])
                    for i in range(vdata.shape[0]): #dim 'nlevsoi' in axis 0, lat/lon in the last 2
                        if interp_soilorg:
                            finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='cubic')
                        else:
                            # nearest
                            finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='linear')
                        
                        if UNSTRUCTURED_DOMAIN: 
                            temp =[finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                            temp = np.squeeze(np.asarray(temp))
                        else:
                            temp = finterp_src1(lon_new, lat_new)
                        vdata_new[i,] = np.copy(temp)
                        del temp
                    #
                    vdata_new[np.where(vdata_new<vdata_min)] = vdata_min
                    vdata_new[np.where(vdata_new>vdata_max)] = vdata_max
                    dst[vname][...] = np.copy(vdata_new)
                    del vdata, vdata_new, vlat, vlon
            fdata_src1.close()
            
#------ soil texture (clay, sand, maybe gravel) 
        if os.path.isfile(fsrfnc_soiltexture):
            fdata_src1 = Dataset(fsrfnc_soiltexture)
            for vname in fdata_src1.variables.keys():
                if vname in src2.variables.keys() and vname not in vnames and \
                  ('PCT_SAND' in vname or 'PCT_CLAY' in vname):
                    print ('variable: ', vname)
                    vdim = fdata_src1.variables[vname].dimensions
                    if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                        vdim_new=[]
                        for i in range(len(vdim)):
                            if vdim[i]=='lsmlon':
                                vdim_new.append('gridcell')
                            elif vdim[i]!='lsmlat':
                                vdim_new.append(vdim[i])
                        vdim = vdim_new

                    vtype = src2.variables[vname].datatype
                    dst.createVariable(vname, vtype, vdim)
                    dst[vname].setncatts(src2[vname].__dict__)
                    vnames.append(vname)
                    
                    vdata = np.asarray(fdata_src1.variables[vname])
                    try:
                        idx_void = np.where(vdata==fdata_src1.variables[vname]._FillValue)
                        vdata_min = np.min(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata_max = np.max(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata[idx_void] = vdata_min  # nan is bad for interpolating
                    except:
                        idx_void = np.where((np.isnan(vdata) | ~np.isfinite(vdata)))
                        if (len(idx_void[0])>0):
                            vdata_min = np.min(vdata[np.where(~(np.isnan(vdata) | ~np.isfinite(vdata)))])
                            vdata_max = np.max(vdata[np.where(~(np.isnan(vdata) | ~np.isfinite(vdata)))])
                            vdata[idx_void] = vdata_min  # nan/inf is bad for interpolating
                        else:
                            vdata_min = 0.0
                            vdata_max = 100.0    # 100% is the max. in ELM model
                    vdata[vdata<vdata_min] = vdata_min
                    vdata[vdata>vdata_max] = vdata_max
                    vlat = np.asarray(fdata_src1.variables['lsmlat'])
                    vlon = np.asarray(fdata_src1.variables['lsmlon'])
                    vdata_new = np.asarray(dst.variables[vname][...])
                    for i in range(vdata.shape[0]): #dim 'nlevsoi' in axis 0, lat/lon in the last 2
                        if interp_soiltexture:
                            finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='cubic')
                        else:
                            # nearest
                            finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='linear')
                        
                        if UNSTRUCTURED_DOMAIN: 
                            temp =[finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                            temp = np.squeeze(np.asarray(temp))
                        else:
                            temp = finterp_src1(lon_new, lat_new)
                        vdata_new[i,] = np.copy(temp)
                        del temp
                    #
                    vdata_new[np.where(vdata_new<vdata_min)] = vdata_min
                    vdata_new[np.where(vdata_new>vdata_max)] = vdata_max
                    dst[vname][...] = np.copy(vdata_new)
                    del vdata, vdata_new, vlat, vlon
            fdata_src1.close()
            
#---
#---natveg pft data
        if os.path.isfile(fsrfnc_natveg_pft):
            fdata_src1 = Dataset(fsrfnc_natveg_pft)
            for vname in fdata_src1.variables.keys():
                if vname in src2.variables.keys() and vname not in vnames and \
                  ('PCT_PFT' in vname or 'PCT_NAT_PFT' in vname):
                    
                    if 'natpft' not in vnames:
                        print('variable: ','natpft')
                        if 'natpft' in fdata_src1.variables.keys():
                            vdim = fdata_src1.variables['natpft'].dimensions
                        else:
                            vdim = src2.variables['natpft'].dimensions
                        vtype = src2.variables['natpft'].datatype
                        dst.createVariable('natpft', vtype, vdim)
                        dst['natpft'].setncatts(src2['natpft'].__dict__)
                        dst['natpft'][...] = np.asarray(range(len(natpft_new)))
                        vnames.append('natpft')
                    
                    print ('variable: ', vname)
                    vdim = fdata_src1.variables[vname].dimensions
                    if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                        vdim_new=[]
                        for i in range(len(vdim)):
                            if vdim[i]=='lsmlon':
                                vdim_new.append('gridcell')
                            elif vdim[i]!='lsmlat':
                                vdim_new.append(vdim[i])
                        vdim = vdim_new
                        
                        
                    vtype = src2.variables[vname].datatype
                    dst.createVariable(vname, vtype, vdim)
                    dst[vname].setncatts(src2[vname].__dict__)
                    vnames.append(vname)
                    
                    vdata = np.asarray(fdata_src1.variables[vname])
                    try:
                        idx_void = np.where(vdata==fdata_src1.variables[vname]._FillValue)
                        vdata_min = np.min(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata_max = np.max(vdata[np.where(vdata!=fdata_src1.variables[vname]._FillValue)])
                        vdata[idx_void] = vdata_min  # nan is bad for interpolating
                    except:
                        idx_void = np.where((np.isnan(vdata) | ~np.isfinite(vdata)))
                        if (len(idx_void[0])>0):
                            vdata_min = np.min(vdata[np.where(~(np.isnan(vdata) | ~np.isfinite(vdata)))])
                            vdata_max = np.max(vdata[np.where(~(np.isnan(vdata) | ~np.isfinite(vdata)))])
                            vdata[idx_void] = vdata_min  # nan/inf is bad for interpolating
                        else:
                            vdata_min = 0.0
                            vdata_max = 100.0
                    vdata[vdata<vdata_min] = vdata_min
                    vdata[vdata>vdata_max] = vdata_max
                    datain2D = False
                    if 'lsmlat' in fdata_src1.variables.keys() and \
                       'lsmlat' in fdata_src1.variables.keys():
                        datain2D=True
                        vlat = np.asarray(fdata_src1.variables['lsmlat'])
                        vlon = np.asarray(fdata_src1.variables['lsmlon'])
                    elif 'LATIXY' in fdata_src1.variables.keys() and \
                       'LONGXY' in fdata_src1.variables.keys():
                        datain2D=False
                        vlat = np.asarray(fdata_src1.variables['LATIXY'])
                        vlon = np.asarray(fdata_src1.variables['LONGXY'])
                    else:
                        print('natpft data dimensions must be either in 2D (lsmlat,lsmlont) or 1D non-regular (gridcell) !')
                    
                    vdata_new = np.asarray(dst.variables[vname][...])
                    for i in range(vdata.shape[0]): #dim 'pft' in axis 0, lat/lon in the last 2
                        #print('point: ',i, datain2D, vdata.shape)
                        if datain2D:
                            if interp_pft:
                                finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='cubic')
                            else:
                                # nearest
                                finterp_src1 =interpolate.interp2d(vlon, vlat, vdata[i,...], kind='linear')
                            
                            if UNSTRUCTURED_DOMAIN: 
                                temp =[finterp_src1(xx,yy) for xx,yy in zip(lon_new,lat_new)]
                                temp = np.squeeze(np.asarray(temp))
                            else:
                                temp = finterp_src1(lon_new, lat_new)
                        else:
                            # nearest search
                            len_new = len(lon_new)
                            allpoints={}
                            allpoints['Latitude']=np.append(lat_new, vlat)
                            allpoints['Longitude']=np.append(lon_new, vlon)
                            if any(allpoints['Longitude']<0.0):
                                idx2=np.where(allpoints['Longitude']<0.0)
                                allpoints['Longitude'][idx2] = allpoints['Longitude'][idx2]+360.0
                            
                            # because mixed points, must make sure nearest points are NONE of lat_new/lon_new
                            AllNearest=False
                            Kpts = 2
                            idx_void = None
                            while not AllNearest:
                                #print('kpts: ', Kpts, i)
                                dist,idx=nearest_using_kdtree(allpoints, latname='Latitude', lonname='Longitude',kpt=Kpts)
                                if idx_void is None:
                                    idx_new = idx[0:len_new]-len_new # for vlat/vlon real index, must adjust 'idx' in mixed points 
                                else:
                                    # fill in the negative idx only
                                    idx_new[idx_void] = (idx[...,0:len_new]-len_new)[idx_void]
                                idx_void= np.where(idx_new<0)
                                if (len(idx_void[0])>0):
                                    Kpts = Kpts+1
                                    # this will re-do the nearest search but with 1 more Kpts
                                else:
                                    AllNearest=True # exit 'while' block
                            # by now, idx_new should be all non-negative
                            temp = vdata[i,idx_new.flatten()]
                            temp[temp<0.0001]=0.0
                        
                        vdata_new[i,] = np.copy(temp)
                        del temp

                    # if 'sum' is less than 100, the residue is actually 'water' in original data
                    vdata_sum = np.sum(vdata_new, axis=0)
                    idx=np.where(vdata_sum<100.0) # excluding 100 and above
                    if len(idx[0])>0:
                        vdata_sumpft = np.ones_like(vdata_sum)*100.0
                        vdata_sumpft[idx]=vdata_sum[idx]
                        # 'water' in original data shall be called 'lake' in ELM land units
                        vname = 'PCT_LAKE'
                        if vname not in vnames: 
                            vnames.append(vname)
                            print ('variable: ', vname)
                            vtype = src2.variables[vname].datatype
                            vdim = src2.variables[vname].dimensions
                            if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                                vdim_new=[]
                                for i in range(len(vdim)):
                                    if vdim[i]=='lsmlon':
                                        vdim_new.append('gridcell')
                                    elif vdim[i]!='lsmlat':
                                        vdim_new.append(vdim[i])
                                vdim = vdim_new
                            dst.createVariable(vname, vtype, vdim)
                            dst[vname].setncatts(src2[vname].__dict__)
                        temp=100.0-vdata_sumpft
                        temp[np.where(temp<0.0001)]=0.0
                        dst[vname][...] = temp

                    #
                    # in original arctic pft data, if non-vegetated includes glacier, urban, wetland 
                    vdata_sumpft = np.ones_like(vdata_sum)*100.0
                    if 'PCT_URBAN' in dst.variables.keys():
                        temp = dst['PCT_URBAN'][...]
                        vdata_new[0,:] = np.maximum(0.0, (vdata_new[0,:] - np.round(np.sum(temp, axis=0),3)))
                        vdata_sumpft = np.maximum(0.0, (vdata_sumpft - np.round(np.sum(temp, axis=0),3)))
                    if 'PCT_LAKE' in dst.variables.keys():
                        temp = dst['PCT_LAKE'][...]
                        vdata_new[0,:] =  np.maximum(0.0, (vdata_new[0,:]-np.round(temp,3)))
                        vdata_sumpft =  np.maximum(0.0, (vdata_sumpft-np.round(temp,3)))
                    if 'PCT_WETLAND' in dst.variables.keys():
                        temp = dst['PCT_WETLAND'][...]
                        vdata_new[0,:] = np.maximum(0.0, (vdata_new[0,:]-np.round(temp,3)))
                        vdata_sumpft = np.maximum(0.0, (vdata_sumpft-np.round(temp,3)))
                    if 'PCT_GLACIER' in dst.variables.keys():
                        temp = dst['PCT_GLACIER'][...]
                        vdata_new[0,:] = np.maximum(0.0, (vdata_new[0,:]-np.round(temp,3)) )
                        vdata_sumpft = np.maximum(0.0, (vdata_sumpft-np.round(temp,3)) )
                    # 'pct_natveg' as residue of any of above
                    vname = 'PCT_NATVEG'
                    if vname not in vnames: 
                        vnames.append(vname)
                        print ('variable: ', vname)
                        vdim = src2.variables[vname].dimensions
                        if 'gridcell' not in vdim and UNSTRUCTURED_DOMAIN:
                            vdim_new=[]
                            for i in range(len(vdim)):
                                if vdim[i]=='lsmlon':
                                    vdim_new.append('gridcell')
                                elif vdim[i]!='lsmlat':
                                    vdim_new.append(vdim[i])
                            vdim = vdim_new
                        vtype = src2.variables[vname].datatype
                        dst.createVariable(vname, vtype, vdim)
                        dst[vname].setncatts(src2[vname].__dict__)
                    dst[vname][...] = vdata_sumpft
                    
 
                    # now, 'vdata_new' is the nat_pft_veg absolute values 
                    # need to make sure all natural pfts with relative values and summed to 100%
                    vname = 'PCT_NAT_PFT'
                    vdata_new[np.where(vdata_new<vdata_min)] = vdata_min
                    vdata_new[np.where(vdata_new>vdata_max)] = vdata_max
                    
                    vdata_sum = np.sum(vdata_new, axis=0)
                    idx_zero = np.where(vdata_sum<=0.0) # ELM will check even for non-natveg land unit
                    vdata_new[0,idx_zero] = 100.0       # So, put 100.0 as non_vegetated
                    idx_nonzero = np.where(vdata_sum>0.0)
                    for i in range(vdata_new.shape[0]):
                        vdata_new[i,idx_nonzero] = vdata_new[i,idx_nonzero]/vdata_sum[idx_nonzero]*100.0
                    # round up to 3 digits to avoid roundup issue upon machine
                    vdata_new = np.round(vdata_new, 3)
                    vdata_new[np.where(vdata_new<0.0)]=0.0; vdata_new[np.where(vdata_new>100.0)]=100.0  #just in case
                    idx_pftmax = np.argmax(vdata_new,0) # put roundup error into max. pft fraction
                    vdata_new[idx_pftmax,np.arange(vdata_new.shape[1])] = \
                         vdata_new[idx_pftmax,np.arange(vdata_new.shape[1])] + (100.0 - np.sum(vdata_new,axis=0))
                    
                    dst[vname][...] = np.copy(vdata_new)
                    
                    #
                    del vdata, vdata_new, vlat, vlon, vdata_sum, vdata_sumpft,temp
            fdata_src1.close()
        #
        
#-----------------------------------------------
#------ rest of surface data for ELM
        for vname in src2.variables.keys():
            if vname not in vnames:
                vdim = src2.variables[vname].dimensions
                vtype = src2.variables[vname].datatype
                dst.createVariable(vname, vtype, vdim)
                dst[vname].setncatts(src2[vname].__dict__)
                
                if 'lsmpft' in vdim and os.path.isfile(fsrfnc_natveg_pft):
                    # this is for a few SP variables - need more work (TODO)
                    vdata_new = np.asarray(dst.variables[vname][...])
                    temp = np.asarray(src2.variables[vname])
                    lsmpft_new = min(vdata_new.shape[1],temp.shape[1])# need to be more general (TODO)
                    vdata_new[:,0:lsmpft_new,...] = temp[:,0:lsmpft_new,...]
                    dst[vname][...] = np.copy(vdata_new)
                else:
                    dst[vname][...] = np.asarray(src2.variables[vname])
                
                    
                vnames.append(vname)
            #
#-----------------------------------------------
#------ final checking of PCT of land units
        sum_wtlunit = dst['PCT_NATVEG'][...] + \
                      np.sum(dst['PCT_URBAN'][...],axis=0) + \
                      dst['PCT_LAKE'][...] + \
                      dst['PCT_GLACIER'][...] + \
                      dst['PCT_WETLAND'][...]
        idx_adj = np.where(sum_wtlunit<=0)
        if len(idx_adj[0])>0: # just in case
            dst['PCT_NATVEG'][...]  = 100.0  # 100% as non-vegetated natural land
            dst['PCT_NAT_VEG'][...] = 0.0
            dst['PCT_NAT_VEG'][0,...] = 100.0
            dst['PCT_URBAN'][...]   = 0.0
            dst['PCT_LAKE'][...]    = 0.0
            dst['PCT_GLACIER'][...] = 0.0
            dst['PCT_WETLAND'][...] = 0.0
            
        idx_adj = np.where(sum_wtlunit!=100.0)
        if len(idx_adj[0])>0:
            dst['PCT_NATVEG'][...]  = dst['PCT_NATVEG'][...]/sum_wtlunit*100.0
            dst['PCT_URBAN'][...]   = dst['PCT_URBAN'][...]/sum_wtlunit*100.0
            dst['PCT_LAKE'][...]    = dst['PCT_LAKE'][...]/sum_wtlunit*100.0
            dst['PCT_GLACIER'][...] = dst['PCT_GLACIER'][...]/sum_wtlunit*100.0
            dst['PCT_WETLAND'][...] = dst['PCT_WETLAND'][...]/sum_wtlunit*100.0

        sum_wtlunit = dst['PCT_NATVEG'][...] + \
                      np.sum(dst['PCT_URBAN'][...],axis=0) + \
                      dst['PCT_LAKE'][...] + \
                      dst['PCT_GLACIER'][...] + \
                      dst['PCT_WETLAND'][...]
        idx_adj = np.where(np.abs(sum_wtlunit-100.0)>1.0E-8)
        if len(idx_adj[0])>0: print('LUNITs wt summed not equal to 100% in ',len(idx_adj[0]), 'grids!')


        #
    #with both 'src2' (original) and 'dst' ncfiles open
    
    if HAS_MPI4PY and myrank==0:
        # merge 'dst' ncfiles
        ncfileall = sorted(glob.glob("*%s" % fsurf_new))
        for i in range(len(ncfileall)):
            if i==0: 
                os.system('cp '+ncfileall[i]+' temp0.nc')
            else:
                ncmod.mergefilesby1dim('temp0.nc', ncfileall[i], fsurf_new, 'gridcell')
                os.system('cp '+fsurf_new+' temp0.nc')
                
    #

# end if true for doing merging all new high-res surfdata 



