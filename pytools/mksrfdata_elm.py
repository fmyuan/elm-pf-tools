#!/usr/bin/env python
import os, sys, csv, time, math
from optparse import OptionParser
import numpy as np
from scipy import interpolate
from netCDF4 import Dataset
import netcdf_modules as ncmod

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
               vname in ['lsmlat', 'lsmlon','LONGXY','LATIXY']:
                
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
                # checking if urban var
                variable1=src1.variables[vname1]
                if 'density_class' not in variable1.dimensions and \
                  vname not in ['lsmlat', 'lsmlon','LONGXY','LATIXY']:
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
    
def mksrfdata_soildtb(fsurfnc_all, fmksrfnc_soildtb, fmksrfnc_soildtb2='', fmksrfnc_soildtb_wt='', redo_grid=False):
    print('#--------------------------------------------------#')
    print("Creating surface data  - soil thickness: 'depth to bedrock' ...")
    fsurf_soildtb ='./surfdata_soildtb.nc'
    
    #--------------------------------------
    # converting from partial global to global
    dnames_elm=['lsmlat', 'lsmlon']
    dnames_dtb=['lat','lon']
    if (redo_grid or os.path.isfile(fsurf_soildtb)==False):
        if fmksrfnc_soildtb2!='' and fmksrfnc_soildtb_wt!='':
            # if provided the second DTB datasets and weight for first-sets of data
            # NOTE: this is for original data sets, which includes 2 sets, one for upland_hillslope and another for Upland_valley&lowland
            f2=Dataset(fmksrfnc_soildtb2,'r')
            f3=Dataset(fmksrfnc_soildtb_wt,'r')
            dtb2 = np.asarray(f2.variables['Band1']).astype('f4')
            dtb2[np.where(dtb2==f2.variables['Band1']._FillValue)]=np.NaN
            dtb_wt = np.asarray(f3.variables['Band1']).astype('f4')
            dtb_wt[np.where(dtb_wt==f2.variables['Band1']._FillValue)]=np.NaN
            f2.close()
            f3.close()
        
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
            laty[:] = np.copy(src1.variables['lat'][...])
            
            vname = 'lsmlon'
            vdim = ('lsmlon')
            vtype = src1.variables['lon'].datatype
            lonx=dst.createVariable(vname, vtype, vdim)
            lonx.units = 'degrees_east'
            lonx.standard_name = 'longitude'
            lonx[:] = np.copy(src1.variables['lon'][...])
            
            vname = 'LATIXY'
            vdim = ('lsmlat','lsmlon')
            vtype = src1.variables['lat'].datatype
            lat2d=dst.createVariable(vname, vtype, vdim, fill_value=-999.)
            lat2d.units = 'degrees_north'
            lat2d.standard_name = 'latitude'
            vals = dst.variables[vname][...]
            vals = np.moveaxis(vals,0,1)
            vals[:,...] = np.copy(src1.variables['lat'][...])
            lat2d[:,:] = np.moveaxis(vals,1,0)
            
            vname = 'LONGXY'
            vdim = ('lsmlat','lsmlon')
            vtype = src1.variables['lon'].datatype
            lon2d=dst.createVariable(vname, vtype, vdim, fill_value=-999.)
            lon2d.units = 'degrees_east'
            lon2d.standard_name = 'longitude'
            vals = dst.variables[vname][...]
            vals[:,...] = np.copy(src1.variables['lon'][...])
            lon2d[:,:] = np.copy(vals)
            
            
            vname = 'aveDTB'
            vdim = ('lsmlat','lsmlon')
            var2d=dst.createVariable(vname, 'f4', vdim, fill_value=-999.)
            var2d.units = 'meters below surface'
            var2d.standard_name = 'aveDTB'
            var2d.long_name = 'mean soil depth to bedrock'
            dtb1 = np.asarray(src1.variables['Band1']).astype('f4')
            src1_fillvalue = src1.variables['Band1']._FillValue
            void1 = np.where(dtb1==src1_fillvalue)
            if fmksrfnc_soildtb2!='' and fmksrfnc_soildtb_wt!='':
                varvals = dtb1*dtb_wt + dtb2*(1.0-dtb_wt)
            # fill in void cell with 'dtb2'
            varvals[void1]=dtb2[void1]
            # still void cells?
            varvals[np.where(varvals==src1_fillvalue)]=dst.variables[vname]._FillValue
            varvals[np.where(varvals==np.NaN)]=dst.variables[vname]._FillValue
            varvals[np.where(varvals==255)]=dst.variables[vname]._FillValue
            var2d[:,:] = varvals
 

#--------------------------------------------------------------------




parser = OptionParser()

parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../../../../ccsm_inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--point_list", dest="point_list", default='', \
                  help = 'File containing list of points to run (unstructured)')
parser.add_option("--point_area_kmxkm", dest="point_area_km2", default=None, \
                  help = 'user-specific area in km2 of each point in point list (unstructured')
parser.add_option("--point_area_degxdeg", dest="point_area_deg2", default=None, \
                  help = 'user-specific area in degreeXdegree of each point in point list (unstructured')
parser.add_option("--keep_duplicates", dest="keep_duplicates", default=False, \
                  help = 'Keep duplicate points', action='store_true')
parser.add_option("--usersurfnc", dest="usersurfnc", default="none", \
                  help = 'User-provided surface data nc file, with one or more variable(s) as defined')
parser.add_option("--usersurfvar", dest="usersurfvar", default="none", \
                  help = 'variable name(s) in User-provided surface data nc file, separated by ","')
parser.add_option("--nco_path", dest="nco_path", default="", \
                     help = 'NCO bin PATH, default "" ')
(options, args) = parser.parse_args()


ccsm_input = os.path.abspath(options.ccsm_input)

#------------------- get site information ----------------------------------
domainfile_orig = ccsm_input+'/share/domains/domain.clm/domain.lnd.360x720_cruncep.c20190221.nc'
surffile_orig = ccsm_input+'/lnd/clm2/surfdata_map/surfdata_360x720cru_simyr1850_c180216.nc'
#pftdyn_orig = ccsm_input+'/lnd/clm2/surfdata_map/landuse.timeseries_360x720cru_hist_simyr1850-2015_c180220.nc'


#get grid cells
longxy = np.asarray(Dataset(surffile_orig).variables['LONGXY'])
latixy = np.asarray(Dataset(surffile_orig).variables['LATIXY'])

#
# urban data --> surfdata ELM standard format
fmksrfnc_urban = ccsm_input+'/lnd/clm2/surfdata_map/high_res/mksrf_urban_0.05x0.05_simyr2000.c220127.nc'
#mksrfdata_urban(surffile_orig, fmksrfnc_urban)                   # from raw data --> grided --> surfdata
#mksrfdata_urban(surffile_orig, fmksrfnc_urban, redo_grid=False)  # from grided --> surfdata

# soil thickness data (30 secs resolution) --> surfdata ELM standard format
#fmksrfnc_soildtb = ccsm_input+'/lnd/clm2/surfdata_map/high_res/average_soil_and_sedimentary-deposit_thickness.nc'  # this datasets NOT really averaged one
fmksrfnc_soildtb = ccsm_input+'/lnd/clm2/surfdata_map/high_res/test1.nc'#upland_hill-slope_soil_thickness.nc'
fmksrfnc_soildtb2 = ccsm_input+'/lnd/clm2/surfdata_map/high_res/test2.nc'#upland_valley-bottom_and_lowland_sedimentary_deposit_thickness.nc'
fmksrfnc_soildtb_wt = ccsm_input+'/lnd/clm2/surfdata_map/high_res/test3.nc'#hill-slope_valley-bottom.nc'
mksrfdata_soildtb(surffile_orig, fmksrfnc_soildtb, \
                  fmksrfnc_soildtb2=fmksrfnc_soildtb2, fmksrfnc_soildtb_wt=fmksrfnc_soildtb_wt,redo_grid=True)              # from raw data --> grided --> surfdata
#mksrfdata_urban(surffile_orig, fmksrfnc_soildtb, redo_grid=False)  # from grided --> surfdata

# land cover type included in soil thickness data (30 secs resolution) --> surfdata ELM standard format
# including: ocean(0), upland(1), lowland(2), lake(3), perennial ice (4), assuming single coverage per grid
# SO, here we aggregate them into 0.05-deg resolution, 
fmksrftiff_landcover = ccsm_input+'/lnd/clm2/surfdata_map/high_res/mksrf_soilthk_land_cover_mask.nc'


