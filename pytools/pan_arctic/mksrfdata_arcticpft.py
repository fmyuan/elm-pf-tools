#!/usr/bin/env python
import os
import numpy as np
from netCDF4 import Dataset
from cmath import nan, inf
from pyproj import Transformer
from pyproj import CRS

from pytools.commons_utils.gridlocator import grids_nearest_or_within

try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False

#---
def arctic_veged_fromraster_jkumaretal(inputpath='./', \
                                    rasterfiles=['./graminoid_toolik.tif'], 
                                    normallized=False):
    #
    from pytools.commons_utils.Modules_netcdf import geotiff2nc
    
    # the following data are not normalized over land area, and only really-vegetated
    rasterfiles = ['lichen_toolik.tif', 
                   'bryophyte_toolik.tif', 
                   'forb_toolik.tif', 
                   'graminoid_toolik.tif', 
                   'evergreen_shrub_toolik.tif', 
                   'deciduous_shrub_toolik.tif']

    # the following are, for data above (exactly matching up with), standared PFT names and order in 14 arcticpft physiology parameters
    # a full list is: (0)not_vegetated, (1)arctic_lichen, (2)arctic_bryophyte,
    #                 (3)arctic_needleleaf_tree, (4)arctic_broadleaf_tree,
    #                 (5)arctic_evergreen_shrub_dwarf, (6)arctic_evergreen_shrub_tall,
    #                 (7)arctic_deciduous_shrub_dwarf, (8)arctic_deciduous_srhub_low, (9)arctic_deciduous_shrub_tall, (10)arctic_deciduous_shrub_alder,
    #                 (11)arctic_forb, (12)arctic_dry_graminoid,(13)arctic_wet_graminoid
    bandinfos={'bands':["arctic_lichen",
                    "arctic_bryophyte",
                    "arctic_forb",
                    "arctic_graminoid",
                    "arctic_evergreen_shrub",
                    "arctic_deciduous_shrub"],
                'pftnum': [1,2,11,12,5,8]
               };
    
    alldata = {}
    for i in range(len(bandinfos['pftnum'])):
        file = inputpath+rasterfiles[i]
        xx, yy, crs_res, crs_wkt, data = geotiff2nc(file, outdata=True)
        alldata[bandinfos['pftnum'][i]] = data[0] # only 1 banded data
    
    if not normallized:
        for i in range(len(bandinfos['pftnum'])):
            if i==0:
                sumall = alldata[bandinfos['pftnum'][i]]
            else:
                sumall = sumall+alldata[bandinfos['pftnum'][i]]
        zerosum_indx = np.where(((sumall==0.0) & (~sumall.mask)))
        nonzero_indx = np.where(((sumall>0.0) & (~sumall.mask)))
        for i in bandinfos['pftnum']: 
            alldata[i][nonzero_indx] = alldata[i][nonzero_indx] \
                                        /sumall[nonzero_indx]
            alldata[i][zerosum_indx] = 0.0  # those probably is fraction of water bodies
    
    alldata['bandinfos'] = bandinfos
    alldata['xx'] = xx
    alldata['yy'] = yy
    alldata['crs_res'] = crs_res
    alldata['crs_wkt'] = crs_wkt
    
    return alldata
    

# 


def arctic_landunit_natpft_fromcavm_jkumaretal(cavm_domainnc='./domain.lnd.pan-arctic_CAVM.1km.1d.c241018.nc', \
                                            vegedonly_pftdata={}, \
                                            outnc_surf='', \
                                            outdata=False):
    #
    # CAVM provides freshwater, barens fractions
    cavm_alldata = Dataset(cavm_domainnc,'r')
    landcov = cavm_alldata['grid_code'][...]
    xc = cavm_alldata['xc'][...]
    yc = cavm_alldata['yc'][...]
    xv = cavm_alldata['xv'][...]
    yv = cavm_alldata['yv'][...]
    area_km = cavm_alldata['area_LAEA'][...]  # appears I got something wrong in orginal data
    area_km[...] = 1.0
    
    cavm_alldata.close()

    # fresh-water (open water), as PCT_LAKE in surfdata
    pct_water = np.zeros_like(landcov)
    pct_water[np.where(landcov==91.0)] = 1.0

    # glaciers, as PCT_GLACIER in surfdata
    pct_glacier = np.zeros_like(landcov)
    pct_glacier[np.where(landcov==93.0)] = 1.0
    
    # partial barren, assumed as category 0 in natural PFTs, i.e. not_vegetated
    # (NOTE: this will have to merged into veged only, and re-normalled to 100%)    
    pct_barren = np.zeros_like(landcov)
    #codes are: 1-B1, 2-B2a, 3-B3, 4-B4, 5-B2b, 
    # assuming 25% veged for B1, and 50% for other barren, if non-zeros existed in veged data
    pct_barren[np.where(landcov==1)] = 0.25
    pct_barren[np.where(landcov==2)] = pct_barren[np.where(landcov==2)]+0.50
    pct_barren[np.where(landcov==3)] = pct_barren[np.where(landcov==3)]+0.50
    pct_barren[np.where(landcov==4)] = pct_barren[np.where(landcov==4)]+0.50
    pct_barren[np.where(landcov==5)] = pct_barren[np.where(landcov==5)]+0.50

    pct_vegedonly = np.ones_like(landcov)
    pct_vegedonly = pct_vegedonly - pct_water - pct_glacier - pct_barren

    # pct_barren + pct_vegedonly, assumed to be PCT_NATVEG in surfdata
    pct_natveg = pct_barren + pct_vegedonly
    pct_barren = np.zeros_like(pct_natveg)
    pct_barren[np.where(pct_natveg>0.0)] = pct_barren[np.where(pct_natveg>0.0)] \
                                            /pct_natveg[np.where(pct_natveg>0.0)]
    pct_vegedonly[np.where(pct_natveg>0.0)] = pct_vegedonly[np.where(pct_natveg>0.0)] \
                                            /pct_natveg[np.where(pct_natveg>0.0)]

    # grids containing veged_natpft data (ormalize PFT fractions) from JKumar & TQ Zhang et al.
    # note that the resolution is about 20m, but in regular lat/lon mesh
    masked_pts = {}
    if 'xc' in vegedonly_pftdata.keys() and 'yc' in vegedonly_pftdata.keys():
        # fully paired centroid xc/yc, unstructured or flatten 
        masked_pts['xc'] = vegedonly_pftdata['xc']
        masked_pts['yc'] = vegedonly_pftdata['yc']
    elif 'xx' in vegedonly_pftdata.keys() and 'yy' in vegedonly_pftdata.keys():
        # xx/yy meshed, need to flatten for easier point search
        x_axis=vegedonly_pftdata['xx']
        y_axis=vegedonly_pftdata['yy']
        yy,xx=np.meshgrid(y_axis, x_axis, indexing='ij')
        masked_data = vegedonly_pftdata[vegedonly_pftdata['bandinfos']['pftnum'][0]]
        masked_pts['xc']=xx[~masked_data.mask]
        masked_pts['yc']=yy[~masked_data.mask]
    #        
    cavm_grids={}
    cavm_grids['xc'] = xc
    cavm_grids['yc'] = yc
    cavm_grids['xv'] = xv
    cavm_grids['yv'] = yv
    cavm_grids_sub, grids_uid, grids_maskptsid = grids_nearest_or_within( \
                                            src_grids=cavm_grids, masked_pts=masked_pts, \
                                            remask_approach = 'within', \
                                            keep_duplicated=False)
    # only need trunked data
    xc = xc.flatten()[grids_uid]
    yc = yc.flatten()[grids_uid]
    xv = xv.flatten()[grids_uid]
    yv = yv.flatten()[grids_uid]
    area_km = area_km.flatten()[grids_uid]
    pct_water = pct_water.flatten()[grids_uid]
    pct_barren = pct_barren.flatten()[grids_uid]
    pct_natveg = pct_natveg.flatten()[grids_uid]
    pct_glacier = pct_glacier.flatten()[grids_uid]
    
    pct_vegedonly = pct_vegedonly.flatten()[grids_uid]
    indx_vegedonly = np.where(pct_vegedonly>0.0)    
    # the following are, for data above (exactly matching up with), standared PFT names and order in 14 arcticpft physiology parameters
    # a full list is: (0)not_vegetated, (1)arctic_lichen, (2)arctic_bryophyte,
    #                 (3)arctic_needleleaf_tree, (4)arctic_broadleaf_tree,
    #                 (5)arctic_evergreen_shrub_dwarf, (6)arctic_evergreen_shrub_tall,
    #                 (7)arctic_deciduous_shrub_dwarf, (8)arctic_deciduous_srhub_low, (9)arctic_deciduous_shrub_tall, (10)arctic_deciduous_shrub_alder,
    #                 (11)arctic_forb, (12)arctic_dry_graminoid,(13)arctic_wet_graminoid

    pct_nat_pft = np.zeros(np.append(14, pct_natveg.shape), dtype=float)
    pct_nat_pft[0,...] = pct_barren
    for iv in range(1,14):
        pct_iv = np.zeros_like(pct_vegedonly)
        
        # get real veged PFT fraction from  JKumar & TQ Zhang et al.
        # vegedonly_pftdata, from JKumar and T-Q Zhang etal
        bandinfos = vegedonly_pftdata['bandinfos']
        if iv in bandinfos['pftnum']:
            newpftdata = vegedonly_pftdata[iv]
            newpftdata = newpftdata[~newpftdata.mask]
            for i in range(len(grids_uid)):
                igrid = grids_uid[i]
                print(i, iv) # for checking
                pct_iv[i] = np.nanmean(newpftdata[grids_maskptsid[igrid]])
            #
            # after all grids_uid done            
            pct_iv[indx_vegedonly] = pct_iv[indx_vegedonly]/pct_vegedonly[indx_vegedonly]
        
        #
        pct_nat_pft[iv,...] = pct_iv
    

    if outnc_surf!='':
        print('write to land unit and natpft to surfdata: ', outnc_surf)

    if outdata:
        surfdata = {}
        surfdata['natpft_info'] = bandinfos
        surfdata['LONGXY'] = xc
        surfdata['LATIXY'] = yc
        surfdata['AREA']   = area_km
        surfdata['PCT_GLACIER'] = pct_glacier
        surfdata['PCT_LAKE']    = pct_water
        surfdata['PCT_NATVEG']  = pct_natveg
        surfdata['PCT_NAT_PFT'] = pct_nat_pft
        surfdata['PCT_CROP'] = np.zeros_like(pct_natveg)
        surfdata['PCT_URBAN'] = np.zeros_like(pct_natveg)
        surfdata['PCT_WETLAND'] = np.zeros_like(pct_natveg)
        
        return surfdata
    

# 
def mksrfdata_updatevals(fsurfnc_in, user_srf_data={}, user_srfnc_file='', user_srf_vars='', OriginPFTclass=True):
    
    print('#--------------------------------------------------#')
    print("Replacing values in surface data by merging user-provided dataset")
    fsurfnc_out ='./'+fsurfnc_in.split('/')[-1]+'-merged'
    
    
    #--------------------------------------
    #
    # PFT classes in B. Sulman et al (2021) paper: 12 arctic PFTs + 2 additional tree PFTs
    
    user_pfts={'pftname':[
                    "non_vegetated",
                    "arctic_lichen",
                    "arctic_bryophyte",
                    "arctic_needleleaf_tree",
                    "arctic_broadleaf_tree",
                    "arctic_evergreen_shrub",
                    "arctic_evergreen_tall_shrub",
                    "arctic_deciduous_dwarf_shrub",
                    "arctic_deciduous_low_shrub",
                    "arctic_low_to_tall_willowbirch_shrub",
                    "arctic_low_to_tall_alder_shrub",
                    "arctic_forb",
                    "arctic_dry_graminoid",
                    "arctic_wet_graminoid"
                    ],
                'pftnum': [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
               };
    
    
    if OriginPFTclass:
        # lichen as not_vegetated (0), moss/forb/graminoids as c3 arctic grass (12),
        # evergreen shrub(9), deci. boreal_shrub(11),
        # evergreen boreal tree(2), deci boreal tree (3)
        user_pfts['pftnum'] = [0,0,12,2,3,9,9,11,11,11,11,12,12,12]
        natpft = np.asarray(range(17))
    else:
        natpft = np.asarray(range(max(user_pfts['pftnum'])+1)) # this is the real arcticpft order number 
 
    #--------------------------------------
    UNSTRUCTURED = False    
    if not user_srfnc_file !='':
        print('read data from: ', user_srfnc_file)
        f=Dataset(user_srfnc_file)
        if 'gridcell' in f.dimensions.items(): UNSTRUCTURED = True
    elif not len(user_srf_data)<=0:
        if len(np.squeeze(user_srf_data['LATIXY']).shape())==1:
            UNSTRUCTURED = True
            
        if user_srf_vars=='':
            user_vname = user_srf_vars.keys()
        else:
            user_vname = user_srf_vars.split(',')
        user_srf = user_srf_data
        
    #--------------------------------------
    #                    
    # write into nc file
    with Dataset(fsurfnc_in,'r') as src, Dataset(fsurfnc_out, "w") as dst:
            
        # new surfdata dimensions
        for dname, dimension in src.dimensions.items():
            if dname == 'natpft':
                len_dimension = len(natpft)            # dim length from new data
            elif dname == 'gridcell':
                len_dimension = user_srf['LATIXY'].flatten().size()
            elif dname in ['lsmlat','lat']:
                len_dimension = user_srf['LATIXY'].shape()[0]
            elif dname in ['lsmlon', 'lon']:
                len_dimension = user_srf['LONGXY'].shape()[1]                
            else:
                len_dimension = len(dimension)
            dst.createDimension(dname, len_dimension if not dimension.isunlimited() else None)
            #
            
        # create variables and write to dst
        for vname, variable in src.variables.items():

            if UNSTRUCTURED:
                vdim = variable.dimensions
                if 'gridcell' not in vdim and \
                    ('lsmlat' in vdim and 'lsmlon' in vdim):
                    vdim = vdim.replace('lsmlon', 'gridcell')
                    vdim = vdim.remove('lsmlat')
            else:
                vdim = variable.dimensions
    
            # create variables, but will update its values later 
            # NOTE: here the variable value size updated due to newly-created dimensions above
            dst.createVariable(vname, variable.datatype, vdim)
            # copy variable attributes all at once via dictionary after created
            dst[vname].setncatts(src[vname].__dict__)
                  
            # values
            src_vals = src[vname][...]

            # dimension length may change, so need to 
            if vname in user_vname:
                varvals = src[vname][...]
                #
            else:
                varvals = src[vname][...]
                #
            
            #                                
            dst[vname][...] = varvals
                    
        # end of variable-loop        
                
        
        print('user surfdata merged and nc file created and written successfully!')
        
    #
#

#--------------------------------------------------------------------
def main():

    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """  
    input_path  = './'
    output_path = './'

    # read-in geotiff data from Zhang and Jitu etal
    vegedonly_pftdata_toolik=arctic_veged_fromraster_jkumaretal(\
        inputpath='/Users/f9y/Desktop/NGEE-P4_sites/ToolikFieldStation_TFS/toolik_clip_2024_11_14/')

    surf_fromcavm_jk = arctic_landunit_natpft_fromcavm_jkumaretal( \
        cavm_domainnc='../domain.lnd.pan-arctic_CAVM.1km.1d.c241018.nc', \
        vegedonly_pftdata=vegedonly_pftdata_toolik, \
        outdata=True)
    
    mksrfdata_updatevals('./surfdata_0.5x0.5_simyr1850_c240308_TOP_cavm1d_TFSarcticpfts.nc', \
                         user_srf_data=surf_fromcavm_jk, user_srf_vars='', OriginPFTclass=False)
    
      
if __name__ == '__main__':
    main()




