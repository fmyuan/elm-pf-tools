#!/usr/bin/env python
import numpy as np

try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False

arctic_pfts={'pftname':[
                    "not_vegetated                           ",
                    "arctic_lichen                           ",
                    "arctic_bryophyte                        ",
                    "arctic_needleleaf_tree                  ",
                    "arctic_broadleaf_tree                   ",
                    "arctic_evergreen_shrub_dwarf            ",
                    "arctic_evergreen_shrub_tall             ",
                    "arctic_deciduous_shrub_dwarf            ",
                    "arctic_deciduous_srhub_low              ",
                    "arctic_deciduous_shrub_tall             ",
                    "arctic_deciduous_shrub_alder            ",
                    "arctic_forb                             ",
                    "arctic_dry_graminoid                    ",
                    "arctic_wet_graminoid                    "],
             'pftnum': [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
            };

#---
def arctic_veged_fromraster_jkumaretal(inputpath='./', \
                                    rasterfiles = ['lichen.tif', 
                                                   'bryophyte.tif', 
                                                   'evergreen_shrub.tif', 
                                                   'deciduous_shrub.tif',
                                                   'forb.tif', 
                                                   'graminoid.tif'], 
                                    bandinfos={'bands':["arctic_lichen",
                                                        "arctic_bryophyte",
                                                        "arctic_evergreen_shrub",
                                                        "arctic_deciduous_shrub",
                                                        "arctic_forb",
                                                        "arctic_graminoid"],
                                                'pftnum': [1,2,5,8,11,12]},
                                    wetness=None,
                                    pctpft_normalized=False):
    #
    import rioxarray
    
    # the following data are not normalized over land area, and only really-vegetated
    # pctpft_normalized = False

    # the following are, for data above (exactly matching up with), standared PFT names and order in 14 arcticpft physiology parameters
    # a full list is: (0)not_vegetated, (1)arctic_lichen, (2)arctic_bryophyte,
    #                 (3)arctic_needleleaf_tree, (4)arctic_broadleaf_tree,
    #                 (5)arctic_evergreen_shrub_dwarf, (6)arctic_evergreen_shrub_tall,
    #                 (7)arctic_deciduous_shrub_dwarf, (8)arctic_deciduous_srhub_low, (9)arctic_deciduous_shrub_tall, (10)arctic_deciduous_shrub_alder,
    #                 (11)arctic_forb, (12)arctic_dry_graminoid,(13)arctic_wet_graminoid
 
    # rasterfiles = ['1.tif', '2.tif']
    
    # read tif data    
    alldata = {}
    for i in range(len(bandinfos['pftnum'])):
        ifile = inputpath+rasterfiles[i]
        iband = bandinfos['bands'][i]
        ipft = bandinfos['pftnum'][i]
        with rioxarray.open_rasterio(ifile, mode='r') as f:
            # geox/y coordinates in 1-D (axis), centroid
            xx = f.x.data
            yy = f.y.data
            crs_res = f.rio.resolution()
            crs = f.rio.crs
            rdata=f[0].data 
            mdata = np.ma.array(rdata,mask=(rdata==f.rio.nodata))
            
            
            #
            if ipft in alldata.keys():
                alldata[ipft] = alldata[ipft] + mdata # add or merge
            else:
                if ipft not in arctic_pfts['pftnum']:
                    # use original 'band name' not numbered
                    alldata[iband] = mdata               
                elif ipft==0:
                    # barren or bareground, if input as PFT0, not included when normalized below                   
                    alldata['barren'] = mdata                  
                else:
                    # really vegetated
                    alldata[ipft] = mdata # only 1 banded data
            
            # wet or dry graminoid, if given info
            if iband=='arctic_graminoid' and not wetness==None:
                if ipft==12: #assumed as arctic_dry_graminoid
                    alldata[ipft] = mdata*(1.0-wetness)
                    alldata[13] = mdata*wetness
                elif ipft==13: #assumed as arctic_wet_graminoid
                    alldata[ipft] = mdata*wetness
                    alldata[12] = mdata*(1.0-wetness)
            
    
    #
    # user must know if fractions are normalized over whole land or vegetated unit only
    # because Zhang's tif files sometime include bare_ground, water, or others (like litter, unknown, ...)
    # normalization here is only for so-called vegetated area, i.e., 
    #     excluding not_vegetated (i.e. bareground/barren et al). 
    # If bareground/barren is provided, called it as 'barren' in data's dict.
    
    if not pctpft_normalized:
        #
        sumall = np.zeros_like(alldata[bandinfos['pftnum'][0]])
        for iv in alldata.keys():
            # only sum real vegs
            if iv not in arctic_pfts['pftnum'] or iv=='0': continue
            sumall = sumall+alldata[iv]
        
        zerosum_indx = np.where(((sumall<=0.0) & (~sumall.mask)))
        nonzero_indx = np.where(((sumall>0.0) & (~sumall.mask)))
        if len(zerosum_indx[0]>0):
            # those probably is fraction of water bodies or other non-vegetated and can be assumed as non-vegetated
            if '0' not in alldata.keys():
                alldata[0] = 1.0 - sumall
            else:
                alldata[0][zerosum_indx] = 1.0
               
        for iv in alldata.keys():
            if iv not in arctic_pfts['pftnum'] or iv=='0': continue
            # normalized non-zero sum only needed
            alldata[iv][nonzero_indx] = alldata[iv][nonzero_indx] \
                                        /sumall[nonzero_indx]
                
    alldata['bandinfos'] = bandinfos
    alldata['xx'] = xx
    alldata['yy'] = yy
    alldata['crs_res'] = crs_res
    alldata['crs_wkt'] = crs.wkt
    alldata['sum_real_pft'] = sumall
    
    return alldata

#
            
#

def arctic_lupft_fromcavm_jkumaretal(cavm_within_domainnc='domain.lnd.pan-arctic_CAVM.0.01d.1D.c250623.nc', \
                                            lu_vegonlypft={}, \
                                            outnc_surf='', \
                                            outdata=False):
    '''
        There are TWO steps to generate ELM surfdata of LandUnit and PFTs for its naturally-vegetated landunit
        LandUnits, other than naturally-vegetated, will be from CAVM categories: 
           (1) PCT_LAKE    -    landcov_code 91 (open water, incl. some wider river)
               PCT_GLACIER -    landcov_code 93  
    '''
    
    
    from netCDF4 import Dataset
    from pytools.commons_utils.gridlocator import grids_nearest_or_within
    from pytools.pan_arctic.arctic_domain import DIN_LOC_ROOT
    
    # the default domain nc file for CAVM extent, contained CAVM 'landcov_code'
    # 'domain.lnd.pan-arctic_CAVM.0.01d.1D.c250623.nc'
    
    if cavm_within_domainnc.startswith('/') \
        or cavm_within_domainnc.startswith('./'):
        # already in full path
        fnc = cavm_within_domainnc
    else:
        fnc = DIN_LOC_ROOT+'/share/domains/domain.clm/'+cavm_within_domainnc
    nonvegs_alldata = Dataset(fnc,'r')

    if 'xc' in nonvegs_alldata.variables.keys():
        xc = nonvegs_alldata['xc'][...]
    else:
        xc = nonvegs_alldata['LONGXY'][...]
    if 'yc' in nonvegs_alldata.variables.keys():
        yc = nonvegs_alldata['yc'][...]
    else:
        yc = nonvegs_alldata['LATIXY'][...]
            
    if 'xv' in nonvegs_alldata.variables.keys(): xv = nonvegs_alldata['xv'][...]
    if 'yv' in nonvegs_alldata.variables.keys(): yv = nonvegs_alldata['yv'][...]
    if 'area_LAEA' in nonvegs_alldata.variables.keys():
        area_km = nonvegs_alldata['area_LAEA'][...]
    elif 'area_km' in nonvegs_alldata.variables.keys():
        area_km = nonvegs_alldata['area_km'][...]
    elif 'AREA' in nonvegs_alldata.variables.keys():
        area_km = nonvegs_alldata['AREA'][...]

    # LU - ELM LandUnit of landcover
    NONVEG_FROM_CAVM = False
    if 'landcov_code' in nonvegs_alldata.variables.keys():
        
        landcov = nonvegs_alldata['landcov_code'][...]
        if landcov.any()>0.0: 
            # CAVM provides freshwater, glacier, barens fractions, if used
            NONVEG_FROM_CAVM = True
        
        if NONVEG_FROM_CAVM:    
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
        
            pct_vegedonly = np.ones_like(landcov, dtype=float)
            pct_vegedonly = pct_vegedonly - pct_water - pct_glacier - pct_barren
        
            # pct_barren + pct_vegedonly, assumed to be PCT_NATVEG in surfdata
            pct_natveg = pct_barren + pct_vegedonly
            
            pct_barren[np.where(pct_natveg>0.0)] = pct_barren[np.where(pct_natveg>0.0)] \
                                                    /pct_natveg[np.where(pct_natveg>0.0)]
            pct_vegedonly[np.where(pct_natveg>0.0)] = pct_vegedonly[np.where(pct_natveg>0.0)] \
                                                    /pct_natveg[np.where(pct_natveg>0.0)]
    if not NONVEG_FROM_CAVM:
        landcov = np.zeros_like(xc, dtype=float)-999.0
        # no data sources for non-vegetated, assuming as 0
        pct_water   = np.zeros_like(xc)
        pct_glacier = np.zeros_like(xc)
        pct_barren  = np.zeros_like(xc)                
        pct_vegedonly = np.ones_like(xc, dtype=float)
        pct_natveg  = pct_barren + pct_vegedonly
    #
    nonvegs_alldata.close()
        

    # grids containing veged_natpft data (normalize PFT fractions) from JKumar & TQ Zhang et al.
    # note that the resolution is about 20m, but in regular lat/lon mesh
    masked_pts = {}
    if 'xc' in lu_vegonlypft.keys() and 'yc' in lu_vegonlypft.keys():
        # fully paired centroid xc/yc, unstructured or flatten 
        masked_pts['xc'] = lu_vegonlypft['xc']
        masked_pts['yc'] = lu_vegonlypft['yc']
    elif 'xx' in lu_vegonlypft.keys() and 'yy' in lu_vegonlypft.keys():
        # xx/yy meshed, need to flatten for easier point search
        x_axis=lu_vegonlypft['xx']
        y_axis=lu_vegonlypft['yy']
        yy,xx=np.meshgrid(y_axis, x_axis, indexing='ij')
        masked_data = lu_vegonlypft[lu_vegonlypft['bandinfos']['pftnum'][0]]
        masked_pts['xc']=xx[~masked_data.mask]
        masked_pts['yc']=yy[~masked_data.mask]
    
    #
    # targetted grids
    #        
    cavm_grids={}
    cavm_grids['xc'] = xc
    cavm_grids['yc'] = yc
    cavm_grids['xv'] = xv
    cavm_grids['yv'] = yv
    
    # the following is not good for large amount of grids or meshes
    cavm_grids_sub, boxed_idx, grids_uid, grids_maskptsid = \
                        grids_nearest_or_within( \
                                                src_grids=cavm_grids, \
                                                masked_pts=masked_pts, \
                                                remask_approach = 'within', \
                                                keep_duplicated=False)   
    #
    #                    
    if not NONVEG_FROM_CAVM:
        for i in range(len(grids_uid)):
            igrid = grids_uid[i]    
            # water/barren included in Zhang/Jitu's data
            if 'water' in lu_vegonlypft.keys():
                newludata = lu_vegonlypft['water']
                newludata = newludata[~newludata.mask]
                #print(i) # for checking
                pct_water.ravel()[i] = np.nanmean(newludata[grids_maskptsid[igrid]])
        pct_natveg = pct_natveg - pct_water
        pct_vegedonly = pct_natveg - pct_barren
            
    # the following are, for data above (exactly matching up with), standared PFT names and order in 14 arcticpft physiology parameters
    # a full list is: (0)not_vegetated, (1)arctic_lichen, (2)arctic_bryophyte,
    #                 (3)arctic_needleleaf_tree, (4)arctic_broadleaf_tree,
    #                 (5)arctic_evergreen_shrub_dwarf, (6)arctic_evergreen_shrub_tall,
    #                 (7)arctic_deciduous_shrub_dwarf, (8)arctic_deciduous_srhub_low, (9)arctic_deciduous_shrub_tall, (10)arctic_deciduous_shrub_alder,
    #                 (11)arctic_forb, (12)arctic_dry_graminoid,(13)arctic_wet_graminoid

    pct_nat_pft = np.zeros(np.append(14, landcov.shape), dtype=float) 
    # barren/not_vegetated
    if NONVEG_FROM_CAVM:
        pct_nat_pft[0,...] = pct_barren
    else:
        pct_iv = pct_nat_pft[0]
        try:
            if '0' in lu_vegonlypft.keys():
                newpftdata = lu_vegonlypft[0]
            elif 'barren' in lu_vegonlypft.keys():
                newpftdata = lu_vegonlypft['barren']            
            newpftdata = newpftdata[~newpftdata.mask]
            for i in range(len(grids_uid)):
                igrid = grids_uid[i]
                pct_iv.ravel()[i] = np.nanmean(newpftdata[grids_maskptsid[igrid]])
            pct_vegedonly = pct_natveg - pct_iv
            pct_nat_pft[0] = pct_iv
        except:
            print('No barren data')  
    # really veged
    for iv in range(1,14):              
        # get real veged PFT fraction from  JKumar & TQ Zhang et al.
        # vegedonly_pftdata, from JKumar and T-Q Zhang etal
        bandinfos = lu_vegonlypft['bandinfos']
        pct_iv = pct_nat_pft[iv]
        if iv in bandinfos['pftnum']:
            newpftdata = lu_vegonlypft[iv]
            newpftdata = newpftdata[~newpftdata.mask]
            for i in range(len(grids_uid)):
                igrid = grids_uid[i]
                #print(i, iv) # for checking
                pct_iv.ravel()[i] = np.nanmean(newpftdata[grids_maskptsid[igrid]])
            #
            pct_nat_pft[iv,...] = pct_iv*pct_vegedonly
    
    #
    # checking if summed to 100%
    sumallveg = np.sum(pct_nat_pft[0:,...], axis=0)
    if ((sumallveg-1.0)<0.0).any():
        #still having some all-pft pixels of zero???. Put those gaps into not_vegetated
        idx = np.where((sumallveg-1.0)<0.0)
        pct_nat_pft[0][idx] = pct_nat_pft[0][idx]+(1.0 - sumallveg[idx])
    sumallveg = np.sum(pct_nat_pft[0:,...], axis=0)
    if ((sumallveg-1.0)>1.e-7).any():
        #still tiny over 100%?
        idx = np.where((sumallveg-1.0)>1.e-7)
        for iv in range(0,14):
            pct_nat_pft[iv,...][idx] = pct_nat_pft[iv,...][idx]/sumallveg[idx]
    

    
    
    
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


    if outnc_surf!='':
        print('write to land unit and natpft to surfdata: ', outnc_surf)

    if outdata:
        surfdata = {}
        surfdata['natpft_info'] = bandinfos
        surfdata['LONGXY'] = xc
        surfdata['LATIXY'] = yc
        surfdata['AREA']   = area_km
        surfdata['PCT_GLACIER'] = pct_glacier*100.0
        surfdata['PCT_LAKE']    = pct_water*100.0
        surfdata['PCT_NATVEG']  = pct_natveg*100.0
        surfdata['PCT_NAT_PFT'] = pct_nat_pft*100.0
        surfdata['PCT_CROP'] = np.zeros_like(pct_natveg)
        surfdata['PCT_URBAN'] = np.zeros((3,)+pct_natveg.shape)
        surfdata['PCT_WETLAND'] = np.zeros_like(pct_natveg)
        
        return surfdata
            
    #
#
#--------------------------------------------------------------------------------------------------------------------
def arctic_lu_barren_fromcavm(cavm_tiff='./raster_cavm_v1_01d.tif', bbox_lb_ru=None, bbox_crs="EPSG:4326", \
                                            outdata_tiff='', \
                                            outdata=True):
    '''
        There are TWO steps to generate ELM surfdata of LandUnit and PFTs for its naturally-vegetated landunit
        LandUnits, other than naturally-vegetated, will be from CAVM categories: 
           (1) PCT_LAKE    -    landcov_code 91 (open water, incl. some wider river)
               PCT_GLACIER -    landcov_code 93
           (2) PCT_NAT_PFT[0] - (Barren) landcov_code: 1-B1, 2-B2a, 3-B3, 4-B4, 5-B2b, all or portion of whose will be not_vegetaged in PFT_NAT_PFT
    '''

    from rasterio.warp import transform_bounds
    import rioxarray
    import xarray
    from pytools.pan_arctic.arctic_domain import elm_km2_from_arcradians2
    
    # rio reading dem
    with rioxarray.open_rasterio(cavm_tiff, masked=True) as src:
    
        if bbox_crs != src.rio.crs and not bbox_lb_ru is None:
            bounds = transform_bounds(bbox_crs, src.rio.crs, *bbox_lb_ru)
            xr = src.rio.clip_box(bounds[0],bounds[1],bounds[2],bounds[3])
        else:
            xr = src[0]
        
        rio_nodata = src.rio.nodata
        rio_res = np.asarray(src.rio.resolution())   
    del src
            
    # calculate a few variables from 'landcov'
    landcov = np.ma.masked_values(xr.data, rio_nodata)
    landcov = np.ma.masked_values(landcov, 99) # masked out non-CAVM grids
    landcov = np.ma.masked_values(landcov, 92) # masked out ocean grids
    
    # LU - ELM LandUnit of landcover
    NONVEG_FROM_CAVM = False
    if landcov.any()>0.0: 
        # CAVM provides freshwater, glacier, barens fractions, if used
        NONVEG_FROM_CAVM = True
        
    if NONVEG_FROM_CAVM:    
        # fresh-water (open water), as PCT_LAKE in surfdata
        pct_water = np.ma.zeros_like(landcov)        
        pct_water[np.ma.where(landcov==91)] = 1.0
        
        # glaciers, as PCT_GLACIER in surfdata
        pct_glacier = np.ma.zeros_like(landcov)
        pct_glacier[np.ma.where(landcov==93)] = 1.0
            
        # partial barren, assumed as category 0 in natural PFTs, i.e. not_vegetated
        # (NOTE: this will have to merged into veged only, and re-normalled to 100%)    
        pct_barren = np.ma.zeros_like(landcov, dtype=float)
        #codes are: 1-B1, 2-B2a, 3-B3, 4-B4, 5-B2b, 
        # assuming 25% veged for B1, and 50% for other barren, if non-zeros existed in veged data
        pct_barren[np.ma.where(landcov==1)] = 0.25
        pct_barren[np.ma.where(landcov==2)] = pct_barren[np.ma.where(landcov==2)]+0.50
        pct_barren[np.ma.where(landcov==3)] = pct_barren[np.ma.where(landcov==3)]+0.50
        pct_barren[np.ma.where(landcov==4)] = pct_barren[np.ma.where(landcov==4)]+0.50
        pct_barren[np.ma.where(landcov==5)] = pct_barren[np.ma.where(landcov==5)]+0.50
        
        pct_vegedonly = np.ma.ones_like(landcov, dtype=float)
        pct_vegedonly = pct_vegedonly - pct_water - pct_glacier - pct_barren
        
        # pct_barren + pct_vegedonly, assumed to be PCT_NATVEG in surfdata. Redo it as PCT_NAT_PFT0
        pct_natveg = pct_barren + pct_vegedonly
        pct_barren[np.ma.where(pct_natveg>0.0)] = pct_barren[np.ma.where(pct_natveg>0.0)] \
                                                /pct_natveg[np.ma.where(pct_natveg>0.0)]
        pct_vegedonly[np.ma.where(pct_natveg>0.0)] = pct_vegedonly[np.ma.where(pct_natveg>0.0)] \
                                                /pct_natveg[np.ma.where(pct_natveg>0.0)]
                                                
        # area in arc-deg2 to that in km2
        area=np.ma.zeros_like(landcov)
        area[np.ma.where(landcov)]=np.abs(np.deg2rad(rio_res[0])*np.deg2rad(rio_res[1]))
        lat2d=np.transpose(np.tile(xr.y.data, (xr.x.size,1)))
        lat2d=np.ma.masked_array(lat2d, mask=np.ma.getmask(landcov))
        area_km2=elm_km2_from_arcradians2(area, latitude=lat2d)
        
    #    
    #

    if outdata or outdata_tiff!='':
        xr.rio.write_nodata(np.nan, inplace=True) # original 'nodata' is an unint
        
        newdata = xr.copy()
        newdata.values = area_km2
        xr = xarray.concat([xr, 
                                newdata.reset_coords('band', drop=True).expand_dims(band=[2])
                                ], dim='band')
    
        newdata.values = pct_glacier*100.0
        xr = xarray.concat([xr, 
                                newdata.reset_coords('band', drop=True).expand_dims(band=[3])
                                ], dim='band')

        newdata.values = pct_water*100.0
        xr = xarray.concat([xr, 
                                newdata.reset_coords('band', drop=True).expand_dims(band=[4])
                                ], dim='band')

        newdata.values = pct_natveg*100.0
        xr = xarray.concat([xr, 
                                newdata.reset_coords('band', drop=True).expand_dims(band=[5])
                                ], dim='band')

        newdata.values = pct_barren*100.0
        xr = xarray.concat([xr, 
                                newdata.reset_coords('band', drop=True).expand_dims(band=[6])
                                ], dim='band')
        

        
        # band naming
        xr = xr.assign_coords(band=["landcov","AREA","PCT_GLACIER","PCT_WATER","PCT_NATVEG","PCT_NAT_PFT0"])
        xr = xr.transpose('band','y','x')
            
        if outdata_tiff!='':
            xr.rio.to_raster(outdata_tiff)
        
        if outdata: 
            return xr

def arctic_vegedonly_fromraster_jkumaretal(inputpath='./', \
                                    rasterfiles = ['lichen.tif', 
                                                   'bryophyte.tif', 
                                                   'evergreen_shrub.tif', 
                                                   'deciduous_shrub.tif',
                                                   'forb.tif', 
                                                   'graminoid.tif'], 
                                    bandinfos={'bands':["arctic_lichen",
                                                        "arctic_bryophyte",
                                                        "arctic_evergreen_shrub",
                                                        "arctic_deciduous_shrub",
                                                        "arctic_forb",
                                                        "arctic_graminoid"],
                                                'pftnum': [1,2,5,8,11,12]},
                                    wetness=None, pctvpft_normalized=False,
                                    bbox_lb_ru=None, bbox_crs="EPSG:4326", 
                                    outdata_tiff=''):
    #
    from rasterio.warp import transform_bounds
    import rioxarray
    import xarray
    from pytools.pan_arctic.arctic_domain import elm_km2_from_arcradians2
       
    # the following data are not normalized over land area, and only really-vegetated
    # pctpft_normalized = False

    # the following are, for data above (exactly matching up with), standared PFT names and order in 14 arcticpft physiology parameters
    # a full list is: (0)not_vegetated, (1)arctic_lichen, (2)arctic_bryophyte,
    #                 (3)arctic_needleleaf_tree, (4)arctic_broadleaf_tree,
    #                 (5)arctic_evergreen_shrub_dwarf, (6)arctic_evergreen_shrub_tall,
    #                 (7)arctic_deciduous_shrub_dwarf, (8)arctic_deciduous_srhub_low, (9)arctic_deciduous_shrub_tall, (10)arctic_deciduous_shrub_alder,
    #                 (11)arctic_forb, (12)arctic_dry_graminoid,(13)arctic_wet_graminoid
 
    # rasterfiles = ['1.tif', '2.tif']
    
    # read tif data    
    alldata = {}
    for i in range(len(bandinfos['pftnum'])):
        ifile = inputpath+rasterfiles[i]
        iband = bandinfos['bands'][i]
        ipft = bandinfos['pftnum'][i]
        with rioxarray.open_rasterio(ifile, mode='r') as f:
            #
            if bbox_crs != f.rio.crs and not bbox_lb_ru is None:
                bounds = transform_bounds(bbox_crs, f.rio.crs, *bbox_lb_ru)
                xr = f.rio.clip_box(bounds[0],bounds[1],bounds[2],bounds[3])
            else:
                xr = f[0]
            xr_nodata = xr.rio.nodata
            xr_res = np.asarray(xr.rio.resolution()) 
            #
            mdata = np.ma.masked_values(xr.data, xr_nodata)

            # wet or dry graminoid, if given info
            if iband=='arctic_graminoid' and not wetness==None:
                if ipft==12: #assumed as arctic_dry_graminoid
                    alldata[ipft] = mdata*(1.0-wetness)
                    alldata[13] = mdata*wetness
                elif ipft==13: #assumed as arctic_wet_graminoid
                    alldata[ipft] = mdata*wetness
                    alldata[12] = mdata*(1.0-wetness)
           
            if ipft in alldata.keys():
                alldata[ipft] = alldata[ipft] + mdata # add or merge
            elif ipft in arctic_pfts['pftnum']:
                alldata[ipft] = mdata # only 1 banded data
                
    #
    # user must know if fractions are normalized over whole land or vegetated area only
    # because Zhang's tif files sometime include bare_ground, water, or others (like litter, unknown, ...)
    # implying it's fraction over whole land grid.
    
    #
    sum_real_pft = np.ma.zeros_like(mdata)
    for iv in alldata.keys():
        # only sum real vegs
        if iv not in arctic_pfts['pftnum'] or iv=='0': continue
        sum_real_pft = sum_real_pft+alldata[iv]        
    zerosum_indx = np.where(((sum_real_pft<=0.0) & (~sum_real_pft.mask)))
    nonzero_indx = np.where(((sum_real_pft>0.0) & (~sum_real_pft.mask)))
    
    if len(zerosum_indx[0]>0):
        # those probably is fraction of water bodies or other non-vegetated 
        # and can be assumed as non-vegetated
        if '0' not in alldata.keys():
            alldata[0] = 1.0 - sum_real_pft
        else:
            alldata[0][zerosum_indx] = 1.0
               
    # normalization here is only for so-called vegetated area, i.e., 
    #     excluding not_vegetated (i.e. bareground/barren et al).     
    if not pctvpft_normalized:
        for iv in alldata.keys():
            if iv not in arctic_pfts['pftnum'] or iv=='0': continue
            # normalized non-zero sum only needed
            alldata[iv][nonzero_indx] = alldata[iv][nonzero_indx] \
                                        /sum_real_pft[nonzero_indx]
    
    # save all data into rioxarray. NOTE: all real vegPFT are normalized over veged-only portion, 
    # but not-vegetated PCT are saved as it is, if raster data provided.
    xr.rio.write_nodata(np.nan, inplace=True) # original 'nodata' is an unint, and 'xr' here is the last raster file read-in already
    newdata = xr.copy()                       # make a copy of xr, to be a template for 'band' value updating and concat.
    for i in range(len(arctic_pfts['pftnum'])):
        ipft = arctic_pfts['pftnum'][i]
        if ipft in alldata.keys():
            newdata.values = alldata[ipft]
        else:
            newdata.values[...] = xr_nodata
        
        if i==0:
            xr.values = newdata.values
        else:
            xr = xarray.concat([xr, 
                                newdata.reset_coords('band', drop=True).expand_dims(band=[i+1])
                                ], dim='band')
    # better to know grid area in unit of km2
    '''
    area=np.ma.zeros_like(mdata)
    area[np.ma.where(mdata)]=np.abs(np.deg2rad(xr_res[0])*np.deg2rad(xr_res[1]))
    lat2d=np.transpose(np.tile(xr.y.data, (xr.x.size,1)))
    lat2d=np.ma.masked_array(lat2d, mask=np.ma.getmask(mdata))
    area_km2=elm_km2_from_arcradians2(area, latitude=lat2d)
    newdata.values = area_km2
    xr = xarray.concat([xr, 
                        newdata.reset_coords('band', drop=True).expand_dims(band=[i+1])
                        ], dim='band')
    '''
    # band index named
    xr = xr.assign_coords(band=[s.strip() for s in arctic_pfts['pftname']])
    xr = xr.transpose('band','y','x')
    
    
    # write into provided tiff file.
    if outdata_tiff!='':
        xr.rio.to_raster(outdata_tiff)
        
    #
    return xr

#
            
#

#--------------------------------------------------------------------

def arctic_lupft_fromcavm_jkumaretal2(elm_domainnc='./domain.lnd.r0125_IcoswISC30E3r5.250918_cavm2d.nc',
                                    veggedonlypft=None,
                                    lupft0=None,
                                    outdata_tiff=''):
    '''
        There are TWO steps to generate ELM surfdata of LandUnit and PFTs for its naturally-vegetated landunit
        LandUnits, other than naturally-vegetated, will be from CAVM categories: 
           (1) PCT_LAKE    -    landcov_code 91 (open water, incl. some wide rivers)
               PCT_GLACIER -    landcov_code 93  
    '''
    
    
    from pytools.commons_utils.gridlocator import elmdomain_xrio
    from pytools.pan_arctic.arctic_domain import DIN_LOC_ROOT
    from pytools.pan_arctic.arctic_domain import elm_km2_from_arcradians2 
    
    from rasterio.enums import Resampling
    
    # the default domain nc file for CAVM extent, contained CAVM 'landcov_code'
    # 'domain.lnd.pan-arctic_CAVM.0.01d.1D.c250623.nc'
    
    if elm_domainnc.startswith('/') \
        or elm_domainnc.startswith('./'):
        # already in full path
        fnc = elm_domainnc
    else:
        fnc = DIN_LOC_ROOT+'/share/domains/domain.clm/'+elm_domainnc

    xrio_elmdomain = elmdomain_xrio(fnc)
    bmask = np.where(xrio_elmdomain.band=='mask')[0][0]
    xrio_elmdomain_mask=xrio_elmdomain[bmask]
    xrio_elmdomain_mask.rio.write_nodata(0.0, inplace=True)
    
                
    #xarray rio operations - this is much faster than geopd sjoin
    
    NONVEG_FROM_CAVM = False
    if not lupft0 is None:
        NONVEG_FROM_CAVM = True
        # the following appears not excludes nodata for multi-banded rxio 
        #lupft0_resampled = lupft0.rio.reproject_match(
        #            xrio_elmdomain_mask,
        #            resampling=Resampling.average)
        
    #
    
    if not NONVEG_FROM_CAVM:
        if 'water' in veggedonlypft.band:
            ib = np.where(veggedonlypft.band=='water')[0][0]
            ib_resampled = lupft0[ib].rio.reproject_match(
                                xrio_elmdomain_mask,
                                resampling=Resampling.average,
                                nodata=np.nan)                
            pct_water = ib_resampled.to_numpy() 
        
        if 'barren' in veggedonlypft.band:
            ib = np.where(veggedonlypft.band=='barren')[0][0]
            ib_resampled = lupft0[ib].rio.reproject_match(
                                xrio_elmdomain_mask,
                                resampling=Resampling.average,
                                nodata=np.nan)                
            pct_barren = ib_resampled.to_numpy() 
    
           
    # the following are, for data above (exactly matching up with), standared PFT names and order in 14 arcticpft physiology parameters
    # a full list is: (0)not_vegetated, (1)arctic_lichen, (2)arctic_bryophyte,
    #                 (3)arctic_needleleaf_tree, (4)arctic_broadleaf_tree,
    #                 (5)arctic_evergreen_shrub_dwarf, (6)arctic_evergreen_shrub_tall,
    #                 (7)arctic_deciduous_shrub_dwarf, (8)arctic_deciduous_srhub_low, (9)arctic_deciduous_shrub_tall, (10)arctic_deciduous_shrub_alder,
    #                 (11)arctic_forb, (12)arctic_dry_graminoid,(13)arctic_wet_graminoid

    pct_nat_pft = np.zeros(np.append(len(arctic_pfts['pftnum']), xrio_elmdomain_mask.rio.shape), dtype=float) 
    # barren/not_vegetated
    if not NONVEG_FROM_CAVM:
        pct_nat_pft[0,...] = pct_barren
    else:
        try:
            if 'PCT_NAT_PFT0' in lupft0.band:
                ib = np.where(lupft0.band=='PCT_NAT_PFT0')[0][0]
                ib_resampled = lupft0[ib].rio.reproject_match(xrio_elmdomain_mask,
                    resampling=Resampling.average, 
                    nodata=np.nan)                
                pct_nat_pft[0] = ib_resampled.to_numpy()
        except:
            print('No barren data')
        
        # other LULC
        if 'PCT_GLACIER' in lupft0.band:
            ib = np.where(lupft0.band=='PCT_GLACIER')[0][0]
            ib_resampled = lupft0[ib].rio.reproject_match(
                                xrio_elmdomain_mask,
                                resampling=Resampling.average,
                                nodata=np.nan)                
            pct_glacier = ib_resampled.to_numpy()
        if 'PCT_WATER' in lupft0.band:
            ib = np.where(lupft0.band=='PCT_WATER')[0][0]
            ib_resampled = lupft0[ib].rio.reproject_match(
                                xrio_elmdomain_mask,
                                resampling=Resampling.average,
                                nodata=np.nan)                
            pct_water = ib_resampled.to_numpy()
        if 'PCT_NATVEG' in lupft0.band:
            ib = np.where(lupft0.band=='PCT_NATVEG')[0][0]
            ib_resampled = lupft0[ib].rio.reproject_match(
                                xrio_elmdomain_mask,
                                resampling=Resampling.average,
                                nodata=np.nan)                
            pct_natveg = ib_resampled.to_numpy()
            
             
    # really veged
    if not veggedonlypft is None:
        for iv in arctic_pfts['pftnum']:
            if iv==0: continue    
                 
            # get real veged PFT fraction from  JKumar & TQ Zhang et al.
            # vegedonly_pftdata, from JKumar and T-Q Zhang etal
            v=arctic_pfts['pftname'][iv].strip()
            if  v in veggedonlypft.band:
                ib = np.where(veggedonlypft.band==v)[0][0]
                ib_resampled = veggedonlypft[ib].rio.reproject_match(
                                    xrio_elmdomain_mask,
                                    resampling=Resampling.average,
                                    nodata=np.nan)                            
                pct_nat_pft[iv,...] = ib_resampled.to_numpy()*100.0  # fraction --> PCT
    
    #
    # normalizing real pft to PCT_NATVEG-lu only
    sumallveg = 100.0 - pct_nat_pft[0,...] # pct_nat_pft[0] already normalized over 'pct_natveg' above
    for iv in arctic_pfts['pftnum']:
        if iv==0: continue
        pct_nat_pft[iv,...] = pct_nat_pft[iv,...]*(sumallveg/100.0)
    
    # re-checking 100% of all summed pfts, inc.
    sum_err = np.zeros_like(sumallveg)
    idx = np.where(~np.isnan(sumallveg))
    sum_err[idx] = (100.0 - np.nansum(pct_nat_pft, axis=0))[idx]
    if (np.abs(sum_err)>1.e-7).any(): 
        idx = np.where(np.abs(sum_err)>1.e-7)
        pct_nat_pft[0,...][idx] = pct_nat_pft[0,...][idx]+sum_err[idx]
    for iv in arctic_pfts['pftnum']:
        if iv==0: continue
        # appears there are data-missing in pan-arctic Arctic-pft datasets after aggregation
        idx = np.where( (~np.isnan(sumallveg)) & (np.isnan(pct_nat_pft[iv,...])) )
        if len(idx[0])>0: pct_nat_pft[iv,...][idx] = 0.0
        
    
    
    # only need trunked data
    idx = np.where(xrio_elmdomain_mask.to_numpy()>0.0)
    ib = np.where(xrio_elmdomain.band=='xc')
    xc = (xrio_elmdomain[ib][0].to_numpy())[idx]
    ib = np.where(xrio_elmdomain.band=='yc')
    yc = (xrio_elmdomain[ib][0].to_numpy())[idx]
    #xv = xrio_elmdomain[idx]
    #yv = xrio_elmdomain[idx]
    ib = np.where(xrio_elmdomain.band=='area')
    area = (xrio_elmdomain[ib][0].to_numpy())[idx]
    area_km2 = elm_km2_from_arcradians2(area, latitude=yc) 
    pct_water = pct_water[idx]
    pct_natveg = pct_natveg[idx]
    pct_glacier = pct_glacier[idx]
    pct_nat_pft = [pct_nat_pft[iv,...][idx] for iv in arctic_pfts['pftnum']]
    
    #
    surfdata = {}
    surfdata['LONGXY'] = xc
    surfdata['LATIXY'] = yc
    surfdata['AREA']   = area_km2
    surfdata['PCT_GLACIER'] = pct_glacier
    surfdata['PCT_LAKE']    = pct_water
    surfdata['PCT_NATVEG']  = pct_natveg
    surfdata['PCT_NAT_PFT'] = pct_nat_pft
    surfdata['PCT_CROP'] = np.zeros_like(pct_natveg)          # not yet
    surfdata['PCT_URBAN'] = np.zeros((3,)+pct_natveg.shape)   # not yet
    surfdata['PCT_WETLAND'] = np.zeros_like(pct_natveg)       # not-yet
        
    return surfdata
    
    #

def test_hires():
    
    import pytools.pan_arctic.mksrfdata_updatevals as srfupdate
    import pytools.pan_arctic.arctic_domain as elm_domain

    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """  

    elm_domain.set_e3sm_input('/Users/f9y/project_e3sm/e3sm_inputdata')
    
    # read-in geotiff data from Zhang and Jitu etal
    '''   
    inputpath='/Users/f9y/Desktop/NGEE-P4/Sites/ToolikFieldStation_TFS/toolik_clip_2024_11_14/'
    rasterfiles = ['lichen_toolik.tif', 
                    'bryophyte_toolik.tif', 
                    'evergreen_shrub_toolik.tif', 
                    'deciduous_shrub_toolik.tif',
                    'forb_toolik.tif', 
                    'graminoid_toolik.tif']
    bandinfos={'bands':["arctic_lichen",
                        "arctic_bryophyte",
                        'arctic_needleleaf_tree', 
                        'arctic_broadleaf_tree',
                        "arctic_evergreen_shrub",
                        "arctic_deciduous_shrub",
                        "arctic_forb",
                        "arctic_graminoid"],
                        'pftnum': [1,2,5,8,11,12]}


    '''
    # canopy height informed shrub PFTs
    inputpath='/Users/f9y/Desktop/NGEE-P4/Sites/ToolikFieldStation_TFS/PFTs_top_height_category/'
    rasterfiles = ['bare_ground_top_20m.tif', 'water_top_20m.tif',
                    'lichen_top_20m.tif', 
                    'bryophyte_top_20m.tif', 
                    'evergreen_shrub_dwarf_top_20m.tif', 
                    'evergreen_shrub_low_top_20m.tif', 
                    'evergreen_shrub_tall_top_20m.tif', 
                    'deciduous_shrub_dwarf_top_20m.tif',
                    'deciduous_shrub_low_top_20m.tif',
                    'deciduous_shrub_tall_top_20m.tif',
                    'forb_top_20m.tif', 
                    'graminoid_top_20m.tif']
    bandinfos={'bands':["barren","water",
                        "arctic_lichen",
                        "arctic_bryophyte",
                        "arctic_evergreen_dwarf_shrub",
                        "arctic_evergreen_dwarf_shrub",  # current no evergreen_low_shrub PFT, so merge into 'dwarf'
                        "arctic_evergreen_tall_shrub",
                        "arctic_deciduous_dwarf_shrub",
                        "arctic_deciduous_low_shrub",
                        "arctic_deciduous_tall_shrub",
                        "arctic_forb",
                        "arctic_graminoid"],
                        'pftnum': [-2,-1,1,2,5,5,6,7,8,9,11,12]}
    
    
    vegedonly_pftdata_toolik=arctic_veged_fromraster_jkumaretal(\
        inputpath=inputpath,
        rasterfiles=rasterfiles,
        bandinfos=bandinfos)

    surf_fromcavm_jk = arctic_lupft_fromcavm_jkumaretal( \
        cavm_within_domainnc='./domain.lnd.0.0025deg.1D.c250624_TFSarcticpfts.nc', \
        lu_vegonlypft=vegedonly_pftdata_toolik, \
        outdata=True)
    
    #srfupdate.updatevals('./surfdata_0.0025deg.1D_simyr1850_c240308_TOP_TFSarcticpfts-default.nc', \
                        #user_srf_data=surf_fromcavm_jk, user_srf_vars='PCT_NAT_VEG', OriginalPFTclass=True)
    srfupdate.updatevals('./surfdata_0.0025deg.1D_simyr1850_c240308_TOP_TFSarcticpfts-sub12+.nc', \
                        user_srf_data=surf_fromcavm_jk, user_srf_vars='PCT_NAT_VEG', OriginalPFTclass=False)
    
#--------------------------------------------------------------------

#---
def test():
    
    import pytools.pan_arctic.mksrfdata_updatevals as srfupdate
    import pytools.pan_arctic.arctic_domain as elm_domain

    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """  

    elm_domain.set_e3sm_input('/Users/f9y/project_e3sm/e3sm_inputdata')
    
    # read-in geotiff data from Zhang and Jitu etal
      
    inputpath='/Users/f9y/Desktop/NGEE-P4/Pan-arctic/'
    rasterfiles = [
                    # barren/water appears not fully covered the CAVM pixels
                    #'cavm_pavc_upscaling_v2_rf_v2_Bareground_fCover_k10fold_mean.tif', 
                    #'cavm_pavc_upscaling_v2_rf_v2_Water_fCover_k10fold_mean.tif',
                    'cavm_pavc_upscaling_v2_rf_v2_Lichen_fCover_k10fold_mean.tif',
                    'cavm_pavc_upscaling_v2_rf_v2_Bryophyte_fCover_k10fold_mean.tif',
                    'cavm_pavc_upscaling_v2_rf_v2_Evergreen_Tree_fCover_k10fold_mean.tif',
                    'cavm_pavc_upscaling_v2_rf_v2_Deciduous_Tree_fCover_k10fold_mean.tif',
                    'cavm_pavc_upscaling_v2_rf_v2_Evergreen_Shrub_fCover_k10fold_mean.tif',
                    'cavm_pavc_upscaling_v2_rf_v2_Deciduous_Shrub_fCover_k10fold_mean.tif',
                    'cavm_pavc_upscaling_v2_rf_v2_Forb_fCover_k10fold_mean.tif',
                    'cavm_pavc_upscaling_v2_rf_v2_Graminoid_fCover_k10fold_mean.tif'
                    ]
    bandinfos={'bands':[
                        #"barren",
                        #"water",
                        "arctic_lichen",
                        "arctic_bryophyte",
                        'arctic_needleleaf_tree', 
                        'arctic_broadleaf_tree',
                        "arctic_evergreen_shrub",
                        "arctic_deciduous_shrub",
                        "arctic_forb",
                        "arctic_graminoid"],
                        #'pftnum': [-2,-1,1,2,3,4,5,8,11,12]}
                        'pftnum': [1,2,3,4,5,8,11,12]}

    
    # vegetated only PFTs from 0.00833 deg pan-Arctic raster data
    # note: data in rioxarray, bandded with Arctic PFT numbering order, but normalized over vegetated area ONLY.
    xrio_vegedonly_pftdata_cavm = arctic_vegedonly_fromraster_jkumaretal(\
                                inputpath=inputpath,
                                rasterfiles=rasterfiles,
                                bandinfos=bandinfos)
    
    # glacier, water (lake), PFT0 (not_vegetated) from CAVM dataset. 
    # note: didn't use what from PFT works, because it appears not fully finished (e.g. glacier not there)
    xrio_lupft0 = arctic_lu_barren_fromcavm(cavm_tiff='./raster_cavm_v1_01d.tif', \
                                bbox_lb_ru=None, bbox_crs="EPSG:4326", \
                                outdata_tiff='', \
                                outdata=True)
    
    # spatially merging above 2 datasets for LandUnits (glacier, water, naturally-veged, and PFTs)     
    srfdata_fromcavm_jk = arctic_lupft_fromcavm_jkumaretal2( \
                                elm_domainnc='./domain.lnd.r05_RRSwISC6to18E3r5.240328_cavm2d.nc', \
                                #elm_domainnc='./domain.lnd.r0125_IcoswISC30E3r5.250918_cavm2d.nc', \
                                veggedonlypft=xrio_vegedonly_pftdata_cavm, \
                                lupft0=xrio_lupft0, \
                                outdata_tiff='')
        
    
    
    #------------------------------------------------------------------------------------    
    # TOPUNITs masked grid
    '''
    IFTGU=True   
    args={}
    args.dem = 'hyd_glo_dem_15s.tif'
    args.bbox = [-180.0, 50, 180, 85]
    args.lonlat_res = [0.125, 0.125] 
    args.group_delta_elev=50.0    
    
    # STEP 1: topographic units from DEM geotiff
    import pytools.pan_arctic.mksrfdata_topounits as mksrfdata_topounits 
 
    if not IFTGU:
        # if don't do topounit dividing    
        xrio, gridded_stats = mksrfdata_topounits.grid_stats_topounits(
                        args, percentiles=np.asarray([100]), verbose=0, 
                        elm_domainnc='domain_test.nc')
    else:
        xrio, gridded_stats = mksrfdata_topounits.grid_stats_topounits(args, verbose=1, 
                        elm_domainnc='domain_test.nc')
                        #elm_domainnc='domain.lnd.r0125_IcoswISC30E3r5.250918_cavm1d.nc')
    
    
    '''
    
    #    
    srfupdate.updatevals('./surfdata_0.5x0.5_simyr1850_c240308_TOP_cavm1d.nc', \
                        user_srf_data=srfdata_fromcavm_jk, \
                        user_srf_vars='PCT_NATVEG, PCT_GLACIER, PCT_LAKE, PCT_WETLAND, PCT_CROP, PCT_URBAN, PCT_NAT_PFT', \
                        OriginalPFTclass=True)
                        #OriginalPFTclass=False)
    
      
if __name__ == '__main__':
    test()




