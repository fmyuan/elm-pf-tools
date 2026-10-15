#!/usr/bin/env python
import numpy as np

try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False

#
# Arctic PFT classes in B. Sulman et al (2021) paper: 12 arctic PFTs + 2 additional tree PFTs
    # when aggregating into 17 ELM default PFTs, it shall be like following
        # lichen as not_vegetated (0), moss/forb/graminoids as c3 arctic grass (12),
        # evergreen shrub(9), deci. boreal_shrub(11),
        # evergreen boreal tree(2), deci boreal tree (3)   
arctic_pfts={'pftname':[
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
            'pftnum': [0,1,2,3,4,5,6,7,8,9,10,11,12,13],
            'pftnum_default17':[0,0,12,2,3,9,9,11,11,11,11,12,12,12]
            };

#---
            
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
                                newdata.assign_coords(band=[2])
                                ], dim='band')
    
        newdata.values = pct_glacier*100.0
        xr = xarray.concat([xr, 
                                newdata.assign_coords(band=[3])
                                ], dim='band')

        newdata.values = pct_water*100.0
        xr = xarray.concat([xr, 
                                newdata.assign_coords(band=[4])
                                ], dim='band')

        newdata.values = pct_natveg*100.0
        xr = xarray.concat([xr, 
                                newdata.assign_coords(band=[5])
                                ], dim='band')

        newdata.values = pct_barren*100.0
        xr = xarray.concat([xr, 
                                newdata.assign_coords(band=[6])
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
                                    rasterunit = 'PCT',
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
    
    # raster could in unit of Percentage, PCT, or Fraction of land coverage
    # BUT, output must be in PCT.
    UNIT2PCT=1.0        # by default raster in unit of percentage, otherwise in fraction
    if rasterunit!='PCT': UNIT2PCT=100.0
    
    # read tif data    
    alldata = {}
    for i in range(len(bandinfos['pftnum'])):
        ifile = inputpath+rasterfiles[i]
        iband = bandinfos['bands'][i]
        ipft = bandinfos['pftnum'][i]
        with rioxarray.open_rasterio(ifile, mode='r') as f:
            #
            if not bbox_lb_ru is None:
                if bbox_crs != f.rio.crs:
                    bounds = transform_bounds(bbox_crs, f.rio.crs, *bbox_lb_ru)
                else:
                    bounds = bbox_lb_ru
                xr = f.rio.clip_box(bounds[0],bounds[1],bounds[2],bounds[3])
            else:
                xr = f[0]
            if 'int' in str(xr.dtype):
                # image may be in integer, but not precise enough for calculations
                xr = xr.astype('float32')
                xr.data[np.where(xr.data==xr.rio.nodata)]=np.nan
                xr.rio.write_nodata(np.nan, inplace=True)
            xr_nodata = xr.rio.nodata
            xr_res = np.asarray(xr.rio.resolution())             

            #
            mdata = np.ma.masked_values(xr.data, xr_nodata)*UNIT2PCT

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
            elif ipft in arctic_pfts['pftnum'] and ipft!=0:
                alldata[ipft] = mdata # only 1 banded data
            else:
                alldata[iband] = mdata # iband should be a string or number not pftnum SO THAT not confused
                
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
        if 'barren' not in alldata.keys():
            alldata['barren'] = 100.0 - sum_real_pft
        else:
            alldata['barren'][zerosum_indx] = 100.0
               
    # normalization here is only for so-called vegetated area, i.e., 
    #     excluding not_vegetated (i.e. bareground/barren et al).     
    if not pctvpft_normalized:
        for iv in alldata.keys():
            if iv not in arctic_pfts['pftnum'] or iv=='0': continue
            # normalized non-zero sum only needed
            alldata[iv][nonzero_indx] = 100.0*alldata[iv][nonzero_indx] \
                                            /sum_real_pft[nonzero_indx]
    
    # save all data into rioxarray. NOTE: all real vegPFT are normalized over veged-only portion, 
    # but not-vegetated PCT are saved as it is, if raster data provided.
    newdata = xr.copy()
    idx_data= np.where(newdata!=xr_nodata)                      
    # make a copy of xr, to be a template for 'band' value updating and concat below
    for i in range(1, len(arctic_pfts['pftnum'])): #excluding pft0
        ipft = arctic_pfts['pftnum'][i]
        if ipft in alldata.keys():
            newdata.values = alldata[ipft]
        else:
            newdata.values[idx_data] = 0.0  #
        
        if i==1:
            xr.values = newdata.values
            xr_band = [arctic_pfts['pftname'][i]]
        else:
            xr = xarray.concat([xr, 
                                newdata.assign_coords(band=[i])
                                ], dim='band')
            xr_band.append(arctic_pfts['pftname'][i])
            
    
    # NON-vegetation data, including LU and PFT0 (barren)
    for v in alldata.keys():
        if not v in arctic_pfts['pftnum']:
            newdata.values = alldata[v]
            i=i+1   # continue from last PFT number
            xr = xarray.concat([xr, 
                                newdata.assign_coords(band=[i])
                                ], dim='band')
            xr_band.append(v)
    
    '''
    # better to know grid area in unit of km2
    area=np.ma.zeros_like(mdata)
    area[np.ma.where(mdata)]=np.abs(np.deg2rad(xr_res[0])*np.deg2rad(xr_res[1]))
    lat2d=np.transpose(np.tile(xr.y.data, (xr.x.size,1)))
    lat2d=np.ma.masked_array(lat2d, mask=np.ma.getmask(mdata))
    area_km2=elm_km2_from_arcradians2(area, latitude=lat2d)
    newdata.values = area_km2
    xr = xarray.concat([xr, 
                        newdata.assign_coords(band=[i+1])
                        ], dim='band')
    '''
    # band index named
    xr = xr.assign_coords(band=[str(s).strip() for s in xr_band])
    
    
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
        LandUnits, other than naturally-vegetated, will be from CAVM categories, by providing 'lupft0' DataArray 
    '''
    
    
    from pytools.commons_utils.gridlocator import elmdomain_xrio
    from pytools.pan_arctic.arctic_domain import DIN_LOC_ROOT
    from pytools.pan_arctic.arctic_domain import elm_km2_from_arcradians2 
    
    from rasterio.enums import Resampling
    import xarray
    import gc
    
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
    #xrio_elmdomain_gid=xrio_elmdomain[0]  # appears this is a messy id dataset in domain.nc created in other scripts (TODO: fix it asap)
    # instead, re-do it from 'mask' data structure
    xrio_elmdomain_gid = xrio_elmdomain_mask.copy()
    gid = np.indices((xrio_elmdomain_mask.to_numpy()).flatten().shape)
    xrio_elmdomain_gid.data = gid.reshape(xrio_elmdomain_mask.shape)
 
    
                
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
        ib_gridIDed = xrio_elmdomain_gid.rio.reproject_match(
                                veggedonlypft[0],
                                resampling=Resampling.nearest)
        ib_gridID = ib_gridIDed.drop_vars(['band','spatial_ref'])
        ib_gridID.name = 'gridID'
        
        if 'water' in veggedonlypft.band:
            ib = np.where(veggedonlypft.band=='water')[0][0]
            ib_data = veggedonlypft[ib].drop_vars(['band','spatial_ref'])
            ib_data.name = str(veggedonlypft.band[ib].data)
            ib_data_df = xarray.merge([ib_data, ib_gridID]).to_dataframe()
            aggr_data = ib_data_df.groupby('gridID')[ib_data.name].mean().to_numpy()
            pct_water = aggr_data.reshape(xrio_elmdomain_mask.rio.shape)
            del ib, ib_data, ib_data_df, aggr_data
            gc.collect()
        else:
            pct_water = np.zeros(xrio_elmdomain_mask.rio.shape, dtype=float)
        
        if 'barren' in veggedonlypft.band:
            ib = np.where(veggedonlypft.band=='barren')[0][0]
            ib_data = veggedonlypft[ib].drop_vars(['band','spatial_ref'])
            ib_data.name = str(veggedonlypft.band[ib].data)
            ib_data_df = xarray.merge([ib_data, ib_gridID]).to_dataframe()
            aggr_data = ib_data_df.groupby('gridID')[ib_data.name].mean().to_numpy()
            pct_barren = aggr_data.reshape(xrio_elmdomain_mask.rio.shape)
            del ib, ib_data, ib_data_df, aggr_data
            gc.collect()
        else:
            pct_barren = np.zeros(xrio_elmdomain_mask.rio.shape, dtype=float)
        
        # the rests
        pct_glacier = np.zeros(xrio_elmdomain_mask.rio.shape, dtype=float)
        pct_natveg = 100.0 - pct_water - pct_glacier
        
    
           
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
            
            ib_gridIDed = xrio_elmdomain_gid.rio.reproject_match(
                                lupft0[0],
                                resampling=Resampling.nearest)
            ib_gridID = ib_gridIDed.drop_vars(['band','spatial_ref'])
            ib_gridID.name = 'gridID'
            
            if 'PCT_NAT_PFT0' in lupft0.band:
                ib = np.where(lupft0.band=='PCT_NAT_PFT0')[0][0]
                ib_data = lupft0[ib].drop_vars(['band','spatial_ref'])
                ib_data.name = str(lupft0.band[ib].data)
                ib_data_df = xarray.merge([ib_data, ib_gridID]).to_dataframe()
                aggr_data = ib_data_df.groupby('gridID')[ib_data.name].mean().to_numpy()
                pct_nat_pft[0] = aggr_data.reshape(xrio_elmdomain_mask.rio.shape)
                del ib, ib_data, ib_data_df, aggr_data  # since memory use could be problematic
                gc.collect() 
        
            # other LULC
            if 'PCT_GLACIER' in lupft0.band:
                ib = np.where(lupft0.band=='PCT_GLACIER')[0][0]
           
                # the following 'Resampling.average' actually did NOT exclude nodata=np.nan
                #ib_resampled = lupft0[ib].rio.reproject_match(
                #                    xrio_elmdomain_mask,
                #                    resampling=Resampling.average,
                #                    nodata=np.nan)             
                #pct_glacier = ib_resampled.to_numpy()
                ib_data = lupft0[ib].drop_vars(['band','spatial_ref'])
                ib_data.name = str(lupft0.band[ib].data)
                ib_data_df = xarray.merge([ib_data, ib_gridID]).to_dataframe()
                aggr_data = ib_data_df.groupby('gridID')[ib_data.name].mean().to_numpy()
                pct_glacier = aggr_data.reshape(xrio_elmdomain_mask.rio.shape)
                del ib, ib_data, ib_data_df, aggr_data
                gc.collect() 
                
            if 'PCT_WATER' in lupft0.band:
                ib = np.where(lupft0.band=='PCT_WATER')[0][0]
                ib_data = lupft0[ib].drop_vars(['band','spatial_ref'])
                ib_data.name = str(lupft0.band[ib].data)
                ib_data_df = xarray.merge([ib_data, ib_gridID]).to_dataframe()
                aggr_data = ib_data_df.groupby('gridID')[ib_data.name].mean().to_numpy()
                pct_water = aggr_data.reshape(xrio_elmdomain_mask.rio.shape)
                del ib, ib_data, ib_data_df, aggr_data
                gc.collect() 

            if 'PCT_NATVEG' in lupft0.band:
                ib = np.where(lupft0.band=='PCT_NATVEG')[0][0]
                ib_data = lupft0[ib].drop_vars(['band','spatial_ref'])
                ib_data.name = str(lupft0.band[ib].data)
                ib_data_df = xarray.merge([ib_data, ib_gridID]).to_dataframe()
                aggr_data = ib_data_df.groupby('gridID')[ib_data.name].mean().to_numpy()
                pct_natveg = aggr_data.reshape(xrio_elmdomain_mask.rio.shape)
                del ib, ib_data, ib_data_df, aggr_data
                gc.collect() 
            
        except ValueError:
            print(ValueError)
             
    # really veged
    if not veggedonlypft is None:
        ib_gridIDed = xrio_elmdomain_gid.rio.reproject_match(
                                veggedonlypft[0],
                                resampling=Resampling.nearest)
        ib_gridID = ib_gridIDed.drop_vars(['band','spatial_ref'])
        ib_gridID.name = 'gridID'
        
        for iv in arctic_pfts['pftnum']:
            if iv==0: continue    
                 
            # get real veged PFT fraction from  JKumar & TQ Zhang et al.
            # vegedonly_pftdata, from JKumar and T-Q Zhang etal
            v=arctic_pfts['pftname'][iv].strip()
            if  v in veggedonlypft.band:
                ib = np.where(veggedonlypft.band==v)[0][0]
                
                # the following 'Resampling.average' actually did NOT exclude nodata=np.nan
                #ib_resampled = veggedonlypft[ib].rio.reproject_match(
                #                    xrio_elmdomain_mask,
                #                    resampling=Resampling.average,
                #                    nodata=np.nan)                            
                #pct_nat_pft[iv,...] = ib_resampled.to_numpy()*100.0  # fraction --> PCT
                
                ib_data = veggedonlypft[ib].drop_vars(['band','spatial_ref'])
                ib_data.name = str(veggedonlypft.band[ib].data)
                ib_data_df = xarray.merge([ib_data, ib_gridID]).to_dataframe()
                aggr_data = ib_data_df.groupby('gridID')[ib_data.name].mean().to_numpy()
                pct_nat_pft[iv] = aggr_data.reshape(xrio_elmdomain_mask.rio.shape)
                pct_nat_pft[iv] = pct_nat_pft[iv]
                del ib, ib_data, ib_data_df, aggr_data
                gc.collect() 
    
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
    gidxy = np.where(xrio_elmdomain_mask.to_numpy()>0.0)
    ib = np.where(xrio_elmdomain.band=='xc')
    xc = (xrio_elmdomain[ib][0].to_numpy())[gidxy]
    ib = np.where(xrio_elmdomain.band=='yc')
    yc = (xrio_elmdomain[ib][0].to_numpy())[gidxy]
    #xv = xrio_elmdomain[idx]
    #yv = xrio_elmdomain[idx]
    ib = np.where(xrio_elmdomain.band=='area')
    area = (xrio_elmdomain[ib][0].to_numpy())[gidxy]
    area_km2 = elm_km2_from_arcradians2(area, latitude=yc) 
    pct_water = pct_water[gidxy]
    pct_natveg = pct_natveg[gidxy]
    pct_glacier = pct_glacier[gidxy]
    pct_nat_pft = [pct_nat_pft[iv,...][gidxy] for iv in arctic_pfts['pftnum']]
    
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
#--------------------------------------------------------------------

def topounits_lupft_fromcavm_jkumaretal(topounit_args, user_domainnc_2d,
                                    xrio_veggedonlypft=None,
                                    xrio_lupft0=None,
                                    outdata_tiff=''):
    
    from rasterio.enums import Resampling

    # STEP 1: topographic units from DEM geotiff
    import pytools.pan_arctic.mksrfdata_topounits as mksrfdata_topounits
    from pytools.pan_arctic.mksrfdata_updatevals import topovar_topounits_stats
    from pytools.pan_arctic.mksrfdata_updatevals import maxtopounits
 
    if xrio_lupft0 is None and xrio_veggedonlypft is None:
        print()
        return {}
    
    xrio_gid, xrio_tid, xrio_elmdomain = \
        mksrfdata_topounits.grid_topounits(topounit_args, 
                                            elm_domainnc=user_domainnc_2d)
    
    # STEP 2: grid-level topounits statistics
    bgid = np.where(xrio_elmdomain.band=='GridID')[0][0]
    xrio_elmdomain_gid=xrio_elmdomain[bgid]

    xmask=xrio_elmdomain[np.where(xrio_elmdomain.band=='mask')[0][0]].to_numpy()
    xmask_gid = xrio_elmdomain_gid.to_numpy()
    xmask_gid[np.where(xmask<=0.0)] = np.nan
    xrio_elmdomain_gid.data = xmask_gid # this will exclude non-landmask gridcell below SO the data is to be unstructured
    
    # STEP 3: resolve LU-PFT data to topounits-Grids
    srfdata_fromcavm_jk2 = {}
    grd_size = np.sum(~np.isnan(xmask_gid))
    zeros_srfdata = np.zeros_like(xmask_gid[~np.isnan(xmask_gid)])           # this is gridcell-wised surfdata variable to be like
    zeros_srfdata = np.zeros((maxtopounits,)+zeros_srfdata.shape)            # currently max. 12 topounits per grid
    pct_nat_pft = np.zeros(np.append(len(arctic_pfts['pftnum']), zeros_srfdata.shape), dtype=float)

    #-----------------------------------------------
    def dfstats_topounit_to_grid(xrio_banded, var):
        ib = np.where(xrio_banded.band==var)[0][0]
        xrio_topovar = xrio_banded[ib]
        xrio_topovar = xrio_topovar.rio.reproject_match(
                                xrio_tid,
                                resampling=Resampling.average,
                                nodata=np.nan)
        df_topovar   = topovar_topounits_stats(
                                    xrio_topovar, xrio_gid, xrio_tid, topovar=var, 
                                    xrio_elmdomain_gid=xrio_elmdomain_gid)
        topo_fraction  = df_topovar["Area_fraction"].drop(columns=['gridID'])
        topo_id   = [np.where(~np.isnan(topo_fraction.iloc[i].values))[0] \
                         for i in range(grd_size)]       
        df_topovar   = df_topovar['Mean_'+var].drop(columns=['gridID'])
            
        #return np.array(df_topovar.T)
        # padding NaN in ending - this is how ELM (topounits,grid) data organized
        padding_tpovar = \
            [np.pad(df_topovar.iloc[i],
                    (0, maxtopounits-topo_id[i].size), 
                    mode='constant',constant_values=np.nan) 
                    for i in range(grd_size)]
        return padding_tpovar
    #-----------------------------------------------
      
    #LU
    if not xrio_lupft0 is None:
        srfdata_fromcavm_jk2['PCT_CROP']    = zeros_srfdata                      # not yet
        srfdata_fromcavm_jk2['PCT_URBAN']   = np.zeros((3,)+zeros_srfdata.shape) # not yet
        srfdata_fromcavm_jk2['PCT_WETLAND'] = zeros_srfdata                      # not-yet
        if 'PCT_GLACIER' in xrio_lupft0.band:
            var = 'PCT_GLACIER'
            srfdata_fromcavm_jk2['PCT_GLACIER'] = dfstats_topounit_to_grid(xrio_lupft0, var)
        else:
            srfdata_fromcavm_jk2['PCT_GLACIER'] = zeros_srfdata

        if 'PCT_WATER' in xrio_lupft0.band:
            var = 'PCT_WATER'
            srfdata_fromcavm_jk2['PCT_LAKE'] = dfstats_topounit_to_grid(xrio_lupft0, var)
        else:
            srfdata_fromcavm_jk2['PCT_LAKE'] = zeros_srfdata

        if 'PCT_NATVEG' in xrio_lupft0.band:
            var = 'PCT_NATVEG'
            srfdata_fromcavm_jk2['PCT_NATVEG'] = dfstats_topounit_to_grid(xrio_lupft0, var)
        else:
            srfdata_fromcavm_jk2['PCT_NATVEG'] = 100.0 - srfdata_fromcavm_jk2['PCT_WETLAND'] \
                                                     - srfdata_fromcavm_jk2['PCT_LAKE'] \
                                                     - srfdata_fromcavm_jk2['PCT_GLACIER'] \
                                                     - srfdata_fromcavm_jk2['PCT_CROP'] \
                                                     - np.sum(srfdata_fromcavm_jk2['PCT_URBAN'], axis=0)
        
        # PCT_NAT_PFT (summed to 100% and occupys 'PCT_NATVEG')
        # barren or bare-ground is from cavm, but will be normalized by 'PCT_NATVEG' as its PFT0
        if 'PCT_NAT_PFT0' in xrio_lupft0.band:
            var = 'PCT_NAT_PFT0'
            pct_nat_pft0 = dfstats_topounit_to_grid(xrio_lupft0, var)
            pct_nat_pft0 = pct_nat_pft0/srfdata_fromcavm_jk2['PCT_NATVEG']*100.0
            if not xrio_veggedonlypft is None:
                srfdata_fromcavm_jk2['PCT_NAT_PFT0'] = pct_nat_pft0
            else:
                pct_nat_pft[0,] = pct_nat_pft0
        
    # PCT_NAT_PFT, with all real PFTs or including 'water', 'barren'
    if not xrio_veggedonlypft is None:
        for iv in arctic_pfts['pftnum']:
            if iv==0: continue    
                 
            # get real veged PFT fraction from  JKumar & TQ Zhang et al.
            # vegedonly_pftdata, from JKumar and T-Q Zhang etal
            v=arctic_pfts['pftname'][iv].strip()
            if  v in xrio_veggedonlypft.band:
                ib = np.where(xrio_veggedonlypft.band==v)[0][0]
                xrio_topovar = xrio_lupft0[ib]
                xrio_topovar = xrio_topovar.rio.reproject_match(
                                    xrio_tid,
                                    resampling=Resampling.average,
                                    nodata=np.nan)            
                df_topovar = topovar_topounits_stats(
                                        xrio_topovar, xrio_gid, xrio_tid, topovar=v, 
                                        xrio_elmdomain_gid=xrio_elmdomain_gid)
                pct_nat_pft[iv,] = df_topovar.to_numpy()
        #
        
                
        
        srfdata_fromcavm_jk2['PCT_NAT_PFT'] = pct_nat_pft
        #
    
    #
    return srfdata_fromcavm_jk2
#--------------------------------------------------------------------



def test_hires():
    ## NOTE: tiff resolution is higher than 1km or so, like 0.0083deg, 
    ## SO CAVM 1km data not good enough for anything finer than 1km res.
    ## in this case, JK et al's landcov rasters include at least 2 additional datasets: barren, water. 
    
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
    ## TFS regional, 10km radius around the Station. NO height categories
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
    '''
    ## TFS regional
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
                        'pftnum': [0,-1,1,2,5,5,6,7,8,9,11,12]}
    '''
    ## AK within CAVM extent
    inputpath='/Users/f9y/Desktop/NGEE-P4/normalized_arcticpft_alaska/'
    rasterfiles = ['normalized_bare_ground_alaska.tif', 'normalized_water_alaska.tif',
                    'normalized_lichen_alaska.tif', 
                    'normalized_bryophyte_alaska.tif', 
                    'normalized_evergreen_shrub_dwarf_alaska.tif', 
                    'normalized_evergreen_shrub_low_alaska.tif', 
                    'normalized_evergreen_shrub_tall_alaska.tif', 
                    'normalized_deciduous_shrub_dwarf_alaska.tif',
                    'normalized_deciduous_shrub_low_alaska.tif',
                    'normalized_deciduous_shrub_tall_alaska.tif',
                    'normalized_forb_alaska.tif', 
                    'normalized_graminoid_alaska.tif']
    bandinfos={'bands':["barren","water",                # those two are separated from real vegs, later barren as PFT0 while water as PCT_LAKE
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
    
    
    vegedonly_pftdata_hires=arctic_vegedonly_fromraster_jkumaretal(\
        inputpath=inputpath,
        rasterfiles=rasterfiles,
        bandinfos=bandinfos,
        # raster unit: PCT (0-100.0) or FRACtion (0-1.0)
        rasterunit = 'PCT',
        bbox_lb_ru=tuple([-150.0, 68.0, -149.0, 69.0]),
        )

    #surf_fromcavm_jk = arctic_lupft_fromcavm_jkumaretal2( \
    #    elm_domainnc='./domain.lnd.0.0025deg.1D.c250624_TFSarcticpfts.nc', \
    #    veggedonlypft=vegedonly_pftdata_hires, \
    #    lupft0=None, \
    #    outdata_tiff='')    
    
    #user_domainnc_2d='./domain.lnd.r0005_IcoswISC30E3r5.250918_upkup_watershed2d.nc'
    user_domainnc_2d='./domain.lnd.r0005_IcoswISC30E3r5.250918_imn_watershed2d.nc'
    #user_domainnc_1d='./domain.lnd.r0005_IcoswISC30E3r5.250918_upkup_watershed1d.nc'  # this NOT YET works
    
    srfdata_fromcavm_jk = arctic_lupft_fromcavm_jkumaretal2( \
                                elm_domainnc=user_domainnc_2d, \
                                veggedonlypft=vegedonly_pftdata_hires, \
                                lupft0=None, \
                                outdata_tiff='')
        
    #user_surfdata_ncf='./surfdata_0.005x0.005_simyr1850_c250910_TOP_upkup_watershed1d-newfrac.nc'
    user_surfdata_ncf='./surfdata_0.005x0.005_simyr1850_c250910_TOP_imn_watershed1d-newfrac.nc'
    srfupdate.updatevals(
                        user_surfdata_ncf, \
                        user_srf_data=srfdata_fromcavm_jk, \
                        user_srf_vars='PCT_NATVEG, PCT_GLACIER, PCT_LAKE, PCT_WETLAND, PCT_CROP, PCT_URBAN, PCT_NAT_PFT', \
                        OriginalPFTclass=True, user_pftname_num=arctic_pfts)

    #user_surfdata_ncf='./surfdata_0.005x0.005_simyr1850_c250910_TOP_upkup_watershed1d-sub12+.nc'
    user_surfdata_ncf='./surfdata_0.005x0.005_simyr1850_c250910_TOP_imn_watershed1d-sub12+.nc'
    srfupdate.updatevals(
                        user_surfdata_ncf, \
                        user_srf_data=srfdata_fromcavm_jk, \
                        user_srf_vars='PCT_NATVEG, PCT_GLACIER, PCT_LAKE, PCT_WETLAND, PCT_CROP, PCT_URBAN, PCT_NAT_PFT', \
                        OriginalPFTclass=False, user_pftname_num=arctic_pfts)
    
    
#-----------------------------------------------------------------------------------------------------------------------

#---
def test_panarctic():
    ## NOTE: tiff resolution is about 1km or so, like 0.00833deg, 
    ## SO CAVM 1km res. data will provide 3 LUs and PFT0: GLACIER, LAKE (water), NATVEG (the rest), and non-vegetated (barren) 
    ## in this case, JK et al's landcov rasters will exclude any NON-vegged ones.
   
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
    xrio_vegedonly_pftdata = arctic_vegedonly_fromraster_jkumaretal(\
                                inputpath=inputpath,
                                bandinfos=bandinfos,
                                # raster unit: PCT (0-100.0) or FRACtion (0-1.0)
                                rasterunit = 'FRAC',
                                bbox_lb_ru=tuple([-180.0, 45.0, 180.0, 89.0]),
                                )
    
    # glacier, water (lake), PFT0 (not_vegetated) from CAVM dataset. 
    # note: didn't use what from PFT works, because it appears not fully finished (e.g. glacier not there)
    xrio_lupft0 = arctic_lu_barren_fromcavm(cavm_tiff='./raster_cavm_v1_01d.tif', \
                                bbox_lb_ru=None, bbox_crs="EPSG:4326", \
                                outdata_tiff='', \
                                outdata=True)
    
    #------------------------------------------------------------------------------------    
    IFTGU=True   

    #user_domainnc_2d = 'domain.lnd.r05_RRSwISC6to18E3r5.240328_cavm2d.nc'
    #user_domainnc_1d = 'domain.lnd.r05_RRSwISC6to18E3r5.240328_cavm1d.nc'
    user_domainnc_2d = 'domain.lnd.r0125_IcoswISC30E3r5.250918_cavm2d_TVC.nc'
    #user_domainnc_1d = 'domain.lnd.r0125_IcoswISC30E3r5.250918_cavm1d_TVC.nc'
    
    #user_surfdata_ncf = './surfdata_0.05x0.05_simyr1850_c240308_TOP_cavm1d.nc'
    user_surfdata_ncf = './topounit_surfdata_0.05x0.05_simyr1850_c240308_TOP_cavm1d.nc'

    if not IFTGU:
        # spatially merging above 2 datasets for LandUnits (glacier, water, naturally-veged, and PFTs)
        # into elm_domain grids, without topounits    
    
        srfdata_fromcavm_jk = arctic_lupft_fromcavm_jkumaretal2( \
                                elm_domainnc=user_domainnc_2d, \
                                veggedonlypft=xrio_vegedonly_pftdata, \
                                lupft0=xrio_lupft0, \
                                outdata_tiff='')
        
        srfupdate.updatevals(
                        user_surfdata_ncf, \
                        user_srf_data=srfdata_fromcavm_jk, \
                        user_srf_vars='PCT_NATVEG, PCT_GLACIER, PCT_LAKE, PCT_WETLAND, PCT_CROP, PCT_URBAN, PCT_NAT_PFT', \
                        OriginalPFTclass=True, user_pftname_num=arctic_pfts)
                        #OriginalPFTclass=False)
    
    
    #------------------------------------------------------------------------------------    
    else:
        # TOPUNITs masked grid
    
        from types import SimpleNamespace

        args = SimpleNamespace( \
                                dem = 'hyd_glo_dem_15s.tif',
                                bbox = [-180.0, 45, 180, 89],
                                lonlat_res = [0.125, 0.125], 
                                group_delta_elev=100.0,
                                bbox_outfile=''
                            )   
    
        srfdata_fromcavm_jk2 = topounits_lupft_fromcavm_jkumaretal(
                        args, 
                        user_domainnc_2d=user_domainnc_2d, 
                        xrio_veggedonlypft=xrio_vegedonly_pftdata, 
                        xrio_lupft0=xrio_lupft0)
        #    
        srfupdate.updatevals(
                        user_surfdata_ncf, \
                        user_srf_data=srfdata_fromcavm_jk2, \
                        user_srf_vars='PCT_NATVEG, PCT_GLACIER, PCT_LAKE, PCT_WETLAND, PCT_CROP, PCT_URBAN, PCT_NAT_PFT', \
                        OriginalPFTclass=False, user_pftname_num=arctic_pfts)
                        #OriginalPFTclass=False)
                            
      
if __name__ == '__main__':
    test_hires()




