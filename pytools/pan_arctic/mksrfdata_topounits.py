#!/usr/bin/env python3
"""
Jitendra Kumar (kumarj@ornl.gov) 

Elevation percentile based topounits:
 - Area-weighted fractions
 - Recursive grouping of consecutive bins using area-weighted mean
 - Group-level mean and std

This follows the topunit scheme from: 
Tesfa, T. K., Leung, L. R., Thornton, P. E., Brunke, 
M. A., & Duan, Z. (2024). Impacts of Topography‐Based 
Subgrid Scheme and Downscaling of Atmospheric Forcing 
on Modeling Land Surface Processes in the Conterminous US. 
Journal of Advances in Modeling Earth Systems, 16(8). 
https://doi.org/10.1029/2023ms004064

Fengming Yuan (yuanf@ornl.gov)
Land surface properties assigned into ELevation percetile based topounits:
 - slope, aspect, skyview factor

"""

import argparse
import numpy as np
import xarray
#import pandas as pd
#import math
import os


#import pytools.commons_utils.gridlocator as gridlocator

# ---------------------------------------------------------------------
# topographic features: slope, aspect, skyview factor
# ---------------------------------------------------------------------
def read_dem_aoi_xrio(dem_file, bbox_lb_ru, bbox_crs="EPSG:4326", bbox_resxy=None, bbox_outfile=''):
    from rasterio.warp import transform_bounds
    import rioxarray
    from rasterio.enums import Resampling
    from rvt import vis as rvtvis
    
    # rio reading dem
    with rioxarray.open_rasterio(dem_file) as src:
    
        if bbox_crs != src.rio.crs:
            bounds = transform_bounds(bbox_crs, src.rio.crs, *bbox_lb_ru)
        else:
            bounds = bbox_lb_ru
        elv_xr = src.rio.clip_box(bounds[0],bounds[1],bounds[2],bounds[3])
        rio_nodata = src.rio.nodata
        rio_res = np.asarray(src.rio.resolution())
        if not bbox_resxy is None:
            rio_res= np.asanyarray(bbox_resxy)
            elv_xr = elv_xr.rio.reproject(
                elv_xr.rio.crs,
                resolution=bbox_resxy,
                resampling=Resampling.bilinear,
                nodata=rio_nodata
                )
    del src
    
    # slope, aspect, and sky view factor (svf) derived from dem
    if 'int' in str(elv_xr.dtype):
        # image may be in integer, but not good for topographic calculation below
        elv_xr = elv_xr.astype('float32')
        elv_xr.data[np.where(elv_xr.data==rio_nodata)]=np.nan
        elv_xr.rio.write_nodata(np.nan, inplace=True)
        rio_nodata = elv_xr.rio.nodata
        
    elev = elv_xr[0].data
    dem_mask = (elev==rio_nodata)
    sa = rvtvis.slope_aspect(elev, output_units="degree", no_data=rio_nodata)
    slope = np.ma.array(sa['slope'], mask=dem_mask)
    aspect = np.ma.array(sa['aspect'], mask=dem_mask)
    svf = rvtvis.sky_view_factor(elev, resolution=abs(rio_res[0]), no_data=rio_nodata)
    svf = np.ma.array(svf['svf'], mask=dem_mask)
    
    # multi-banded rioxarray
    topo_xr = xarray.concat([elv_xr, elv_xr, elv_xr, elv_xr],
                            dim='band')
    topo_xr = topo_xr.assign_coords(band=[1,2,3,4])
    topo_xr.loc[dict(band=2)] = slope
    topo_xr.loc[dict(band=3)] = aspect
    topo_xr.loc[dict(band=4)] = svf
    #topo_xr.attrs['long_name']=['elevation','slope','aspect','svf']
    
    
    # output raster
    if bbox_outfile != '':
        # the following will only save elevation data to raster. 
        # there is an issue of converting all other data into integer, which not right in raster image.
        topo_xr.loc[dict(band=1)].rio.to_raster(bbox_outfile, dtype=topo_xr.dtype, nodata=float('nan'))
                         
    return topo_xr

# ---------------------------------------------------------------------
# xrio/geopandas used

def percentile_bins_xriolabled(topo_element_xrio, percentiles, topo_xdelta=None):
    xrio_data = topo_element_xrio.to_numpy()
    xmask = ((xrio_data!=topo_element_xrio.rio.nodata) & (~np.isnan(xrio_data)))
    
    bins_labeledid = topo_element_xrio.copy()
    bins_labeledid.data[...] = np.nan
     
    # Nan excluded data bins and indexing
    xrio_data = xrio_data[xmask]
    thresholds = np.percentile(xrio_data, percentiles)
    # duplicates of threshods, if any
    if not np.unique(thresholds).size==thresholds.size:
        idx = np.searchsorted(np.unique(xrio_data), thresholds)
        for i in range(1, thresholds.size):
            if thresholds[i]==thresholds[i-1]: 
                if i<thresholds.size-1:
                    thresholds[i]=np.unique(xrio_data)[idx[i]+1] 
                else:
                    thresholds[i-1]=np.unique(xrio_data)[idx[i]-1]
    
    # if smaller of delta of thresholds than provided-criterial
    if not topo_xdelta is None:
        n = thresholds.size
        if abs(topo_xdelta)>0.0:
            #min. thresholds number
            n = max(1, np.round(abs(np.max(thresholds)-np.min(thresholds))/topo_xdelta))
        
        if n > 1: 
            def recursive_array1d_diff_remv(array1d, xdiff=0.0):                  
                ndiff = np.abs(np.diff(array1d))
                x = np.nanmin(ndiff)
                i = np.nanargmin(ndiff)
                if x<xdiff and ndiff.size>n:
                    if i==0:
                        array1d_new = np.delete(array1d, i+1)
                    elif i==ndiff.size-1:
                        array1d_new = np.delete(array1d, i)
                    else:
                        if ndiff[i-1]>ndiff[i+1]:
                            array1d_new = np.delete(array1d, i+1)
                        else:
                            array1d_new = np.delete(array1d, i)
                
                    return recursive_array1d_diff_remv(array1d_new, xdiff)
                else:
                    return array1d
        
            thresholds = recursive_array1d_diff_remv(thresholds, xdiff=topo_xdelta)
                
    # ID bins
    bins_labeledid.data[xmask] = thresholds_bins_indexing(xrio_data, thresholds)
       
    return bins_labeledid

def thresholds_bins_indexing(topo_element, thresholds_element):
    # note: topo_element must be Nan removed.    
    bins_indx = np.digitize(topo_element, thresholds_element, right=True)
    return bins_indx

def pd_df_groupby_weighted(group, valcol_name="elevation", wtcol_name="area", ret='mean'):
    # Extract values and weights for the current group
    vals = group[valcol_name]
    weights = group[wtcol_name]

    # 1. Calculate the weighted average
    weighted_avg = np.average(vals, weights=weights)
    if ret=='mean':
        return weighted_avg

    # 2. Calculate the weighted variance
    if ret=='std':
        weighted_var = np.average((vals - weighted_avg) ** 2, weights=weights)
        return np.sqrt(weighted_var)

def stats_gridbins_area_weighted(topo_xrio, topo_gid, topo_binid, aggr_to_gid=None, 
                                 topo_area_km2=None, topo_var="elevation", out_gxy=False, out_2d=False):
    # GeoArea-weighted bin statistics
    import pandas
    
    topo_xrio = topo_xrio.drop_vars(['band','spatial_ref'])
    topo_xrio.name = topo_var
    
    topo_binid = topo_binid.drop_vars(['band','spatial_ref'])
    topo_binid.name = 'topounitID'
    uniq_binid = np.unique(topo_binid.to_numpy())
    
    topo_gid = topo_gid.drop_vars(['band','spatial_ref'])
    topo_gid.name = 'gridID'
    #uniq_gid = np.unique(topo_gid.to_numpy())
    
    
    if not topo_area_km2 is None:
        # for calculating area-weights within individual grid
        topo_area_km2 = topo_area_km2.drop_vars(['band','spatial_ref'])
        topo_area_km2.name = 'area_km2'
        topo_df = xarray.merge([topo_xrio,topo_area_km2,topo_binid,topo_gid]).to_dataframe()
    else:
        topo_df = xarray.merge([topo_xrio, topo_gid, topo_binid]).to_dataframe()
    
    aggr_df = {}
    # gridIDs to be aggregated into
    if aggr_to_gid is None:
        aggr_df['x']= topo_gid.x.data
        aggr_df['y']= topo_gid.y.data
        aggr_df['gridID'] = topo_df.groupby('gridID')['gridID'].mean().astype(int)
    else:
        # everything below will be aggregated into this gridID system
        aggr_df['x']= aggr_to_gid.x.data
        aggr_df['y']= aggr_to_gid.y.data
        
        aggr_to_gid = aggr_to_gid.drop_vars(['band','spatial_ref'])
        aggr_to_gid.name = 'gridID'
        aggr_to_gid = aggr_to_gid.to_dataframe()
        aggr_df['gridID'] = aggr_to_gid.groupby('gridID')['gridID'].mean().astype(int)
    
    # need grid's xy coordinates, which are all-pixels' median
    # this is an expensive loop, and may not need by default
    if out_gxy:
        gridxy = topo_df.groupby('gridID')['gridID'].groups
        gridxy = [np.nanmedian([np.array(gridxy[i][j]) for j in range(gridxy[i].size)], axis=0) for i in gridxy.keys()]  # could be long loops - needs a better method    

        gridx  = topo_df.groupby('gridID')['gridID'].mean().astype(float) # as a df as gridID, but renamed and assigned values below 
        gridx.name = 'gridX'
        gridx.values[...]=np.array(gridxy)[:,0]   
        aggr_df['gridX'] = pandas.merge(aggr_df['gridID'], gridx, left_index=True, right_index=True, how='left')

        gridy  = topo_df.groupby('gridID')['gridID'].mean().astype(float) # as a df as gridID, but renamed and assigned values below 
        gridy.name = 'gridY'
        gridy.values[...]=np.array(gridxy)[:,1]   
        aggr_df['gridY'] = pandas.merge(aggr_df['gridID'], gridy, left_index=True, right_index=True, how='left')
    
    # pixel area if any
    if not topo_area_km2 is None:
        area_km2 = topo_df.groupby('gridID')['area_km2'].sum().astype(float)
        area_km2.name = 'area_km2'
        aggr_df['area_km2'] = pandas.merge(aggr_df['gridID'], area_km2, left_index=True, right_index=True, how='left')
    
    # aggregating data bin-by-bin
    #nanbin_group = topo_df.groupby('topounitID', dropna=False).get_group(np.nan) # this didn't work to exclude nan-topoid'ed pixel area
    idx=np.where(np.isnan(topo_df['topounitID'].to_numpy()))
    topo_df['gridID'].values[idx]=np.nan
    for i in range(uniq_binid.size):
        ibin=uniq_binid[i]
        if np.isnan(ibin): continue
        bin_group = topo_df.groupby('topounitID').get_group(ibin)
        if topo_area_km2 is None:
            bin_fraction = bin_group.groupby('gridID')[topo_var].count()/ \
                                topo_df.groupby('gridID')[topo_var].count()    # here, assumming equal area of grid
            bin_mean = bin_group.groupby('gridID')[topo_var].mean()
            bin_std = bin_group.groupby('gridID')[topo_var].std()
        else:
            bin_fraction = bin_group.groupby('gridID')['area_km2'].sum()/ \
                                topo_df.groupby('gridID')['area_km2'].sum()
            bin_mean = bin_group.groupby('gridID'). \
                apply(pd_df_groupby_weighted,
                      valcol_name=topo_var, wtcol_name='area_km2', 
                      ret='mean', include_groups=False)

            bin_std = bin_group.groupby('gridID'). \
                apply(pd_df_groupby_weighted, 
                      valcol_name=topo_var, wtcol_name='area_km2', 
                      ret='std', include_groups=False)
            
        bin_fraction.name = int(ibin)
        bin_mean.name = int(ibin)
        bin_std.name = int(ibin)
        if i==0:
            aggr_df['Area_fraction']  = pandas.merge(aggr_df['gridID'], bin_fraction, left_index=True, right_index=True, how='left')
            aggr_df['Mean_'+topo_var] = pandas.merge(aggr_df['gridID'], bin_mean, left_index=True, right_index=True, how='left')
            aggr_df['Std_'+topo_var]  = pandas.merge(aggr_df['gridID'], bin_std, left_index=True, right_index=True, how='left')
        
        else:
            aggr_df['Area_fraction']  = pandas.merge(aggr_df['Area_fraction'], bin_fraction, left_index=True, right_index=True, how='left')
            aggr_df['Mean_'+topo_var] = pandas.merge(aggr_df['Mean_'+topo_var], bin_mean, left_index=True, right_index=True, how='left')
            aggr_df['Std_'+topo_var] = pandas.merge(aggr_df['Std_'+topo_var], bin_std, left_index=True, right_index=True, how='left')

    return aggr_df

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def grid_topounits(args, 
                     percentiles = np.asarray([10,20,30,40,50,60,70,80,85,90,95,100]), 
                     elm_domainnc='',
                     xrio_topo_out=False):

    from pytools.commons_utils.gridlocator import elmdomain_xrio
    
    from rasterio.enums import Resampling
    #import geopandas
 
    #percentiles = np.asarray([10,20,30,40,50,60,70,80,85,90,95,100])  # total 12 bins

    # boxed boundary from domain.nc, if any
    if not elm_domainnc=='':
        xrio_elmdomain = elmdomain_xrio(elm_domainnc, gridID_redo=True)
        
        bgid = np.where(xrio_elmdomain.band=='GridID')[0][0]
        xrio_elmdomain_gid=xrio_elmdomain[bgid]
        
        # gid as float, so that nan can be used to mask. 
        # and aslo, if gid=0, it will be lost below when do xrio reproject_match
        gid = xrio_elmdomain_gid.to_numpy()+1.0 
        gid[np.where(gid<=0)]=np.nan #mask out 
        xrio_elmdomain_gid.data = gid

        #domain_xr = xarray.open_dataset(elm_domainnc, engine='netcdf4')
        res_x = abs(xrio_elmdomain.rio.resolution()[0])
        res_y = abs(xrio_elmdomain.rio.resolution()[1])
        xlmt = [xrio_elmdomain.rio.bounds()[0],xrio_elmdomain.rio.bounds()[2]]
        ylmt = [xrio_elmdomain.rio.bounds()[1],xrio_elmdomain.rio.bounds()[3]]  
        args.bbox = [np.nanmin(xlmt)-res_x,np.nanmin(ylmt)-res_y,
                     np.nanmax(xlmt)+res_x,np.nanmax(ylmt)+res_y] # extra 1 cell beyond bounds to make sure all grids included
    else:
        xrio_elmdomain = None    
    
    # by rioxarray to read/truncate a geotiff of DEM
    resampleXY=None
    if not args.lonlat_res is None: resampleXY=tuple(args.lonlat_res)
        
    xrio_topo = read_dem_aoi_xrio(
                    args.dem, 
                    tuple(args.bbox), 
                    bbox_resxy=resampleXY, 
                    bbox_outfile=args.bbox_outfile
                    )
        
    # masking topounit ID (labeled xrio arrays)
    elv_xrio  = xrio_topo[0]
    topo_binid = percentile_bins_xriolabled(elv_xrio, percentiles, topo_xdelta=args.group_delta_elev)
    
    #slp_xrio  = xrio_topo[1]
    #slp_binid = percentile_bins_xriolabled(slp_xrio, np.asarray([100]))
    #asp_xrio  = xrio_topo[2]
    #asp_binid = percentile_bins_xriolabled(asp_xrio, np.asarray([50.0,100.0])) # could be merged with elevation-classes (contour bins)
    #svf_xrio  = xrio_topo[3]
    #svf_binid = percentile_bins_xriolabled(svf_xrio, np.asarray([100]))
    
    # assigning grid id to DEM raster cells
    if not elm_domainnc=='':
        topo_gid = xrio_elmdomain_gid.rio.reproject_match(
                                elv_xrio,
                                resampling=Resampling.nearest,
                                nodata=np.nan)
        #
        topo_gid.values = topo_gid.values - 1.0
        if xrio_topo_out:
            return xrio_topo, topo_gid, topo_binid, xrio_elmdomain
        else:
            return None, topo_gid, topo_binid, xrio_elmdomain
    else:
        topo_gid = elv_xrio.copy()
        topo_gid.values[...] = np.indices(elv_xrio.to_numpy().flatten().shape).reshape(topo_gid.shape)
        #
        if xrio_topo_out:
            return xrio_topo, topo_gid, topo_binid, None
        else:
            return None, topo_gid, topo_binid, None

    #

#

def grid_topounits_stats(xrio_topo, xrio_topo_gid, xrio_topo_binid,  
                     xrio_elmdomain_gid=None):

    #
    from pytools.pan_arctic.arctic_domain import elm_km2_from_arcradians2
    from pytools.pan_arctic.arctic_domain import elm_arcradians2_from_km2
    from cmath import pi


    if not xrio_elmdomain_gid is None:
        aggr_to_gid = xrio_elmdomain_gid
    else:
        aggr_to_gid = xrio_topo_gid
        


    # summarizing by bin ids following grid ids, via pandas' xarray dataframe.groupby

    # pixel area in unit of km2, calculated from lat/lon axis. for weighted-averaging
    xrio_topoarea = xrio_topo[0].copy()
    res_x = np.mean(np.diff(xrio_topoarea.x.to_numpy()))
    res_y = np.mean(np.diff(xrio_topoarea.y.to_numpy()))
    res_rad2 = abs(res_x*res_y/180.0/180.0*pi*pi)
    yc = np.meshgrid(xrio_topoarea.x,xrio_topoarea.y)[1]
    xrio_topoarea.data = elm_km2_from_arcradians2(res_rad2, latitude=yc)
    

    elv_df = stats_gridbins_area_weighted(xrio_topo[0], xrio_topo_gid, xrio_topo_binid, aggr_to_gid = aggr_to_gid, \
                                 topo_area_km2=xrio_topoarea, topo_var="elevation", out_gxy=True)
    slp_df = stats_gridbins_area_weighted(xrio_topo[1], xrio_topo_gid, xrio_topo_binid, aggr_to_gid = aggr_to_gid, \
                                 topo_area_km2=xrio_topoarea, topo_var="slope")
    asp_df = stats_gridbins_area_weighted(xrio_topo[2], xrio_topo_gid, xrio_topo_binid, aggr_to_gid = aggr_to_gid, \
                                 topo_area_km2=xrio_topoarea, topo_var="aspect")
    svf_df = stats_gridbins_area_weighted(xrio_topo[3], xrio_topo_gid, xrio_topo_binid, aggr_to_gid = aggr_to_gid, \
                                 topo_area_km2=xrio_topoarea, topo_var="svf")
    
    '''
                "gridID"                : elv_df['gridID'],           # single-column pandas-dataframe, with indexing of lat/lon coordinates
                "Topounit_area_fraction": elv_df['Area_fraction'],    # topounit-columnned pandas-dataframe
                "Topounit_elevation"    : elv_df['Mean_elevation'],
                "Topounit_elevation_std": elv_df['Std_elevation'],
                "Topounit_slope":         slp_df['Mean_slope'],
                "Topounit_slope_std":     slp_df['Std_slope'],
                "Topounit_aspect":        asp_df['Mean_aspect'],
                "Topounit_aspect_std":    asp_df['Std_aspect'],
                "Topounit_skyview":       svf_df['Mean_svf'],
                "Topounit_skyview_std":   svf_df['Std_svf']
    '''
        
    df_topo = {}
    df_topo['gridID'] = elv_df['gridID']
    xc = elv_df['gridX'].drop(columns=['gridID'])
    df_topo['grid_xc'] = xc
    yc = elv_df['gridY'].drop(columns=['gridID'])
    df_topo['grid_yc'] = yc
    if 'bounds' in elv_df.keys():
        df_topo['bounds'] = elv_df['bounds']
    else:
        #bounds are in order of box of LowBot(X,Y), RightUp(X,Y)
        df_topo['bounds'] = [xc-res_x, yc-res_y, xc+res_x, yc+res_y]
    if 'area_km2' in elv_df.keys():
        area_km2 = elv_df['area_km2'].drop(columns=['gridID'])
        df_topo['area_km2'] = area_km2
        df_topo['area'] = elm_arcradians2_from_km2(area_km2)
    if 'x' in elv_df.keys():
        df_topo['x'] = elv_df['x']
    if 'y' in elv_df.keys():
        df_topo['y'] = elv_df['y']
        
    df_topo['Topounit_area_fraction'] = elv_df['Area_fraction'].drop(columns=['gridID'])   
    df_topo['Topounit_elevation']     = elv_df['Mean_elevation'].drop(columns=['gridID'])
    df_topo['Topounit_elevation_std'] = elv_df['Std_elevation'].drop(columns=['gridID']) 
    df_topo['Topounit_slope']     = slp_df['Mean_slope'].drop(columns=['gridID'])
    df_topo['Topounit_slope_std'] = slp_df['Std_slope'].drop(columns=['gridID'])
    df_topo['Topounit_aspect']     = asp_df['Mean_aspect'].drop(columns=['gridID'])
    df_topo['Topounit_aspect_std'] = asp_df['Std_aspect'].drop(columns=['gridID'])
    df_topo['Topounit_skyview']     = svf_df['Mean_svf'].drop(columns=['gridID'])
    df_topo['Topounit_skyview_std'] = svf_df['Std_svf'].drop(columns=['gridID'])
    
    return df_topo


def srf_topo_from_gridstats(pd_topo_data,
                            TOPOUNIT = False, outdata=True, outnc_surf=''):
    
    #from pytools.commons_utils.gridlocator import grids_nearest_or_within
    import pytools.pan_arctic.arctic_domain as elm_domain
    from pytools.pan_arctic.mksrfdata_updatevals import maxtopounits
    from cmath import pi
    
    # directly from grided topo feature dataset    
    srf_grids_sub = pd_topo_data
    
    grd_size = srf_grids_sub['gridID'].values.size
    topounit_size = maxtopounits
    
    # grid-level      
    if 'x' in srf_grids_sub.keys():
        x = srf_grids_sub['x']
    else:
        x = None
    if 'y' in srf_grids_sub.keys():
        y = srf_grids_sub['y']
    else:
        y = None
    if 'grid_xc' in srf_grids_sub.keys():
        xc = srf_grids_sub['grid_xc'].values.squeeze()  # need to squeeze any single-dim
    else:
        xc = None
    if 'grid_yc' in srf_grids_sub.keys():
        yc = srf_grids_sub['grid_yc'].values.squeeze()
    else:
        yc = None
    if 'bounds' in srf_grids_sub.keys():
        nb = len(srf_grids_sub["bounds"])
        bounds = [srf_grids_sub["bounds"][i].values.squeeze() for i in range(nb)]
        bounds = np.transpose(np.asarray(bounds))
    else:
        bounds = None       
    if 'area' in srf_grids_sub.keys():
        area_rad2 = srf_grids_sub['area'].values.squeeze()/180/180*pi*pi
        area_km2 = elm_domain.elm_km2_from_arcradians2(area_rad2, latitude=yc)
        area_rad2 = elm_domain.elm_arcradians2_from_km2(area_km2)
    elif 'area_km2' in srf_grids_sub.keys():
        area_km2 = srf_grids_sub['area_km2'].values.squeeze()
    else:
        area_km2 = None
          

    # topounits-level
    topounit_id   = [np.where(~np.isnan(srf_grids_sub["Topounit_area_fraction"].iloc[i].values))[0] for i in range(grd_size)]
    area_fraction = [srf_grids_sub["Topounit_area_fraction"].iloc[i].values.astype(float) for i in range(grd_size)]  # to np float datatype for convenience
    elevation     = [srf_grids_sub["Topounit_elevation"].iloc[i].values.astype(float) for i in range(grd_size)]    
    elevation_std = [srf_grids_sub["Topounit_elevation_std"].iloc[i].values.astype(float) for i in range(grd_size)]
    slope         = [srf_grids_sub["Topounit_slope"].iloc[i].values.astype(float) for i in range(grd_size)]
    aspect        = [srf_grids_sub["Topounit_aspect"].iloc[i].values.astype(float) for i in range(grd_size)]
    sky_view      = [srf_grids_sub["Topounit_skyview"].iloc[i].values.astype(float) for i in range(grd_size)]               
    
    # need to move all nans to be in trail while reals ahead,
    # and assume reals are continuously ordered with starting-position of topounit_id[:][0]
    area_fraction = [area_fraction[i][~np.isnan(area_fraction[i])] for i in range(grd_size)] # this is not good for ndarray ops
    area_fraction = [np.pad(area_fraction[i],(0, topounit_size-topounit_id[i].size), mode='constant',constant_values=np.nan) for i in range(grd_size)]
    
    elevation     = [elevation[i][~np.isnan(elevation[i])] for i in range(grd_size)]
    elevation     = [np.pad(elevation[i],(0, topounit_size-topounit_id[i].size), mode='constant',constant_values=np.nan) for i in range(grd_size)]
    
    elevation_std = [elevation_std[i][~np.isnan(elevation_std[i])] for i in range(grd_size)]
    elevation_std = [np.pad(elevation_std[i],(0, topounit_size-topounit_id[i].size), mode='constant',constant_values=np.nan) for i in range(grd_size)]
    
    slope         = [slope[i][~np.isnan(slope[i])] for i in range(grd_size)]
    slope         = [np.pad(slope[i],(0, topounit_size-topounit_id[i].size), mode='constant',constant_values=np.nan) for i in range(grd_size)]
    
    aspect        = [aspect[i][~np.isnan(aspect[i])] for i in range(grd_size)]
    aspect        = [np.pad(aspect[i],(0, topounit_size-topounit_id[i].size), mode='constant',constant_values=np.nan) for i in range(grd_size)]
    
    sky_view      = [sky_view[i][~np.isnan(sky_view[i])] for i in range(grd_size)]
    sky_view      = [np.pad(sky_view[i],(0, topounit_size-topounit_id[i].size), mode='constant',constant_values=np.nan) for i in range(grd_size)]
    
    # tranpose to (topounit,gridcell)
    area_fraction = np.transpose(area_fraction)
    elevation     = np.transpose(elevation)
    elevation_std = np.transpose(elevation_std)
    slope         = np.transpose(slope)
    aspect        = np.transpose(aspect)
    sky_view      = np.transpose(sky_view)
    
    #--------------------------------------------------------------------------------------------------------
    # standardard surfdata relevant to topo-features or topounits
    # note: dim 'gridcell', unstructed, can be 2D-meshed as '(lsmlat,lsmlon)'
    
    # by default, the following are included in surfdata, although not really used. 'STD_ELEV' may be named as 'STDEV_ELEV'
    '''
    double SLOPE(gridcell) ;
        SLOPE:_FillValue = NaN ;
        SLOPE:long_name = "mean topographic slope" ;
        SLOPE:units = "degrees" ;
    double STD_ELEV(gridcell) ;
        STD_ELEV:_FillValue = NaN ;
        STD_ELEV:long_name = "standard deviation of elevation" ;
        STD_ELEV:units = "m" ;  
    double TOPO(gridcell) ;
        TOPO:_FillValue = NaN ;
        TOPO:long_name = "mean elevation on land" ;
        TOPO:units = "m" ;
    '''
    
    # topounit on dim 0, averaged for grid-level
    topo1     = np.nanmean(elevation, axis=0, keepdims=True)
    slope1    = np.nanmean(slope, axis=0, keepdims=True)
    topo1_std = np.nanmean(elevation_std, axis=0, keepdims=True)
    aspect1   = np.nanmean(aspect, axis=0, keepdims=True)
    sky_view1 = np.nanmean(sky_view, axis=0, keepdims=True)
        
        
    
    #------------------
    # the following is for downscaling solar radiation 
    '''    
    double SINSL_COSAS(topounit, gridcell) ;
        SINSL_COSAS:_FillValue = NaN ;
        SINSL_COSAS:long_name = "sin(slope) * cos(aspect)" ;
        SINSL_COSAS:units = "unitless" ;
    double SINSL_SINAS(topounit, gridcell) ;
        SINSL_SINAS:_FillValue = NaN ;
        SINSL_SINAS:long_name = "sin(slope) * sin(aspect)" ;
        SINSL_SINAS:units = "unitless" ;
    double SKY_VIEW(topounit, gridcell) ;
        SKY_VIEW:_FillValue = NaN ;
        SKY_VIEW:long_name = "sky view factor" ;
        SKY_VIEW:units = "unitless" ;
    double TERRAIN_CONFIG(topounit, gridcell) ;
        TERRAIN_CONFIG:_FillValue = NaN ;
        TERRAIN_CONFIG:long_name = "terrain configuration factor" ;
        TERRAIN_CONFIG:units = "unitless" ;
    '''

    if TOPOUNIT:
        sinsl_cosas = np.sin(slope)*np.cos(aspect)
        sinsl_sinas = np.sin(slope)*np.sin(aspect)
        # (1+cos(slope))/2-SKYVIEW
        terrain_config = (1.0+np.cos(slope))/2.0 - sky_view

    else:
        sinsl_cosas = np.sin(slope1)*np.cos(aspect1)
        sinsl_sinas = np.sin(slope1)*np.sin(aspect1)

        #TERRAIN_CONFIG: may be estimated, according to Lee et al. (2011) (Eq. (4)), as following: 
        # (1+cos(slope))/2-SKYVIEW
        # Wei-Liang Lee, K.N. Liou, and Alex Hall, 2011. Parameterization of solar fluxes over mountain surfaces for application to climate models. JGR,116, D01101
        terrain_config = (1.0+np.cos(slope1))/2.0 - sky_view1
    
    #------------------
    # the following is for topounit greater than 1.
    # NOTE: all variables are in (topounit, gridcell) or (topounit, lsmlat, lsmlon),
    #       except for 'AREA', 'LATIXY', 'LONGXY', 'TOPO', 'SLOPE', 'STD_ELEV' (or 'STDEV_ELEV'), scalars,
    #                and of course 'MaxTopounitElv', 'TOPO2', 'topoPerGrid'
    '''
    dimensions:
        topounit = 12 ;
    variables:
        float MaxTopounitElv(gridcell) ;
            MaxTopounitElv:_FillValue = NaNf ;
            MaxTopounitElv:descriptions = "Maximum topounits elevation in a gridcell" ;
        float TOPO2(gridcell) ;
            TOPO2:_FillValue = NaNf ;
            TOPO2:descriptions = "Weighted average of topounits elevation in a gridcell" ;
        double TopounitAveElv(topounit, gridcell) ;
            TopounitAveElv:_FillValue = -999. ;
            TopounitAveElv:descriptions = "Average elevation of subgrid" ;
            TopounitAveElv:units = "meters asl" ;
        double TopounitFracArea(topounit, gridcell) ;
            TopounitFracArea:_FillValue = -999. ;
            TopounitFracArea:descriptions = "Fractional area of the subgrid of the land fraction of grid" ;
        int topoPerGrid(gridcell) ;
            topoPerGrid:descriptions = "Number of topounit for each grid" ;    
    '''
    if TOPOUNIT:
        topounit_elev = elevation
        topounit_frac = area_fraction
        
        topo2     = np.nansum(elevation*area_fraction, axis=0)
        topo2_max = np.nanmax(elevation, axis=0)
        topo2_n   = [topounit_id[i].size for i in range(len(topounit_id))]  #topounit_id still a list of grid, with varying id length
        #topo2_n0  = [topounit_id[i][0] for i in range(len(topounit_id))]
        
        
    #--------------------------------------------------------------------------------------------------------

    if outnc_surf!='':
        
        print('write to topounits and/or topo-features to surfdata: ', outnc_surf)
        print('TODO!')
        

    if outdata:
        surfdata = {}
        if not x is None: surfdata['x']   = x # x-axis node
        if not y is None: surfdata['y']   = y # y-axis node
        if not xc is None: surfdata['LONGXY']   = xc
        if not yc is None: surfdata['LATIXY']   = yc
        if not bounds is None: surfdata['bounds']   = bounds
        if not area_km2 is None: 
            surfdata['AREA'] = area_km2
        surfdata['TOPO']     = topo1
        surfdata['STD_ELEV'] = topo1_std
        surfdata['SLOPE']    = slope1
        
        surfdata['SINSL_COSAS']    = sinsl_cosas
        surfdata['SINSL_SINAS']    = sinsl_sinas
        surfdata['SKY_VIEW']       = sky_view1
        surfdata['TERRAIN_CONFIG'] = terrain_config
        
        if TOPOUNIT:
            surfdata['MaxTopounitElv']    = topo2_max
            surfdata['TOPO2']             = topo2
            surfdata['TopounitAveElv']    = topounit_elev
            surfdata['TopounitFracArea']  = topounit_frac
            surfdata['topoPerGrid']       = topo2_n
            #surfdata['topounit0PerGrid']  = topo2_n0 # this is useful to allow lateral connection
        
        
        return surfdata
            
    #
#
# ---------------------------------------------------------------------

def test():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dem", default="global_dem_3s.vrt")
    parser.add_argument("--bbox", nargs=4, type=float, required=True,
                        help="left bottom right top bounding box")
    parser.add_argument("--lonlat_res", nargs=2, type=float, default=None,
                        help="resoultion along longitude and latitude axis within bounding box, if resampling needed")
    parser.add_argument("--sitename", required=True)
    parser.add_argument("--group-delta-elev", type=float, default=100.0,
                        help="Max area-weighted mean elevation difference for recursive grouping")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--bbox_outfile", default="")
    
    args = parser.parse_args()


    #------------------------------------------------------------------------------------    
    IFTGU=False                  # topounit sub-grided or not
    DOMAIN_Extraction = False    # True: extract an area or point; False: generate new domain from template
    DOMAIN_2D = True            # domain/surfdata are in 2D or 1D-unstructured (i.e. drop masked-off grid)
    
    #
    #user_domainnc_2d = 'domain.lnd.r05_RRSwISC6to18E3r5.240328_cavm2d.nc'
    #user_domainnc_1d = 'domain.lnd.r05_RRSwISC6to18E3r5.240328_cavm1d.nc'
    user_domainnc_2d = ''#'domain.lnd.r0125_IcoswISC30E3r5.250918_cavm2d.nc'   
    user_domainnc_1d = ''#'domain.lnd.r0125_IcoswISC30E3r5.250918_cavm1d.nc'
    #user_domainnc='domain.lnd.r05_RRSwISC6to18E3r5.240328_AnakRBurns.nc'
    
    # STEP 1: topographic units from DEM geotiff, and pixel-wide gridID and topounitID
 
    if not IFTGU:
        # if don't do topounit dividing
        xrio, xrio_gid, xrio_tid, xrio_elmdomain = grid_topounits(args, percentiles=np.asarray([100]), 
                        elm_domainnc=user_domainnc_2d,
                        xrio_topo_out=True)
            
    else:
        xrio, xrio_gid, xrio_tid, xrio_elmdomain = grid_topounits(args, 
                        elm_domainnc=user_domainnc_2d,
                        xrio_topo_out=True)
    
    # grid-level topounits statistics
    if not xrio_elmdomain is None:
        bgid = np.where(xrio_elmdomain.band=='GridID')[0][0]
        xrio_elmdomain_gid=xrio_elmdomain[bgid]

        xmask=xrio_elmdomain[np.where(xrio_elmdomain.band=='mask')[0][0]].to_numpy()
        xmask_gid = xrio_elmdomain_gid.to_numpy()
        xmask_gid[np.where(xmask<=0.0)] = np.nan
        xrio_elmdomain_gid.data = xmask_gid # this will exclude non-landmask gridcell below SO the data is to be unstructured
    else:
        if not DOMAIN_2D:
            # removal of NAN data grid-id
            xrio_elmdomain_gid=xrio_gid.where(xrio[0].notnull())
        else:
            # as it is
            xrio_elmdomain_gid=None
        
    gridded_stats = grid_topounits_stats(xrio, xrio_gid, xrio_tid, xrio_elmdomain_gid=xrio_elmdomain_gid)
    
    #------------------------------------------------------------------------------------
    # STEP 2: write surfdata with TOPO features 
    import pytools.pan_arctic.arctic_domain as elm_domain
    import pytools.pan_arctic.mksrfdata_updatevals as srfupdate

    elm_domain.set_e3sm_input('/Users/f9y/project_e3sm/e3sm_inputdata')
 
    
    # example 1: generating template domain and surfdata, from DEM only or user_domainnc, and then filling in new TOPO
    surfdata_topo = srf_topo_from_gridstats(gridded_stats, TOPOUNIT=IFTGU)
        
    xdomain = {}
    xdomain['xc'] = surfdata_topo['LONGXY']
    xdomain['yc'] = surfdata_topo['LATIXY']
    if 'bounds' in surfdata_topo.keys():
        bxy0 = surfdata_topo['bounds'][:,0]
        bxy1 = surfdata_topo['bounds'][:,1]
        bxy2 = surfdata_topo['bounds'][:,2]
        bxy3 = surfdata_topo['bounds'][:,3]
        xdomain['xv'] = np.stack([bxy0, bxy2, bxy2, bxy0], axis=1)
        xdomain['yv'] = np.stack([bxy1, bxy1, bxy3, bxy3], axis=1)
    if 'x' in surfdata_topo.keys():
        xdomain['X_axis'] = surfdata_topo['x']
    if 'y' in surfdata_topo.keys():
        xdomain['Y_axis'] = surfdata_topo['y']
        
    xdomain['frac'] = np.ones_like(xdomain['xc'])
    xdomain['mask'] = np.ones_like(xdomain['xc'])
    if 'TOPO' in surfdata_topo.keys():
        # nan MASKED as 0. NOTE - TOPO is assumed 2-D with dims (TOPOUNITs, GRIDS)
        idx=np.where(np.isnan(surfdata_topo['TOPO'][0,:]))   # assuming at least 1 topounit
        xdomain['mask'][idx]=0
        xdomain['frac'][idx]=0.0    
    if 'AREA' in surfdata_topo.keys():
        xdomain['area_km2'] = surfdata_topo['AREA']
        xdomain['area'] = elm_domain.elm_arcradians2_from_km2(xdomain['area_km2'])    
    
    if user_domainnc_2d=='' and user_domainnc_1d=='': 
        #from DEM
        elm_domain.domain_ncwrite(xdomain, WRITE2D=DOMAIN_2D, ncfile='xdomain.nc', coord_system=True)
        xdomain_ncf='xdomain.nc'
        surfdata_topo['PFTDATA_MASK'] = xdomain['mask']                        
        surfdata_topo['LANDFRAC_PFT'] = xdomain['frac']                        

        if DOMAIN_2D:
            # above xdomain.nc is in 2D. For convenience it's needed for a 1d-unstructured version.
            elm_domain.domain_ncwrite(xdomain, WRITE2D=False, ncfile='xdomain1d.nc', coord_system=True)
                      
    else:
        #from user domain
        xdomain_ncf=user_domainnc_2d
    
        #xc,yc,area directly from user_domain
        # because those in 'surfdata_topo' are DEM pixels aggregated data, which excluded non-values pixels
        bgx = np.where(xrio_elmdomain.band=='xc')[0][0]
        xrio_bgx=xrio_elmdomain[bgx].to_numpy()
        surfdata_topo['LONGXY'] = xrio_bgx[np.where(xmask>0.0)].reshape(surfdata_topo['LONGXY'].shape) # this will overide
        
        bgx = np.where(xrio_elmdomain.band=='yc')[0][0]
        xrio_bgx=xrio_elmdomain[bgx].to_numpy()
        surfdata_topo['LATIXY'] = xrio_bgx[np.where(xmask>0.0)].reshape(surfdata_topo['LATIXY'].shape) # this will overide
        
        bgx = np.where(xrio_elmdomain.band=='area')[0][0]
        xrio_bgx=xrio_elmdomain[bgx].to_numpy()
        area_km2 = xrio_bgx[np.where(xmask>0.0)].reshape(surfdata_topo['LATIXY'].shape)
        area_km2 = elm_domain.elm_km2_from_arcradians2(area_km2, latitude=surfdata_topo['LATIXY'])
        surfdata_topo['AREA'] = area_km2
        
        surfdata_topo['PFTDATA_MASK'] = xmask.reshape(surfdata_topo['LATIXY'].shape)       
    
    # generate surfdata template, and domain if required (e.g. different resolution from user-provided domain)
    #DOMAIN_Extraction = True # if want to overide
    if IFTGU:
        allout = elm_domain.refine_surfdata( \
                lnd_domain_file=elm_domain.DIN_LOC_ROOT+'/share/domains/domain.clm/domain.lnd.r05_RRSwISC6to18E3r5_N45.240328.nc', \
                fsurdat=elm_domain.DIN_LOC_ROOT+'/lnd/clm2/surfdata_map/topounit_surfdata_0.5x0.5_simyr1850_N45.c20220204_TOP.nc', \
                flanduse_timeseries=None, \
                userdomain=xdomain_ncf, \
                SUBDOMAIN_REORDER=(not DOMAIN_Extraction), 
                KEEP_DUPLICATED=(not DOMAIN_Extraction), \
                domain_included=DOMAIN_Extraction)
    else:
        allout = elm_domain.refine_surfdata( \
                lnd_domain_file=elm_domain.DIN_LOC_ROOT+'/share/domains/domain.clm/domain.lnd.r0125_IcoswISC30E3r5.250918_cavm1d.nc', \
                fsurdat=elm_domain.DIN_LOC_ROOT+'/lnd/clm2/surfdata_map/surfdata_0.125x0.125_simyr1850_c250910_TOP_cavm1d.nc', \
                flanduse_timeseries=None, \
                userdomain=xdomain_ncf, \
                SUBDOMAIN_REORDER=(not DOMAIN_Extraction), 
                KEEP_DUPLICATED=(not DOMAIN_Extraction), \
                domain_included=DOMAIN_Extraction)

    
    # the surfdata is named as '*_subset.nc' from above
    os.system('mv '+allout[0]+' surfdata.nc')
   
    
    if DOMAIN_Extraction:
        srf_topo_vars = ['TOPO','STD_ELEV','SLOPE', 'TOPO2','MaxTopounitElv', 'TopounitAveElv','TopounitFracArea','topoPerGrid', \
                    'SINSL_COSAS', 'SINSL_SINAS', 'SKY_VIEW', 'TERRAIN_CONFIG', 'PFTDATA_MASK','LANDFRAC_PFT']
    else:
        srf_topo_vars = ['LONGXY','LATIXY', 'AREA', \
                    'TOPO','STD_ELEV','SLOPE', 'TOPO2','MaxTopounitElv', 'TopounitAveElv','TopounitFracArea','topoPerGrid', \
                    'SINSL_COSAS', 'SINSL_SINAS', 'SKY_VIEW', 'TERRAIN_CONFIG', 'PFTDATA_MASK','LANDFRAC_PFT']
        
    srfupdate.updatevals('surfdata.nc', \
                        user_srf_data=surfdata_topo, user_srf_vars=srf_topo_vars)
    
    print('DONE!')
    
    
   
if __name__ == "__main__":
    test()

