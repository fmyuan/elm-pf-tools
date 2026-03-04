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
import pandas as pd
import math


# ---------------------------------------------------------------------
# topographic features: slope, aspect, skyview factor
# ---------------------------------------------------------------------
def read_dem_aoi_xrio(dem_file, bbox_lb_ru, bbox_crs="EPSG:4326", bbox_outfile=''):
    from rasterio.warp import transform_bounds
    import rioxarray
    import xarray
    from rvt import vis as rvtvis
    
    # rio reading dem
    with rioxarray.open_rasterio(dem_file, mask_and_scale=True) as src:
    
        if bbox_crs != src.rio.crs:
            bounds = transform_bounds(bbox_crs, src.rio.crs, *bbox_lb_ru)
        else:
            bounds = bbox_lb_ru
        elv_xr = src.rio.clip_box(bounds[0],bounds[1],bounds[2],bounds[3])
        rio_nodata = src.rio.nodata
        rio_res = np.asarray([np.mean(np.diff(elv_xr.y)), np.mean(np.diff(elv_xr.x))])    
    del src
    
    # slope, aspect, and sky view factor (svf) derived from dem
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
# Percentile (thresholds) binning
# ---------------------------------------------------------------------
def percentile_bins(topo_element, xlon, ylat, percentiles):
    assert topo_element.shape == xlon.shape == ylat.shape
    thresholds = np.percentile(topo_element, percentiles)    
    return thresholds_bins(topo_element, xlon, ylat, thresholds, percentiles)

def thresholds_bins(topo_element, xlon, ylat, thresholds_element, thresholds_label):    
    assert topo_element.shape == xlon.shape == ylat.shape
    assert thresholds_element.shape == thresholds_label.shape
    
    bins_indx = np.digitize(topo_element, thresholds_element, right=True)
    bins = {}
    for ip in range(thresholds_label.size):
        m = (bins_indx==ip)
        bins[thresholds_label[ip]] = {"masked": m, 
                                      "topo_element": topo_element[m], 
                                      "lon":xlon[m] ,"lat": ylat[m]}
        
    return bins

# ---------------------------------------------------------------------
# GeoArea-weighted bin statistics
# ---------------------------------------------------------------------
def summarize_bins_area_weighted(bins, topo="elev"):
    records = []
    total_weight = sum(np.cos(np.deg2rad(b["lat"])).sum() for b in bins.values())
    for p, b in bins.items():
        elev = b[topo]
        lat = b["lat"]
        if elev.size == 0:
            records.append({
                "Percentile_bin": f"≤{p}",
                "Mean_"+topo: np.nan,
                "Std_"+topo: np.nan,
                "Area_fraction": 0.0
            })
            continue
        weight = np.cos(np.deg2rad(lat))
        mean_elev = np.sum(elev * weight) / np.sum(weight)
        std_elev = np.sqrt(np.sum(weight * (elev - mean_elev)**2) / np.sum(weight))
        area_frac = np.sum(weight) / total_weight
        records.append({
            "Percentile_bin": f"≤{p}",
            "Mean_"+topo: mean_elev,
            "Std_"+topo: std_elev,
            "Area_fraction": area_frac
        })
    return pd.DataFrame(records)

# ---------------------------------------------------------------------
# Recursive grouping with area-weighted mean
# ---------------------------------------------------------------------
def recursive_group_bins_area_weighted(df, delta_elev):
    group_ids = np.zeros(len(df), dtype=int)
    current_group = 1
    start_idx = 0

    #print(df)
    # initialize first bin
    merged_mean = df["Mean_elev"].iloc[0]
    merged_area = df["Area_fraction"].iloc[0]

    for i in range(1, len(df)):
        # weighted mean including current bin
        w1 = merged_area
        w2 = df["Area_fraction"].iloc[i]
        candidate_mean = (merged_mean * w1 + df["Mean_elev"].iloc[i] * w2) / (w1 + w2)
        #print(f'{i}: df {df["Mean_elev"].iloc[i]} cand {candidate_mean} merged_mean {merged_mean}')

        #if abs(candidate_mean - merged_mean) <= delta_elev:
        if abs(candidate_mean - df["Mean_elev"].iloc[i]) <= delta_elev:
            merged_mean = candidate_mean
            merged_area += w2
        else:
            group_ids[start_idx:i] = current_group
            current_group += 1
            start_idx = i
            merged_mean = df["Mean_elev"].iloc[i]
            merged_area = df["Area_fraction"].iloc[i]

    group_ids[start_idx:] = current_group
    df["Group_ID"] = group_ids

    # Generate group tuples for plotting
    groups = []
    for gid in range(1, current_group + 1):
        idxs = np.where(group_ids == gid)[0]
        groups.append((gid, idxs[0], idxs[-1]))

    # Compute group-level mean, std, area_fraction
    group_stats = []
    for gid, start, end in groups:
        bins_in_group = df.iloc[start:end+1]
        # area-weighted mean
        weights = bins_in_group["Area_fraction"].values
        if np.sum(weights)<=0: continue
        
        means = bins_in_group["Mean_elev"].values
        stds = bins_in_group["Std_elev"].values
        group_mean = np.sum(weights * means) / np.sum(weights)
        # area-weighted std including within-bin variance
        group_std = np.sqrt(np.sum(weights * (stds**2 + (means - group_mean)**2)) / np.sum(weights))
        group_area = np.sum(weights)
                
        group_stats.append({
            "Group_ID": gid,
            "Group_Mean": group_mean,
            "Group_Std": group_std,
            "Group_Area_fraction": group_area            
        })
    group_stats_df = pd.DataFrame(group_stats)
    return df, groups, group_stats_df


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def grided_topounits(args, versbose=0):
    
    percentiles = np.asarray([10,20,30,40,50,60,70,80,85,90,95,100])

    # by rioxarray to read/truncate a geotiff of DEM
    xrio_topo = read_dem_aoi_xrio(args.dem, tuple(args.bbox), bbox_outfile=args.bbox_outfile)
    elev = xrio_topo[0].data
    xy = np.meshgrid(xrio_topo.x.data, xrio_topo.y.data)
    lon = xy[0][elev!=xrio_topo.rio.nodata]
    lat = xy[1][elev!=xrio_topo.rio.nodata]
    elev = elev[elev!=xrio_topo.rio.nodata]
    
    # grouping for user-defined resolution grids within bounding box.
    # (default: one grid within bounding box)
    ylat = np.asarray([np.min(lat), np.max(lat)])
    xlon = np.asarray([np.min(lon), np.max(lon)])
    if args.lonlat_res!=None:
        dy = tuple(args.lonlat_res)[1]
        ny = math.ceil((np.max(lat)-np.min(lat))/dy)  # math.ceil() will be over bbox, but this what we want.
        ylat = np.min(lat)+np.asarray(range(ny+1))*dy

        dx = tuple(args.lonlat_res)[0]
        nx = math.ceil((np.max(lon)-np.min(lon))/dx)
        xlon = np.min(lon)+np.asarray(range(nx+1))*dx
    
    # contour-binning for whole bounding-box
    contour_bins = percentile_bins(elev, lon, lat, percentiles)
    
    # grided summarizing
    grid_stats = []
    for iy in range(ylat.size-1):
        for ix in range(xlon.size-1):
            g_mask = ((lat>=ylat[iy]) & (lat<=ylat[iy+1])) & \
                ((lon>=xlon[ix]) & (lon<=xlon[ix+1]))
            if not any(g_mask) or sum(g_mask)/g_mask.size<0.001: continue #skip if none or less than 0.01%
            
            bins_masked = []
            bins = {}
            for p in contour_bins.keys():
                contour_mask = contour_bins[p]['masked']
                assert g_mask.shape == contour_mask.shape
                m = (g_mask & contour_mask)
                bins_masked.append(m)
                bins[p]={'elev': elev[m],
                         'lat': lat[m], 'lon':lon[m]}
                
            df = summarize_bins_area_weighted(bins)
            df, groups, group_stats_df = recursive_group_bins_area_weighted(df, args.group_delta_elev)
            
            # need to merge bins' pixel mask after grouping
            # NOTE: this logical mask is useful to retrieve pixel position of group contained. 
            #       e.g. for each group, we could get its pixel location by lat[m]/lon[m] pairs
            bins_masked = np.asarray(bins_masked)
            group_pixels_masked = {}
            for grpid, start, end in groups:
                pmasked = np.any(bins_masked[start:end+1,:], axis=0)
                if np.any(pmasked): group_pixels_masked[grpid] = pmasked
            
            #
            grid_stats.append({
                "Xindex": ix,
                "Yindex": iy,
                "Longitude_range": np.asarray([xlon[ix],xlon[ix+1]]),
                "Latitude_range": np.asarray([ylat[iy],ylat[iy+1]]),
                "Topounit_pixels_mask": group_pixels_masked,
                "Topounit_id:": group_stats_df['Group_ID'],
                "Topounit_elevation": group_stats_df['Group_Mean'],
                "Topounit_elevation_std": group_stats_df['Group_Std'],
                "Topounit_area_fraction": group_stats_df['Group_Area_fraction']
                })
            
            # Print group stats
            print("Group-level statistics for grid: ", iy, ix, ylat[iy],'~~',ylat[iy+1], xlon[ix],'~~',xlon[ix+1])
            if versbose>0:
                print(group_stats_df,'\n')
    
    #
    return xrio_topo, grid_stats

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
    
    xrio, gridded_stats = grided_topounits(args, versbose=0)


if __name__ == "__main__":
    test()

