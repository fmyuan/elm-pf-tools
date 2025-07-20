#!/usr/bin/env python3
"""
Jitendra Kumar (kumarj@ornl.gov) 

Elevation percentile based topounits:
 - Area-weighted fractions
 - Recursive grouping of consecutive bins using area-weighted mean
 - Group-level mean and std
 - Group annotation on plots
 - CSV exports (summary + group membership)

This follows the topunit scheme from: 
Tesfa, T. K., Leung, L. R., Thornton, P. E., Brunke, 
M. A., & Duan, Z. (2024). Impacts of Topography‐Based 
Subgrid Scheme and Downscaling of Atmospheric Forcing 
on Modeling Land Surface Processes in the Conterminous US. 
Journal of Advances in Modeling Earth Systems, 16(8). 
https://doi.org/10.1029/2023ms004064

"""

import argparse
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# Read DEM for AOI
# ---------------------------------------------------------------------
def read_dem_aoi(dem_file, bbox, bbox_crs="EPSG:4326"):
    with rasterio.open(dem_file) as src:
        if bbox_crs != src.crs:
            bounds = transform_bounds(bbox_crs, src.crs, *bbox)
        else:
            bounds = bbox
        window = from_bounds(*bounds, transform=src.transform)
        data = src.read(1, window=window, masked=True)
        transform = src.window_transform(window)
        rows, _ = np.indices(data.shape)
        lat2d = transform.f + rows * transform.e
        valid = ~data.mask
        elev = data.data[valid]
        lat = lat2d[valid]
    assert elev.shape == lat.shape
    return elev, lat

# ---------------------------------------------------------------------
# Percentile binning
# ---------------------------------------------------------------------
def elevation_percentile_bins(elev, lat, percentiles):
    thresholds = np.percentile(elev, percentiles)
    bins = {}
    prev = -np.inf
    for p, thr in zip(percentiles, thresholds):
        m = (elev > prev) & (elev <= thr)
        bins[p] = {"elev": elev[m], "lat": lat[m]}
        prev = thr
    return bins

# ---------------------------------------------------------------------
# Area-weighted bin statistics
# ---------------------------------------------------------------------
def summarize_bins_area_weighted(bins):
    records = []
    total_weight = sum(np.cos(np.deg2rad(b["lat"])).sum() for b in bins.values())
    for p, b in bins.items():
        elev = b["elev"]
        lat = b["lat"]
        if elev.size == 0:
            records.append({
                "Percentile_bin": f"≤{p}",
                "Mean_elev": np.nan,
                "Std_elev": np.nan,
                "Area_fraction": 0.0
            })
            continue
        weight = np.cos(np.deg2rad(lat))
        mean_elev = np.sum(elev * weight) / np.sum(weight)
        std_elev = np.sqrt(np.sum(weight * (elev - mean_elev)**2) / np.sum(weight))
        area_frac = np.sum(weight) / total_weight
        records.append({
            "Percentile_bin": f"≤{p}",
            "Mean_elev": mean_elev,
            "Std_elev": std_elev,
            "Area_fraction": area_frac
        })
    return pd.DataFrame(records)

# ---------------------------------------------------------------------
# Recursive grouping with area-weighted mean
# ---------------------------------------------------------------------
def recursive_group_bins_area_weighted(df, bins, delta_elev):
    group_ids = np.zeros(len(df), dtype=int)
    current_group = 1
    start_idx = 0

    print(df)
    # initialize first bin
    merged_mean = df["Mean_elev"].iloc[0]
    merged_area = df["Area_fraction"].iloc[0]

    for i in range(1, len(df)):
        # weighted mean including current bin
        w1 = merged_area
        w2 = df["Area_fraction"].iloc[i]
        candidate_mean = (merged_mean * w1 + df["Mean_elev"].iloc[i] * w2) / (w1 + w2)
        print(f'{i}: df {df["Mean_elev"].iloc[i]} cand {candidate_mean} merged_mean {merged_mean}')

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
# Plot helpers
# ---------------------------------------------------------------------
def draw_group_blocks(ax, groups, alpha=0.15):
    colors = plt.cm.tab20.colors
    for i, (gid, start, end) in enumerate(groups):
        ax.axvspan(start-0.5, end+0.5, color=colors[i % len(colors)], alpha=alpha, zorder=0)

def annotate_groups(ax, groups, y):
    for gid, start, end in groups:
        x = (start + end)/2
        ax.text(x, y, f"G{gid}", ha="center", va="bottom", fontsize=10, fontweight="bold")

# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
def plot_boxplot(bins, df, groups, sitename, bbox, outdir):
    labels = df["Percentile_bin"]
    data = [bins[int(p[1:])]["elev"] for p in labels]
    bbox_str = f"AOI: [{bbox[0]}, {bbox[1]}, {bbox[2]}, {bbox[3]}]"
    fig, ax = plt.subplots(figsize=(12,6))
    draw_group_blocks(ax, groups)
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1,len(labels)+1))
    ax.set_xticklabels(labels, rotation=45)
    ax.set_ylabel("Elevation")
    ax.set_xlabel("Elevation percentile bins")
    ax.set_title(f"{sitename}\n{bbox_str}")
    annotate_groups(ax, groups, ax.get_ylim()[1])
    plt.tight_layout()
    fig.savefig(f"{outdir}/{sitename}_elevation_boxplot.png", dpi=300)
    plt.close(fig)

def plot_summary(df, groups, sitename, bbox, outdir):
    x = np.arange(len(df))
    bbox_str = f"AOI: [{bbox[0]}, {bbox[1]}, {bbox[2]}, {bbox[3]}]"
    fig, ax1 = plt.subplots(figsize=(12,6))
    draw_group_blocks(ax1, groups)
    ax1.errorbar(x, df["Mean_elev"], yerr=df["Std_elev"], fmt="o", capsize=4, label="Mean ± Std")
    ax1.set_ylabel("Elevation")
    ax1.set_xticks(x)
    ax1.set_xticklabels(df["Percentile_bin"], rotation=45)
    ax2 = ax1.twinx()
    ax2.bar(x, df["Area_fraction"], alpha=0.35, width=0.6, label="Area-weighted fraction")
    ax2.set_ylabel("Fraction of AOI area")
    annotate_groups(ax1, groups, ax1.get_ylim()[1])
    h1,l1 = ax1.get_legend_handles_labels()
    h2,l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc="upper left")
    ax1.set_title(f"{sitename}\n{bbox_str}")
    plt.tight_layout()
    fig.savefig(f"{outdir}/{sitename}_elevation_summary.png", dpi=300)
    plt.close(fig)

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def test():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dem", default="global_dem_3s.vrt")
    parser.add_argument("--bbox", nargs=4, type=float, required=True)
    parser.add_argument("--sitename", required=True)
    parser.add_argument("--group-delta-elev", type=float, default=100.0,
                        help="Max area-weighted mean elevation difference for recursive grouping")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    percentiles = [10,20,30,40,50,60,70,80,85,90,95,100]

    elev, lat = read_dem_aoi(args.dem, tuple(args.bbox))
    bins = elevation_percentile_bins(elev, lat, percentiles)

    df = summarize_bins_area_weighted(bins)
    df, groups, group_stats_df = recursive_group_bins_area_weighted(df, bins, args.group_delta_elev)

    # Save CSVs
    df.to_csv(f"{args.outdir}/{args.sitename}_elevation_summary.csv", index=False)
    df[["Percentile_bin","Group_ID"]].to_csv(f"{args.outdir}/{args.sitename}_bin_groups.csv", index=False)
    group_stats_df.to_csv(f"{args.outdir}/{args.sitename}_group_stats.csv", index=False)

    # Print group stats
    print("Group-level statistics:")
    print(group_stats_df)

    plot_boxplot(bins, df, groups, args.sitename, args.bbox, args.outdir)
    plot_summary(df, groups, args.sitename, args.bbox, args.outdir)


if __name__ == "__main__":
    test()

