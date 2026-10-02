# Seasonal Steric Sea-Level Components with Objective Mapping
#
# Columns:
#   1 = Total steric height
#   2 = Thermosteric height
#   3 = Halosteric height
#
# Rows:
#   Spring, Summer, Autumn, Winter
#
# Method:
#   - Seasonal climatology from Argo-derived steric NetCDF files
#   - Local objective mapping / optimal interpolation
#   - Small black dots show original valid Argo grid cells
#
# Output:
# /content/drive/MyDrive/SAM_Thesis/Fig/Obj2Fig15_SO_ObjectiveMapped_corrected.png


# Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')


# 0) Install
!pip -q install xarray netCDF4 h5netcdf dask cartopy scipy

# 1) Imports
import os
import gc
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker

from scipy.spatial import cKDTree

import cartopy.crs as ccrs
import cartopy.feature as cfeature


# 2) Input / output paths
total_steric_file = "/content/drive/MyDrive/SAM_Thesis/Processed/total_steric_height_monthly_2008_2025_SO.nc"
thermo_steric_file = "/content/drive/MyDrive/SAM_Thesis/Processed/thermosteric_height_monthly_2008_2025_SO.nc"
halo_steric_file = "/content/drive/MyDrive/SAM_Thesis/Processed/halosteric_height_monthly_2008_2025_SO.nc"

out_dir = "/content/drive/MyDrive/SAM_Thesis/Fig"
os.makedirs(out_dir, exist_ok=True)

OUT_PATH = os.path.join(out_dir, "Obj2Fig15_SO_ObjectiveMapped_corrected.png")

for f in [total_steric_file, thermo_steric_file, halo_steric_file]:
    if not os.path.exists(f):
        raise FileNotFoundError(f"File not found: {f}")

print("All input steric files found.")


# 3) Settings
CLIM_START_YEAR = 2008
CLIM_END_YEAR = 2025

SEASONS = {
    "Spring": [9, 10, 11],   # SON
    "Summer": [12, 1, 2],    # DJF
    "Autumn": [3, 4, 5],     # MAM
    "Winter": [6, 7, 8],     # JJA
}

ROW_ORDER = ["Spring", "Summer", "Autumn", "Winter"]

# Convert steric height from metre to cm
STERIC_TO_CM = 100.0

# Objective mapping settings — corrected, less over-smoothed
LENGTH_SCALE_KM = 350.0
SEARCH_RADIUS_KM = 800.0
MIN_NEIGHBORS = 3
MAX_NEIGHBORS = 20
NOISE_TO_SIGNAL = 0.25

# Original-data dot overlay
SHOW_VALID_DOTS = True
DOT_SIZE = 0.7
DOT_ALPHA = 0.30
DOT_COLOR = "black"

# Colorbar settings
# Use fixed ranges to avoid extreme percentile stretching
USE_PERCENTILE_LIMITS = False

TOTAL_VMIN, TOTAL_VMAX = -30, 30
THERMO_VMIN, THERMO_VMAX = -25, 25
HALO_VMIN, HALO_VMAX = -15, 15

TOTAL_TICKS = [-30, -15, 0, 15, 30]
THERMO_TICKS = [-25, -12.5, 0, 12.5, 25]
HALO_TICKS = [-15, -7.5, 0, 7.5, 15]

TOTAL_CMAP = "RdBu_r"
THERMO_CMAP = "RdBu_r"
HALO_CMAP = "RdBu_r"

# Figure layout
FIGSIZE = (14.2, 15.5)

LEFT = 0.10
RIGHT = 0.97
TOP = 0.95
BOTTOM = 0.04
WSPACE = 0.15
HSPACE = 0.17

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "savefig.dpi": 300
})


# 4) Helper functions
def standardize_coords(ds):
    rename_dict = {}

    if "latitude" in ds.coords:
        rename_dict["latitude"] = "lat"

    if "longitude" in ds.coords:
        rename_dict["longitude"] = "lon"

    if "valid_time" in ds.coords:
        rename_dict["valid_time"] = "time"

    ds = ds.rename(rename_dict)

    if "lon" in ds.coords and float(ds["lon"].max()) > 180:
        ds = ds.assign_coords(
            lon=(((ds["lon"] + 180) % 360) - 180)
        )

    if "lat" in ds.coords:
        ds = ds.sortby("lat")

    if "lon" in ds.coords:
        ds = ds.sortby("lon")

    if "time" in ds.coords:
        ds = ds.sortby("time")

    return ds


def find_var(ds, possible_names):
    for name in possible_names:
        if name in ds.data_vars:
            return name

    vars_list = [v for v in ds.data_vars if v != "nprof"]

    if len(vars_list) == 1:
        print(f"Using only variable found: {vars_list[0]}")
        return vars_list[0]

    raise ValueError(
        f"Could not find variable from {possible_names}. "
        f"Available variables: {list(ds.data_vars)}"
    )


def get_season_year(time_da, months):
    """
    DJF rule:
    Dec 2024 + Jan 2025 + Feb 2025 = Summer 2025
    """
    year = time_da.dt.year

    if 12 in months and 1 in months:
        year = xr.where(time_da.dt.month == 12, year + 1, year)

    return year.rename("season_year")


def seasonal_mean_by_year(da, months, min_months=2):
    da_season = da.where(da["time"].dt.month.isin(months), drop=True)

    if da_season.sizes.get("time", 0) == 0:
        return None

    sy = get_season_year(da_season["time"], months)
    da_season = da_season.assign_coords(season_year=sy)

    mean_da = da_season.groupby("season_year").mean("time", skipna=True)

    years = sy.values
    unique_years = np.unique(years)

    valid_years = []

    for y in unique_years:
        n_months = np.sum(years == y)
        if n_months >= min_months:
            valid_years.append(int(y))

    if len(valid_years) == 0:
        return None

    mean_da = mean_da.sel(season_year=valid_years)

    return mean_da


def seasonal_climatology(da, months):
    yearly = seasonal_mean_by_year(da, months, min_months=2)

    if yearly is None:
        raise ValueError("No valid seasonal data found.")

    clim = yearly.sel(
        season_year=slice(CLIM_START_YEAR, CLIM_END_YEAR)
    ).mean("season_year", skipna=True)

    return clim.compute()


def get_lonlat_mesh(ds):
    lon_vals = ds["lon"].values
    lat_vals = ds["lat"].values
    return np.meshgrid(lon_vals, lat_vals)


def lonlat_to_unit_xyz(lon_deg, lat_deg):
    lon_rad = np.deg2rad(lon_deg)
    lat_rad = np.deg2rad(lat_deg)

    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)

    return np.column_stack([x, y, z])


def chord_radius_from_km(radius_km):
    earth_radius_km = 6371.0
    angle = radius_km / earth_radius_km
    return 2.0 * np.sin(angle / 2.0)


def great_circle_km_from_chord(chord):
    earth_radius_km = 6371.0
    chord = np.clip(chord, 0.0, 2.0)
    angle = 2.0 * np.arcsin(chord / 2.0)
    return earth_radius_km * angle


def objective_map_field(raw_da, lon2, lat2,
                        length_scale_km=350.0,
                        search_radius_km=800.0,
                        min_neighbors=3,
                        max_neighbors=20,
                        noise_to_signal=0.25):
    """
    Local objective mapping / optimal interpolation for sparse Argo-derived fields.

    This is for visualization. Original valid cells are shown by dots.
    """

    raw = raw_da.values.astype("float64")
    valid = np.isfinite(raw)

    mapped = np.full_like(raw, np.nan, dtype="float64")

    if np.sum(valid) < min_neighbors:
        return xr.DataArray(
            mapped,
            coords=raw_da.coords,
            dims=raw_da.dims,
            attrs=raw_da.attrs
        )

    obs_lon = lon2[valid]
    obs_lat = lat2[valid]
    obs_val = raw[valid]

    bg = np.nanmean(obs_val)
    obs_anom = obs_val - bg

    target_lon = lon2.ravel()
    target_lat = lat2.ravel()

    obs_xyz = lonlat_to_unit_xyz(obs_lon, obs_lat)
    target_xyz = lonlat_to_unit_xyz(target_lon, target_lat)

    tree = cKDTree(obs_xyz)
    search_chord = chord_radius_from_km(search_radius_km)

    sigma2 = np.nanvar(obs_anom)

    if not np.isfinite(sigma2) or sigma2 <= 0:
        sigma2 = 1.0

    mapped_flat = mapped.ravel()

    for ti, target_point in enumerate(target_xyz):

        neighbor_idx = tree.query_ball_point(target_point, r=search_chord)

        if len(neighbor_idx) < min_neighbors:
            continue

        neighbor_idx = np.asarray(neighbor_idx, dtype=int)

        neighbor_xyz = obs_xyz[neighbor_idx]
        chord_to_target = np.linalg.norm(
            neighbor_xyz - target_point[None, :],
            axis=1
        )

        order = np.argsort(chord_to_target)

        if len(order) > max_neighbors:
            order = order[:max_neighbors]

        neighbor_idx = neighbor_idx[order]
        neighbor_xyz = obs_xyz[neighbor_idx]
        obs_local = obs_anom[neighbor_idx]

        n = len(neighbor_idx)

        if n < min_neighbors:
            continue

        diff = neighbor_xyz[:, None, :] - neighbor_xyz[None, :, :]
        chord_oo = np.linalg.norm(diff, axis=2)
        dist_oo = great_circle_km_from_chord(chord_oo)

        chord_ot = np.linalg.norm(neighbor_xyz - target_point[None, :], axis=1)
        dist_ot = great_circle_km_from_chord(chord_ot)

        C_oo = sigma2 * np.exp(-0.5 * (dist_oo / length_scale_km) ** 2)
        C_ot = sigma2 * np.exp(-0.5 * (dist_ot / length_scale_km) ** 2)

        C_oo = C_oo + np.eye(n) * (noise_to_signal * sigma2)

        try:
            weights = np.linalg.solve(C_oo, C_ot)
            estimate = bg + np.dot(weights, obs_local)
            mapped_flat[ti] = estimate

        except np.linalg.LinAlgError:
            w = np.exp(-0.5 * (dist_ot / length_scale_km) ** 2)
            if np.sum(w) > 0:
                mapped_flat[ti] = bg + np.sum(w * obs_local) / np.sum(w)

    mapped = mapped_flat.reshape(raw.shape)

    mapped_da = xr.DataArray(
        mapped.astype("float32"),
        coords=raw_da.coords,
        dims=raw_da.dims,
        attrs=raw_da.attrs
    )

    mapped_da.attrs["mapping_method"] = "local objective mapping / optimal interpolation"
    mapped_da.attrs["length_scale_km"] = length_scale_km
    mapped_da.attrs["search_radius_km"] = search_radius_km
    mapped_da.attrs["noise_to_signal"] = noise_to_signal

    return mapped_da


def symmetric_percentile_limit(list_of_maps, percentile=98):
    vals = []

    for da in list_of_maps:
        arr = da.values
        arr = arr[np.isfinite(arr)]
        if arr.size > 0:
            vals.append(arr)

    if len(vals) == 0:
        return -1, 1, [-1, 0, 1]

    vals = np.concatenate(vals)
    lim = np.nanpercentile(np.abs(vals), percentile)

    if lim <= 5:
        lim = np.ceil(lim)
    elif lim <= 20:
        lim = np.ceil(lim / 2) * 2
    else:
        lim = np.ceil(lim / 5) * 5

    ticks = np.linspace(-lim, lim, 5)

    return -lim, lim, ticks


def polar_ax(fig, nrows, ncols, idx):
    proj = ccrs.SouthPolarStereo()
    ax = fig.add_subplot(nrows, ncols, idx, projection=proj)

    theta = np.linspace(0, 2 * np.pi, 240)
    center = [0.5, 0.5]
    radius = 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())

    ax.add_feature(
        cfeature.LAND,
        facecolor="0.88",
        edgecolor="black",
        linewidth=0.45,
        zorder=3
    )

    ax.coastlines(linewidth=0.55, zorder=4)

    lon_grid = [-180, -150, -120, -90, -60, -30, 0, 30, 60, 120, 150]
    lat_grid = [-60, -70, -80]

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        linewidth=0.35,
        linestyle=":",
        color="0.55",
        alpha=0.60,
        zorder=5
    )

    gl.xlocator = mticker.FixedLocator(lon_grid)
    gl.ylocator = mticker.FixedLocator(lat_grid)

    edge_lat = -57.5

    for lo in lon_grid:
        if lo < 0:
            label = f"{abs(lo)}°W"
        elif lo > 0:
            label = f"{lo}°E"
        else:
            label = "0°"

        ax.text(
            lo,
            edge_lat,
            label,
            transform=ccrs.PlateCarree(),
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            zorder=10
        )

    for la in [-70, -80]:
        ax.text(
            0,
            la,
            f"{abs(la)}°S",
            transform=ccrs.PlateCarree(),
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            zorder=10
        )

    return ax


def add_valid_dots(ax, lon2, lat2, raw_da):
    valid = np.isfinite(raw_da.values)

    if np.any(valid):
        ax.scatter(
            lon2[valid],
            lat2[valid],
            s=DOT_SIZE,
            c=DOT_COLOR,
            alpha=DOT_ALPHA,
            transform=ccrs.PlateCarree(),
            linewidths=0,
            zorder=6
        )


# 5) Open datasets
ds_total = xr.open_dataset(
    total_steric_file,
    decode_times=True,
    chunks={"time": 12}
)

ds_thermo = xr.open_dataset(
    thermo_steric_file,
    decode_times=True,
    chunks={"time": 12}
)

ds_halo = xr.open_dataset(
    halo_steric_file,
    decode_times=True,
    chunks={"time": 12}
)

ds_total = standardize_coords(ds_total)
ds_thermo = standardize_coords(ds_thermo)
ds_halo = standardize_coords(ds_halo)

print("\nTotal steric dataset:")
print(ds_total)

print("\nThermosteric dataset:")
print(ds_thermo)

print("\nHalosteric dataset:")
print(ds_halo)


# 6) Detect variables
total_var = find_var(
    ds_total,
    ["total_steric_height", "steric_height", "total_steric", "steric"]
)

thermo_var = find_var(
    ds_thermo,
    ["thermosteric_height", "thermosteric", "thermo_steric_height", "thermo_steric"]
)

halo_var = find_var(
    ds_halo,
    ["halosteric_height", "halosteric", "halo_steric_height", "halo_steric"]
)

print("\nDetected variables:")
print("Total steric :", total_var)
print("Thermosteric :", thermo_var)
print("Halosteric   :", halo_var)


# 7) Prepare variables
total_steric = ds_total[total_var].astype("float32") * STERIC_TO_CM
thermo_steric = ds_thermo[thermo_var].astype("float32") * STERIC_TO_CM
halo_steric = ds_halo[halo_var].astype("float32") * STERIC_TO_CM

total_steric.attrs["units"] = "cm"
thermo_steric.attrs["units"] = "cm"
halo_steric.attrs["units"] = "cm"

lon2, lat2 = get_lonlat_mesh(ds_total)


# 8) Calculate seasonal raw climatology maps
raw_total = {}
raw_thermo = {}
raw_halo = {}

for season_name, months in SEASONS.items():

    print(f"\nCalculating raw seasonal climatology: {season_name}")

    raw_total[season_name] = seasonal_climatology(total_steric, months)
    raw_thermo[season_name] = seasonal_climatology(thermo_steric, months)
    raw_halo[season_name] = seasonal_climatology(halo_steric, months)

    print("  Total valid cells:", int(np.isfinite(raw_total[season_name]).sum()))
    print("  Thermo valid cells:", int(np.isfinite(raw_thermo[season_name]).sum()))
    print("  Halo valid cells:", int(np.isfinite(raw_halo[season_name]).sum()))

    gc.collect()


# 9) Objective mapping
maps_total = {}
maps_thermo = {}
maps_halo = {}

for season_name in ROW_ORDER:

    print(f"\nObjective mapping: {season_name}")

    maps_total[season_name] = objective_map_field(
        raw_total[season_name],
        lon2,
        lat2,
        length_scale_km=LENGTH_SCALE_KM,
        search_radius_km=SEARCH_RADIUS_KM,
        min_neighbors=MIN_NEIGHBORS,
        max_neighbors=MAX_NEIGHBORS,
        noise_to_signal=NOISE_TO_SIGNAL
    )

    maps_thermo[season_name] = objective_map_field(
        raw_thermo[season_name],
        lon2,
        lat2,
        length_scale_km=LENGTH_SCALE_KM,
        search_radius_km=SEARCH_RADIUS_KM,
        min_neighbors=MIN_NEIGHBORS,
        max_neighbors=MAX_NEIGHBORS,
        noise_to_signal=NOISE_TO_SIGNAL
    )

    maps_halo[season_name] = objective_map_field(
        raw_halo[season_name],
        lon2,
        lat2,
        length_scale_km=LENGTH_SCALE_KM,
        search_radius_km=SEARCH_RADIUS_KM,
        min_neighbors=MIN_NEIGHBORS,
        max_neighbors=MAX_NEIGHBORS,
        noise_to_signal=NOISE_TO_SIGNAL
    )

    print("  Total mapped cells:", int(np.isfinite(maps_total[season_name]).sum()))
    print("  Thermo mapped cells:", int(np.isfinite(maps_thermo[season_name]).sum()))
    print("  Halo mapped cells:", int(np.isfinite(maps_halo[season_name]).sum()))

    gc.collect()

print("\nObjective mapping finished.")


# 10) Color limits
if USE_PERCENTILE_LIMITS:
    TOTAL_VMIN, TOTAL_VMAX, TOTAL_TICKS = symmetric_percentile_limit(
        [maps_total[s] for s in ROW_ORDER],
        PERCENTILE_LIMIT
    )

    THERMO_VMIN, THERMO_VMAX, THERMO_TICKS = symmetric_percentile_limit(
        [maps_thermo[s] for s in ROW_ORDER],
        PERCENTILE_LIMIT
    )

    HALO_VMIN, HALO_VMAX, HALO_TICKS = symmetric_percentile_limit(
        [maps_halo[s] for s in ROW_ORDER],
        PERCENTILE_LIMIT
    )

else:
    TOTAL_TICKS = TOTAL_TICKS
    THERMO_TICKS = THERMO_TICKS
    HALO_TICKS = HALO_TICKS

print("\nColor limits used:")
print("Total:", TOTAL_VMIN, TOTAL_VMAX)
print("Thermo:", THERMO_VMIN, THERMO_VMAX)
print("Halo:", HALO_VMIN, HALO_VMAX)


# 11) Plot figure
fig = plt.figure(figsize=FIGSIZE)

letters = list("abcdefghijkl")
letter_i = 0

axs = np.empty((4, 3), dtype=object)

for r, season_name in enumerate(ROW_ORDER):

    # Column 1: Total steric height
    ax1 = polar_ax(fig, 4, 3, r * 3 + 1)
    axs[r, 0] = ax1

    pcm1 = ax1.pcolormesh(
        lon2,
        lat2,
        maps_total[season_name],
        transform=ccrs.PlateCarree(),
        cmap=TOTAL_CMAP,
        vmin=TOTAL_VMIN,
        vmax=TOTAL_VMAX,
        shading="auto",
        zorder=1
    )

    if SHOW_VALID_DOTS:
        add_valid_dots(ax1, lon2, lat2, raw_total[season_name])

    ax1.text(
        -0.08,
        1.02,
        f"({letters[letter_i]})",
        transform=ax1.transAxes,
        ha="left",
        va="bottom",
        fontsize=16,
        fontweight="bold"
    )

    letter_i += 1

    cb1 = fig.colorbar(
        pcm1,
        ax=ax1,
        shrink=0.82,
        pad=0.02,
        ticks=TOTAL_TICKS
    )

    cb1.set_label("cm", fontsize=12, fontweight="bold")
    cb1.ax.tick_params(labelsize=9)

    for t in cb1.ax.get_yticklabels():
        t.set_fontweight("bold")

    # Column 2: Thermosteric height
    ax2 = polar_ax(fig, 4, 3, r * 3 + 2)
    axs[r, 1] = ax2

    pcm2 = ax2.pcolormesh(
        lon2,
        lat2,
        maps_thermo[season_name],
        transform=ccrs.PlateCarree(),
        cmap=THERMO_CMAP,
        vmin=THERMO_VMIN,
        vmax=THERMO_VMAX,
        shading="auto",
        zorder=1
    )

    if SHOW_VALID_DOTS:
        add_valid_dots(ax2, lon2, lat2, raw_thermo[season_name])

    ax2.text(
        -0.08,
        1.02,
        f"({letters[letter_i]})",
        transform=ax2.transAxes,
        ha="left",
        va="bottom",
        fontsize=16,
        fontweight="bold"
    )

    letter_i += 1

    cb2 = fig.colorbar(
        pcm2,
        ax=ax2,
        shrink=0.82,
        pad=0.02,
        ticks=THERMO_TICKS
    )

    cb2.set_label("cm", fontsize=11.3, fontweight="bold")
    cb2.ax.tick_params(labelsize=9)

    for t in cb2.ax.get_yticklabels():
        t.set_fontweight("bold")

    # Column 3: Halosteric height
    ax3 = polar_ax(fig, 4, 3, r * 3 + 3)
    axs[r, 2] = ax3

    pcm3 = ax3.pcolormesh(
        lon2,
        lat2,
        maps_halo[season_name],
        transform=ccrs.PlateCarree(),
        cmap=HALO_CMAP,
        vmin=HALO_VMIN,
        vmax=HALO_VMAX,
        shading="auto",
        zorder=1
    )

    if SHOW_VALID_DOTS:
        add_valid_dots(ax3, lon2, lat2, raw_halo[season_name])

    ax3.text(
        -0.08,
        1.02,
        f"({letters[letter_i]})",
        transform=ax3.transAxes,
        ha="left",
        va="bottom",
        fontsize=16,
        fontweight="bold"
    )

    letter_i += 1

    cb3 = fig.colorbar(
        pcm3,
        ax=ax3,
        shrink=0.82,
        pad=0.02,
        ticks=HALO_TICKS
    )

    cb3.set_label("cm", fontsize=12, fontweight="bold")
    cb3.ax.tick_params(labelsize=9)

    for t in cb3.ax.get_yticklabels():
        t.set_fontweight("bold")


# 12) Layout labels
plt.subplots_adjust(
    left=LEFT,
    right=RIGHT,
    top=TOP,
    bottom=BOTTOM,
    wspace=WSPACE,
    hspace=HSPACE
)

column_titles = [
    "Total steric height",
    "Thermosteric height",
    "Halosteric height"
]

for j, title in enumerate(column_titles):
    pos = axs[0, j].get_position()
    x_center = 0.5 * (pos.x0 + pos.x1)
    y_top = pos.y1

    fig.text(
        x_center,
        y_top + 0.015,
        title,
        ha="center",
        va="bottom",
        fontsize=18,
        fontweight="bold"
    )

for r, season_name in enumerate(ROW_ORDER):
    pos0 = axs[r, 0].get_position()
    pos2 = axs[r, 2].get_position()

    y_center = 0.5 * (min(pos0.y0, pos2.y0) + max(pos0.y1, pos2.y1))
    x_left = pos0.x0

    fig.text(
        x_left - 0.07,
        y_center,
        season_name,
        rotation=90,
        ha="center",
        va="center",
        fontsize=22,
        fontweight="bold"
    )


# 13) Save output
plt.savefig(
    OUT_PATH,
    dpi=1080,
    bbox_inches="tight"
)

plt.show()

print("\nSaved figure:")
print(OUT_PATH)

ds_total.close()
ds_thermo.close()
ds_halo.close()
