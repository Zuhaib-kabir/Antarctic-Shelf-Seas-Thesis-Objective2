"""
FIGURE 10 — FULL ORGANIZED LATEST WORKFLOW


FIGURE

Sea-wise Monthly Climatological Cycles
Southern Ocean / Antarctic Shelf Seas, 2008-2025

PANELS

(a) Sea Level Anomaly (SLA)
(b) EN4 total steric height
(c) EN4 thermosteric height
(d) EN4 halosteric height

SOURCE

Extracted from:
    so_sealevel_paper_fig.py

VERSION SELECTION
The supplied source contains one Figure 10 implementation. Therefore this
standalone script uses that complete and latest Figure 10 workflow directly;
there is no older duplicated Figure 10 renderer to retain.

ORGANIZED EXECUTION ORDER

1. Mount Google Drive.
2. Install missing Python packages.
3. Import all required libraries.
4. Define input/output paths.
5. Define analysis settings.
6. Define the 13 Antarctic shelf-sea sectors.
7. Define dataset, mask, weighting and plotting helper functions.
8. Open SLA, total steric, thermosteric and halosteric NetCDF files.
9. Restrict all datasets to the 2008-2025 analysis period.
10. Interpolate SLA to the EN4 grid.
11. Build one fixed SLA-EN4 common mask.
12. Calculate area-weighted monthly time series for each shelf sea.
13. Calculate Jan-Dec climatological cycles.
14. Express each monthly value as departure from that sea's annual mean.
15. Track monthly SLA coverage and identify low-coverage cells.
16. Save the complete analysis table to CSV.
17. Calculate one shared symmetric colour limit for all four panels.
18. Draw four vertically stacked heatmaps.
19. Mark low-coverage SLA cells in panel (a).
20. Save the 1080-dpi PNG and vector PDF.
21. Save the Figure 10 processing summary text.
22. Close all open datasets.

INPUT FILES

1. SLA:
   /content/drive/MyDrive/SAM_Thesis/Data/
   SLA_Antarctic_monthly_2008_2025.nc

2. EN4 total steric:
   /content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
   EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc

3. EN4 thermosteric:
   /content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
   EN4_thermosteric_0_1000m_monthly_2008_2025_SO.nc

4. EN4 halosteric:
   /content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
   EN4_halosteric_0_1000m_monthly_2008_2025_SO.nc

FINAL OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/

    Figure10_sea_wise_monthly_climatological_cycles_1080dpi.png
    Figure10_sea_wise_monthly_climatological_cycles.pdf
    Fig10_sea_wise_monthly_climatological_cycles.csv
    Fig10_processing_summary.txt

ANALYSIS SETTINGS PRESERVED

Study period:
    2008-01-01 through 2025-12-31

SLA core-validity threshold:
    0.70

Low-coverage SLA threshold:
    0.70

Minimum valid years contributing to a monthly climatology:
    5

Output units:
    cm

Heatmap definition:
    Each cell is the monthly departure from that sea's annual mean.

IMPORTANT SOURCE NOTE

The original Figure 10 code itself contains a note saying that the SEA_SECTORS
dictionary should be replaced with the exact sector dictionary used in
Figures 4 and 5 if those earlier analyses used slightly different boundaries.
That source warning is preserved in the workflow below; this organization step
does not silently change those scientific sector definitions.

ORGANIZATION CHANGES
Only structural organization was changed:
    - Google Drive mounting is moved before installation/imports.
    - The complete original Figure 10 processing and plotting workflow follows.
    - Figure 11 and all later code are excluded.

No Figure 10 thresholds, calculations, filenames, masks, sector definitions,
statistical requirements or plotting logic were intentionally changed.

"""

# FIGURE 10 — SEA-WISE MONTHLY CLIMATOLOGICAL CYCLES
#
# 10a: SLA
# 10b: EN4 total steric
# 10c: EN4 thermosteric
# 10d: EN4 halosteric
#
# Output:
#   Four vertically stacked heatmaps
#   X-axis: Jan–Dec
#   Y-axis: 13 Antarctic seas
#   Colour: monthly departure from each sea's annual mean
#   Unit: cm
#
# Required:
#   - common SLA–EN4 mask
#   - mark low-coverage SLA cells
#
# Results section:
#   3.8 Sea-wise monthly cycles


# 1. MOUNT GOOGLE DRIVE

try:
    from google.colab import drive
    drive.mount("/content/drive")
except Exception:
    print("Google Drive mount skipped (not running in Colab).")
# 2. INSTALL MISSING PACKAGES

import sys
import subprocess
import importlib.util

REQUIRED_PACKAGES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "xarray": "xarray",
    "matplotlib": "matplotlib",
    "netCDF4": "netCDF4",
    "h5netcdf": "h5netcdf",
}

missing_packages = [
    pip_name
    for import_name, pip_name in REQUIRED_PACKAGES.items()
    if importlib.util.find_spec(import_name) is None
]

if missing_packages:
    print("Installing missing packages:", ", ".join(missing_packages))
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "-q", *missing_packages]
    )
    print("Packages installed.")
else:
    print("All required packages are already installed.")
# 3. IMPORT LIBRARIES

import os
import warnings
from collections import OrderedDict

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

warnings.filterwarnings("ignore", category=RuntimeWarning)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.weight": "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
})
# 4. COMPLETE FIGURE 10 PROCESSING, OUTPUT CREATION AND PLOTTING

# 3. INPUT / OUTPUT PATHS
SLA_FILE = "/content/drive/MyDrive/SAM_Thesis/Data/SLA_Antarctic_monthly_2008_2025.nc"
TOTAL_FILE = "/content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc"
THERMO_FILE = "/content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/EN4_thermosteric_0_1000m_monthly_2008_2025_SO.nc"
HALO_FILE = "/content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/EN4_halosteric_0_1000m_monthly_2008_2025_SO.nc"

OUTPUT_DIR = "/content/drive/MyDrive/SAM_Thesis/paper2"
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUT_PNG = os.path.join(
    OUTPUT_DIR,
    "Figure10_sea_wise_monthly_climatological_cycles_1080dpi.png",
)
OUT_PDF = os.path.join(
    OUTPUT_DIR,
    "Figure10_sea_wise_monthly_climatological_cycles.pdf",
)
OUT_CSV = os.path.join(
    OUTPUT_DIR,
    "Fig10_sea_wise_monthly_climatological_cycles.csv",
)
OUT_TXT = os.path.join(
    OUTPUT_DIR,
    "Fig10_processing_summary.txt",
)


# 4. USER SETTINGS
ANALYSIS_START = "2008-01-01"
ANALYSIS_END   = "2025-12-31"

# SLA core-validity threshold for the common mask
SLA_CORE_VALID_FRACTION = 0.70

# A sea-month cell in the SLA panel is marked as "low coverage"
# if the valid coverage fraction for that sea and month is below this.
LOW_COVERAGE_THRESHOLD = 0.70

# Minimum number of valid years contributing to a month climatology
MIN_VALID_YEARS_PER_MONTH = 5

SAVE_DPI = 1080

MONTH_NAMES = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
]


# 5. SEA DEFINITIONS
# IMPORTANT:
# Replace these sector boundaries with the EXACT same sea-sector
# dictionary used in your Figures 4 and 5 if your earlier code uses
# slightly different boundaries.
#
# lon_min and lon_max are in degrees east, range [-180, 180].
# If lon_min > lon_max, the sector crosses the dateline.
SEA_SECTORS = OrderedDict([
    ("WED", {"name": "Weddell Sea",           "lon_min": -60,  "lon_max": -20}),
    ("KHV", {"name": "King Haakon VII Sea",   "lon_min": -20,  "lon_max":  10}),
    ("RLS", {"name": "Riiser-Larsen Sea",     "lon_min":  10,  "lon_max":  35}),
    ("LAZ", {"name": "Lazarev Sea",           "lon_min":  35,  "lon_max":  60}),
    ("COS", {"name": "Cosmonauts Sea",        "lon_min":  60,  "lon_max":  90}),
    ("COO", {"name": "Cooperation Sea",       "lon_min":  90,  "lon_max": 115}),
    ("DAV", {"name": "Davis Sea",             "lon_min": 115,  "lon_max": 130}),
    ("MAW", {"name": "Mawson Sea",            "lon_min": 130,  "lon_max": 150}),
    ("DUR", {"name": "D'Urville Sea",         "lon_min": 150,  "lon_max": 170}),
    ("SOM", {"name": "Somov Sea",             "lon_min": 170,  "lon_max": -160}),
    ("ROS", {"name": "Ross Sea",              "lon_min": -160, "lon_max": -130}),
    ("AMU", {"name": "Amundsen Sea",          "lon_min": -130, "lon_max": -100}),
    ("BEL", {"name": "Bellingshausen Sea",    "lon_min": -100, "lon_max":  -60}),
])

SEA_CODES = list(SEA_SECTORS.keys())
N_SEAS = len(SEA_CODES)


# 6. HELPER FUNCTIONS
def open_dataset_safely(file_path):
    attempts = []
    for engine in [None, "netcdf4", "h5netcdf", "scipy"]:
        try:
            kwargs = dict(decode_times=True, mask_and_scale=True)
            if engine is not None:
                kwargs["engine"] = engine
            ds = xr.open_dataset(file_path, **kwargs)
            engine_name = "xarray-default" if engine is None else engine
            print(f"Opened {os.path.basename(file_path)} with engine={engine_name}")
            return ds
        except Exception as e:
            attempts.append(f"{engine}: {e}")

    raise RuntimeError(
        f"Could not open {file_path}\n" + "\n".join(attempts)
    )


def detect_coord_name(ds, kind):
    aliases = {
        "time": ["time"],
        "lat": ["lat", "latitude", "nav_lat", "y"],
        "lon": ["lon", "longitude", "nav_lon", "x"],
    }

    for name in list(ds.coords) + list(ds.variables):
        lname = name.lower()
        var = ds[name]

        std = str(var.attrs.get("standard_name", "")).lower()
        axis = str(var.attrs.get("axis", "")).upper()
        units = str(var.attrs.get("units", "")).lower()

        if kind == "time":
            if lname in aliases["time"] or std == "time" or axis == "T":
                return name

        if kind == "lat":
            if lname in aliases["lat"] or std == "latitude" or axis == "Y" or "degrees_north" in units:
                return name

        if kind == "lon":
            if lname in aliases["lon"] or std == "longitude" or axis == "X" or "degrees_east" in units:
                return name

    raise KeyError(f"Could not detect {kind} coordinate.")


def standardize_dataset(ds):
    time_name = detect_coord_name(ds, "time")
    lat_name = detect_coord_name(ds, "lat")
    lon_name = detect_coord_name(ds, "lon")

    rename_map = {}
    if time_name != "time":
        rename_map[time_name] = "time"
    if lat_name != "lat":
        rename_map[lat_name] = "lat"
    if lon_name != "lon":
        rename_map[lon_name] = "lon"

    if rename_map:
        ds = ds.rename(rename_map)

    # normalize longitude to [-180, 180)
    ds = ds.assign_coords(
        lon=(((ds["lon"].astype(float) + 180) % 360) - 180)
    )
    ds = ds.sortby("lon")
    ds = ds.sortby("lat")

    return ds


def choose_variable(ds, exact_names, contains_names=None):
    available = list(ds.data_vars)

    for name in exact_names:
        if name in available:
            return name

    if contains_names is not None:
        for token in contains_names:
            for name in available:
                if token.lower() in name.lower():
                    return name

    raise KeyError(
        f"Could not detect science variable. Available variables: {available}"
    )


def convert_to_cm(da):
    units = str(da.attrs.get("units", "")).strip().lower()

    if units in ["m", "meter", "metre", "meters", "metres"]:
        out = da * 100.0
        out.attrs["units"] = "cm"
        return out

    if units in ["cm", "centimeter", "centimetre", "centimeters", "centimetres"]:
        out = da.copy()
        out.attrs["units"] = "cm"
        return out

    # If no units are given, assume SLA and steric are already in meters only if values are small.
    vmax = float(np.nanmax(np.abs(da.values)))
    if vmax < 2.0:
        out = da * 100.0
        out.attrs["units"] = "cm"
        return out

    out = da.copy()
    out.attrs["units"] = "cm (assumed)"
    return out


def lon_sector_mask(lon_1d, lon_min, lon_max):
    lon_vals = np.asarray(lon_1d.values)

    if lon_min <= lon_max:
        mask = (lon_vals >= lon_min) & (lon_vals < lon_max)
    else:
        # dateline crossing sector
        mask = (lon_vals >= lon_min) | (lon_vals < lon_max)

    return xr.DataArray(mask, coords={"lon": lon_1d}, dims=["lon"])


def make_2d_sector_mask(lat_1d, lon_1d, lon_min, lon_max):
    lon_mask = lon_sector_mask(lon_1d, lon_min, lon_max)
    lat_mask = xr.DataArray(
        np.ones(lat_1d.size, dtype=bool),
        coords={"lat": lat_1d},
        dims=["lat"],
    )
    return lat_mask & lon_mask


def area_weights_2d(lat_1d, lon_1d):
    w_lat = np.cos(np.deg2rad(lat_1d))
    w2d = xr.DataArray(
        np.repeat(w_lat.values[:, None], lon_1d.size, axis=1),
        coords={"lat": lat_1d, "lon": lon_1d},
        dims=["lat", "lon"],
    )
    return w2d


def weighted_spatial_mean(data_3d, weights_2d, spatial_mask_2d):
    """
    data_3d: (time, lat, lon)
    weights_2d: (lat, lon)
    spatial_mask_2d: boolean (lat, lon)
    returns: time series (time)
    """
    masked = data_3d.where(spatial_mask_2d)
    valid_weights = weights_2d.where(spatial_mask_2d)

    numerator = (masked * valid_weights).sum(dim=("lat", "lon"), skipna=True)
    denominator = valid_weights.where(np.isfinite(masked)).sum(dim=("lat", "lon"), skipna=True)

    return numerator / denominator


def coverage_fraction(data_3d, weights_2d, spatial_mask_2d):
    """
    Coverage fraction for each time step inside a sea mask,
    relative to the total available weight of that sea in the fixed common mask.
    """
    sea_weights_total = weights_2d.where(spatial_mask_2d).sum(dim=("lat", "lon"), skipna=True)

    valid_weight = weights_2d.where(spatial_mask_2d & np.isfinite(data_3d))
    valid_weight_sum = valid_weight.sum(dim=("lat", "lon"), skipna=True)

    return valid_weight_sum / sea_weights_total


def monthly_climatology_and_anomaly(ts):
    """
    ts: time series (time)
    Returns:
      monthly climatology (month=1..12)
      anomaly from the sea's annual mean
      monthly sample count
    """
    monthly_clim = ts.groupby("time.month").mean(dim="time", skipna=True)
    monthly_count = ts.groupby("time.month").count(dim="time")

    annual_mean = monthly_clim.mean(dim="month", skipna=True)
    monthly_anom = monthly_clim - annual_mean

    return monthly_clim, monthly_anom, monthly_count


def nice_symmetric_limit(arrays, percentile=99.0):
    vals = []
    for arr in arrays:
        a = np.asarray(arr, dtype=float)
        a = a[np.isfinite(a)]
        if a.size > 0:
            vals.append(a)

    if not vals:
        return 1.0

    allvals = np.concatenate(vals)
    vmax = float(np.nanpercentile(np.abs(allvals), percentile))

    if vmax <= 1:
        candidates = [0.5, 1.0, 1.5, 2.0]
    elif vmax <= 5:
        candidates = [2.0, 3.0, 4.0, 5.0, 6.0]
    elif vmax <= 10:
        candidates = [6.0, 8.0, 10.0, 12.0]
    elif vmax <= 20:
        candidates = [12.0, 15.0, 20.0]
    else:
        candidates = [20.0, 25.0, 30.0, 40.0]

    for c in candidates:
        if vmax <= c:
            return c

    return np.ceil(vmax)


def draw_low_coverage_cells(ax, low_cov_mask_2d, nrows, ncols):
    """
    Overlay grey hatching on low-coverage SLA cells.
    low_cov_mask_2d shape: (nrows, ncols), True = low coverage
    """
    for i in range(nrows):
        for j in range(ncols):
            if low_cov_mask_2d[i, j]:
                rect = Rectangle(
                    (j - 0.5, i - 0.5),
                    1.0,
                    1.0,
                    facecolor=(0.75, 0.75, 0.75, 0.55),
                    edgecolor="none",
                    hatch=None,
                    zorder=5,
                )
                ax.add_patch(rect)


def style_heatmap_axis(ax, title, ylabels):
    ax.set_title(title, fontsize=16, fontweight="bold", pad=10)

    ax.set_xticks(np.arange(12))
    ax.set_xticklabels(MONTH_NAMES, fontsize=10, fontweight="bold")

    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=11, fontweight="bold")

    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)

    ax.set_xlim(-0.5, 11.5)
    ax.set_ylim(len(ylabels) - 0.5, -0.5)

    # cell grid
    for x in np.arange(-0.5, 12.5, 1):
        ax.axvline(x, color="white", linewidth=1.2, zorder=6)

    for y in np.arange(-0.5, len(ylabels) + 0.5, 1):
        ax.axhline(y, color="white", linewidth=1.2, zorder=6)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


# 7. OPEN DATASETS
ds_sla = standardize_dataset(open_dataset_safely(SLA_FILE))
ds_tot = standardize_dataset(open_dataset_safely(TOTAL_FILE))
ds_the = standardize_dataset(open_dataset_safely(THERMO_FILE))
ds_hal = standardize_dataset(open_dataset_safely(HALO_FILE))

sla_var = choose_variable(ds_sla, exact_names=["sla"], contains_names=["sla"])
tot_var = choose_variable(ds_tot, exact_names=["total_steric_height"], contains_names=["total_steric"])
the_var = choose_variable(ds_the, exact_names=["thermosteric_height"], contains_names=["thermosteric"])
hal_var = choose_variable(ds_hal, exact_names=["halosteric_height"], contains_names=["halosteric"])

print("\nDetected variables:")
print("SLA         :", sla_var)
print("Total       :", tot_var)
print("Thermosteric:", the_var)
print("Halosteric  :", hal_var)


# 8. TIME SUBSET AND COMMON TIME PERIOD
sla = ds_sla[sla_var].sel(time=slice(ANALYSIS_START, ANALYSIS_END))
tot = ds_tot[tot_var].sel(time=slice(ANALYSIS_START, ANALYSIS_END))
the = ds_the[the_var].sel(time=slice(ANALYSIS_START, ANALYSIS_END))
hal = ds_hal[hal_var].sel(time=slice(ANALYSIS_START, ANALYSIS_END))

common_time = np.intersect1d(
    np.intersect1d(sla["time"].values, tot["time"].values),
    np.intersect1d(the["time"].values, hal["time"].values),
)

sla = sla.sel(time=common_time)
tot = tot.sel(time=common_time)
the = the.sel(time=common_time)
hal = hal.sel(time=common_time)

print("\nCommon analysis months:", sla.sizes["time"])


# 9. USE EN4 GRID AS THE COMMON GRID
# Restrict SLA to EN4 latitude range before interpolation
sla = sla.sel(
    lat=slice(float(tot["lat"].min()), float(tot["lat"].max()))
)

# Interpolate SLA to EN4 grid
sla_on_en4 = sla.interp(
    lat=tot["lat"],
    lon=tot["lon"],
    method="linear",
)

# Convert all to cm
sla_on_en4 = convert_to_cm(sla_on_en4)
tot = convert_to_cm(tot)
the = convert_to_cm(the)
hal = convert_to_cm(hal)

print("Units after conversion:")
print("SLA         :", sla_on_en4.attrs.get("units", ""))
print("Total       :", tot.attrs.get("units", ""))
print("Thermosteric:", the.attrs.get("units", ""))
print("Halosteric  :", hal.attrs.get("units", ""))


# 10. BUILD FIXED COMMON SLA–EN4 MASK
nmonths = int(sla_on_en4.sizes["time"])
min_valid_months = int(np.ceil(SLA_CORE_VALID_FRACTION * nmonths))

sla_core_mask = np.isfinite(sla_on_en4).sum(dim="time") >= min_valid_months

# EN4 common validity across all three steric variables
en4_core_mask = (
    (np.isfinite(tot).sum(dim="time") >= min_valid_months) &
    (np.isfinite(the).sum(dim="time") >= min_valid_months) &
    (np.isfinite(hal).sum(dim="time") >= min_valid_months)
)

common_mask = sla_core_mask & en4_core_mask

print("\nCommon fixed mask summary:")
print("Total grid cells on EN4 grid:", common_mask.size)
print("Common valid cells         :", int(common_mask.sum().values))


# 11. SEA-WISE MONTHLY SERIES
weights2d = area_weights_2d(tot["lat"], tot["lon"])

sla_matrix = np.full((N_SEAS, 12), np.nan)
tot_matrix = np.full((N_SEAS, 12), np.nan)
the_matrix = np.full((N_SEAS, 12), np.nan)
hal_matrix = np.full((N_SEAS, 12), np.nan)

sla_monthly_count_matrix = np.full((N_SEAS, 12), np.nan)
low_coverage_matrix = np.zeros((N_SEAS, 12), dtype=bool)

output_rows = []

for i, sea_code in enumerate(SEA_CODES):
    sector = SEA_SECTORS[sea_code]
    sea_mask = make_2d_sector_mask(
        tot["lat"],
        tot["lon"],
        sector["lon_min"],
        sector["lon_max"],
    )

    # final sea mask = sector ∩ common mask
    sea_common_mask = sea_mask & common_mask

    # if a sea has zero valid cells, skip
    if int(sea_common_mask.sum().values) == 0:
        print(f"WARNING: {sea_code} has zero common valid cells.")
        continue

    # weighted monthly time series for each variable
    sla_ts = weighted_spatial_mean(sla_on_en4, weights2d, sea_common_mask)
    tot_ts = weighted_spatial_mean(tot,      weights2d, sea_common_mask)
    the_ts = weighted_spatial_mean(the,      weights2d, sea_common_mask)
    hal_ts = weighted_spatial_mean(hal,      weights2d, sea_common_mask)

    # SLA coverage fraction per time step inside this sea
    sla_cov_ts = coverage_fraction(sla_on_en4, weights2d, sea_common_mask)

    # monthly climatology and anomaly from annual mean
    sla_clim, sla_anom, sla_count = monthly_climatology_and_anomaly(sla_ts)
    tot_clim, tot_anom, tot_count = monthly_climatology_and_anomaly(tot_ts)
    the_clim, the_anom, the_count = monthly_climatology_and_anomaly(the_ts)
    hal_clim, hal_anom, hal_count = monthly_climatology_and_anomaly(hal_ts)

    sla_cov_clim = sla_cov_ts.groupby("time.month").mean(dim="time", skipna=True)

    for m in range(1, 13):
        if int(sla_count.sel(month=m).values) >= MIN_VALID_YEARS_PER_MONTH:
            sla_matrix[i, m-1] = float(sla_anom.sel(month=m).values)
            sla_monthly_count_matrix[i, m-1] = float(sla_count.sel(month=m).values)
        else:
            sla_matrix[i, m-1] = np.nan

        if int(tot_count.sel(month=m).values) >= MIN_VALID_YEARS_PER_MONTH:
            tot_matrix[i, m-1] = float(tot_anom.sel(month=m).values)

        if int(the_count.sel(month=m).values) >= MIN_VALID_YEARS_PER_MONTH:
            the_matrix[i, m-1] = float(the_anom.sel(month=m).values)

        if int(hal_count.sel(month=m).values) >= MIN_VALID_YEARS_PER_MONTH:
            hal_matrix[i, m-1] = float(hal_anom.sel(month=m).values)

        cov_val = float(sla_cov_clim.sel(month=m).values)
        low_coverage_matrix[i, m-1] = np.isfinite(cov_val) and (cov_val < LOW_COVERAGE_THRESHOLD)

        output_rows.append({
            "sea_code": sea_code,
            "sea_name": sector["name"],
            "month": m,
            "month_name": MONTH_NAMES[m-1],
            "sla_anomaly_cm": sla_matrix[i, m-1],
            "total_steric_anomaly_cm": tot_matrix[i, m-1],
            "thermosteric_anomaly_cm": the_matrix[i, m-1],
            "halosteric_anomaly_cm": hal_matrix[i, m-1],
            "sla_valid_years": sla_monthly_count_matrix[i, m-1],
            "sla_low_coverage_flag": int(low_coverage_matrix[i, m-1]),
        })


# 12. SAVE OUTPUT TABLE
df_out = pd.DataFrame(output_rows)
df_out.to_csv(OUT_CSV, index=False)
print("\nSaved processed table:")
print(OUT_CSV)


# 13. SHARED COLOR LIMIT FOR THE FOUR PANELS
shared_limit = nice_symmetric_limit(
    [sla_matrix, tot_matrix, the_matrix, hal_matrix],
    percentile=99.0
)

print("\nShared symmetric colour limit (cm):", shared_limit)


# 14. PLOT FOUR VERTICALLY STACKED HEATMAPS
fig = plt.figure(figsize=(10, 14))

gs = fig.add_gridspec(
    nrows=4,
    ncols=2,
    width_ratios=[24, 1.2],
    height_ratios=[1, 1, 1, 1],
    left=0.10,
    right=0.93,
    bottom=0.06,
    top=0.97,
    hspace=0.24,
    wspace=0.05,
)

heatmap_data = [
    sla_matrix,
    tot_matrix,
    the_matrix,
    hal_matrix,
]

panel_titles = [
    "(a) SLA",
    "(b) Total steric",
    "(c) Thermosteric",
    "(d) Halosteric",
]

axes = []
caxes = []

for r in range(4):
    ax = fig.add_subplot(gs[r, 0])
    cax = fig.add_subplot(gs[r, 1])
    axes.append(ax)
    caxes.append(cax)

cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad(color="0.78")

images = []

for idx, (ax, cax, data, title) in enumerate(zip(axes, caxes, heatmap_data, panel_titles)):
    im = ax.imshow(
        data,
        aspect="auto",
        cmap=cmap,
        vmin=-shared_limit,
        vmax=shared_limit,
        interpolation="nearest",
        origin="upper",
    )
    images.append(im)

    style_heatmap_axis(ax, title, SEA_CODES)

    # only panel (a): mark low-coverage SLA cells
    if idx == 0:
        draw_low_coverage_cells(
            ax,
            low_coverage_matrix,
            nrows=N_SEAS,
            ncols=12,
        )

    cb = plt.colorbar(im, cax=cax)
    cb.set_label("cm", fontsize=12, fontweight="bold", rotation=0, labelpad=8)
    cb.ax.tick_params(labelsize=10, width=1.0, length=3)
    cb.formatter = FuncFormatter(lambda x, pos: f"{x:.0f}")
    cb.update_ticks()

    for tick in cb.ax.get_yticklabels():
        tick.set_fontweight("bold")



# Common x label only for the bottom panel
axes[-1].set_xlabel("Month", fontsize=12, fontweight="bold", labelpad=8)

# Figure-level y label
fig.text(
    0.03,
    0.5,
    "Antarctic shelf sea",
    rotation=90,
    va="center",
    ha="center",
    fontsize=13,
    fontweight="bold",
)

# Save
fig.savefig(OUT_PNG, dpi=SAVE_DPI, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("\nSaved figure:")
print(OUT_PNG)
print(OUT_PDF)


# 15. SAVE PROCESSING SUMMARY
summary_lines = [
    "=" * 80,
    "FIGURE 10 — SEA-WISE MONTHLY CLIMATOLOGICAL CYCLES",
    "=" * 80,
    "",
    f"SLA file         : {SLA_FILE}",
    f"Total steric     : {TOTAL_FILE}",
    f"Thermosteric     : {THERMO_FILE}",
    f"Halosteric       : {HALO_FILE}",
    "",
    f"SLA variable     : {sla_var}",
    f"Total variable   : {tot_var}",
    f"Thermo variable  : {the_var}",
    f"Halo variable    : {hal_var}",
    "",
    f"Common months    : {nmonths}",
    f"SLA core fraction threshold : {SLA_CORE_VALID_FRACTION}",
    f"Low-coverage SLA threshold : {LOW_COVERAGE_THRESHOLD}",
    f"Min valid years per month  : {MIN_VALID_YEARS_PER_MONTH}",
    f"Shared colour limit (cm)   : ±{shared_limit}",
    "",
    f"Common valid grid cells    : {int(common_mask.sum().values)}",
    "",
    f"Output CSV : {OUT_CSV}",
    f"Output PNG : {OUT_PNG}",
    f"Output PDF : {OUT_PDF}",
    "",
    "IMPORTANT:",
    "If your Figures 4 and 5 used a different sector dictionary,",
    "replace SEA_SECTORS in this script with that exact same definition.",
]

with open(OUT_TXT, "w", encoding="utf-8") as f:
    f.write("\n".join(summary_lines))

print("\n".join(summary_lines))


# 16. CLOSE DATASETS
ds_sla.close()
ds_tot.close()
ds_the.close()
ds_hal.close()

print("\nFIGURE 10 COMPLETED SUCCESSFULLY.")
