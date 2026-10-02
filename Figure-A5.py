"""

SUPPLEMENTARY FIGURE S1 — FULL ORGANIZED LATEST WORKFLOW

FIGURE

Argo Coverage, Southern Ocean
2001-2025

SOURCE

Extracted from:
    so_sealevel_paper_fig.py

VERSION SELECTION
The supplied source contains one Supplementary Figure S1 implementation.
That implementation is already explicitly described as:

Therefore, there is no older duplicated S1 plotting block to retain or merge.
This standalone script preserves the complete corrected S1 workflow and only
reorders the startup so that Google Drive mounting appears first.

PANEL STRUCTURE
---------------
(a) Spring profile-count map
(b) Summer profile-count map
(c) Autumn profile-count map
(d) Winter profile-count map
(e) Sea x season profile-count heatmap
(f) Percentage of profiles reaching 1000 and 2000 dbar

IMPORTANT DEPTH-SUPPORT CORRECTION
The Argo source dataset contains:

    PRES_GRID (N_LEVELS)
    TEMP      (N_PROF, N_LEVELS)
    PSAL      (N_PROF, N_LEVELS)

The corrected S1 workflow therefore calculates maximum pressure separately for
each profile as the deepest PRES_GRID level at which BOTH temperature and
salinity are finite.

A profile is classified as reaching 1000 or 2000 dbar when at least one common
valid temperature-salinity sample exists at or below the corresponding pressure.

INPUT FILE
/content/drive/MyDrive/SAM_Thesis/Data/
    argo_SO_profiles_2001_2025_cleaned_gridded.nc

FINAL OUTPUT DIRECTORY
/content/drive/MyDrive/SAM_Thesis/paper2/

FINAL OUTPUTS
Figure:
    Figure_S1_Argo_coverage_2001_2025_FINAL_LAYOUT_1080dpi.png
    Figure_S1_Argo_coverage_2001_2025_FINAL_LAYOUT.pdf

Supporting data:
    Figure_S1_Argo_coverage_2001_2025_FINAL_LAYOUT_sea_season_counts.csv
    Figure_S1_Argo_coverage_2001_2025_FINAL_LAYOUT_depth_threshold_percentages.csv
    Figure_S1_Argo_coverage_2001_2025_FINAL_LAYOUT_seasonal_profile_count_maps.nc

Optional profile-level output:
    Figure_S1_Argo_coverage_2001_2025_FINAL_LAYOUT_profile_summary.csv.gz

Processing record:
    Figure_S1_Argo_coverage_2001_2025_FINAL_LAYOUT_processing_log.txt

ORGANIZED EXECUTION ORDER
1. Mount Google Drive.
2. Install required Python packages when missing.
3. Import all required libraries.
4. Define Argo input and output paths.
5. Define Southern Ocean region and 1-degree map grid.
6. Define Southern Hemisphere seasons.
7. Define the 13 Antarctic shelf-sea longitude sectors.
8. Define figure/layout settings.
9. Define safe NetCDF opening and coordinate/time helpers.
10. Convert Argo JULD to timestamps.
11. Normalize longitude to [-180, 180).
12. Calculate profile depth support from PRES_GRID/TEMP/PSAL.
13. Build the filtered profile-level summary.
14. Create seasonal 1-degree profile-count maps.
15. Create the 13-sea x 4-season profile-count table.
16. Calculate sea-wise percentages reaching 1000 and 2000 dbar.
17. Save the sea-season CSV.
18. Save the depth-threshold CSV.
19. Save the seasonal-count NetCDF.
20. Optionally save the compressed profile-level summary.
21. Determine the shared logarithmic map scale.
22. Build the corrected six-panel S1 layout.
23. Keep a dedicated spacer between panels (e) and (f).
24. Place panel (f) title on its own two-line title axis.
25. Save the vector PDF.
26. Save the fixed-canvas 1080-dpi PNG.
27. Save the processing log.
28. Close/release large data resources.

LATEST LAYOUT CORRECTIONS PRESERVED
-----------------------------------
The corrected source layout retains:

    - no overall figure title;
    - no bottom explanatory note;
    - a dedicated spacer between panel (e)'s colorbar and panel (f)'s y labels;
    - a dedicated two-line title axis for panel (f);
    - four seasonal polar coverage maps;
    - shared logarithmic profile-count scaling;
    - fixed-canvas 1080-dpi output without tight-bbox expansion.

SOURCE SETTINGS PRESERVED
Region:
    latitude  -90 to -60 degrees
    longitude -180 to 180 degrees

Map grid:
    1 degree

Profile chunk size:
    5000

Profile-level summary:
    enabled by default in the source

Seasons:
    Spring = SON
    Summer = DJF
    Autumn = MAM
    Winter = JJA

The original 13 shelf-sea longitude sectors are retained exactly as written in
the S1 source. The source itself notes that these should be kept identical to
Figures 4 and 5 if those figures use slightly different sector limits.

ORGANIZATION CHANGES
Only structural organization has been changed:

    - Google Drive mounting is moved to the beginning.
    - Package installation follows the mount.
    - Library imports follow package installation.
    - The complete original corrected S1 workflow then follows.
    - Figure S3 and all later source code are excluded.

No S1 pressure-support rule, seasons, sea-sector definitions, thresholds,
gridding, output filenames, plotting logic, or layout corrections are
intentionally changed.

"""

# SUPPLEMENTARY FIGURE S1 — ARGO COVERAGE
# Full corrected code for:
#   argo_SO_profiles_2001_2025_cleaned_gridded.nc
#
# Panels:
#   (a) Spring profile-count map
#   (b) Summer profile-count map
#   (c) Autumn profile-count map
#   (d) Winter profile-count map
#   (e) Sea × season profile-count heatmap
#   (f) Percentage reaching 1000 and 2000 dbar
#
# IMPORTANT CORRECTION
# --------------------
# This dataset has:
#   PRES_GRID (N_LEVELS)
#   TEMP      (N_PROF, N_LEVELS)
#   PSAL      (N_PROF, N_LEVELS)
#
# Therefore, maximum pressure is calculated profile-by-profile as the deepest
# PRES_GRID level at which BOTH temperature and salinity are finite.
#
# No overall figure title is added.
# Output is saved at 1080 dpi.

# 1. MOUNT GOOGLE DRIVE
try:
    from google.colab import drive
    drive.mount("/content/drive", force_remount=False)
except Exception:
    print("Google Drive mounting skipped because this is not Google Colab.")


# 2. INSTALL REQUIRED PACKAGES WHEN NEEDED
import sys
import subprocess
import importlib.util

REQUIRED_PACKAGES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "xarray": "xarray",
    "matplotlib": "matplotlib",
    "cartopy": "cartopy",
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
else:
    print("All required packages are already installed.")
# 3. IMPORT LIBRARIES

import os
import gc
import warnings
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
from matplotlib.path import Path as MplPath

import cartopy.crs as ccrs
import cartopy.feature as cfeature

warnings.filterwarnings("ignore", category=RuntimeWarning)


# 4. COMPLETE CORRECTED SUPPLEMENTARY FIGURE S1 WORKFLOW

# 3. USER SETTINGS
ARGO_FILE = Path(
    "/content/drive/MyDrive/SAM_Thesis/Data/"
    "argo_SO_profiles_2001_2025_cleaned_gridded.nc"
)

OUTPUT_DIR = Path("/content/drive/MyDrive/SAM_Thesis/paper2")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_STEM = "Figure_S1_Argo_coverage_2001_2025_FINAL_LAYOUT"

OUTPUT_PNG = OUTPUT_DIR / f"{OUTPUT_STEM}_1080dpi.png"
OUTPUT_PDF = OUTPUT_DIR / f"{OUTPUT_STEM}.pdf"
OUTPUT_HEATMAP_CSV = OUTPUT_DIR / f"{OUTPUT_STEM}_sea_season_counts.csv"
OUTPUT_DEPTH_CSV = OUTPUT_DIR / f"{OUTPUT_STEM}_depth_threshold_percentages.csv"
OUTPUT_SEASONAL_NC = OUTPUT_DIR / f"{OUTPUT_STEM}_seasonal_profile_count_maps.nc"
OUTPUT_PROFILE_SUMMARY = OUTPUT_DIR / f"{OUTPUT_STEM}_profile_summary.csv.gz"
OUTPUT_LOG = OUTPUT_DIR / f"{OUTPUT_STEM}_processing_log.txt"

# Region and gridding
LAT_MIN = -90.0
LAT_MAX = -60.0
LON_MIN = -180.0
LON_MAX = 180.0

MAP_GRID_RESOLUTION_DEG = 1.0
PROFILE_CHUNK_SIZE = 5000

# Set False after the first successful run if you want to avoid writing the
# large profile-level compressed CSV. Aggregated CSV/NetCDF outputs are always
# written.
SAVE_PROFILE_LEVEL_SUMMARY = True

# Southern Hemisphere seasons
SEASONS = OrderedDict([
    ("Spring", [9, 10, 11]),
    ("Summer", [12, 1, 2]),
    ("Autumn", [3, 4, 5]),
    ("Winter", [6, 7, 8]),
])

# 13 Antarctic shelf-sea longitude sectors.
# These are configurable. Keep them identical to the boundaries used in your
# main-paper Figures 4 and 5 if those figures use slightly different limits.
# A wrapped sector is represented by lon_min > lon_max.
SEA_BOUNDS = OrderedDict([
    ("WED", (-60.0,  -20.0)),
    ("KHV", (-20.0,   10.0)),
    ("RLS", ( 10.0,   35.0)),
    ("LAZ", ( 35.0,   60.0)),
    ("COS", ( 60.0,   90.0)),
    ("COO", ( 90.0,  115.0)),
    ("DAV", (115.0,  130.0)),
    ("MAW", (130.0,  150.0)),
    ("DUR", (150.0,  170.0)),
    ("SOM", (170.0, -160.0)),
    ("ROS", (-160.0, -130.0)),
    ("AMU", (-130.0, -100.0)),
    ("BEL", (-100.0,  -60.0)),
])

SEA_ORDER = list(SEA_BOUNDS.keys())

# Figure settings
SAVE_DPI = 1080
FIGSIZE = (11.2, 7.8)

FONT_PANEL = 12.2
FONT_MAP_TICK = 6.5
FONT_AXIS = 9.2
FONT_TICK = 8.0
FONT_CELL = 8.3
FONT_LEGEND = 8.2
FONT_NOTE = 8.3

LAND_COLOR = "#d9d9d9"
COAST_COLOR = "#666666"
GRID_COLOR = "#8c8c8c"

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.weight": "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "none",
})

# DATASET
def open_dataset_safely(path):
    attempts = []

    for engine in [None, "netcdf4", "h5netcdf", "scipy"]:
        try:
            kwargs = {
                "decode_times": False,
                "mask_and_scale": True,
                "cache": False,
            }

            if engine is not None:
                kwargs["engine"] = engine

            dataset = xr.open_dataset(path, **kwargs)

            print(
                f"Opened {Path(path).name} with "
                f"engine={engine or 'xarray-default'}"
            )

            return dataset

        except Exception as error:
            attempts.append(f"{engine}: {error}")

    raise RuntimeError(
        f"Could not open file:\n{path}\n\n" + "\n".join(attempts)
    )


def normalize_longitude(longitude):
    longitude = np.asarray(longitude, dtype=np.float64)
    return ((longitude + 180.0) % 360.0) - 180.0


def convert_juld_to_datetime(juld_data_array):
    """
    Convert Argo JULD to monthly-compatible pandas timestamps.

    Standard Argo JULD is days since 1950-01-01 00:00:00 UTC.
    The code first checks the variable's units attribute. If units are absent,
    it uses the standard Argo epoch and prints a warning.
    """
    values = np.asarray(juld_data_array.values)

    if np.issubdtype(values.dtype, np.datetime64):
        return pd.DatetimeIndex(pd.to_datetime(values))

    units = str(juld_data_array.attrs.get("units", "")).strip()
    calendar = str(juld_data_array.attrs.get("calendar", "standard")).strip()

    if "since" in units.lower():
        temporary = xr.Dataset(
            {
                "JULD": xr.DataArray(
                    values,
                    dims=("N_PROF",),
                    attrs={"units": units, "calendar": calendar},
                )
            }
        )

        try:
            decoded = xr.decode_cf(temporary)["JULD"].values
            return pd.DatetimeIndex(pd.to_datetime(decoded))
        except Exception as error:
            print(
                "WARNING: CF decoding of JULD failed; using standard Argo "
                f"epoch instead. Reason: {error}"
            )

    if not units:
        print(
            "WARNING: JULD has no units attribute. The standard Argo epoch "
            "1950-01-01 is being used."
        )
    else:
        print(
            f"WARNING: JULD units {units!r} were not decoded. The standard "
            "Argo epoch 1950-01-01 is being used."
        )

    numeric = pd.to_numeric(
        pd.Series(values.reshape(-1)),
        errors="coerce",
    ).to_numpy(dtype=float)

    result = pd.DatetimeIndex(
        pd.Timestamp("1950-01-01")
        + pd.to_timedelta(numeric, unit="D")
    )

    return result


def month_to_season(month):
    for season_name, months in SEASONS.items():
        if int(month) in months:
            return season_name
    return np.nan


def longitude_in_sector(longitude, lon_min, lon_max):
    if lon_min <= lon_max:
        return (longitude >= lon_min) & (longitude < lon_max)

    return (longitude >= lon_min) | (longitude < lon_max)


def assign_sea_codes(longitudes, latitudes):
    longitudes = np.asarray(longitudes, dtype=float)
    latitudes = np.asarray(latitudes, dtype=float)

    output = np.full(longitudes.shape, None, dtype=object)

    in_latitude_range = (
        np.isfinite(latitudes)
        & (latitudes >= LAT_MIN)
        & (latitudes <= LAT_MAX)
    )

    for sea_code, (lon_min, lon_max) in SEA_BOUNDS.items():
        mask = (
            in_latitude_range
            & longitude_in_sector(longitudes, lon_min, lon_max)
        )

        output[mask] = sea_code

    return output


def power_of_ten_ceiling(value):
    if not np.isfinite(value) or value <= 1.0:
        return 10.0

    return float(10.0 ** np.ceil(np.log10(value)))


def logarithmic_ticks(vmax):
    maximum_power = int(np.ceil(np.log10(max(vmax, 10.0))))
    return [10.0 ** power for power in range(0, maximum_power + 1)]


def format_log_tick(value, _position=None):
    if value >= 1000:
        return f"{int(round(value)):,}"
    return f"{int(round(value))}"


def circular_map_boundary():
    theta = np.linspace(0.0, 2.0 * np.pi, 361)

    vertices = (
        np.vstack([np.sin(theta), np.cos(theta)]).T
        * 0.5
        + np.array([0.5, 0.5])
    )

    return MplPath(vertices)


def add_manual_polar_labels(axis):
    # Six longitude labels keep four small maps readable and non-overlapping.
    longitude_labels = [
        (0.0, "0°"),
        (60.0, "60°E"),
        (120.0, "120°E"),
        (180.0, "180°"),
        (-120.0, "120°W"),
        (-60.0, "60°W"),
    ]

    for longitude, text in longitude_labels:
        axis.text(
            longitude,
            -58.7,
            text,
            transform=ccrs.PlateCarree(),
            ha="center",
            va="center",
            fontsize=FONT_MAP_TICK,
            fontweight="bold",
            clip_on=False,
            zorder=20,
        )

    latitude_labels = [
        (-60.0, "60°S"),
        (-70.0, "70°S"),
        (-80.0, "80°S"),
    ]

    for latitude, text in latitude_labels:
        axis.text(
            0.0,
            latitude,
            text,
            transform=ccrs.PlateCarree(),
            ha="center",
            va="bottom",
            fontsize=FONT_MAP_TICK,
            fontweight="bold",
            zorder=20,
        )

    axis.text(
        0.5,
        0.5,
        "90°S",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=FONT_MAP_TICK,
        fontweight="bold",
        zorder=20,
    )


def configure_polar_axis(axis):
    axis.set_extent(
        [LON_MIN, LON_MAX, LAT_MIN, LAT_MAX],
        crs=ccrs.PlateCarree(),
    )

    axis.set_boundary(
        circular_map_boundary(),
        transform=axis.transAxes,
    )

    axis.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        xlocs=np.arange(-180.0, 181.0, 30.0),
        ylocs=[-60.0, -70.0, -80.0],
        linewidth=0.42,
        color=GRID_COLOR,
        alpha=0.62,
        linestyle=":",
        zorder=2,
    )

    # Cartopy Natural Earth land is used only as a visual background.
    axis.add_feature(
        cfeature.LAND,
        facecolor=LAND_COLOR,
        edgecolor=COAST_COLOR,
        linewidth=0.42,
        zorder=8,
    )

    axis.coastlines(
        resolution="110m",
        color=COAST_COLOR,
        linewidth=0.42,
        zorder=9,
    )

    add_manual_polar_labels(axis)


def add_panel_title(axis, title):
    axis.text(
        0.01,
        1.06,
        title,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=FONT_PANEL,
        fontweight="bold",
        clip_on=False,
    )


# 5. PROFILE DEPTH CALCULATION USING PRES_GRID

def calculate_profile_depth_support(dataset):
    """
    Return three one-dimensional arrays with length N_PROF:
      maximum_common_ts_pressure
      reaches_1000_dbar
      reaches_2000_dbar

    A pressure level is accepted only when BOTH TEMP and PSAL are finite.
    Processing is chunked along N_PROF to avoid loading both full 132 MB
    variables into memory simultaneously.
    """
    required_variables = ["PRES_GRID", "TEMP", "PSAL"]

    missing = [
        variable_name
        for variable_name in required_variables
        if variable_name not in dataset.variables
    ]

    if missing:
        raise KeyError(
            "The following required variables are missing: "
            + ", ".join(missing)
        )

    pressure_grid = np.asarray(
        dataset["PRES_GRID"].values,
        dtype=np.float32,
    ).reshape(-1)

    if pressure_grid.size != dataset.sizes["N_LEVELS"]:
        raise ValueError(
            "PRES_GRID length does not match the N_LEVELS dimension."
        )

    if not np.isfinite(pressure_grid).any():
        raise ValueError("PRES_GRID contains no finite values.")

    number_of_profiles = int(dataset.sizes["N_PROF"])

    maximum_pressure = np.full(
        number_of_profiles,
        np.nan,
        dtype=np.float32,
    )

    reaches_1000 = np.zeros(number_of_profiles, dtype=bool)
    reaches_2000 = np.zeros(number_of_profiles, dtype=bool)

    pressure_matrix = pressure_grid[None, :]
    threshold_1000 = pressure_grid >= 1000.0
    threshold_2000 = pressure_grid >= 2000.0

    print("\nCalculating profile depth support from PRES_GRID + valid TEMP/PSAL...")
    print(f"Profiles       : {number_of_profiles:,}")
    print(f"Pressure levels: {pressure_grid.size}")
    print(
        f"Pressure range : {np.nanmin(pressure_grid):.2f} to "
        f"{np.nanmax(pressure_grid):.2f} dbar"
    )
    print(f"Chunk size     : {PROFILE_CHUNK_SIZE:,} profiles")

    for start in range(0, number_of_profiles, PROFILE_CHUNK_SIZE):
        stop = min(start + PROFILE_CHUNK_SIZE, number_of_profiles)

        temperature = np.asarray(
            dataset["TEMP"].isel(N_PROF=slice(start, stop)).values,
            dtype=np.float32,
        )

        salinity = np.asarray(
            dataset["PSAL"].isel(N_PROF=slice(start, stop)).values,
            dtype=np.float32,
        )

        if temperature.shape != salinity.shape:
            raise ValueError(
                "TEMP and PSAL chunk shapes do not match: "
                f"{temperature.shape} versus {salinity.shape}"
            )

        valid_common_ts = (
            np.isfinite(temperature)
            & np.isfinite(salinity)
        )

        deepest = np.max(
            np.where(valid_common_ts, pressure_matrix, -np.inf),
            axis=1,
        )

        deepest[~np.isfinite(deepest)] = np.nan

        maximum_pressure[start:stop] = deepest.astype(np.float32)

        if threshold_1000.any():
            reaches_1000[start:stop] = np.any(
                valid_common_ts[:, threshold_1000],
                axis=1,
            )

        if threshold_2000.any():
            reaches_2000[start:stop] = np.any(
                valid_common_ts[:, threshold_2000],
                axis=1,
            )

        del temperature, salinity, valid_common_ts, deepest
        gc.collect()

        print(
            f"  [{stop:>7,}/{number_of_profiles:,}] "
            f"profiles processed"
        )

    return maximum_pressure, reaches_1000, reaches_2000


# 6. BUILD PROFILE-LEVEL SUMMARY

def build_profile_summary(dataset):
    required_variables = ["JULD", "LATITUDE", "LONGITUDE"]

    missing = [
        variable_name
        for variable_name in required_variables
        if variable_name not in dataset.variables
    ]

    if missing:
        raise KeyError(
            "Required profile coordinates are missing: "
            + ", ".join(missing)
        )

    number_of_profiles = int(dataset.sizes["N_PROF"])

    latitude = np.asarray(
        dataset["LATITUDE"].values,
        dtype=np.float32,
    ).reshape(-1)

    longitude = normalize_longitude(
        dataset["LONGITUDE"].values
    ).astype(np.float32)

    dates = convert_juld_to_datetime(dataset["JULD"])

    if not (
        len(latitude)
        == len(longitude)
        == len(dates)
        == number_of_profiles
    ):
        raise ValueError(
            "Profile coordinate lengths do not match N_PROF."
        )

    maximum_pressure, reaches_1000, reaches_2000 = (
        calculate_profile_depth_support(dataset)
    )

    summary = pd.DataFrame({
        "profile_index": np.arange(number_of_profiles, dtype=np.int64),
        "time": dates,
        "latitude": latitude,
        "longitude": longitude,
        "maximum_common_ts_pressure_dbar": maximum_pressure,
        "reaches_1000_dbar": reaches_1000,
        "reaches_2000_dbar": reaches_2000,
    })

    summary = summary.replace([np.inf, -np.inf], np.nan)

    summary = summary.dropna(
        subset=["time", "latitude", "longitude"]
    ).copy()

    summary = summary[
        (summary["latitude"] >= LAT_MIN)
        & (summary["latitude"] <= LAT_MAX)
    ].copy()

    summary["month"] = summary["time"].dt.month.astype(int)
    summary["season"] = summary["month"].map(month_to_season)

    summary["sea"] = assign_sea_codes(
        summary["longitude"].to_numpy(),
        summary["latitude"].to_numpy(),
    )

    summary = summary.dropna(
        subset=["season", "sea"]
    ).copy()

    summary["sea"] = pd.Categorical(
        summary["sea"],
        categories=SEA_ORDER,
        ordered=True,
    )

    summary["season"] = pd.Categorical(
        summary["season"],
        categories=list(SEASONS.keys()),
        ordered=True,
    )

    summary = summary.sort_values(
        ["time", "profile_index"]
    ).reset_index(drop=True)

    print("\nProfile summary created:")
    print(f"Retained profiles : {len(summary):,}")
    print(f"First date        : {summary['time'].min()}")
    print(f"Last date         : {summary['time'].max()}")
    print(
        "Profiles with valid common T–S depth: "
        f"{summary['maximum_common_ts_pressure_dbar'].notna().sum():,}"
    )

    return summary


# 7. AGGREGATED STATISTICS
def create_seasonal_profile_count_maps(profile_summary):
    longitude_edges = np.arange(
        LON_MIN,
        LON_MAX + MAP_GRID_RESOLUTION_DEG,
        MAP_GRID_RESOLUTION_DEG,
    )

    latitude_edges = np.arange(
        LAT_MIN,
        LAT_MAX + MAP_GRID_RESOLUTION_DEG,
        MAP_GRID_RESOLUTION_DEG,
    )

    longitude_centres = (
        longitude_edges[:-1]
        + 0.5 * MAP_GRID_RESOLUTION_DEG
    )

    latitude_centres = (
        latitude_edges[:-1]
        + 0.5 * MAP_GRID_RESOLUTION_DEG
    )

    data_arrays = []

    for season_name in SEASONS.keys():
        subset = profile_summary[
            profile_summary["season"] == season_name
        ]

        counts, _, _ = np.histogram2d(
            subset["latitude"].to_numpy(dtype=float),
            subset["longitude"].to_numpy(dtype=float),
            bins=[latitude_edges, longitude_edges],
        )

        data_arrays.append(
            xr.DataArray(
                counts.astype(np.int32),
                coords={
                    "lat": latitude_centres,
                    "lon": longitude_centres,
                },
                dims=("lat", "lon"),
                name=season_name,
            )
        )

    seasonal = xr.concat(
        data_arrays,
        dim=pd.Index(list(SEASONS.keys()), name="season"),
    )

    seasonal.name = "profile_count"

    seasonal.attrs.update({
        "long_name": "Seasonal Argo profile count per map grid cell",
        "grid_resolution_degrees": MAP_GRID_RESOLUTION_DEG,
        "profile_definition": "One count per Argo profile",
        "period": "2001-2025",
    })

    return seasonal


def create_sea_season_count_table(profile_summary):
    table = pd.crosstab(
        profile_summary["sea"],
        profile_summary["season"],
        dropna=False,
    )

    table = table.reindex(
        index=SEA_ORDER,
        columns=list(SEASONS.keys()),
        fill_value=0,
    )

    return table.astype(int)


def create_depth_threshold_table(profile_summary):
    rows = []

    for sea_code in SEA_ORDER:
        subset = profile_summary[
            profile_summary["sea"] == sea_code
        ]

        number_of_profiles = len(subset)

        if number_of_profiles == 0:
            percentage_1000 = np.nan
            percentage_2000 = np.nan
            count_1000 = 0
            count_2000 = 0
        else:
            count_1000 = int(subset["reaches_1000_dbar"].sum())
            count_2000 = int(subset["reaches_2000_dbar"].sum())

            percentage_1000 = (
                100.0 * count_1000 / number_of_profiles
            )

            percentage_2000 = (
                100.0 * count_2000 / number_of_profiles
            )

        rows.append({
            "sea": sea_code,
            "total_profiles": number_of_profiles,
            "profiles_reaching_1000_dbar": count_1000,
            "percentage_reaching_1000_dbar": percentage_1000,
            "profiles_reaching_2000_dbar": count_2000,
            "percentage_reaching_2000_dbar": percentage_2000,
        })

    return pd.DataFrame(rows).set_index("sea")


# 8. PLOTTING HELPERS
def draw_seasonal_map(
    figure,
    map_axis,
    colorbar_axis,
    seasonal_counts,
    season_name,
    panel_letter,
    map_norm,
    map_ticks,
):
    configure_polar_axis(map_axis)

    field = seasonal_counts.sel(
        season=season_name
    ).values.astype(float)

    field[field <= 0.0] = np.nan

    longitude_2d, latitude_2d = np.meshgrid(
        seasonal_counts["lon"].values,
        seasonal_counts["lat"].values,
    )

    image = map_axis.pcolormesh(
        longitude_2d,
        latitude_2d,
        field,
        transform=ccrs.PlateCarree(),
        cmap="Blues",
        norm=map_norm,
        shading="auto",
        zorder=4,
    )

    add_panel_title(
        map_axis,
        f"({panel_letter}) {season_name}",
    )

    colorbar = figure.colorbar(
        image,
        cax=colorbar_axis,
        orientation="horizontal",
    )

    colorbar.set_ticks(map_ticks)
    colorbar.ax.xaxis.set_major_formatter(
        FuncFormatter(format_log_tick)
    )

    colorbar.ax.tick_params(
        labelsize=FONT_TICK,
        length=2.5,
        pad=1.0,
    )

    for label in colorbar.ax.get_xticklabels():
        label.set_fontweight("bold")

    colorbar.set_label(
        "Profile count",
        fontsize=FONT_AXIS,
        fontweight="bold",
        labelpad=2.0,
    )


def draw_sea_season_heatmap(
    figure,
    axis,
    colorbar_axis,
    count_table,
):
    values = count_table.to_numpy(dtype=float)
    plot_values = values.copy()
    plot_values[plot_values <= 0.0] = np.nan

    positive = plot_values[np.isfinite(plot_values)]
    vmax = power_of_ten_ceiling(np.nanmax(positive)) if positive.size else 10.0

    norm = LogNorm(vmin=1.0, vmax=vmax)

    image = axis.imshow(
        plot_values,
        cmap="Blues",
        norm=norm,
        aspect="auto",
        interpolation="nearest",
    )

    axis.set_xticks(np.arange(len(SEASONS)))
    axis.set_xticklabels(
        list(SEASONS.keys()),
        fontsize=FONT_TICK,
        fontweight="bold",
    )

    axis.xaxis.tick_top()
    axis.tick_params(
        axis="x",
        top=True,
        labeltop=True,
        bottom=False,
        labelbottom=False,
        length=0,
        pad=2.5,
    )

    axis.set_yticks(np.arange(len(SEA_ORDER)))
    axis.set_yticklabels(
        SEA_ORDER,
        fontsize=FONT_TICK,
        fontweight="bold",
    )

    axis.tick_params(axis="y", length=0, pad=2.0)

    axis.set_xticks(
        np.arange(-0.5, len(SEASONS), 1),
        minor=True,
    )

    axis.set_yticks(
        np.arange(-0.5, len(SEA_ORDER), 1),
        minor=True,
    )

    axis.grid(
        which="minor",
        color="white",
        linewidth=0.85,
    )

    axis.tick_params(
        which="minor",
        bottom=False,
        left=False,
    )

    maximum_count = float(np.nanmax(values)) if np.isfinite(values).any() else 1.0

    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]

            if not np.isfinite(value):
                continue

            text_color = (
                "white"
                if value >= 0.42 * maximum_count
                else "black"
            )

            axis.text(
                column,
                row,
                f"{int(round(value)):,}",
                ha="center",
                va="center",
                fontsize=FONT_CELL,
                fontweight="bold",
                color=text_color,
            )

    axis.set_title(
        "(e) Sea × season profile-count heatmap",
        loc="left",
        fontsize=FONT_PANEL,
        fontweight="bold",
        pad=9.0,
    )

    ticks = logarithmic_ticks(vmax)

    colorbar = figure.colorbar(
        image,
        cax=colorbar_axis,
        orientation="vertical",
    )

    colorbar.set_ticks(ticks)
    colorbar.ax.yaxis.set_major_formatter(
        FuncFormatter(format_log_tick)
    )

    colorbar.ax.tick_params(
        axis="y",
        labelsize=FONT_TICK,
        length=2.5,
        pad=1.0,
        labelright=True,
        labelleft=False,
    )

    for label in colorbar.ax.get_yticklabels():
        label.set_fontweight("bold")

    # Keep the numerical tick labels on the right side of the colorbar.
    colorbar.ax.yaxis.set_ticks_position("right")

    # Move only the vertical colorbar label to the left so that it cannot
    # overlap panel (f).
    colorbar.ax.yaxis.set_label_position("left")

    colorbar.set_label(
        "Profile count",
        fontsize=FONT_AXIS,
        fontweight="bold",
        rotation=90,
        labelpad=7.0,
    )


def draw_depth_threshold_bars(axis, depth_table):
    values_1000 = depth_table.loc[
        SEA_ORDER,
        "percentage_reaching_1000_dbar",
    ].to_numpy(dtype=float)

    values_2000 = depth_table.loc[
        SEA_ORDER,
        "percentage_reaching_2000_dbar",
    ].to_numpy(dtype=float)

    positions = np.arange(len(SEA_ORDER))
    bar_height = 0.34

    axis.barh(
        positions - 0.5 * bar_height,
        values_1000,
        height=bar_height,
        color="#555555",
        edgecolor="#333333",
        linewidth=0.6,
        label="% reaching 1000 dbar",
        zorder=3,
    )

    axis.barh(
        positions + 0.5 * bar_height,
        values_2000,
        height=bar_height,
        color="#eeeeee",
        edgecolor="#8c8c8c",
        linewidth=0.6,
        hatch="////",
        label="% reaching 2000 dbar",
        zorder=3,
    )

    axis.set_yticks(positions)
    axis.set_yticklabels(
        SEA_ORDER,
        fontsize=FONT_TICK,
        fontweight="bold",
    )

    axis.invert_yaxis()
    axis.set_xlim(0.0, 108.0)
    axis.set_xticks(np.arange(0.0, 101.0, 20.0))

    axis.set_xlabel(
        "Percentage (%)",
        fontsize=FONT_AXIS,
        fontweight="bold",
        labelpad=3.0,
    )

    axis.tick_params(
        axis="x",
        labelsize=FONT_TICK,
        pad=2.0,
    )

    for label in axis.get_xticklabels():
        label.set_fontweight("bold")

    axis.tick_params(axis="y", length=0, pad=2.0)

    axis.grid(
        axis="x",
        linestyle=":",
        linewidth=0.55,
        color="0.70",
        alpha=0.80,
        zorder=0,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    for index, value in enumerate(values_1000):
        if np.isfinite(value):
            axis.text(
                min(value + 1.4, 104.5),
                positions[index] - 0.5 * bar_height,
                f"{value:.1f}",
                ha="left",
                va="center",
                fontsize=FONT_CELL,
                fontweight="bold",
            )

    for index, value in enumerate(values_2000):
        if np.isfinite(value):
            axis.text(
                min(value + 1.4, 104.5),
                positions[index] + 0.5 * bar_height,
                f"{value:.1f}",
                ha="left",
                va="center",
                fontsize=FONT_CELL,
                fontweight="bold",
            )

    legend = axis.legend(
        loc="lower right",
        fontsize=FONT_LEGEND,
        frameon=True,
        framealpha=0.95,
        borderpad=0.45,
        labelspacing=0.55,
        handlelength=1.4,
    )

    for label in legend.get_texts():
        label.set_fontweight("bold")



# 9. MAIN PROCESSING
print("=" * 78)
print("LOADING ARGO DATA")
print("=" * 78)

if not ARGO_FILE.exists():
    raise FileNotFoundError(f"Argo file not found:\n{ARGO_FILE}")

argo_dataset = open_dataset_safely(ARGO_FILE)
print(argo_dataset)

required_dimensions = ["N_PROF", "N_LEVELS"]
missing_dimensions = [
    dimension
    for dimension in required_dimensions
    if dimension not in argo_dataset.dims
]

if missing_dimensions:
    raise ValueError(
        "The input file does not have the required dimensions: "
        + ", ".join(missing_dimensions)
    )

profile_summary = build_profile_summary(argo_dataset)

# Close the large NetCDF before plotting and aggregation.
argo_dataset.close()
del argo_dataset
gc.collect()

seasonal_counts = create_seasonal_profile_count_maps(profile_summary)
sea_season_counts = create_sea_season_count_table(profile_summary)
depth_thresholds = create_depth_threshold_table(profile_summary)

print("\nSea × season profile counts:")
print(sea_season_counts)

print("\nDepth-threshold percentages:")
print(depth_thresholds)

# Save supporting outputs before plotting.
sea_season_counts.to_csv(OUTPUT_HEATMAP_CSV)
depth_thresholds.to_csv(OUTPUT_DEPTH_CSV)

seasonal_counts.to_dataset().to_netcdf(
    OUTPUT_SEASONAL_NC,
    encoding={
        "profile_count": {
            "zlib": True,
            "complevel": 4,
            "dtype": "int32",
        }
    },
)

if SAVE_PROFILE_LEVEL_SUMMARY:
    profile_summary.to_csv(
        OUTPUT_PROFILE_SUMMARY,
        index=False,
        compression="gzip",
    )


# 10. SHARED MAP SCALE
all_map_values = seasonal_counts.values.astype(float)
positive_map_values = all_map_values[
    np.isfinite(all_map_values)
    & (all_map_values > 0.0)
]

if positive_map_values.size:
    map_vmax = power_of_ten_ceiling(
        np.nanmax(positive_map_values)
    )
else:
    map_vmax = 10.0

map_norm = LogNorm(vmin=1.0, vmax=map_vmax)
map_ticks = logarithmic_ticks(map_vmax)


# 11. FIGURE LAYOUT — NO OVERALL TITLE
figure = plt.figure(figsize=FIGSIZE)

main_grid = figure.add_gridspec(
    2,
    4,
    left=0.045,
    right=0.985,
    top=0.955,
    bottom=0.060,
    height_ratios=[1.04, 0.96],
    wspace=0.19,
    hspace=0.32,
)

map_axes = []
map_colorbar_axes = []

for column in range(4):
    map_grid = main_grid[0, column].subgridspec(
        2,
        1,
        height_ratios=[19.0, 1.05],
        hspace=0.08,
    )

    map_axes.append(
        figure.add_subplot(
            map_grid[0, 0],
            projection=ccrs.SouthPolarStereo(),
        )
    )

    map_colorbar_axes.append(
        figure.add_subplot(map_grid[1, 0])
    )

# Bottom row uses an explicit spacer between panel (e)'s vertical colorbar
# and panel (f)'s y-axis labels. This prevents colorbar ticks from touching
# or overlapping WED/KHV/... labels in panel (f).
bottom_grid = main_grid[1, :].subgridspec(
    1,
    4,
    width_ratios=[20.5, 1.0, 2.8, 25.0],
    wspace=0.05,
)

heatmap_axis = figure.add_subplot(bottom_grid[0, 0])
heatmap_colorbar_axis = figure.add_subplot(bottom_grid[0, 1])

# Dedicated blank spacer column.
spacer_axis = figure.add_subplot(bottom_grid[0, 2])
spacer_axis.set_axis_off()

# Panel (f) receives its own title row so the complete title cannot be clipped.
depth_grid = bottom_grid[0, 3].subgridspec(
    2,
    1,
    height_ratios=[2.0, 18.0],
    hspace=0.04,
)

depth_title_axis = figure.add_subplot(depth_grid[0, 0])
depth_title_axis.set_axis_off()
depth_title_axis.text(
    0.5,
    0.54,
    "(f) Percentage of profiles reaching\n1000 and 2000 dbar",
    ha="center",
    va="center",
    fontsize=FONT_PANEL,
    fontweight="bold",
    linespacing=1.0,
    clip_on=False,
)

depth_axis = figure.add_subplot(depth_grid[1, 0])

panel_letters = ["a", "b", "c", "d"]

for axis, colorbar_axis, season_name, panel_letter in zip(
    map_axes,
    map_colorbar_axes,
    SEASONS.keys(),
    panel_letters,
):
    draw_seasonal_map(
        figure,
        axis,
        colorbar_axis,
        seasonal_counts,
        season_name,
        panel_letter,
        map_norm,
        map_ticks,
    )


draw_sea_season_heatmap(
    figure,
    heatmap_axis,
    heatmap_colorbar_axis,
    sea_season_counts,
)


draw_depth_threshold_bars(
    depth_axis,
    depth_thresholds,
)




# 12. SAVE AT 1080 DPI
print("\nSaving vector PDF...")
figure.savefig(
    OUTPUT_PDF,
    facecolor="white",
    edgecolor="none",
    bbox_inches=None,
)

gc.collect()

print("Saving fixed-canvas PNG at 1080 dpi...")
figure.savefig(
    OUTPUT_PNG,
    dpi=SAVE_DPI,
    facecolor="white",
    edgecolor="none",
    bbox_inches=None,
    pil_kwargs={"compress_level": 6},
)

plt.close(figure)
gc.collect()


# 13. PROCESSING LOG
log_lines = [
    "=" * 92,
    "SUPPLEMENTARY FIGURE S1 — ARGO COVERAGE",
    "=" * 92,
    "",
    f"Input file: {ARGO_FILE}",
    f"Profiles in original file: {len(profile_summary):,} retained after filtering",
    f"Date range: {profile_summary['time'].min()} to {profile_summary['time'].max()}",
    f"Map grid resolution: {MAP_GRID_RESOLUTION_DEG} degree",
    "",
    "Layout corrections:",
    "Panel-e colorbar separated from panel-f y-axis labels by a dedicated spacer.",
    "Panel-f title placed on a dedicated two-line title axis.",
    "No bottom explanatory text is plotted.",
    "",
    "Depth rule:",
    "Maximum depth is the deepest PRES_GRID level where TEMP and PSAL are both finite.",
    "A profile reaches 1000/2000 dbar when at least one common valid T-S sample exists at or below that pressure.",
    "",
    "Outputs:",
    f"PNG: {OUTPUT_PNG}",
    f"PDF: {OUTPUT_PDF}",
    f"Sea-season CSV: {OUTPUT_HEATMAP_CSV}",
    f"Depth-threshold CSV: {OUTPUT_DEPTH_CSV}",
    f"Seasonal-map NetCDF: {OUTPUT_SEASONAL_NC}",
]

if SAVE_PROFILE_LEVEL_SUMMARY:
    log_lines.append(
        f"Profile summary: {OUTPUT_PROFILE_SUMMARY}"
    )

OUTPUT_LOG.write_text(
    "\n".join(log_lines),
    encoding="utf-8",
)

print("\n".join(log_lines))
print("\nFigure S1 completed successfully without an overall title or bottom note.")
