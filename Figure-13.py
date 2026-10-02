"""
FIGURE 7 — FULL ORGANIZED LATEST WORKFLOW

FIGURE

Observed Sea Level Anomaly (SLA) versus EN4 total steric height
Southern Ocean, 2008-2025.

SOURCE

Extracted from the latest Figure 7 implementation in:
    so_sealevel_paper_fig.py

WHY THIS VERSION

The source file contains an earlier Figure 7 based on Argo objective mapping
and a later Figure 7 based on EN4 total steric height on the native EN4 grid.

This standalone script uses ONLY the later/latest Figure 7 methodology:

    - observed satellite SLA in the left column;
    - EN4 total steric-height anomaly integrated over 0-1000 m in the right;
    - both fields represented on the native EN4 1-degree grid;
    - one fixed common mask for direct comparison;
    - no Argo objective mapping;
    - no Argo grid-cell dots;
    - Argo-EN4 validation left for supplementary analysis.

ORGANIZED EXECUTION ORDER
1. Mount Google Drive
2. Install missing Python packages
3. Import all libraries
4. Define input/output paths
5. Define analysis and figure settings
6. Define dataset/coordinate/statistical helper functions
7. Verify all required input files
8. Open and identify SLA, EN4, SIC and GEBCO variables
9. Standardize coordinates and longitude
10. Prepare science variables
11. Align the complete 2008-2025 monthly period
12. Create the OSTIA ocean mask
13. Create the GEBCO >=1000 m depth mask on the EN4 grid
14. Build the fixed SLA/EN4 common mask
15. Calculate seasonal climatological anomalies
16. Calculate seasonal spatial statistics
17. Save fixed/common masks to NetCDF
18. Save seasonal SLA and EN4 anomaly fields to NetCDF
19. Determine common plotting limits and no-data masks
20. Prepare the bathymetric contour
21. Draw the final four-season x two-variable polar figure
22. Save final 1080-dpi PNG and vector PDF
23. Print the complete output report and close all datasets

INPUT FILES

1. SLA:
   /content/drive/MyDrive/SAM_Thesis/Data/
   SLA_Antarctic_monthly_2008_2025.nc

2. EN4 total steric height:
   /content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
   EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc

3. OSTIA sea-ice fraction:
   /content/drive/MyDrive/SAM_Thesis/Data/
   OSTIA_sea_ice_fraction_monthly_2008_2025_SO.nc

4. GEBCO bathymetry:
   GEBCO_2024_CEC.nc or GEBCO_2024_CF.nc

FINAL OUTPUTS

/content/drive/MyDrive/SAM_Thesis/paper2/

    Figure07_SLA_vs_EN4_total_steric_1080dpi.png
    Figure07_SLA_vs_EN4_total_steric.pdf
    Fig07_SLA_EN4_total_steric_statistics.csv
    Fig07_fixed_common_mask_EN4_grid.nc
    Fig07_SLA_EN4_seasonal_anomalies_common_mask.nc

SCIENTIFIC DESIGN PRESERVED

The latest source methodology is retained:
    - study period 2008-2025;
    - domain south of 60 degrees S;
    - SLA minimum valid fraction = 0.70;
    - EN4 minimum valid fraction = 1.00;
    - minimum water depth = 1000 m;
    - Spring SON, Summer DJF, Autumn MAM, Winter JJA;
    - seasonal anomaly = seasonal climatology - all-month mean;
    - seasonal spatial Pearson correlation, RMSD, bias and coverage;
    - bias = EN4 total steric - observed SLA.

CORRECTION APPLIED
One standard longitude-normalization expression in the supplied latest Figure 7
code had `% 360.0` accidentally commented out. It is restored so longitude is
correctly converted with:

    ((longitude + 180.0) % 360.0) - 180.0

No other scientific thresholds, file names, analysis definitions or output
products are intentionally changed.

"""



# ORIGINAL LATEST FIGURE 7 SCIENTIFIC DESCRIPTION

# FIGURE 7 — OBSERVED SLA VERSUS EN4 TOTAL STERIC HEIGHT
# Southern Ocean, 2008–2025
#
# Scientific design
# • Eight panels arranged as four seasons × two variables.
# • Left column  : satellite-observed SLA seasonal anomaly.
# • Right column : EN4 total steric-height seasonal anomaly,
#                  integrated over 0–1000 m.
# • Both columns are placed on the native EN4 1° grid.
# • The same fixed common mask is used in both columns.
# • No Argo objective mapping and no Argo grid-cell dots.
# • Argo–EN4 validation belongs in supplementary material.
#
# Seasonal anomalies
# SLA'(season) = SLA seasonal climatology - SLA all-month mean
#
# EN4 steric'(season) =
#     EN4 seasonal climatology - EN4 all-month mean
#
# Fixed common mask
# A cell is eligible only when:
#   1. OSTIA identifies it as an ocean cell.
#   2. GEBCO water depth is at least 1000 m.
#   3. SLA is available during at least 70% of the 216 months.
#   4. EN4 valid_layer_mask is valid during every month.
#
# For each seasonal row, an exact finite intersection is then
# applied to both maps. Therefore the two maps in a row always
# display precisely the same cells.
#
# Statistics reported for every season
# • Spatial Pearson correlation
# • RMSD
# • Mean bias = EN4 total steric - observed SLA
# • Number of common cells
# • Spatial coverage = seasonal common cells / fixed-mask cells × 100
#
# Statistics are calculated on the common native EN4 grid.
#
# Output directory
# /content/drive/MyDrive/SAM_Thesis/paper2/




# 1. MOUNT GOOGLE DRIVE

try:
    from google.colab import drive
    drive.mount("/content/drive")
except ImportError:
    print("Google Drive mounting skipped because this is not Colab.")


# 2. INSTALL MISSING PACKAGES

import sys
import subprocess
import importlib.util

REQUIRED_PACKAGES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "xarray": "xarray",
    "netCDF4": "netCDF4",
    "h5netcdf": "h5netcdf",
    "cftime": "cftime",
    "dask": "dask[array]",
    "matplotlib": "matplotlib",
    "cartopy": "cartopy",
}

missing_packages = [
    pip_name
    for import_name, pip_name in REQUIRED_PACKAGES.items()
    if importlib.util.find_spec(import_name) is None
]

if missing_packages:
    print("Installing missing packages:", ", ".join(missing_packages))
    subprocess.check_call(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "-q",
            *missing_packages,
        ]
    )
    print("Package installation completed.")
else:
    print("All required packages are already installed.")


# 3. IMPORT LIBRARIES

import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

warnings.filterwarnings("ignore", category=RuntimeWarning)


# 4. COMPLETE LATEST FIGURE 7 PROCESSING AND PLOTTING WORKFLOW

# 3. INPUT AND OUTPUT PATHS
SLA_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Data/"
    "SLA_Antarctic_monthly_2008_2025.nc"
)

EN4_STERIC_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Processed/"
    "EN4_NetCDF_inventory/"
    "EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc"
)

SIC_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Data/"
    "OSTIA_sea_ice_fraction_monthly_2008_2025_SO.nc"
)

GEBCO_CANDIDATES = [
    (
        "/content/drive/MyDrive/SAM_Thesis/Data/"
        "GEBCO_2024_CEC.nc"
    ),
    (
        "/content/drive/MyDrive/SAM_Thesis/Data/"
        "GEBCO_2024_CF.nc"
    ),
]

GEBCO_FILE = next(
    (
        candidate
        for candidate in GEBCO_CANDIDATES
        if os.path.exists(candidate)
    ),
    GEBCO_CANDIDATES[0],
)

OUTPUT_DIRECTORY = Path(
    "/content/drive/MyDrive/SAM_Thesis/paper2"
)
OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

OUT_PNG = (
    OUTPUT_DIRECTORY
    / "Figure07_SLA_vs_EN4_total_steric_1080dpi.png"
)

OUT_PDF = (
    OUTPUT_DIRECTORY
    / "Figure07_SLA_vs_EN4_total_steric.pdf"
)

OUT_STATS_CSV = (
    OUTPUT_DIRECTORY
    / "Fig07_SLA_EN4_total_steric_statistics.csv"
)

OUT_MASK_NC = (
    OUTPUT_DIRECTORY
    / "Fig07_fixed_common_mask_EN4_grid.nc"
)

OUT_SEASONAL_NC = (
    OUTPUT_DIRECTORY
    / "Fig07_SLA_EN4_seasonal_anomalies_common_mask.nc"
)


# 4. ANALYSIS SETTINGS
LAT_MIN = -90.0
LAT_MAX = -60.0

YEAR_START = 2008
YEAR_END = 2025

MINIMUM_SLA_VALID_FRACTION = 0.70

# EN4 must be valid in all months because the product was generated
# using a complete-layer 0–1000 m requirement.
MINIMUM_EN4_VALID_FRACTION = 1.00

MINIMUM_WATER_DEPTH_M = 1000.0

SEASON_ORDER = [
    "Spring",
    "Summer",
    "Autumn",
    "Winter",
]

SEASON_MONTHS = {
    "Spring": [9, 10, 11],
    "Summer": [12, 1, 2],
    "Autumn": [3, 4, 5],
    "Winter": [6, 7, 8],
}

SEASON_ROW_LABEL = {
    "Spring": "Spring (SON)",
    "Summer": "Summer (DJF)",
    "Autumn": "Autumn (MAM)",
    "Winter": "Winter (JJA)",
}

MERIDIANS = np.arange(-180, 181, 30)
PARALLELS = [-60, -70, -80]

DRAW_1000_M_CONTOUR = True
BATHYMETRY_CONTOUR_LEVEL = -1000.0

# Use one shared symmetric scale across both columns for direct
# visual comparison. Set False to use separate column scales.
USE_SHARED_COLOR_LIMIT = True

SAVE_DPI = 1080


# 5. FIGURE SETTINGS
FIGSIZE = (14.8, 20.8)

FONT_PANEL_TITLE = 16.4
FONT_ROW_LABEL = 18.5
FONT_GEO_LABEL = 9.6
FONT_STATISTICS = 11.8
FONT_COLORBAR = 12.4
FONT_COLORBAR_TICK = 10.8
FONT_LEGEND = 10.8

NO_DATA_COLOR = "0.84"
LAND_COLOR = "0.78"
COAST_COLOR = "0.35"
GRID_COLOR = "0.56"

# Longitude labels are positioned just outside the circular border.
RIM_LABEL_LATITUDE = -58.25

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.weight": "bold",
    "axes.titleweight": "bold",
    "axes.labelweight": "bold",
    "axes.linewidth": 1.1,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "mathtext.default": "regular",
})


# 6. DATASET HELPERS
def safe_is_datetime(dtype):
    try:
        return np.issubdtype(dtype, np.datetime64)
    except TypeError:
        return False


def open_dataset_safely(file_path, chunks=None):
    """
    Open a NetCDF file using several xarray engines.
    """
    attempts = []

    for engine in [
        None,
        "netcdf4",
        "h5netcdf",
        "scipy",
    ]:

        for decode_times in [True, False]:

            try:

                kwargs = {
                    "decode_times": decode_times,
                    "mask_and_scale": True,
                }

                if engine is not None:
                    kwargs["engine"] = engine

                if chunks is not None:
                    kwargs["chunks"] = chunks

                dataset = xr.open_dataset(
                    file_path,
                    **kwargs,
                )

                engine_name = (
                    "xarray-default"
                    if engine is None
                    else engine
                )

                print(
                    f"Opened {os.path.basename(file_path)} "
                    f"with engine={engine_name}, "
                    f"decode_times={decode_times}"
                )

                return dataset

            except Exception as error:

                attempts.append(
                    f"engine={engine}, "
                    f"decode_times={decode_times}: {error}"
                )

    raise RuntimeError(
        f"Could not open:\n{file_path}\n\n"
        + "\n".join(attempts)
    )


def detect_coordinate(dataset, coordinate_type):
    """
    Detect time, latitude, and longitude coordinate names.
    """
    names = list(dataset.coords)

    names += [
        name
        for name in dataset.variables
        if name not in dataset.coords
    ]

    aliases = {
        "time": [
            "time",
            "date",
            "datetime",
            "month",
            "juld",
            "t",
        ],
        "lat": [
            "lat",
            "latitude",
            "nav_lat",
            "y",
        ],
        "lon": [
            "lon",
            "longitude",
            "nav_lon",
            "x",
        ],
    }

    for name in names:

        variable = dataset[name]

        lower_name = name.lower()

        standard_name = str(
            variable.attrs.get(
                "standard_name",
                "",
            )
        ).lower()

        axis = str(
            variable.attrs.get(
                "axis",
                "",
            )
        ).upper()

        units = str(
            variable.attrs.get(
                "units",
                "",
            )
        ).lower()

        if coordinate_type == "time":

            if (
                lower_name in aliases["time"]
                or standard_name == "time"
                or axis == "T"
                or " since " in units
                or safe_is_datetime(variable.dtype)
            ):
                return name

        elif coordinate_type == "lat":

            if (
                lower_name in aliases["lat"]
                or standard_name == "latitude"
                or axis == "Y"
                or "degree_north" in units
                or "degrees_north" in units
            ):
                return name

        elif coordinate_type == "lon":

            if (
                lower_name in aliases["lon"]
                or standard_name == "longitude"
                or axis == "X"
                or "degree_east" in units
                or "degrees_east" in units
            ):
                return name

    raise KeyError(
        f"Could not identify the {coordinate_type} coordinate."
    )


def choose_data_variable(
    dataset,
    preferred_names,
    excluded_names=None,
):
    """
    Select the intended science variable.
    """
    if excluded_names is None:
        excluded_names = []

    variables = [
        name
        for name in dataset.data_vars
        if name not in excluded_names
    ]

    for preferred_name in preferred_names:

        for variable in variables:

            if preferred_name.lower() in variable.lower():
                return variable

    if not variables:
        raise ValueError(
            "No suitable data variable was found."
        )

    return variables[0]


def standardize_dataset_coordinates(
    dataset,
    time_name=None,
    lat_name=None,
    lon_name=None,
):
    """
    Rename coordinates to canonical time, lat, and lon names before
    any interpolation.
    """
    rename_mapping = {}

    if (
        time_name is not None
        and time_name != "time"
    ):
        rename_mapping[time_name] = "time"

    if (
        lat_name is not None
        and lat_name != "lat"
    ):
        rename_mapping[lat_name] = "lat"

    if (
        lon_name is not None
        and lon_name != "lon"
    ):
        rename_mapping[lon_name] = "lon"

    if rename_mapping:
        dataset = dataset.rename(rename_mapping)

    return dataset


def normalize_and_sort_longitude(dataset):
    """
    Normalize longitude to -180 ... 180 and remove duplicate endpoints.
    """
    if "lon" not in dataset.coords:
        return dataset

    normalized_lon = (
        (
            dataset["lon"].astype(float)
            + 180.0
        )
         % 360.0
    ) - 180.0

    dataset = dataset.assign_coords(
        lon=normalized_lon
    )

    longitude_values = dataset["lon"].values

    _, unique_indices = np.unique(
        longitude_values,
        return_index=True,
    )

    dataset = dataset.isel(
        lon=np.sort(unique_indices)
    )

    return dataset.sortby("lon")


def subset_antarctic_latitudes(
    dataset,
    southern_limit=LAT_MIN,
    northern_limit=LAT_MAX,
):
    """
    Subset and sort a one-dimensional latitude coordinate.
    """
    latitude = dataset["lat"].values

    if latitude[0] <= latitude[-1]:

        subset = dataset.sel(
            lat=slice(
                southern_limit,
                northern_limit,
            )
        )

    else:

        subset = dataset.sel(
            lat=slice(
                northern_limit,
                southern_limit,
            )
        )

    return subset.sortby("lat")


def restrict_time_period(data_array):
    """
    Restrict a time-dependent variable to 2008–2025.
    """
    return data_array.sel(
        time=slice(
            f"{YEAR_START}-01-01",
            f"{YEAR_END}-12-31",
        )
    )


def regrid_boolean_mask(
    boolean_mask,
    target_latitude,
    target_longitude,
):
    """
    Regrid a Boolean mask using numeric nearest-neighbour interpolation.
    """
    numeric_mask = boolean_mask.astype(np.float32)

    interpolated = numeric_mask.interp(
        lat=target_latitude,
        lon=target_longitude,
        method="nearest",
    )

    return interpolated.fillna(0.0) >= 0.5


def convert_height_to_cm(data_array):
    """
    Convert metres to centimetres when required.
    """
    units = str(
        data_array.attrs.get(
            "units",
            "",
        )
    ).strip().lower()

    if units in [
        "cm",
        "centimeter",
        "centimeters",
        "centimetre",
        "centimetres",
    ]:

        output = data_array.copy()
        output.attrs["units"] = "cm"
        return output

    if units in [
        "m",
        "meter",
        "meters",
        "metre",
        "metres",
    ]:

        output = data_array * 100.0
        output.attrs.update(data_array.attrs)
        output.attrs["units"] = "cm"
        return output

    # Conservative magnitude-based fallback.
    sample = data_array.isel(
        time=slice(
            0,
            min(
                12,
                data_array.sizes.get("time", 1),
            ),
        )
    )

    values = np.asarray(
        sample.values,
        dtype=float,
    )

    values = values[np.isfinite(values)]

    if (
        values.size > 0
        and np.nanpercentile(
            np.abs(values),
            99,
        ) < 1.0
    ):

        output = data_array * 100.0
        output.attrs.update(data_array.attrs)
        output.attrs["units"] = "cm"
        return output

    output = data_array.copy()
    output.attrs["units"] = units if units else "unknown"
    return output


def seasonal_climatology(data_array):
    """
    Calculate austral seasonal climatological means.
    """
    output = {}

    for season in SEASON_ORDER:

        selected = data_array.where(
            data_array["time"].dt.month.isin(
                SEASON_MONTHS[season]
            ),
            drop=True,
        )

        output[season] = selected.mean(
            dim="time",
            skipna=True,
        )

    return output


def symmetric_limit(fields):
    """
    Calculate a robust, rounded symmetric color limit.
    """
    collected = []

    for field in fields:

        values = np.asarray(
            field.values,
            dtype=float,
        )

        values = values[np.isfinite(values)]

        if values.size > 0:
            collected.append(values)

    if not collected:
        return 1.0

    all_values = np.concatenate(collected)

    robust_maximum = float(
        np.nanpercentile(
            np.abs(all_values),
            99,
        )
    )

    if robust_maximum <= 2.0:
        interval = 0.5

    elif robust_maximum <= 5.0:
        interval = 1.0

    elif robust_maximum <= 15.0:
        interval = 2.0

    else:
        interval = 5.0

    return float(
        interval
        * np.ceil(
            robust_maximum / interval
        )
    )


def coordinate_mesh(data_array):
    """
    Return two-dimensional longitude and latitude arrays.
    """
    return np.meshgrid(
        data_array["lon"].values,
        data_array["lat"].values,
    )


def format_longitude(longitude):
    """
    Format longitude labels for the polar-map rim.
    """
    value = int(longitude)

    if value == 0:
        return "0°"

    if abs(value) == 180:
        return "180°"

    if value < 0:
        return f"{abs(value)}°W"

    return f"{value}°E"


# 7. STATISTICAL HELPERS
def unweighted_spatial_statistics(
    observed_sla,
    en4_steric,
    fixed_mask_count,
):
    """
    Calculate requested statistics on the exact common finite cells.

    Bias sign:
        EN4 total steric minus observed SLA
    """
    observed_values = np.asarray(
        observed_sla.values,
        dtype=float,
    )

    steric_values = np.asarray(
        en4_steric.values,
        dtype=float,
    )

    valid = (
        np.isfinite(observed_values)
        & np.isfinite(steric_values)
    )

    n_common = int(valid.sum())

    if n_common < 2:

        return {
            "spatial_correlation_r": np.nan,
            "RMSD_cm": np.nan,
            "mean_bias_EN4_minus_SLA_cm": np.nan,
            "n_common_cells": n_common,
            "spatial_coverage_percent": (
                100.0 * n_common / fixed_mask_count
                if fixed_mask_count > 0
                else np.nan
            ),
        }

    observed = observed_values[valid]
    modeled = steric_values[valid]

    difference = modeled - observed

    return {
        "spatial_correlation_r": float(
            np.corrcoef(
                observed,
                modeled,
            )[0, 1]
        ),
        "RMSD_cm": float(
            np.sqrt(
                np.mean(
                    difference ** 2
                )
            )
        ),
        "mean_bias_EN4_minus_SLA_cm": float(
            np.mean(difference)
        ),
        "n_common_cells": n_common,
        "spatial_coverage_percent": float(
            100.0
            * n_common
            / fixed_mask_count
        ) if fixed_mask_count > 0 else np.nan,
    }


# 8. VERIFY INPUT FILES
for file_path in [
    SLA_FILE,
    EN4_STERIC_FILE,
    SIC_FILE,
    GEBCO_FILE,
]:

    if not os.path.exists(file_path):

        raise FileNotFoundError(
            f"Input file not found:\n{file_path}"
        )

print("\nUsing EN4 steric file:")
print(EN4_STERIC_FILE)

print("\nUsing GEBCO file:")
print(GEBCO_FILE)


# 9. OPEN INPUT DATASETS
ds_sla = open_dataset_safely(
    SLA_FILE,
    chunks="auto",
)

ds_en4 = open_dataset_safely(
    EN4_STERIC_FILE,
    chunks="auto",
)

ds_sic = open_dataset_safely(
    SIC_FILE,
    chunks="auto",
)

ds_gebco = open_dataset_safely(
    GEBCO_FILE,
    chunks="auto",
)


# 10. DETECT VARIABLES AND COORDINATES
sla_time_name = detect_coordinate(
    ds_sla,
    "time",
)
sla_lat_name = detect_coordinate(
    ds_sla,
    "lat",
)
sla_lon_name = detect_coordinate(
    ds_sla,
    "lon",
)

en4_time_name = detect_coordinate(
    ds_en4,
    "time",
)
en4_lat_name = detect_coordinate(
    ds_en4,
    "lat",
)
en4_lon_name = detect_coordinate(
    ds_en4,
    "lon",
)

sic_time_name = detect_coordinate(
    ds_sic,
    "time",
)
sic_lat_name = detect_coordinate(
    ds_sic,
    "lat",
)
sic_lon_name = detect_coordinate(
    ds_sic,
    "lon",
)

gebco_lat_name = detect_coordinate(
    ds_gebco,
    "lat",
)
gebco_lon_name = detect_coordinate(
    ds_gebco,
    "lon",
)

sla_variable_name = choose_data_variable(
    ds_sla,
    preferred_names=[
        "sla",
        "sea_level",
        "adt",
        "ssh",
    ],
)

en4_steric_variable_name = choose_data_variable(
    ds_en4,
    preferred_names=[
        "total_steric_height",
        "total_steric",
        "steric",
    ],
    excluded_names=[
        "valid_layer_mask",
        "gebco_water_depth",
    ],
)

sic_variable_name = choose_data_variable(
    ds_sic,
    preferred_names=[
        "SIF",
        "sea_ice_fraction",
        "ice_fraction",
        "sic",
    ],
)

gebco_variable_name = choose_data_variable(
    ds_gebco,
    preferred_names=[
        "elevation",
        "bathymetry",
        "depth",
        "z",
    ],
)

print("\nDetected variables:")
print(f"Observed SLA       : {sla_variable_name}")
print(f"EN4 total steric   : {en4_steric_variable_name}")
print(f"OSTIA sea ice      : {sic_variable_name}")
print(f"GEBCO bathymetry   : {gebco_variable_name}")


# 11. STANDARDIZE ALL COORDINATES
ds_sla = standardize_dataset_coordinates(
    ds_sla,
    time_name=sla_time_name,
    lat_name=sla_lat_name,
    lon_name=sla_lon_name,
)

ds_en4 = standardize_dataset_coordinates(
    ds_en4,
    time_name=en4_time_name,
    lat_name=en4_lat_name,
    lon_name=en4_lon_name,
)

ds_sic = standardize_dataset_coordinates(
    ds_sic,
    time_name=sic_time_name,
    lat_name=sic_lat_name,
    lon_name=sic_lon_name,
)

ds_gebco = standardize_dataset_coordinates(
    ds_gebco,
    lat_name=gebco_lat_name,
    lon_name=gebco_lon_name,
)

ds_sla = subset_antarctic_latitudes(
    normalize_and_sort_longitude(ds_sla)
)

ds_en4 = subset_antarctic_latitudes(
    normalize_and_sort_longitude(ds_en4)
)

ds_sic = subset_antarctic_latitudes(
    normalize_and_sort_longitude(ds_sic)
)

ds_gebco = subset_antarctic_latitudes(
    normalize_and_sort_longitude(ds_gebco),
    southern_limit=-90.0,
    northern_limit=-55.0,
)


# 12. PREPARE SCIENCE VARIABLES
sla = convert_height_to_cm(
    restrict_time_period(
        ds_sla[sla_variable_name]
    )
)

en4_steric = convert_height_to_cm(
    restrict_time_period(
        ds_en4[en4_steric_variable_name]
    )
)

sic = restrict_time_period(
    ds_sic[sic_variable_name]
)

gebco = ds_gebco[gebco_variable_name]

if "valid_layer_mask" in ds_en4.data_vars:

    en4_valid_layer_mask = restrict_time_period(
        ds_en4["valid_layer_mask"]
    )

else:

    en4_valid_layer_mask = xr.where(
        np.isfinite(en4_steric),
        1,
        0,
    )

print("\nUnits after conversion:")
print("Observed SLA     :", sla.attrs.get("units"))
print("EN4 total steric:", en4_steric.attrs.get("units"))


# 13. ALIGN THE COMMON MONTHLY PERIOD
common_times = np.intersect1d(
    pd.to_datetime(
        sla["time"].values
    ),
    pd.to_datetime(
        en4_steric["time"].values
    ),
)

if common_times.size == 0:

    raise ValueError(
        "No common SLA and EN4 steric months were found."
    )

sla = sla.sel(time=common_times)
en4_steric = en4_steric.sel(time=common_times)
en4_valid_layer_mask = en4_valid_layer_mask.sel(
    time=common_times
)

number_of_months = len(common_times)

print("\nCommon analysis period:")
print(pd.to_datetime(common_times[0]))
print("to")
print(pd.to_datetime(common_times[-1]))
print(f"Months = {number_of_months}")


# 14. OSTIA STATIC OCEAN MASK
sic_reference = sic.isel(time=0)

sic_ocean_native = xr.where(
    np.isfinite(sic_reference),
    True,
    False,
)

sic_ocean_native = (
    sic_ocean_native
    .astype(np.float32)
    .compute()
    >= 0.5
)

sic_ocean_on_sla = regrid_boolean_mask(
    boolean_mask=sic_ocean_native,
    target_latitude=sla["lat"],
    target_longitude=sla["lon"],
).compute()

sic_ocean_on_en4 = regrid_boolean_mask(
    boolean_mask=sic_ocean_native,
    target_latitude=en4_steric["lat"],
    target_longitude=en4_steric["lon"],
).compute()

print("\nOSTIA static ocean masks created.")


# 15. GEBCO DEPTH MASK ON THE EN4 GRID
gebco_on_en4 = gebco.interp(
    lat=en4_steric["lat"],
    lon=en4_steric["lon"],
    method="nearest",
).compute()

gebco_values_on_en4 = np.asarray(
    gebco_on_en4.values,
    dtype=float,
)

# GEBCO elevation is negative below sea level.
gebco_water_depth_on_en4 = xr.DataArray(
    np.where(
        gebco_values_on_en4 < 0.0,
        -gebco_values_on_en4,
        np.nan,
    ).astype(np.float32),
    dims=("lat", "lon"),
    coords={
        "lat": en4_steric["lat"],
        "lon": en4_steric["lon"],
    },
    name="gebco_water_depth",
)

deep_ocean_mask_en4 = (
    np.isfinite(gebco_water_depth_on_en4)
    & (
        gebco_water_depth_on_en4
        >= MINIMUM_WATER_DEPTH_M
    )
)


# 16. BUILD THE FIXED COMMON MASK
# SLA core validity on its native high-resolution grid.
sla_valid_fraction_native = (
    sla.notnull().sum(dim="time")
    / number_of_months
)

sla_core_mask_native = (
    sic_ocean_on_sla
    & (
        sla_valid_fraction_native
        >= MINIMUM_SLA_VALID_FRACTION
    )
)

# Transfer SLA validity to the EN4 grid.
sla_core_mask_on_en4 = regrid_boolean_mask(
    boolean_mask=sla_core_mask_native,
    target_latitude=en4_steric["lat"],
    target_longitude=en4_steric["lon"],
).compute()

# EN4 valid-layer fraction on the native EN4 grid.
en4_valid_fraction = (
    (
        en4_valid_layer_mask == 1
    )
    .sum(dim="time")
    / number_of_months
)

en4_fixed_valid_mask = (
    en4_valid_fraction
    >= (
        MINIMUM_EN4_VALID_FRACTION
        - 1.0e-10
    )
)

fixed_common_mask = (
    sic_ocean_on_en4
    & deep_ocean_mask_en4
    & sla_core_mask_on_en4
    & en4_fixed_valid_mask
)

fixed_common_mask = fixed_common_mask.compute()

fixed_common_cell_count = int(
    fixed_common_mask.sum().values
)

if fixed_common_cell_count == 0:

    raise RuntimeError(
        "The fixed SLA–EN4 common mask contains no cells."
    )

print("\nFixed common mask:")
print(
    f"SLA threshold = "
    f"{MINIMUM_SLA_VALID_FRACTION * 100:.0f}%"
)
print(
    f"EN4 threshold = "
    f"{MINIMUM_EN4_VALID_FRACTION * 100:.0f}%"
)
print(
    f"Minimum water depth = "
    f"{MINIMUM_WATER_DEPTH_M:.0f} m"
)
print(
    f"Fixed common cells = "
    f"{fixed_common_cell_count:,}"
)


# 17. CALCULATE SEASONAL ANOMALIES
# SLA anomalies are first calculated on the native SLA grid.
sla_all_month_mean_native = sla.mean(
    dim="time",
    skipna=True,
)

sla_seasonal_mean_native = seasonal_climatology(
    sla
)

# EN4 anomalies remain on the native EN4 grid.
en4_all_month_mean = en4_steric.mean(
    dim="time",
    skipna=True,
)

en4_seasonal_mean = seasonal_climatology(
    en4_steric
)

sla_anomaly_common = {}
en4_anomaly_common = {}
seasonal_common_mask = {}
statistics_rows = []

for season in SEASON_ORDER:

    sla_anomaly_native = (
        sla_seasonal_mean_native[season]
        - sla_all_month_mean_native
    ).where(
        sla_core_mask_native
    )

    # Regrid the already calculated SLA seasonal anomaly to EN4.
    sla_anomaly_en4 = sla_anomaly_native.interp(
        lat=en4_steric["lat"],
        lon=en4_steric["lon"],
        method="linear",
    )

    en4_anomaly = (
        en4_seasonal_mean[season]
        - en4_all_month_mean
    )

    # Exact finite seasonal intersection, applied to both columns.
    row_common_mask = (
        fixed_common_mask
        & np.isfinite(sla_anomaly_en4)
        & np.isfinite(en4_anomaly)
    )

    row_common_mask = row_common_mask.compute()

    sla_field = sla_anomaly_en4.where(
        row_common_mask
    ).compute()

    en4_field = en4_anomaly.where(
        row_common_mask
    ).compute()

    sla_anomaly_common[season] = sla_field
    en4_anomaly_common[season] = en4_field
    seasonal_common_mask[season] = row_common_mask

    statistics = unweighted_spatial_statistics(
        observed_sla=sla_field,
        en4_steric=en4_field,
        fixed_mask_count=fixed_common_cell_count,
    )

    statistics_rows.append(
        {
            "season": season,
            **statistics,
            "bias_definition": (
                "EN4 total steric minus observed SLA"
            ),
            "statistics_grid": (
                "Native EN4 grid"
            ),
            "fixed_common_mask_cells": (
                fixed_common_cell_count
            ),
        }
    )

statistics_table = pd.DataFrame(
    statistics_rows
)

statistics_table.to_csv(
    OUT_STATS_CSV,
    index=False,
)

print("\nSeasonal statistics:")
print(
    statistics_table.to_string(
        index=False
    )
)


# 18. SAVE THE FIXED AND SEASONAL COMMON MASKS
mask_output = xr.Dataset(
    {
        "fixed_common_mask": (
            fixed_common_mask.astype(np.int8)
        ),
        "sla_valid_fraction": (
            sla_valid_fraction_native.interp(
                lat=en4_steric["lat"],
                lon=en4_steric["lon"],
                method="linear",
            )
        ),
        "en4_valid_fraction": (
            en4_valid_fraction
        ),
        "gebco_water_depth": (
            gebco_water_depth_on_en4
        ),
    }
)

mask_output["fixed_common_mask"].attrs.update({
    "long_name": "fixed common SLA–EN4 analysis mask",
    "flag_values": np.array([0, 1], dtype=np.int8),
    "flag_meanings": "excluded included",
})

mask_output.attrs.update({
    "title": "Figure 7 fixed SLA–EN4 common mask",
    "SLA_valid_fraction_threshold": (
        MINIMUM_SLA_VALID_FRACTION
    ),
    "EN4_valid_fraction_threshold": (
        MINIMUM_EN4_VALID_FRACTION
    ),
    "minimum_GEBCO_water_depth_m": (
        MINIMUM_WATER_DEPTH_M
    ),
})

mask_output.rename(
    {
        "lat": "latitude",
        "lon": "longitude",
    }
).to_netcdf(
    OUT_MASK_NC,
    encoding={
        "fixed_common_mask": {
            "zlib": True,
            "complevel": 4,
        },
        "sla_valid_fraction": {
            "zlib": True,
            "complevel": 4,
        },
        "en4_valid_fraction": {
            "zlib": True,
            "complevel": 4,
        },
        "gebco_water_depth": {
            "zlib": True,
            "complevel": 4,
        },
    },
)


# 19. SAVE SEASONAL FIELDS USED IN THE FIGURE
season_coordinate = xr.DataArray(
    SEASON_ORDER,
    dims="season",
    name="season",
)

sla_stack = xr.concat(
    [
        sla_anomaly_common[season]
        for season in SEASON_ORDER
    ],
    dim=season_coordinate,
)

en4_stack = xr.concat(
    [
        en4_anomaly_common[season]
        for season in SEASON_ORDER
    ],
    dim=season_coordinate,
)

seasonal_mask_stack = xr.concat(
    [
        seasonal_common_mask[season].astype(np.int8)
        for season in SEASON_ORDER
    ],
    dim=season_coordinate,
)

seasonal_output = xr.Dataset(
    {
        "observed_SLA_anomaly": sla_stack,
        "EN4_total_steric_anomaly": en4_stack,
        "seasonal_common_mask": seasonal_mask_stack,
    }
)

seasonal_output[
    "observed_SLA_anomaly"
].attrs.update({
    "long_name": (
        "seasonal observed sea-level anomaly relative "
        "to the 2008-2025 all-month mean"
    ),
    "units": "cm",
})

seasonal_output[
    "EN4_total_steric_anomaly"
].attrs.update({
    "long_name": (
        "seasonal EN4 0-1000 m total steric-height anomaly "
        "relative to the 2008-2025 all-month mean"
    ),
    "units": "cm",
})

seasonal_output[
    "seasonal_common_mask"
].attrs.update({
    "long_name": (
        "season-specific exact common finite mask"
    ),
    "flag_values": np.array([0, 1], dtype=np.int8),
    "flag_meanings": "excluded included",
})

seasonal_output.attrs.update({
    "title": (
        "Figure 7 observed SLA and EN4 total steric "
        "seasonal anomalies"
    ),
    "bias_definition": (
        "EN4 total steric minus observed SLA"
    ),
})

seasonal_output.rename(
    {
        "lat": "latitude",
        "lon": "longitude",
    }
).to_netcdf(
    OUT_SEASONAL_NC,
    encoding={
        "observed_SLA_anomaly": {
            "zlib": True,
            "complevel": 4,
        },
        "EN4_total_steric_anomaly": {
            "zlib": True,
            "complevel": 4,
        },
        "seasonal_common_mask": {
            "zlib": True,
            "complevel": 4,
        },
    },
)


# 20. BUILD NO-DATA MASKS AND COLOR LIMITS
no_data_mask = {}

for season in SEASON_ORDER:

    no_data_mask[season] = xr.where(
        seasonal_common_mask[season],
        np.nan,
        1.0,
    )

if USE_SHARED_COLOR_LIMIT:

    shared_limit = symmetric_limit(
        [
            sla_anomaly_common[season]
            for season in SEASON_ORDER
        ]
        + [
            en4_anomaly_common[season]
            for season in SEASON_ORDER
        ]
    )

    sla_color_limit = shared_limit
    en4_color_limit = shared_limit

else:

    sla_color_limit = symmetric_limit(
        [
            sla_anomaly_common[season]
            for season in SEASON_ORDER
        ]
    )

    en4_color_limit = symmetric_limit(
        [
            en4_anomaly_common[season]
            for season in SEASON_ORDER
        ]
    )

print("\nColor limits:")
print(f"Observed SLA: ±{sla_color_limit:.2f} cm")
print(f"EN4 steric : ±{en4_color_limit:.2f} cm")


# 21. PREPARE GEBCO FOR THE 1000 M CONTOUR
if (
    gebco.ndim == 2
    and gebco.size > 1_200_000
):

    stride = int(
        np.ceil(
            np.sqrt(
                gebco.size
                / 1_200_000
            )
        )
    )

    gebco = gebco.isel(
        lat=slice(None, None, stride),
        lon=slice(None, None, stride),
    )

gebco = gebco.compute()

gebco_lon2d, gebco_lat2d = coordinate_mesh(
    gebco
)


# 22. MAP-DRAWING HELPERS
def add_circular_boundary(axis):
    """
    Apply a circular boundary to a polar map.
    """
    theta = np.linspace(
        0,
        2 * np.pi,
        500,
    )

    vertices = np.vstack(
        [
            np.sin(theta),
            np.cos(theta),
        ]
    ).T

    circle = mpath.Path(
        vertices * 0.5
        + np.array([0.5, 0.5])
    )

    axis.set_boundary(
        circle,
        transform=axis.transAxes,
    )


def add_geographic_labels(axis):
    """
    Draw bold longitude labels outside the circular border.
    """
    for longitude in MERIDIANS:

        axis.text(
            longitude,
            RIM_LABEL_LATITUDE,
            format_longitude(longitude),
            transform=ccrs.PlateCarree(),
            ha="center",
            va="center",
            fontsize=FONT_GEO_LABEL,
            fontweight="bold",
            color="black",
            clip_on=False,
            zorder=30,
        )

    axis.text(
        0.0,
        -62.0,
        "60°S",
        transform=ccrs.PlateCarree(),
        ha="center",
        va="center",
        fontsize=FONT_GEO_LABEL,
        fontweight="bold",
        color="black",
        zorder=30,
    )

    axis.text(
        0.50,
        0.50,
        "90°S",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=FONT_GEO_LABEL + 1.0,
        fontweight="bold",
        color="black",
        zorder=30,
    )


def style_map_axis(axis):
    """
    Apply the common South Polar map style.
    """
    axis.set_extent(
        [
            -180,
            180,
            LAT_MIN,
            LAT_MAX,
        ],
        crs=ccrs.PlateCarree(),
    )

    add_circular_boundary(axis)

    axis.add_feature(
        cfeature.LAND.with_scale("110m"),
        facecolor=LAND_COLOR,
        edgecolor=COAST_COLOR,
        linewidth=0.65,
        zorder=6,
    )

    axis.coastlines(
        resolution="110m",
        color=COAST_COLOR,
        linewidth=0.55,
        zorder=7,
    )

    # Geographic gridlines only; no sea-sector division lines.
    axis.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        xlocs=MERIDIANS,
        ylocs=PARALLELS,
        linewidth=0.48,
        linestyle=":",
        color=GRID_COLOR,
        alpha=0.88,
        zorder=3,
    )

    add_geographic_labels(axis)

    for spine in axis.spines.values():

        spine.set_linewidth(1.0)
        spine.set_color("black")


def plot_map_field(
    axis,
    data_array,
    common_no_data_mask,
    color_limit,
    panel_title,
):
    """
    Plot one seasonal field with the fixed/common mask.
    """
    style_map_axis(axis)

    longitude_2d, latitude_2d = coordinate_mesh(
        data_array
    )

    axis.pcolormesh(
        longitude_2d,
        latitude_2d,
        common_no_data_mask.values,
        transform=ccrs.PlateCarree(),
        cmap=ListedColormap([NO_DATA_COLOR]),
        vmin=0.0,
        vmax=1.0,
        shading="auto",
        zorder=1,
    )

    image = axis.pcolormesh(
        longitude_2d,
        latitude_2d,
        data_array.values,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        vmin=-color_limit,
        vmax=color_limit,
        shading="auto",
        rasterized=True,
        zorder=4,
    )

    if DRAW_1000_M_CONTOUR:

        axis.contour(
            gebco_lon2d,
            gebco_lat2d,
            gebco.values,
            levels=[
                BATHYMETRY_CONTOUR_LEVEL
            ],
            colors="black",
            linewidths=0.68,
            linestyles="-",
            alpha=0.92,
            transform=ccrs.PlateCarree(),
            zorder=9,
        )

    axis.set_title(
        panel_title,
        fontsize=FONT_PANEL_TITLE,
        fontweight="bold",
        pad=22,
    )

    return image


# 23. BUILD THE 4 × 2 FIGURE
projection = ccrs.SouthPolarStereo(
    central_longitude=0.0
)

fig = plt.figure(
    figsize=FIGSIZE
)

# The central column holds one vertical multi-line statistics box
# per seasonal row.
grid = gridspec.GridSpec(
    nrows=4,
    ncols=3,
    figure=fig,
    width_ratios=[
        1.0,
        0.145,
        1.0,
    ],
    height_ratios=[
        1,
        1,
        1,
        1,
    ],
    left=0.110,
    right=0.985,
    top=0.975,
    bottom=0.105,
    wspace=0.020,
    hspace=0.165,
)

map_axes = np.empty(
    (4, 2),
    dtype=object,
)

panel_title_left = {
    "Spring": "(a) Spring — Observed SLA",
    "Summer": "(c) Summer — Observed SLA",
    "Autumn": "(e) Autumn — Observed SLA",
    "Winter": "(g) Winter — Observed SLA",
}

panel_title_right = {
    "Spring": "(b) Spring — EN4 total steric",
    "Summer": "(d) Summer — EN4 total steric",
    "Autumn": "(f) Autumn — EN4 total steric",
    "Winter": "(h) Winter — EN4 total steric",
}

last_sla_image = None
last_en4_image = None

for row_index, season in enumerate(SEASON_ORDER):

    left_axis = fig.add_subplot(
        grid[row_index, 0],
        projection=projection,
    )

    middle_axis = fig.add_subplot(
        grid[row_index, 1]
    )

    right_axis = fig.add_subplot(
        grid[row_index, 2],
        projection=projection,
    )

    middle_axis.set_axis_off()

    map_axes[row_index, 0] = left_axis
    map_axes[row_index, 1] = right_axis

    last_sla_image = plot_map_field(
        axis=left_axis,
        data_array=sla_anomaly_common[season],
        common_no_data_mask=no_data_mask[season],
        color_limit=sla_color_limit,
        panel_title=panel_title_left[season],
    )

    last_en4_image = plot_map_field(
        axis=right_axis,
        data_array=en4_anomaly_common[season],
        common_no_data_mask=no_data_mask[season],
        color_limit=en4_color_limit,
        panel_title=panel_title_right[season],
    )

    statistics_row = (
        statistics_table.loc[
            statistics_table["season"]
            == season
        ]
        .iloc[0]
    )

    statistics_text = (
        f"r$_{{spatial}}$ = "
        f"{statistics_row['spatial_correlation_r']:.2f}\n"
        f"RMSD = "
        f"{statistics_row['RMSD_cm']:.2f} cm\n"
        f"Bias = "
        f"{statistics_row['mean_bias_EN4_minus_SLA_cm']:+.2f} cm\n"
        f"n = "
        f"{int(statistics_row['n_common_cells']):,}\n"
        f"Coverage = "
        f"{statistics_row['spatial_coverage_percent']:.1f}%"
    )

    middle_axis.text(
        0.50,
        0.50,
        statistics_text,
        rotation=90,
        rotation_mode="anchor",
        transform=middle_axis.transAxes,
        ha="center",
        va="center",
        fontsize=FONT_STATISTICS,
        fontweight="bold",
        color="black",
        linespacing=1.18,
        bbox={
            "boxstyle": "round,pad=0.55",
            "facecolor": "white",
            "edgecolor": "0.22",
            "linewidth": 1.05,
            "alpha": 1.0,
        },
    )


# 24. LEFT-SIDE SEASON LABELS
for row_index, season in enumerate(SEASON_ORDER):

    position = map_axes[
        row_index,
        0
    ].get_position()

    fig.text(
        position.x0 - 0.050,
        position.y0
        + position.height / 2,
        SEASON_ROW_LABEL[season],
        rotation=90,
        rotation_mode="anchor",
        ha="center",
        va="center",
        fontsize=FONT_ROW_LABEL,
        fontweight="bold",
        color="black",
    )


# 25. BOTTOM COLORBARS
left_bottom_position = map_axes[
    3,
    0
].get_position()

right_bottom_position = map_axes[
    3,
    1
].get_position()

colorbar_y = 0.062
colorbar_height = 0.014

left_colorbar_axis = fig.add_axes(
    [
        left_bottom_position.x0 + 0.012,
        colorbar_y,
        left_bottom_position.width - 0.024,
        colorbar_height,
    ]
)

right_colorbar_axis = fig.add_axes(
    [
        right_bottom_position.x0 + 0.012,
        colorbar_y,
        right_bottom_position.width - 0.024,
        colorbar_height,
    ]
)

sla_colorbar = fig.colorbar(
    last_sla_image,
    cax=left_colorbar_axis,
    orientation="horizontal",
    extend="both",
)

sla_colorbar.set_label(
    "Observed SLA anomaly (cm)",
    fontsize=FONT_COLORBAR,
    fontweight="bold",
    labelpad=4,
)

sla_colorbar.locator = MaxNLocator(
    nbins=6
)
sla_colorbar.update_ticks()

sla_colorbar.ax.tick_params(
    labelsize=FONT_COLORBAR_TICK,
    width=1.0,
    length=4,
)

for label in sla_colorbar.ax.get_xticklabels():
    label.set_fontweight("bold")

en4_colorbar = fig.colorbar(
    last_en4_image,
    cax=right_colorbar_axis,
    orientation="horizontal",
    extend="both",
)

en4_colorbar.set_label(
    "EN4 total steric-height anomaly, 0–1000 m (cm)",
    fontsize=FONT_COLORBAR,
    fontweight="bold",
    labelpad=4,
)

en4_colorbar.locator = MaxNLocator(
    nbins=6
)
en4_colorbar.update_ticks()

en4_colorbar.ax.tick_params(
    labelsize=FONT_COLORBAR_TICK,
    width=1.0,
    length=4,
)

for label in en4_colorbar.ax.get_xticklabels():
    label.set_fontweight("bold")


# 26. BOTTOM LEGEND
legend_handles = [
    Line2D(
        [0],
        [0],
        color="black",
        linewidth=1.6,
        label="1000 m bathymetry",
    ),

    mpatches.Patch(
        facecolor=NO_DATA_COLOR,
        edgecolor="0.45",
        linewidth=0.8,
        label="Outside fixed/common valid mask",
    ),
]

bottom_legend = fig.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(
        0.50,
        0.018,
    ),
    ncol=2,
    frameon=False,
    fontsize=FONT_LEGEND,
    handlelength=2.0,
    handletextpad=0.7,
    columnspacing=1.8,
)

for text in bottom_legend.get_texts():
    text.set_fontweight("bold")


# 27. SAVE THE FIGURE
fig.savefig(
    OUT_PNG,
    dpi=SAVE_DPI,
    bbox_inches="tight",
    pad_inches=0.10,
    facecolor="white",
    edgecolor="none",
)

fig.savefig(
    OUT_PDF,
    bbox_inches="tight",
    pad_inches=0.10,
    facecolor="white",
    edgecolor="none",
)

plt.show()
plt.close(fig)


# 28. FINAL REPORT
print("\n" + "=" * 100)
print("FIGURE 7 — SLA VERSUS EN4 TOTAL STERIC COMPLETED")
print("=" * 100)

print("\n1080-dpi PNG:")
print(OUT_PNG)

print("\nVector PDF:")
print(OUT_PDF)

print("\nSeasonal statistics:")
print(OUT_STATS_CSV)

print("\nFixed common mask:")
print(OUT_MASK_NC)

print("\nSeasonal anomaly fields:")
print(OUT_SEASONAL_NC)

print("\nMethod summary:")
print(
    "Both columns use the same EN4-grid common mask."
)
print(
    "Bias is defined as EN4 total steric minus observed SLA."
)
print(
    "No Argo objective mapping or Argo grid-cell dots were used."
)
print(
    "Argo validation should be presented separately in the "
    "supplementary material."
)

ds_sla.close()
ds_en4.close()
ds_sic.close()
ds_gebco.close()
