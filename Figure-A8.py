"""
SUPPLEMENTARY FIGURE S4 — FULL ORGANIZED LATEST WORKFLOW

FIGURE
EN4 Integration-Depth Sensitivity
Southern Ocean, 2008-2025

SOURCE
Extracted from:
    so_sealevel_paper_fig.py

VERSION SELECTION
The supplied source contains one Figure S4 implementation. It is already a
complete processing-and-plotting workflow, so there is no older duplicated
S4 renderer to retain or merge.

This standalone script preserves that full Figure S4 workflow and only
reorganizes the startup order to:

    Google Drive mount
        -> package installation
        -> imports
        -> full S4 processing
        -> final outputs

FIGURE STRUCTURE
                    0-1000 m       0-2000 m       Difference

Total steric            (a)             (b)             (c)
Thermosteric            (d)             (e)             (f)
Halosteric              (g)             (h)             (i)

COMMON SPATIAL RULE

Every map uses the SAME fixed spatial mask:

    - GEBCO water depth >= 2000 m;
    - finite 0-1000 m total / thermosteric / halosteric values;
    - finite 0-2000 m total / thermosteric / halosteric values;
    - required valid-data fraction through the complete study period.

The 0-1000 m maps are therefore not allowed to display cells that are
unavailable in the 0-2000 m products.

DIFFERENCE DEFINITION

    Difference = 0-2000 m field - 0-1000 m field

DEFAULT DISPLAYED DIAGNOSTIC

The source default is:

    COMPARISON_MODE = "winter_minus_summer"

defined as:

    Winter seasonal climatology - Summer seasonal climatology

This signed seasonal contrast is retained because the Figure S4 layout contains
one map for each depth/component combination and uses a diverging color scale.

The source also supports:

    "spring_anomaly"
    "summer_anomaly"
    "autumn_anomaly"
    "winter_anomaly"
    "seasonal_peak_to_peak"
    "seasonal_rms"

0-2000 M PRODUCT CREATION
If all three 0-2000 m products already exist, the workflow can reuse them.

Otherwise, the script creates them month-by-month from the raw EN4 temperature
and salinity dataset using the same fixed 2008-2025 grid-cell reference-state
method used by the corrected 0-1000 m products.

The three 0-2000 m products are:

    EN4_total_steric_0_2000m_monthly_2008_2025_SO.nc
    EN4_thermosteric_0_2000m_monthly_2008_2025_SO.nc
    EN4_halosteric_0_2000m_monthly_2008_2025_SO.nc

MAIN INPUTS

Raw EN4:
    /content/drive/MyDrive/SAM_Thesis/Data/EN4/final_monthly/
    EN4.2.2_g10_Antarctic_monthly_2008_2025.nc

Existing corrected 0-1000 m products:
    EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc
    EN4_thermosteric_0_1000m_monthly_2008_2025_SO.nc
    EN4_halosteric_0_1000m_monthly_2008_2025_SO.nc

FINAL FIGURE / ANALYSIS OUTPUTS

/content/drive/MyDrive/SAM_Thesis/paper2/

    Figure_S4_EN4_depth_sensitivity_1080dpi.png
    Figure_S4_EN4_depth_sensitivity.pdf
    Figure_S4_EN4_depth_sensitivity_fields.nc
    Figure_S4_EN4_depth_sensitivity_summary.csv
    Figure_S4_EN4_depth_sensitivity_summary.txt

ORGANIZED EXECUTION ORDER
1. Mount Google Drive.
2. Install all required packages.
3. Import all required libraries.
4. Define raw EN4, derived-product and output paths.
5. Define analysis years, latitude domain and integration depths.
6. Configure reuse/local staging for 0-2000 m product creation.
7. Define the strict common valid-data fraction.
8. Define comparison mode and Southern Hemisphere seasons.
9. Define plotting/layout settings.
10. Define robust NetCDF, coordinate, masking and interpolation helpers.
11. Verify all required 0-1000 m and raw EN4 inputs.
12. Create or reuse the three monthly 0-2000 m products.
13. Open all six monthly steric products.
14. Align them in time and on the same grid.
15. Create the fixed common mask requiring GEBCO depth >=2000 m.
16. Require the configured valid-data fraction in all six products.
17. Calculate the requested depth-sensitivity diagnostic.
18. Calculate 0-2000 m minus 0-1000 m difference fields.
19. Save the processed 3-component S4 fields to NetCDF.
20. Calculate and save the S4 summary table to CSV.
21. Calculate common color limits.
22. Draw the 3 x 3 polar-map figure.
23. Add final column headings and component labels.
24. Add bottom colorbars.
25. Save 1080-dpi PNG and vector PDF.
26. Save the text summary.
27. Close all six opened product datasets.

SOURCE SETTINGS PRESERVED
Study period:
    2008-2025

Latitude domain:
    -90 to -60 degrees

Integration depths:
    1000 m and 2000 m

Reuse existing 0-2000 m products:
    True

Stage the large raw EN4 file locally:
    True

Local cache:
    /content/Figure_S4_EN4_cache

Fixed common valid-data fraction:
    1.00

This means all 216 study months must be valid in all six comparison products
for a grid cell to enter the common comparison mask.

SCIENTIFIC RULES PRESERVED
The source final summary states that:

    - both depth columns use exactly the same common >2000 m cells;
    - differences are calculated after applying the common mask;
    - 0-2000 m products use one fixed 2008-2025 grid-cell reference profile;
    - monthly T/S must cover the complete 0-2000 m layer;
    - shallow cells are excluded rather than being integrated only to local
      seabed depth;
    - halosteric height is total minus corrected thermosteric height.

ORGANIZATION CHANGES
Only structural organization has been changed:

    - Google Drive mounting is moved to the beginning.
    - Package installation follows the mount.
    - Library imports follow installation.
    - The complete original S4 workflow follows unchanged.
    - Table S4 and every later source block are excluded.

No S4 comparison definition, integration method, common-mask rule, valid-data
threshold, reference-state method, product filenames, summary products, or
plotting methodology are intentionally changed.

"""

# FIGURE S4 — EN4 INTEGRATION-DEPTH SENSITIVITY
# Southern Ocean, 2008–2025
#
# Structure
#                     0–1000 m       0–2000 m       Difference
# Total steric            (a)             (b)             (c)
# Thermosteric            (d)             (e)             (f)
# Halosteric              (g)             (h)             (i)
#
# Spatial rule
# Every map uses the SAME fixed spatial mask:
#   • GEBCO water depth >= 2000 m
#   • finite 0–1000 m total / thermo / halo values
#   • finite 0–2000 m total / thermo / halo values
#   • required valid-data fraction through the study period
#
# The 0–1000 m maps are therefore NOT allowed to show cells that are
# unavailable in the 0–2000 m products.
#
# Difference
# Difference = 0–2000 m field − 0–1000 m field
#
# Default displayed diagnostic
# WINTER_MINUS_SUMMER:
#   Winter seasonal climatology − Summer seasonal climatology
#
# This signed seasonal contrast is used because Figure S4 contains only
# one map for each depth/component combination. It retains direction and
# is compatible with a diverging red–blue colour scale.
#
# Other supported diagnostics can be selected below:
#   "winter_minus_summer"
#   "spring_anomaly"
#   "summer_anomaly"
#   "autumn_anomaly"
#   "winter_anomaly"
#   "seasonal_peak_to_peak"
#   "seasonal_rms"
#
# 0–2000 m files
# If the three 0–2000 m monthly products already exist, the script reuses
# them. Otherwise, it creates them month-by-month from the raw EN4 T–S
# file using the same fixed 2008–2025 grid-cell reference-state method
# used for the corrected 0–1000 m products.
#
# No overall figure title is added.


# 1. MOUNT GOOGLE DRIVE
try:
    from google.colab import drive
    drive.mount("/content/drive", force_remount=False)
except Exception:
    print("Google Drive mounting skipped because this is not Colab.")
# 2. INSTALL REQUIRED PACKAGES
import sys
import subprocess
import importlib.util

REQUIRED_PACKAGES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "xarray": "xarray",
    "h5netcdf": "h5netcdf",
    "netCDF4": "netCDF4",
    "cftime": "cftime",
    "dask": "dask[array]",
    "scipy": "scipy",
    "gsw": "gsw",
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
        [sys.executable, "-m", "pip", "install", "-q", *missing_packages]
    )
    print("Package installation completed.")
else:
    print("All required packages are already installed.")
# 3. IMPORT LIBRARIES

import os
import gc
import shutil
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import gsw

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from scipy.interpolate import interp1d

import cartopy.crs as ccrs
import cartopy.feature as cfeature

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
# 4. COMPLETE SUPPLEMENTARY FIGURE S4 WORKFLOW

# 3. INPUT AND OUTPUT PATHS
BASE_DIR = Path("/content/drive/MyDrive/SAM_Thesis")

RAW_EN4_FILE = (
    BASE_DIR
    / "Data/EN4/final_monthly/"
    / "EN4.2.2_g10_Antarctic_monthly_2008_2025.nc"
)

DERIVED_DIR = (
    BASE_DIR
    / "Processed/EN4_NetCDF_inventory"
)

OUTPUT_DIR = (
    BASE_DIR
    / "paper2"
)

DERIVED_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Existing corrected 0–1000 m products
FILES_1000 = {
    "total": (
        DERIVED_DIR
        / "EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc"
    ),
    "thermo": (
        DERIVED_DIR
        / "EN4_thermosteric_0_1000m_monthly_2008_2025_SO.nc"
    ),
    "halo": (
        DERIVED_DIR
        / "EN4_halosteric_0_1000m_monthly_2008_2025_SO.nc"
    ),
}

# 0–2000 m products created/reused by this script
FILES_2000 = {
    "total": (
        DERIVED_DIR
        / "EN4_total_steric_0_2000m_monthly_2008_2025_SO.nc"
    ),
    "thermo": (
        DERIVED_DIR
        / "EN4_thermosteric_0_2000m_monthly_2008_2025_SO.nc"
    ),
    "halo": (
        DERIVED_DIR
        / "EN4_halosteric_0_2000m_monthly_2008_2025_SO.nc"
    ),
}

OUT_PNG = (
    OUTPUT_DIR
    / "Figure_S4_EN4_depth_sensitivity_1080dpi.png"
)

OUT_PDF = (
    OUTPUT_DIR
    / "Figure_S4_EN4_depth_sensitivity.pdf"
)

OUT_FIELDS_NC = (
    OUTPUT_DIR
    / "Figure_S4_EN4_depth_sensitivity_fields.nc"
)

OUT_SUMMARY_CSV = (
    OUTPUT_DIR
    / "Figure_S4_EN4_depth_sensitivity_summary.csv"
)

OUT_SUMMARY_TXT = (
    OUTPUT_DIR
    / "Figure_S4_EN4_depth_sensitivity_summary.txt"
)


# 4. ANALYSIS SETTINGS
YEAR_START = 2008
YEAR_END = 2025

LAT_MIN = -90.0
LAT_MAX = -60.0

DEPTH_1000 = 1000.0
DEPTH_2000 = 2000.0

# Reuse existing 0–2000 products when all three exist.
REUSE_EXISTING_2000_PRODUCTS = True

# Copy the approximately 2 GB raw EN4 file from Drive to local Colab
# storage before computing 0–2000 m products. This avoids many HDF5
# read failures caused by repeatedly reading a large file from Drive.
STAGE_RAW_EN4_TO_LOCAL = True
LOCAL_CACHE_DIR = Path("/content/Figure_S4_EN4_cache")

# Monthly support required for the fixed comparison mask.
# 1.00 = all 216 months must be valid in all six products.
COMMON_VALID_FRACTION = 1.00

# Default signed diagnostic displayed in the nine maps.
COMPARISON_MODE = "winter_minus_summer"

SUPPORTED_COMPARISON_MODES = {
    "winter_minus_summer",
    "spring_anomaly",
    "summer_anomaly",
    "autumn_anomaly",
    "winter_anomaly",
    "seasonal_peak_to_peak",
    "seasonal_rms",
}

if COMPARISON_MODE not in SUPPORTED_COMPARISON_MODES:
    raise ValueError(
        f"Unsupported COMPARISON_MODE={COMPARISON_MODE!r}. "
        f"Choose one of {sorted(SUPPORTED_COMPARISON_MODES)}"
    )

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

COMPONENT_ORDER = [
    "total",
    "thermo",
    "halo",
]

COMPONENT_LABELS = {
    "total": "Total steric",
    "thermo": "Thermosteric",
    "halo": "Halosteric",
}

COLUMN_LABELS = [
    "0–1000 m",
    "0–2000 m",
    "Difference",
]

PANEL_LETTERS = [
    ["(a)", "(b)", "(c)"],
    ["(d)", "(e)", "(f)"],
    ["(g)", "(h)", "(i)"],
]

DIFFERENCE_DEFINITION = "0–2000 m minus 0–1000 m"

MERIDIANS = np.arange(-180, 181, 30)
PARALLELS = [-60, -70, -80]

DRAW_2000_M_CONTOUR = True
BATHYMETRY_CONTOUR_LEVEL = 2000.0

SAVE_DPI = 1080


# 5. FIGURE SETTINGS
FIGSIZE = (11.6, 11.2)

FONT_COLUMN_TITLE = 18.0
FONT_PANEL_LETTER = 20.0
FONT_ROW_LABEL = 18.0
FONT_GEO_LABEL = 9.0
FONT_COLORBAR = 12.0
FONT_COLORBAR_TICK = 9.0

NO_DATA_COLOR = "0.84"
LAND_COLOR = "0.78"
COAST_COLOR = "0.30"
GRID_COLOR = "0.56"

RIM_LABEL_LATITUDE = -58.2

# Panel letters are placed outside the 30°W rim label.
PANEL_LETTER_X = -0.155
PANEL_LETTER_Y = 1.055

COLUMN_TITLE_OFFSET = 0.050

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "axes.linewidth": 1.0,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "mathtext.default": "regular",
    }
)


# 6. GENERIC HELPERS
def print_header(text):
    print("\n" + "=" * 92)
    print(text)
    print("=" * 92)


def safe_is_datetime(dtype):
    try:
        return np.issubdtype(dtype, np.datetime64)
    except TypeError:
        return False


def open_dataset_safely(
    file_path,
    chunks=None,
    decode_times=True,
):
    """
    Open a NetCDF file using several engines.
    h5netcdf is attempted first because it is generally reliable for
    large HDF5-backed EN4 files.
    """
    attempts = []

    for engine in [
        "h5netcdf",
        "netcdf4",
        None,
        "scipy",
    ]:
        for decode in (
            [decode_times]
            if decode_times is False
            else [True, False]
        ):
            try:
                kwargs = {
                    "decode_times": decode,
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
                    f"Opened {Path(file_path).name} "
                    f"with engine={engine_name}, "
                    f"decode_times={decode}"
                )

                return dataset

            except Exception as error:
                attempts.append(
                    f"engine={engine}, decode_times={decode}: {error}"
                )

    raise RuntimeError(
        f"Could not open:\n{file_path}\n\n"
        + "\n".join(attempts)
    )


def stage_file_to_local(source_path, cache_directory):
    """
    Copy a large Drive file to local Colab storage when required.
    """
    source_path = Path(source_path)
    cache_directory.mkdir(parents=True, exist_ok=True)

    destination = cache_directory / source_path.name

    if destination.exists():
        source_size = source_path.stat().st_size
        destination_size = destination.stat().st_size

        if source_size == destination_size:
            print("Using existing local staged file:")
            print(destination)
            return destination

        destination.unlink()

    print("Copying large EN4 input from Drive to local runtime...")
    print("Source     :", source_path)
    print("Destination:", destination)

    shutil.copy2(source_path, destination)

    print("Local staging completed.")
    return destination


def detect_coordinate(dataset, coordinate_type):
    names = list(dataset.coords) + [
        name
        for name in dataset.variables
        if name not in dataset.coords
    ]

    aliases = {
        "time": ["time", "valid_time", "date", "datetime", "month", "t"],
        "depth": ["depth", "lev", "level", "pressure", "pres", "z"],
        "lat": ["lat", "latitude", "nav_lat", "y"],
        "lon": ["lon", "longitude", "nav_lon", "x"],
    }

    for name in names:
        variable = dataset[name]

        lower_name = name.lower()
        standard_name = str(
            variable.attrs.get("standard_name", "")
        ).lower()
        axis = str(
            variable.attrs.get("axis", "")
        ).upper()
        units = str(
            variable.attrs.get("units", "")
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

        elif coordinate_type == "depth":
            if (
                lower_name in aliases["depth"]
                or standard_name in {
                    "depth",
                    "sea_water_pressure",
                }
                or axis == "Z"
                or units in {
                    "m",
                    "meter",
                    "metre",
                    "dbar",
                }
            ):
                return name

        elif coordinate_type == "lat":
            if (
                lower_name in aliases["lat"]
                or standard_name == "latitude"
                or axis == "Y"
                or "degree_north" in units
            ):
                return name

        elif coordinate_type == "lon":
            if (
                lower_name in aliases["lon"]
                or standard_name == "longitude"
                or axis == "X"
                or "degree_east" in units
            ):
                return name

    raise KeyError(
        f"Could not identify the {coordinate_type} coordinate. "
        f"Available coordinates/variables: {names}"
    )


def standardize_coordinates(
    dataset,
    time_name=None,
    depth_name=None,
    lat_name=None,
    lon_name=None,
):
    rename_mapping = {}

    if time_name is not None and time_name != "time":
        rename_mapping[time_name] = "time"

    if depth_name is not None and depth_name != "depth":
        rename_mapping[depth_name] = "depth"

    if lat_name is not None and lat_name != "lat":
        rename_mapping[lat_name] = "lat"

    if lon_name is not None and lon_name != "lon":
        rename_mapping[lon_name] = "lon"

    if rename_mapping:
        dataset = dataset.rename(rename_mapping)

    return dataset


def normalize_and_sort_longitude(dataset):
    if "lon" not in dataset.coords:
        return dataset

    if dataset["lon"].ndim != 1:
        raise ValueError("Only one-dimensional longitude is supported.")

    longitude = (
        (dataset["lon"].astype(float) + 180.0) % 360.0
    ) - 180.0

    dataset = dataset.assign_coords(lon=longitude)

    _, unique_indices = np.unique(
        dataset["lon"].values,
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
    latitude = dataset["lat"].values

    if latitude[0] <= latitude[-1]:
        output = dataset.sel(
            lat=slice(
                southern_limit,
                northern_limit,
            )
        )
    else:
        output = dataset.sel(
            lat=slice(
                northern_limit,
                southern_limit,
            )
        )

    return output.sortby("lat")


def choose_data_variable(
    dataset,
    preferred_names,
    excluded_names=None,
):
    if excluded_names is None:
        excluded_names = []

    variables = [
        name
        for name in dataset.data_vars
        if name not in excluded_names
    ]

    for preferred_name in preferred_names:
        for variable_name in variables:
            if preferred_name.lower() in variable_name.lower():
                return variable_name

    if not variables:
        raise ValueError(
            "No suitable science variable was found."
        )

    return variables[0]


def find_exact_or_partial_variable(dataset, names):
    for name in names:
        if name in dataset.data_vars:
            return name

    for name in names:
        for candidate in dataset.data_vars:
            if name.lower() in candidate.lower():
                return candidate

    raise KeyError(
        f"Could not find any of {names}. "
        f"Available variables: {list(dataset.data_vars)}"
    )


def target_depth_grid(native_depth, maximum_depth):
    """
    Construct a depth grid containing:
      • every native EN4 level inside the target layer
      • an explicit 0 m upper boundary
      • an explicit maximum-depth lower boundary
    """
    native_depth = np.asarray(
        native_depth,
        dtype=float,
    )

    native_inside = native_depth[
        (native_depth > 0.0)
        & (native_depth < maximum_depth)
    ]

    mandatory = np.array(
        [0.0, maximum_depth],
        dtype=float,
    )

    target = np.unique(
        np.round(
            np.concatenate(
                [native_inside, mandatory]
            ),
            6,
        )
    )

    target.sort()
    return target


def interpolate_depth(
    values,
    native_depth,
    target_depth,
):
    """
    Linear interpolation along the first (depth) axis.
    The approximately 5 m shallowest EN4 value represents the 0 m
    boundary, matching the existing corrected 0–1000 m workflow.
    """
    interpolator = interp1d(
        np.asarray(native_depth, dtype=float),
        np.asarray(values, dtype=float),
        axis=0,
        kind="linear",
        bounds_error=False,
        fill_value=np.nan,
        assume_sorted=True,
    )

    output = interpolator(
        np.asarray(target_depth, dtype=float)
    )

    zero_index = int(
        np.argmin(
            np.abs(
                np.asarray(target_depth, dtype=float)
                - 0.0
            )
        )
    )

    output[zero_index] = np.asarray(
        values,
        dtype=float,
    )[0]

    return output


def trapezoid(values, coordinate, axis=0):
    function = (
        np.trapezoid
        if hasattr(np, "trapezoid")
        else np.trapz
    )

    return function(
        values,
        x=coordinate,
        axis=axis,
    )


def integrate_steric(
    delta_specific_volume,
    pressure_dbar,
    gravity,
):
    """
    1 dbar = 1e4 Pa. The result is metres.
    """
    return trapezoid(
        delta_specific_volume
        * 1.0e4
        / gravity,
        pressure_dbar,
        axis=0,
    )


def convert_height_to_cm(data_array):
    units = str(
        data_array.attrs.get("units", "")
    ).strip().lower()

    if units in {
        "cm",
        "centimeter",
        "centimeters",
        "centimetre",
        "centimetres",
    }:
        output = data_array.copy()
        output.attrs["units"] = "cm"
        return output

    if units in {
        "m",
        "meter",
        "meters",
        "metre",
        "metres",
    }:
        output = data_array * 100.0
        output.attrs.update(data_array.attrs)
        output.attrs["units"] = "cm"
        return output

    sample = np.asarray(
        data_array.isel(
            time=slice(
                0,
                min(
                    12,
                    data_array.sizes.get("time", 1),
                ),
            )
        ).values,
        dtype=float,
    )

    sample = sample[np.isfinite(sample)]

    if (
        sample.size > 0
        and np.nanpercentile(
            np.abs(sample),
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


def robust_symmetric_limit(
    fields,
    percentile=99.0,
):
    values = []

    for field in fields:
        array = np.asarray(
            field.values,
            dtype=float,
        )

        array = array[np.isfinite(array)]

        if array.size:
            values.append(array)

    if not values:
        return 1.0

    combined = np.concatenate(values)

    maximum = float(
        np.nanpercentile(
            np.abs(combined),
            percentile,
        )
    )

    if maximum <= 1.0:
        interval = 0.25
    elif maximum <= 2.0:
        interval = 0.5
    elif maximum <= 5.0:
        interval = 1.0
    elif maximum <= 10.0:
        interval = 2.0
    elif maximum <= 20.0:
        interval = 5.0
    else:
        interval = 10.0

    return max(
        interval,
        float(
            interval
            * np.ceil(
                maximum / interval
            )
        ),
    )


def robust_positive_limit(
    fields,
    percentile=99.0,
):
    values = []

    for field in fields:
        array = np.asarray(
            field.values,
            dtype=float,
        )

        array = array[np.isfinite(array)]

        if array.size:
            values.append(array)

    if not values:
        return 1.0

    combined = np.concatenate(values)

    maximum = float(
        np.nanpercentile(
            combined,
            percentile,
        )
    )

    if maximum <= 1.0:
        interval = 0.25
    elif maximum <= 2.0:
        interval = 0.5
    elif maximum <= 5.0:
        interval = 1.0
    elif maximum <= 10.0:
        interval = 2.0
    elif maximum <= 20.0:
        interval = 5.0
    else:
        interval = 10.0

    return max(
        interval,
        float(
            interval
            * np.ceil(
                maximum / interval
            )
        ),
    )


def coordinate_mesh(data_array):
    return np.meshgrid(
        data_array["lon"].values,
        data_array["lat"].values,
    )


def format_longitude(longitude):
    value = int(longitude)

    if value == 0:
        return "0°"

    if abs(value) == 180:
        return "180°"

    if value < 0:
        return f"{abs(value)}°W"

    return f"{value}°E"


# 7. VERIFY REQUIRED INPUTS
for component, file_path in FILES_1000.items():
    if not file_path.exists():
        raise FileNotFoundError(
            f"Missing corrected 0–1000 m {component} product:\n"
            f"{file_path}"
        )

have_all_2000_files = all(
    file_path.exists()
    for file_path in FILES_2000.values()
)

if not have_all_2000_files:
    if not RAW_EN4_FILE.exists():
        raise FileNotFoundError(
            "The 0–2000 m products do not exist and the raw EN4 file "
            f"was not found:\n{RAW_EN4_FILE}"
        )


# 8. CREATE 0–2000 M MONTHLY PRODUCTS WHEN NEEDED
def create_2000m_products():
    print_header(
        "CREATING CORRECTED EN4 0–2000 M MONTHLY STERIC PRODUCTS"
    )

    raw_input = RAW_EN4_FILE

    if STAGE_RAW_EN4_TO_LOCAL:
        raw_input = stage_file_to_local(
            RAW_EN4_FILE,
            LOCAL_CACHE_DIR,
        )

    raw = open_dataset_safely(
        raw_input,
        chunks={"time": 1},
    )

    time_name = detect_coordinate(raw, "time")
    depth_name = detect_coordinate(raw, "depth")
    lat_name = detect_coordinate(raw, "lat")
    lon_name = detect_coordinate(raw, "lon")

    raw = standardize_coordinates(
        raw,
        time_name=time_name,
        depth_name=depth_name,
        lat_name=lat_name,
        lon_name=lon_name,
    )

    raw = subset_antarctic_latitudes(
        normalize_and_sort_longitude(raw)
    )

    temperature_name = find_exact_or_partial_variable(
        raw,
        [
            "temperature",
            "potential_temperature",
            "temp",
        ],
    )

    salinity_name = find_exact_or_partial_variable(
        raw,
        [
            "salinity",
            "practical_salinity",
            "psal",
        ],
    )

    temperature = raw[temperature_name].sel(
        time=slice(
            f"{YEAR_START}-01-01",
            f"{YEAR_END}-12-31",
        )
    )

    salinity = raw[salinity_name].sel(
        time=temperature["time"]
    )

    dates = pd.to_datetime(
        temperature["time"].values
    )

    native_depth = np.asarray(
        raw["depth"].values,
        dtype=float,
    )

    latitude = np.asarray(
        raw["lat"].values,
        dtype=float,
    )

    longitude = np.asarray(
        raw["lon"].values,
        dtype=float,
    )

    number_of_times = len(dates)
    number_of_latitudes = len(latitude)
    number_of_longitudes = len(longitude)

    print(
        f"Raw EN4 dimensions: time={number_of_times}, "
        f"depth={len(native_depth)}, "
        f"lat={number_of_latitudes}, "
        f"lon={number_of_longitudes}"
    )

    analysis_depth = target_depth_grid(
        native_depth,
        DEPTH_2000,
    )

    print(
        "0–2000 m analysis levels:",
        len(analysis_depth),
    )

    # The existing corrected 0–1000 m file contains GEBCO water depth
    # already interpolated onto the EN4 grid.
    mask_source = open_dataset_safely(
        FILES_1000["total"],
        chunks=None,
    )

    mask_source = standardize_coordinates(
        mask_source,
        time_name=detect_coordinate(mask_source, "time"),
        lat_name=detect_coordinate(mask_source, "lat"),
        lon_name=detect_coordinate(mask_source, "lon"),
    )

    mask_source = subset_antarctic_latitudes(
        normalize_and_sort_longitude(mask_source)
    )

    if "gebco_water_depth" not in mask_source.data_vars:
        raise KeyError(
            "The corrected 0–1000 m file does not contain "
            "'gebco_water_depth'."
        )

    water_depth = (
        mask_source["gebco_water_depth"]
        .sel(
            lat=latitude,
            lon=longitude,
        )
        .load()
        .values
        .astype(float)
    )

    bathymetry_mask = (
        np.isfinite(water_depth)
        & (water_depth >= DEPTH_2000)
    )

    print(
        "Cells with GEBCO depth >=2000 m:",
        f"{int(bathymetry_mask.sum()):,}",
    )

    # Fixed grid-cell all-month reference profiles.
    print(
        "Calculating fixed 2008–2025 all-month reference profiles..."
    )

    reference_temperature_native = (
        temperature.mean(
            dim="time",
            skipna=True,
        )
        .load()
        .values
        .astype(float)
    )

    reference_salinity_native = (
        salinity.mean(
            dim="time",
            skipna=True,
        )
        .load()
        .values
        .astype(float)
    )

    reference_temperature_k = interpolate_depth(
        reference_temperature_native,
        native_depth,
        analysis_depth,
    )

    reference_salinity = interpolate_depth(
        reference_salinity_native,
        native_depth,
        analysis_depth,
    )

    # EN4 temperature is treated as potential temperature in kelvin.
    reference_potential_temperature_c = (
        reference_temperature_k
        - 273.15
    )

    reference_complete = (
        np.isfinite(
            reference_potential_temperature_c
        ).all(axis=0)
        & np.isfinite(
            reference_salinity
        ).all(axis=0)
    )

    pressure_depth_lat = gsw.p_from_z(
        -analysis_depth[:, None],
        latitude[None, :],
    )

    pressure_3d = np.broadcast_to(
        pressure_depth_lat[:, :, None],
        (
            len(analysis_depth),
            number_of_latitudes,
            number_of_longitudes,
        ),
    )

    latitude_3d = latitude[None, :, None]
    longitude_3d = longitude[None, None, :]

    gravity_3d = gsw.grav(
        latitude_3d,
        pressure_3d,
    )

    reference_absolute_salinity = gsw.SA_from_SP(
        reference_salinity,
        pressure_3d,
        longitude_3d,
        latitude_3d,
    )

    reference_conservative_temperature = gsw.CT_from_pt(
        reference_absolute_salinity,
        reference_potential_temperature_c,
    )

    reference_specific_volume = gsw.specvol(
        reference_absolute_salinity,
        reference_conservative_temperature,
        pressure_3d,
    )

    total_output = np.full(
        (
            number_of_times,
            number_of_latitudes,
            number_of_longitudes,
        ),
        np.nan,
        dtype=np.float32,
    )

    thermo_output = np.full_like(
        total_output,
        np.nan,
    )

    halo_output = np.full_like(
        total_output,
        np.nan,
    )

    valid_output = np.zeros(
        total_output.shape,
        dtype=np.int8,
    )

    for time_index in range(number_of_times):
        date = pd.Timestamp(dates[time_index])

        print(
            f"[{time_index + 1:03d}/{number_of_times:03d}] "
            f"{date:%Y-%m}"
        )

        monthly_temperature_native = (
            temperature
            .isel(time=time_index)
            .load()
            .values
            .astype(float)
        )

        monthly_salinity_native = (
            salinity
            .isel(time=time_index)
            .load()
            .values
            .astype(float)
        )

        monthly_temperature_k = interpolate_depth(
            monthly_temperature_native,
            native_depth,
            analysis_depth,
        )

        monthly_salinity = interpolate_depth(
            monthly_salinity_native,
            native_depth,
            analysis_depth,
        )

        monthly_potential_temperature_c = (
            monthly_temperature_k
            - 273.15
        )

        monthly_complete = (
            np.isfinite(
                monthly_potential_temperature_c
            ).all(axis=0)
            & np.isfinite(
                monthly_salinity
            ).all(axis=0)
        )

        valid = (
            bathymetry_mask
            & reference_complete
            & monthly_complete
        )

        monthly_absolute_salinity = gsw.SA_from_SP(
            monthly_salinity,
            pressure_3d,
            longitude_3d,
            latitude_3d,
        )

        monthly_conservative_temperature = gsw.CT_from_pt(
            monthly_absolute_salinity,
            monthly_potential_temperature_c,
        )

        monthly_specific_volume = gsw.specvol(
            monthly_absolute_salinity,
            monthly_conservative_temperature,
            pressure_3d,
        )

        total_height = integrate_steric(
            monthly_specific_volume
            - reference_specific_volume,
            pressure_3d,
            gravity_3d,
        )

        thermosteric_conservative_temperature = gsw.CT_from_pt(
            reference_absolute_salinity,
            monthly_potential_temperature_c,
        )

        thermosteric_specific_volume = gsw.specvol(
            reference_absolute_salinity,
            thermosteric_conservative_temperature,
            pressure_3d,
        )

        thermosteric_height = integrate_steric(
            thermosteric_specific_volume
            - reference_specific_volume,
            pressure_3d,
            gravity_3d,
        )

        halosteric_height = (
            total_height
            - thermosteric_height
        )

        total_output[time_index] = np.where(
            valid,
            total_height,
            np.nan,
        ).astype(np.float32)

        thermo_output[time_index] = np.where(
            valid,
            thermosteric_height,
            np.nan,
        ).astype(np.float32)

        halo_output[time_index] = np.where(
            valid,
            halosteric_height,
            np.nan,
        ).astype(np.float32)

        valid_output[time_index] = valid.astype(
            np.int8
        )

        del (
            monthly_temperature_native,
            monthly_salinity_native,
            monthly_temperature_k,
            monthly_salinity,
            monthly_potential_temperature_c,
            monthly_absolute_salinity,
            monthly_conservative_temperature,
            monthly_specific_volume,
            thermosteric_conservative_temperature,
            thermosteric_specific_volume,
            total_height,
            thermosteric_height,
            halosteric_height,
        )

        gc.collect()

    coordinates = {
        "time": dates,
        "latitude": latitude.astype(np.float32),
        "longitude": longitude.astype(np.float32),
    }

    common_dataset_attributes = {
        "Conventions": "CF-1.8",
        "source_en4_file": str(RAW_EN4_FILE),
        "period": f"{YEAR_START}-01 through {YEAR_END}-12",
        "region": "Southern Ocean EN4 Antarctic domain",
        "integration_layer": "0–2000 m",
        "reference_state": (
            "One fixed grid-cell, depth-dependent all-month mean "
            "potential-temperature and practical-salinity profile "
            "for 2008–2025."
        ),
        "depth_mask_rule": (
            "GEBCO water depth >=2000 m and complete monthly/reference "
            "temperature and salinity through the complete 0–2000 m layer."
        ),
        "seasonal_cycle_preserved": (
            "Yes; no calendar-month climatology was removed."
        ),
        "temperature_interpretation": (
            "EN4 potential temperature in kelvin."
        ),
        "salinity_interpretation": (
            "EN4 Practical Salinity."
        ),
        "processing_software": (
            "Python, NumPy, Xarray, SciPy and GSW-Python"
        ),
        "processing_time_utc": (
            datetime.now(timezone.utc).isoformat()
        ),
    }

    output_specs = {
        "total": {
            "variable": "total_steric_height",
            "long_name": (
                "EN4 total steric height anomaly "
                "integrated over 0–2000 m"
            ),
            "values": total_output,
            "calculation": (
                "Integral of monthly minus fixed-reference "
                "specific volume over pressure divided by local gravity."
            ),
        },
        "thermo": {
            "variable": "thermosteric_height",
            "long_name": (
                "EN4 thermosteric height anomaly "
                "integrated over 0–2000 m"
            ),
            "values": thermo_output,
            "calculation": (
                "Monthly potential-temperature effect with Absolute "
                "Salinity held at the fixed reference profile."
            ),
        },
        "halo": {
            "variable": "halosteric_height",
            "long_name": (
                "EN4 halosteric height anomaly "
                "integrated over 0–2000 m"
            ),
            "values": halo_output,
            "calculation": (
                "Total steric height minus corrected thermosteric height."
            ),
        },
    }

    for component, specification in output_specs.items():
        output_dataset = xr.Dataset(
            {
                specification["variable"]: (
                    (
                        "time",
                        "latitude",
                        "longitude",
                    ),
                    specification["values"],
                ),
                "valid_layer_mask": (
                    (
                        "time",
                        "latitude",
                        "longitude",
                    ),
                    valid_output,
                ),
                "gebco_water_depth": (
                    (
                        "latitude",
                        "longitude",
                    ),
                    water_depth.astype(np.float32),
                ),
            },
            coords=coordinates,
        )

        output_dataset[
            specification["variable"]
        ].attrs.update(
            {
                "long_name": specification["long_name"],
                "units": "m",
                "calculation": specification["calculation"],
            }
        )

        output_dataset["valid_layer_mask"].attrs.update(
            {
                "long_name": "monthly valid 0–2000 m layer mask",
                "flag_values": np.array(
                    [0, 1],
                    dtype=np.int8,
                ),
                "flag_meanings": "invalid valid",
            }
        )

        output_dataset["gebco_water_depth"].attrs.update(
            {
                "long_name": "GEBCO ocean water depth",
                "units": "m",
                "positive": "down",
            }
        )

        output_dataset.attrs.update(
            common_dataset_attributes
        )

        encoding = {
            specification["variable"]: {
                "zlib": True,
                "complevel": 4,
                "shuffle": True,
                "dtype": "float32",
                "_FillValue": np.float32(9.96921e36),
            },
            "valid_layer_mask": {
                "zlib": True,
                "complevel": 4,
                "shuffle": True,
                "dtype": "int8",
                "_FillValue": np.int8(-127),
            },
            "gebco_water_depth": {
                "zlib": True,
                "complevel": 4,
                "shuffle": True,
                "dtype": "float32",
                "_FillValue": np.float32(9.96921e36),
            },
        }

        if FILES_2000[component].exists():
            FILES_2000[component].unlink()

        output_dataset.to_netcdf(
            FILES_2000[component],
            engine="h5netcdf",
            encoding=encoding,
        )

        print(
            f"Saved {component} 0–2000 m product:\n"
            f"{FILES_2000[component]}"
        )

        output_dataset.close()

    mask_source.close()
    raw.close()

    del (
        total_output,
        thermo_output,
        halo_output,
        valid_output,
        reference_temperature_native,
        reference_salinity_native,
        reference_temperature_k,
        reference_salinity,
        reference_potential_temperature_c,
        reference_absolute_salinity,
        reference_conservative_temperature,
        reference_specific_volume,
    )

    gc.collect()


if (
    REUSE_EXISTING_2000_PRODUCTS
    and have_all_2000_files
):
    print_header(
        "USING EXISTING CORRECTED EN4 0–2000 M PRODUCTS"
    )

    for component in COMPONENT_ORDER:
        print(
            f"{component:>7}: {FILES_2000[component]}"
        )

else:
    create_2000m_products()


# 9. OPEN THE SIX MONTHLY PRODUCTS
print_header(
    "OPENING 0–1000 M AND 0–2000 M MONTHLY PRODUCTS"
)

PREFERRED_VARIABLES = {
    "total": [
        "total_steric_height",
        "total_steric",
    ],
    "thermo": [
        "thermosteric_height",
        "thermosteric",
    ],
    "halo": [
        "halosteric_height",
        "halosteric",
    ],
}


def load_steric_product(file_path, component):
    dataset = open_dataset_safely(
        file_path,
        chunks={"time": 12},
    )

    dataset = standardize_coordinates(
        dataset,
        time_name=detect_coordinate(dataset, "time"),
        lat_name=detect_coordinate(dataset, "lat"),
        lon_name=detect_coordinate(dataset, "lon"),
    )

    dataset = subset_antarctic_latitudes(
        normalize_and_sort_longitude(dataset)
    )

    variable_name = choose_data_variable(
        dataset,
        preferred_names=PREFERRED_VARIABLES[component],
        excluded_names=[
            "valid_layer_mask",
            "gebco_water_depth",
        ],
    )

    data = convert_height_to_cm(
        dataset[variable_name].sel(
            time=slice(
                f"{YEAR_START}-01-01",
                f"{YEAR_END}-12-31",
            )
        )
    )

    return dataset, data


datasets_1000 = {}
datasets_2000 = {}
data_1000 = {}
data_2000 = {}

for component in COMPONENT_ORDER:
    datasets_1000[component], data_1000[component] = (
        load_steric_product(
            FILES_1000[component],
            component,
        )
    )

    datasets_2000[component], data_2000[component] = (
        load_steric_product(
            FILES_2000[component],
            component,
        )
    )


# 10. ALIGN TIME AND GRID
alignment_list = []

for component in COMPONENT_ORDER:
    alignment_list.append(data_1000[component])
    alignment_list.append(data_2000[component])

aligned = xr.align(
    *alignment_list,
    join="inner",
)

for component_index, component in enumerate(COMPONENT_ORDER):
    data_1000[component] = aligned[2 * component_index]
    data_2000[component] = aligned[2 * component_index + 1]

common_time = data_1000["total"]["time"]
common_latitude = data_1000["total"]["lat"]
common_longitude = data_1000["total"]["lon"]

number_of_months = len(common_time)

print(
    "Common period:",
    pd.Timestamp(common_time.values[0]),
    "to",
    pd.Timestamp(common_time.values[-1]),
)

print("Common months:", number_of_months)
print(
    "Common grid:",
    len(common_latitude),
    "latitudes ×",
    len(common_longitude),
    "longitudes",
)


# 11. FIXED COMMON MASK: ONLY CELLS DEEPER THAN 2000 M
if "gebco_water_depth" not in datasets_1000["total"].data_vars:
    raise KeyError(
        "The 0–1000 m total product does not contain "
        "'gebco_water_depth'."
    )

water_depth = (
    datasets_1000["total"]["gebco_water_depth"]
    .rename(
        {
            name: replacement
            for name, replacement in {
                "latitude": "lat",
                "longitude": "lon",
            }.items()
            if (
                name in datasets_1000["total"][
                    "gebco_water_depth"
                ].dims
                and replacement not in datasets_1000[
                    "total"
                ]["gebco_water_depth"].dims
            )
        }
    )
    .sel(
        lat=common_latitude,
        lon=common_longitude,
    )
    .load()
)

deep_water_mask = (
    np.isfinite(water_depth)
    & (water_depth >= DEPTH_2000)
)

common_mask = deep_water_mask.copy()

for component in COMPONENT_ORDER:
    valid_fraction_1000 = (
        data_1000[component]
        .notnull()
        .mean(dim="time")
    )

    valid_fraction_2000 = (
        data_2000[component]
        .notnull()
        .mean(dim="time")
    )

    common_mask = (
        common_mask
        & (
            valid_fraction_1000
            >= COMMON_VALID_FRACTION - 1.0e-10
        )
        & (
            valid_fraction_2000
            >= COMMON_VALID_FRACTION - 1.0e-10
        )
    )

common_mask = common_mask.compute()

common_cell_count = int(
    common_mask.sum().values
)

if common_cell_count == 0:
    raise RuntimeError(
        "The fixed common >2000 m mask contains no cells. "
        "Reduce COMMON_VALID_FRACTION only if scientifically justified."
    )

print_header(
    "FIXED COMMON MASK"
)

print(
    "GEBCO depth threshold     : >=2000 m"
)

print(
    "Monthly valid fraction    :",
    f"{COMMON_VALID_FRACTION:.2f}",
)

print(
    "Common spatial cells      :",
    f"{common_cell_count:,}",
)


# 12. CALCULATE THE DEPTH-SENSITIVITY FIELDS
def seasonal_anomalies(data_array):
    all_month_mean = data_array.mean(
        dim="time",
        skipna=True,
    )

    output = {}

    for season in SEASON_ORDER:
        seasonal_mean = (
            data_array
            .where(
                data_array["time"].dt.month.isin(
                    SEASON_MONTHS[season]
                ),
                drop=True,
            )
            .mean(
                dim="time",
                skipna=True,
            )
        )

        output[season] = (
            seasonal_mean
            - all_month_mean
        )

    return output


def select_comparison_field(data_array, mode):
    anomalies = seasonal_anomalies(data_array)

    if mode == "winter_minus_summer":
        field = (
            anomalies["Winter"]
            - anomalies["Summer"]
        )

    elif mode == "spring_anomaly":
        field = anomalies["Spring"]

    elif mode == "summer_anomaly":
        field = anomalies["Summer"]

    elif mode == "autumn_anomaly":
        field = anomalies["Autumn"]

    elif mode == "winter_anomaly":
        field = anomalies["Winter"]

    elif mode == "seasonal_peak_to_peak":
        stack = xr.concat(
            [
                anomalies[season]
                for season in SEASON_ORDER
            ],
            dim="season",
        )

        field = (
            stack.max(
                dim="season",
                skipna=True,
            )
            - stack.min(
                dim="season",
                skipna=True,
            )
        )

    elif mode == "seasonal_rms":
        stack = xr.concat(
            [
                anomalies[season]
                for season in SEASON_ORDER
            ],
            dim="season",
        )

        field = np.sqrt(
            (
                stack ** 2
            ).mean(
                dim="season",
                skipna=True,
            )
        )

    else:
        raise ValueError(mode)

    return field


fields_1000 = {}
fields_2000 = {}
difference_fields = {}

for component in COMPONENT_ORDER:
    fields_1000[component] = (
        select_comparison_field(
            data_1000[component],
            COMPARISON_MODE,
        )
        .where(common_mask)
        .compute()
    )

    fields_2000[component] = (
        select_comparison_field(
            data_2000[component],
            COMPARISON_MODE,
        )
        .where(common_mask)
        .compute()
    )

    difference_fields[component] = (
        fields_2000[component]
        - fields_1000[component]
    ).where(common_mask)


# 13. SAVE PROCESSED FIGURE S4 FIELDS
component_coordinate = xr.DataArray(
    COMPONENT_ORDER,
    dims="component",
    name="component",
)

fields_1000_stack = xr.concat(
    [
        fields_1000[component]
        for component in COMPONENT_ORDER
    ],
    dim=component_coordinate,
)

fields_2000_stack = xr.concat(
    [
        fields_2000[component]
        for component in COMPONENT_ORDER
    ],
    dim=component_coordinate,
)

difference_stack = xr.concat(
    [
        difference_fields[component]
        for component in COMPONENT_ORDER
    ],
    dim=component_coordinate,
)

output_fields = xr.Dataset(
    {
        "field_0_1000m": fields_1000_stack,
        "field_0_2000m": fields_2000_stack,
        "difference_0_2000m_minus_0_1000m": difference_stack,
        "common_deeper_than_2000m_mask": common_mask.astype(np.int8),
        "gebco_water_depth": water_depth,
    }
)

for variable_name in [
    "field_0_1000m",
    "field_0_2000m",
    "difference_0_2000m_minus_0_1000m",
]:
    output_fields[variable_name].attrs["units"] = "cm"

output_fields[
    "field_0_1000m"
].attrs["long_name"] = (
    f"0–1000 m EN4 {COMPARISON_MODE.replace('_', ' ')} field"
)

output_fields[
    "field_0_2000m"
].attrs["long_name"] = (
    f"0–2000 m EN4 {COMPARISON_MODE.replace('_', ' ')} field"
)

output_fields[
    "difference_0_2000m_minus_0_1000m"
].attrs["long_name"] = (
    f"0–2000 m minus 0–1000 m difference in "
    f"EN4 {COMPARISON_MODE.replace('_', ' ')}"
)

output_fields[
    "common_deeper_than_2000m_mask"
].attrs.update(
    {
        "long_name": (
            "fixed common cells deeper than 2000 m"
        ),
        "flag_values": np.array(
            [0, 1],
            dtype=np.int8,
        ),
        "flag_meanings": "excluded included",
    }
)

output_fields.attrs.update(
    {
        "title": "Figure S4 EN4 integration-depth sensitivity",
        "comparison_mode": COMPARISON_MODE,
        "difference_definition": DIFFERENCE_DEFINITION,
        "common_valid_fraction": COMMON_VALID_FRACTION,
        "common_cell_count": common_cell_count,
        "period": f"{YEAR_START}-01 through {YEAR_END}-12",
    }
)

output_fields.rename(
    {
        "lat": "latitude",
        "lon": "longitude",
    }
).to_netcdf(
    OUT_FIELDS_NC,
    engine="h5netcdf",
    encoding={
        "field_0_1000m": {
            "zlib": True,
            "complevel": 4,
        },
        "field_0_2000m": {
            "zlib": True,
            "complevel": 4,
        },
        "difference_0_2000m_minus_0_1000m": {
            "zlib": True,
            "complevel": 4,
        },
        "common_deeper_than_2000m_mask": {
            "zlib": True,
            "complevel": 4,
        },
        "gebco_water_depth": {
            "zlib": True,
            "complevel": 4,
        },
    },
)


# 14. SUMMARY TABLE
summary_rows = []

for component in COMPONENT_ORDER:
    for case_name, field in [
        ("0–1000 m", fields_1000[component]),
        ("0–2000 m", fields_2000[component]),
        ("Difference", difference_fields[component]),
    ]:
        values = np.asarray(
            field.values,
            dtype=float,
        )

        values = values[np.isfinite(values)]

        summary_rows.append(
            {
                "component": COMPONENT_LABELS[component],
                "depth_case": case_name,
                "comparison_mode": COMPARISON_MODE,
                "n_common_cells": int(values.size),
                "mean_cm": (
                    float(np.mean(values))
                    if values.size
                    else np.nan
                ),
                "median_cm": (
                    float(np.median(values))
                    if values.size
                    else np.nan
                ),
                "standard_deviation_cm": (
                    float(np.std(values))
                    if values.size
                    else np.nan
                ),
                "minimum_cm": (
                    float(np.min(values))
                    if values.size
                    else np.nan
                ),
                "maximum_cm": (
                    float(np.max(values))
                    if values.size
                    else np.nan
                ),
                "RMSE_cm": (
                    float(
                        np.sqrt(
                            np.mean(
                                values ** 2
                            )
                        )
                    )
                    if values.size
                    else np.nan
                ),
            }
        )

summary_table = pd.DataFrame(
    summary_rows
)

summary_table.to_csv(
    OUT_SUMMARY_CSV,
    index=False,
)


# 15. COLOUR LIMITS
signed_mode = COMPARISON_MODE not in {
    "seasonal_peak_to_peak",
    "seasonal_rms",
}

all_depth_fields = [
    fields_1000[component]
    for component in COMPONENT_ORDER
] + [
    fields_2000[component]
    for component in COMPONENT_ORDER
]

all_difference_fields = [
    difference_fields[component]
    for component in COMPONENT_ORDER
]

if signed_mode:
    depth_limit = robust_symmetric_limit(
        all_depth_fields
    )

    difference_limit = robust_symmetric_limit(
        all_difference_fields
    )

    depth_cmap = "RdBu_r"
    difference_cmap = "RdBu_r"

else:
    depth_limit = robust_positive_limit(
        all_depth_fields
    )

    difference_limit = robust_symmetric_limit(
        all_difference_fields
    )

    depth_cmap = "YlGnBu"
    difference_cmap = "RdBu_r"

print_header(
    "FIGURE COLOUR LIMITS"
)

if signed_mode:
    print(
        f"0–1000 / 0–2000 maps : ±{depth_limit:.2f} cm"
    )
else:
    print(
        f"0–1000 / 0–2000 maps : 0 to {depth_limit:.2f} cm"
    )

print(
    f"Difference maps        : ±{difference_limit:.2f} cm"
)


# 16. MAP-DRAWING HELPERS
def add_circular_boundary(axis):
    theta = np.linspace(
        0.0,
        2.0 * np.pi,
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
        + np.array(
            [0.5, 0.5]
        )
    )

    axis.set_boundary(
        circle,
        transform=axis.transAxes,
    )


def add_geographic_labels(axis):
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
    axis.set_extent(
        [
            -180.0,
            180.0,
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
        zorder=8,
    )

    axis.coastlines(
        resolution="110m",
        color=COAST_COLOR,
        linewidth=0.55,
        zorder=9,
    )

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
    panel_letter,
    cmap,
    value_minimum,
    value_maximum,
):
    style_map_axis(axis)

    longitude_2d, latitude_2d = coordinate_mesh(
        data_array
    )

    no_data = xr.where(
        common_mask,
        np.nan,
        1.0,
    )

    axis.pcolormesh(
        longitude_2d,
        latitude_2d,
        no_data.values,
        transform=ccrs.PlateCarree(),
        cmap=ListedColormap(
            [NO_DATA_COLOR]
        ),
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
        cmap=cmap,
        vmin=value_minimum,
        vmax=value_maximum,
        shading="auto",
        rasterized=True,
        zorder=4,
    )

    if DRAW_2000_M_CONTOUR:
        axis.contour(
            longitude_2d,
            latitude_2d,
            water_depth.values,
            levels=[
                BATHYMETRY_CONTOUR_LEVEL
            ],
            colors="black",
            linewidths=0.62,
            linestyles="-",
            alpha=0.86,
            transform=ccrs.PlateCarree(),
            zorder=10,
        )

    axis.text(
        PANEL_LETTER_X,
        PANEL_LETTER_Y,
        panel_letter,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=FONT_PANEL_LETTER,
        fontweight="bold",
        color="black",
        clip_on=False,
        zorder=50,
    )

    return image


# 17. BUILD THE 3 × 3 FIGURE
projection = ccrs.SouthPolarStereo(
    central_longitude=0.0
)

figure = plt.figure(
    figsize=FIGSIZE
)

grid = gridspec.GridSpec(
    nrows=3,
    ncols=3,
    figure=figure,
    left=0.085,
    right=0.985,
    top=0.920,
    bottom=0.150,
    wspace=0.055,
    hspace=0.205,
)

map_axes = np.empty(
    (3, 3),
    dtype=object,
)

last_images = [
    None,
    None,
    None,
]

for row_index, component in enumerate(
    COMPONENT_ORDER
):
    row_fields = [
        fields_1000[component],
        fields_2000[component],
        difference_fields[component],
    ]

    for column_index, field in enumerate(
        row_fields
    ):
        axis = figure.add_subplot(
            grid[row_index, column_index],
            projection=projection,
        )

        map_axes[
            row_index,
            column_index,
        ] = axis

        if column_index < 2:
            if signed_mode:
                value_minimum = -depth_limit
                value_maximum = depth_limit
            else:
                value_minimum = 0.0
                value_maximum = depth_limit

            cmap = depth_cmap

        else:
            value_minimum = -difference_limit
            value_maximum = difference_limit
            cmap = difference_cmap

        last_images[column_index] = plot_map_field(
            axis=axis,
            data_array=field,
            panel_letter=PANEL_LETTERS[
                row_index
            ][column_index],
            cmap=cmap,
            value_minimum=value_minimum,
            value_maximum=value_maximum,
        )


# 18. COLUMN HEADINGS
for column_index in range(3):
    position = map_axes[
        0,
        column_index,
    ].get_position()

    centre_x = (
        position.x0
        + position.x1
    ) / 2.0

    title_y = min(
        0.985,
        position.y1
        + COLUMN_TITLE_OFFSET,
    )

    figure.text(
        centre_x,
        title_y,
        COLUMN_LABELS[column_index],
        ha="center",
        va="top",
        fontsize=FONT_COLUMN_TITLE,
        fontweight="bold",
        color="black",
    )


# 19. LEFT-SIDE COMPONENT LABELS
for row_index, component in enumerate(
    COMPONENT_ORDER
):
    position = map_axes[
        row_index,
        0,
    ].get_position()

    figure.text(
        position.x0 - 0.045,
        position.y0
        + position.height / 2.0,
        COMPONENT_LABELS[component],
        rotation=90,
        rotation_mode="anchor",
        ha="center",
        va="center",
        fontsize=FONT_ROW_LABEL,
        fontweight="bold",
        color="black",
    )


# 20. BOTTOM COLOURBARS
bottom_positions = [
    map_axes[2, 0].get_position(),
    map_axes[2, 1].get_position(),
    map_axes[2, 2].get_position(),
]

colourbar_y = 0.083
colourbar_height = 0.014

colourbar_axes = [
    figure.add_axes(
        [
            bottom_positions[index].x0
            + 0.020,
            colourbar_y,
            bottom_positions[index].width
            - 0.040,
            colourbar_height,
        ]
    )
    for index in range(3)
]

for index, colourbar_axis in enumerate(
    colourbar_axes
):
    colourbar = figure.colorbar(
        last_images[index],
        cax=colourbar_axis,
        orientation="horizontal",
        extend="both",
    )

    colourbar.set_label(
        "(cm)",
        fontsize=FONT_COLORBAR,
        fontweight="bold",
        labelpad=3.0,
    )

    if index < 2:
        if signed_mode:
            ticks = np.linspace(
                -depth_limit,
                depth_limit,
                5,
            )
        else:
            ticks = np.linspace(
                0.0,
                depth_limit,
                5,
            )

    else:
        ticks = np.linspace(
            -difference_limit,
            difference_limit,
            5,
        )

    colourbar.set_ticks(ticks)

    colourbar.ax.tick_params(
        labelsize=FONT_COLORBAR_TICK,
        width=1.0,
        length=4.0,
        pad=2.0,
    )

    for label in colourbar.ax.get_xticklabels():
        label.set_fontweight("bold")


# 21. SAVE FIGURE AND TEXT SUMMARY
figure.savefig(
    OUT_PNG,
    dpi=SAVE_DPI,
    bbox_inches="tight",
    pad_inches=0.08,
    facecolor="white",
    edgecolor="none",
)

figure.savefig(
    OUT_PDF,
    bbox_inches="tight",
    pad_inches=0.08,
    facecolor="white",
    edgecolor="none",
)

plt.show()
plt.close(figure)

summary_lines = [
    "=" * 92,
    "FIGURE S4 — EN4 INTEGRATION-DEPTH SENSITIVITY",
    "=" * 92,
    "",
    f"Comparison mode : {COMPARISON_MODE}",
    f"Difference      : {DIFFERENCE_DEFINITION}",
    f"Common mask     : GEBCO depth >=2000 m",
    f"Valid fraction  : {COMMON_VALID_FRACTION:.2f}",
    f"Common cells    : {common_cell_count:,}",
    f"Study period    : {YEAR_START}-01 through {YEAR_END}-12",
    "",
    f"PNG             : {OUT_PNG}",
    f"PDF             : {OUT_PDF}",
    f"Processed fields: {OUT_FIELDS_NC}",
    f"Summary CSV     : {OUT_SUMMARY_CSV}",
    "",
    "Colour limits:",
    (
        f"0–1000 / 0–2000 maps: "
        + (
            f"±{depth_limit:.2f} cm"
            if signed_mode
            else f"0 to {depth_limit:.2f} cm"
        )
    ),
    f"Difference maps       : ±{difference_limit:.2f} cm",
    "",
    "Scientific rules:",
    "• Both depth columns use exactly the same common >2000 m cells.",
    "• Difference is calculated after applying the common mask.",
    "• 0–2000 m products use one fixed 2008–2025 grid-cell reference profile.",
    "• Monthly temperature and salinity must cover the full 0–2000 m layer.",
    "• Shallow cells are excluded and are never integrated only to the local seabed.",
    "• Halosteric height is total minus corrected thermosteric height.",
]

OUT_SUMMARY_TXT.write_text(
    "\n".join(summary_lines),
    encoding="utf-8",
)

print("\n".join(summary_lines))


# 22. CLOSE DATASETS
for dataset in datasets_1000.values():
    dataset.close()

for dataset in datasets_2000.values():
    dataset.close()

gc.collect()
