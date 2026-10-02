
"""
FIGURE 11 — FULL ORGANIZED LATEST WORKFLOW

FIGURE

Study-Period Trends During 2008-2025
Southern Ocean / 13 Antarctic Shelf Seas

SOURCE

Extracted and organized from:
    so_sealevel_paper_fig.py

VERSION SELECTION
The source contains two Figure 11 stages:

1. A complete trend-analysis workflow that calculates:
       - deseasonalized annual anomalies;
       - grid-cell Sen slopes;
       - Hamed-Rao modified Mann-Kendall significance;
       - Benjamini-Hochberg FDR correction;
       - 13-sea Sen slopes;
       - 3-year circular moving-block bootstrap confidence intervals.

   It saves:
       Fig11_gridcell_trend_statistics_2008_2025.nc
       Fig11_13sea_Sen_slopes_bootstrap_CI.csv
       Fig11_deseasonalized_annual_anomalies_common_grid.nc

2. A later FINAL CORRECTED plot-only workflow that reads those precomputed
   products and creates:
       Figure11_study_period_trends_FINAL_CORRECTED_1080dpi.png
       Figure11_study_period_trends_FINAL_CORRECTED.pdf
       Figure11_FINAL_CORRECTED_plot_summary.txt

This standalone file combines those two stages correctly:

    full trend computation
        -> saved analysis products
        -> latest FINAL CORRECTED plotting

The older first Figure 11 plotting block is intentionally excluded.

FIGURE STRUCTURE
Top row:
    (a) SLA trend map
    (b) Total steric trend map
    (c) Thermosteric trend map
    (d) Halosteric trend map

Bottom row:
    (e) 13-sea SLA Sen slopes
    (f) 13-sea total steric Sen slopes
    (g) 13-sea thermosteric Sen slopes
    (h) 13-sea halosteric Sen slopes

TREND METHOD
------------
1. Interpolate monthly SLA to the EN4 grid.
2. Apply one fixed SLA-EN4 common mask.
3. Remove each grid cell's 2008-2025 calendar-month climatology.
4. Average deseasonalized monthly anomalies by calendar year.
5. Estimate Sen's slope from annual anomalies.
6. Test significance with the Hamed-Rao modified Mann-Kendall test.
7. Correct grid-cell and sector p values using Benjamini-Hochberg FDR.
8. Estimate sector 95% confidence intervals using a 3-year circular
   moving-block bootstrap of residuals.

UNITS

cm decade^-1

ORGANIZED EXECUTION ORDER

1. Mount Google Drive.
2. Install required Python packages.
3. Import all required libraries.
4. Define source/output paths.
5. Define trend-analysis settings.
6. Define the 13 Antarctic sea sectors.
7. Validate and standardize input datasets.
8. Load EN4 total, thermosteric and halosteric monthly fields.
9. Interpolate monthly SLA to the EN4 grid.
10. Build the fixed SLA-EN4 common mask.
11. Remove monthly climatology and calculate annual anomalies.
12. Calculate grid-cell Sen slopes.
13. Apply Hamed-Rao modified Mann-Kendall tests.
14. Apply Benjamini-Hochberg FDR correction.
15. Calculate sector trends and bootstrap confidence intervals.
16. Save the 13-sea sector CSV.
17. Save the grid-cell trend-statistics NetCDF.
18. Save the deseasonalized annual-anomaly NetCDF.
19. Open the precomputed products with the latest final plot workflow.
20. Build the final corrected 2 x 4 Figure 11 layout.
21. Plot FDR-significant trend dots.
22. Plot 13-sea Sen slopes and confidence intervals.
23. Save the final vector PDF.
24. Save the final 1080-dpi PNG.
25. Save the final corrected plotting summary.

PRIMARY ANALYSIS OUTPUTS

/content/drive/MyDrive/SAM_Thesis/paper2/

    Fig11_gridcell_trend_statistics_2008_2025.nc
    Fig11_13sea_Sen_slopes_bootstrap_CI.csv
    Fig11_deseasonalized_annual_anomalies_common_grid.nc

LATEST FINAL FIGURE OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/

    Figure11_study_period_trends_FINAL_CORRECTED_1080dpi.png
    Figure11_study_period_trends_FINAL_CORRECTED.pdf
    Figure11_FINAL_CORRECTED_plot_summary.txt

FINAL-PLOT CORRECTIONS PRESERVED
The later source version explicitly includes:
    - no overall figure title;
    - separate clean map-colorbar row;
    - removal of the long overlapping statistics box;
    - separate map and sector explanation legends;
    - larger visible FDR-significance dots;
    - two-line bottom-panel titles;
    - stronger sector-point and confidence-interval colors;
    - dedicated non-overlapping space for labels and legends;
    - 1080-dpi PNG plus vector PDF output.

CODE CORRECTION APPLIED
Two standard longitude-normalization expressions in the supplied Figure 11
source had the `% 360.0` term accidentally commented out. They are restored
so longitude normalization correctly uses:

    ((longitude + 180.0) % 360.0) - 180.0

No Figure 11 scientific thresholds, trend definitions, statistical tests,
bootstrap settings, FDR settings, sector definitions, file names, or final
plot methodology are intentionally changed.

"""

# 1. MOUNT GOOGLE DRIVE

try:
    from google.colab import drive

    drive.mount("/content/drive")

except ImportError:
    print(
        "Google Drive mounting skipped because "
        "this is not Google Colab."
    )


# 2. INSTALL REQUIRED PACKAGES
import importlib.util
import subprocess
import sys

REQUIRED_PACKAGES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "xarray": "xarray",
    "netCDF4": "netCDF4",
    "h5netcdf": "h5netcdf",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "cartopy": "cartopy",
}

missing_packages = [
    pip_name
    for import_name, pip_name in REQUIRED_PACKAGES.items()
    if importlib.util.find_spec(import_name) is None
]

if missing_packages:
    print(
        "Installing missing packages:",
        ", ".join(missing_packages),
    )

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
    print(
        "All required packages are already installed."
    )


# 3. IMPORT LIBRARIES

import gc
import os
import warnings

from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from scipy.stats import norm, rankdata

import matplotlib

# Non-interactive backend is safer for high-resolution Colab export.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.colors as mcolors

from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from cartopy.util import add_cyclic_point

import cartopy.crs as ccrs
import cartopy.feature as cfeature

warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
)


# PART A — COMPLETE FIGURE 11 TREND COMPUTATION
# This section creates the precomputed grid, sector and annual-anomaly products
# required by the latest/final corrected Figure 11 renderer.

# FIGURE 11 — STUDY-PERIOD TRENDS DURING 2008–2025
#
# Structure
# Top row:
#   (a) SLA trend map
#   (b) Total steric trend map
#   (c) Thermosteric trend map
#   (d) Halosteric trend map
#
# Bottom row:
#   (e) 13-sea SLA Sen slopes
#   (f) 13-sea total steric Sen slopes
#   (g) 13-sea thermosteric Sen slopes
#   (h) 13-sea halosteric Sen slopes
#
# Trend method
# 1. Interpolate monthly SLA to the EN4 grid.
# 2. Apply one fixed SLA–EN4 common mask.
# 3. Remove each grid cell's 2008–2025 calendar-month climatology.
# 4. Average the deseasonalized monthly anomalies by calendar year.
# 5. Estimate Sen's slope from the annual anomalies.
# 6. Test significance with the Hamed–Rao modified Mann–Kendall test.
# 7. Correct grid-cell and sector p values using Benjamini–Hochberg FDR.
# 8. Estimate sector 95% confidence intervals using a 3-year
#    circular moving-block bootstrap of residuals.
#
# Units
# cm decade^-1
#
# Results section
# 3.9 Study-period trends and regional contrasts
#
# Output directory
# /content/drive/MyDrive/SAM_Thesis/paper2/


# 3. INPUT AND OUTPUT PATHS
SLA_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Data/"
    "SLA_Antarctic_monthly_2008_2025.nc"
)

TOTAL_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Processed/"
    "EN4_NetCDF_inventory/"
    "EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc"
)

THERMO_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Processed/"
    "EN4_NetCDF_inventory/"
    "EN4_thermosteric_0_1000m_monthly_2008_2025_SO.nc"
)

HALO_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Processed/"
    "EN4_NetCDF_inventory/"
    "EN4_halosteric_0_1000m_monthly_2008_2025_SO.nc"
)

OUTPUT_DIRECTORY = Path(
    "/content/drive/MyDrive/SAM_Thesis/paper2"
)

OUTPUT_DIRECTORY.mkdir(
    parents=True,
    exist_ok=True,
)

OUT_GRID_NC = (
    OUTPUT_DIRECTORY
    / "Fig11_gridcell_trend_statistics_2008_2025.nc"
)

OUT_SECTOR_CSV = (
    OUTPUT_DIRECTORY
    / "Fig11_13sea_Sen_slopes_bootstrap_CI.csv"
)

OUT_ANNUAL_NC = (
    OUTPUT_DIRECTORY
    / "Fig11_deseasonalized_annual_anomalies_common_grid.nc"
)

# 4. ANALYSIS SETTINGS
YEAR_START = 2008
YEAR_END = 2025

ANALYSIS_START = f"{YEAR_START}-01-01"
ANALYSIS_END = f"{YEAR_END}-12-31"

LAT_MIN = -90.0
LAT_MAX = -60.0

SLA_CORE_VALID_FRACTION = 0.70
EN4_CORE_VALID_FRACTION = 1.00

# At least this many valid monthly anomalies are required to form
# a grid-cell annual anomaly.
MIN_VALID_MONTHS_PER_YEAR = 8

# At least this many valid annual anomalies are required for trend.
MIN_VALID_YEARS = 12

FDR_ALPHA = 0.05
MK_ALPHA = 0.05

BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_BLOCK_LENGTH_YEARS = 3
BOOTSTRAP_SEED = 4217

DRAW_SECTOR_BOUNDARIES = True

# A compact physical figure permits 1080-dpi export without the
# enormous canvas created by a 15–18 inch figure.
FIGSIZE = (
    10.2,
    7.0,
)

SAVE_DPI = 1080


# 5. 13 ANTARCTIC SEA-SECTOR DEFINITIONS
# Use the same sector order and longitude limits used in the
# preceding sea-wise analysis figures.
#
# Longitudes use the -180 to 180 convention. A sector with
# lon_min > lon_max crosses the dateline.
SEA_SECTORS = OrderedDict([
    (
        "WED",
        {
            "name": "Weddell Sea",
            "lon_min": -60.0,
            "lon_max": -20.0,
        },
    ),
    (
        "KHV",
        {
            "name": "King Haakon VII Sea",
            "lon_min": -20.0,
            "lon_max": 10.0,
        },
    ),
    (
        "RLS",
        {
            "name": "Riiser-Larsen Sea",
            "lon_min": 10.0,
            "lon_max": 35.0,
        },
    ),
    (
        "LAZ",
        {
            "name": "Lazarev Sea",
            "lon_min": 35.0,
            "lon_max": 60.0,
        },
    ),
    (
        "COS",
        {
            "name": "Cosmonauts Sea",
            "lon_min": 60.0,
            "lon_max": 90.0,
        },
    ),
    (
        "COO",
        {
            "name": "Cooperation Sea",
            "lon_min": 90.0,
            "lon_max": 115.0,
        },
    ),
    (
        "DAV",
        {
            "name": "Davis Sea",
            "lon_min": 115.0,
            "lon_max": 130.0,
        },
    ),
    (
        "MAW",
        {
            "name": "Mawson Sea",
            "lon_min": 130.0,
            "lon_max": 150.0,
        },
    ),
    (
        "DUR",
        {
            "name": "D'Urville Sea",
            "lon_min": 150.0,
            "lon_max": 170.0,
        },
    ),
    (
        "SOM",
        {
            "name": "Somov Sea",
            "lon_min": 170.0,
            "lon_max": -160.0,
        },
    ),
    (
        "ROS",
        {
            "name": "Ross Sea",
            "lon_min": -160.0,
            "lon_max": -130.0,
        },
    ),
    (
        "AMU",
        {
            "name": "Amundsen Sea",
            "lon_min": -130.0,
            "lon_max": -100.0,
        },
    ),
    (
        "BEL",
        {
            "name": "Bellingshausen Sea",
            "lon_min": -100.0,
            "lon_max": -60.0,
        },
    ),
])

SEA_CODES = list(
    SEA_SECTORS.keys()
)

SEA_NAMES = [
    SEA_SECTORS[
        code
    ][
        "name"
    ]
    for code in SEA_CODES
]


# 6. FIGURE STYLE
FONT_FIGURE_TITLE = 14.0
FONT_MAP_TITLE = 10.8
FONT_FOREST_TITLE = 10.5
FONT_GEO_LABEL = 7.0
FONT_FOREST_TICK = 7.6
FONT_AXIS_LABEL = 8.8
FONT_COLORBAR = 8.4
FONT_NOTE = 7.6
FONT_LEGEND = 7.5

LAND_COLOR = "0.84"
COAST_COLOR = "0.25"
NO_DATA_COLOR = "0.86"
SECTOR_LINE_COLOR = "white"

MERIDIANS_FOR_GRID = [
    -180,
    -90,
    0,
    90,
    180,
]

PARALLELS_FOR_GRID = [
    -60,
    -70,
    -80,
]

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "mathtext.default": "regular",
    }
)


# 7. FILE AND DATASET HELPERS
for required_file in [
    SLA_FILE,
    TOTAL_FILE,
    THERMO_FILE,
    HALO_FILE,
]:

    if not os.path.exists(
        required_file
    ):

        raise FileNotFoundError(
            f"Required file was not found:\n"
            f"{required_file}"
        )


def safe_is_datetime(
    dtype,
):

    try:

        return np.issubdtype(
            dtype,
            np.datetime64,
        )

    except TypeError:

        return False


def open_dataset_safely(
    file_path,
):
    """
    Open without Dask. The SLA file is read month-by-month, while
    the EN4 fields are already small on the 24 × 360 target grid.
    """
    attempts = []

    for engine in [
        None,
        "netcdf4",
        "h5netcdf",
        "scipy",
    ]:

        try:

            kwargs = {
                "decode_times": True,
                "mask_and_scale": True,
                "cache": False,
            }

            if engine is not None:

                kwargs[
                    "engine"
                ] = engine

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
                f"with engine={engine_name}"
            )

            return dataset

        except Exception as error:

            attempts.append(
                f"engine={engine}: {error}"
            )

    raise RuntimeError(
        f"Could not open:\n"
        f"{file_path}\n\n"
        + "\n".join(
            attempts
        )
    )


def detect_coordinate(
    dataset,
    coordinate_type,
):
    """
    Detect time, latitude and longitude coordinate names.
    """
    aliases = {
        "time": [
            "time",
            "date",
            "datetime",
            "valid_time",
            "month",
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

    candidates = list(
        dataset.coords
    ) + [
        name
        for name in dataset.variables
        if name not in dataset.coords
    ]

    for name in candidates:

        variable = dataset[
            name
        ]

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
                lower_name
                in aliases["time"]
                or standard_name
                == "time"
                or axis == "T"
                or " since "
                in units
                or safe_is_datetime(
                    variable.dtype
                )
            ):

                return name

        elif coordinate_type == "lat":

            if (
                lower_name
                in aliases["lat"]
                or standard_name
                == "latitude"
                or axis == "Y"
                or "degrees_north"
                in units
            ):

                return name

        elif coordinate_type == "lon":

            if (
                lower_name
                in aliases["lon"]
                or standard_name
                == "longitude"
                or axis == "X"
                or "degrees_east"
                in units
            ):

                return name

    raise KeyError(
        f"Could not detect "
        f"{coordinate_type}."
    )


def standardize_dataset(
    dataset,
):
    """
    Standardize coordinate names and longitude convention.
    """
    time_name = detect_coordinate(
        dataset,
        "time",
    )

    lat_name = detect_coordinate(
        dataset,
        "lat",
    )

    lon_name = detect_coordinate(
        dataset,
        "lon",
    )

    rename_mapping = {}

    if time_name != "time":

        rename_mapping[
            time_name
        ] = "time"

    if lat_name != "lat":

        rename_mapping[
            lat_name
        ] = "lat"

    if lon_name != "lon":

        rename_mapping[
            lon_name
        ] = "lon"

    if rename_mapping:

        dataset = dataset.rename(
            rename_mapping
        )

    normalized_longitude = (
        (
            dataset[
                "lon"
            ].astype(float)
            + 180.0
        )
         % 360.0
    ) - 180.0

    dataset = dataset.assign_coords(
        lon=normalized_longitude
    )

    longitude_values = np.asarray(
        dataset[
            "lon"
        ].values
    )

    _, unique_indices = np.unique(
        longitude_values,
        return_index=True,
    )

    dataset = dataset.isel(
        lon=np.sort(
            unique_indices
        )
    )

    dataset = dataset.sortby(
        "lon"
    )

    dataset = dataset.sortby(
        "lat"
    )

    dataset = dataset.sel(
        time=slice(
            ANALYSIS_START,
            ANALYSIS_END,
        )
    )

    return dataset


def choose_variable(
    dataset,
    exact_names,
    contains_names,
):
    """
    Select the requested science variable.
    """
    excluded_names = {
        "valid_layer_mask",
        "gebco_water_depth",
        "surface_salinity_reference",
        "time_bnds",
        "depth_bnds",
    }

    available = [
        variable_name
        for variable_name
        in dataset.data_vars
        if variable_name
        not in excluded_names
    ]

    for candidate in exact_names:

        if candidate in available:

            return candidate

    for token in contains_names:

        for variable_name in available:

            if (
                token.lower()
                in variable_name.lower()
            ):

                return variable_name

    raise KeyError(
        "Could not find the requested science variable.\n"
        f"Available variables: {available}"
    )


def unit_factor_to_cm(
    data_array,
):
    """
    Determine a multiplication factor without loading the full array.
    """
    units = str(
        data_array.attrs.get(
            "units",
            "",
        )
    ).strip().lower()

    if units in [
        "m",
        "meter",
        "meters",
        "metre",
        "metres",
    ]:

        return 100.0

    if units in [
        "cm",
        "centimeter",
        "centimeters",
        "centimetre",
        "centimetres",
    ]:

        return 1.0

    sample = np.asarray(
        data_array.isel(
            time=0
        ).values,
        dtype=np.float32,
    )

    finite = sample[
        np.isfinite(
            sample
        )
    ]

    if (
        finite.size > 0
        and np.nanpercentile(
            np.abs(
                finite
            ),
            99.0,
        )
        < 1.0
    ):

        return 100.0

    return 1.0


def prepare_monthly_field(
    field,
):
    """
    Return a 2-D lat × lon monthly field.
    """
    field = field.squeeze(
        drop=True
    )

    extra_dimensions = [
        dimension
        for dimension
        in field.dims
        if dimension
        not in [
            "lat",
            "lon",
        ]
    ]

    if extra_dimensions:

        raise ValueError(
            "Unexpected dimensions remain in a monthly field: "
            f"{extra_dimensions}"
        )

    return field.transpose(
        "lat",
        "lon",
    )


# 8. LOAD EN4 TARGET GRID AND MONTHLY VARIABLES
total_dataset = standardize_dataset(
    open_dataset_safely(
        TOTAL_FILE
    )
)

total_variable = choose_variable(
    total_dataset,
    exact_names=[
        "total_steric_height",
    ],
    contains_names=[
        "total_steric",
    ],
)

target_latitude = np.asarray(
    total_dataset[
        "lat"
    ].values,
    dtype=np.float64,
)

target_longitude = np.asarray(
    total_dataset[
        "lon"
    ].values,
    dtype=np.float64,
)

target_times = pd.DatetimeIndex(
    pd.to_datetime(
        total_dataset[
            "time"
        ].values
    )
)

total_factor = unit_factor_to_cm(
    total_dataset[
        total_variable
    ]
)

total_monthly = (
    np.asarray(
        total_dataset[
            total_variable
        ].transpose(
            "time",
            "lat",
            "lon",
        ).values,
        dtype=np.float32,
    )
    * np.float32(
        total_factor
    )
)

if "valid_layer_mask" in total_dataset:

    total_layer_mask = np.asarray(
        total_dataset[
            "valid_layer_mask"
        ].transpose(
            "time",
            "lat",
            "lon",
        ).values,
        dtype=np.int8,
    )

else:

    total_layer_mask = np.isfinite(
        total_monthly
    ).astype(
        np.int8
    )

total_dataset.close()

del total_dataset

gc.collect()


def load_en4_variable(
    file_path,
    exact_names,
    contains_names,
    target_times,
    target_latitude,
    target_longitude,
):
    """
    Load one EN4 variable on the target grid and common monthly period.
    """
    dataset = standardize_dataset(
        open_dataset_safely(
            file_path
        )
    )

    variable_name = choose_variable(
        dataset,
        exact_names=exact_names,
        contains_names=contains_names,
    )

    factor = unit_factor_to_cm(
        dataset[
            variable_name
        ]
    )

    selected = dataset[
        variable_name
    ].sel(
        time=target_times
    )

    grids_match = (
        selected.sizes.get(
            "lat"
        )
        == target_latitude.size
        and selected.sizes.get(
            "lon"
        )
        == target_longitude.size
        and np.allclose(
            selected[
                "lat"
            ].values,
            target_latitude,
        )
        and np.allclose(
            selected[
                "lon"
            ].values,
            target_longitude,
        )
    )

    if not grids_match:

        selected = selected.interp(
            lat=xr.DataArray(
                target_latitude,
                dims="lat",
                coords={
                    "lat": target_latitude,
                },
            ),
            lon=xr.DataArray(
                target_longitude,
                dims="lon",
                coords={
                    "lon": target_longitude,
                },
            ),
            method="linear",
        )

    values = (
        np.asarray(
            selected.transpose(
                "time",
                "lat",
                "lon",
            ).values,
            dtype=np.float32,
        )
        * np.float32(
            factor
        )
    )

    dataset.close()

    del dataset
    del selected

    gc.collect()

    return (
        values,
        variable_name,
    )


thermosteric_monthly, thermo_variable = (
    load_en4_variable(
        THERMO_FILE,
        exact_names=[
            "thermosteric_height",
        ],
        contains_names=[
            "thermosteric",
        ],
        target_times=target_times,
        target_latitude=target_latitude,
        target_longitude=target_longitude,
    )
)

halosteric_monthly, halo_variable = (
    load_en4_variable(
        HALO_FILE,
        exact_names=[
            "halosteric_height",
        ],
        contains_names=[
            "halosteric",
        ],
        target_times=target_times,
        target_latitude=target_latitude,
        target_longitude=target_longitude,
    )
)


# 9. INTERPOLATE SLA MONTH-BY-MONTH TO THE EN4 GRID
sla_dataset = standardize_dataset(
    open_dataset_safely(
        SLA_FILE
    )
)

sla_variable = choose_variable(
    sla_dataset,
    exact_names=[
        "sla",
    ],
    contains_names=[
        "sea_level_anomaly",
        "sla",
    ],
)

sla_factor = unit_factor_to_cm(
    sla_dataset[
        sla_variable
    ]
)

sla_times = pd.DatetimeIndex(
    pd.to_datetime(
        sla_dataset[
            "time"
        ].values
    )
)

common_times = target_times.intersection(
    sla_times
)

if len(
    common_times
) != len(
    target_times
):

    raise ValueError(
        "The SLA and EN4 monthly sequences do not match exactly.\n"
        f"EN4 months: {len(target_times)}\n"
        f"Common months: {len(common_times)}"
    )

sla_monthly = np.full(
    (
        len(target_times),
        target_latitude.size,
        target_longitude.size,
    ),
    np.nan,
    dtype=np.float32,
)

print(
    "\nInterpolating SLA month-by-month to the EN4 grid..."
)

for time_index, timestamp in enumerate(
    target_times
):

    monthly_field = prepare_monthly_field(
        sla_dataset[
            sla_variable
        ].sel(
            time=timestamp
        )
    )

    monthly_field = monthly_field.interp(
        lat=xr.DataArray(
            target_latitude,
            dims="lat",
            coords={
                "lat": target_latitude,
            },
        ),
        lon=xr.DataArray(
            target_longitude,
            dims="lon",
            coords={
                "lon": target_longitude,
            },
        ),
        method="linear",
    )

    sla_monthly[
        time_index
    ] = (
        np.asarray(
            monthly_field.values,
            dtype=np.float32,
        )
        * np.float32(
            sla_factor
        )
    )

    del monthly_field

    if (
        (time_index + 1) % 12 == 0
        or time_index
        == len(target_times) - 1
    ):

        print(
            f"  [{time_index + 1:03d}/"
            f"{len(target_times):03d}] "
            f"{timestamp:%Y-%m}"
        )

        gc.collect()

sla_dataset.close()

del sla_dataset

gc.collect()


# 10. FIXED COMMON SLA–EN4 MASK
number_of_months = len(
    target_times
)

sla_valid_fraction = (
    np.isfinite(
        sla_monthly
    ).sum(
        axis=0
    )
    / number_of_months
)

total_valid_fraction = (
    (
        np.isfinite(
            total_monthly
        )
        & (
            total_layer_mask
            == 1
        )
    ).sum(
        axis=0
    )
    / number_of_months
)

thermo_valid_fraction = (
    np.isfinite(
        thermosteric_monthly
    ).sum(
        axis=0
    )
    / number_of_months
)

halo_valid_fraction = (
    np.isfinite(
        halosteric_monthly
    ).sum(
        axis=0
    )
    / number_of_months
)

common_mask = (
    (
        sla_valid_fraction
        >= SLA_CORE_VALID_FRACTION
    )
    & (
        total_valid_fraction
        >= (
            EN4_CORE_VALID_FRACTION
            - 1.0e-10
        )
    )
    & (
        thermo_valid_fraction
        >= (
            EN4_CORE_VALID_FRACTION
            - 1.0e-10
        )
    )
    & (
        halo_valid_fraction
        >= (
            EN4_CORE_VALID_FRACTION
            - 1.0e-10
        )
    )
)

common_cell_count = int(
    common_mask.sum()
)

if common_cell_count == 0:

    raise RuntimeError(
        "The fixed SLA–EN4 common mask contains no cells."
    )

print(
    "\nFixed common mask cells:",
    f"{common_cell_count:,}",
)


# 11. REMOVE MONTHLY CLIMATOLOGY AND FORM ANNUAL ANOMALIES
calendar_months = np.asarray(
    target_times.month,
    dtype=np.int16,
)

calendar_years = np.asarray(
    target_times.year,
    dtype=np.int16,
)

analysis_years = np.arange(
    YEAR_START,
    YEAR_END + 1,
    dtype=np.int16,
)


def deseasonalize_and_annualize(
    monthly_values,
    calendar_months,
    calendar_years,
    analysis_years,
    fixed_mask,
):
    """
    Remove the grid-cell calendar-month climatology, then average
    the deseasonalized monthly anomalies by calendar year.
    """
    monthly_values = np.asarray(
        monthly_values,
        dtype=np.float32,
    )

    monthly_climatology = np.full(
        (
            12,
            monthly_values.shape[1],
            monthly_values.shape[2],
        ),
        np.nan,
        dtype=np.float32,
    )

    for month in range(
        1,
        13,
    ):

        month_indices = (
            calendar_months
            == month
        )

        monthly_climatology[
            month - 1
        ] = np.nanmean(
            monthly_values[
                month_indices
            ],
            axis=0,
        ).astype(
            np.float32
        )

    deseasonalized = np.full_like(
        monthly_values,
        np.nan,
        dtype=np.float32,
    )

    for time_index, month in enumerate(
        calendar_months
    ):

        deseasonalized[
            time_index
        ] = (
            monthly_values[
                time_index
            ]
            - monthly_climatology[
                month - 1
            ]
        )

    deseasonalized[
        :,
        ~fixed_mask,
    ] = np.nan

    annual_values = np.full(
        (
            analysis_years.size,
            monthly_values.shape[1],
            monthly_values.shape[2],
        ),
        np.nan,
        dtype=np.float32,
    )

    for year_index, year in enumerate(
        analysis_years
    ):

        year_indices = (
            calendar_years
            == year
        )

        year_block = deseasonalized[
            year_indices
        ]

        valid_month_count = np.isfinite(
            year_block
        ).sum(
            axis=0
        )

        year_mean = np.nanmean(
            year_block,
            axis=0,
        ).astype(
            np.float32
        )

        year_mean[
            valid_month_count
            < MIN_VALID_MONTHS_PER_YEAR
        ] = np.nan

        annual_values[
            year_index
        ] = year_mean

    return (
        annual_values,
        monthly_climatology,
    )


print(
    "\nRemoving monthly climatology and forming annual anomalies..."
)

sla_annual, sla_monthly_climatology = (
    deseasonalize_and_annualize(
        sla_monthly,
        calendar_months,
        calendar_years,
        analysis_years,
        common_mask,
    )
)

total_annual, total_monthly_climatology = (
    deseasonalize_and_annualize(
        total_monthly,
        calendar_months,
        calendar_years,
        analysis_years,
        common_mask,
    )
)

thermo_annual, thermo_monthly_climatology = (
    deseasonalize_and_annualize(
        thermosteric_monthly,
        calendar_months,
        calendar_years,
        analysis_years,
        common_mask,
    )
)

halo_annual, halo_monthly_climatology = (
    deseasonalize_and_annualize(
        halosteric_monthly,
        calendar_months,
        calendar_years,
        analysis_years,
        common_mask,
    )
)

# Release monthly arrays before the trend calculations.
del sla_monthly
del total_monthly
del thermosteric_monthly
del halosteric_monthly
del total_layer_mask

gc.collect()


# 12. SEN SLOPE AND HAMED–RAO MODIFIED MANN–KENDALL
def sen_slope(
    values,
    times,
):
    """
    Median pairwise Sen slope.
    """
    values = np.asarray(
        values,
        dtype=np.float64,
    )

    times = np.asarray(
        times,
        dtype=np.float64,
    )

    valid = (
        np.isfinite(
            values
        )
        & np.isfinite(
            times
        )
    )

    values = values[
        valid
    ]

    times = times[
        valid
    ]

    number_of_values = values.size

    if number_of_values < 2:

        return np.nan

    upper_i, upper_j = np.triu_indices(
        number_of_values,
        k=1,
    )

    time_difference = (
        times[
            upper_j
        ]
        - times[
            upper_i
        ]
    )

    value_difference = (
        values[
            upper_j
        ]
        - values[
            upper_i
        ]
    )

    valid_pairs = (
        time_difference
        != 0.0
    )

    slopes = (
        value_difference[
            valid_pairs
        ]
        / time_difference[
            valid_pairs
        ]
    )

    if slopes.size == 0:

        return np.nan

    return float(
        np.nanmedian(
            slopes
        )
    )


def mann_kendall_score_and_variance(
    values,
):
    """
    Mann–Kendall S score and tie-corrected variance.
    """
    values = np.asarray(
        values,
        dtype=np.float64,
    )

    number_of_values = values.size

    score = 0.0

    for index in range(
        number_of_values - 1
    ):

        score += np.sign(
            values[
                index + 1:
            ]
            - values[
                index
            ]
        ).sum()

    _, tie_counts = np.unique(
        values,
        return_counts=True,
    )

    tie_term = np.sum(
        tie_counts
        * (
            tie_counts
            - 1
        )
        * (
            2
            * tie_counts
            + 5
        )
    )

    variance = (
        number_of_values
        * (
            number_of_values
            - 1
        )
        * (
            2
            * number_of_values
            + 5
        )
        - tie_term
    ) / 18.0

    return (
        float(
            score
        ),
        float(
            variance
        ),
    )


def autocorrelation(
    values,
    lag,
):
    """
    Pearson autocorrelation at one lag.
    """
    values = np.asarray(
        values,
        dtype=np.float64,
    )

    if lag <= 0 or lag >= values.size:

        return np.nan

    first = values[
        :-lag
    ]

    second = values[
        lag:
    ]

    first = first - np.mean(
        first
    )

    second = second - np.mean(
        second
    )

    denominator = np.sqrt(
        np.sum(
            first ** 2
        )
        * np.sum(
            second ** 2
        )
    )

    if denominator <= 0.0:

        return 0.0

    return float(
        np.sum(
            first
            * second
        )
        / denominator
    )


def hamed_rao_modified_mk(
    values,
    times,
    alpha=MK_ALPHA,
):
    """
    Hamed–Rao variance-corrected Mann–Kendall test.

    The autocorrelation correction is calculated from ranks of the
    Sen-slope-detrended series. Only rank autocorrelations outside
    the approximate 95% white-noise bounds contribute to the
    correction factor.
    """
    values = np.asarray(
        values,
        dtype=np.float64,
    )

    times = np.asarray(
        times,
        dtype=np.float64,
    )

    valid = (
        np.isfinite(
            values
        )
        & np.isfinite(
            times
        )
    )

    values = values[
        valid
    ]

    times = times[
        valid
    ]

    number_of_values = values.size

    if number_of_values < 4:

        return (
            np.nan,
            np.nan,
            np.nan,
        )

    score, variance = (
        mann_kendall_score_and_variance(
            values
        )
    )

    if variance <= 0.0:

        return (
            0.0,
            1.0,
            1.0,
        )

    slope = sen_slope(
        values,
        times,
    )

    intercept = np.nanmedian(
        values
        - slope
        * times
    )

    detrended = (
        values
        - (
            intercept
            + slope
            * times
        )
    )

    ranked = rankdata(
        detrended,
        method="average",
    )

    confidence_limit = (
        norm.ppf(
            1.0
            - alpha
            / 2.0
        )
        / np.sqrt(
            number_of_values
        )
    )

    weighted_autocorrelation_sum = 0.0

    for lag in range(
        1,
        number_of_values,
    ):

        rank_autocorrelation = (
            autocorrelation(
                ranked,
                lag,
            )
        )

        if (
            np.isfinite(
                rank_autocorrelation
            )
            and abs(
                rank_autocorrelation
            )
            > confidence_limit
        ):

            weighted_autocorrelation_sum += (
                (
                    number_of_values
                    - lag
                )
                * (
                    number_of_values
                    - lag
                    - 1
                )
                * (
                    number_of_values
                    - lag
                    - 2
                )
                * rank_autocorrelation
            )

    denominator = (
        number_of_values
        * (
            number_of_values
            - 1
        )
        * (
            number_of_values
            - 2
        )
    )

    correction_factor = (
        1.0
        + (
            2.0
            * weighted_autocorrelation_sum
            / denominator
        )
    )

    correction_factor = max(
        correction_factor,
        1.0e-6,
    )

    modified_variance = (
        variance
        * correction_factor
    )

    if score > 0.0:

        z_score = (
            score
            - 1.0
        ) / np.sqrt(
            modified_variance
        )

    elif score < 0.0:

        z_score = (
            score
            + 1.0
        ) / np.sqrt(
            modified_variance
        )

    else:

        z_score = 0.0

    p_value = (
        2.0
        * (
            1.0
            - norm.cdf(
                abs(
                    z_score
                )
            )
        )
    )

    return (
        float(
            z_score
        ),
        float(
            p_value
        ),
        float(
            correction_factor
        ),
    )


def benjamini_hochberg_qvalues(
    p_values,
):
    """
    Benjamini–Hochberg FDR-adjusted q values.
    """
    p_values = np.asarray(
        p_values,
        dtype=np.float64,
    )

    q_values = np.full_like(
        p_values,
        np.nan,
        dtype=np.float64,
    )

    finite_mask = np.isfinite(
        p_values
    )

    finite_p = p_values[
        finite_mask
    ]

    number_of_tests = finite_p.size

    if number_of_tests == 0:

        return q_values

    order = np.argsort(
        finite_p
    )

    sorted_p = finite_p[
        order
    ]

    ranks = np.arange(
        1,
        number_of_tests + 1,
        dtype=np.float64,
    )

    sorted_q = (
        sorted_p
        * number_of_tests
        / ranks
    )

    sorted_q = np.minimum.accumulate(
        sorted_q[
            ::-1
        ]
    )[
        ::-1
    ]

    sorted_q = np.clip(
        sorted_q,
        0.0,
        1.0,
    )

    finite_q = np.empty_like(
        sorted_q
    )

    finite_q[
        order
    ] = sorted_q

    q_values[
        finite_mask
    ] = finite_q

    return q_values


def grid_trend_statistics(
    annual_values,
    years,
    fixed_mask,
    variable_label,
):
    """
    Calculate grid-cell Sen slope and Hamed–Rao p values.

    Slopes are returned in cm decade^-1.
    """
    latitude_count = annual_values.shape[1]
    longitude_count = annual_values.shape[2]

    slope_map = np.full(
        (
            latitude_count,
            longitude_count,
        ),
        np.nan,
        dtype=np.float32,
    )

    p_map = np.full(
        (
            latitude_count,
            longitude_count,
        ),
        np.nan,
        dtype=np.float32,
    )

    correction_map = np.full(
        (
            latitude_count,
            longitude_count,
        ),
        np.nan,
        dtype=np.float32,
    )

    valid_positions = np.argwhere(
        fixed_mask
    )

    print(
        f"\nCalculating {variable_label} grid trends "
        f"for {valid_positions.shape[0]:,} cells..."
    )

    for position_index, (
        latitude_index,
        longitude_index,
    ) in enumerate(
        valid_positions
    ):

        series = annual_values[
            :,
            latitude_index,
            longitude_index,
        ]

        valid_years = np.isfinite(
            series
        )

        if valid_years.sum() < MIN_VALID_YEARS:

            continue

        selected_values = series[
            valid_years
        ]

        selected_years = years[
            valid_years
        ].astype(
            np.float64
        )

        slope_per_year = sen_slope(
            selected_values,
            selected_years,
        )

        _, p_value, correction_factor = (
            hamed_rao_modified_mk(
                selected_values,
                selected_years,
            )
        )

        slope_map[
            latitude_index,
            longitude_index,
        ] = np.float32(
            slope_per_year
            * 10.0
        )

        p_map[
            latitude_index,
            longitude_index,
        ] = np.float32(
            p_value
        )

        correction_map[
            latitude_index,
            longitude_index,
        ] = np.float32(
            correction_factor
        )

        if (
            (position_index + 1) % 500 == 0
            or position_index
            == valid_positions.shape[0] - 1
        ):

            print(
                f"  [{position_index + 1:4d}/"
                f"{valid_positions.shape[0]:4d}]"
            )

    q_map = np.full_like(
        p_map,
        np.nan,
        dtype=np.float32,
    )

    q_values = benjamini_hochberg_qvalues(
        p_map[
            fixed_mask
        ]
    )

    q_map[
        fixed_mask
    ] = q_values.astype(
        np.float32
    )

    significant_map = (
        np.isfinite(
            q_map
        )
        & (
            q_map
            < FDR_ALPHA
        )
    )

    return {
        "slope": slope_map,
        "p": p_map,
        "q": q_map,
        "significant": significant_map,
        "hamed_rao_factor": correction_map,
    }


sla_grid_statistics = grid_trend_statistics(
    sla_annual,
    analysis_years.astype(
        np.float64
    ),
    common_mask,
    "SLA",
)

total_grid_statistics = grid_trend_statistics(
    total_annual,
    analysis_years.astype(
        np.float64
    ),
    common_mask,
    "total steric",
)

thermo_grid_statistics = grid_trend_statistics(
    thermo_annual,
    analysis_years.astype(
        np.float64
    ),
    common_mask,
    "thermosteric",
)

halo_grid_statistics = grid_trend_statistics(
    halo_annual,
    analysis_years.astype(
        np.float64
    ),
    common_mask,
    "halosteric",
)


# 13. SECTOR MASKS, WEIGHTS AND BLOCK-BOOTSTRAP CIs
longitude_2d, latitude_2d = np.meshgrid(
    target_longitude,
    target_latitude,
)

area_weights = np.cos(
    np.deg2rad(
        latitude_2d
    )
).astype(
    np.float64
)


def longitude_sector_mask(
    longitude_values,
    minimum_longitude,
    maximum_longitude,
):
    """
    One-dimensional sector mask, including dateline-crossing sectors.
    """
    if (
        minimum_longitude
        <= maximum_longitude
    ):

        return (
            (
                longitude_values
                >= minimum_longitude
            )
            & (
                longitude_values
                < maximum_longitude
            )
        )

    return (
        (
            longitude_values
            >= minimum_longitude
        )
        | (
            longitude_values
            < maximum_longitude
        )
    )


def sector_mask(
    sector_information,
):
    """
    Two-dimensional sector mask on the common EN4 grid.
    """
    longitude_mask = longitude_sector_mask(
        target_longitude,
        sector_information[
            "lon_min"
        ],
        sector_information[
            "lon_max"
        ],
    )

    return (
        np.broadcast_to(
            longitude_mask[
                None,
                :
            ],
            common_mask.shape,
        )
        & common_mask
    )


def area_weighted_sector_series(
    annual_values,
    mask,
):
    """
    Area-weighted annual sector mean.
    """
    number_of_years = annual_values.shape[0]

    output = np.full(
        number_of_years,
        np.nan,
        dtype=np.float64,
    )

    fixed_weights = np.where(
        mask,
        area_weights,
        np.nan,
    )

    for year_index in range(
        number_of_years
    ):

        field = annual_values[
            year_index
        ]

        valid = (
            mask
            & np.isfinite(
                field
            )
        )

        if not np.any(
            valid
        ):

            continue

        denominator = np.nansum(
            fixed_weights[
                valid
            ]
        )

        if denominator <= 0.0:

            continue

        output[
            year_index
        ] = (
            np.nansum(
                field[
                    valid
                ]
                * fixed_weights[
                    valid
                ]
            )
            / denominator
        )

    return output


def moving_block_bootstrap_ci(
    annual_values,
    years,
    block_length=BOOTSTRAP_BLOCK_LENGTH_YEARS,
    number_of_bootstraps=BOOTSTRAP_REPLICATES,
    random_seed=BOOTSTRAP_SEED,
):
    """
    Residual circular moving-block bootstrap confidence interval.

    A Sen trend is fitted first. Residuals are sampled in contiguous
    circular blocks and added to the fitted trend. Each bootstrap
    realization is assigned the original sequential year positions.
    """
    annual_values = np.asarray(
        annual_values,
        dtype=np.float64,
    )

    years = np.asarray(
        years,
        dtype=np.float64,
    )

    valid = (
        np.isfinite(
            annual_values
        )
        & np.isfinite(
            years
        )
    )

    annual_values = annual_values[
        valid
    ]

    years = years[
        valid
    ]

    number_of_years = annual_values.size

    if number_of_years < MIN_VALID_YEARS:

        return (
            np.nan,
            np.nan,
        )

    point_slope = sen_slope(
        annual_values,
        years,
    )

    intercept = np.nanmedian(
        annual_values
        - point_slope
        * years
    )

    fitted = (
        intercept
        + point_slope
        * years
    )

    residuals = (
        annual_values
        - fitted
    )

    random_generator = np.random.default_rng(
        random_seed
    )

    bootstrap_slopes = np.full(
        number_of_bootstraps,
        np.nan,
        dtype=np.float64,
    )

    relative_time = np.arange(
        number_of_years,
        dtype=np.float64,
    )

    for bootstrap_index in range(
        number_of_bootstraps
    ):

        sampled_indices = []

        while len(
            sampled_indices
        ) < number_of_years:

            block_start = int(
                random_generator.integers(
                    0,
                    number_of_years,
                )
            )

            sampled_indices.extend(
                [
                    (
                        block_start
                        + offset
                    )
#                     % number_of_years
                    for offset in range(
                        block_length
                    )
                ]
            )

        sampled_indices = np.asarray(
            sampled_indices[
                :number_of_years
            ],
            dtype=np.int64,
        )

        bootstrap_values = (
            fitted
            + residuals[
                sampled_indices
            ]
        )

        bootstrap_slopes[
            bootstrap_index
        ] = (
            sen_slope(
                bootstrap_values,
                relative_time,
            )
            * 10.0
        )

    finite_slopes = bootstrap_slopes[
        np.isfinite(
            bootstrap_slopes
        )
    ]

    if finite_slopes.size == 0:

        return (
            np.nan,
            np.nan,
        )

    return (
        float(
            np.nanpercentile(
                finite_slopes,
                2.5,
            )
        ),
        float(
            np.nanpercentile(
                finite_slopes,
                97.5,
            )
        ),
    )


annual_field_dictionary = OrderedDict([
    (
        "SLA",
        sla_annual,
    ),
    (
        "Total steric",
        total_annual,
    ),
    (
        "Thermosteric",
        thermo_annual,
    ),
    (
        "Halosteric",
        halo_annual,
    ),
])

sector_rows = []

print(
    "\nCalculating 13-sea trends and 3-year block-bootstrap CIs..."
)

for variable_index, (
    variable_label,
    annual_values,
) in enumerate(
    annual_field_dictionary.items()
):

    variable_rows = []

    for sea_index, sea_code in enumerate(
        SEA_CODES
    ):

        sea_information = SEA_SECTORS[
            sea_code
        ]

        mask = sector_mask(
            sea_information
        )

        annual_series = (
            area_weighted_sector_series(
                annual_values,
                mask,
            )
        )

        valid_year_mask = np.isfinite(
            annual_series
        )

        selected_values = annual_series[
            valid_year_mask
        ]

        selected_years = analysis_years[
            valid_year_mask
        ].astype(
            np.float64
        )

        if selected_values.size < MIN_VALID_YEARS:

            slope_decade = np.nan
            p_value = np.nan
            correction_factor = np.nan
            ci_lower = np.nan
            ci_upper = np.nan

        else:

            slope_decade = (
                sen_slope(
                    selected_values,
                    selected_years,
                )
                * 10.0
            )

            _, p_value, correction_factor = (
                hamed_rao_modified_mk(
                    selected_values,
                    selected_years,
                )
            )

            ci_lower, ci_upper = (
                moving_block_bootstrap_ci(
                    selected_values,
                    selected_years,
                    block_length=(
                        BOOTSTRAP_BLOCK_LENGTH_YEARS
                    ),
                    number_of_bootstraps=(
                        BOOTSTRAP_REPLICATES
                    ),
                    random_seed=(
                        BOOTSTRAP_SEED
                        + variable_index
                        * 100
                        + sea_index
                    ),
                )
            )

        variable_rows.append(
            {
                "variable": variable_label,
                "sea_code": sea_code,
                "sea_name": sea_information[
                    "name"
                ],
                "sen_slope_cm_decade": slope_decade,
                "ci_lower_cm_decade": ci_lower,
                "ci_upper_cm_decade": ci_upper,
                "hamed_rao_p": p_value,
                "hamed_rao_factor": correction_factor,
                "valid_years": int(
                    selected_values.size
                ),
                "common_grid_cells": int(
                    mask.sum()
                ),
            }
        )

    variable_p_values = np.asarray(
        [
            row[
                "hamed_rao_p"
            ]
            for row in variable_rows
        ],
        dtype=np.float64,
    )

    variable_q_values = (
        benjamini_hochberg_qvalues(
            variable_p_values
        )
    )

    for row, q_value in zip(
        variable_rows,
        variable_q_values,
    ):

        row[
            "fdr_q"
        ] = q_value

        row[
            "fdr_significant"
        ] = bool(
            np.isfinite(
                q_value
            )
            and q_value
            < FDR_ALPHA
        )

        sector_rows.append(
            row
        )

sector_table = pd.DataFrame(
    sector_rows
)

sector_table.to_csv(
    OUT_SECTOR_CSV,
    index=False,
)

print(
    "Saved sector statistics:",
    OUT_SECTOR_CSV,
)


# 14. SAVE GRID AND ANNUAL OUTPUTS
grid_output = xr.Dataset(
    coords={
        "latitude": (
            "latitude",
            target_latitude.astype(
                np.float32
            ),
        ),
        "longitude": (
            "longitude",
            target_longitude.astype(
                np.float32
            ),
        ),
    }
)

grid_statistics_dictionary = OrderedDict([
    (
        "sla",
        sla_grid_statistics,
    ),
    (
        "total_steric",
        total_grid_statistics,
    ),
    (
        "thermosteric",
        thermo_grid_statistics,
    ),
    (
        "halosteric",
        halo_grid_statistics,
    ),
])

for short_name, statistics in (
    grid_statistics_dictionary.items()
):

    grid_output[
        f"{short_name}_sen_slope"
    ] = (
        (
            "latitude",
            "longitude",
        ),
        statistics[
            "slope"
        ],
    )

    grid_output[
        f"{short_name}_hamed_rao_p"
    ] = (
        (
            "latitude",
            "longitude",
        ),
        statistics[
            "p"
        ],
    )

    grid_output[
        f"{short_name}_fdr_q"
    ] = (
        (
            "latitude",
            "longitude",
        ),
        statistics[
            "q"
        ],
    )

    grid_output[
        f"{short_name}_fdr_significant"
    ] = (
        (
            "latitude",
            "longitude",
        ),
        statistics[
            "significant"
        ].astype(
            np.int8
        ),
    )

    grid_output[
        f"{short_name}_hamed_rao_factor"
    ] = (
        (
            "latitude",
            "longitude",
        ),
        statistics[
            "hamed_rao_factor"
        ],
    )

    grid_output[
        f"{short_name}_sen_slope"
    ].attrs[
        "units"
    ] = "cm decade-1"

grid_output[
    "fixed_common_mask"
] = (
    (
        "latitude",
        "longitude",
    ),
    common_mask.astype(
        np.int8
    ),
)

grid_output.attrs.update(
    {
        "title": (
            "Figure 11 study-period grid-cell trends"
        ),
        "study_period": (
            "2008-2025"
        ),
        "trend_phrase": (
            "Study-period trend during 2008–2025"
        ),
        "deseasonalization": (
            "Calendar-month climatology removed before annual averaging"
        ),
        "trend_estimator": (
            "Sen slope"
        ),
        "significance_test": (
            "Hamed-Rao modified Mann-Kendall"
        ),
        "multiple_testing": (
            "Benjamini-Hochberg FDR"
        ),
        "FDR_alpha": FDR_ALPHA,
    }
)

grid_encoding = {
    variable_name: {
        "zlib": True,
        "complevel": 4,
    }
    for variable_name
    in grid_output.data_vars
}

grid_output.to_netcdf(
    OUT_GRID_NC,
    encoding=grid_encoding,
)

grid_output.close()

del grid_output

annual_output = xr.Dataset(
    {
        "sla_annual_anomaly": (
            (
                "year",
                "latitude",
                "longitude",
            ),
            sla_annual,
        ),
        "total_steric_annual_anomaly": (
            (
                "year",
                "latitude",
                "longitude",
            ),
            total_annual,
        ),
        "thermosteric_annual_anomaly": (
            (
                "year",
                "latitude",
                "longitude",
            ),
            thermo_annual,
        ),
        "halosteric_annual_anomaly": (
            (
                "year",
                "latitude",
                "longitude",
            ),
            halo_annual,
        ),
        "fixed_common_mask": (
            (
                "latitude",
                "longitude",
            ),
            common_mask.astype(
                np.int8
            ),
        ),
    },
    coords={
        "year": analysis_years,
        "latitude": target_latitude.astype(
            np.float32
        ),
        "longitude": target_longitude.astype(
            np.float32
        ),
    },
)

for variable_name in [
    "sla_annual_anomaly",
    "total_steric_annual_anomaly",
    "thermosteric_annual_anomaly",
    "halosteric_annual_anomaly",
]:

    annual_output[
        variable_name
    ].attrs[
        "units"
    ] = "cm"

annual_output.attrs.update(
    {
        "title": (
            "Deseasonalized annual anomalies used in Figure 11"
        ),
        "method": (
            "Grid-cell monthly climatology removed, then annual mean calculated"
        ),
    }
)

annual_encoding = {
    variable_name: {
        "dtype": "float32",
        "zlib": True,
        "complevel": 4,
    }
    for variable_name
    in annual_output.data_vars
    if variable_name
    != "fixed_common_mask"
}

annual_encoding[
    "fixed_common_mask"
] = {
    "dtype": "int8",
    "zlib": True,
    "complevel": 4,
}

annual_output.to_netcdf(
    OUT_ANNUAL_NC,
    encoding=annual_encoding,
)

annual_output.close()

del annual_output

gc.collect()



# PART B — LATEST / FINAL CORRECTED FIGURE 11 PLOTTING
# This later source revision reads the products created in Part A.
# Trend and significance calculations are not repeated here.

# FIGURE 11 — STUDY-PERIOD TRENDS
# FINAL PLOT-ONLY WORKFLOW
#
# IMPORTANT
# This code DOES NOT recompute:
#   • deseasonalized annual anomalies;
#   • Sen slopes;
#   • Hamed–Rao modified Mann–Kendall tests;
#   • FDR correction;
#   • block-bootstrap confidence intervals.
#
# It reads the already calculated outputs:
#
# 1. Fig11_gridcell_trend_statistics_2008_2025.nc
# 2. Fig11_13sea_Sen_slopes_bootstrap_CI.csv
#
# The annual-anomaly NetCDF is checked for completeness but is not
# loaded because it is not required to redraw Figure 11.


# 3. INPUT AND OUTPUT PATHS
GRID_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/paper2/"
    "Fig11_gridcell_trend_statistics_2008_2025.nc"
)

ANNUAL_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/paper2/"
    "Fig11_deseasonalized_annual_anomalies_common_grid.nc"
)

SECTOR_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/paper2/"
    "Fig11_13sea_Sen_slopes_bootstrap_CI.csv"
)

OUTPUT_DIRECTORY = Path(
    "/content/drive/MyDrive/SAM_Thesis/paper2"
)

OUTPUT_DIRECTORY.mkdir(
    parents=True,
    exist_ok=True,
)

OUT_PNG = (
    OUTPUT_DIRECTORY
    / "Figure11_study_period_trends_FINAL_CORRECTED_1080dpi.png"
)

OUT_PDF = (
    OUTPUT_DIRECTORY
    / "Figure11_study_period_trends_FINAL_CORRECTED.pdf"
)

OUT_SUMMARY = (
    OUTPUT_DIRECTORY
    / "Figure11_FINAL_CORRECTED_plot_summary.txt"
)


# 4. DISPLAY SETTINGS
LAT_MIN = -90.0
LAT_MAX = -60.0

FDR_ALPHA = 0.05
SAVE_DPI = 1080

DRAW_SECTOR_BOUNDARIES = True

# FDR-significance stippling size requested for the final map.
SIGNIFICANCE_DOT_SIZE = 3.0
SIGNIFICANCE_DOT_ALPHA = 0.95

# A compact physical canvas prevents excessive RAM use at 1080 dpi.
FIGSIZE = (
    10.4,
    7.25,
)

FONT_MAP_TITLE = 10.4
FONT_FOREST_TITLE = 10.0
FONT_GEO_LABEL = 7.0
FONT_FOREST_TICK = 7.3
FONT_AXIS_LABEL = 8.3
FONT_COLORBAR = 7.8
FONT_LEGEND = 7.3
FONT_FOOTNOTE = 7.1

LAND_COLOR = "0.84"
COAST_COLOR = "0.25"
NO_DATA_COLOR = "0.86"
SECTOR_LINE_COLOR = "white"

# Radial meridian gridlines every 30 degrees.
MERIDIANS_FOR_GRID = np.arange(
    -180,
    181,
    30,
)

# Longitude labels around the circular border every 30 degrees,
# excluding 0 degrees as requested.
LONGITUDE_LABELS = [
    -150,
    -120,
    -90,
    -60,
    -30,
    30,
    60,
    90,
    120,
    150,
    180,
]

PARALLELS_FOR_GRID = [
    -60,
    -70,
    -80,
]

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "mathtext.default": "regular",
    }
)


# 5. 13 ANTARCTIC SEA-SECTOR DEFINITIONS
SEA_SECTORS = OrderedDict([
    (
        "WED",
        {
            "name": "Weddell Sea",
            "lon_min": -60.0,
            "lon_max": -20.0,
        },
    ),
    (
        "KHV",
        {
            "name": "King Haakon VII Sea",
            "lon_min": -20.0,
            "lon_max": 10.0,
        },
    ),
    (
        "RLS",
        {
            "name": "Riiser-Larsen Sea",
            "lon_min": 10.0,
            "lon_max": 35.0,
        },
    ),
    (
        "LAZ",
        {
            "name": "Lazarev Sea",
            "lon_min": 35.0,
            "lon_max": 60.0,
        },
    ),
    (
        "COS",
        {
            "name": "Cosmonauts Sea",
            "lon_min": 60.0,
            "lon_max": 90.0,
        },
    ),
    (
        "COO",
        {
            "name": "Cooperation Sea",
            "lon_min": 90.0,
            "lon_max": 115.0,
        },
    ),
    (
        "DAV",
        {
            "name": "Davis Sea",
            "lon_min": 115.0,
            "lon_max": 130.0,
        },
    ),
    (
        "MAW",
        {
            "name": "Mawson Sea",
            "lon_min": 130.0,
            "lon_max": 150.0,
        },
    ),
    (
        "DUR",
        {
            "name": "D'Urville Sea",
            "lon_min": 150.0,
            "lon_max": 170.0,
        },
    ),
    (
        "SOM",
        {
            "name": "Somov Sea",
            "lon_min": 170.0,
            "lon_max": -160.0,
        },
    ),
    (
        "ROS",
        {
            "name": "Ross Sea",
            "lon_min": -160.0,
            "lon_max": -130.0,
        },
    ),
    (
        "AMU",
        {
            "name": "Amundsen Sea",
            "lon_min": -130.0,
            "lon_max": -100.0,
        },
    ),
    (
        "BEL",
        {
            "name": "Bellingshausen Sea",
            "lon_min": -100.0,
            "lon_max": -60.0,
        },
    ),
])

SEA_CODES = list(
    SEA_SECTORS.keys()
)


# 6. INPUT VERIFICATION
for required_file in [
    GRID_FILE,
    SECTOR_FILE,
]:

    if not os.path.exists(
        required_file
    ):

        raise FileNotFoundError(
            "Required precomputed output was not found:\n"
            f"{required_file}"
        )

if os.path.exists(
    ANNUAL_FILE
):

    print(
        "Annual anomaly file found and retained without loading:"
    )

    print(
        ANNUAL_FILE
    )

else:

    print(
        "WARNING: annual anomaly file was not found, "
        "but it is not required for this plot-only workflow."
    )


# 7. OPEN PRECOMPUTED GRID STATISTICS
def open_dataset_safely(
    file_path,
):
    """
    Open the small precomputed trend-statistics NetCDF.
    """
    attempts = []

    for engine in [
        None,
        "netcdf4",
        "h5netcdf",
        "scipy",
    ]:

        try:

            kwargs = {
                "decode_times": True,
                "mask_and_scale": True,
                "cache": False,
            }

            if engine is not None:

                kwargs[
                    "engine"
                ] = engine

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
                f"with engine={engine_name}"
            )

            return dataset

        except Exception as error:

            attempts.append(
                f"engine={engine}: {error}"
            )

    raise RuntimeError(
        f"Could not open:\n"
        f"{file_path}\n\n"
        + "\n".join(
            attempts
        )
    )


def detect_coordinate(
    dataset,
    coordinate_type,
):
    """
    Detect latitude and longitude coordinate names.
    """
    aliases = {
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

    candidates = list(
        dataset.coords
    ) + list(
        dataset.variables
    )

    for name in candidates:

        variable = dataset[
            name
        ]

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

        if coordinate_type == "lat":

            if (
                lower_name
                in aliases["lat"]
                or standard_name
                == "latitude"
                or axis == "Y"
                or "degrees_north"
                in units
            ):

                return name

        elif coordinate_type == "lon":

            if (
                lower_name
                in aliases["lon"]
                or standard_name
                == "longitude"
                or axis == "X"
                or "degrees_east"
                in units
            ):

                return name

    raise KeyError(
        f"Could not detect the "
        f"{coordinate_type} coordinate."
    )


def standardize_grid_dataset(
    dataset,
):
    """
    Rename coordinates and normalize longitude to -180...180.
    """
    latitude_name = detect_coordinate(
        dataset,
        "lat",
    )

    longitude_name = detect_coordinate(
        dataset,
        "lon",
    )

    rename_mapping = {}

    if latitude_name != "lat":

        rename_mapping[
            latitude_name
        ] = "lat"

    if longitude_name != "lon":

        rename_mapping[
            longitude_name
        ] = "lon"

    if rename_mapping:

        dataset = dataset.rename(
            rename_mapping
        )

    normalized_longitude = (
        (
            dataset[
                "lon"
            ].astype(float)
            + 180.0
        )
         % 360.0
    ) - 180.0

    dataset = dataset.assign_coords(
        lon=normalized_longitude
    )

    longitude_values = np.asarray(
        dataset[
            "lon"
        ].values
    )

    _, unique_indices = np.unique(
        longitude_values,
        return_index=True,
    )

    dataset = dataset.isel(
        lon=np.sort(
            unique_indices
        )
    )

    dataset = dataset.sortby(
        "lon"
    )

    dataset = dataset.sortby(
        "lat"
    )

    return dataset


grid_dataset = standardize_grid_dataset(
    open_dataset_safely(
        GRID_FILE
    )
)

required_grid_variables = [
    "sla_sen_slope",
    "sla_fdr_significant",
    "total_steric_sen_slope",
    "total_steric_fdr_significant",
    "thermosteric_sen_slope",
    "thermosteric_fdr_significant",
    "halosteric_sen_slope",
    "halosteric_fdr_significant",
]

missing_grid_variables = [
    variable_name
    for variable_name
    in required_grid_variables
    if variable_name
    not in grid_dataset.data_vars
]

if missing_grid_variables:

    raise KeyError(
        "The grid NetCDF is missing:\n"
        + "\n".join(
            missing_grid_variables
        )
        + "\n\nAvailable variables:\n"
        + "\n".join(
            grid_dataset.data_vars
        )
    )

target_latitude = np.asarray(
    grid_dataset[
        "lat"
    ].values,
    dtype=np.float64,
)

target_longitude = np.asarray(
    grid_dataset[
        "lon"
    ].values,
    dtype=np.float64,
)

grid_statistics = OrderedDict([
    (
        "SLA",
        {
            "slope": np.asarray(
                grid_dataset[
                    "sla_sen_slope"
                ].values,
                dtype=np.float32,
            ),
            "significant": (
                np.asarray(
                    grid_dataset[
                        "sla_fdr_significant"
                    ].values
                )
                > 0
            ),
        },
    ),
    (
        "Total steric",
        {
            "slope": np.asarray(
                grid_dataset[
                    "total_steric_sen_slope"
                ].values,
                dtype=np.float32,
            ),
            "significant": (
                np.asarray(
                    grid_dataset[
                        "total_steric_fdr_significant"
                    ].values
                )
                > 0
            ),
        },
    ),
    (
        "Thermosteric",
        {
            "slope": np.asarray(
                grid_dataset[
                    "thermosteric_sen_slope"
                ].values,
                dtype=np.float32,
            ),
            "significant": (
                np.asarray(
                    grid_dataset[
                        "thermosteric_fdr_significant"
                    ].values
                )
                > 0
            ),
        },
    ),
    (
        "Halosteric",
        {
            "slope": np.asarray(
                grid_dataset[
                    "halosteric_sen_slope"
                ].values,
                dtype=np.float32,
            ),
            "significant": (
                np.asarray(
                    grid_dataset[
                        "halosteric_fdr_significant"
                    ].values
                )
                > 0
            ),
        },
    ),
])

grid_dataset.close()

del grid_dataset

gc.collect()


# 8. OPEN PRECOMPUTED SECTOR STATISTICS
sector_table = pd.read_csv(
    SECTOR_FILE
)

required_sector_columns = [
    "variable",
    "sea_code",
    "sen_slope_cm_decade",
    "ci_lower_cm_decade",
    "ci_upper_cm_decade",
    "fdr_significant",
]

missing_sector_columns = [
    column_name
    for column_name
    in required_sector_columns
    if column_name
    not in sector_table.columns
]

if missing_sector_columns:

    raise KeyError(
        "The sector CSV is missing:\n"
        + "\n".join(
            missing_sector_columns
        )
        + "\n\nAvailable columns:\n"
        + "\n".join(
            sector_table.columns
        )
    )


def to_boolean(
    value,
):
    """
    Robust conversion of CSV boolean fields.
    """
    if isinstance(
        value,
        bool,
    ):

        return value

    if pd.isna(
        value
    ):

        return False

    text = str(
        value
    ).strip().lower()

    return text in {
        "true",
        "1",
        "yes",
        "y",
        "t",
    }


sector_table[
    "fdr_significant"
] = sector_table[
    "fdr_significant"
].map(
    to_boolean
)


# 9. FIGURE LIMITS
def collect_finite_values(
    arrays,
):
    """
    Combine finite values from small precomputed arrays.
    """
    collected = []

    for array in arrays:

        values = np.asarray(
            array,
            dtype=np.float64,
        )

        finite = values[
            np.isfinite(
                values
            )
        ]

        if finite.size > 0:

            collected.append(
                finite
            )

    if not collected:

        return np.array(
            [
                -1.0,
                1.0,
            ],
            dtype=np.float64,
        )

    return np.concatenate(
        collected
    )


def nice_symmetric_limit(
    value,
):
    """
    Round upward to a readable symmetric endpoint.
    """
    value = abs(
        float(
            value
        )
    )

    if (
        not np.isfinite(
            value
        )
        or value <= 0.0
    ):

        return 1.0

    exponent = np.floor(
        np.log10(
            value
        )
    )

    scale = 10.0 ** exponent

    normalized = value / scale

    for candidate in [
        1.0,
        1.5,
        2.0,
        2.5,
        3.0,
        4.0,
        5.0,
        6.0,
        8.0,
        10.0,
    ]:

        if normalized <= candidate:

            return candidate * scale

    return 10.0 * scale


map_values = collect_finite_values(
    [
        statistics[
            "slope"
        ]
        for statistics
        in grid_statistics.values()
    ]
)

MAP_COLOR_LIMIT = nice_symmetric_limit(
    np.nanpercentile(
        np.abs(
            map_values
        ),
        99.0,
    )
)

sector_values = collect_finite_values(
    [
        sector_table[
            "sen_slope_cm_decade"
        ].values,
        sector_table[
            "ci_lower_cm_decade"
        ].values,
        sector_table[
            "ci_upper_cm_decade"
        ].values,
    ]
)

SECTOR_X_LIMIT = nice_symmetric_limit(
    np.nanpercentile(
        np.abs(
            sector_values
        ),
        99.5,
    )
)

print(
    "\nFigure limits:"
)

print(
    f"Map color scale: ±"
    f"{MAP_COLOR_LIMIT:.4g} cm decade-1"
)

print(
    f"Sector x-axis : ±"
    f"{SECTOR_X_LIMIT:.4g} cm decade-1"
)


# ================================================================
# 10. MAP HELPERS
# ================================================================
longitude_2d, latitude_2d = np.meshgrid(
    target_longitude,
    target_latitude,
)


def circular_boundary():
    """
    Circular clipping path for South Polar maps.
    """
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

    return mpath.Path(
        vertices * 0.5
        + np.array(
            [
                0.5,
                0.5,
            ]
        )
    )


def format_longitude_label(
    longitude,
):
    """
    Format a longitude label without displaying 0 degrees.
    """
    longitude = int(
        longitude
    )

    if abs(
        longitude
    ) == 180:

        return "180°"

    if longitude < 0:

        return f"{abs(longitude)}°W"

    return f"{longitude}°E"


def add_geographic_labels(
    axis,
):
    """
    Add longitude labels every 30 degrees except 0 degrees.

    The labels are placed slightly outside the circular map boundary.
    Their reduced font size prevents overlap among the four maps.
    """
    for longitude in LONGITUDE_LABELS:

        axis.text(
            float(
                longitude
            ),
            -58.35,
            format_longitude_label(
                longitude
            ),
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
        fontsize=FONT_GEO_LABEL + 0.5,
        fontweight="bold",
        color="black",
        zorder=30,
    )


def style_map_axis(
    axis,
):
    """
    Common South Polar map styling.
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

    axis.set_boundary(
        circular_boundary(),
        transform=axis.transAxes,
    )

    axis.add_feature(
        cfeature.LAND.with_scale(
            "110m"
        ),
        facecolor=LAND_COLOR,
        edgecolor=COAST_COLOR,
        linewidth=0.50,
        zorder=8,
    )

    axis.coastlines(
        resolution="110m",
        color=COAST_COLOR,
        linewidth=0.50,
        zorder=9,
    )

    axis.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        xlocs=MERIDIANS_FOR_GRID,
        ylocs=PARALLELS_FOR_GRID,
        linewidth=0.50,
        linestyle=":",
        color="0.38",
        alpha=0.90,
        zorder=3,
    )

    add_geographic_labels(
        axis
    )

    if DRAW_SECTOR_BOUNDARIES:

        boundary_longitudes = sorted(
            {
                float(
                    sea_information[
                        "lon_min"
                    ]
                )
                for sea_information
                in SEA_SECTORS.values()
            }
        )

        southern_latitude = max(
            -85.0,
            float(
                np.nanmin(
                    target_latitude
                )
            ),
        )

        for boundary_longitude in boundary_longitudes:

            axis.plot(
                [
                    boundary_longitude,
                    boundary_longitude,
                ],
                [
                    LAT_MAX,
                    southern_latitude,
                ],
                transform=ccrs.PlateCarree(),
                color=SECTOR_LINE_COLOR,
                linewidth=0.70,
                alpha=0.98,
                zorder=12,
            )


def add_cyclic_map_field(
    values,
):
    """
    Remove the longitude seam.
    """
    cyclic_values, cyclic_longitude = add_cyclic_point(
        np.asarray(
            values,
            dtype=np.float32,
        ),
        coord=target_longitude,
        axis=-1,
    )

    return (
        cyclic_values,
        cyclic_longitude,
    )


def plot_trend_map(
    axis,
    slope_map,
    significant_map,
    title,
):
    """
    Plot one Sen-slope map and larger FDR-significance dots.
    """
    style_map_axis(
        axis
    )

    cyclic_slope, cyclic_longitude = add_cyclic_map_field(
        slope_map
    )

    no_data = np.where(
        np.isfinite(
            cyclic_slope
        ),
        np.nan,
        1.0,
    )

    axis.pcolormesh(
        cyclic_longitude,
        target_latitude,
        no_data,
        transform=ccrs.PlateCarree(),
        cmap=ListedColormap(
            [
                NO_DATA_COLOR
            ]
        ),
        vmin=0.0,
        vmax=1.0,
        shading="auto",
        zorder=1,
    )

    image = axis.pcolormesh(
        cyclic_longitude,
        target_latitude,
        cyclic_slope,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        vmin=-MAP_COLOR_LIMIT,
        vmax=MAP_COLOR_LIMIT,
        shading="auto",
        rasterized=True,
        zorder=4,
    )

    stipple_mask = (
        np.asarray(
            significant_map,
            dtype=bool,
        )
        & np.isfinite(
            slope_map
        )
    )

    if np.any(
        stipple_mask
    ):

        axis.scatter(
            longitude_2d[
                stipple_mask
            ],
            latitude_2d[
                stipple_mask
            ],
            transform=ccrs.PlateCarree(),
            s=SIGNIFICANCE_DOT_SIZE,
            marker="o",
            facecolor="black",
            edgecolor="none",
            linewidths=0.0,
            alpha=SIGNIFICANCE_DOT_ALPHA,
            zorder=20,
        )

    axis.set_title(
        title,
        fontsize=FONT_MAP_TITLE,
        fontweight="bold",
        pad=8,
    )

    return image


# 11. STRONG FOREST-PLOT COLOR SYSTEM
# The original RdBu center is nearly white and made small trends
# difficult to see. This custom map uses a dark-grey midpoint and
# saturated blue/red sides so every estimate remains visible.
STRONG_TREND_CMAP = (
    mcolors.LinearSegmentedColormap.from_list(
        "StrongTrend",
        [
            (
                0.00,
                "#08306B",
            ),
            (
                0.22,
                "#2171B5",
            ),
            (
                0.42,
                "#6BAED6",
            ),
            (
                0.50,
                "#4A4A4A",
            ),
            (
                0.58,
                "#F08080",
            ),
            (
                0.78,
                "#D7301F",
            ),
            (
                1.00,
                "#67000D",
            ),
        ],
    )
)

SECTOR_NORMALIZATION = mcolors.TwoSlopeNorm(
    vmin=-SECTOR_X_LIMIT,
    vcenter=0.0,
    vmax=SECTOR_X_LIMIT,
)


def style_forest_axis(
    axis,
    title,
):
    """
    Compact forest-plot styling with a two-line title.
    """
    y_positions = np.arange(
        len(
            SEA_CODES
        )
    )

    axis.axvline(
        0.0,
        color="0.35",
        linestyle="--",
        linewidth=0.85,
        zorder=1,
    )

    axis.set_yticks(
        y_positions
    )

    axis.set_yticklabels(
        SEA_CODES,
        fontsize=FONT_FOREST_TICK,
        fontweight="bold",
    )

    axis.set_ylim(
        len(
            SEA_CODES
        )
        - 0.5,
        -0.5,
    )

    axis.set_xlim(
        -SECTOR_X_LIMIT,
        SECTOR_X_LIMIT,
    )

    axis.set_xlabel(
        "Trend (cm decade$^{-1}$)",
        fontsize=FONT_AXIS_LABEL,
        fontweight="bold",
        labelpad=3,
    )

    axis.set_title(
        title,
        fontsize=FONT_FOREST_TITLE,
        fontweight="bold",
        linespacing=1.05,
        pad=7,
    )

    axis.tick_params(
        axis="x",
        labelsize=FONT_FOREST_TICK,
        width=0.8,
        length=3,
        pad=2,
    )

    axis.tick_params(
        axis="y",
        width=0.0,
        length=0,
        pad=2,
    )

    axis.spines[
        "top"
    ].set_visible(
        False
    )

    axis.spines[
        "right"
    ].set_visible(
        False
    )

    axis.spines[
        "left"
    ].set_linewidth(
        0.8
    )

    axis.spines[
        "bottom"
    ].set_linewidth(
        0.8
    )


def plot_sector_forest(
    axis,
    variable_label,
    title,
):
    """
    Plot strong-color Sen slopes and bootstrap confidence intervals.
    """
    subset = (
        sector_table[
            sector_table[
                "variable"
            ]
            == variable_label
        ]
        .set_index(
            "sea_code"
        )
        .reindex(
            SEA_CODES
        )
        .reset_index()
    )

    style_forest_axis(
        axis,
        title,
    )

    for sea_index, row in subset.iterrows():

        slope = row[
            "sen_slope_cm_decade"
        ]

        lower = row[
            "ci_lower_cm_decade"
        ]

        upper = row[
            "ci_upper_cm_decade"
        ]

        significant = bool(
            row[
                "fdr_significant"
            ]
        )

        if not np.isfinite(
            slope
        ):

            continue

        trend_color = STRONG_TREND_CMAP(
            SECTOR_NORMALIZATION(
                slope
            )
        )

        if (
            np.isfinite(
                lower
            )
            and np.isfinite(
                upper
            )
        ):

            x_error = np.array(
                [
                    [
                        max(
                            0.0,
                            slope - lower,
                        )
                    ],
                    [
                        max(
                            0.0,
                            upper - slope,
                        )
                    ],
                ]
            )

        else:

            x_error = None

        marker_facecolor = (
            trend_color
            if significant
            else "white"
        )

        axis.errorbar(
            slope,
            sea_index,
            xerr=x_error,
            fmt="o",
            markersize=5.4,
            markerfacecolor=marker_facecolor,
            markeredgecolor=trend_color,
            markeredgewidth=1.25,
            ecolor=trend_color,
            elinewidth=1.15,
            capsize=2.8,
            capthick=1.05,
            zorder=5,
        )


# 12. BUILD THE FINAL NON-OVERLAPPING FIGURE
projection = ccrs.SouthPolarStereo(
    central_longitude=0.0
)

fig = plt.figure(
    figsize=FIGSIZE
)


# Top row: maps
MAP_HEIGHT = 0.266

MAP_WIDTH = (
    MAP_HEIGHT
    * FIGSIZE[1]
    / FIGSIZE[0]
)

COLUMN_CENTERS = [
    0.135,
    0.380,
    0.625,
    0.870,
]

MAP_BOTTOM = 0.688

map_lefts = [
    center
    - MAP_WIDTH / 2.0
    for center in COLUMN_CENTERS
]

map_titles = [
    "(a) SLA trend map",
    "(b) Total steric trend map",
    "(c) Thermosteric trend map",
    "(d) Halosteric trend map",
]

map_keys = [
    "SLA",
    "Total steric",
    "Thermosteric",
    "Halosteric",
]

map_images = []

for column_index, map_key in enumerate(
    map_keys
):

    map_axis = fig.add_axes(
        [
            map_lefts[
                column_index
            ],
            MAP_BOTTOM,
            MAP_WIDTH,
            MAP_HEIGHT,
        ],
        projection=projection,
    )

    map_image = plot_trend_map(
        axis=map_axis,
        slope_map=grid_statistics[
            map_key
        ][
            "slope"
        ],
        significant_map=grid_statistics[
            map_key
        ][
            "significant"
        ],
        title=map_titles[
            column_index
        ],
    )

    map_images.append(
        map_image
    )


# Dedicated map-colorbar row
MAP_COLORBAR_Y = 0.634
MAP_COLORBAR_HEIGHT = 0.0105
MAP_COLORBAR_WIDTH = MAP_WIDTH * 0.91

for column_index, map_image in enumerate(
    map_images
):

    colorbar_axis = fig.add_axes(
        [
            COLUMN_CENTERS[
                column_index
            ]
            - MAP_COLORBAR_WIDTH / 2.0,
            MAP_COLORBAR_Y,
            MAP_COLORBAR_WIDTH,
            MAP_COLORBAR_HEIGHT,
        ]
    )

    colorbar = fig.colorbar(
        map_image,
        cax=colorbar_axis,
        orientation="horizontal",
        extend="both",
    )

    colorbar.set_ticks(
        np.linspace(
            -MAP_COLOR_LIMIT,
            MAP_COLOR_LIMIT,
            5,
        )
    )

    colorbar.set_label(
        "cm decade$^{-1}$",
        fontsize=FONT_COLORBAR,
        fontweight="bold",
        labelpad=1.5,
    )

    colorbar.ax.tick_params(
        labelsize=FONT_FOREST_TICK,
        width=0.7,
        length=2.5,
        pad=1,
    )

    for tick_label in colorbar.ax.get_xticklabels():

        tick_label.set_fontweight(
            "bold"
        )


# Separate clean information row
map_significance_handle = Line2D(
    [0],
    [0],
    marker="o",
    linestyle="None",
    markerfacecolor="black",
    markeredgecolor="black",
    markersize=4.2,
    label="FDR-significant grid cell (q < 0.05)",
)

map_legend = fig.legend(
    handles=[
        map_significance_handle
    ],
    loc="center",
    bbox_to_anchor=(
        0.285,
        0.570,
    ),
    frameon=True,
    fontsize=FONT_LEGEND,
    borderpad=0.50,
    handlelength=1.0,
    handletextpad=0.55,
)

map_legend.get_frame().set_edgecolor(
    "0.25"
)

map_legend.get_frame().set_linewidth(
    0.8
)

for text in map_legend.get_texts():

    text.set_fontweight(
        "bold"
    )


sector_method_handle = Line2D(
    [0],
    [0],
    color="black",
    linewidth=1.2,
    marker="|",
    markersize=8,
    markeredgewidth=1.0,
    label=(
        "Sen slope with 3-year block-bootstrap 95% CI"
    ),
)

sector_method_legend = fig.legend(
    handles=[
        sector_method_handle
    ],
    loc="center",
    bbox_to_anchor=(
        0.720,
        0.570,
    ),
    frameon=True,
    fontsize=FONT_LEGEND,
    borderpad=0.50,
    handlelength=2.0,
    handletextpad=0.65,
)

sector_method_legend.get_frame().set_edgecolor(
    "0.25"
)

sector_method_legend.get_frame().set_linewidth(
    0.8
)

for text in sector_method_legend.get_texts():

    text.set_fontweight(
        "bold"
    )


# Bottom row: forest plots
FOREST_BOTTOM = 0.205
FOREST_HEIGHT = 0.290
FOREST_WIDTH = 0.205

forest_lefts = [
    center
    - FOREST_WIDTH / 2.0
    for center in COLUMN_CENTERS
]

forest_titles = [
    "(e) 13-sea SLA\nSen slopes",
    "(f) 13-sea total steric\nSen slopes",
    "(g) 13-sea thermosteric\nSen slopes",
    "(h) 13-sea halosteric\nSen slopes",
]

forest_variable_labels = [
    "SLA",
    "Total steric",
    "Thermosteric",
    "Halosteric",
]

for column_index in range(
    4
):

    forest_axis = fig.add_axes(
        [
            forest_lefts[
                column_index
            ],
            FOREST_BOTTOM,
            FOREST_WIDTH,
            FOREST_HEIGHT,
        ]
    )

    plot_sector_forest(
        axis=forest_axis,
        variable_label=forest_variable_labels[
            column_index
        ],
        title=forest_titles[
            column_index
        ],
    )


# 13. BOTTOM LEGENDS AND SHARED SECTOR COLORBAR
significance_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor="black",
        markeredgecolor="black",
        markersize=5.5,
        label="FDR significant (q < 0.05)",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor="white",
        markeredgecolor="black",
        markersize=5.5,
        label="Not significant (q ≥ 0.05)",
    ),
]

sector_significance_legend = fig.legend(
    handles=significance_handles,
    loc="lower left",
    bbox_to_anchor=(
        0.045,
        0.055,
    ),
    frameon=True,
    fontsize=FONT_LEGEND,
    labelspacing=0.42,
    handlelength=1.1,
    handletextpad=0.50,
    borderpad=0.55,
)

sector_significance_legend.get_frame().set_edgecolor(
    "0.25"
)

sector_significance_legend.get_frame().set_linewidth(
    0.8
)

for text in sector_significance_legend.get_texts():

    text.set_fontweight(
        "bold"
    )


# Strong-color shared slope colorbar.
sector_colorbar_axis = fig.add_axes(
    [
        0.365,
        0.077,
        0.270,
        0.011,
    ]
)

sector_scalar_mappable = plt.cm.ScalarMappable(
    norm=SECTOR_NORMALIZATION,
    cmap=STRONG_TREND_CMAP,
)

sector_colorbar = fig.colorbar(
    sector_scalar_mappable,
    cax=sector_colorbar_axis,
    orientation="horizontal",
    extend="both",
)

sector_colorbar.set_ticks(
    np.linspace(
        -SECTOR_X_LIMIT,
        SECTOR_X_LIMIT,
        5,
    )
)

sector_colorbar.set_label(
    "Sen’s slope (cm decade$^{-1}$)",
    fontsize=FONT_COLORBAR,
    fontweight="bold",
    labelpad=2,
)

sector_colorbar.ax.tick_params(
    labelsize=FONT_FOREST_TICK,
    width=0.7,
    length=2.5,
    pad=1,
)

for tick_label in sector_colorbar.ax.get_xticklabels():

    tick_label.set_fontweight(
        "bold"
    )


# Statistical-method footnote, separate from all legends/colorbars.
fig.text(
    0.805,
    0.076,
    (
        "Hamed–Rao modified Mann–Kendall\n"
        "Benjamini–Hochberg FDR correction"
    ),
    ha="center",
    va="center",
    fontsize=FONT_FOOTNOTE,
    fontweight="bold",
    linespacing=1.20,
)


# 14. SAVE FIGURE
print(
    "\nSaving vector PDF..."
)

fig.savefig(
    OUT_PDF,
    facecolor="white",
    edgecolor="none",
    bbox_inches=None,
)

gc.collect()

print(
    "Saving compact 1080-dpi PNG..."
)

fig.savefig(
    OUT_PNG,
    dpi=SAVE_DPI,
    facecolor="white",
    edgecolor="none",
    bbox_inches=None,
)

plt.close(
    fig
)

gc.collect()


# 15. FINAL SUMMARY
summary_lines = [
    "=" * 94,
    "FIGURE 11 — FINAL CORRECTED PLOT-ONLY WORKFLOW",
    "=" * 94,
    "",
    "Precomputed inputs:",
    f"Grid statistics : {GRID_FILE}",
    f"Sector table    : {SECTOR_FILE}",
    f"Annual fields   : {ANNUAL_FILE}",
    "",
    "No trend or significance calculation was repeated.",
    "",
    "Corrections:",
    "• Overall figure title removed.",
    "• Map colorbars placed in an independent row.",
    "• Overlapping statistics text replaced by two separate legends.",
    f"• Map FDR-dot size set to {SIGNIFICANCE_DOT_SIZE}.",
    "• Radial meridian gridlines drawn every 30 degrees.",
    "• Longitude labels shown every 30 degrees except 0 degrees.",
    "• Forest-plot titles split across two lines.",
    "• Strong custom trend colors used for sector estimates and CIs.",
    "• Bottom significance legend, colorbar and method note separated.",
    "",
    "Output:",
    f"1080-dpi PNG : {OUT_PNG}",
    f"Vector PDF   : {OUT_PDF}",
    "",
    "Phrase retained in the processing record:",
    "Study-period trend during 2008–2025",
]

OUT_SUMMARY.write_text(
    "\n".join(
        summary_lines
    ),
    encoding="utf-8",
)

print(
    "\n".join(
        summary_lines
    )
)

print(
    "\nFIGURE 11 FINAL PLOT COMPLETED SUCCESSFULLY."
)
