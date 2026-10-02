"""
FIGURE 9 — FULL ORGANIZED LATEST WORKFLOW

FIGURE

EN4 Freshwater Storage and Upper-Ocean Stability
Southern Ocean, 2008-2025

SOURCE

Extracted and organized from:
    so_sealevel_paper_fig.py

VERSION SELECTION
The source contains two Figure 9 stages:

1. A complete LOW-RAM / COLAB-SAFE monthly-processing workflow that creates:
       Fig09_EN4_freshwater_stratification_seasonal_fields.nc

2. A later FINAL plot-only workflow that reads that seasonal NetCDF and creates:
       Figure09_EN4_freshwater_stratification_FINAL_1080dpi.png
       Figure09_EN4_freshwater_stratification_FINAL.pdf

This standalone file combines the required processing stage with only the
latest/final plotting stage. The older 600-dpi Figure 9 renderer is excluded.

ORGANIZED ORDER
1. Mount Google Drive.
2. Install required packages.
3. Import all libraries.
4. Define monthly source files and paths.
5. Define analysis settings.
6. Validate files and standardize datasets.
7. Run the low-RAM sequential monthly aggregation.
8. Create seasonal FWC, SSS anomaly, stratification, MLD and SIC fields.
9. Save the precomputed seasonal Figure 9 NetCDF.
10. Read that seasonal NetCDF with the latest final plotting workflow.
11. Calculate robust final display limits.
12. Draw the final 4 x 3 polar figure.
13. Overlay solid MLD contours and the dashed 15% SIC contour.
14. Save the final vector PDF.
15. Save the final 1080-dpi PNG.
16. Save the latest final summary text.

PRIMARY DATA OUTPUT
/content/drive/MyDrive/SAM_Thesis/paper2/
    Fig09_EN4_freshwater_stratification_seasonal_fields.nc

LATEST FINAL OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/
    Figure09_EN4_freshwater_stratification_FINAL_1080dpi.png
    Figure09_EN4_freshwater_stratification_FINAL.pdf
    Figure09_EN4_freshwater_stratification_FINAL_summary.txt

SCIENTIFIC PRESERVATION
The source-defined Southern Ocean domain, 2008-2025 study period, seasonal
definitions, FWC/SSS/stratification calculations, sequential low-RAM strategy,
MLD contours, SIC 15% contour, regridding choices and latest plotting logic
are retained. Only organization and removal of the obsolete earlier renderer
were performed.

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
    "cftime": "cftime",
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
    print("All required packages are already installed.")


# 3. IMPORT LIBRARIES

import gc
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.path as mpath

from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from cartopy.util import add_cyclic_point

import cartopy.crs as ccrs
import cartopy.feature as cfeature

warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
)


# PART A — LOW-RAM FIGURE 9 DATA GENERATION
# Creates the seasonal NetCDF required by the latest Figure 9 renderer.

# FIGURE 9 — EN4 FRESHWATER AND STRATIFICATION
# LOW-RAM / COLAB-SAFE WORKFLOW
# Southern Ocean, 2008–2025
#
# The script reads only ONE MONTH at a time.
# It never loads complete 216-month MLD or SIC arrays into memory.
#
# Columns
# 1. FWC seasonal climatology
# 2. SSS seasonal departure from the 2008–2025 all-month mean
# 3. Stratification seasonal climatology
#       + solid MLD contours
#       + dashed SIC 15% contour
#
# Output directory
# /content/drive/MyDrive/SAM_Thesis/paper2/


# 3. INPUT AND OUTPUT PATHS
FWC_CANDIDATES = [
    (
        "/content/drive/MyDrive/SAM_Thesis/Processed/"
        "EN4_FWC_0_200m_monthly_2008_2025_SO.nc"
    ),
    (
        "/content/drive/MyDrive/SAM_Thesis/Processed/"
        "EN4_NetCDF_inventory/"
        "EN4_freshwater_content_0_200m_monthly_2008_2025_SO.nc"
    ),
]

SSS_CANDIDATES = [
    (
        "/content/drive/MyDrive/SAM_Thesis/Processed/"
        "EN4_SSS_0_30m_monthly_2008_2025_SO.nc"
    ),
    (
        "/content/drive/MyDrive/SAM_Thesis/Processed/"
        "EN4_NetCDF_inventory/"
        "EN4_surface_salinity_0_30m_monthly_2008_2025_SO.nc"
    ),
]

STRATIFICATION_CANDIDATES = [
    (
        "/content/drive/MyDrive/SAM_Thesis/Processed/"
        "EN4_stratification_20_200m_monthly_2008_2025_SO.nc"
    ),
    (
        "/content/drive/MyDrive/SAM_Thesis/Processed/"
        "EN4_NetCDF_inventory/"
        "EN4_stratification_20_200m_monthly_2008_2025_SO.nc"
    ),
]

MLD_CANDIDATES = [
    (
        "/content/drive/MyDrive/SAM_Thesis/Processed/"
        "MLD_monthly_2008_2025_SO.nc"
    ),
    (
        "/content/drive/MyDrive/SAM_Thesis/Data/"
        "MLD_monthly_2008_2025_SO.nc"
    ),
]

SIC_CANDIDATES = [
    (
        "/content/drive/MyDrive/SAM_Thesis/Data/"
        "OSTIA_sea_ice_fraction_monthly_2008_2025_SO.nc"
    ),
]


def first_existing_path(candidates):

    for candidate in candidates:

        if os.path.exists(
            candidate
        ):

            return candidate

    return candidates[0]


FWC_FILE = first_existing_path(
    FWC_CANDIDATES
)

SSS_FILE = first_existing_path(
    SSS_CANDIDATES
)

STRATIFICATION_FILE = first_existing_path(
    STRATIFICATION_CANDIDATES
)

MLD_FILE = first_existing_path(
    MLD_CANDIDATES
)

SIC_FILE = first_existing_path(
    SIC_CANDIDATES
)


OUTPUT_DIRECTORY = Path(
    "/content/drive/MyDrive/SAM_Thesis/paper2"
)

OUTPUT_DIRECTORY.mkdir(
    parents=True,
    exist_ok=True,
)

OUT_SEASONAL_NC = (
    OUTPUT_DIRECTORY
    / "Fig09_EN4_freshwater_stratification_seasonal_fields.nc"
)


# 4. ANALYSIS SETTINGS
YEAR_START = 2008
YEAR_END = 2025

LAT_MIN = -90.0
LAT_MAX = -60.0

SEASON_ORDER = [
    "Spring",
    "Summer",
    "Autumn",
    "Winter",
]

SEASON_MONTHS = {
    "Spring": [
        9,
        10,
        11,
    ],
    "Summer": [
        12,
        1,
        2,
    ],
    "Autumn": [
        3,
        4,
        5,
    ],
    "Winter": [
        6,
        7,
        8,
    ],
}

SEASON_ROW_LABEL = {
    "Spring": "Spring",
    "Summer": "Summer",
    "Autumn": "Autumn",
    "Winter": "Winter",
}

COLUMN_LABELS = [
    "FWC\nclimatology",
    "SSS\nanomaly",
    "Stratification + MLD",
]

PANEL_LETTERS = [
    [
        "(a)",
        "(b)",
        "(c)",
    ],
    [
        "(d)",
        "(e)",
        "(f)",
    ],
    [
        "(g)",
        "(h)",
        "(i)",
    ],
    [
        "(j)",
        "(k)",
        "(l)",
    ],
]

MERIDIANS = np.arange(
    -180,
    181,
    30,
)

PARALLELS = [
    -60,
    -70,
    -80,
]

SIC_CONTOUR_PERCENT = 15.0

# Only levels falling within a seasonal field are plotted.
PREFERRED_MLD_LEVELS_M = [
    50,
    100,
    200,
    500,
    1000,
]

# The PDF is vector and should be used for publication.
# A 600-dpi PNG is the stable low-RAM option.
SAVE_PNG_DPI = 600

# Setting this to True may require substantially more RAM.
SAVE_OPTIONAL_1080_DPI = False

# 5. FIGURE SETTINGS
FIGSIZE = (
    9.2,
    8.4,
)

FONT_COLUMN_TITLE = 14.0
FONT_PANEL_LETTER = 17.0
FONT_ROW_LABEL = 15.5
FONT_GEO_LABEL = 7.6
FONT_CONTOUR_LABEL = 7.0
FONT_COLORBAR = 10.5
FONT_COLORBAR_TICK = 8.0
FONT_LEGEND = 8.8

LAND_COLOR = "0.82"
COAST_COLOR = "0.30"
GRID_COLOR = "0.62"
NO_DATA_COLOR = "0.86"

PANEL_LETTER_X = -0.12
PANEL_LETTER_Y = 1.025

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


# 6. INPUT CHECK
for required_file in [
    FWC_FILE,
    SSS_FILE,
    STRATIFICATION_FILE,
    MLD_FILE,
    SIC_FILE,
]:

    if not os.path.exists(
        required_file
    ):

        raise FileNotFoundError(
            f"Required file was not found:\n"
            f"{required_file}"
        )


print(
    "\nUsing files:"
)

print(
    "FWC  :",
    FWC_FILE,
)

print(
    "SSS  :",
    SSS_FILE,
)

print(
    "STRAT:",
    STRATIFICATION_FILE,
)

print(
    "MLD  :",
    MLD_FILE,
)

print(
    "SIC  :",
    SIC_FILE,
)


# 7. DATASET HELPERS
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
    Open without Dask chunks.

    This is intentional: monthly slices are loaded explicitly one
    at a time, preventing Dask from building a large task graph.
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
                f"{engine}: {error}"
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
    Detect time, latitude and longitude names.
    """
    aliases = {
        "time": [
            "time",
            "date",
            "datetime",
            "month",
            "valid_time",
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


def choose_data_variable(
    dataset,
    preferred_names,
):
    """
    Prefer named science variables and avoid known auxiliary fields.
    """
    excluded = {
        "valid_layer_mask",
        "gebco_water_depth",
        "time_bnds",
        "depth_bnds",
        "lat_bnds",
        "lon_bnds",
        "surface_salinity_reference",
    }

    available = [
        variable_name
        for variable_name
        in dataset.data_vars
        if variable_name
        not in excluded
    ]

    for preferred_name in preferred_names:

        for variable_name in available:

            if (
                preferred_name.lower()
                in variable_name.lower()
            ):

                return variable_name

    if not available:

        raise ValueError(
            "No suitable science variable was found."
        )

    return available[0]


def standardize_dataset(
    dataset,
):
    """
    Standardize names, longitude, latitude and time range lazily.
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
        lat=slice(
            LAT_MIN,
            LAT_MAX,
        )
    )

    dataset = dataset.sel(
        time=slice(
            f"{YEAR_START}-01-01",
            f"{YEAR_END}-12-31",
        )
    )

    return dataset


def prepare_monthly_field(
    field,
):
    """
    Return one 2-D lat × lon monthly field.
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
            "A monthly field still contains "
            f"unexpected dimensions: "
            f"{extra_dimensions}"
        )

    return field.transpose(
        "lat",
        "lon",
    )


def infer_multiplication_factor(
    data_array,
    role,
):
    """
    Infer only from metadata and ONE monthly slice.
    Never load the complete time series.
    """
    units = str(
        data_array.attrs.get(
            "units",
            "",
        )
    ).lower()

    if role == "mld":

        if (
            "kilomet" in units
            or units.strip()
            in [
                "km",
                "kilometer",
                "kilometers",
            ]
        ):

            return 1000.0

        return 1.0

    if role == "sic":

        if "%" in units:

            return 1.0

        first_slice = prepare_monthly_field(
            data_array.isel(
                time=0
            ).load()
        )

        sample_values = np.asarray(
            first_slice.values,
            dtype=np.float32,
        )

        finite_values = sample_values[
            np.isfinite(
                sample_values
            )
        ]

        del first_slice
        del sample_values

        gc.collect()

        if (
            finite_values.size > 0
            and np.nanmax(
                finite_values
            ) <= 1.01
        ):

            return 100.0

        return 1.0

    return 1.0


def grids_match(
    data_array,
    target_latitude,
    target_longitude,
):
    """
    Check whether a monthly field is already on the target grid.
    """
    return (
        data_array.sizes.get(
            "lat"
        )
        == target_latitude.size
        and data_array.sizes.get(
            "lon"
        )
        == target_longitude.size
        and np.allclose(
            data_array[
                "lat"
            ].values,
            target_latitude,
        )
        and np.allclose(
            data_array[
                "lon"
            ].values,
            target_longitude,
        )
    )


def monthly_field_to_target_numpy(
    monthly_field,
    target_latitude,
    target_longitude,
    multiplication_factor=1.0,
):
    """
    Load and, when required, interpolate ONE monthly field only.
    """
    monthly_field = prepare_monthly_field(
        monthly_field
    )

    if not grids_match(
        monthly_field,
        target_latitude,
        target_longitude,
    ):

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

    values = np.asarray(
        monthly_field.values,
        dtype=np.float32,
    )

    if multiplication_factor != 1.0:

        values = (
            values
            * np.float32(
                multiplication_factor
            )
        )

    return values


def month_to_season(
    month,
):

    for season_name, season_months in (
        SEASON_MONTHS.items()
    ):

        if month in season_months:

            return season_name

    raise ValueError(
        f"Month {month} is not assigned "
        "to an austral season."
    )


def sequential_seasonal_aggregation(
    file_path,
    preferred_variables,
    target_latitude,
    target_longitude,
    role,
    calculate_all_month_mean=False,
):
    """
    Process one source file at a time and one month at a time.

    Memory at any instant contains:
      • one monthly source field;
      • one monthly target-grid field;
      • small target-grid accumulators.
    """
    dataset = open_dataset_safely(
        file_path
    )

    dataset = standardize_dataset(
        dataset
    )

    variable_name = choose_data_variable(
        dataset,
        preferred_variables,
    )

    data_array = dataset[
        variable_name
    ]

    multiplication_factor = (
        infer_multiplication_factor(
            data_array,
            role,
        )
    )

    times = pd.DatetimeIndex(
        pd.to_datetime(
            data_array[
                "time"
            ].values
        )
    )

    target_shape = (
        target_latitude.size,
        target_longitude.size,
    )

    seasonal_sum = {
        season: np.zeros(
            target_shape,
            dtype=np.float64,
        )
        for season in SEASON_ORDER
    }

    seasonal_count = {
        season: np.zeros(
            target_shape,
            dtype=np.uint16,
        )
        for season in SEASON_ORDER
    }

    if calculate_all_month_mean:

        all_month_sum = np.zeros(
            target_shape,
            dtype=np.float64,
        )

        all_month_count = np.zeros(
            target_shape,
            dtype=np.uint16,
        )

    else:

        all_month_sum = None
        all_month_count = None

    number_of_times = len(
        times
    )

    print(
        f"\nProcessing {role}: "
        f"{number_of_times} monthly records"
    )

    for time_index, timestamp in enumerate(
        times
    ):

        season = month_to_season(
            int(
                timestamp.month
            )
        )

        monthly_values = (
            monthly_field_to_target_numpy(
                data_array.isel(
                    time=time_index
                ),
                target_latitude,
                target_longitude,
                multiplication_factor,
            )
        )

        finite = np.isfinite(
            monthly_values
        )

        seasonal_sum[
            season
        ][finite] += monthly_values[
            finite
        ].astype(
            np.float64
        )

        seasonal_count[
            season
        ][finite] += 1

        if calculate_all_month_mean:

            all_month_sum[
                finite
            ] += monthly_values[
                finite
            ].astype(
                np.float64
            )

            all_month_count[
                finite
            ] += 1

        del monthly_values
        del finite

        if (
            (time_index + 1) % 12 == 0
            or time_index
            == number_of_times - 1
        ):

            print(
                f"  [{time_index + 1:03d}/"
                f"{number_of_times:03d}] "
                f"{timestamp:%Y-%m}"
            )

            gc.collect()

    seasonal_mean = {}

    for season in SEASON_ORDER:

        with np.errstate(
            invalid="ignore",
            divide="ignore",
        ):

            mean_values = np.divide(
                seasonal_sum[
                    season
                ],
                seasonal_count[
                    season
                ],
                out=np.full(
                    target_shape,
                    np.nan,
                    dtype=np.float64,
                ),
                where=(
                    seasonal_count[
                        season
                    ]
                    > 0
                ),
            )

        seasonal_mean[
            season
        ] = mean_values.astype(
            np.float32
        )

    if calculate_all_month_mean:

        with np.errstate(
            invalid="ignore",
            divide="ignore",
        ):

            all_month_mean = np.divide(
                all_month_sum,
                all_month_count,
                out=np.full(
                    target_shape,
                    np.nan,
                    dtype=np.float64,
                ),
                where=(
                    all_month_count
                    > 0
                ),
            ).astype(
                np.float32
            )

    else:

        all_month_mean = None

    units = str(
        data_array.attrs.get(
            "units",
            "",
        )
    )

    dataset.close()

    del dataset
    del data_array
    del seasonal_sum
    del seasonal_count
    del all_month_sum
    del all_month_count

    gc.collect()

    return {
        "variable_name": variable_name,
        "units": units,
        "seasonal_mean": seasonal_mean,
        "all_month_mean": all_month_mean,
        "multiplication_factor": (
            multiplication_factor
        ),
        "first_time": times[0],
        "last_time": times[-1],
        "number_of_months": len(times),
    }


def data_array_from_numpy(
    values,
    latitude,
    longitude,
    name,
    units,
):
    """
    Construct a small target-grid DataArray.
    """
    output = xr.DataArray(
        values,
        dims=(
            "lat",
            "lon",
        ),
        coords={
            "lat": latitude,
            "lon": longitude,
        },
        name=name,
    )

    output.attrs[
        "units"
    ] = units

    return output


# 8. READ TARGET EN4 GRID ONLY
target_dataset = open_dataset_safely(
    STRATIFICATION_FILE
)

target_dataset = standardize_dataset(
    target_dataset
)

TARGET_LATITUDE = np.asarray(
    target_dataset[
        "lat"
    ].values,
    dtype=np.float64,
)

TARGET_LONGITUDE = np.asarray(
    target_dataset[
        "lon"
    ].values,
    dtype=np.float64,
)

target_dataset.close()

del target_dataset

gc.collect()

print(
    "\nTarget EN4 grid:"
)

print(
    f"Latitude cells : "
    f"{TARGET_LATITUDE.size}"
)

print(
    f"Longitude cells: "
    f"{TARGET_LONGITUDE.size}"
)


# 9. SEQUENTIAL LOW-RAM AGGREGATION
fwc_result = sequential_seasonal_aggregation(
    file_path=FWC_FILE,
    preferred_variables=[
        "freshwater_content",
        "fwc",
    ],
    target_latitude=TARGET_LATITUDE,
    target_longitude=TARGET_LONGITUDE,
    role="fwc",
    calculate_all_month_mean=False,
)

sss_result = sequential_seasonal_aggregation(
    file_path=SSS_FILE,
    preferred_variables=[
        "surface_salinity",
        "sss",
        "salinity",
    ],
    target_latitude=TARGET_LATITUDE,
    target_longitude=TARGET_LONGITUDE,
    role="sss",
    calculate_all_month_mean=True,
)

stratification_result = (
    sequential_seasonal_aggregation(
        file_path=STRATIFICATION_FILE,
        preferred_variables=[
            "stratification_20_200m",
            "stratification",
            "density_difference",
        ],
        target_latitude=TARGET_LATITUDE,
        target_longitude=TARGET_LONGITUDE,
        role="stratification",
        calculate_all_month_mean=False,
    )
)

mld_result = sequential_seasonal_aggregation(
    file_path=MLD_FILE,
    preferred_variables=[
        "MLD",
        "mld",
        "mixed_layer_depth",
        "mlotst",
    ],
    target_latitude=TARGET_LATITUDE,
    target_longitude=TARGET_LONGITUDE,
    role="mld",
    calculate_all_month_mean=False,
)

sic_result = sequential_seasonal_aggregation(
    file_path=SIC_FILE,
    preferred_variables=[
        "SIF",
        "sif",
        "sea_ice_fraction",
        "sic",
    ],
    target_latitude=TARGET_LATITUDE,
    target_longitude=TARGET_LONGITUDE,
    role="sic",
    calculate_all_month_mean=False,
)


# 10. BUILD SMALL SEASONAL DATA ARRAYS
fwc_seasonal = {}
sss_anomaly_seasonal = {}
stratification_seasonal = {}
mld_seasonal = {}
sic_seasonal = {}

for season in SEASON_ORDER:

    fwc_seasonal[
        season
    ] = data_array_from_numpy(
        fwc_result[
            "seasonal_mean"
        ][
            season
        ],
        TARGET_LATITUDE,
        TARGET_LONGITUDE,
        "freshwater_content",
        "m",
    )

    sss_anomaly_values = (
        sss_result[
            "seasonal_mean"
        ][
            season
        ]
        - sss_result[
            "all_month_mean"
        ]
    ).astype(
        np.float32
    )

    sss_anomaly_seasonal[
        season
    ] = data_array_from_numpy(
        sss_anomaly_values,
        TARGET_LATITUDE,
        TARGET_LONGITUDE,
        "surface_salinity_anomaly",
        "1",
    )

    stratification_seasonal[
        season
    ] = data_array_from_numpy(
        stratification_result[
            "seasonal_mean"
        ][
            season
        ],
        TARGET_LATITUDE,
        TARGET_LONGITUDE,
        "stratification_20_200m",
        "kg m-3",
    )

    mld_seasonal[
        season
    ] = data_array_from_numpy(
        mld_result[
            "seasonal_mean"
        ][
            season
        ],
        TARGET_LATITUDE,
        TARGET_LONGITUDE,
        "mixed_layer_depth",
        "m",
    )

    sic_seasonal[
        season
    ] = data_array_from_numpy(
        sic_result[
            "seasonal_mean"
        ][
            season
        ],
        TARGET_LATITUDE,
        TARGET_LONGITUDE,
        "sea_ice_fraction",
        "%",
    )


# Release the aggregation dictionaries.
del fwc_result
del sss_result
del stratification_result
del mld_result
del sic_result

gc.collect()


# 11. COLOUR LIMIT HELPERS
def collect_finite_values(
    data_dictionary,
):

    collected = []

    for season in SEASON_ORDER:

        values = np.asarray(
            data_dictionary[
                season
            ].values,
            dtype=np.float32,
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
                0.0,
                1.0,
            ],
            dtype=np.float32,
        )

    return np.concatenate(
        collected
    )


def nice_upper_limit(
    value,
):

    if (
        not np.isfinite(
            value
        )
        or value <= 0
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


fwc_values = collect_finite_values(
    fwc_seasonal
)

sss_values = collect_finite_values(
    sss_anomaly_seasonal
)

stratification_values = collect_finite_values(
    stratification_seasonal
)

# FWC can be slightly negative under the selected reference
# salinity. Preserve negative values when present.
fwc_lower_percentile = float(
    np.nanpercentile(
        fwc_values,
        1.0,
    )
)

fwc_upper_percentile = float(
    np.nanpercentile(
        fwc_values,
        99.0,
    )
)

FWC_VMIN = min(
    0.0,
    fwc_lower_percentile,
)

FWC_VMAX = nice_upper_limit(
    fwc_upper_percentile
)

SSS_LIMIT = nice_upper_limit(
    float(
        np.nanpercentile(
            np.abs(
                sss_values
            ),
            99.0,
        )
    )
)

STRATIFICATION_VMIN = min(
    0.0,
    float(
        np.nanpercentile(
            stratification_values,
            1.0,
        )
    ),
)

STRATIFICATION_VMAX = nice_upper_limit(
    float(
        np.nanpercentile(
            stratification_values,
            99.0,
        )
    )
)

print(
    "\nFinal colour limits:"
)

print(
    f"FWC            : "
    f"{FWC_VMIN:.4g} to "
    f"{FWC_VMAX:.4g} m"
)

print(
    f"SSS anomaly    : "
    f"±{SSS_LIMIT:.4g}"
)

print(
    f"Stratification : "
    f"{STRATIFICATION_VMIN:.4g} to "
    f"{STRATIFICATION_VMAX:.4g} "
    f"kg m-3"
)


# 12. SAVE SMALL SEASONAL NETCDF
season_coordinate = xr.DataArray(
    SEASON_ORDER,
    dims="season",
    name="season",
)

seasonal_output = xr.Dataset(
    {
        "freshwater_content_climatology": xr.concat(
            [
                fwc_seasonal[
                    season
                ]
                for season in SEASON_ORDER
            ],
            dim=season_coordinate,
        ),
        "surface_salinity_anomaly": xr.concat(
            [
                sss_anomaly_seasonal[
                    season
                ]
                for season in SEASON_ORDER
            ],
            dim=season_coordinate,
        ),
        "stratification_climatology": xr.concat(
            [
                stratification_seasonal[
                    season
                ]
                for season in SEASON_ORDER
            ],
            dim=season_coordinate,
        ),
        "mixed_layer_depth_climatology": xr.concat(
            [
                mld_seasonal[
                    season
                ]
                for season in SEASON_ORDER
            ],
            dim=season_coordinate,
        ),
        "sea_ice_fraction_climatology": xr.concat(
            [
                sic_seasonal[
                    season
                ]
                for season in SEASON_ORDER
            ],
            dim=season_coordinate,
        ),
    }
)

seasonal_output.attrs.update(
    {
        "title": (
            "Figure 9 EN4 freshwater and "
            "stratification seasonal fields"
        ),
        "period": (
            "2008-01 to 2025-12"
        ),
        "SSS_anomaly_reference": (
            "2008-2025 all-month grid-cell mean"
        ),
        "SIC_contour_level_percent": (
            SIC_CONTOUR_PERCENT
        ),
    }
)

encoding = {
    variable_name: {
        "dtype": "float32",
        "zlib": True,
        "complevel": 4,
        "_FillValue": np.float32(
            9.96921e36
        ),
    }
    for variable_name
    in seasonal_output.data_vars
}

seasonal_output.rename(
    {
        "lat": "latitude",
        "lon": "longitude",
    }
).to_netcdf(
    OUT_SEASONAL_NC,
    encoding=encoding,
)

seasonal_output.close()

del seasonal_output

gc.collect()


# PART B — LATEST / FINAL FIGURE 9 PLOTTING
# Reads the seasonal NetCDF created in Part A. Monthly fields are not
# recomputed in this section.

# FIGURE 9 — EN4 FRESHWATER STORAGE AND UPPER-OCEAN STABILITY
# PLOT-ONLY WORKFLOW USING THE PRECOMPUTED SEASONAL NETCDF
#
# No monthly fields are recomputed.
#
# Input
# /content/drive/MyDrive/SAM_Thesis/paper2/
# Fig09_EN4_freshwater_stratification_seasonal_fields.nc
#
# Figure structure
# Rows:
#   Spring, Summer, Autumn, Winter
#
# Columns:
#   1. FWC seasonal climatology
#   2. SSS seasonal anomaly
#   3. Stratification seasonal climatology
#      + MLD solid contours
#      + SIC 15% dashed contour
#
# Outputs
# 1080-dpi PNG
# Vector PDF
#
# Results section
# 3.7 Freshwater storage and upper-ocean stability


# 3. INPUT AND OUTPUT PATHS
SEASONAL_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/paper2/"
    "Fig09_EN4_freshwater_stratification_seasonal_fields.nc"
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
    / "Figure09_EN4_freshwater_stratification_"
      "FINAL_1080dpi.png"
)

OUT_PDF = (
    OUTPUT_DIRECTORY
    / "Figure09_EN4_freshwater_stratification_FINAL.pdf"
)

OUT_SUMMARY = (
    OUTPUT_DIRECTORY
    / "Figure09_EN4_freshwater_stratification_"
      "FINAL_summary.txt"
)


# 4. ANALYSIS AND DISPLAY SETTINGS
LAT_MIN = -90.0
LAT_MAX = -60.0

SEASON_ORDER = [
    "Spring",
    "Summer",
    "Autumn",
    "Winter",
]

SEASON_LABELS = {
    "Spring": "Spring",
    "Summer": "Summer",
    "Autumn": "Autumn",
    "Winter": "Winter",
}

COLUMN_TITLES = [
    "FWC\nclimatology",
    "SSS\nanomaly",
    "Stratification + MLD",
]

PANEL_LETTERS = [
    ["(a)", "(b)", "(c)"],
    ["(d)", "(e)", "(f)"],
    ["(g)", "(h)", "(i)"],
    ["(j)", "(k)", "(l)"],
]

MERIDIANS = np.arange(
    -180,
    181,
    30,
)

PARALLELS = [
    -60,
    -70,
    -80,
]

SIC_CONTOUR_PERCENT = 15.0

# Scientifically interpretable MLD contours.
# Only levels within a seasonal field's actual range are plotted.
PREFERRED_MLD_LEVELS_M = [
    50,
    100,
    200,
    500,
    1000,
]

SAVE_DPI = 1080


# 5. COMPACT 1080-DPI FIGURE SETTINGS
# The compact physical canvas keeps the 1080-dpi raster below the
# very large memory footprint produced by a 12–16 inch-wide figure.
FIGSIZE = (
    8.8,
    8.4,
)

FONT_COLUMN_TITLE = 14.5
FONT_PANEL_LETTER = 17.5
FONT_ROW_LABEL = 15.5
FONT_GEO_LABEL = 5.8
FONT_CONTOUR_LABEL = 6.8
FONT_COLORBAR_LABEL = 10.7
FONT_COLORBAR_TICK = 8.0
FONT_LEGEND = 9.0

LAND_COLOR = "0.80"
COAST_COLOR = "0.30"
GRID_COLOR = "0.62"
NO_DATA_COLOR = "0.86"

# The panel letter sits left of the circular map and no longer
# touches the 30°W longitude label.
PANEL_LETTER_X = -0.19
PANEL_LETTER_Y = 1.025

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


# 6. DATASET HELPERS
def open_dataset_safely(
    file_path,
):
    """
    Open the small seasonal file without Dask.
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
                kwargs["engine"] = engine

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
        f"Could not open:\n{file_path}\n\n"
        + "\n".join(attempts)
    )


def detect_coordinate(
    dataset,
    coordinate_type,
):
    """
    Detect season, latitude and longitude coordinates.
    """
    aliases = {
        "season": [
            "season",
            "seasons",
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

        if coordinate_type == "season":
            if lower_name in aliases["season"]:
                return name

        elif coordinate_type == "lat":
            if (
                lower_name in aliases["lat"]
                or standard_name == "latitude"
                or axis == "Y"
                or "degrees_north" in units
            ):
                return name

        elif coordinate_type == "lon":
            if (
                lower_name in aliases["lon"]
                or standard_name == "longitude"
                or axis == "X"
                or "degrees_east" in units
            ):
                return name

    raise KeyError(
        f"Could not detect the {coordinate_type} coordinate."
    )


def choose_variable(
    dataset,
    exact_candidates,
    contains_candidates,
):
    """
    Choose a variable robustly from the seasonal output.
    """
    available = list(
        dataset.data_vars
    )

    for candidate in exact_candidates:
        if candidate in available:
            return candidate

    for candidate in contains_candidates:
        for variable_name in available:
            if candidate.lower() in variable_name.lower():
                return variable_name

    raise KeyError(
        "Could not find a required variable.\n"
        f"Available variables: {available}"
    )


def standardize_dataset(
    dataset,
):
    """
    Rename season/latitude/longitude and normalize longitude.
    """
    season_name = detect_coordinate(
        dataset,
        "season",
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

    if season_name != "season":
        rename_mapping[season_name] = "season"

    if lat_name != "lat":
        rename_mapping[lat_name] = "lat"

    if lon_name != "lon":
        rename_mapping[lon_name] = "lon"

    if rename_mapping:
        dataset = dataset.rename(
            rename_mapping
        )

    normalized_longitude = (
        (
            dataset["lon"].astype(float)
            + 180.0
        )
         % 360.0
    ) - 180.0

    dataset = dataset.assign_coords(
        lon=normalized_longitude
    )

    longitude_values = np.asarray(
        dataset["lon"].values
    )

    _, unique_indices = np.unique(
        longitude_values,
        return_index=True,
    )

    dataset = dataset.isel(
        lon=np.sort(unique_indices)
    )

    dataset = dataset.sortby(
        "lon"
    )

    dataset = dataset.sortby(
        "lat"
    )

    dataset = dataset.sel(
        lat=slice(
            LAT_MIN,
            LAT_MAX,
        )
    )

    # Normalize season coordinate values to title case strings.
    normalized_seasons = [
        str(value).strip().title()
        for value in dataset["season"].values
    ]

    dataset = dataset.assign_coords(
        season=normalized_seasons
    )

    return dataset


def extract_season(
    data_array,
    season,
):
    """
    Select and fully load one small 2-D seasonal field.
    """
    field = data_array.sel(
        season=season
    ).squeeze(
        drop=True
    )

    extra_dimensions = [
        dimension
        for dimension in field.dims
        if dimension not in [
            "lat",
            "lon",
        ]
    ]

    if extra_dimensions:
        raise ValueError(
            f"{data_array.name} contains unexpected "
            f"dimensions after season selection: "
            f"{extra_dimensions}"
        )

    return field.transpose(
        "lat",
        "lon",
    ).load()


def ensure_sic_percent(
    sic_dictionary,
):
    """
    Convert SIC from 0–1 to percent when necessary.
    """
    all_finite = []

    for season in SEASON_ORDER:
        values = np.asarray(
            sic_dictionary[season].values,
            dtype=np.float32,
        )

        finite = values[
            np.isfinite(values)
        ]

        if finite.size > 0:
            all_finite.append(
                finite
            )

    if not all_finite:
        return sic_dictionary

    maximum = float(
        np.nanmax(
            np.concatenate(all_finite)
        )
    )

    if maximum <= 1.01:
        print(
            "SIC detected as fraction; converting to percent."
        )

        return {
            season: (
                sic_dictionary[season]
                * 100.0
            )
            for season in SEASON_ORDER
        }

    return sic_dictionary

# 7. OPEN AND LOAD ONLY THE PRECOMPUTED SEASONAL FILE
if not os.path.exists(
    SEASONAL_FILE
):
    raise FileNotFoundError(
        "The precomputed seasonal NetCDF was not found:\n"
        f"{SEASONAL_FILE}"
    )

dataset = open_dataset_safely(
    SEASONAL_FILE
)

dataset = standardize_dataset(
    dataset
)

fwc_variable = choose_variable(
    dataset,
    exact_candidates=[
        "freshwater_content_climatology",
    ],
    contains_candidates=[
        "freshwater_content",
        "fwc",
    ],
)

sss_variable = choose_variable(
    dataset,
    exact_candidates=[
        "surface_salinity_anomaly",
    ],
    contains_candidates=[
        "salinity_anomaly",
        "sss_anomaly",
    ],
)

stratification_variable = choose_variable(
    dataset,
    exact_candidates=[
        "stratification_climatology",
    ],
    contains_candidates=[
        "stratification",
    ],
)

mld_variable = choose_variable(
    dataset,
    exact_candidates=[
        "mixed_layer_depth_climatology",
    ],
    contains_candidates=[
        "mixed_layer_depth",
        "mld",
    ],
)

sic_variable = choose_variable(
    dataset,
    exact_candidates=[
        "sea_ice_fraction_climatology",
    ],
    contains_candidates=[
        "sea_ice_fraction",
        "sic",
    ],
)

print(
    "\nDetected seasonal variables:"
)

print(
    "FWC            :",
    fwc_variable,
)

print(
    "SSS anomaly    :",
    sss_variable,
)

print(
    "Stratification :",
    stratification_variable,
)

print(
    "MLD            :",
    mld_variable,
)

print(
    "SIC            :",
    sic_variable,
)

available_seasons = [
    str(value).strip().title()
    for value in dataset["season"].values
]

missing_seasons = [
    season
    for season in SEASON_ORDER
    if season not in available_seasons
]

if missing_seasons:
    raise ValueError(
        "The seasonal file does not contain: "
        + ", ".join(missing_seasons)
    )

fwc = {
    season: extract_season(
        dataset[fwc_variable],
        season,
    )
    for season in SEASON_ORDER
}

sss_anomaly = {
    season: extract_season(
        dataset[sss_variable],
        season,
    )
    for season in SEASON_ORDER
}

stratification = {
    season: extract_season(
        dataset[stratification_variable],
        season,
    )
    for season in SEASON_ORDER
}

mld = {
    season: extract_season(
        dataset[mld_variable],
        season,
    )
    for season in SEASON_ORDER
}

sic = {
    season: extract_season(
        dataset[sic_variable],
        season,
    )
    for season in SEASON_ORDER
}

sic = ensure_sic_percent(
    sic
)

dataset.close()

del dataset

gc.collect()


# 8. ROBUST COLUMN-SPECIFIC COLOR LIMITS
def combined_finite_values(
    seasonal_dictionary,
):
    """
    Combine finite values from only four small seasonal fields.
    """
    collected = []

    for season in SEASON_ORDER:
        values = np.asarray(
            seasonal_dictionary[season].values,
            dtype=np.float32,
        )

        finite = values[
            np.isfinite(values)
        ]

        if finite.size > 0:
            collected.append(
                finite
            )

    if not collected:
        raise ValueError(
            "No finite seasonal values were available."
        )

    return np.concatenate(
        collected
    )


def nice_upper_bound(
    value,
):
    """
    Round upward to a readable colorbar endpoint.
    """
    if not np.isfinite(value) or value <= 0:
        return 1.0

    exponent = np.floor(
        np.log10(value)
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


fwc_values = combined_finite_values(
    fwc
)

sss_values = combined_finite_values(
    sss_anomaly
)

stratification_values = combined_finite_values(
    stratification
)

FWC_VMIN = min(
    0.0,
    float(
        np.nanpercentile(
            fwc_values,
            1.0,
        )
    ),
)

FWC_VMAX = nice_upper_bound(
    float(
        np.nanpercentile(
            fwc_values,
            99.0,
        )
    )
)

SSS_LIMIT = nice_upper_bound(
    float(
        np.nanpercentile(
            np.abs(sss_values),
            99.0,
        )
    )
)

STRATIFICATION_VMIN = min(
    0.0,
    float(
        np.nanpercentile(
            stratification_values,
            1.0,
        )
    ),
)

STRATIFICATION_VMAX = nice_upper_bound(
    float(
        np.nanpercentile(
            stratification_values,
            99.0,
        )
    )
)

print(
    "\nColor limits:"
)

print(
    f"FWC            : "
    f"{FWC_VMIN:.4g} to {FWC_VMAX:.4g} m"
)

print(
    f"SSS anomaly    : "
    f"±{SSS_LIMIT:.4g}"
)

print(
    f"Stratification : "
    f"{STRATIFICATION_VMIN:.4g} to "
    f"{STRATIFICATION_VMAX:.4g} kg m-3"
)


# 9. MAP HELPERS
def circular_boundary():
    """
    Circular clipping path for the South Polar map.
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


def format_longitude(
    longitude,
):
    """
    Format polar-rim longitude labels.
    """
    value = int(
        longitude
    )

    if value == 0:
        return "0°"

    if abs(value) == 180:
        return "180°"

    if value < 0:
        return f"{abs(value)}°W"

    return f"{value}°E"


def style_polar_axis(
    axis,
):
    """
    Apply the common professional South Polar map style.
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
        linewidth=0.55,
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
        linewidth=0.42,
        linestyle=":",
        color=GRID_COLOR,
        alpha=0.78,
        zorder=3,
    )

    # Rim labels are outside the circle but columns are spaced far
    # enough that neighboring labels cannot collide.
    for longitude in MERIDIANS:
        axis.text(
            longitude,
            -58.1,
            format_longitude(
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
        fontsize=FONT_GEO_LABEL + 1.0,
        fontweight="bold",
        color="black",
        zorder=30,
    )


def cyclic_field(
    field,
):
    """
    Add the cyclic longitude column to remove the map seam.
    """
    cyclic_values, cyclic_longitude = add_cyclic_point(
        np.asarray(
            field.values,
            dtype=np.float32,
        ),
        coord=np.asarray(
            field["lon"].values,
            dtype=np.float64,
        ),
        axis=-1,
    )

    return (
        cyclic_values,
        cyclic_longitude,
        np.asarray(
            field["lat"].values,
            dtype=np.float64,
        ),
    )


def available_mld_levels(
    mld_field,
):
    """
    Keep only preferred MLD levels inside the seasonal data range.
    """
    values = np.asarray(
        mld_field.values,
        dtype=np.float32,
    )

    finite = values[
        np.isfinite(values)
    ]

    if finite.size == 0:
        return []

    minimum = float(
        np.nanmin(finite)
    )

    maximum = float(
        np.nanmax(finite)
    )

    return [
        level
        for level in PREFERRED_MLD_LEVELS_M
        if minimum <= level <= maximum
    ]


def plot_panel(
    axis,
    field,
    panel_letter,
    cmap,
    vmin,
    vmax,
    mld_field=None,
    sic_field=None,
):
    """
    Plot one seasonal panel.
    """
    style_polar_axis(
        axis
    )

    field_values, cyclic_lon, latitude = cyclic_field(
        field
    )

    no_data_values = np.where(
        np.isfinite(field_values),
        np.nan,
        1.0,
    )

    axis.pcolormesh(
        cyclic_lon,
        latitude,
        no_data_values,
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
        cyclic_lon,
        latitude,
        field_values,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="auto",
        rasterized=True,
        zorder=4,
    )

    # Overlay contours only inside valid stratification cells.
    if mld_field is not None:
        masked_mld = mld_field.where(
            np.isfinite(field)
        )

        levels = available_mld_levels(
            masked_mld
        )

        if levels:
            mld_values, mld_lon, mld_lat = cyclic_field(
                masked_mld
            )

            mld_contours = axis.contour(
                mld_lon,
                mld_lat,
                mld_values,
                levels=levels,
                colors="black",
                linewidths=0.82,
                linestyles="-",
                transform=ccrs.PlateCarree(),
                zorder=14,
            )

            # Compact labels reduce the heavy clutter visible in the
            # earlier version.
            axis.clabel(
                mld_contours,
                inline=True,
                inline_spacing=2,
                fontsize=FONT_CONTOUR_LABEL,
                fmt=lambda value: (
                    f"{int(round(value))}"
                ),
                rightside_up=True,
                use_clabeltext=True,
            )

    if sic_field is not None:
        masked_sic = sic_field.where(
            np.isfinite(field)
        )

        sic_values, sic_lon, sic_lat = cyclic_field(
            masked_sic
        )

        finite_sic = sic_values[
            np.isfinite(sic_values)
        ]

        if (
            finite_sic.size > 0
            and np.nanmin(finite_sic)
            <= SIC_CONTOUR_PERCENT
            <= np.nanmax(finite_sic)
        ):
            axis.contour(
                sic_lon,
                sic_lat,
                sic_values,
                levels=[
                    SIC_CONTOUR_PERCENT
                ],
                colors="0.18",
                linewidths=0.95,
                linestyles="--",
                transform=ccrs.PlateCarree(),
                zorder=15,
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


# 10. BUILD THE PROFESSIONAL 4 × 3 LAYOUT
projection = ccrs.SouthPolarStereo(
    central_longitude=0.0
)

fig = plt.figure(
    figsize=FIGSIZE
)

# Explicit physical-square axes preserve round polar maps while
# allowing wide, controlled separation among all three columns.
MAP_HEIGHT = 0.180

MAP_WIDTH = (
    MAP_HEIGHT
    * FIGSIZE[1]
    / FIGSIZE[0]
)

# Fixed column centers give substantially more separation than the
# previous low-RAM figure and eliminate longitude-label collisions.
MAP_COLUMN_CENTERS = [
    0.205,
    0.500,
    0.795,
]

MAP_COLUMN_LEFTS = [
    center
    - MAP_WIDTH / 2.0
    for center in MAP_COLUMN_CENTERS
]

MAP_ROW_BOTTOMS = [
    0.744,
    0.539,
    0.334,
    0.129,
]

map_axes = np.empty(
    (
        4,
        3,
    ),
    dtype=object,
)

last_images = [
    None,
    None,
    None,
]

field_dictionaries = [
    fwc,
    sss_anomaly,
    stratification,
]

column_cmaps = [
    plt.get_cmap(
        "Blues"
    ),
    plt.get_cmap(
        "RdBu_r"
    ),
    plt.get_cmap(
        "turbo"
    ),
]

column_vmins = [
    FWC_VMIN,
    -SSS_LIMIT,
    STRATIFICATION_VMIN,
]

column_vmaxs = [
    FWC_VMAX,
    SSS_LIMIT,
    STRATIFICATION_VMAX,
]

for row_index, season in enumerate(
    SEASON_ORDER
):
    for column_index in range(
        3
    ):
        axis = fig.add_axes(
            [
                MAP_COLUMN_LEFTS[
                    column_index
                ],
                MAP_ROW_BOTTOMS[
                    row_index
                ],
                MAP_WIDTH,
                MAP_HEIGHT,
            ],
            projection=projection,
        )

        map_axes[
            row_index,
            column_index,
        ] = axis

        if column_index == 2:
            mld_overlay = mld[
                season
            ]

            sic_overlay = sic[
                season
            ]

        else:
            mld_overlay = None
            sic_overlay = None

        last_images[
            column_index
        ] = plot_panel(
            axis=axis,
            field=field_dictionaries[
                column_index
            ][
                season
            ],
            panel_letter=PANEL_LETTERS[
                row_index
            ][
                column_index
            ],
            cmap=column_cmaps[
                column_index
            ],
            vmin=column_vmins[
                column_index
            ],
            vmax=column_vmaxs[
                column_index
            ],
            mld_field=mld_overlay,
            sic_field=sic_overlay,
        )


# 11. COLUMN TITLES
COLUMN_TITLE_Y = 0.988

for column_index in range(
    3
):
    fig.text(
        MAP_COLUMN_CENTERS[
            column_index
        ],
        COLUMN_TITLE_Y,
        COLUMN_TITLES[
            column_index
        ],
        ha="center",
        va="top",
        fontsize=FONT_COLUMN_TITLE,
        fontweight="bold",
        linespacing=0.98,
        color="black",
    )


# 12. SEASON LABELS
SEASON_LABEL_X = 0.040

for row_index, season in enumerate(
    SEASON_ORDER
):
    fig.text(
        SEASON_LABEL_X,
        (
            MAP_ROW_BOTTOMS[
                row_index
            ]
            + MAP_HEIGHT / 2.0
        ),
        SEASON_LABELS[
            season
        ],
        rotation=90,
        rotation_mode="anchor",
        ha="center",
        va="center",
        fontsize=FONT_ROW_LABEL,
        fontweight="bold",
        color="black",
    )


# 13. COLORBAR HELPERS
def decimal_formatter(
    maximum_absolute_value,
):
    """
    Choose a compact colorbar label format.
    """
    if maximum_absolute_value < 0.1:
        return FuncFormatter(
            lambda value, position: (
                f"{value:.3f}"
            )
        )

    if maximum_absolute_value < 1.0:
        return FuncFormatter(
            lambda value, position: (
                f"{value:.2f}"
            )
        )

    if maximum_absolute_value < 10.0:
        return FuncFormatter(
            lambda value, position: (
                f"{value:.1f}"
            )
        )

    return FuncFormatter(
        lambda value, position: (
            f"{value:.0f}"
        )
    )


# 14. BOTTOM COLORBARS
COLORBAR_Y = 0.061
COLORBAR_HEIGHT = 0.012
COLORBAR_WIDTH = MAP_WIDTH * 0.94

colorbar_axes = [
    fig.add_axes(
        [
            MAP_COLUMN_CENTERS[
                column_index
            ]
            - COLORBAR_WIDTH / 2.0,
            COLORBAR_Y,
            COLORBAR_WIDTH,
            COLORBAR_HEIGHT,
        ]
    )
    for column_index in range(3)
]

colorbar_labels = [
    "FWC (m)",
    "SSS anomaly",
    "Stratification (kg m$^{-3}$)",
]

colorbar_ticks = [
    np.linspace(
        FWC_VMIN,
        FWC_VMAX,
        5,
    ),
    np.linspace(
        -SSS_LIMIT,
        SSS_LIMIT,
        5,
    ),
    np.linspace(
        STRATIFICATION_VMIN,
        STRATIFICATION_VMAX,
        5,
    ),
]

colorbar_formatters = [
    decimal_formatter(
        max(
            abs(FWC_VMIN),
            abs(FWC_VMAX),
        )
    ),
    decimal_formatter(
        SSS_LIMIT
    ),
    decimal_formatter(
        max(
            abs(STRATIFICATION_VMIN),
            abs(STRATIFICATION_VMAX),
        )
    ),
]

for column_index, colorbar_axis in enumerate(
    colorbar_axes
):
    if column_index == 1:
        extend = "both"

    elif (
        column_index == 0
        and FWC_VMIN < 0.0
    ):
        extend = "both"

    elif (
        column_index == 2
        and STRATIFICATION_VMIN < 0.0
    ):
        extend = "both"

    else:
        extend = "max"

    colorbar = fig.colorbar(
        last_images[
            column_index
        ],
        cax=colorbar_axis,
        orientation="horizontal",
        extend=extend,
    )

    colorbar.set_ticks(
        colorbar_ticks[
            column_index
        ]
    )

    colorbar.ax.xaxis.set_major_formatter(
        colorbar_formatters[
            column_index
        ]
    )

    colorbar.set_label(
        colorbar_labels[
            column_index
        ],
        fontsize=FONT_COLORBAR_LABEL,
        fontweight="bold",
        labelpad=4,
    )

    colorbar.ax.tick_params(
        labelsize=FONT_COLORBAR_TICK,
        width=0.9,
        length=3,
        pad=2,
    )

    for tick_label in (
        colorbar.ax.get_xticklabels()
    ):
        tick_label.set_fontweight(
            "bold"
        )


# 15. MLD AND SIC LEGEND
legend_handles = [
    Line2D(
        [0],
        [0],
        color="black",
        linewidth=1.25,
        linestyle="-",
        label="MLD (m)",
    ),
    Line2D(
        [0],
        [0],
        color="0.18",
        linewidth=1.25,
        linestyle="--",
        label="SIC 15%",
    ),
]

legend = fig.legend(
    handles=legend_handles,
    loc="upper left",
    bbox_to_anchor=(
        0.895,
        0.074,
    ),
    frameon=False,
    fontsize=FONT_LEGEND,
    handlelength=1.8,
    handletextpad=0.5,
    labelspacing=0.45,
    borderaxespad=0.0,
)

for legend_text in (
    legend.get_texts()
):
    legend_text.set_fontweight(
        "bold"
    )


# 16. SAVE THE FIGURE
# Save vector output first.
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

# Compact physical dimensions make the requested 1080-dpi raster
# considerably safer than the previous large figure.
print(
    "Saving 1080-dpi PNG..."
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


# 17. FINAL SUMMARY
summary_lines = [
    "=" * 94,
    "FIGURE 9 — EN4 FRESHWATER STORAGE AND UPPER-OCEAN STABILITY",
    "PLOT-ONLY FINAL WORKFLOW",
    "=" * 94,
    "",
    f"Seasonal input: {SEASONAL_FILE}",
    "",
    "Variables:",
    f"FWC            : {fwc_variable}",
    f"SSS anomaly    : {sss_variable}",
    f"Stratification : {stratification_variable}",
    f"MLD            : {mld_variable}",
    f"SIC            : {sic_variable}",
    "",
    "No monthly calculations were repeated.",
    "",
    "Figure design:",
    "• 4 seasonal rows × 3 variable columns.",
    "• FWC seasonal climatology.",
    "• SSS seasonal anomaly.",
    "• Stratification seasonal climatology.",
    "• Solid MLD contours.",
    "• Dashed SIC 15% contour.",
    "• Wider map-column separation prevents longitude-label overlap.",
    "• Compact contour labels reduce clutter.",
    "",
    "Outputs:",
    f"1080-dpi PNG: {OUT_PNG}",
    f"Vector PDF  : {OUT_PDF}",
    "",
    "Color limits:",
    (
        f"FWC: {FWC_VMIN:.6g} to "
        f"{FWC_VMAX:.6g} m"
    ),
    (
        f"SSS anomaly: ±"
        f"{SSS_LIMIT:.6g}"
    ),
    (
        f"Stratification: "
        f"{STRATIFICATION_VMIN:.6g} to "
        f"{STRATIFICATION_VMAX:.6g} kg m-3"
    ),
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
    "\nFIGURE 9 COMPLETED SUCCESSFULLY."
)
