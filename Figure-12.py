"""
FIGURE 6 — FULL ORGANIZED FINAL WORKFLOW

Title
-----
Quantitative integration of T-S-derived hydrographic properties and
steric-height variability across 13 Antarctic shelf seas.

Main period

2008-2025

Purpose

This standalone script contains the complete Figure 6 workflow extracted from
`so_sealevel_paper_fig.py`, organized from first to last:

    1. Mount Google Drive
    2. Install required packages
    3. Import libraries
    4. Define all input/output paths
    5. Define analysis settings
    6. Define Antarctic shelf-sea sectors and seasons
    7. Read and quality-control Argo profile data
    8. Calculate sea-season hydrographic metrics
    9. Read thermosteric, halosteric and freshwater NetCDF products
   10. Calculate seasonal sea-wise gridded metrics
   11. Merge all 13 seas x 4 seasons = 52 records
   12. Save analysis-ready CSV and XLSX outputs
   13. Calculate regression/correlation statistics
   14. Save regression statistics
   15. Render ONLY the latest corrected Figure 6 plotting version
   16. Save final 1080-dpi PNG and vector PDF

LATEST PLOT VERSION USED
The source file contains repeated Figure 6 plotting revisions. This organized
script intentionally removes the older repeated versions and uses only the last
Figure 6 plotting revision in the source:

    Figure06_TS_steric_13seas_clean_trajectories_final_1080dpi.png
    Figure06_TS_steric_13seas_clean_trajectories_final.pdf

The latest plotting version:
    - reads the computed Figure 6 metrics/regression tables,
    - checks all 13 seas x 4 seasons,
    - uses clean seasonal trajectories,
    - contains no overall figure title,
    - preserves the final panel/legend/profile-count layout.

INPUT FILES

1. /content/drive/MyDrive/SAM_Thesis/Data/
   argo_SO_profiles_2001_2025_cleaned_gridded.nc

2. /content/drive/MyDrive/SAM_Thesis/Processed/
   thermosteric_height_monthly_2008_2025_SO.nc

3. /content/drive/MyDrive/SAM_Thesis/Processed/
   halosteric_height_monthly_2008_2025_SO.nc

4. /content/drive/MyDrive/SAM_Thesis/Processed/
   freshwater_content_monthly_2008_2025_SO.nc

ANALYSIS OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/
    Fig06_TS_steric_metrics_13seas_seasonal.csv
    Fig06_TS_steric_metrics_13seas_seasonal.xlsx
    Fig06_regression_statistics.csv

FINAL FIGURE OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/
    Figure06_TS_steric_13seas_clean_trajectories_final_1080dpi.png
    Figure06_TS_steric_13seas_clean_trajectories_final.pdf

SCIENTIFIC PRESERVATION

The original Figure 6 scientific settings, water-mass criteria, sea sectors,
season definitions, profile-support thresholds, seasonal calculations,
regression methods, and saved analysis tables are preserved.

CORRECTION APPLIED

The supplied source contained an accidentally commented `% 360.0` operation in
longitude normalization. It is restored here so longitude conversion correctly
uses:

    ((longitude + 180.0) % 360.0) - 180.0


"""


# 1. MOUNT GOOGLE DRIVE

try:
    from google.colab import drive

    drive.mount("/content/drive")

except ImportError:
    print(
        "Google Drive mounting skipped because this environment "
        "is not Google Colab."
    )


# 2. INSTALL REQUIRED PACKAGES
# Installs only packages missing from the current Python environment.
# This is valid .py syntax and does not use Colab-only `!pip` commands.

import importlib.util
import subprocess
import sys


REQUIRED_PACKAGES = {
    "netCDF4": "netCDF4",
    "cftime": "cftime",
    "xarray": "xarray",
    "scipy": "scipy",
    "openpyxl": "openpyxl",
    "matplotlib": "matplotlib",
    "pandas": "pandas",
    "numpy": "numpy",
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

from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from netCDF4 import Dataset
from scipy import stats


warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
)


# 4. INPUT / OUTPUT FILES

ARGO_NC = (
    "/content/drive/MyDrive/SAM_Thesis/Data/"
    "argo_SO_profiles_2001_2025_cleaned_gridded.nc"
)

THERMOSTERIC_NC = (
    "/content/drive/MyDrive/SAM_Thesis/Processed/"
    "thermosteric_height_monthly_2008_2025_SO.nc"
)

HALOSTERIC_NC = (
    "/content/drive/MyDrive/SAM_Thesis/Processed/"
    "halosteric_height_monthly_2008_2025_SO.nc"
)

FRESHWATER_NC = (
    "/content/drive/MyDrive/SAM_Thesis/Processed/"
    "freshwater_content_monthly_2008_2025_SO.nc"
)


OUTPUT_DIRECTORY = (
    "/content/drive/MyDrive/SAM_Thesis/paper2"
)

os.makedirs(
    OUTPUT_DIRECTORY,
    exist_ok=True,
)


OUTPUT_METRICS_CSV = os.path.join(
    OUTPUT_DIRECTORY,
    "Fig06_TS_steric_metrics_13seas_seasonal.csv",
)

OUTPUT_METRICS_XLSX = os.path.join(
    OUTPUT_DIRECTORY,
    "Fig06_TS_steric_metrics_13seas_seasonal.xlsx",
)

OUTPUT_REGRESSION_CSV = os.path.join(
    OUTPUT_DIRECTORY,
    "Fig06_regression_statistics.csv",
)


# Final/latest plotting outputs.
OUTPUT_PNG = os.path.join(
    OUTPUT_DIRECTORY,
    "Figure06_TS_steric_13seas_clean_trajectories_final_1080dpi.png",
)

OUTPUT_PDF = os.path.join(
    OUTPUT_DIRECTORY,
    "Figure06_TS_steric_13seas_clean_trajectories_final.pdf",
)


# 3. ANALYSIS SETTINGS

YEAR_START = 2008
YEAR_END = 2025

LAT_MIN = -90.0
LAT_MAX = -60.0

# Argo quality-control ranges
MIN_TEMPERATURE = -3.0
MAX_TEMPERATURE = 20.0

MIN_SALINITY = 0.0
MAX_SALINITY = 42.0

# A profile is considered usable when it has at least this many
# valid paired temperature-salinity observations between 0–2000 dbar.
MIN_VALID_TS_LEVELS = 5

# Minimum profile support requested for each sea-season estimate.
MIN_PROFILES_PER_SEA_SEASON = 10

# Minimum number of valid T-S levels required for a profile-level
# diagnostic in each search interval.
MIN_LEVELS_WW_SEARCH = 3
MIN_LEVELS_SUBSURFACE_SEARCH = 5

# Gridded seasonal values are retained only when at least this many
# monthly sector means are available.
MIN_VALID_GRIDDED_MONTHS = 3

# Depth/pressure ranges
PROFILE_PRESSURE_MIN = 0.0
PROFILE_PRESSURE_MAX = 2000.0

WW_SEARCH_MIN = 50.0
WW_SEARCH_MAX = 300.0

SUBSURFACE_MIN = 200.0
SUBSURFACE_MAX = 1000.0

# Save resolution requested by the user
SAVE_DPI = 1080

# A moderate physical figure size prevents excessive memory use
# during the 1080-dpi PNG render.
FIGURE_SIZE_INCHES = (
    12.5,
    9.5
)


# 4. SEAS AND NON-OVERLAPPING LONGITUDE SECTORS

# Each interval is interpreted as [minimum, maximum), so profiles
# located exactly on a common boundary are not counted twice.
SEA_INFORMATION = [
    (
        "WED",
        "Weddell Sea",
        -60.0,
        -20.0
    ),
    (
        "KHV",
        "King Haakon VII Sea",
        -20.0,
        0.0
    ),
    (
        "RLS",
        "Riiser-Larsen Sea",
        0.0,
        10.0
    ),
    (
        "LAZ",
        "Lazarev Sea",
        10.0,
        30.0
    ),
    (
        "COS",
        "Cosmonauts Sea",
        30.0,
        50.0
    ),
    (
        "COO",
        "Cooperation Sea",
        50.0,
        70.0
    ),
    (
        "DAV",
        "Davis Sea",
        70.0,
        90.0
    ),
    (
        "MAW",
        "Mawson Sea",
        90.0,
        130.0
    ),
    (
        "DUR",
        "D'Urville Sea",
        130.0,
        150.0
    ),
    (
        "SOM",
        "Somov Sea",
        150.0,
        170.0
    ),
    (
        "ROS",
        "Ross Sea",
        170.0,
        -130.0
    ),
    (
        "AMU",
        "Amundsen Sea",
        -130.0,
        -100.0
    ),
    (
        "BEL",
        "Bellingshausen Sea",
        -100.0,
        -60.0
    ),
]


SEA_CODES = [
    item[0]
    for item in SEA_INFORMATION
]

SEA_NAMES = {
    item[0]: item[1]
    for item in SEA_INFORMATION
}

SEA_RANGES = {
    item[0]: (
        item[2],
        item[3]
    )
    for item in SEA_INFORMATION
}


# 5. AUSTRAL SEASONS AND COLOURS

SEASON_ORDER = [
    "Spring",
    "Summer",
    "Autumn",
    "Winter",
]

SEASON_CODE = {
    "Spring": "SON",
    "Summer": "DJF",
    "Autumn": "MAM",
    "Winter": "JJA",
}

SEASON_MONTHS = {
    "Spring": [
        9,
        10,
        11
    ],
    "Summer": [
        12,
        1,
        2
    ],
    "Autumn": [
        3,
        4,
        5
    ],
    "Winter": [
        6,
        7,
        8
    ],
}

# Same colours are used in all four panels.
SEASON_COLOUR = {
    "Spring": "#D62728",
    "Summer": "#1F77B4",
    "Autumn": "#FF7F0E",
    "Winter": "#2CA02C",
}

# Small vertical offsets separate the four seasonal dots in panels a-b.
SEASON_Y_OFFSET = {
    "Spring": -0.24,
    "Summer": -0.08,
    "Autumn": 0.08,
    "Winter": 0.24,
}


# 6. GENERAL HELPER FUNCTIONS

def check_input_files(
    file_paths
):
    """
    Stop immediately and list any missing input files.
    """

    missing_files = [
        path
        for path in file_paths
        if not os.path.exists(
            path
        )
    ]

    if missing_files:

        missing_text = "\n".join(
            missing_files
        )

        raise FileNotFoundError(
            "The following required files were not found:\n"
            f"{missing_text}"
        )


def normalize_longitude(
    longitude
):
    """
    Convert valid longitude values to the interval [-180, 180).
    """

    longitude = np.asarray(
        longitude,
        dtype=np.float64
    )

    return (
        (
            longitude
            + 180.0
        )
         % 360.0
    ) - 180.0


def longitude_sector_mask(
    longitude,
    minimum_longitude,
    maximum_longitude
):
    """
    Create a non-overlapping longitude mask.

    Normal sector:
        minimum <= longitude < maximum

    Dateline-crossing sector:
        longitude >= minimum OR longitude < maximum
    """

    longitude = np.asarray(
        longitude
    )

    if minimum_longitude < maximum_longitude:

        return (
            (
                longitude
                >= minimum_longitude
            )
            &
            (
                longitude
                < maximum_longitude
            )
        )

    return (
        (
            longitude
            >= minimum_longitude
        )
        |
        (
            longitude
            < maximum_longitude
        )
    )


def month_to_season_name(
    month
):
    """
    Convert a month number to an austral season name.
    """

    for season_name, season_month_list in (
        SEASON_MONTHS.items()
    ):

        if month in season_month_list:

            return season_name

    return None


def safe_float(
    value
):
    """
    Convert a scalar to float, returning NaN when conversion fails.
    """

    try:

        return float(
            value
        )

    except (
        TypeError,
        ValueError
    ):

        return np.nan


def nanmedian_or_nan(
    values
):
    """
    Return the median of finite values or NaN.
    """

    values = np.asarray(
        values,
        dtype=np.float64
    )

    values = values[
        np.isfinite(
            values
        )
    ]

    if values.size == 0:

        return np.nan

    return float(
        np.median(
            values
        )
    )


def percentile_or_nan(
    values,
    percentile
):
    """
    Return a percentile of finite values or NaN.
    """

    values = np.asarray(
        values,
        dtype=np.float64
    )

    values = values[
        np.isfinite(
            values
        )
    ]

    if values.size == 0:

        return np.nan

    return float(
        np.percentile(
            values,
            percentile
        )
    )


# 7. EXACT WATER-MASS RULES USED IN FIGURES 4 AND 5

def classify_watermass_existing_criteria(
    salinity,
    temperature,
    pressure
):
    """
    Apply exactly the simplified water-mass hierarchy used in the
    supplied Figure 4 and Figure 5 scripts.

    Note:
    Salinity is retained as an input because this is a T-S diagnostic,
    although the supplied simplified thresholds do not impose a
    numerical salinity boundary.
    """

    if not (
        np.isfinite(
            salinity
        )
        and np.isfinite(
            temperature
        )
        and np.isfinite(
            pressure
        )
    ):

        return None

    # Antarctic Bottom Water
    if (
        pressure
        >= 1500.0
        and temperature
        <= 0.5
    ):

        return "AABW"

    # Surface waters
    if pressure < 200.0:

        if temperature <= -0.5:

            return "WW"

        return "AASW"

    # Subsurface and deep waters
    if pressure >= 200.0:

        if temperature >= 1.0:

            return "CDW"

        return "mCDW"

    return "mCDW"


def existing_mcdw_cdw_mask(
    temperature,
    pressure,
    valid_ts_mask
):
    """
    Vectorized equivalent of the existing classification for the
    combined mCDW/CDW categories.

    The AABW condition has priority and is therefore excluded.
    """

    pressure_2d = np.broadcast_to(
        pressure[
            None,
            :
        ],
        temperature.shape
    )

    aabw_mask = (
        valid_ts_mask
        &
        (
            pressure_2d
            >= 1500.0
        )
        &
        (
            temperature
            <= 0.5
        )
    )

    return (
        valid_ts_mask
        &
        (
            pressure_2d
            >= 200.0
        )
        &
        ~aabw_mask
    )


# 8. READ AND DECODE ARGO PROFILE DATA

def decode_argo_juld(
    julian_days
):
    """
    Decode Argo JULD.

    The source scripts use the standard Argo origin:
        1950-01-01 00:00:00
    """

    julian_days = np.asarray(
        julian_days,
        dtype=np.float64
    )

    output_time = np.full(
        julian_days.shape,
        np.datetime64(
            "NaT"
        ),
        dtype="datetime64[ns]"
    )

    valid_juld = (
        np.isfinite(
            julian_days
        )
        &
        (
            julian_days
            > 0.0
        )
        &
        (
            julian_days
            < 50000.0
        )
    )

    if np.any(
        valid_juld
    ):

        base_date = datetime(
            1950,
            1,
            1
        )

        output_time[
            valid_juld
        ] = np.array(
            [
                np.datetime64(
                    base_date
                    + timedelta(
                        days=float(
                            day
                        )
                    )
                )
                for day in julian_days[
                    valid_juld
                ]
            ],
            dtype="datetime64[ns]"
        )

    return output_time


def read_argo_profiles():
    """
    Read the Argo profile file and apply explicit physical-range masks.
    """

    print(
        "\nReading Argo profile data..."
    )

    with Dataset(
        ARGO_NC,
        mode="r"
    ) as dataset:

        pressure = np.ma.filled(
            dataset.variables[
                "PRES_GRID"
            ][
                :
            ],
            np.nan
        ).astype(
            np.float64
        )

        juld = np.ma.filled(
            dataset.variables[
                "JULD"
            ][
                :
            ],
            np.nan
        ).astype(
            np.float64
        )

        latitude = np.ma.filled(
            dataset.variables[
                "LATITUDE"
            ][
                :
            ],
            np.nan
        ).astype(
            np.float64
        )

        longitude_raw = np.ma.filled(
            dataset.variables[
                "LONGITUDE"
            ][
                :
            ],
            np.nan
        ).astype(
            np.float64
        )

        temperature = np.ma.filled(
            dataset.variables[
                "TEMP"
            ][
                :
            ],
            np.nan
        ).astype(
            np.float32
        )

        salinity = np.ma.filled(
            dataset.variables[
                "PSAL"
            ][
                :
            ],
            np.nan
        ).astype(
            np.float32
        )

    # Remove invalid pressure values.
    pressure[
        ~np.isfinite(
            pressure
        )
        |
        (
            pressure
            < PROFILE_PRESSURE_MIN
        )
        |
        (
            pressure
            > PROFILE_PRESSURE_MAX
        )
    ] = np.nan

    # Remove invalid positions before longitude normalization.
    latitude[
        ~np.isfinite(
            latitude
        )
        |
        (
            latitude
            < -90.0
        )
        |
        (
            latitude
            > 90.0
        )
    ] = np.nan

    longitude_valid = (
        np.isfinite(
            longitude_raw
        )
        &
        (
            longitude_raw
            >= -360.0
        )
        &
        (
            longitude_raw
            <= 360.0
        )
    )

    longitude = np.full(
        longitude_raw.shape,
        np.nan,
        dtype=np.float64
    )

    longitude[
        longitude_valid
    ] = normalize_longitude(
        longitude_raw[
            longitude_valid
        ]
    )

    # Explicitly remove 99999-style values and other physical outliers.
    temperature[
        ~np.isfinite(
            temperature
        )
        |
        (
            temperature
            <= MIN_TEMPERATURE
        )
        |
        (
            temperature
            >= MAX_TEMPERATURE
        )
    ] = np.nan

    salinity[
        ~np.isfinite(
            salinity
        )
        |
        (
            salinity
            <= MIN_SALINITY
        )
        |
        (
            salinity
            >= MAX_SALINITY
        )
    ] = np.nan

    time = decode_argo_juld(
        juld
    )

    time_index = pd.DatetimeIndex(
        pd.to_datetime(
            time
        )
    )

    year = time_index.year.to_numpy()
    month = time_index.month.to_numpy()

    season = np.array(
        [
            (
                month_to_season_name(
                    int(
                        month_value
                    )
                )
                if pd.notna(
                    month_value
                )
                else None
            )
            for month_value in month
        ],
        dtype=object
    )

    print(
        "Argo arrays loaded:",
        f"{temperature.shape[0]:,} profiles × "
        f"{temperature.shape[1]:,} pressure levels"
    )

    return {
        "pressure": pressure,
        "latitude": latitude,
        "longitude": longitude,
        "time": time,
        "year": year,
        "month": month,
        "season": season,
        "temperature": temperature,
        "salinity": salinity,
    }


# 9. CALCULATE PROFILE-BASED SEA-SEASON METRICS

def calculate_profile_metrics(
    argo_data
):
    """
    Calculate:
      - valid profile count
      - profile occurrence of mCDW/CDW
      - median WW core temperature and pressure
      - median profile-wise maximum temperature from 200–1000 dbar
    """

    pressure = argo_data[
        "pressure"
    ]

    latitude = argo_data[
        "latitude"
    ]

    longitude = argo_data[
        "longitude"
    ]

    year = argo_data[
        "year"
    ]

    season_array = argo_data[
        "season"
    ]

    temperature = argo_data[
        "temperature"
    ]

    salinity = argo_data[
        "salinity"
    ]

    pressure_in_profile_range = (
        np.isfinite(
            pressure
        )
        &
        (
            pressure
            >= PROFILE_PRESSURE_MIN
        )
        &
        (
            pressure
            <= PROFILE_PRESSURE_MAX
        )
    )

    ww_pressure_indices = np.where(
        np.isfinite(
            pressure
        )
        &
        (
            pressure
            >= WW_SEARCH_MIN
        )
        &
        (
            pressure
            <= WW_SEARCH_MAX
        )
    )[0]

    subsurface_pressure_indices = np.where(
        np.isfinite(
            pressure
        )
        &
        (
            pressure
            >= SUBSURFACE_MIN
        )
        &
        (
            pressure
            <= SUBSURFACE_MAX
        )
    )[0]

    if ww_pressure_indices.size == 0:

        raise ValueError(
            "No pressure levels were found in the WW search range."
        )

    if subsurface_pressure_indices.size == 0:

        raise ValueError(
            "No pressure levels were found in the 200–1000 dbar range."
        )

    base_profile_mask = (
        np.isfinite(
            latitude
        )
        &
        np.isfinite(
            longitude
        )
        &
        (
            latitude
            >= LAT_MIN
        )
        &
        (
            latitude
            <= LAT_MAX
        )
        &
        (
            year
            >= YEAR_START
        )
        &
        (
            year
            <= YEAR_END
        )
        &
        pd.notna(
            season_array
        )
    )

    print(
        "Profiles within 2008–2025 and 60–90°S:",
        f"{int(base_profile_mask.sum()):,}"
    )

    output_rows = []

    for sea_code in SEA_CODES:

        longitude_minimum, longitude_maximum = (
            SEA_RANGES[
                sea_code
            ]
        )

        sea_mask = longitude_sector_mask(
            longitude,
            longitude_minimum,
            longitude_maximum
        )

        for season_name in SEASON_ORDER:

            profile_mask = (
                base_profile_mask
                &
                sea_mask
                &
                (
                    season_array
                    == season_name
                )
            )

            profile_indices = np.where(
                profile_mask
            )[0]

            selected_profile_count = int(
                profile_indices.size
            )

            row = {
                "sea_abbreviation": sea_code,
                "sea_name": SEA_NAMES[
                    sea_code
                ],
                "season": season_name,
                "season_code": SEASON_CODE[
                    season_name
                ],
                "selected_profile_count": selected_profile_count,
                "valid_profile_count": 0,
                "mCDW_CDW_profile_count": 0,
                "mCDW_CDW_occurrence_raw_percent": np.nan,
                "mCDW_CDW_occurrence_percent": np.nan,
                "WW_core_profile_count": 0,
                "WW_core_temperature_median_degC": np.nan,
                "WW_core_temperature_p25_degC": np.nan,
                "WW_core_temperature_p75_degC": np.nan,
                "WW_core_pressure_median_dbar": np.nan,
                "WW_core_pressure_p25_dbar": np.nan,
                "WW_core_pressure_p75_dbar": np.nan,
                "subsurface_maximum_profile_count": 0,
                (
                    "maximum_subsurface_temperature_"
                    "200_1000dbar_median_degC"
                ): np.nan,
                (
                    "maximum_subsurface_temperature_"
                    "200_1000dbar_p25_degC"
                ): np.nan,
                (
                    "maximum_subsurface_temperature_"
                    "200_1000dbar_p75_degC"
                ): np.nan,
                "profile_support_status": "insufficient",
            }

            if selected_profile_count == 0:

                output_rows.append(
                    row
                )

                continue

            temperature_subset = temperature[
                profile_indices,
                :
            ].astype(
                np.float64,
                copy=False
            )

            salinity_subset = salinity[
                profile_indices,
                :
            ].astype(
                np.float64,
                copy=False
            )

            valid_ts = (
                np.isfinite(
                    temperature_subset
                )
                &
                np.isfinite(
                    salinity_subset
                )
                &
                pressure_in_profile_range[
                    None,
                    :
                ]
            )

            valid_level_count = valid_ts.sum(
                axis=1
            )

            valid_profile_mask = (
                valid_level_count
                >= MIN_VALID_TS_LEVELS
            )

            valid_profile_count = int(
                valid_profile_mask.sum()
            )

            row[
                "valid_profile_count"
            ] = valid_profile_count

            # mCDW/CDW profile occurrence

            mcdw_cdw_samples = existing_mcdw_cdw_mask(
                temperature_subset,
                pressure,
                valid_ts
            )

            profile_contains_mcdw_cdw = (
                valid_profile_mask
                &
                np.any(
                    mcdw_cdw_samples,
                    axis=1
                )
            )

            mcdw_cdw_profile_count = int(
                profile_contains_mcdw_cdw.sum()
            )

            row[
                "mCDW_CDW_profile_count"
            ] = mcdw_cdw_profile_count

            if valid_profile_count > 0:

                raw_occurrence = (
                    100.0
                    * mcdw_cdw_profile_count
                    / valid_profile_count
                )

                row[
                    "mCDW_CDW_occurrence_raw_percent"
                ] = float(
                    raw_occurrence
                )

                if (
                    valid_profile_count
                    >= MIN_PROFILES_PER_SEA_SEASON
                ):

                    row[
                        "mCDW_CDW_occurrence_percent"
                    ] = float(
                        raw_occurrence
                    )

            # Winter Water core

            ww_temperature = temperature_subset[
                :,
                ww_pressure_indices
            ]

            ww_salinity = salinity_subset[
                :,
                ww_pressure_indices
            ]

            ww_pressure = pressure[
                ww_pressure_indices
            ]

            ww_valid = (
                np.isfinite(
                    ww_temperature
                )
                &
                np.isfinite(
                    ww_salinity
                )
                &
                valid_profile_mask[
                    :,
                    None
                ]
            )

            ww_valid_level_count = ww_valid.sum(
                axis=1
            )

            ww_search_profile_mask = (
                ww_valid_level_count
                >= MIN_LEVELS_WW_SEARCH
            )

            ww_core_temperatures = []
            ww_core_pressures = []

            for local_profile_index in np.where(
                ww_search_profile_mask
            )[0]:

                profile_temperature = np.where(
                    ww_valid[
                        local_profile_index
                    ],
                    ww_temperature[
                        local_profile_index
                    ],
                    np.nan
                )

                minimum_local_index = int(
                    np.nanargmin(
                        profile_temperature
                    )
                )

                core_temperature = float(
                    ww_temperature[
                        local_profile_index,
                        minimum_local_index
                    ]
                )

                core_salinity = float(
                    ww_salinity[
                        local_profile_index,
                        minimum_local_index
                    ]
                )

                core_pressure = float(
                    ww_pressure[
                        minimum_local_index
                    ]
                )

                core_classification = (
                    classify_watermass_existing_criteria(
                        core_salinity,
                        core_temperature,
                        core_pressure
                    )
                )

                if core_classification == "WW":

                    ww_core_temperatures.append(
                        core_temperature
                    )

                    ww_core_pressures.append(
                        core_pressure
                    )

            ww_core_profile_count = len(
                ww_core_temperatures
            )

            row[
                "WW_core_profile_count"
            ] = ww_core_profile_count

            if (
                ww_core_profile_count
                >= MIN_PROFILES_PER_SEA_SEASON
            ):

                row[
                    "WW_core_temperature_median_degC"
                ] = nanmedian_or_nan(
                    ww_core_temperatures
                )

                row[
                    "WW_core_temperature_p25_degC"
                ] = percentile_or_nan(
                    ww_core_temperatures,
                    25
                )

                row[
                    "WW_core_temperature_p75_degC"
                ] = percentile_or_nan(
                    ww_core_temperatures,
                    75
                )

                row[
                    "WW_core_pressure_median_dbar"
                ] = nanmedian_or_nan(
                    ww_core_pressures
                )

                row[
                    "WW_core_pressure_p25_dbar"
                ] = percentile_or_nan(
                    ww_core_pressures,
                    25
                )

                row[
                    "WW_core_pressure_p75_dbar"
                ] = percentile_or_nan(
                    ww_core_pressures,
                    75
                )

            # Profile-wise maximum subsurface temperature

            subsurface_temperature = temperature_subset[
                :,
                subsurface_pressure_indices
            ]

            subsurface_salinity = salinity_subset[
                :,
                subsurface_pressure_indices
            ]

            subsurface_valid = (
                np.isfinite(
                    subsurface_temperature
                )
                &
                np.isfinite(
                    subsurface_salinity
                )
                &
                valid_profile_mask[
                    :,
                    None
                ]
            )

            subsurface_level_count = subsurface_valid.sum(
                axis=1
            )

            subsurface_profile_mask = (
                subsurface_level_count
                >= MIN_LEVELS_SUBSURFACE_SEARCH
            )

            profile_maximum_temperature = np.full(
                selected_profile_count,
                np.nan,
                dtype=np.float64
            )

            if np.any(
                subsurface_profile_mask
            ):

                valid_subsurface_temperature = np.where(
                    subsurface_valid[
                        subsurface_profile_mask
                    ],
                    subsurface_temperature[
                        subsurface_profile_mask
                    ],
                    np.nan
                )

                profile_maximum_temperature[
                    subsurface_profile_mask
                ] = np.nanmax(
                    valid_subsurface_temperature,
                    axis=1
                )

            finite_profile_maximum = (
                profile_maximum_temperature[
                    np.isfinite(
                        profile_maximum_temperature
                    )
                ]
            )

            subsurface_profile_count = int(
                finite_profile_maximum.size
            )

            row[
                "subsurface_maximum_profile_count"
            ] = subsurface_profile_count

            if (
                subsurface_profile_count
                >= MIN_PROFILES_PER_SEA_SEASON
            ):

                row[
                    (
                        "maximum_subsurface_temperature_"
                        "200_1000dbar_median_degC"
                    )
                ] = nanmedian_or_nan(
                    finite_profile_maximum
                )

                row[
                    (
                        "maximum_subsurface_temperature_"
                        "200_1000dbar_p25_degC"
                    )
                ] = percentile_or_nan(
                    finite_profile_maximum,
                    25
                )

                row[
                    (
                        "maximum_subsurface_temperature_"
                        "200_1000dbar_p75_degC"
                    )
                ] = percentile_or_nan(
                    finite_profile_maximum,
                    75
                )

            if (
                valid_profile_count
                >= MIN_PROFILES_PER_SEA_SEASON
            ):

                row[
                    "profile_support_status"
                ] = "sufficient"

            output_rows.append(
                row
            )

            del (
                temperature_subset,
                salinity_subset,
                valid_ts,
                mcdw_cdw_samples,
                ww_temperature,
                ww_salinity,
                subsurface_temperature,
                subsurface_salinity
            )

    output_dataframe = pd.DataFrame(
        output_rows
    )

    return output_dataframe


# 10. GRIDDED PRODUCT HELPERS
def mask_fill_values(
    data_array
):
    """
    Convert encoded fill values to NaN without assuming that xarray
    decoded every source file identically.
    """

    cleaned = data_array.astype(
        np.float64
    )

    fill_candidates = [
        data_array.attrs.get(
            "_FillValue"
        ),
        data_array.attrs.get(
            "missing_value"
        ),
        data_array.encoding.get(
            "_FillValue"
        ),
    ]

    for fill_value in fill_candidates:

        if fill_value is None:

            continue

        try:

            cleaned = cleaned.where(
                ~np.isclose(
                    cleaned,
                    float(
                        fill_value
                    )
                )
            )

        except (
            TypeError,
            ValueError
        ):

            pass

    cleaned = cleaned.where(
        np.isfinite(
            cleaned
        )
    )

    # The source inventory identifies -9999 as the fill value.
    cleaned = cleaned.where(
        cleaned
        > -9000.0
    )

    return cleaned


def prepare_gridded_dataset(
    netcdf_path,
    variable_name
):
    """
    Open one monthly gridded product, standardize longitude and mask
    fill values.
    """

    dataset = xr.open_dataset(
        netcdf_path,
        decode_times=True,
        mask_and_scale=True
    )

    required_names = {
        "time",
        "lat",
        "lon",
        variable_name,
    }

    missing_names = [
        name
        for name in required_names
        if name not in dataset
        and name not in dataset.coords
    ]

    if missing_names:

        dataset.close()

        raise KeyError(
            f"{netcdf_path} is missing: {missing_names}"
        )

    data_array = mask_fill_values(
        dataset[
            variable_name
        ]
    )

    longitude_normalized = normalize_longitude(
        dataset[
            "lon"
        ].values
    )

    data_array = data_array.assign_coords(
        lon=(
            "lon",
            longitude_normalized
        )
    ).sortby(
        "lon"
    )

    if "nprof" in dataset:

        nprof = dataset[
            "nprof"
        ].astype(
            np.float64
        )

        nprof = nprof.assign_coords(
            lon=(
                "lon",
                longitude_normalized
            )
        ).sortby(
            "lon"
        )

        data_array = data_array.where(
            nprof
            > 0.0
        )

    else:

        nprof = None

    data_array = data_array.sel(
        time=slice(
            f"{YEAR_START}-01-01",
            f"{YEAR_END}-12-31"
        )
    )

    if nprof is not None:

        nprof = nprof.sel(
            time=data_array[
                "time"
            ]
        )

    return (
        dataset,
        data_array,
        nprof
    )


def sector_monthly_area_weighted_mean(
    data_array,
    nprof,
    sea_code,
    season_name
):
    """
    Calculate monthly sector means using cosine-latitude weighting,
    then return the mean of the available monthly sector means.

    The nprof field is used as a validity mask and support diagnostic,
    not as a spatial weight. This avoids biasing the regional mean
    toward the most heavily sampled cells.
    """

    longitude_minimum, longitude_maximum = (
        SEA_RANGES[
            sea_code
        ]
    )

    longitude_values = data_array[
        "lon"
    ].values

    longitude_mask_values = longitude_sector_mask(
        longitude_values,
        longitude_minimum,
        longitude_maximum
    )

    longitude_mask = xr.DataArray(
        longitude_mask_values,
        dims=[
            "lon"
        ],
        coords={
            "lon": data_array[
                "lon"
            ]
        }
    )

    month_mask_values = np.isin(
        data_array[
            "time"
        ].dt.month.values,
        SEASON_MONTHS[
            season_name
        ]
    )

    month_mask = xr.DataArray(
        month_mask_values,
        dims=[
            "time"
        ],
        coords={
            "time": data_array[
                "time"
            ]
        }
    )

    selected_data = data_array.where(
        longitude_mask,
        drop=True
    ).where(
        month_mask,
        drop=True
    )

    if selected_data.sizes.get(
        "time",
        0
    ) == 0:

        return {
            "seasonal_mean": np.nan,
            "seasonal_standard_deviation": np.nan,
            "valid_month_count": 0,
            "valid_gridcell_month_count": 0,
            "product_profile_count": 0,
        }

    latitude_weight = xr.DataArray(
        np.cos(
            np.deg2rad(
                selected_data[
                    "lat"
                ].values
            )
        ),
        dims=[
            "lat"
        ],
        coords={
            "lat": selected_data[
                "lat"
            ]
        }
    )

    broadcast_weight = latitude_weight.broadcast_like(
        selected_data
    )

    valid_data = selected_data.notnull()

    weighted_numerator = (
        selected_data
        * broadcast_weight
    ).sum(
        dim=[
            "lat",
            "lon"
        ],
        skipna=True
    )

    weighted_denominator = broadcast_weight.where(
        valid_data
    ).sum(
        dim=[
            "lat",
            "lon"
        ],
        skipna=True
    )

    monthly_sector_mean = (
        weighted_numerator
        / weighted_denominator
    ).where(
        weighted_denominator
        > 0.0
    )

    monthly_values = np.asarray(
        monthly_sector_mean.values,
        dtype=np.float64
    )

    monthly_values = monthly_values[
        np.isfinite(
            monthly_values
        )
    ]

    valid_month_count = int(
        monthly_values.size
    )

    valid_gridcell_month_count = int(
        valid_data.sum().values
    )

    if nprof is not None:

        selected_nprof = nprof.where(
            longitude_mask,
            drop=True
        ).where(
            month_mask,
            drop=True
        )

        product_profile_count = int(
            np.nansum(
                selected_nprof.values
            )
        )

    else:

        product_profile_count = 0

    if (
        valid_month_count
        < MIN_VALID_GRIDDED_MONTHS
    ):

        seasonal_mean = np.nan
        seasonal_standard_deviation = np.nan

    else:

        seasonal_mean = float(
            np.mean(
                monthly_values
            )
        )

        seasonal_standard_deviation = float(
            np.std(
                monthly_values,
                ddof=1
            )
        ) if valid_month_count > 1 else 0.0

    return {
        "seasonal_mean": seasonal_mean,
        "seasonal_standard_deviation": (
            seasonal_standard_deviation
        ),
        "valid_month_count": valid_month_count,
        "valid_gridcell_month_count": (
            valid_gridcell_month_count
        ),
        "product_profile_count": product_profile_count,
    }


def calculate_gridded_product_metrics(
    netcdf_path,
    variable_name,
    output_prefix,
    conversion_factor=1.0
):
    """
    Calculate one seasonal metric for every sea-season combination.
    """

    print(
        f"\nReading gridded product: {variable_name}"
    )

    (
        dataset,
        data_array,
        nprof
    ) = prepare_gridded_dataset(
        netcdf_path,
        variable_name
    )

    output_rows = []

    for sea_code in SEA_CODES:

        for season_name in SEASON_ORDER:

            summary = sector_monthly_area_weighted_mean(
                data_array,
                nprof,
                sea_code,
                season_name
            )

            seasonal_mean = summary[
                "seasonal_mean"
            ]

            seasonal_standard_deviation = summary[
                "seasonal_standard_deviation"
            ]

            if np.isfinite(
                seasonal_mean
            ):

                seasonal_mean = (
                    seasonal_mean
                    * conversion_factor
                )

            if np.isfinite(
                seasonal_standard_deviation
            ):

                seasonal_standard_deviation = (
                    seasonal_standard_deviation
                    * conversion_factor
                )

            output_rows.append(
                {
                    "sea_abbreviation": sea_code,
                    "season": season_name,
                    output_prefix: seasonal_mean,
                    (
                        output_prefix
                        + "_monthly_sd"
                    ): seasonal_standard_deviation,
                    (
                        output_prefix
                        + "_valid_month_count"
                    ): summary[
                        "valid_month_count"
                    ],
                    (
                        output_prefix
                        + "_valid_gridcell_month_count"
                    ): summary[
                        "valid_gridcell_month_count"
                    ],
                    (
                        output_prefix
                        + "_product_profile_count"
                    ): summary[
                        "product_profile_count"
                    ],
                }
            )

    dataset.close()

    return pd.DataFrame(
        output_rows
    )


# 11. REGRESSION AND CONFIDENCE-INTERVAL FUNCTIONS

def calculate_regression_statistics(
    dataframe,
    x_column,
    y_column,
    panel_name,
    slope_unit
):
    """
    Calculate Pearson correlation, Spearman correlation and ordinary
    least-squares slope with a two-sided 95% confidence interval.
    """

    valid_data = dataframe[
        [
            x_column,
            y_column
        ]
    ].replace(
        [
            np.inf,
            -np.inf
        ],
        np.nan
    ).dropna()

    x = valid_data[
        x_column
    ].to_numpy(
        dtype=np.float64
    )

    y = valid_data[
        y_column
    ].to_numpy(
        dtype=np.float64
    )

    n_observations = int(
        x.size
    )

    output = {
        "panel": panel_name,
        "x_variable": x_column,
        "y_variable": y_column,
        "n": n_observations,
        "pearson_r": np.nan,
        "pearson_p": np.nan,
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
        "regression_slope": np.nan,
        "slope_unit": slope_unit,
        "slope_ci95_lower": np.nan,
        "slope_ci95_upper": np.nan,
        "intercept": np.nan,
        "regression_p": np.nan,
        "r_squared": np.nan,
    }

    if (
        n_observations
        < 3
        or np.nanstd(
            x
        )
        == 0.0
        or np.nanstd(
            y
        )
        == 0.0
    ):

        return output

    pearson_result = stats.pearsonr(
        x,
        y
    )

    spearman_result = stats.spearmanr(
        x,
        y
    )

    regression_result = stats.linregress(
        x,
        y
    )

    degrees_of_freedom = (
        n_observations
        - 2
    )

    t_critical = stats.t.ppf(
        0.975,
        degrees_of_freedom
    )

    slope_margin = (
        t_critical
        * regression_result.stderr
    )

    output.update(
        {
            "pearson_r": float(
                pearson_result.statistic
            ),
            "pearson_p": float(
                pearson_result.pvalue
            ),
            "spearman_rho": float(
                spearman_result.statistic
            ),
            "spearman_p": float(
                spearman_result.pvalue
            ),
            "regression_slope": float(
                regression_result.slope
            ),
            "slope_ci95_lower": float(
                regression_result.slope
                - slope_margin
            ),
            "slope_ci95_upper": float(
                regression_result.slope
                + slope_margin
            ),
            "intercept": float(
                regression_result.intercept
            ),
            "regression_p": float(
                regression_result.pvalue
            ),
            "r_squared": float(
                regression_result.rvalue
                ** 2
            ),
        }
    )

    return output


def regression_line_and_mean_ci(
    x,
    y,
    number_of_line_points=250
):
    """
    Return fitted values and the 95% confidence interval for the mean
    regression response.
    """

    x = np.asarray(
        x,
        dtype=np.float64
    )

    y = np.asarray(
        y,
        dtype=np.float64
    )

    valid = (
        np.isfinite(
            x
        )
        &
        np.isfinite(
            y
        )
    )

    x = x[
        valid
    ]

    y = y[
        valid
    ]

    if (
        x.size
        < 3
        or np.std(
            x
        )
        == 0.0
    ):

        return None

    regression_result = stats.linregress(
        x,
        y
    )

    x_line = np.linspace(
        np.min(
            x
        ),
        np.max(
            x
        ),
        number_of_line_points
    )

    y_line = (
        regression_result.intercept
        + regression_result.slope
        * x_line
    )

    fitted_observations = (
        regression_result.intercept
        + regression_result.slope
        * x
    )

    residuals = (
        y
        - fitted_observations
    )

    degrees_of_freedom = (
        x.size
        - 2
    )

    residual_standard_error = np.sqrt(
        np.sum(
            residuals
            ** 2
        )
        / degrees_of_freedom
    )

    x_mean = np.mean(
        x
    )

    sum_squared_x = np.sum(
        (
            x
            - x_mean
        )
        ** 2
    )

    t_critical = stats.t.ppf(
        0.975,
        degrees_of_freedom
    )

    confidence_half_width = (
        t_critical
        * residual_standard_error
        * np.sqrt(
            (
                1.0
                / x.size
            )
            +
            (
                (
                    x_line
                    - x_mean
                )
                ** 2
                / sum_squared_x
            )
        )
    )

    return {
        "x_line": x_line,
        "y_line": y_line,
        "ci_lower": (
            y_line
            - confidence_half_width
        ),
        "ci_upper": (
            y_line
            + confidence_half_width
        ),
    }


# 13. RUN THE COMPLETE ANALYSIS

check_input_files(
    [
        ARGO_NC,
        THERMOSTERIC_NC,
        HALOSTERIC_NC,
        FRESHWATER_NC,
    ]
)


# 13.1 Argo metrics

argo_data = read_argo_profiles()

profile_metrics = calculate_profile_metrics(
    argo_data
)


# Free the two large profile matrices before the high-resolution render.
del argo_data

gc.collect()


# 13.2 Thermosteric seasonal sector means

thermosteric_metrics = (
    calculate_gridded_product_metrics(
        THERMOSTERIC_NC,
        "thermosteric_height",
        "thermosteric_height_anomaly_cm",
        conversion_factor=100.0
    )
)


# 13.3 Halosteric seasonal sector means

halosteric_metrics = (
    calculate_gridded_product_metrics(
        HALOSTERIC_NC,
        "halosteric_height",
        "halosteric_height_anomaly_cm",
        conversion_factor=100.0
    )
)


# 13.4 Freshwater seasonal sector means

freshwater_metrics = (
    calculate_gridded_product_metrics(
        FRESHWATER_NC,
        "freshwater_content",
        "freshwater_content_0_200dbar_m",
        conversion_factor=1.0
    )
)


# 13.5 Merge all sea-season metrics

metrics = profile_metrics.merge(
    thermosteric_metrics,
    on=[
        "sea_abbreviation",
        "season"
    ],
    how="left"
)

metrics = metrics.merge(
    halosteric_metrics,
    on=[
        "sea_abbreviation",
        "season"
    ],
    how="left"
)

metrics = metrics.merge(
    freshwater_metrics,
    on=[
        "sea_abbreviation",
        "season"
    ],
    how="left"
)


sea_order_map = {
    sea_code: index
    for index, sea_code in enumerate(
        SEA_CODES
    )
}

season_order_map = {
    season_name: index
    for index, season_name in enumerate(
        SEASON_ORDER
    )
}

metrics[
    "_sea_order"
] = metrics[
    "sea_abbreviation"
].map(
    sea_order_map
)

metrics[
    "_season_order"
] = metrics[
    "season"
].map(
    season_order_map
)

metrics = metrics.sort_values(
    [
        "_sea_order",
        "_season_order"
    ]
).drop(
    columns=[
        "_sea_order",
        "_season_order"
    ]
).reset_index(
    drop=True
)


# 13.6 Place the required columns first

required_columns_first = [
    "sea_name",
    "sea_abbreviation",
    "season",
    "season_code",
    "valid_profile_count",
    "mCDW_CDW_occurrence_percent",
    "WW_core_temperature_median_degC",
    "WW_core_pressure_median_dbar",
    (
        "maximum_subsurface_temperature_"
        "200_1000dbar_median_degC"
    ),
    "thermosteric_height_anomaly_cm",
    "freshwater_content_0_200dbar_m",
    "halosteric_height_anomaly_cm",
]

remaining_columns = [
    column
    for column in metrics.columns
    if column not in required_columns_first
]

metrics = metrics[
    required_columns_first
    + remaining_columns
]


# 13.7 Save supporting tables

metrics.to_csv(
    OUTPUT_METRICS_CSV,
    index=False,
    float_format="%.6f"
)

metrics.to_excel(
    OUTPUT_METRICS_XLSX,
    index=False
)


# 14. STATISTICAL ANALYSIS FOR PANELS C AND D

PANEL_C_X = (
    "maximum_subsurface_temperature_"
    "200_1000dbar_median_degC"
)

PANEL_C_Y = (
    "thermosteric_height_anomaly_cm"
)

PANEL_D_X = (
    "freshwater_content_0_200dbar_m"
)

PANEL_D_Y = (
    "halosteric_height_anomaly_cm"
)


panel_c_statistics = calculate_regression_statistics(
    metrics,
    PANEL_C_X,
    PANEL_C_Y,
    panel_name="Figure 6c",
    slope_unit="cm per °C"
)

panel_d_statistics = calculate_regression_statistics(
    metrics,
    PANEL_D_X,
    PANEL_D_Y,
    panel_name="Figure 6d",
    slope_unit="cm per m FWC"
)


regression_statistics = pd.DataFrame(
    [
        panel_c_statistics,
        panel_d_statistics,
    ]
)

regression_statistics.to_csv(
    OUTPUT_REGRESSION_CSV,
    index=False,
    float_format="%.8f"
)


# 15. DATA-SUPPORT DIAGNOSTICS

print(
    "\n"
    + "=" * 78
)

print(
    "FIGURE 6 SUPPORT SUMMARY"
)

print(
    "=" * 78
)

print(
    "\nValid sea-season mCDW/CDW points:",
    int(
        metrics[
            "mCDW_CDW_occurrence_percent"
        ].notna().sum()
    ),
    "/ 52"
)

print(
    "Valid sea-season WW-core points:",
    int(
        metrics[
            "WW_core_temperature_median_degC"
        ].notna().sum()
    ),
    "/ 52"
)

print(
    "Valid Figure 6c pairs:",
    panel_c_statistics[
        "n"
    ]
)

print(
    "Valid Figure 6d pairs:",
    panel_d_statistics[
        "n"
    ]
)


finite_occurrence = metrics[
    "mCDW_CDW_occurrence_percent"
].dropna()

if (
    not finite_occurrence.empty
    and np.nanmin(
        finite_occurrence
    )
    >= 99.0
):

    print(
        "\nIMPORTANT THRESHOLD DIAGNOSTIC:"
    )

    print(
        "All supported mCDW/CDW occurrence values are approximately "
        "100%. This follows directly from the existing Figure 4/5 "
        "rules, because almost every valid non-AABW observation at "
        "pressure >= 200 dbar is classified as mCDW or CDW."
    )

    print(
        "No alternative threshold has been introduced in this code."
    )




# 16. LATEST CORRECTED FIGURE 6 PLOTTING VERSION
# Older repeated Figure 6 plotting blocks from the source file are intentionally
# omitted. The code below is the last/latest Figure 6 plotting revision.

# 2. INPUT AND OUTPUT PATHS
OUTPUT_DIRECTORY = "/content/drive/MyDrive/SAM_Thesis/paper2"

METRICS_CSV = os.path.join(
    OUTPUT_DIRECTORY,
    "Fig06_TS_steric_metrics_13seas_seasonal.csv",
)

REGRESSION_CSV = os.path.join(
    OUTPUT_DIRECTORY,
    "Fig06_regression_statistics.csv",
)

OUTPUT_PNG = os.path.join(
    OUTPUT_DIRECTORY,
    "Figure06_TS_steric_13seas_clean_trajectories_final_1080dpi.png",
)

OUTPUT_PDF = os.path.join(
    OUTPUT_DIRECTORY,
    "Figure06_TS_steric_13seas_clean_trajectories_final.pdf",
)

os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)


# 3. USER-CONTROLLABLE PLOT SETTINGS
SAVE_DPI = 1080

# Compact enough to avoid an unnecessarily huge PNG while retaining
# publication-quality text and panel spacing.
FIGURE_SIZE = (14.8, 10.8)

FONT_PANEL_TITLE = 13.5
FONT_AXIS_LABEL = 12.0
FONT_TICK = 10.0
FONT_LEGEND = 10.0
FONT_SEA_LABEL = 8.0
FONT_STATISTICS = 9.0

DOT_SIZE_TOP = 54
EDGE_WIDTH = 0.75

# Marker-size limits for panels c and d.
MIN_SCATTER_SIZE = 42
MAX_SCATTER_SIZE = 125

# The four seasonal points in panels a and b are vertically offset
# around each sea row so that overlapping values remain visible.
SEASON_VERTICAL_OFFSET = {
    "Spring": -0.18,
    "Summer": -0.06,
    "Autumn":  0.06,
    "Winter":  0.18,
}


# 4. SEA AND SEASON ORDER
SEA_ORDER = [
    "WED",
    "KHV",
    "RLS",
    "LAZ",
    "COS",
    "COO",
    "DAV",
    "MAW",
    "DUR",
    "SOM",
    "ROS",
    "AMU",
    "BEL",
]

SEA_FULL_NAME = {
    "WED": "Weddell Sea",
    "KHV": "King Haakon VII Sea",
    "RLS": "Riiser-Larsen Sea",
    "LAZ": "Lazarev Sea",
    "COS": "Cosmonauts Sea",
    "COO": "Cooperation Sea",
    "DAV": "Davis Sea",
    "MAW": "Mawson Sea",
    "DUR": "D'Urville Sea",
    "SOM": "Somov Sea",
    "ROS": "Ross Sea",
    "AMU": "Amundsen Sea",
    "BEL": "Bellingshausen Sea",
}

SEASON_ORDER = [
    "Spring",
    "Summer",
    "Autumn",
    "Winter",
]

SEASON_COLOR = {
    "Spring": "#d62728",
    "Summer": "#1f77b4",
    "Autumn": "#ff7f0e",
    "Winter": "#2ca02c",
}

# A separate marker shape identifies every sea in panels c and d.
# Colour still identifies season.
SEA_MARKER = {
    "WED": "o",
    "KHV": "s",
    "RLS": "^",
    "LAZ": "v",
    "COS": "D",
    "COO": "P",
    "DAV": "X",
    "MAW": "<",
    "DUR": ">",
    "SOM": "*",
    "ROS": "h",
    "AMU": "p",
    "BEL": "8",
}


# 5. COLUMN NAMES FROM THE EXISTING METRICS TABLE
COL_SEA = "sea_abbreviation"
COL_SEASON = "season"
COL_PROFILE_COUNT = "valid_profile_count"

COL_OCCURRENCE = "mCDW_CDW_occurrence_percent"
COL_WW_TEMPERATURE = "WW_core_temperature_median_degC"

COL_SUBSURFACE_TEMPERATURE = (
    "maximum_subsurface_temperature_200_1000dbar_median_degC"
)
COL_THERMOSTERIC = "thermosteric_height_anomaly_cm"

COL_FRESHWATER = "freshwater_content_0_200dbar_m"
COL_HALOSTERIC = "halosteric_height_anomaly_cm"


# 6. READ THE EXISTING TABLES — NO METRIC RECOMPUTATION
if not os.path.exists(METRICS_CSV):
    raise FileNotFoundError(
        "Metrics CSV was not found:\n"
        f"{METRICS_CSV}\n\n"
        "Place the previously generated CSV in the paper2 folder."
    )

if not os.path.exists(REGRESSION_CSV):
    raise FileNotFoundError(
        "Regression-statistics CSV was not found:\n"
        f"{REGRESSION_CSV}\n\n"
        "Place the previously generated regression CSV in the paper2 folder."
    )

metrics = pd.read_csv(METRICS_CSV)
regression = pd.read_csv(REGRESSION_CSV)


# 7. STANDARDIZE LABELS AND VALIDATE THE 52 COMBINATIONS
metrics[COL_SEA] = (
    metrics[COL_SEA]
    .astype(str)
    .str.strip()
    .replace({
        "KHVII": "KHV",
        "KH VII": "KHV",
    })
)

metrics[COL_SEASON] = (
    metrics[COL_SEASON]
    .astype(str)
    .str.strip()
    .str.title()
)

required_columns = [
    COL_SEA,
    COL_SEASON,
    COL_PROFILE_COUNT,
    COL_OCCURRENCE,
    COL_WW_TEMPERATURE,
    COL_SUBSURFACE_TEMPERATURE,
    COL_THERMOSTERIC,
    COL_FRESHWATER,
    COL_HALOSTERIC,
]

missing_columns = [
    column
    for column in required_columns
    if column not in metrics.columns
]

if missing_columns:
    raise KeyError(
        "The following required columns are missing from the metrics CSV:\n"
        + "\n".join(missing_columns)
    )

# Keep the intended 13 seas and four seasons only.
metrics = metrics[
    metrics[COL_SEA].isin(SEA_ORDER)
    & metrics[COL_SEASON].isin(SEASON_ORDER)
].copy()

# Detect duplicate sea-season records.
duplicate_mask = metrics.duplicated(
    subset=[COL_SEA, COL_SEASON],
    keep=False,
)

if duplicate_mask.any():
    duplicate_rows = metrics.loc[
        duplicate_mask,
        [COL_SEA, COL_SEASON],
    ]

    raise ValueError(
        "Duplicate sea-season rows were found:\n"
        f"{duplicate_rows.to_string(index=False)}"
    )

expected_combinations = pd.MultiIndex.from_product(
    [SEA_ORDER, SEASON_ORDER],
    names=[COL_SEA, COL_SEASON],
)

actual_combinations = pd.MultiIndex.from_frame(
    metrics[[COL_SEA, COL_SEASON]]
)

missing_combinations = expected_combinations.difference(
    actual_combinations
)

if len(missing_combinations) > 0:
    raise ValueError(
        "The metrics table does not contain all 52 sea-season combinations.\n"
        "Missing combinations:\n"
        + "\n".join(
            f"{sea} — {season}"
            for sea, season in missing_combinations
        )
    )

if len(metrics) != 52:
    raise ValueError(
        f"Expected 52 rows, but found {len(metrics)} rows."
    )

# Sort in the exact study order.
metrics[COL_SEA] = pd.Categorical(
    metrics[COL_SEA],
    categories=SEA_ORDER,
    ordered=True,
)

metrics[COL_SEASON] = pd.Categorical(
    metrics[COL_SEASON],
    categories=SEASON_ORDER,
    ordered=True,
)

metrics = (
    metrics
    .sort_values([COL_SEA, COL_SEASON])
    .reset_index(drop=True)
)


# 8. VERIFY THAT PANEL VARIABLES CONTAIN ALL 13 SEAS
panel_columns = {
    "Figure 6a": [COL_OCCURRENCE],
    "Figure 6b": [COL_WW_TEMPERATURE],
    "Figure 6c": [
        COL_SUBSURFACE_TEMPERATURE,
        COL_THERMOSTERIC,
    ],
    "Figure 6d": [
        COL_FRESHWATER,
        COL_HALOSTERIC,
    ],
}

for panel_name, columns in panel_columns.items():

    valid_rows = metrics.dropna(
        subset=columns
    )

    valid_seas = set(
        valid_rows[COL_SEA].astype(str)
    )

    absent_seas = [
        sea
        for sea in SEA_ORDER
        if sea not in valid_seas
    ]

    if absent_seas:
        raise ValueError(
            f"{panel_name} is missing data for: "
            + ", ".join(absent_seas)
        )


# 9. SUPPORT FUNCTIONS
def padded_limits(
    values,
    lower_fraction=0.06,
    upper_fraction=0.06,
    hard_lower=None,
    hard_upper=None,
):
    """
    Return axis limits with padding so edge points are not clipped.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    minimum = float(np.min(values))
    maximum = float(np.max(values))

    value_range = maximum - minimum

    if value_range == 0:
        value_range = max(abs(maximum), 1.0)

    lower = minimum - lower_fraction * value_range
    upper = maximum + upper_fraction * value_range

    if hard_lower is not None:
        lower = min(lower, hard_lower)

    if hard_upper is not None:
        upper = max(upper, hard_upper)

    return lower, upper


def scale_profile_count(
    profile_count,
    minimum_size=MIN_SCATTER_SIZE,
    maximum_size=MAX_SCATTER_SIZE,
):
    """
    Scale marker area by the square root of profile count.
    This preserves the profile-count information without allowing
    the largest Ross Sea values to dominate the figure.
    """
    counts = np.asarray(
        metrics[COL_PROFILE_COUNT],
        dtype=float,
    )

    transformed_counts = np.sqrt(counts)
    transformed_value = np.sqrt(float(profile_count))

    minimum_count = np.nanmin(transformed_counts)
    maximum_count = np.nanmax(transformed_counts)

    if np.isclose(minimum_count, maximum_count):
        return (minimum_size + maximum_size) / 2

    scaled = (
        minimum_size
        + (
            (transformed_value - minimum_count)
            / (maximum_count - minimum_count)
        )
        * (maximum_size - minimum_size)
    )

    return float(scaled)


def add_panel_letter_title(
    axis,
    panel_letter,
    title,
):
    """
    Add a left-aligned panel title without an overall figure title.
    """
    axis.set_title(
        f"({panel_letter}) {title}",
        loc="left",
        fontsize=FONT_PANEL_TITLE,
        fontweight="bold",
        pad=9,
    )


def style_axis(axis):
    """
    Apply a consistent publication-style appearance.
    """
    axis.grid(
        True,
        which="major",
        color="0.83",
        linewidth=0.65,
        alpha=0.72,
        zorder=0,
    )

    axis.set_axisbelow(True)

    axis.tick_params(
        axis="both",
        which="major",
        labelsize=FONT_TICK,
        width=0.9,
        length=4.5,
        direction="out",
    )

    # Make all major tick labels bold.
    for tick in axis.get_xticklabels() + axis.get_yticklabels():
        tick.set_fontweight("bold")

    for spine in axis.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("0.15")


def get_regression_row(panel_name):
    """
    Retrieve stored regression statistics without recalculating them.
    """
    selected = regression.loc[
        regression["panel"].astype(str).str.strip() == panel_name
    ]

    if selected.empty:
        raise ValueError(
            f"No stored regression row was found for {panel_name}."
        )

    return selected.iloc[0]


def add_regression_line_and_statistics(
    axis,
    panel_name,
    x_values,
):
    """
    Draw the regression line and statistics using the existing
    Fig06_regression_statistics.csv file.
    """
    result = get_regression_row(panel_name)

    x_values = np.asarray(x_values, dtype=float)
    x_values = x_values[np.isfinite(x_values)]

    x_line = np.linspace(
        np.min(x_values),
        np.max(x_values),
        300,
    )

    slope = float(result["regression_slope"])
    intercept = float(result["intercept"])

    y_line = intercept + slope * x_line

    axis.plot(
        x_line,
        y_line,
        color="0.15",
        linewidth=1.8,
        linestyle="--",
        zorder=2,
        label="Stored linear regression",
    )

    if panel_name == "Figure 6c":
        displayed_slope_unit = "cm °C⁻¹"
    elif panel_name == "Figure 6d":
        displayed_slope_unit = "cm m⁻¹"
    else:
        displayed_slope_unit = str(result["slope_unit"])

    statistics_text = (
        f"Pearson r = {float(result['pearson_r']):.2f}\n"
        f"Spearman ρ = {float(result['spearman_rho']):.2f}\n"
        f"Slope = {slope:.2f} {displayed_slope_unit}\n"
        f"95% CI = [{float(result['slope_ci95_lower']):.2f}, "
        f"{float(result['slope_ci95_upper']):.2f}]\n"
        f"n = {int(result['n'])}"
    )

    axis.text(
        0.020,
        0.975,
        statistics_text,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=FONT_STATISTICS,
        fontweight="bold",
        linespacing=1.25,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "0.45",
            "linewidth": 0.8,
            "alpha": 0.96,
        },
        zorder=30,
    )


def repel_labels(
    axis,
    texts,
    anchor_points,
    max_iterations=180,
    padding_pixels=2.5,
):
    """
    Separate overlapping direct labels using only Matplotlib.

    The method works in display coordinates, then converts the adjusted
    positions back to data coordinates. It avoids requiring adjustText or
    any additional plotting package.
    """
    figure = axis.figure
    figure.canvas.draw()

    for iteration in range(max_iterations):

        renderer = figure.canvas.get_renderer()
        bounding_boxes = [
            text.get_window_extent(renderer=renderer).expanded(1.08, 1.18)
            for text in texts
        ]

        moved = False

        for first_index in range(len(texts)):
            for second_index in range(first_index + 1, len(texts)):

                first_box = bounding_boxes[first_index]
                second_box = bounding_boxes[second_index]

                if not first_box.overlaps(second_box):
                    continue

                moved = True

                first_center = np.array([
                    (first_box.x0 + first_box.x1) / 2,
                    (first_box.y0 + first_box.y1) / 2,
                ])

                second_center = np.array([
                    (second_box.x0 + second_box.x1) / 2,
                    (second_box.y0 + second_box.y1) / 2,
                ])

                direction = first_center - second_center

                if np.allclose(direction, 0):
                    direction = np.array([1.0, -1.0])

                direction = direction / np.linalg.norm(direction)

                # More vertical than horizontal movement gives cleaner labels.
                displacement = np.array([
                    direction[0] * 1.3,
                    direction[1] * 2.5,
                ]) * padding_pixels

                for text_index, sign in [
                    (first_index, 1.0),
                    (second_index, -1.0),
                ]:
                    current_data = np.array(texts[text_index].get_position())
                    current_display = axis.transData.transform(current_data)
                    new_display = current_display + sign * displacement

                    # Keep the label centre inside the plotting area.
                    new_display[0] = np.clip(
                        new_display[0],
                        axis.bbox.x0 + 12,
                        axis.bbox.x1 - 12,
                    )
                    new_display[1] = np.clip(
                        new_display[1],
                        axis.bbox.y0 + 10,
                        axis.bbox.y1 - 10,
                    )

                    new_data = axis.transData.inverted().transform(new_display)
                    texts[text_index].set_position(new_data)

        if not moved:
            break

        figure.canvas.draw()

    # Draw short leader lines from each label to its sea-track anchor.
    for text, anchor in zip(texts, anchor_points):
        label_position = text.get_position()

        axis.plot(
            [anchor[0], label_position[0]],
            [anchor[1], label_position[1]],
            color="0.55",
            linewidth=0.60,
            alpha=0.75,
            zorder=6,
        )


def draw_scatter_panel(
    axis,
    x_column,
    y_column,
    panel_letter,
    panel_title,
    x_label,
    y_label,
    regression_panel_name,
):
    """
    Draw a clean sea-trajectory scatter-regression panel.

    Encoding:
      • colour = austral season
      • marker area = valid profile count
      • all markers are circles
      • a faint line joins the four seasons belonging to one sea
      • one repelled direct label identifies each sea trajectory

    This keeps all 52 observations while avoiding 13 marker shapes and
    the large sea-marker legend used in the previous version.
    """
    add_panel_letter_title(
        axis,
        panel_letter,
        panel_title,
    )

    x_all = metrics[x_column].to_numpy(dtype=float)
    y_all = metrics[y_column].to_numpy(dtype=float)

    x_lower, x_upper = padded_limits(
        x_all,
        lower_fraction=0.10,
        upper_fraction=0.10,
    )

    y_lower, y_upper = padded_limits(
        y_all,
        lower_fraction=0.13,
        upper_fraction=0.13,
    )

    axis.set_xlim(x_lower, x_upper)
    axis.set_ylim(y_lower, y_upper)

    x_center = float(np.nanmedian(x_all))
    y_center = float(np.nanmedian(y_all))
    x_span = max(x_upper - x_lower, 1e-12)
    y_span = max(y_upper - y_lower, 1e-12)

    label_texts = []
    label_anchors = []

    for sea_index, sea in enumerate(SEA_ORDER):

        sea_data = metrics.loc[
            metrics[COL_SEA].astype(str) == sea
        ].copy()

        sea_data[COL_SEASON] = pd.Categorical(
            sea_data[COL_SEASON].astype(str),
            categories=SEASON_ORDER,
            ordered=True,
        )

        sea_data = sea_data.sort_values(COL_SEASON)

        x_track = sea_data[x_column].to_numpy(dtype=float)
        y_track = sea_data[y_column].to_numpy(dtype=float)

        axis.plot(
            x_track,
            y_track,
            color="0.60",
            linewidth=0.85,
            alpha=0.60,
            zorder=1,
        )

        for season in SEASON_ORDER:

            row = sea_data.loc[
                sea_data[COL_SEASON].astype(str) == season
            ].iloc[0]

            axis.scatter(
                float(row[x_column]),
                float(row[y_column]),
                s=scale_profile_count(row[COL_PROFILE_COUNT]),
                marker="o",
                facecolor=SEASON_COLOR[season],
                edgecolor="white",
                linewidth=0.75,
                alpha=0.95,
                zorder=4,
            )

        # Use the seasonal point farthest from the overall cloud centre as
        # the sea-label anchor. This places labels around the outside of
        # dense clusters more effectively than labelling each centroid.
        normalized_distance = np.sqrt(
            ((x_track - x_center) / x_span) ** 2
            + ((y_track - y_center) / y_span) ** 2
        )

        anchor_index = int(np.nanargmax(normalized_distance))
        anchor_x = float(x_track[anchor_index])
        anchor_y = float(y_track[anchor_index])

        outward_x = (anchor_x - x_center) / x_span
        outward_y = (anchor_y - y_center) / y_span
        outward_norm = np.hypot(outward_x, outward_y)

        if outward_norm < 1e-8:
            angle = 2 * np.pi * sea_index / len(SEA_ORDER)
            outward_x = np.cos(angle)
            outward_y = np.sin(angle)
            outward_norm = 1.0

        outward_x /= outward_norm
        outward_y /= outward_norm

        # Stagger initial label distances slightly to reduce collisions.
        stagger = 1.0 + 0.16 * (sea_index % 3)

        label_x = anchor_x + outward_x * 0.030 * x_span * stagger
        label_y = anchor_y + outward_y * 0.045 * y_span * stagger

        label_texts.append(
            axis.text(
                label_x,
                label_y,
                sea,
                ha="center",
                va="center",
                fontsize=8.5,
                fontweight="bold",
                color="0.15",
                bbox={
                    "boxstyle": "round,pad=0.10",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "linewidth": 0.0,
                    "alpha": 0.70,
                },
                clip_on=True,
                zorder=8,
            )
        )

        label_anchors.append((anchor_x, anchor_y))

    add_regression_line_and_statistics(
        axis,
        regression_panel_name,
        x_all,
    )

    axis.axhline(
        0,
        color="0.30",
        linewidth=0.85,
        linestyle=":",
        zorder=1,
    )

    axis.set_xlabel(
        x_label,
        fontsize=FONT_AXIS_LABEL,
        fontweight="bold",
        labelpad=8,
    )

    axis.set_ylabel(
        y_label,
        fontsize=FONT_AXIS_LABEL,
        fontweight="bold",
        labelpad=8,
    )

    axis.xaxis.set_major_locator(MaxNLocator(nbins=6))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=6))

    style_axis(axis)

    repel_labels(
        axis,
        label_texts,
        label_anchors,
    )


# 10. CREATE THE 2 × 2 FIGURE
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.weight": "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.linewidth": 1.0,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "mathtext.default": "regular",
})

fig, axes = plt.subplots(
    nrows=2,
    ncols=2,
    figsize=FIGURE_SIZE,
    squeeze=False,
)

ax_a = axes[0, 0]
ax_b = axes[0, 1]
ax_c = axes[1, 0]
ax_d = axes[1, 1]

# Reserve space above for the shared seasonal legend and below for a short explanatory note.
fig.subplots_adjust(
    left=0.075,
    right=0.985,
    top=0.895,
    bottom=0.095,
    wspace=0.18,
    hspace=0.29,
)


# 11. PANEL A — mCDW/CDW PROFILE OCCURRENCE
add_panel_letter_title(
    ax_a,
    "a",
    "Seasonal mCDW/CDW profile occurrence",
)

sea_y_position = {
    sea: index
    for index, sea in enumerate(SEA_ORDER)
}

for sea in SEA_ORDER:

    sea_data = metrics.loc[
        metrics[COL_SEA].astype(str) == sea
    ]

    valid_values = sea_data[COL_OCCURRENCE].dropna()

    if not valid_values.empty:
        ax_a.hlines(
            y=sea_y_position[sea],
            xmin=float(valid_values.min()),
            xmax=float(valid_values.max()),
            color="0.58",
            linewidth=0.85,
            zorder=1,
        )

    for season in SEASON_ORDER:

        row = sea_data.loc[
            sea_data[COL_SEASON].astype(str) == season
        ].iloc[0]

        value = float(row[COL_OCCURRENCE])

        ax_a.scatter(
            value,
            sea_y_position[sea]
            + SEASON_VERTICAL_OFFSET[season],
            s=DOT_SIZE_TOP,
            marker="o",
            facecolor=SEASON_COLOR[season],
            edgecolor="white",
            linewidth=EDGE_WIDTH,
            zorder=4,
            clip_on=False,
        )

occurrence_minimum = float(
    metrics[COL_OCCURRENCE].min()
)

# The upper limit deliberately exceeds 100 so dots at exactly 100%
# remain fully visible instead of being clipped by the panel border.
ax_a.set_xlim(
    occurrence_minimum
    - max(0.8, 0.06 * (100.0 - occurrence_minimum)),
    100.55,
)

ax_a.set_ylim(
    len(SEA_ORDER) - 0.50,
    -0.50,
)

ax_a.set_yticks(
    range(len(SEA_ORDER))
)

ax_a.set_yticklabels(
    SEA_ORDER,
    fontsize=FONT_TICK,
)

ax_a.set_xlabel(
    "mCDW/CDW profile occurrence (%)",
    fontsize=FONT_AXIS_LABEL,
    fontweight="bold",
    labelpad=8,
)

ax_a.set_ylabel(
    "Antarctic shelf sea",
    fontsize=FONT_AXIS_LABEL,
    fontweight="bold",
    labelpad=8,
)

ax_a.xaxis.set_major_locator(
    MaxNLocator(nbins=6)
)

style_axis(ax_a)


# 12. PANEL B — WINTER WATER CORE TEMPERATURE
add_panel_letter_title(
    ax_b,
    "b",
    "Seasonal Winter Water core temperature",
)

for sea in SEA_ORDER:

    sea_data = metrics.loc[
        metrics[COL_SEA].astype(str) == sea
    ]

    valid_values = sea_data[
        COL_WW_TEMPERATURE
    ].dropna()

    if not valid_values.empty:
        ax_b.hlines(
            y=sea_y_position[sea],
            xmin=float(valid_values.min()),
            xmax=float(valid_values.max()),
            color="0.58",
            linewidth=0.85,
            zorder=1,
        )

    for season in SEASON_ORDER:

        row = sea_data.loc[
            sea_data[COL_SEASON].astype(str) == season
        ].iloc[0]

        value = float(
            row[COL_WW_TEMPERATURE]
        )

        ax_b.scatter(
            value,
            sea_y_position[sea]
            + SEASON_VERTICAL_OFFSET[season],
            s=DOT_SIZE_TOP,
            marker="o",
            facecolor=SEASON_COLOR[season],
            edgecolor="white",
            linewidth=EDGE_WIDTH,
            zorder=4,
            clip_on=False,
        )

ww_lower, ww_upper = padded_limits(
    metrics[COL_WW_TEMPERATURE],
    lower_fraction=0.05,
    upper_fraction=0.07,
)

ax_b.set_xlim(
    ww_lower,
    ww_upper,
)

ax_b.set_ylim(
    len(SEA_ORDER) - 0.50,
    -0.50,
)

ax_b.set_yticks(
    range(len(SEA_ORDER))
)

ax_b.set_yticklabels(
    SEA_ORDER,
    fontsize=FONT_TICK,
)

ax_b.set_xlabel(
    "Winter Water core temperature (°C)",
    fontsize=FONT_AXIS_LABEL,
    fontweight="bold",
    labelpad=8,
)

ax_b.set_ylabel(
    "Antarctic shelf sea",
    fontsize=FONT_AXIS_LABEL,
    fontweight="bold",
    labelpad=8,
)

ax_b.xaxis.set_major_locator(
    MaxNLocator(nbins=6)
)

style_axis(ax_b)


# 13. PANEL C — SUBSURFACE TEMPERATURE VS THERMOSTERIC HEIGHT
draw_scatter_panel(
    axis=ax_c,
    x_column=COL_SUBSURFACE_TEMPERATURE,
    y_column=COL_THERMOSTERIC,
    panel_letter="c",
    panel_title=(
        "Subsurface thermal structure versus thermosteric height"
    ),
    x_label=(
        "Median profile-wise maximum temperature, "
        "200–1000 dbar (°C)"
    ),
    y_label=(
        "Seasonal thermosteric-height anomaly (cm)"
    ),
    regression_panel_name="Figure 6c",
)


# 14. PANEL D — FRESHWATER CONTENT VS HALOSTERIC HEIGHT
draw_scatter_panel(
    axis=ax_d,
    x_column=COL_FRESHWATER,
    y_column=COL_HALOSTERIC,
    panel_letter="d",
    panel_title=(
        "Freshwater storage versus halosteric height"
    ),
    x_label=(
        "Freshwater content, 0–200 dbar (m)"
    ),
    y_label=(
        "Seasonal halosteric-height anomaly (cm)"
    ),
    regression_panel_name="Figure 6d",
)


# 15. SHARED SEASON LEGEND — ABOVE PANELS, NO OVERLAP
season_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markersize=7.8,
        markerfacecolor=SEASON_COLOR[season],
        markeredgecolor="white",
        markeredgewidth=0.8,
        label=season,
    )
    for season in SEASON_ORDER
]

season_legend = fig.legend(
    handles=season_handles,
    title="Austral season",
    loc="upper center",
    bbox_to_anchor=(0.50, 0.982),
    ncol=4,
    frameon=False,
    fontsize=FONT_LEGEND,
    title_fontsize=FONT_LEGEND,
    handletextpad=0.45,
    columnspacing=1.35,
)

season_legend.get_title().set_fontweight("bold")

for legend_text in season_legend.get_texts():
    legend_text.set_color("0.10")
    legend_text.set_fontweight("bold")


# 16. EXPLANATORY NOTE FOR PANELS C AND D
fig.text(
    0.50,
    0.025,
    (
        "Panels (c) and (d): colours denote austral season; "
        "grey trajectories connect the four seasonal observations for each sea, "
        "and labels identify individual seas."
    ),
    ha="center",
    va="center",
    fontsize=8.8,
    fontweight="bold",
    color="0.25",
)


# 17. PROFILE-COUNT SIZE KEY — PANEL C
size_key_values = [
    int(
        np.nanpercentile(
            metrics[COL_PROFILE_COUNT],
            15,
        )
    ),
    int(
        np.nanmedian(
            metrics[COL_PROFILE_COUNT]
        )
    ),
    int(
        np.nanpercentile(
            metrics[COL_PROFILE_COUNT],
            90,
        )
    ),
]

# Remove accidental duplicate values while preserving order.
size_key_values = list(
    dict.fromkeys(size_key_values)
)

size_handles = [
    ax_c.scatter(
        [],
        [],
        s=scale_profile_count(value),
        marker="o",
        facecolor="0.72",
        edgecolor="white",
        linewidth=0.8,
        label=f"{value:,}",
    )
    for value in size_key_values
]

profile_legend = ax_c.legend(
    handles=size_handles,
    title="Valid profiles",
    loc="lower right",
    frameon=True,
    fancybox=False,
    facecolor="white",
    edgecolor="0.55",
    framealpha=0.95,
    fontsize=8.0,
    title_fontsize=8.4,
    borderpad=0.45,
    labelspacing=0.45,
    handletextpad=0.55,
)

profile_legend.get_title().set_fontweight("bold")

for txt in profile_legend.get_texts():
    txt.set_fontweight("bold")

profile_legend.get_frame().set_linewidth(0.75)


# 18. SAVE AND DISPLAY
fig.savefig(
    OUTPUT_PNG,
    dpi=SAVE_DPI,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
    pad_inches=0.12,
)

fig.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
    pad_inches=0.12,
)

plt.show()
plt.close(fig)


# 19. FINAL CHECKS
presence_table = pd.crosstab(
    metrics[COL_SEA],
    metrics[COL_SEASON],
).reindex(
    index=SEA_ORDER,
    columns=SEASON_ORDER,
)

print("\nSea-season presence table:")
print(presence_table)

print("\nAll required records are present:")
print(
    f"{metrics[COL_SEA].nunique()} seas × "
    f"{metrics[COL_SEASON].nunique()} seasons = "
    f"{len(metrics)} points per paired analysis."
)

print("\nClean trajectory figure saved successfully:")
print(OUTPUT_PNG)

print("\nVector PDF saved successfully:")
print(OUTPUT_PDF)
