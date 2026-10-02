"""
FIGURE 8 — FULL ORGANIZED LATEST WORKFLOW

FIGURE

EN4 Seasonal Steric Decomposition
Southern Ocean, 2008-2025

SOURCE

Extracted from:
    so_sealevel_paper_fig.py

VERSION SELECTION
The source file contains one Figure 8 implementation. It is already the
layout-corrected/final Figure 8 version, identified by these figure outputs:

    Figure08_EN4_steric_decomposition_LAYOUT_CORRECTED_1080dpi.png
    Figure08_EN4_steric_decomposition_LAYOUT_CORRECTED.pdf

Therefore, no older Figure 8 plotting version is retained in this standalone
script.

PURPOSE

This script performs the complete Figure 8 analysis from first to last:

    1. Mount Google Drive.
    2. Install missing Python packages.
    3. Import all required libraries.
    4. Define input and output files.
    5. Define study period, seasons and plotting settings.
    6. Verify required NetCDF inputs.
    7. Open total, thermosteric and halosteric EN4 products.
    8. Standardize coordinates and longitude.
    9. Restrict data to the Southern Ocean study period/domain.
   10. Calculate seasonal climatological anomalies.
   11. Apply common finite masks to total/thermo/halo fields.
   12. Calculate steric closure:
           residual = total - (thermosteric + halosteric)
   13. Calculate seasonal and overall closure statistics.
   14. Save closure metrics to CSV.
   15. Save seasonal anomaly fields to NetCDF.
   16. Prepare GEBCO bathymetry and the 1000 m contour.
   17. Create the final 4 x 3 polar-map Figure 8.
   18. Save the final 1080-dpi PNG and vector PDF.
   19. Save a text processing/closure summary.
   20. Close all opened datasets.

INPUT FILES
1. EN4 total steric height:
   /content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
   EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc

2. EN4 thermosteric height:
   /content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
   EN4_thermosteric_0_1000m_monthly_2008_2025_SO.nc

3. EN4 halosteric height:
   /content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
   EN4_halosteric_0_1000m_monthly_2008_2025_SO.nc

4. GEBCO bathymetry:
   /content/drive/MyDrive/SAM_Thesis/Data/GEBCO_2024_CEC.nc
   or
   /content/drive/MyDrive/SAM_Thesis/Data/GEBCO_2024_CF.nc

FINAL OUTPUTS

/content/drive/MyDrive/SAM_Thesis/paper2/

    Figure08_EN4_steric_decomposition_LAYOUT_CORRECTED_1080dpi.png
    Figure08_EN4_steric_decomposition_LAYOUT_CORRECTED.pdf
    Fig08_EN4_steric_decomposition_closure_metrics.csv
    Fig08_EN4_steric_decomposition_seasonal_anomalies.nc
    Fig08_EN4_steric_decomposition_summary.txt

SCIENTIFIC DESIGN PRESERVED
The original Figure 8 methodology is retained:

    Study period:
        2008-2025

    Domain:
        Southern Ocean / Antarctic region south of 60 degrees S

    Seasons:
        Spring = SON
        Summer = DJF
        Autumn = MAM
        Winter = JJA

    Seasonal anomaly:
        seasonal climatology - full-period all-month mean

    Steric closure residual:
        total - (thermosteric + halosteric)

    Figure columns:
        Column 1 = total steric-height anomaly
        Column 2 = thermosteric-height anomaly
        Column 3 = halosteric-height anomaly

    Closure tolerance:
        ±1 cm

    Bathymetric reference:
        1000 m contour

ORGANIZATION CHANGES

Only code organization has been changed:
    - Google Drive mounting is moved to the beginning.
    - Package installation follows the Drive mount.
    - Library imports follow installation.
    - The complete original Figure 8 scientific workflow follows unchanged.
    - Figure 9 and all later code are excluded.

No Figure 8 scientific thresholds, equations, file names, seasonal definitions,
output names or plotting logic are intentionally changed.

"""


# ORIGINAL FIGURE 8 SCIENTIFIC DESCRIPTION

# FIGURE 8 — EN4 SEASONAL STERIC DECOMPOSITION
# Southern Ocean, 2008–2025
#
# Inputs
# EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc
# EN4_thermosteric_0_1000m_monthly_2008_2025_SO.nc
# EN4_halosteric_0_1000m_monthly_2008_2025_SO.nc
# GEBCO_2024_CEC.nc
#
# Output
# /content/drive/MyDrive/SAM_Thesis/paper2/
#
# Main calculation
# For each product and each season:
#
# seasonal anomaly = seasonal climatology - all-month mean
#
# where the all-month mean is calculated from the full 2008–2025
# monthly series.
#
# Closure residual
# epsilon = total - (thermo + halo)
#
# Reported closure metrics are calculated from all four seasonal
# residual fields combined, using the common finite mask.
#
# Columns
# Column 1: Total steric height anomaly
# Column 2: Thermosteric height anomaly
# Column 3: Halosteric height anomaly
#
# One symmetric colour scale is used per column.



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


# 4. COMPLETE FIGURE 8 PROCESSING, OUTPUT CREATION AND FINAL PLOTTING

# 3. INPUT AND OUTPUT PATHS
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
    / "Figure08_EN4_steric_decomposition_LAYOUT_CORRECTED_1080dpi.png"
)

OUT_PDF = (
    OUTPUT_DIRECTORY
    / "Figure08_EN4_steric_decomposition_LAYOUT_CORRECTED.pdf"
)

OUT_CLOSURE_CSV = (
    OUTPUT_DIRECTORY
    / "Fig08_EN4_steric_decomposition_closure_metrics.csv"
)

OUT_SEASONAL_NC = (
    OUTPUT_DIRECTORY
    / "Fig08_EN4_steric_decomposition_seasonal_anomalies.nc"
)

OUT_SUMMARY_TXT = (
    OUTPUT_DIRECTORY
    / "Fig08_EN4_steric_decomposition_summary.txt"
)


# 4. ANALYSIS SETTINGS
LAT_MIN = -90.0
LAT_MAX = -60.0

YEAR_START = 2008
YEAR_END = 2025

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
    "Spring": "Spring",
    "Summer": "Summer",
    "Autumn": "Autumn",
    "Winter": "Winter",
}

COLUMN_LABELS = [
    "Total steric height\nanomaly (cm)",
    "Thermosteric height\nanomaly (cm)",
    "Halosteric height\nanomaly (cm)",
]

PANEL_LETTERS = [
    ["(a)", "(b)", "(c)"],
    ["(d)", "(e)", "(f)"],
    ["(g)", "(h)", "(i)"],
    ["(j)", "(k)", "(l)"],
]

MERIDIANS = np.arange(-180, 181, 30)
PARALLELS = [-60, -70, -80]

DRAW_1000_M_CONTOUR = True
BATHYMETRY_CONTOUR_LEVEL = -1000.0

CLOSURE_TOLERANCE_CM = 1.0

SAVE_DPI = 1080


# 5. FIGURE SETTINGS
FIGSIZE = (12.8, 12.0)

FONT_COLUMN_TITLE = 15.2
FONT_PANEL_LETTER = 20.5
FONT_ROW_LABEL = 18.5
FONT_GEO_LABEL = 9.5
FONT_BOX_TITLE = 13.5
FONT_BOX_TEXT = 11.0
FONT_COLORBAR = 12.2
FONT_COLORBAR_TICK = 8.8
FONT_LEGEND = 10.8

NO_DATA_COLOR = "0.84"
LAND_COLOR = "0.78"
COAST_COLOR = "0.35"
GRID_COLOR = "0.56"

RIM_LABEL_LATITUDE = -58.25

# Panel letters are positioned well to the upper-left of each map,
# outside the 30°W rim label.
PANEL_LETTER_X = -0.145
PANEL_LETTER_Y = 1.045

# Extra vertical separation between the top maps and the two-line
# column headings.
COLUMN_TITLE_OFFSET = 0.052

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.weight": "bold",
    "axes.titleweight": "bold",
    "axes.labelweight": "bold",
    "axes.linewidth": 1.1,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "mathtext.default": "regular",
})


# 6. GENERIC HELPERS
def safe_is_datetime(dtype):
    try:
        return np.issubdtype(dtype, np.datetime64)
    except TypeError:
        return False


def open_dataset_safely(file_path, chunks=None):
    """
    Open a NetCDF file using several possible engines.
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
                    f"engine={engine}, decode_times={decode_times}: {error}"
                )

    raise RuntimeError(
        f"Could not open:\n{file_path}\n\n"
        + "\n".join(attempts)
    )


def detect_coordinate(dataset, coordinate_type):
    """
    Detect time, latitude, and longitude coordinate names.
    """
    names = list(dataset.coords) + [
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
            "No suitable science variable was found."
        )

    return variables[0]


def standardize_dataset_coordinates(
    dataset,
    time_name=None,
    lat_name=None,
    lon_name=None,
):
    """
    Rename coordinates to canonical time/lat/lon names.
    """
    rename_mapping = {}

    if time_name is not None and time_name != "time":
        rename_mapping[time_name] = "time"

    if lat_name is not None and lat_name != "lat":
        rename_mapping[lat_name] = "lat"

    if lon_name is not None and lon_name != "lon":
        rename_mapping[lon_name] = "lon"

    if rename_mapping:
        dataset = dataset.rename(rename_mapping)

    return dataset


def normalize_and_sort_longitude(dataset):
    """
    Normalize longitude to -180 ... 180 and remove duplicates.
    """
    if "lon" not in dataset.coords:
        return dataset

    normalized_lon = (
        (dataset["lon"].astype(float) + 180.0) % 360.0
    ) - 180.0

    dataset = dataset.assign_coords(lon=normalized_lon)

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
    Subset a dataset to Antarctic latitudes and sort ascending.
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
    Restrict the analysis period to 2008–2025.
    """
    return data_array.sel(
        time=slice(
            f"{YEAR_START}-01-01",
            f"{YEAR_END}-12-31",
        )
    )


def convert_height_to_cm(data_array):
    """
    Convert metres to centimetres when needed.
    """
    units = str(
        data_array.attrs.get("units", "")
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

    sample = data_array.isel(
        time=slice(
            0,
            min(
                12,
                data_array.sizes.get("time", 1),
            )
        )
    )

    sample_values = np.asarray(
        sample.values,
        dtype=float,
    )
    sample_values = sample_values[
        np.isfinite(sample_values)
    ]

    if (
        sample_values.size > 0
        and np.nanpercentile(
            np.abs(sample_values),
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
    Calculate a robust rounded symmetric colour limit.
    """
    collected = []

    for field in fields:
        values = np.asarray(
            field.values,
            dtype=float,
        )

        values = values[
            np.isfinite(values)
        ]

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
        interval * np.ceil(robust_maximum / interval)
    )


def coordinate_mesh(data_array):
    """
    Return 2-D longitude and latitude arrays.
    """
    return np.meshgrid(
        data_array["lon"].values,
        data_array["lat"].values,
    )


def format_longitude(longitude):
    """
    Format longitude for circular-map rim labels.
    """
    value = int(longitude)

    if value == 0:
        return "0°"

    if abs(value) == 180:
        return "180°"

    if value < 0:
        return f"{abs(value)}°W"

    return f"{value}°E"


# 7. VERIFY INPUT FILES
for file_path in [
    TOTAL_FILE,
    THERMO_FILE,
    HALO_FILE,
    GEBCO_FILE,
]:
    if not os.path.exists(file_path):
        raise FileNotFoundError(
            f"Input file not found:\n{file_path}"
        )

print("\nUsing GEBCO file:")
print(GEBCO_FILE)


# 8. OPEN DATASETS
ds_total = open_dataset_safely(
    TOTAL_FILE,
    chunks="auto",
)

ds_thermo = open_dataset_safely(
    THERMO_FILE,
    chunks="auto",
)

ds_halo = open_dataset_safely(
    HALO_FILE,
    chunks="auto",
)

ds_gebco = open_dataset_safely(
    GEBCO_FILE,
    chunks="auto",
)


# 9. DETECT COORDINATES AND VARIABLES
total_time_name = detect_coordinate(ds_total, "time")
total_lat_name = detect_coordinate(ds_total, "lat")
total_lon_name = detect_coordinate(ds_total, "lon")

thermo_time_name = detect_coordinate(ds_thermo, "time")
thermo_lat_name = detect_coordinate(ds_thermo, "lat")
thermo_lon_name = detect_coordinate(ds_thermo, "lon")

halo_time_name = detect_coordinate(ds_halo, "time")
halo_lat_name = detect_coordinate(ds_halo, "lat")
halo_lon_name = detect_coordinate(ds_halo, "lon")

gebco_lat_name = detect_coordinate(ds_gebco, "lat")
gebco_lon_name = detect_coordinate(ds_gebco, "lon")

total_variable_name = choose_data_variable(
    ds_total,
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

thermo_variable_name = choose_data_variable(
    ds_thermo,
    preferred_names=[
        "thermosteric_height",
        "thermosteric",
        "steric",
    ],
    excluded_names=[
        "valid_layer_mask",
        "gebco_water_depth",
    ],
)

halo_variable_name = choose_data_variable(
    ds_halo,
    preferred_names=[
        "halosteric_height",
        "halosteric",
        "steric",
    ],
    excluded_names=[
        "valid_layer_mask",
        "gebco_water_depth",
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
print(f"Total steric : {total_variable_name}")
print(f"Thermosteric : {thermo_variable_name}")
print(f"Halosteric   : {halo_variable_name}")
print(f"Bathymetry   : {gebco_variable_name}")


# 10. STANDARDIZE COORDINATES
ds_total = standardize_dataset_coordinates(
    ds_total,
    time_name=total_time_name,
    lat_name=total_lat_name,
    lon_name=total_lon_name,
)

ds_thermo = standardize_dataset_coordinates(
    ds_thermo,
    time_name=thermo_time_name,
    lat_name=thermo_lat_name,
    lon_name=thermo_lon_name,
)

ds_halo = standardize_dataset_coordinates(
    ds_halo,
    time_name=halo_time_name,
    lat_name=halo_lat_name,
    lon_name=halo_lon_name,
)

ds_gebco = standardize_dataset_coordinates(
    ds_gebco,
    lat_name=gebco_lat_name,
    lon_name=gebco_lon_name,
)

ds_total = subset_antarctic_latitudes(
    normalize_and_sort_longitude(ds_total)
)

ds_thermo = subset_antarctic_latitudes(
    normalize_and_sort_longitude(ds_thermo)
)

ds_halo = subset_antarctic_latitudes(
    normalize_and_sort_longitude(ds_halo)
)

ds_gebco = subset_antarctic_latitudes(
    normalize_and_sort_longitude(ds_gebco),
    southern_limit=-90.0,
    northern_limit=-55.0,
)


# 11. PREPARE SCIENCE VARIABLES
total = convert_height_to_cm(
    restrict_time_period(ds_total[total_variable_name])
)

thermo = convert_height_to_cm(
    restrict_time_period(ds_thermo[thermo_variable_name])
)

halo = convert_height_to_cm(
    restrict_time_period(ds_halo[halo_variable_name])
)

gebco = ds_gebco[gebco_variable_name]

print("\nUnits after conversion:")
print("Total steric :", total.attrs.get("units"))
print("Thermosteric :", thermo.attrs.get("units"))
print("Halosteric   :", halo.attrs.get("units"))


# 12. ALIGN THE COMMON MONTHLY PERIOD
common_times = np.intersect1d(
    pd.to_datetime(total["time"].values),
    np.intersect1d(
        pd.to_datetime(thermo["time"].values),
        pd.to_datetime(halo["time"].values),
    ),
)

if common_times.size == 0:
    raise ValueError(
        "No common monthly period among total, thermo, and halo files."
    )

total = total.sel(time=common_times)
thermo = thermo.sel(time=common_times)
halo = halo.sel(time=common_times)

number_of_months = len(common_times)

print("\nCommon analysis period:")
print(pd.to_datetime(common_times[0]))
print("to")
print(pd.to_datetime(common_times[-1]))
print(f"Months = {number_of_months}")


# 13. BUILD A COMMON STATIC MASK
# Use the intersection of finite monthly values across all three
# products during the full period.
total_valid_fraction = (
    total.notnull().sum(dim="time")
    / number_of_months
)

thermo_valid_fraction = (
    thermo.notnull().sum(dim="time")
    / number_of_months
)

halo_valid_fraction = (
    halo.notnull().sum(dim="time")
    / number_of_months
)

static_common_mask = (
    (total_valid_fraction >= 1.0 - 1e-10)
    & (thermo_valid_fraction >= 1.0 - 1e-10)
    & (halo_valid_fraction >= 1.0 - 1e-10)
).compute()

static_common_cell_count = int(
    static_common_mask.sum().values
)

if static_common_cell_count == 0:
    raise RuntimeError(
        "The common static mask contains no valid cells."
    )

print("\nStatic common mask cells:", f"{static_common_cell_count:,}")


# 14. CALCULATE SEASONAL ANOMALIES
total_all_month_mean = total.mean(
    dim="time",
    skipna=True,
)

thermo_all_month_mean = thermo.mean(
    dim="time",
    skipna=True,
)

halo_all_month_mean = halo.mean(
    dim="time",
    skipna=True,
)

total_season_mean = seasonal_climatology(total)
thermo_season_mean = seasonal_climatology(thermo)
halo_season_mean = seasonal_climatology(halo)

total_anomaly = {}
thermo_anomaly = {}
halo_anomaly = {}
seasonal_mask = {}
residual_field = {}
no_data_mask = {}

for season in SEASON_ORDER:

    total_field = (
        total_season_mean[season]
        - total_all_month_mean
    )

    thermo_field = (
        thermo_season_mean[season]
        - thermo_all_month_mean
    )

    halo_field = (
        halo_season_mean[season]
        - halo_all_month_mean
    )

    row_mask = (
        static_common_mask
        & np.isfinite(total_field)
        & np.isfinite(thermo_field)
        & np.isfinite(halo_field)
    ).compute()

    total_anomaly[season] = total_field.where(
        row_mask
    ).compute()

    thermo_anomaly[season] = thermo_field.where(
        row_mask
    ).compute()

    halo_anomaly[season] = halo_field.where(
        row_mask
    ).compute()

    seasonal_mask[season] = row_mask

    residual = (
        total_anomaly[season]
        - (
            thermo_anomaly[season]
            + halo_anomaly[season]
        )
    )

    residual_field[season] = residual.where(row_mask)

    no_data_mask[season] = xr.where(
        row_mask,
        np.nan,
        1.0,
    )


# 15. CLOSURE METRICS
all_residual_values = []

seasonal_closure_rows = []

for season in SEASON_ORDER:

    values = np.asarray(
        residual_field[season].values,
        dtype=float,
    )

    values = values[np.isfinite(values)]

    if values.size > 0:
        all_residual_values.append(values)

        seasonal_closure_rows.append(
            {
                "season": season,
                "mean_residual_cm": float(np.mean(values)),
                "median_absolute_residual_cm": float(
                    np.median(np.abs(values))
                ),
                "RMSE_cm": float(
                    np.sqrt(np.mean(values ** 2))
                ),
                "maximum_absolute_residual_cm": float(
                    np.max(np.abs(values))
                ),
                "cells_within_tolerance_percent": float(
                    100.0
                    * np.mean(
                        np.abs(values)
                        <= CLOSURE_TOLERANCE_CM
                    )
                ),
                "n_cells": int(values.size),
            }
        )

if not all_residual_values:
    raise RuntimeError(
        "No finite residual values were available for closure analysis."
    )

combined_residual_values = np.concatenate(
    all_residual_values
)

overall_closure_metrics = {
    "scope": "All seasons combined",
    "tolerance_cm": CLOSURE_TOLERANCE_CM,
    "mean_residual_cm": float(
        np.mean(combined_residual_values)
    ),
    "median_absolute_residual_cm": float(
        np.median(
            np.abs(combined_residual_values)
        )
    ),
    "RMSE_cm": float(
        np.sqrt(
            np.mean(combined_residual_values ** 2)
        )
    ),
    "maximum_absolute_residual_cm": float(
        np.max(
            np.abs(combined_residual_values)
        )
    ),
    "cells_within_tolerance_percent": float(
        100.0
        * np.mean(
            np.abs(combined_residual_values)
            <= CLOSURE_TOLERANCE_CM
        )
    ),
    "n_cells": int(
        combined_residual_values.size
    ),
}

closure_metrics_table = pd.concat(
    [
        pd.DataFrame([overall_closure_metrics]),
        pd.DataFrame(seasonal_closure_rows),
    ],
    ignore_index=True,
)

closure_metrics_table.to_csv(
    OUT_CLOSURE_CSV,
    index=False,
)

print("\nClosure metrics:")
print(
    closure_metrics_table.to_string(index=False)
)


# 16. SAVE SEASONAL FIELDS
season_coordinate = xr.DataArray(
    SEASON_ORDER,
    dims="season",
    name="season",
)

total_stack = xr.concat(
    [total_anomaly[season] for season in SEASON_ORDER],
    dim=season_coordinate,
)

thermo_stack = xr.concat(
    [thermo_anomaly[season] for season in SEASON_ORDER],
    dim=season_coordinate,
)

halo_stack = xr.concat(
    [halo_anomaly[season] for season in SEASON_ORDER],
    dim=season_coordinate,
)

residual_stack = xr.concat(
    [residual_field[season] for season in SEASON_ORDER],
    dim=season_coordinate,
)

mask_stack = xr.concat(
    [
        seasonal_mask[season].astype(np.int8)
        for season in SEASON_ORDER
    ],
    dim=season_coordinate,
)

seasonal_output = xr.Dataset(
    {
        "total_steric_anomaly": total_stack,
        "thermosteric_anomaly": thermo_stack,
        "halosteric_anomaly": halo_stack,
        "closure_residual": residual_stack,
        "seasonal_common_mask": mask_stack,
    }
)

seasonal_output["total_steric_anomaly"].attrs.update({
    "long_name": (
        "seasonal total steric-height anomaly "
        "relative to the 2008-2025 all-month mean"
    ),
    "units": "cm",
})

seasonal_output["thermosteric_anomaly"].attrs.update({
    "long_name": (
        "seasonal thermosteric-height anomaly "
        "relative to the 2008-2025 all-month mean"
    ),
    "units": "cm",
})

seasonal_output["halosteric_anomaly"].attrs.update({
    "long_name": (
        "seasonal halosteric-height anomaly "
        "relative to the 2008-2025 all-month mean"
    ),
    "units": "cm",
})

seasonal_output["closure_residual"].attrs.update({
    "long_name": (
        "closure residual: total - (thermo + halo)"
    ),
    "units": "cm",
})

seasonal_output["seasonal_common_mask"].attrs.update({
    "long_name": "season-specific common finite mask",
    "flag_values": np.array([0, 1], dtype=np.int8),
    "flag_meanings": "excluded included",
})

seasonal_output.attrs.update({
    "title": "Figure 8 EN4 seasonal steric decomposition",
    "closure_definition": "total - (thermo + halo)",
    "time_period": f"{YEAR_START}-01 to {YEAR_END}-12",
})

seasonal_output.rename(
    {
        "lat": "latitude",
        "lon": "longitude",
    }
).to_netcdf(
    OUT_SEASONAL_NC,
    encoding={
        "total_steric_anomaly": {
            "zlib": True,
            "complevel": 4,
        },
        "thermosteric_anomaly": {
            "zlib": True,
            "complevel": 4,
        },
        "halosteric_anomaly": {
            "zlib": True,
            "complevel": 4,
        },
        "closure_residual": {
            "zlib": True,
            "complevel": 4,
        },
        "seasonal_common_mask": {
            "zlib": True,
            "complevel": 4,
        },
    },
)


# 17. COLUMN-SPECIFIC COLOUR LIMITS
total_limit = symmetric_limit(
    [total_anomaly[season] for season in SEASON_ORDER]
)

thermo_limit = symmetric_limit(
    [thermo_anomaly[season] for season in SEASON_ORDER]
)

halo_limit = symmetric_limit(
    [halo_anomaly[season] for season in SEASON_ORDER]
)

print("\nColour limits:")
print(f"Total steric : ±{total_limit:.2f} cm")
print(f"Thermosteric : ±{thermo_limit:.2f} cm")
print(f"Halosteric   : ±{halo_limit:.2f} cm")


# 18. PREPARE GEBCO FOR THE 1000 M CONTOUR
if (
    gebco.ndim == 2
    and gebco.size > 1_200_000
):
    stride = int(
        np.ceil(
            np.sqrt(
                gebco.size / 1_200_000
            )
        )
    )

    gebco = gebco.isel(
        lat=slice(None, None, stride),
        lon=slice(None, None, stride),
    )

gebco = gebco.compute()

gebco_lon2d, gebco_lat2d = coordinate_mesh(gebco)


# 19. MAP-DRAWING HELPERS
def add_circular_boundary(axis):
    """
    Apply a circular boundary to the polar map.
    """
    theta = np.linspace(
        0,
        2 * np.pi,
        500,
    )

    vertices = np.vstack(
        [np.sin(theta), np.cos(theta)]
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
    Draw rim longitude labels and two latitude labels.
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
    Common South Polar map style.
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
    no_data,
    color_limit,
    panel_letter,
):
    """
    Plot one map panel.
    """
    style_map_axis(axis)

    longitude_2d, latitude_2d = coordinate_mesh(data_array)

    axis.pcolormesh(
        longitude_2d,
        latitude_2d,
        no_data.values,
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
            levels=[BATHYMETRY_CONTOUR_LEVEL],
            colors="black",
            linewidths=0.68,
            linestyles="-",
            alpha=0.92,
            transform=ccrs.PlateCarree(),
            zorder=9,
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


# 20. BUILD THE FIGURE
projection = ccrs.SouthPolarStereo(
    central_longitude=0.0
)

fig = plt.figure(figsize=FIGSIZE)

grid = gridspec.GridSpec(
    nrows=4,
    ncols=4,
    figure=fig,
    width_ratios=[1.0, 1.0, 1.0, 0.52],
    height_ratios=[1, 1, 1, 1],
    left=0.080,
    right=0.990,
    top=0.925,
    bottom=0.135,
    wspace=-0.035,   # negative value compensates for Cartopy's internal space
    hspace=0.18,
)

map_axes = np.empty((4, 3), dtype=object)

last_images = [None, None, None]

for row_index, season in enumerate(SEASON_ORDER):

    for column_index, field_dictionary in enumerate(
        [
            total_anomaly,
            thermo_anomaly,
            halo_anomaly,
        ]
    ):

        axis = fig.add_subplot(
            grid[row_index, column_index],
            projection=projection,
        )

        map_axes[row_index, column_index] = axis

        if column_index == 0:
            limit = total_limit
        elif column_index == 1:
            limit = thermo_limit
        else:
            limit = halo_limit

        last_images[column_index] = plot_map_field(
            axis=axis,
            data_array=field_dictionary[season],
            no_data=no_data_mask[season],
            color_limit=limit,
            panel_letter=PANEL_LETTERS[row_index][column_index],
        )

# 21. TWO-LINE COLUMN HEADINGS
# Figure-level titles are used instead of axis titles. Their centers
# follow the actual positions of the three top-row maps, preventing
# overlap between neighbouring column headings.
for column_index in range(3):

    column_position = map_axes[
        0,
        column_index,
    ].get_position()

    column_center_x = (
        column_position.x0
        + column_position.x1
    ) / 2.0

    title_y = min(
        0.986,
        column_position.y1
        + COLUMN_TITLE_OFFSET,
    )

    fig.text(
        column_center_x,
        title_y,
        COLUMN_LABELS[column_index],
        ha="center",
        va="top",
        fontsize=FONT_COLUMN_TITLE,
        fontweight="bold",
        linespacing=1.05,
        color="black",
    )


# 23. LEFT-SIDE SEASON LABELS
for row_index, season in enumerate(SEASON_ORDER):

    position = map_axes[row_index, 0].get_position()

    fig.text(
        position.x0 - 0.040,
        position.y0 + position.height / 2,
        SEASON_ROW_LABEL[season],
        rotation=90,
        rotation_mode="anchor",
        ha="center",
        va="center",
        fontsize=FONT_ROW_LABEL,
        fontweight="bold",
        color="black",
    )


# 24. BOTTOM COLOURBARS
bottom_positions = [
    map_axes[3, 0].get_position(),
    map_axes[3, 1].get_position(),
    map_axes[3, 2].get_position(),
]

colorbar_y = 0.080
colorbar_height = 0.012

colorbar_axes = [
    fig.add_axes(
        [
            bottom_positions[i].x0 + 0.010,
            colorbar_y,
            bottom_positions[i].width - 0.020,
            colorbar_height,
        ]
    )
    for i in range(3)
]

limits = [total_limit, thermo_limit, halo_limit]

for i, cax in enumerate(colorbar_axes):

    colorbar = fig.colorbar(
        last_images[i],
        cax=cax,
        orientation="horizontal",
        extend="both",
    )

    colorbar.set_label(
        "(cm)",
        fontsize=FONT_COLORBAR,
        fontweight="bold",
        labelpad=3,
    )

    # Five evenly spaced ticks prevent labels from touching,
    # especially in the thermosteric colourbar.
    colorbar.set_ticks(
        np.linspace(
            -limits[i],
            limits[i],
            5,
        )
    )

    colorbar.ax.tick_params(
        labelsize=FONT_COLORBAR_TICK,
        width=1.0,
        length=4,
    )

    for label in colorbar.ax.get_xticklabels():
        label.set_fontweight("bold")


# 25. SAVE FIGURE AND SUMMARY
fig.savefig(
    OUT_PNG,
    dpi=SAVE_DPI,
    bbox_inches="tight",
    pad_inches=0.08,
    facecolor="white",
    edgecolor="none",
)

fig.savefig(
    OUT_PDF,
    bbox_inches="tight",
    pad_inches=0.08,
    facecolor="white",
    edgecolor="none",
)

plt.show()
plt.close(fig)

summary_lines = [
    "=" * 90,
    "FIGURE 8 — EN4 SEASONAL STERIC DECOMPOSITION",
    "=" * 90,
    "",
    f"PNG : {OUT_PNG}",
    f"PDF : {OUT_PDF}",
    f"Closure CSV : {OUT_CLOSURE_CSV}",
    f"Seasonal NetCDF : {OUT_SEASONAL_NC}",
    "",
    "Overall closure metrics:",
    f"Mean residual (cm)               : {overall_closure_metrics['mean_residual_cm']:.6f}",
    f"Median absolute residual (cm)    : {overall_closure_metrics['median_absolute_residual_cm']:.6f}",
    f"RMSE (cm)                        : {overall_closure_metrics['RMSE_cm']:.6f}",
    f"Maximum absolute residual (cm)   : {overall_closure_metrics['maximum_absolute_residual_cm']:.6f}",
    f"Cells within ±{CLOSURE_TOLERANCE_CM:.0f} cm (%) : {overall_closure_metrics['cells_within_tolerance_percent']:.3f}",
    f"Combined cells                   : {overall_closure_metrics['n_cells']}",
    "",
    "Notes:",
    "• Seasonal anomalies are departures from the all-month mean.",
    "• One symmetric colour scale is used for each column.",
    "• The same common finite seasonal mask is used among total, thermo, and halo in each row.",
    "• Column headings are centred two-line figure labels.",
    "• Panel letters are offset from the 30°W rim labels.",
    "• The information boxes are positioned close to the third map column.",
]

with open(OUT_SUMMARY_TXT, "w", encoding="utf-8") as handle:
    handle.write("\n".join(summary_lines))

print("\n" + "\n".join(summary_lines))

ds_total.close()
ds_thermo.close()
ds_halo.close()
ds_gebco.close()

# Commented out IPython magic to ensure Python compatibility.
