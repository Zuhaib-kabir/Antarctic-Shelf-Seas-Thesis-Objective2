# -*- coding: utf-8 -*-
"""
SUPPLEMENTARY FIGURE S3 — FULL ORGANIZED LATEST WORKFLOW

FIGURE

Argo–EN4 0–1000 m Validation
Southern Ocean, 2008-2025

SOURCE

Extracted and organized from:
    so_sealevel_paper_fig.py

VERSION SELECTION
The source contains two S3 stages:

1. A complete validation workflow labelled:

       HDF-ERROR FIXED
       LOW-RAM
       CHECKPOINTED
       RESUMABLE
       1080-DPI

   This stage performs the Argo–EN4 matching, TEOS-10 calculations,
   summary statistics, caching and CSV creation.

2. A later plot-only correction using the existing precomputed CSV files.

   Final plotting corrections:
       - panel (d) legend row moved LOWER so it does not overlap the x-axis;
       - panel (e) seasonal labels remain visible;
       - no scientific recomputation;
       - output filenames remain unchanged.

This standalone script therefore runs:

    full low-RAM validation science
        -> precomputed S3 CSV files
        -> latest corrected plot-only renderer

The superseded first S3 renderer is excluded.

PANELS

(a) Argo versus EN4 total steric height
(b) Argo versus EN4 thermosteric height
(c) Argo versus EN4 halosteric height
(d) Sea-wise mean bias: Argo - EN4
(e) Seasonal Pearson correlation
(f) Number of matched Argo profiles

RELIABILITY / LOW-RAM CORRECTIONS PRESERVED
1. Large NetCDF inputs are staged from Google Drive to /content first.
2. h5netcdf is attempted before netCDF4.
3. Argo TEMP/PSAL arrays are never loaded in full.
4. Argo profiles are read in contiguous chunks.
5. The fixed EN4 reference profile is calculated once and cached.
6. Completed Argo chunks are cached as CSV files and reused after interruption.
7. The completed matched-profile CSV is reused unless FORCE_RECOMPUTE=True.
8. The fixed 2008-2025 EN4 grid-cell reference profiles and TEOS-10
   decomposition used by the source are retained.

VALIDATION SETTINGS PRESERVED
Period:
    2008-01-01 through 2025-12-31

Layer:
    0-1000 m

Argo profile chunk size:
    1500 profiles

Maximum vertical interpolation gap:
    100 m

Maximum nearest EN4 grid-cell distance:
    120 km

Maximum accepted top depth:
    15 m

Required bottom support:
    at least 1000 m

Temperature guard:
    -3.5 to 15.0 deg C

Salinity guard:
    20.0 to 40.0

MATCHING METHOD

- same calendar month;
- nearest EN4 grid cell within 120 km;
- Argo T/S interpolated to the EN4 0-1000 m analysis-depth grid;
- complete common T/S coverage from near-surface to 1000 m;
- fixed 2008-2025 EN4 grid-cell reference profiles;
- TEOS-10 total, thermosteric and residual-halosteric decomposition.

PRECOMPUTED SCIENCE OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/

    Figure_S3_Argo_EN4_0_1000_validation_matched_profiles.csv
    Figure_S3_Argo_EN4_0_1000_validation_sea_bias.csv
    Figure_S3_Argo_EN4_0_1000_validation_seasonal_correlations.csv
    Figure_S3_Argo_EN4_0_1000_validation_matched_counts.csv

CACHE / RESUME PRODUCTS
    Figure_S3_EN4_fixed_reference_0_1000m.npz
    Figure_S3_matching_chunks_v2/

Large source files may be locally staged in:
    /content/Figure_S3_input_cache/

LATEST FIGURE OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/

    Figure_S3_Argo_EN4_0_1000_validation_1080dpi.png
    Figure_S3_Argo_EN4_0_1000_validation.pdf
    Figure_S3_Argo_EN4_0_1000_validation_processing_log.txt

ORGANIZED EXECUTION ORDER
1. Mount Google Drive.
2. Install required packages.
3. Import all libraries.
4. Define cache/recompute/low-RAM settings.
5. Define Argo, raw-EN4 and EN4-product paths.
6. Define 13 Antarctic sea sectors and seasonal order.
7. Define robust NetCDF/file/time/coordinate helpers.
8. Discover the actual Argo and raw-EN4 inputs.
9. Stage large inputs locally under /content.
10. Read EN4 monthly total/thermo/halo products.
11. Build or load the fixed EN4 0-1000 m reference cache.
12. Process Argo profiles in low-RAM chunks.
13. Match profiles to same-month nearest EN4 cells.
14. Enforce distance, depth-support, vertical-gap and range rules.
15. Compute Argo TEOS-10 steric decomposition.
16. Save/reuse chunk checkpoints.
17. Save the final matched-profile table.
18. Calculate and save sea-wise biases.
19. Calculate and save seasonal correlations.
20. Calculate and save sea-wise matched counts.
21. Read those precomputed tables with the latest renderer.
22. Plot panels (a)-(c) validation scatterplots.
23. Plot panel (d) with its corrected lower legend row.
24. Plot panel (e) with visible seasonal labels.
25. Plot panel (f) matched-profile counts.
26. Save corrected PDF and 1080-dpi PNG.
27. Save the complete scientific processing log.

ORGANIZATION CHANGES
Only organization/version selection has changed:

    - Google Drive mount moved first.
    - Package installation consolidated.
    - Imports consolidated.
    - Full S3 scientific validation retained.
    - Superseded first S3 plot removed.
    - Later corrected plot-only renderer retained.
    - Scientific processing log retained.
    - Figure S4 and all later code excluded.

No scientific thresholds, matching rules, TEOS-10 calculations, cache
behavior, data-table filenames, or final plotting corrections are
intentionally changed.

"""

# 1. MOUNT GOOGLE DRIVE

try:
    from google.colab import drive

    drive.mount(
        "/content/drive",
        force_remount=False,
    )

except Exception:
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
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "gsw": "gsw",
    "h5netcdf": "h5netcdf",
    "netCDF4": "netCDF4",
    "dask": "dask[array]",
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

print("All required packages are available.")


# 3. IMPORT LIBRARIES

import gc
import os
import shutil
import time
import warnings

from collections import OrderedDict
from pathlib import Path

import dask
import gsw
import numpy as np
import pandas as pd
import xarray as xr

from scipy import stats
from scipy.interpolate import interp1d

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
)


# PART A — HDF-ERROR-FIXED LOW-RAM S3 VALIDATION SCIENCE
# Creates the four precomputed S3 validation CSV files.

# FIGURE S3 — ARGO–EN4 0–1000 m VALIDATION
# HDF-ERROR FIXED • LOW-RAM • CHECKPOINTED • RESUMABLE • 1080-DPI
#
# Panels
# (a) Argo versus EN4 total steric height
# (b) Argo versus EN4 thermosteric height
# (c) Argo versus EN4 halosteric height
# (d) Sea-wise mean bias: Argo − EN4
# (e) Seasonal Pearson correlation
# (f) Number of matched Argo profiles
#
# Main corrections compared with the failed workflow
# 1. Large NetCDF inputs are copied from Google Drive to /content first.
#    This avoids intermittent Google-Drive HDF5 read errors.
# 2. h5netcdf is tried before netCDF4.
# 3. The Argo TEMP/PSAL arrays are never loaded in full. Profiles are read
#    in contiguous chunks.
# 4. The EN4 all-month reference profile is calculated once and cached.
# 5. Completed Argo chunks are cached as CSV files and reused after a crash.
# 6. The final matched-profile CSV is reused unless FORCE_RECOMPUTE=True.
# 7. The scientific calculations use the same fixed 2008–2025 grid-cell
#    reference profile and TEOS-10 decomposition used for the EN4 products.
#
# Outputs are written to:
# /content/drive/MyDrive/SAM_Thesis/paper2

# 3. USER SETTINGS
FORCE_RECOMPUTE = False
REUSE_REFERENCE_CACHE = True
REUSE_CHUNK_CACHE = True
STAGE_INPUTS_TO_LOCAL = True

# Chunk size controls RAM usage. 1000–2500 is safe for standard Colab RAM.
ARGO_PROFILE_CHUNK = 1500

# Plot settings
SAVE_DPI = 1080
FIGSIZE = (9.6, 6.9)
ADD_OVERALL_TITLE = False
MAX_SCATTER_POINTS_DISPLAY = 18000
RANDOM_SEED = 42

# Validation period and layer
START_DATE = pd.Timestamp("2008-01-01")
END_DATE = pd.Timestamp("2025-12-31 23:59:59")
MAX_DEPTH_M = 1000.0
MAX_VERTICAL_GAP_M = 100.0
MAX_NEAREST_CELL_DISTANCE_KM = 120.0

# A profile must have common T/S coverage at the top and down to 1000 m.
MAX_ALLOWED_TOP_DEPTH_M = 15.0
MIN_ALLOWED_BOTTOM_DEPTH_M = 1000.0

# Scientific range guard. Values outside these broad limits are rejected.
TEMP_RANGE_C = (-3.5, 15.0)
SAL_RANGE = (20.0, 40.0)

VERBOSE = True


# 4. PATHS — OUTPUT NAMES REMAIN FIXED
BASE_DIR = Path("/content/drive/MyDrive/SAM_Thesis")
DATA_DIR = BASE_DIR / "Data"
EN4_PRODUCT_DIR = BASE_DIR / "Processed" / "EN4_NetCDF_inventory"
OUTPUT_DIR = BASE_DIR / "paper2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

ARGO_CANDIDATES = [
    DATA_DIR / "argo_SO_profiles_2001_2025_cleaned_gridded.nc",
    DATA_DIR / "argo_SO_profiles_2001_2025_cleaned.nc",
]

EN4_RAW_CANDIDATES = [
    DATA_DIR / "EN4" / "final_monthly" / "EN4.2.2_g10_Antarctic_monthly_2008_2025.nc",
    DATA_DIR / "EN4.2.2_g10_Antarctic_monthly_2008_2025.nc",
]

EN4_TOTAL_FILE = EN4_PRODUCT_DIR / "EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc"
EN4_THERMO_FILE = EN4_PRODUCT_DIR / "EN4_thermosteric_0_1000m_monthly_2008_2025_SO.nc"
EN4_HALO_FILE = EN4_PRODUCT_DIR / "EN4_halosteric_0_1000m_monthly_2008_2025_SO.nc"

MATCHED_CSV = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_matched_profiles.csv"
SEA_BIAS_CSV = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_sea_bias.csv"
SEASONAL_CORR_CSV = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_seasonal_correlations.csv"
SEA_COUNT_CSV = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_matched_counts.csv"
FIG_PNG = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_1080dpi.png"
FIG_PDF = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation.pdf"
PROCESSING_LOG = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_processing_log.txt"

REFERENCE_CACHE = OUTPUT_DIR / "Figure_S3_EN4_fixed_reference_0_1000m.npz"
CHUNK_CACHE_DIR = OUTPUT_DIR / "Figure_S3_matching_chunks_v2"
CHUNK_CACHE_DIR.mkdir(parents=True, exist_ok=True)

LOCAL_STAGE_DIR = Path("/content/Figure_S3_input_cache")
LOCAL_STAGE_DIR.mkdir(parents=True, exist_ok=True)


# 5. ANTARCTIC SEA DEFINITIONS
SEA_SECTORS = OrderedDict([
    ("WED", {"name": "Weddell Sea",         "lon_min": -60.0,  "lon_max": -20.0}),
    ("KHV", {"name": "King Haakon VII Sea", "lon_min": -20.0,  "lon_max":  10.0}),
    ("RLS", {"name": "Riiser-Larsen Sea",   "lon_min":  10.0,  "lon_max":  35.0}),
    ("LAZ", {"name": "Lazarev Sea",         "lon_min":  35.0,  "lon_max":  60.0}),
    ("COS", {"name": "Cosmonauts Sea",      "lon_min":  60.0,  "lon_max":  90.0}),
    ("COO", {"name": "Cooperation Sea",     "lon_min":  90.0,  "lon_max": 115.0}),
    ("DAV", {"name": "Davis Sea",           "lon_min": 115.0,  "lon_max": 130.0}),
    ("MAW", {"name": "Mawson Sea",          "lon_min": 130.0,  "lon_max": 150.0}),
    ("DUR", {"name": "D'Urville Sea",       "lon_min": 150.0,  "lon_max": 170.0}),
    ("SOM", {"name": "Somov Sea",           "lon_min": 170.0,  "lon_max": -160.0}),
    ("ROS", {"name": "Ross Sea",            "lon_min": -160.0, "lon_max": -130.0}),
    ("AMU", {"name": "Amundsen Sea",        "lon_min": -130.0, "lon_max": -100.0}),
    ("BEL", {"name": "Bellingshausen Sea",  "lon_min": -100.0, "lon_max": -60.0}),
])

SEA_ORDER = list(SEA_SECTORS.keys())
SEASON_ORDER = ["Spring", "Summer", "Autumn", "Winter"]


# 6. GENERAL HELPERS
def log(message=""):
    if VERBOSE:
        print(message, flush=True)


def first_existing(candidates, label):
    for candidate in candidates:
        candidate = Path(candidate)
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"No {label} file was found. Checked:\n"
        + "\n".join(str(path) for path in candidates)
    )


def require_file(path, label):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{label} file not found:\n{path}")
    return path


def normalize_longitude(values):
    values = np.asarray(values, dtype=float)
    return ((values + 180.0) % 360.0) - 180.0


def assign_sea(longitude):
    longitude = float(normalize_longitude([longitude])[0])
    for code, information in SEA_SECTORS.items():
        west = information["lon_min"]
        east = information["lon_max"]
        if west <= east:
            inside = west <= longitude < east
        else:
            inside = longitude >= west or longitude < east
        if inside:
            return code
    return None


def month_to_season(month):
    month = int(month)
    if month in (9, 10, 11):
        return "Spring"
    if month in (12, 1, 2):
        return "Summer"
    if month in (3, 4, 5):
        return "Autumn"
    return "Winter"


def open_dataset_robust(
    path,
    *,
    decode_times=True,
    chunks=None,
    drop_variables=None,
):
    """Open a NetCDF with multiple engines; h5netcdf is attempted first."""
    path = Path(path)
    attempts = []

    for engine in ("h5netcdf", "netcdf4", None, "scipy"):
        try:
            kwargs = {
                "decode_times": decode_times,
                "mask_and_scale": True,
                "cache": False,
            }
            if engine is not None:
                kwargs["engine"] = engine
            if chunks is not None:
                kwargs["chunks"] = chunks
            if drop_variables is not None:
                kwargs["drop_variables"] = drop_variables

            dataset = xr.open_dataset(path, **kwargs)
            log(f"Opened {path.name} with engine={engine or 'xarray-default'}")
            return dataset
        except Exception as error:
            attempts.append(f"{engine or 'default'}: {type(error).__name__}: {error}")
            gc.collect()

    raise RuntimeError(
        f"Could not open NetCDF file:\n{path}\n\nAttempts:\n"
        + "\n".join(attempts)
    )


def stage_file_to_local(source_path, label, retries=3):
    """Copy a Drive file to /content and verify byte size."""
    source_path = Path(source_path)

    if not STAGE_INPUTS_TO_LOCAL:
        return source_path

    destination = LOCAL_STAGE_DIR / source_path.name
    source_size = source_path.stat().st_size

    if destination.exists() and destination.stat().st_size == source_size:
        log(f"Using staged local {label}: {destination}")
        return destination

    temporary = destination.with_suffix(destination.suffix + ".partial")

    for attempt in range(1, retries + 1):
        try:
            log(
                f"Staging {label} to local runtime "
                f"({source_size / 1024**2:.1f} MB), attempt {attempt}/{retries} ..."
            )
            if temporary.exists():
                temporary.unlink()

            with source_path.open("rb") as source, temporary.open("wb") as target:
                shutil.copyfileobj(source, target, length=16 * 1024 * 1024)
                target.flush()
                os.fsync(target.fileno())

            if temporary.stat().st_size != source_size:
                raise IOError(
                    f"Staged file size mismatch: "
                    f"{temporary.stat().st_size} != {source_size}"
                )

            temporary.replace(destination)
            log(f"Staged local file: {destination}")
            return destination

        except Exception as error:
            log(f"  staging attempt failed: {error}")
            time.sleep(3 * attempt)

    raise RuntimeError(
        f"Could not stage {label} after {retries} attempts:\n{source_path}"
    )


def standardize_3d_product(dataset, variable_name):
    rename = {}
    for source, target in (
        ("latitude", "lat"),
        ("longitude", "lon"),
        ("valid_time", "time"),
    ):
        if source in dataset.coords or source in dataset.dims:
            if target not in dataset.coords and target not in dataset.dims:
                rename[source] = target
    if rename:
        dataset = dataset.rename(rename)

    data = dataset[variable_name]
    data = data.assign_coords(lon=normalize_longitude(data["lon"].values))
    _, unique_indices = np.unique(data["lon"].values, return_index=True)
    data = data.isel(lon=np.sort(unique_indices)).sortby("lon").sortby("lat")
    return data.transpose("time", "lat", "lon")


def decode_argo_juld(juld):
    values = np.asarray(juld.values)
    if np.issubdtype(values.dtype, np.datetime64):
        return pd.to_datetime(values)

    units = str(juld.attrs.get("units", "")).lower()
    candidate_origins = []

    if "since" in units:
        try:
            origin_text = units.split("since", 1)[1].strip()
            candidate_origins.append(pd.Timestamp(origin_text))
        except Exception:
            pass

    candidate_origins.extend([
        pd.Timestamp("1950-01-01"),
        pd.Timestamp("1970-01-01"),
    ])

    best_dates = None
    best_score = -1

    for origin in candidate_origins:
        dates = origin + pd.to_timedelta(values.astype(float), unit="D")
        valid = pd.DatetimeIndex(dates)
        score = int(((valid.year >= 1990) & (valid.year <= 2035)).sum())
        if score > best_score:
            best_score = score
            best_dates = valid

    if best_dates is None:
        raise ValueError("Could not decode Argo JULD values.")

    return pd.DatetimeIndex(best_dates)


def nearest_regular_index(grid, values):
    grid = np.asarray(grid, dtype=float)
    values = np.asarray(values, dtype=float)
    positions = np.searchsorted(grid, values)
    positions = np.clip(positions, 1, len(grid) - 1)
    left = positions - 1
    right = positions
    choose_right = np.abs(values - grid[right]) < np.abs(values - grid[left])
    return np.where(choose_right, right, left).astype(int)


def nearest_circular_longitude_index(grid, values):
    grid = normalize_longitude(grid)
    values = normalize_longitude(values)

    positions = np.searchsorted(grid, values)
    candidate_a = np.mod(positions - 1, len(grid))
    candidate_b = np.mod(positions, len(grid))

    distance_a = np.abs(
        ((values - grid[candidate_a] + 180.0) % 360.0) - 180.0
    )
    distance_b = np.abs(
        ((values - grid[candidate_b] + 180.0) % 360.0) - 180.0
    )

    return np.where(distance_b < distance_a, candidate_b, candidate_a).astype(int)


def haversine_km(lon1, lat1, lon2, lat2):
    radius_km = 6371.0
    lon1 = np.deg2rad(lon1)
    lat1 = np.deg2rad(lat1)
    lon2 = np.deg2rad(lon2)
    lat2 = np.deg2rad(lat2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 2.0 * radius_km * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def target_depth_grid(native_depth):
    native_depth = np.asarray(native_depth, dtype=float)
    inside = native_depth[(native_depth > 0.0) & (native_depth < MAX_DEPTH_M)]
    mandatory = np.array([0.0, MAX_DEPTH_M], dtype=float)
    depth = np.unique(np.round(np.concatenate([inside, mandatory]), 6))
    depth.sort()
    return depth


def interpolate_depth_axis(values, native_depth, target_depth):
    function = interp1d(
        native_depth,
        values,
        axis=0,
        kind="linear",
        bounds_error=False,
        fill_value=np.nan,
        assume_sorted=True,
    )
    output = function(target_depth)
    # Match the EN4 product construction: represent 0 m with first EN4 level.
    zero_index = int(np.argmin(np.abs(target_depth - 0.0)))
    output[zero_index] = values[0]
    return output


def trapezoid(values, x):
    function = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return function(values, x=x)


def integrate_steric(delta_specific_volume, pressure_dbar, gravity):
    return trapezoid(
        delta_specific_volume * 1.0e4 / gravity,
        pressure_dbar,
    )


def scatter_statistics(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 3:
        return {
            "r": np.nan,
            "slope": np.nan,
            "intercept": np.nan,
            "bias": np.nan,
            "rmse": np.nan,
            "n": len(x),
        }

    regression = stats.linregress(x, y)
    return {
        "r": float(np.corrcoef(x, y)[0, 1]),
        "slope": float(regression.slope),
        "intercept": float(regression.intercept),
        "bias": float(np.mean(y - x)),
        "rmse": float(np.sqrt(np.mean((y - x) ** 2))),
        "n": int(len(x)),
    }


def correlation(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 3:
        return np.nan
    return float(np.corrcoef(x[valid], y[valid])[0, 1])


def nice_symmetric_limit(values, minimum=1.0):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return minimum
    maximum = max(float(np.nanpercentile(np.abs(values), 99.5)), minimum)
    steps = np.array([1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50], dtype=float)
    larger = steps[steps >= maximum]
    return float(larger[0] if larger.size else np.ceil(maximum / 10.0) * 10.0)


# ============================================================================
# 7. INPUT DISCOVERY AND LOCAL STAGING
# ============================================================================
ARGO_SOURCE = first_existing(ARGO_CANDIDATES, "Argo profile")
EN4_RAW_SOURCE = first_existing(EN4_RAW_CANDIDATES, "raw EN4")
require_file(EN4_TOTAL_FILE, "EN4 total steric")
require_file(EN4_THERMO_FILE, "EN4 thermosteric")
require_file(EN4_HALO_FILE, "EN4 halosteric")

log("=" * 88)
log("FIGURE S3 — ARGO–EN4 0–1000 m VALIDATION")
log("=" * 88)
log(f"Argo source      : {ARGO_SOURCE}")
log(f"EN4 raw source   : {EN4_RAW_SOURCE}")
log(f"EN4 total product: {EN4_TOTAL_FILE}")
log(f"EN4 thermo       : {EN4_THERMO_FILE}")
log(f"EN4 halo         : {EN4_HALO_FILE}")

ARGO_FILE = stage_file_to_local(ARGO_SOURCE, "Argo file")
EN4_RAW_FILE = stage_file_to_local(EN4_RAW_SOURCE, "raw EN4 file")
EN4_TOTAL_LOCAL = stage_file_to_local(EN4_TOTAL_FILE, "EN4 total product")
EN4_THERMO_LOCAL = stage_file_to_local(EN4_THERMO_FILE, "EN4 thermosteric product")
EN4_HALO_LOCAL = stage_file_to_local(EN4_HALO_FILE, "EN4 halosteric product")


# 8. READ EN4 MONTHLY PRODUCTS
def load_en4_product(path, variable_name):
    dataset = open_dataset_robust(path, decode_times=True)
    try:
        data = standardize_3d_product(dataset, variable_name).load()
        units = str(data.attrs.get("units", "")).lower()
        if units in {"m", "metre", "meter", "metres", "meters"}:
            data = data * 100.0
        elif "cm" not in units:
            # The verified products are metres; this fallback preserves that behavior.
            data = data * 100.0
        data.attrs["units"] = "cm"
        return data
    finally:
        dataset.close()
        gc.collect()


log("\nLoading compact EN4 monthly steric products ...")
en4_total = load_en4_product(EN4_TOTAL_LOCAL, "total_steric_height")
en4_thermo = load_en4_product(EN4_THERMO_LOCAL, "thermosteric_height")
en4_halo = load_en4_product(EN4_HALO_LOCAL, "halosteric_height")

# Ensure exact common coordinates.
en4_thermo = en4_thermo.sel(time=en4_total.time, lat=en4_total.lat, lon=en4_total.lon)
en4_halo = en4_halo.sel(time=en4_total.time, lat=en4_total.lat, lon=en4_total.lon)

en4_times = pd.DatetimeIndex(pd.to_datetime(en4_total.time.values))
en4_latitudes = en4_total.lat.values.astype(float)
en4_longitudes = en4_total.lon.values.astype(float)
month_to_time_index = {
    pd.Timestamp(date).to_period("M"): index
    for index, date in enumerate(en4_times)
}

log(
    f"EN4 product grid: time={len(en4_times)}, "
    f"lat={len(en4_latitudes)}, lon={len(en4_longitudes)}"
)


# 9. BUILD OR LOAD FIXED EN4 REFERENCE PROFILES
def build_reference_cache():
    source_size = EN4_RAW_SOURCE.stat().st_size

    if REFERENCE_CACHE.exists() and REUSE_REFERENCE_CACHE and not FORCE_RECOMPUTE:
        cached = np.load(REFERENCE_CACHE, allow_pickle=False)
        cached_size = int(cached["source_size_bytes"])
        if cached_size == source_size:
            log(f"\nLoading cached EN4 reference: {REFERENCE_CACHE.name}")
            return {
                "depth": cached["depth"].astype(float),
                "temperature_c": cached["temperature_c"].astype(np.float32),
                "salinity": cached["salinity"].astype(np.float32),
                "lat": cached["lat"].astype(float),
                "lon": cached["lon"].astype(float),
            }
        log("Reference cache source-size mismatch; recalculating.")

    log("\nCalculating fixed 2008–2025 EN4 reference profiles with low RAM ...")

    dataset = open_dataset_robust(
        EN4_RAW_FILE,
        decode_times=False,
        chunks={"time": 12},
    )

    try:
        rename = {}
        for source, target in (
            ("latitude", "lat"),
            ("longitude", "lon"),
        ):
            if source in dataset.coords or source in dataset.dims:
                if target not in dataset.coords and target not in dataset.dims:
                    rename[source] = target
        if rename:
            dataset = dataset.rename(rename)

        native_depth = dataset["depth"].values.astype(float)
        analysis_depth = target_depth_grid(native_depth)

        with dask.config.set(scheduler="single-threaded"):
            temperature_native = (
                dataset["temperature"]
                .mean("time", skipna=True)
                .compute()
                .values
                .astype(np.float32)
            )
            salinity_native = (
                dataset["salinity"]
                .mean("time", skipna=True)
                .compute()
                .values
                .astype(np.float32)
            )

        # EN4 temperature is potential temperature in Kelvin.
        if np.nanmedian(temperature_native) > 100.0:
            temperature_native = temperature_native - 273.15

        temperature_reference = interpolate_depth_axis(
            temperature_native,
            native_depth,
            analysis_depth,
        ).astype(np.float32)

        salinity_reference = interpolate_depth_axis(
            salinity_native,
            native_depth,
            analysis_depth,
        ).astype(np.float32)

        latitudes = dataset["lat"].values.astype(float)
        longitudes = normalize_longitude(dataset["lon"].values.astype(float))

        # Sort longitude and remove duplicates if necessary.
        sort_indices = np.argsort(longitudes)
        longitudes = longitudes[sort_indices]
        temperature_reference = temperature_reference[:, :, sort_indices]
        salinity_reference = salinity_reference[:, :, sort_indices]
        unique_longitudes, unique_indices = np.unique(longitudes, return_index=True)
        longitudes = unique_longitudes
        temperature_reference = temperature_reference[:, :, unique_indices]
        salinity_reference = salinity_reference[:, :, unique_indices]

        np.savez_compressed(
            REFERENCE_CACHE,
            depth=analysis_depth.astype(np.float32),
            temperature_c=temperature_reference,
            salinity=salinity_reference,
            lat=latitudes.astype(np.float32),
            lon=longitudes.astype(np.float32),
            source_size_bytes=np.array(source_size, dtype=np.int64),
        )

        log(f"Saved EN4 reference cache: {REFERENCE_CACHE}")

        return {
            "depth": analysis_depth,
            "temperature_c": temperature_reference,
            "salinity": salinity_reference,
            "lat": latitudes,
            "lon": longitudes,
        }

    finally:
        dataset.close()
        gc.collect()


reference = build_reference_cache()
analysis_depth = reference["depth"]
reference_temperature = reference["temperature_c"]
reference_salinity = reference["salinity"]

# Verify reference and product grids agree.
if not np.allclose(reference["lat"], en4_latitudes, atol=1e-5):
    raise ValueError("Raw EN4 and derived EN4 latitude grids do not agree.")
if not np.allclose(reference["lon"], en4_longitudes, atol=1e-5):
    raise ValueError("Raw EN4 and derived EN4 longitude grids do not agree.")

log(f"Reference depth levels: {len(analysis_depth)}; {analysis_depth[0]:.1f}–{analysis_depth[-1]:.1f} m")


# 10. ARGO CHUNK PROCESSING
def process_argo_chunks():
    if MATCHED_CSV.exists() and not FORCE_RECOMPUTE:
        log(f"\nLoading completed matched table: {MATCHED_CSV.name}")
        return pd.read_csv(MATCHED_CSV, parse_dates=["time"])

    dataset = open_dataset_robust(ARGO_FILE, decode_times=False)

    try:
        required_variables = [
            "PRES_GRID", "JULD", "LATITUDE", "LONGITUDE", "TEMP", "PSAL"
        ]
        missing = [name for name in required_variables if name not in dataset.variables]
        if missing:
            raise KeyError(
                f"Argo file is missing variables {missing}. "
                f"Available variables: {list(dataset.variables)}"
            )

        pressure_grid = dataset["PRES_GRID"].values.astype(float)
        if np.any(np.diff(pressure_grid) <= 0):
            raise ValueError("PRES_GRID must be strictly increasing.")

        argo_times = decode_argo_juld(dataset["JULD"])
        latitudes_all = dataset["LATITUDE"].values.astype(float)
        longitudes_all = normalize_longitude(dataset["LONGITUDE"].values.astype(float))

        n_profiles = int(dataset.sizes["N_PROF"])
        log(f"\nArgo profiles in source: {n_profiles:,}")
        log(f"Processing chunk size: {ARGO_PROFILE_CHUNK:,}")

        all_chunk_files = []

        for start in range(0, n_profiles, ARGO_PROFILE_CHUNK):
            stop = min(start + ARGO_PROFILE_CHUNK, n_profiles)
            chunk_file = CHUNK_CACHE_DIR / f"matched_{start:06d}_{stop:06d}.csv"
            all_chunk_files.append(chunk_file)

            if chunk_file.exists() and REUSE_CHUNK_CACHE and not FORCE_RECOMPUTE:
                log(f"Reusing chunk {start:06d}:{stop:06d}")
                continue

            chunk_times = argo_times[start:stop]
            chunk_latitudes = latitudes_all[start:stop]
            chunk_longitudes = longitudes_all[start:stop]

            preliminary = (
                (chunk_times >= START_DATE)
                & (chunk_times <= END_DATE)
                & np.isfinite(chunk_latitudes)
                & np.isfinite(chunk_longitudes)
                & (chunk_latitudes >= -90.0)
                & (chunk_latitudes <= -50.0)
            )

            if preliminary.sum() == 0:
                pd.DataFrame().to_csv(chunk_file, index=False)
                log(f"[{stop:>6}/{n_profiles}] no profiles in analysis period/domain")
                continue

            # Read only this contiguous block from the large HDF5 arrays.
            temperature_chunk = (
                dataset["TEMP"]
                .isel(N_PROF=slice(start, stop))
                .load()
                .values
                .astype(np.float32)
            )
            salinity_chunk = (
                dataset["PSAL"]
                .isel(N_PROF=slice(start, stop))
                .load()
                .values
                .astype(np.float32)
            )

            if np.nanmedian(temperature_chunk) > 100.0:
                temperature_chunk = temperature_chunk - 273.15

            rows = []
            local_indices = np.where(preliminary)[0]

            for local_index in local_indices:
                profile_time = pd.Timestamp(chunk_times[local_index])
                profile_latitude = float(chunk_latitudes[local_index])
                profile_longitude = float(chunk_longitudes[local_index])
                sea = assign_sea(profile_longitude)
                if sea is None:
                    continue

                month_period = profile_time.to_period("M")
                if month_period not in month_to_time_index:
                    continue
                time_index = month_to_time_index[month_period]

                lat_index = int(nearest_regular_index(en4_latitudes, [profile_latitude])[0])
                lon_index = int(nearest_circular_longitude_index(en4_longitudes, [profile_longitude])[0])

                en4_latitude = float(en4_latitudes[lat_index])
                en4_longitude = float(en4_longitudes[lon_index])
                match_distance = float(
                    haversine_km(
                        profile_longitude,
                        profile_latitude,
                        en4_longitude,
                        en4_latitude,
                    )
                )
                if match_distance > MAX_NEAREST_CELL_DISTANCE_KM:
                    continue

                en4_total_value = float(en4_total.values[time_index, lat_index, lon_index])
                en4_thermo_value = float(en4_thermo.values[time_index, lat_index, lon_index])
                en4_halo_value = float(en4_halo.values[time_index, lat_index, lon_index])

                if not np.all(np.isfinite([
                    en4_total_value,
                    en4_thermo_value,
                    en4_halo_value,
                ])):
                    continue

                temperature_native = temperature_chunk[local_index].astype(float)
                salinity_native = salinity_chunk[local_index].astype(float)

                common = (
                    np.isfinite(pressure_grid)
                    & np.isfinite(temperature_native)
                    & np.isfinite(salinity_native)
                    & (temperature_native >= TEMP_RANGE_C[0])
                    & (temperature_native <= TEMP_RANGE_C[1])
                    & (salinity_native >= SAL_RANGE[0])
                    & (salinity_native <= SAL_RANGE[1])
                )

                if common.sum() < 4:
                    continue

                pressure_valid = pressure_grid[common]
                temperature_valid = temperature_native[common]
                salinity_valid = salinity_native[common]

                # Convert pressure to positive depth at the profile latitude.
                depth_valid = -gsw.z_from_p(pressure_valid, profile_latitude)
                ordering = np.argsort(depth_valid)
                depth_valid = depth_valid[ordering]
                temperature_valid = temperature_valid[ordering]
                salinity_valid = salinity_valid[ordering]

                unique_depth, unique_indices = np.unique(depth_valid, return_index=True)
                depth_valid = unique_depth
                temperature_valid = temperature_valid[unique_indices]
                salinity_valid = salinity_valid[unique_indices]

                if depth_valid[0] > MAX_ALLOWED_TOP_DEPTH_M:
                    continue
                if depth_valid[-1] < MIN_ALLOWED_BOTTOM_DEPTH_M:
                    continue
                if np.nanmax(np.diff(depth_valid)) > MAX_VERTICAL_GAP_M:
                    continue

                temperature_function = interp1d(
                    depth_valid,
                    temperature_valid,
                    kind="linear",
                    bounds_error=False,
                    fill_value=np.nan,
                    assume_sorted=True,
                )
                salinity_function = interp1d(
                    depth_valid,
                    salinity_valid,
                    kind="linear",
                    bounds_error=False,
                    fill_value=np.nan,
                    assume_sorted=True,
                )

                profile_temperature = temperature_function(analysis_depth)
                profile_salinity = salinity_function(analysis_depth)

                # Match the EN4 construction at 0 m: use the shallowest valid value.
                profile_temperature[0] = temperature_valid[0]
                profile_salinity[0] = salinity_valid[0]

                if not (
                    np.isfinite(profile_temperature).all()
                    and np.isfinite(profile_salinity).all()
                ):
                    continue

                reference_temperature_profile = reference_temperature[:, lat_index, lon_index].astype(float)
                reference_salinity_profile = reference_salinity[:, lat_index, lon_index].astype(float)
                if not (
                    np.isfinite(reference_temperature_profile).all()
                    and np.isfinite(reference_salinity_profile).all()
                ):
                    continue

                # Use the matched EN4 grid-cell pressure/gravity to reproduce
                # the EN4 product integration as directly as possible.
                pressure_target = gsw.p_from_z(-analysis_depth, en4_latitude)
                gravity_target = gsw.grav(en4_latitude, pressure_target)

                absolute_salinity_reference = gsw.SA_from_SP(
                    reference_salinity_profile,
                    pressure_target,
                    en4_longitude,
                    en4_latitude,
                )
                conservative_temperature_reference = gsw.CT_from_pt(
                    absolute_salinity_reference,
                    reference_temperature_profile,
                )
                specific_volume_reference = gsw.specvol(
                    absolute_salinity_reference,
                    conservative_temperature_reference,
                    pressure_target,
                )

                absolute_salinity_profile = gsw.SA_from_SP(
                    profile_salinity,
                    pressure_target,
                    en4_longitude,
                    en4_latitude,
                )
                conservative_temperature_profile = gsw.CT_from_t(
                    absolute_salinity_profile,
                    profile_temperature,
                    pressure_target,
                )
                potential_temperature_profile = gsw.pt0_from_t(
                    absolute_salinity_profile,
                    profile_temperature,
                    pressure_target,
                )
                specific_volume_profile = gsw.specvol(
                    absolute_salinity_profile,
                    conservative_temperature_profile,
                    pressure_target,
                )

                argo_total_cm = 100.0 * integrate_steric(
                    specific_volume_profile - specific_volume_reference,
                    pressure_target,
                    gravity_target,
                )

                conservative_temperature_thermo = gsw.CT_from_pt(
                    absolute_salinity_reference,
                    potential_temperature_profile,
                )
                specific_volume_thermo = gsw.specvol(
                    absolute_salinity_reference,
                    conservative_temperature_thermo,
                    pressure_target,
                )

                argo_thermo_cm = 100.0 * integrate_steric(
                    specific_volume_thermo - specific_volume_reference,
                    pressure_target,
                    gravity_target,
                )
                argo_halo_cm = argo_total_cm - argo_thermo_cm

                if not np.all(np.isfinite([
                    argo_total_cm,
                    argo_thermo_cm,
                    argo_halo_cm,
                ])):
                    continue

                rows.append({
                    "profile_index": int(start + local_index),
                    "time": profile_time,
                    "season": month_to_season(profile_time.month),
                    "sea": sea,
                    "latitude": profile_latitude,
                    "longitude": profile_longitude,
                    "en4_latitude": en4_latitude,
                    "en4_longitude": en4_longitude,
                    "match_distance_km": match_distance,
                    "argo_total_cm": float(argo_total_cm),
                    "en4_total_cm": en4_total_value,
                    "argo_thermo_cm": float(argo_thermo_cm),
                    "en4_thermo_cm": en4_thermo_value,
                    "argo_halo_cm": float(argo_halo_cm),
                    "en4_halo_cm": en4_halo_value,
                })

            pd.DataFrame(rows).to_csv(chunk_file, index=False)
            log(
                f"[{stop:>6}/{n_profiles}] saved {len(rows):>4} matches: "
                f"{chunk_file.name}"
            )

            del temperature_chunk, salinity_chunk
            gc.collect()

        log("\nCombining cached matching chunks ...")
        tables = []
        for chunk_file in all_chunk_files:
            if not chunk_file.exists() or chunk_file.stat().st_size == 0:
                continue
            try:
                table = pd.read_csv(chunk_file, parse_dates=["time"])
            except pd.errors.EmptyDataError:
                continue
            if not table.empty:
                tables.append(table)

        if not tables:
            raise RuntimeError(
                "No valid Argo–EN4 profile matches were produced. "
                "Review the pressure coverage and variable ranges."
            )

        matched_table = pd.concat(tables, ignore_index=True)
        matched_table = matched_table.drop_duplicates("profile_index").sort_values("time")
        matched_table.to_csv(MATCHED_CSV, index=False)
        log(f"Saved completed matched table: {MATCHED_CSV}")
        log(f"Matched profiles: {len(matched_table):,}")
        return matched_table

    finally:
        dataset.close()
        gc.collect()


matched = process_argo_chunks()

# Free the reference and monthly product arrays before high-DPI plotting.
del reference_temperature, reference_salinity, reference
try:
    del en4_total, en4_thermo, en4_halo
except Exception:
    pass
gc.collect()


# 11. SUMMARY TABLES
component_definitions = OrderedDict([
    ("total", {
        "argo": "argo_total_cm",
        "en4": "en4_total_cm",
        "label": "Total steric",
        "color": "#1665C1",
    }),
    ("thermo", {
        "argo": "argo_thermo_cm",
        "en4": "en4_thermo_cm",
        "label": "Thermosteric",
        "color": "#229A22",
    }),
    ("halo", {
        "argo": "argo_halo_cm",
        "en4": "en4_halo_cm",
        "label": "Halosteric",
        "color": "#9227A8",
    }),
])

scatter_summary = {
    key: scatter_statistics(matched[info["en4"]], matched[info["argo"]])
    for key, info in component_definitions.items()
}

sea_bias_rows = []
for sea in SEA_ORDER:
    subset = matched[matched["sea"] == sea]
    row = {"sea": sea, "matched_profiles": int(len(subset))}
    for key, info in component_definitions.items():
        row[f"bias_{key}_cm"] = (
            float(np.mean(subset[info["argo"]] - subset[info["en4"]]))
            if len(subset) else np.nan
        )
    sea_bias_rows.append(row)

sea_bias = pd.DataFrame(sea_bias_rows)
sea_bias.to_csv(SEA_BIAS_CSV, index=False)

season_rows = []
for season in SEASON_ORDER:
    subset = matched[matched["season"] == season]
    row = {"season": season, "matched_profiles": int(len(subset))}
    for key, info in component_definitions.items():
        row[f"r_{key}"] = correlation(subset[info["en4"]], subset[info["argo"]])
    season_rows.append(row)

seasonal_correlation = pd.DataFrame(season_rows)
seasonal_correlation.to_csv(SEASONAL_CORR_CSV, index=False)

sea_counts = (
    matched.groupby("sea")
    .size()
    .reindex(SEA_ORDER, fill_value=0)
    .astype(int)
    .rename("matched_profiles")
    .reset_index()
)
sea_counts.to_csv(SEA_COUNT_CSV, index=False)



# PART B — LATEST CORRECTED S3 PLOT-ONLY WORKFLOW
# Reads the precomputed CSVs from Part A. No science is recomputed here.

# Figure S3 — Argo–EN4 0–1000 m validation
# PLOT-ONLY VERSION USING THE EXISTING PRECOMPUTED CSV FILES
#
# Final correction:
#   - Panel (d) legend row moved LOWER so it does not overlap
#     the x-axis label.
#   - Panel (e) seasonal labels remain visible.
#   - No recomputation is done.
#   - Output filenames remain unchanged.
# 1. Paths — keep filenames unchanged
OUTPUT_DIR = Path("/content/drive/MyDrive/SAM_Thesis/paper2")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MATCHED_CSV = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_matched_profiles.csv"
SEA_BIAS_CSV = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_sea_bias.csv"
SEASONAL_CORR_CSV = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_seasonal_correlations.csv"
SEA_COUNTS_CSV = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_matched_counts.csv"

# SAME output names as before
OUTPUT_PNG = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation_1080dpi.png"
OUTPUT_PDF = OUTPUT_DIR / "Figure_S3_Argo_EN4_0_1000_validation.pdf"


# 2. Check input files
required_files = [MATCHED_CSV, SEA_BIAS_CSV, SEASONAL_CORR_CSV, SEA_COUNTS_CSV]
missing_files = [path for path in required_files if not path.exists()]

if missing_files:
    raise FileNotFoundError(
        "The following precomputed Figure S3 files were not found:\n"
        + "\n".join(str(path) for path in missing_files)
    )

print("Using precomputed files:")
print("Matched CSV             :", MATCHED_CSV)
print("Sea bias CSV            :", SEA_BIAS_CSV)
print("Seasonal correlation CSV:", SEASONAL_CORR_CSV)
print("Sea counts CSV          :", SEA_COUNTS_CSV)


# 3. Read the exact existing CSV structures
matched = pd.read_csv(MATCHED_CSV, parse_dates=["time"])
sea_bias = pd.read_csv(SEA_BIAS_CSV)
seasonal_corr = pd.read_csv(SEASONAL_CORR_CSV)
sea_counts = pd.read_csv(SEA_COUNTS_CSV)

required_matched_columns = [
    "sea", "season",
    "argo_total_cm", "en4_total_cm",
    "argo_thermo_cm", "en4_thermo_cm",
    "argo_halo_cm", "en4_halo_cm",
]

required_bias_columns = [
    "sea", "matched_profiles",
    "bias_total_cm", "bias_thermo_cm", "bias_halo_cm",
]

required_season_columns = [
    "season", "matched_profiles",
    "r_total", "r_thermo", "r_halo",
]

required_count_columns = [
    "sea", "matched_profiles",
]

for column in required_matched_columns:
    if column not in matched.columns:
        raise KeyError(
            f"Missing column '{column}' in:\n{MATCHED_CSV}\n"
            f"Available columns: {list(matched.columns)}"
        )

for column in required_bias_columns:
    if column not in sea_bias.columns:
        raise KeyError(
            f"Missing column '{column}' in:\n{SEA_BIAS_CSV}\n"
            f"Available columns: {list(sea_bias.columns)}"
        )

for column in required_season_columns:
    if column not in seasonal_corr.columns:
        raise KeyError(
            f"Missing column '{column}' in:\n{SEASONAL_CORR_CSV}\n"
            f"Available columns: {list(seasonal_corr.columns)}"
        )

for column in required_count_columns:
    if column not in sea_counts.columns:
        raise KeyError(
            f"Missing column '{column}' in:\n{SEA_COUNTS_CSV}\n"
            f"Available columns: {list(sea_counts.columns)}"
        )


# 4. Fixed ordering and colours
SEA_ORDER = [
    "WED", "KHV", "RLS", "LAZ", "COS", "COO", "DAV",
    "MAW", "DUR", "SOM", "ROS", "AMU", "BEL"
]

SEASON_ORDER = ["Spring", "Summer", "Autumn", "Winter"]

COLOR_TOTAL = "#1769E0"
COLOR_THERMO = "#20A226"
COLOR_HALO = "#9C27B0"
COLOR_COUNTS = "#72A6D9"

sea_bias = sea_bias.set_index("sea").reindex(SEA_ORDER).reset_index()
seasonal_corr = seasonal_corr.set_index("season").reindex(SEASON_ORDER).reset_index()
sea_counts = sea_counts.set_index("sea").reindex(SEA_ORDER).reset_index()


# 5. Plot settings
SAVE_DPI = 1080
FIGSIZE = (16.5, 10.8)

FONT_PANEL_TITLE = 15.0
FONT_AXIS_LABEL = 12.5
FONT_TICK = 10.0
FONT_LEGEND = 9.5
FONT_STATS = 10.3
FONT_ANNOTATION = 8.8

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "axes.linewidth": 0.9,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


# 6. Helper functions
def calculate_scatter_statistics(x_values, y_values):
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 2:
        return {
            "x": x, "y": y, "n": len(x),
            "r": np.nan, "slope": np.nan, "intercept": np.nan,
            "bias": np.nan, "rmse": np.nan,
        }

    slope, intercept = np.polyfit(x, y, 1)
    correlation = np.corrcoef(x, y)[0, 1]
    bias = np.mean(y - x)
    rmse = np.sqrt(np.mean((y - x) ** 2))

    return {
        "x": x, "y": y, "n": len(x),
        "r": correlation, "slope": slope, "intercept": intercept,
        "bias": bias, "rmse": rmse,
    }


def common_scatter_limits(x_values, y_values):
    values = np.concatenate([
        np.asarray(x_values, dtype=float),
        np.asarray(y_values, dtype=float),
    ])
    values = values[np.isfinite(values)]

    if values.size == 0:
        return -1.0, 1.0

    lower = float(np.nanmin(values))
    upper = float(np.nanmax(values))

    span = upper - lower
    if span <= 0:
        span = max(abs(lower), 1.0)

    padding = max(0.5, 0.07 * span)
    lower -= padding
    upper += padding

    absolute_limit = max(abs(lower), abs(upper))
    return -absolute_limit, absolute_limit


def format_axis(axis):
    axis.grid(True, linestyle=":", linewidth=0.55, alpha=0.32, zorder=0)
    axis.tick_params(axis="both", labelsize=FONT_TICK, pad=3.0)

    for label in axis.get_xticklabels():
        label.set_fontweight("bold")
    for label in axis.get_yticklabels():
        label.set_fontweight("bold")


def plot_validation_scatter(axis, en4_values, argo_values, colour, title, x_label, y_label):
    statistics = calculate_scatter_statistics(en4_values, argo_values)
    lower, upper = common_scatter_limits(statistics["x"], statistics["y"])
    comparison_line = np.linspace(lower, upper, 300)

    axis.scatter(
        statistics["x"], statistics["y"],
        s=13, facecolor=colour, edgecolor="none", alpha=0.65,
        rasterized=True, zorder=2,
    )

    axis.plot(
        comparison_line, comparison_line,
        color="0.25", linestyle="--", linewidth=1.1, dashes=(4, 3),
        label="1:1 line", zorder=3,
    )

    if np.isfinite(statistics["slope"]):
        regression_line = statistics["intercept"] + statistics["slope"] * comparison_line
        axis.plot(
            comparison_line, regression_line,
            color=colour, linewidth=1.8, label="Regression", zorder=4,
        )

    axis.set_xlim(lower, upper)
    axis.set_ylim(lower, upper)
    axis.set_aspect("equal", adjustable="box")

    axis.set_title(title, fontsize=FONT_PANEL_TITLE, fontweight="bold", pad=8.0)
    axis.set_xlabel(x_label, fontsize=FONT_AXIS_LABEL, fontweight="bold")
    axis.set_ylabel(y_label, fontsize=FONT_AXIS_LABEL, fontweight="bold")

    statistics_text = (
        f"r = {statistics['r']:.2f}\n"
        f"slope = {statistics['slope']:.2f}\n"
        f"bias = {statistics['bias']:.2f} cm\n"
        f"RMSE = {statistics['rmse']:.2f} cm\n"
        f"n = {statistics['n']}"
    )

    axis.text(
        0.04, 0.95, statistics_text,
        transform=axis.transAxes,
        ha="left", va="top",
        fontsize=FONT_STATS, fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="white",
            edgecolor="0.30",
            alpha=0.94,
        ),
        zorder=8,
    )

    axis.legend(
        loc="lower right",
        fontsize=FONT_LEGEND,
        frameon=True,
        facecolor="white",
        framealpha=0.95,
        borderpad=0.35,
        handlelength=1.7,
        labelspacing=0.25,
    )

    format_axis(axis)


# 7. Create the figure
figure = plt.figure(figsize=FIGSIZE, constrained_layout=False)

main_grid = figure.add_gridspec(
    2, 3,
    left=0.060,
    right=0.985,
    top=0.965,
    bottom=0.080,
    width_ratios=[1.0, 1.0, 1.0],
    height_ratios=[1.0, 1.10],
    wspace=0.24,
    hspace=0.35,
)

axis_a = figure.add_subplot(main_grid[0, 0])
axis_b = figure.add_subplot(main_grid[0, 1])
axis_c = figure.add_subplot(main_grid[0, 2])

# PANEL (d): extra space + lower legend row
panel_d_grid = main_grid[1, 0].subgridspec(
    2, 1,
    height_ratios=[17.4, 3.1],   # more height for legend row
    hspace=0.22,                 # larger gap between x-label and legend row
)

axis_d = figure.add_subplot(panel_d_grid[0, 0])
axis_d_legend = figure.add_subplot(panel_d_grid[1, 0])
axis_d_legend.set_axis_off()

# PANEL (e): keep separate legend row
panel_e_grid = main_grid[1, 1].subgridspec(
    2, 1,
    height_ratios=[18.0, 2.4],
    hspace=0.15,
)

axis_e = figure.add_subplot(panel_e_grid[0, 0])
axis_e_legend = figure.add_subplot(panel_e_grid[1, 0])
axis_e_legend.set_axis_off()

axis_f = figure.add_subplot(main_grid[1, 2])


# 8. Panels (a)–(c)
plot_validation_scatter(
    axis_a,
    matched["en4_total_cm"],
    matched["argo_total_cm"],
    COLOR_TOTAL,
    "(a) Argo versus EN4 total steric",
    "EN4 total steric (cm)",
    "Argo total steric (cm)",
)

plot_validation_scatter(
    axis_b,
    matched["en4_thermo_cm"],
    matched["argo_thermo_cm"],
    COLOR_THERMO,
    "(b) Argo versus EN4 thermosteric",
    "EN4 thermosteric (cm)",
    "Argo thermosteric (cm)",
)

plot_validation_scatter(
    axis_c,
    matched["en4_halo_cm"],
    matched["argo_halo_cm"],
    COLOR_HALO,
    "(c) Argo versus EN4 halosteric",
    "EN4 halosteric (cm)",
    "Argo halosteric (cm)",
)


# 9. Panel (d): sea-wise bias
sea_positions = np.arange(len(SEA_ORDER))
bar_height = 0.22

bars_d_total = axis_d.barh(
    sea_positions - bar_height,
    sea_bias["bias_total_cm"].to_numpy(dtype=float),
    height=bar_height,
    color=COLOR_TOTAL,
    label="Total steric",
    zorder=3,
)

bars_d_thermo = axis_d.barh(
    sea_positions,
    sea_bias["bias_thermo_cm"].to_numpy(dtype=float),
    height=bar_height,
    color=COLOR_THERMO,
    label="Thermosteric",
    zorder=3,
)

bars_d_halo = axis_d.barh(
    sea_positions + bar_height,
    sea_bias["bias_halo_cm"].to_numpy(dtype=float),
    height=bar_height,
    color=COLOR_HALO,
    label="Halosteric",
    zorder=3,
)

axis_d.axvline(0.0, color="0.20", linewidth=0.9, zorder=2)

axis_d.set_yticks(sea_positions)
axis_d.set_yticklabels(SEA_ORDER, fontsize=FONT_TICK, fontweight="bold")
axis_d.invert_yaxis()

all_biases = np.concatenate([
    sea_bias["bias_total_cm"].to_numpy(dtype=float),
    sea_bias["bias_thermo_cm"].to_numpy(dtype=float),
    sea_bias["bias_halo_cm"].to_numpy(dtype=float),
])

finite_biases = all_biases[np.isfinite(all_biases)]
if finite_biases.size:
    bias_limit = max(1.0, float(np.nanmax(np.abs(finite_biases))) * 1.18)
else:
    bias_limit = 1.0

axis_d.set_xlim(-bias_limit, bias_limit)

axis_d.set_title("(d) Sea-wise bias", fontsize=FONT_PANEL_TITLE, fontweight="bold", pad=8.0)
axis_d.set_xlabel("Mean bias: Argo − EN4 (cm)", fontsize=FONT_AXIS_LABEL, fontweight="bold", labelpad=6)

format_axis(axis_d)

# Moved LOWER so it no longer overlaps the x-axis label
axis_d_legend.legend(
    handles=[bars_d_total, bars_d_thermo, bars_d_halo],
    labels=["Total steric", "Thermosteric", "Halosteric"],
    loc="lower left",
    bbox_to_anchor=(0.01, 0.02),
    ncol=3,
    fontsize=FONT_LEGEND,
    frameon=True,
    facecolor="white",
    framealpha=0.96,
    borderpad=0.35,
    columnspacing=1.0,
    handlelength=1.3,
)


# 10. Panel (e): seasonal correlation
season_positions = np.arange(len(SEASON_ORDER))
bar_width = 0.23

bars_e_total = axis_e.bar(
    season_positions - bar_width,
    seasonal_corr["r_total"].to_numpy(dtype=float),
    width=bar_width,
    color=COLOR_TOTAL,
    label="Total steric",
    zorder=3,
)

bars_e_thermo = axis_e.bar(
    season_positions,
    seasonal_corr["r_thermo"].to_numpy(dtype=float),
    width=bar_width,
    color=COLOR_THERMO,
    label="Thermosteric",
    zorder=3,
)

bars_e_halo = axis_e.bar(
    season_positions + bar_width,
    seasonal_corr["r_halo"].to_numpy(dtype=float),
    width=bar_width,
    color=COLOR_HALO,
    label="Halosteric",
    zorder=3,
)

axis_e.set_xticks(season_positions)
axis_e.set_xticklabels(SEASON_ORDER, fontsize=FONT_TICK, fontweight="bold")
axis_e.tick_params(axis="x", labelsize=FONT_TICK, pad=5.0)

axis_e.set_ylim(0.0, 1.04)
axis_e.set_ylabel("Pearson correlation (r)", fontsize=FONT_AXIS_LABEL, fontweight="bold")
axis_e.set_title("(e) Seasonal correlation", fontsize=FONT_PANEL_TITLE, fontweight="bold", pad=8.0)

for bar_group in [bars_e_total, bars_e_thermo, bars_e_halo]:
    for bar in bar_group:
        value = float(bar.get_height())
        if not np.isfinite(value):
            continue
        axis_e.text(
            bar.get_x() + bar.get_width() / 2.0,
            value + 0.018,
            f"{value:.2f}",
            ha="center",
            va="bottom",
            fontsize=FONT_ANNOTATION,
            fontweight="bold",
            clip_on=False,
            zorder=5,
        )

format_axis(axis_e)

axis_e_legend.legend(
    handles=[bars_e_total, bars_e_thermo, bars_e_halo],
    labels=["Total steric", "Thermosteric", "Halosteric"],
    loc="center",
    ncol=3,
    fontsize=FONT_LEGEND,
    frameon=True,
    facecolor="white",
    framealpha=0.96,
    borderpad=0.35,
    columnspacing=1.2,
    handlelength=1.4,
)


# 11. Panel (f): number of matched profiles
count_positions = np.arange(len(SEA_ORDER))
count_values = sea_counts["matched_profiles"].fillna(0).to_numpy(dtype=int)

bars_f = axis_f.bar(
    count_positions,
    count_values,
    width=0.68,
    color=COLOR_COUNTS,
    edgecolor="#3F6F9D",
    linewidth=0.7,
    zorder=3,
)

axis_f.set_xticks(count_positions)
axis_f.set_xticklabels(SEA_ORDER, rotation=45, ha="right", fontsize=FONT_TICK, fontweight="bold")
axis_f.set_ylabel("Matched profiles (n)", fontsize=FONT_AXIS_LABEL, fontweight="bold")
axis_f.set_title("(f) Number of matched profiles", fontsize=FONT_PANEL_TITLE, fontweight="bold", pad=8.0)

maximum_count = max(int(np.nanmax(count_values)), 1)
axis_f.set_ylim(0.0, maximum_count * 1.19)

for bar, value in zip(bars_f, count_values):
    axis_f.text(
        bar.get_x() + bar.get_width() / 2.0,
        value + maximum_count * 0.018,
        f"{int(value)}",
        ha="center",
        va="bottom",
        fontsize=FONT_ANNOTATION,
        fontweight="bold",
        zorder=5,
    )

format_axis(axis_f)


# 12. Save — filenames unchanged
print("Saving corrected plot-only Figure S3...")

figure.savefig(
    OUTPUT_PDF,
    facecolor="white",
    edgecolor="none",
    bbox_inches=None,
)

figure.savefig(
    OUTPUT_PNG,
    dpi=SAVE_DPI,
    facecolor="white",
    edgecolor="none",
    bbox_inches=None,
    pil_kwargs={"compress_level": 6},
)

plt.show()
plt.close(figure)

print("\nCompleted successfully.")
print("PNG:", OUTPUT_PNG)
print("PDF:", OUTPUT_PDF)
# FINAL SCIENTIFIC PROCESSING LOG

log_lines = [
    "=" * 88,
    "FIGURE S3 — ARGO–EN4 0–1000 m VALIDATION",
    "=" * 88,
    "",
    f"Argo source: {ARGO_SOURCE}",
    f"EN4 raw source: {EN4_RAW_SOURCE}",
    f"Matched profiles: {len(matched):,}",
    "",
    "Matching method:",
    "• Same calendar month.",
    "• Nearest EN4 grid cell, maximum distance 120 km.",
    "• Argo T/S interpolated to the EN4 0–1000 m analysis-depth grid.",
    "• Complete common T/S coverage from the near-surface to 1000 m.",
    "• Fixed 2008–2025 EN4 grid-cell reference profiles.",
    "• TEOS-10 total, thermosteric and residual halosteric decomposition.",
    "",
    "Reliability corrections:",
    "• Google-Drive inputs staged to /content before HDF5 access.",
    "• h5netcdf attempted before netCDF4.",
    "• Argo processed in contiguous low-RAM chunks.",
    "• EN4 reference and matching chunks cached for resume.",
    "",
    f"Matched CSV: {MATCHED_CSV}",
    f"Sea bias CSV: {SEA_BIAS_CSV}",
    f"Seasonal correlation CSV: {SEASONAL_CORR_CSV}",
    f"Sea counts CSV: {SEA_COUNT_CSV}",
    f"PNG: {FIG_PNG}",
    f"PDF: {FIG_PDF}",
]

PROCESSING_LOG.write_text("\n".join(log_lines), encoding="utf-8")
log("\n".join(log_lines))
log("\nFigure S3 completed successfully.")
