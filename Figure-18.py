"""
FIGURE 12 — FULL ORGANIZED LATEST WORKFLOW

FIGURE

Sector-Based EOF Analysis
Southern Ocean / 13 Antarctic Shelf Seas

SOURCE
Extracted and organized from:
    so_sealevel_paper_fig.py

VERSION SELECTION
The source contains two Figure 12 stages:

1. A complete EOF-analysis workflow that:
       - builds 216-month x 13-sector SLA and total-steric matrices;
       - removes each sector's monthly climatology;
       - detrends each sector;
       - performs covariance EOF analysis as the main analysis;
       - performs correlation EOF analysis as a sensitivity analysis;
       - calculates the first three EOF modes;
       - calculates standardized PCs;
       - calculates 13-month running-mean PCs;
       - saves all numerical products to CSV.

2. A later PLOT-ONLY workflow explicitly labeled:
       FULL CORRECTED — NO SEA LABELS

   It reads the precomputed EOF products and creates the latest corrected
   Figure 12 without repeating the EOF calculations.

This standalone script combines those two stages correctly:

    full EOF computation
        -> saved numerical products
        -> latest FULL CORRECTED plotting

The older first Figure 12 renderer is intentionally excluded.

SCIENTIFIC SETTINGS
Input matrices:
    216 months x 13 Antarctic sectors

Preprocessing:
    - remove monthly climatology;
    - detrend each sector.

Main EOF analysis:
    covariance EOFs

Sensitivity analysis:
    correlation EOFs

Number of displayed modes:
    first 3 modes

Principal components:
    standardized PCs with 13-month running means

MAIN INPUTS
1. SLA:
   /content/drive/MyDrive/SAM_Thesis/Data/
   SLA_Antarctic_monthly_2008_2025.nc

2. EN4 total steric:
   /content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
   EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc

3. GEBCO bathymetry:
   GEBCO_2024_CEC.nc or GEBCO_2024_CF.nc

NUMERICAL OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/

    Fig12_SLA_sector_matrix_216x13.csv
    Fig12_steric_sector_matrix_216x13.csv
    Fig12_EOF_sector_loadings_summary.csv
    Fig12_EOF_PC_timeseries_summary.csv
    Fig12_correlation_EOF_sensitivity_summary.csv

LATEST FINAL FIGURE OUTPUTS
/content/drive/MyDrive/SAM_Thesis/paper2/

    Figure12_sector_based_EOF_analysis_FROM_PRECOMPUTED_outputs_1080dpi.png
    Figure12_sector_based_EOF_analysis_FROM_PRECOMPUTED_outputs.pdf
    Figure12_sector_based_EOF_analysis_FROM_PRECOMPUTED_outputs_log.txt

LATEST PLOT CORRECTIONS PRESERVED
The final source revision explicitly specifies:

    - do not recompute sector matrices or EOF summaries in the plotting stage;
    - read the precomputed Figure 12 CSV outputs;
    - show Antarctica using GEBCO bathymetry;
    - broadcast each sector EOF loading across ocean cells in that sector;
    - no sea-name labels on EOF maps;
    - latitude labels only at 60°S and 90°S;
    - longitude spacing every 30°;
    - omit the 0° longitude label;
    - place longitude labels outside the subplot without overlap;
    - rotate PC-panel x-axis tick labels by 45°;
    - separate the bottom colorbar from the footer box;
    - save a 1080-dpi PNG and vector PDF.

ORGANIZED EXECUTION ORDER
1. Mount Google Drive.
2. Install required packages.
3. Import all libraries.
4. Define input and output paths.
5. Define Figure 12/EOF settings.
6. Define the 13 Antarctic sector boundaries.
7. Define dataset, preprocessing and EOF helper functions.
8. Verify the SLA and steric input files.
9. Read and standardize the monthly datasets.
10. Build the 216 x 13 SLA sector matrix.
11. Build the 216 x 13 steric sector matrix.
12. Remove monthly climatology and detrend each sector.
13. Save both processed sector matrices.
14. Run covariance EOF analysis.
15. Run correlation EOF sensitivity analysis.
16. Calculate 13-month running-mean PCs.
17. Save EOF loading summaries.
18. Save PC time-series summaries.
19. Save correlation-EOF sensitivity summaries.
20. Read those precomputed products using the latest corrected renderer.
21. Prepare GEBCO bathymetry for map display.
22. Broadcast sector EOF loadings to the gridded ocean mesh.
23. Draw the corrected EOF maps and PC panels.
24. Save the final vector PDF.
25. Save the final 1080-dpi PNG.
26. Save the final Figure 12 plotting log.

ORGANIZATION CHANGES
Only workflow organization has been changed:

    - Google Drive mounting is placed first.
    - Package installation is consolidated.
    - Library imports are consolidated.
    - The full numerical EOF workflow is retained.
    - The superseded first Figure 12 renderer is removed.
    - The later FULL CORRECTED plot-only renderer is retained.
    - All content after Figure 12 is excluded.

No EOF methodology, sector definitions, preprocessing rules, file names,
scientific calculations, sensitivity analysis, or final corrected plotting
instructions are intentionally changed.

"""

# 1. MOUNT GOOGLE DRIVE

try:
    from google.colab import drive

    drive.mount("/content/drive")

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
    "netCDF4": "netCDF4",
    "h5netcdf": "h5netcdf",
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

import matplotlib

# Non-interactive backend is safer for high-resolution Colab export.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates

from matplotlib.path import Path as MplPath

import cartopy.crs as ccrs

warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
)


# PART A — COMPLETE FIGURE 12 EOF DATA PROCESSING
# This section creates all precomputed sector matrices, EOF loadings,
# PC time-series products and correlation-EOF sensitivity outputs.

# FIGURE 12 — SECTOR-BASED EOF ANALYSIS (CORRECTED MAP VERSION)
#
# Main changes requested:
#   1) EOF panels are plotted as south-polar maps, not circular ring charts.
#   2) Antarctica / land mask comes from GEBCO bathymetry.
#   3) EOF loadings are displayed on a gridded ocean mesh by broadcasting
#      each sector loading across the bathymetry ocean cells that belong
#      to that sector.
#   4) No sea labels are drawn on the EOF maps.
#   5) PC-panel x ticks are rotated 45 degrees.
#   6) Layout spacing is tightened so labels and text do not overlap.
#
# Scientific settings:
#   • Input matrices: 216 months × 13 sectors
#   • Remove monthly climatology
#   • Detrend each sector
#   • Covariance EOFs = main analysis
#   • Correlation EOFs = supplementary sensitivity output
#   • PCs shown as 13-month running means
#
# Output folder:
#   /content/drive/MyDrive/SAM_Thesis/paper2


# 3. INPUT / OUTPUT PATHS
SLA_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Data/"
    "SLA_Antarctic_monthly_2008_2025.nc"
)

STERIC_FILE = (
    "/content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/"
    "EN4_total_steric_0_1000m_monthly_2008_2025_SO.nc"
)

GEBCO_CANDIDATES = [
    "/content/drive/MyDrive/SAM_Thesis/Data/GEBCO_2024_CEC.nc",
    "/content/drive/MyDrive/SAM_Thesis/Data/GEBCO_2024_CF.nc",
]

OUTPUT_DIR = Path("/content/drive/MyDrive/SAM_Thesis/paper2")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


OUT_SLA_MATRIX = OUTPUT_DIR / "Fig12_SLA_sector_matrix_216x13.csv"
OUT_STERIC_MATRIX = OUTPUT_DIR / "Fig12_steric_sector_matrix_216x13.csv"
OUT_EOF_LOADINGS = OUTPUT_DIR / "Fig12_EOF_sector_loadings_summary.csv"
OUT_EOF_PCS = OUTPUT_DIR / "Fig12_EOF_PC_timeseries_summary.csv"
OUT_SUPP_CORR = OUTPUT_DIR / "Fig12_correlation_EOF_sensitivity_summary.csv"


# 4. DISPLAY SETTINGS
SAVE_DPI = 1080
FIGSIZE = (16.2, 8.8)

LAT_MIN = -90.0
LAT_MAX = -60.0
LON_MIN = -180.0
LON_MAX = 180.0

FONT_PANEL = 13.0
FONT_AXIS = 10.0
FONT_TICK = 9.0
FONT_MAP_LABEL = 8.5
FONT_LEGEND = 9.0
FONT_NOTE = 9.5

PC_COLORS = ["red", "blue", "green"]
EOF_CMAP = plt.get_cmap("coolwarm")
EOF_NORM = mcolors.TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


# 5. 13 SECTOR DEFINITIONS
SEA_SECTORS = OrderedDict([
    ("WED", {"name": "Weddell Sea",         "lon_min":  -60.0, "lon_max":  -20.0}),
    ("KHV", {"name": "King Haakon VII Sea", "lon_min":  -20.0, "lon_max":   10.0}),
    ("RLS", {"name": "Riiser-Larsen Sea",   "lon_min":   10.0, "lon_max":   35.0}),
    ("LAZ", {"name": "Lazarev Sea",         "lon_min":   35.0, "lon_max":   60.0}),
    ("COS", {"name": "Cosmonauts Sea",      "lon_min":   60.0, "lon_max":   90.0}),
    ("COO", {"name": "Cooperation Sea",     "lon_min":   90.0, "lon_max":  115.0}),
    ("DAV", {"name": "Davis Sea",           "lon_min":  115.0, "lon_max":  130.0}),
    ("MAW", {"name": "Mawson Sea",          "lon_min":  130.0, "lon_max":  150.0}),
    ("DUR", {"name": "D'Urville Sea",       "lon_min":  150.0, "lon_max":  170.0}),
    ("SOM", {"name": "Somov Sea",           "lon_min":  170.0, "lon_max": -160.0}),
    ("ROS", {"name": "Ross Sea",            "lon_min": -160.0, "lon_max": -130.0}),
    ("AMU", {"name": "Amundsen Sea",        "lon_min": -130.0, "lon_max": -100.0}),
    ("BEL", {"name": "Bellingshausen Sea",  "lon_min": -100.0, "lon_max":  -60.0}),
])
SEA_CODES = list(SEA_SECTORS.keys())


# 6. HELPERS
def open_dataset_safely(file_path):
    attempts = []
    for engine in [None, "netcdf4", "h5netcdf", "scipy"]:
        try:
            kwargs = {
                "decode_times": True,
                "mask_and_scale": True,
                "cache": False,
            }
            if engine is not None:
                kwargs["engine"] = engine
            ds = xr.open_dataset(file_path, **kwargs)
            engine_name = "xarray-default" if engine is None else engine
            print(f"Opened {os.path.basename(file_path)} with engine={engine_name}")
            return ds
        except Exception as error:
            attempts.append(f"engine={engine}: {error}")

    raise RuntimeError(
        "Could not open dataset:\n"
        f"{file_path}\n\n"
        + "\n".join(attempts)
    )


def detect_coordinate(dataset, coordinate_type):
    aliases = {
        "lat": ["lat", "latitude", "nav_lat", "y"],
        "lon": ["lon", "longitude", "nav_lon", "x"],
        "time": ["time"],
    }

    candidates = list(dataset.coords) + list(dataset.variables)

    for name in candidates:
        lower = name.lower()
        da = dataset[name]
        standard_name = str(da.attrs.get("standard_name", "")).lower()
        axis = str(da.attrs.get("axis", "")).upper()
        units = str(da.attrs.get("units", "")).lower()

        if coordinate_type == "lat":
            if (
                lower in aliases["lat"]
                or standard_name == "latitude"
                or axis == "Y"
                or "degrees_north" in units
            ):
                return name

        elif coordinate_type == "lon":
            if (
                lower in aliases["lon"]
                or standard_name == "longitude"
                or axis == "X"
                or "degrees_east" in units
            ):
                return name

        elif coordinate_type == "time":
            if lower in aliases["time"] or np.issubdtype(da.dtype, np.datetime64):
                return name

    raise KeyError(f"Could not detect coordinate: {coordinate_type}")


def standardize_dataset(dataset):
    rename_map = {}

    lat_name = detect_coordinate(dataset, "lat")
    lon_name = detect_coordinate(dataset, "lon")
    time_name = None
    try:
        time_name = detect_coordinate(dataset, "time")
    except Exception:
        time_name = None

    if lat_name != "lat":
        rename_map[lat_name] = "lat"
    if lon_name != "lon":
        rename_map[lon_name] = "lon"
    if time_name is not None and time_name != "time":
        rename_map[time_name] = "time"

    if rename_map:
        dataset = dataset.rename(rename_map)

    lon_values = (((dataset["lon"].astype(float) + 180.0) % 360.0) - 180.0)
    dataset = dataset.assign_coords(lon=lon_values)

    lon_values = np.asarray(dataset["lon"].values)
    _, unique_indices = np.unique(lon_values, return_index=True)
    dataset = dataset.isel(lon=np.sort(unique_indices))

    dataset = dataset.sortby("lon").sortby("lat")
    if "time" in dataset.coords:
        dataset = dataset.sortby("time")

    return dataset


def detect_main_variable(dataset, preferred_names):
    for name in preferred_names:
        if name in dataset.data_vars:
            return name

    candidates = []
    for name in dataset.data_vars:
        dims = set(dataset[name].dims)
        if {"lat", "lon"}.issubset(dims):
            candidates.append(name)

    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        print("Candidate variables:", candidates)
        return candidates[0]

    raise KeyError("Could not detect main data variable.")


def convert_to_cm(data_array):
    units = str(data_array.attrs.get("units", "")).strip().lower()

    if units in {"cm", "centimeter", "centimeters", "centimetre", "centimetres"}:
        out = data_array.astype(np.float32)
        out.attrs["units"] = "cm"
        return out

    out = (data_array.astype(np.float32) * 100.0)
    out.attrs["units"] = "cm"
    return out


def monthly_anomaly(time_series):
    climatology = time_series.groupby("time.month").mean("time", skipna=True)
    anomaly = time_series.groupby("time.month") - climatology
    return anomaly


def detrend_series_1d(values):
    values = np.asarray(values, dtype=np.float64)
    x = np.arange(values.size, dtype=np.float64)
    valid = np.isfinite(values)
    out = np.full(values.shape, np.nan, dtype=np.float64)

    if valid.sum() < 2:
        return out

    coef = np.polyfit(x[valid], values[valid], deg=1)
    fit = coef[0] * x + coef[1]
    out[valid] = values[valid] - fit[valid]
    return out


def fill_small_time_gaps(series):
    s = pd.Series(series)
    s = s.interpolate(method="linear", limit_direction="both")
    if s.isna().any():
        s = s.fillna(s.mean())
    if s.isna().any():
        s = s.fillna(0.0)
    return s.values.astype(np.float64)


def sector_mask(lon2d, lat2d, lon_min, lon_max):
    lat_cond = (lat2d >= LAT_MIN) & (lat2d <= LAT_MAX)

    if lon_max < lon_min:
        lon_cond = (lon2d >= lon_min) | (lon2d < lon_max)
    else:
        lon_cond = (lon2d >= lon_min) & (lon2d < lon_max)

    return lat_cond & lon_cond


def weighted_sector_mean(field_2d, mask_2d, lat_2d):
    valid = np.isfinite(field_2d) & mask_2d
    if not np.any(valid):
        return np.nan

    weights = np.cos(np.deg2rad(lat_2d))
    weights = np.where(valid, weights, np.nan)

    numerator = np.nansum(field_2d * weights)
    denominator = np.nansum(weights)

    if not np.isfinite(denominator) or denominator == 0.0:
        return np.nan

    return float(numerator / denominator)


def build_sector_matrix(data_array, validity_fraction_threshold=0.70):
    lon2d, lat2d = np.meshgrid(
        data_array["lon"].values.astype(np.float64),
        data_array["lat"].values.astype(np.float64),
    )

    core_validity = np.isfinite(data_array).mean("time") >= validity_fraction_threshold
    core_validity = core_validity.values.astype(bool)

    time_values = pd.to_datetime(data_array["time"].values)
    sector_df = pd.DataFrame(index=time_values, columns=SEA_CODES, dtype=float)

    for sea_code, sea_info in SEA_SECTORS.items():
        mask = sector_mask(
            lon2d,
            lat2d,
            sea_info["lon_min"],
            sea_info["lon_max"],
        ) & core_validity

        for it in range(data_array.sizes["time"]):
            field = data_array.isel(time=it).values.astype(np.float64)
            sector_df.iloc[it, sector_df.columns.get_loc(sea_code)] = (
                weighted_sector_mean(field, mask, lat2d)
            )

    return sector_df


def preprocess_sector_matrix(sector_df):
    da = xr.DataArray(
        sector_df.values.astype(np.float64),
        coords={"time": sector_df.index, "sector": sector_df.columns},
        dims=("time", "sector"),
    )

    anomaly = monthly_anomaly(da)

    out = np.full(anomaly.shape, np.nan, dtype=np.float64)
    for j in range(anomaly.shape[1]):
        out[:, j] = detrend_series_1d(anomaly[:, j].values)
        out[:, j] = fill_small_time_gaps(out[:, j])

    return pd.DataFrame(out, index=sector_df.index, columns=sector_df.columns)


def run_eof(matrix_df, n_modes=3, use_correlation=False):
    X = matrix_df.values.astype(np.float64)

    if use_correlation:
        mean = X.mean(axis=0, keepdims=True)
        std = X.std(axis=0, ddof=1, keepdims=True)
        std = np.where(std == 0.0, 1.0, std)
        X = (X - mean) / std
    else:
        X = X - X.mean(axis=0, keepdims=True)

    U, singular_values, VT = np.linalg.svd(X, full_matrices=False)
    eigenvalues = (singular_values ** 2) / (X.shape[0] - 1)
    explained_variance = 100.0 * eigenvalues / eigenvalues.sum()

    eofs = VT[:n_modes, :].copy()
    pcs = (U[:, :n_modes] * singular_values[:n_modes]).copy()

    for mode in range(n_modes):
        idx = int(np.nanargmax(np.abs(eofs[mode])))
        if eofs[mode, idx] < 0:
            eofs[mode] *= -1.0
            pcs[:, mode] *= -1.0

    pcs_std = pcs.copy()
    for mode in range(n_modes):
        std = pcs_std[:, mode].std(ddof=1)
        if std == 0.0:
            std = 1.0
        pcs_std[:, mode] = (pcs_std[:, mode] - pcs_std[:, mode].mean()) / std

    return {
        "eofs": eofs,
        "pcs": pcs,
        "pcs_std": pcs_std,
        "explained_variance": explained_variance[:n_modes],
    }


def running_mean_13(values_2d, time_index):
    return pd.DataFrame(values_2d, index=time_index).rolling(
        window=13,
        center=True,
        min_periods=1,
    ).mean()


def normalize_loading_vector(vector):
    vector = np.asarray(vector, dtype=np.float64)
    max_abs = np.nanmax(np.abs(vector))
    if not np.isfinite(max_abs) or max_abs == 0.0:
        max_abs = 1.0
    return vector / max_abs


def find_gebco_file():
    for candidate in GEBCO_CANDIDATES:
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(
        "Could not find a GEBCO bathymetry file.\nChecked:\n"
        + "\n".join(GEBCO_CANDIDATES)
    )


def prepare_bathymetry_grid(gebco_file):
    bathy_ds = standardize_dataset(open_dataset_safely(gebco_file))
    bathy_var = detect_main_variable(
        bathy_ds,
        ["elevation", "z", "depth", "Band1"],
    )

    bathy = bathy_ds[bathy_var].sel(lat=slice(LAT_MIN, LAT_MAX))

    # Reduce memory if GEBCO is very high resolution.
    n_lat = bathy.sizes["lat"]
    n_lon = bathy.sizes["lon"]

    lat_step = max(1, int(np.ceil(n_lat / 420)))
    lon_step = max(1, int(np.ceil(n_lon / 720)))

    bathy = bathy.isel(
        lat=slice(None, None, lat_step),
        lon=slice(None, None, lon_step),
    )

    lon2d, lat2d = np.meshgrid(
        bathy["lon"].values.astype(np.float64),
        bathy["lat"].values.astype(np.float64),
    )

    elevation = bathy.values.astype(np.float32)
    ocean_mask = np.isfinite(elevation) & (elevation < 0.0)
    land_mask = np.isfinite(elevation) & (elevation >= 0.0)

    return {
        "dataset": bathy_ds,
        "variable": bathy_var,
        "bathy": bathy,
        "lon2d": lon2d,
        "lat2d": lat2d,
        "elevation": elevation,
        "ocean_mask": ocean_mask,
        "land_mask": land_mask,
    }


def broadcast_sector_loading_to_grid(loading_vector, lon2d, lat2d, ocean_mask):
    field = np.full(lon2d.shape, np.nan, dtype=np.float32)
    norm_load = normalize_loading_vector(loading_vector)

    for idx, sea_code in enumerate(SEA_CODES):
        info = SEA_SECTORS[sea_code]
        mask = sector_mask(
            lon2d,
            lat2d,
            info["lon_min"],
            info["lon_max"],
        ) & ocean_mask
        field[mask] = norm_load[idx]

    return field


def circular_boundary_path():
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    center = np.array([0.5, 0.5])
    radius = 0.5
    vertices = np.vstack(
        [np.sin(theta), np.cos(theta)]
    ).T * radius + center
    return MplPath(vertices)


def add_polar_labels(ax):
    # Longitude labels around outer boundary
    label_specs = [
        (-180, -62.3, "180°"),
        (-120, -61.3, "120°W"),
        (-60, -61.3, "60°W"),
        (0, -61.0, "0°"),
        (60, -61.3, "60°E"),
        (120, -61.3, "120°E"),
        (180, -62.3, "180°"),
    ]

    for lon, lat, txt in label_specs:
        ax.text(
            lon,
            lat,
            txt,
            transform=ccrs.PlateCarree(),
            fontsize=FONT_MAP_LABEL,
            fontweight="bold",
            ha="center",
            va="center",
            clip_on=False,
            zorder=30,
        )

    # Latitude labels
    for lat, lon, txt in [(-60, 0, "60°S"), (-70, 0, "70°S"), (-80, 0, "80°S")]:
        ax.text(
            lon,
            lat,
            txt,
            transform=ccrs.PlateCarree(),
            fontsize=FONT_MAP_LABEL,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=30,
            bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.55),
        )

    ax.text(
        0.50,
        0.50,
        "90°S",
        transform=ax.transAxes,
        fontsize=FONT_MAP_LABEL + 0.2,
        fontweight="bold",
        ha="center",
        va="center",
        zorder=31,
    )


def style_polar_map_axis(ax):
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=ccrs.PlateCarree())
    ax.set_boundary(circular_boundary_path(), transform=ax.transAxes)

    ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        xlocs=np.arange(-180, 181, 60),
        ylocs=np.arange(-80, -59, 10),
        linewidth=0.50,
        color="0.55",
        linestyle="--",
        alpha=0.65,
        zorder=2,
    )

    add_polar_labels(ax)


def plot_eof_map(ax, field, bathy_dict, title):
    lon2d = bathy_dict["lon2d"]
    lat2d = bathy_dict["lat2d"]
    elevation = bathy_dict["elevation"]
    land_mask = bathy_dict["land_mask"]

    style_polar_map_axis(ax)

    mesh = ax.pcolormesh(
        lon2d,
        lat2d,
        field,
        transform=ccrs.PlateCarree(),
        cmap=EOF_CMAP,
        norm=EOF_NORM,
        shading="auto",
        zorder=4,
    )

    # Antarctica / land from bathymetry
    land_binary = np.where(land_mask, 1.0, np.nan)
    ax.contourf(
        lon2d,
        lat2d,
        land_binary,
        levels=[0.5, 1.5],
        colors=["0.85"],
        transform=ccrs.PlateCarree(),
        zorder=8,
    )

    # Coastline and bathymetric context
    ax.contour(
        lon2d,
        lat2d,
        elevation,
        levels=[0.0],
        colors="0.35",
        linewidths=0.7,
        transform=ccrs.PlateCarree(),
        zorder=9,
    )

    ax.contour(
        lon2d,
        lat2d,
        elevation,
        levels=[-1000.0],
        colors="0.65",
        linewidths=0.45,
        linestyles="--",
        transform=ccrs.PlateCarree(),
        zorder=7,
    )

    ax.set_title(title, fontsize=FONT_PANEL, fontweight="bold", pad=9)
    return mesh


def plot_pc_panel(ax, time_index, pcs_running_mean, explained_variance, title):
    years = pd.to_datetime(time_index)

    for mode in range(3):
        ax.plot(
            years,
            pcs_running_mean.iloc[:, mode].values,
            color=PC_COLORS[mode],
            lw=1.25,
            label=f"PC{mode+1} ({explained_variance[mode]:.1f}%)",
        )

    ax.axhline(0.0, color="0.35", lw=0.9, ls="--")
    ax.grid(True, axis="x", linestyle=":", color="0.72", lw=0.8)

    ax.set_ylabel("Amplitude (std. dev.)", fontsize=FONT_AXIS, fontweight="bold")
    ax.set_xlabel("Year", fontsize=FONT_AXIS, fontweight="bold")
    ax.set_title(title, fontsize=FONT_PANEL, fontweight="bold", pad=8)

    ax.tick_params(axis="both", labelsize=FONT_TICK)

    # X-axis spacing and rotation to prevent overlap
    ax.xaxis.set_major_locator(mdates.YearLocator(3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_ha("right")
        tick.set_fontweight("bold")
    for tick in ax.get_yticklabels():
        tick.set_fontweight("bold")

    y_max = float(np.nanmax(np.abs(pcs_running_mean.values)))
    if not np.isfinite(y_max) or y_max == 0:
        y_max = 1.0
    y_lim = max(2.0, np.ceil(y_max * 1.15 * 2) / 2)
    ax.set_ylim(-y_lim, y_lim)

    ax.legend(
        loc="upper right",
        fontsize=FONT_LEGEND,
        frameon=False,
        handlelength=2.4,
        borderaxespad=0.3,
    )


# 7. VERIFY INPUT FILES
for path in [SLA_FILE, STERIC_FILE]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required file not found:\n{path}")

GEBCO_FILE = find_gebco_file()

print("\nUsing input files:")
print("SLA   :", SLA_FILE)
print("Steric:", STERIC_FILE)
print("GEBCO :", GEBCO_FILE)


# 8. READ DATASETS
sla_ds = standardize_dataset(open_dataset_safely(SLA_FILE))
steric_ds = standardize_dataset(open_dataset_safely(STERIC_FILE))

sla_var = detect_main_variable(sla_ds, ["sla", "SLA", "adt", "sla_anomaly"])
steric_var = detect_main_variable(steric_ds, ["total_steric_height", "steric", "eta_total"])

print("\nDetected variables:")
print("SLA   :", sla_var)
print("Steric:", steric_var)

sla = convert_to_cm(sla_ds[sla_var]).sel(lat=slice(LAT_MIN, LAT_MAX))
steric = convert_to_cm(steric_ds[steric_var]).sel(lat=slice(LAT_MIN, LAT_MAX))

common_time = np.intersect1d(
    pd.to_datetime(sla["time"].values).values,
    pd.to_datetime(steric["time"].values).values,
)

sla = sla.sel(time=common_time)
steric = steric.sel(time=common_time)

print("\nCommon analysis period:")
print(pd.to_datetime(sla["time"].values[0]))
print("to")
print(pd.to_datetime(sla["time"].values[-1]))
print("Months =", sla.sizes["time"])


# 9. BUILD 216 × 13 SECTOR MATRICES
print("\nBuilding 216 × 13 SLA sector matrix...")
sla_sector_raw = build_sector_matrix(sla, validity_fraction_threshold=0.70)

print("Building 216 × 13 steric sector matrix...")
steric_sector_raw = build_sector_matrix(steric, validity_fraction_threshold=0.70)

print("Removing monthly climatology and detrending each sector...")
sla_sector_preprocessed = preprocess_sector_matrix(sla_sector_raw)
steric_sector_preprocessed = preprocess_sector_matrix(steric_sector_raw)

sla_sector_preprocessed.to_csv(OUT_SLA_MATRIX, index_label="time")
steric_sector_preprocessed.to_csv(OUT_STERIC_MATRIX, index_label="time")


# 10. EOF ANALYSIS
print("\nRunning covariance EOF analysis...")
sla_eof = run_eof(sla_sector_preprocessed, n_modes=3, use_correlation=False)
steric_eof = run_eof(steric_sector_preprocessed, n_modes=3, use_correlation=False)

print("Running correlation EOF sensitivity analysis...")
sla_eof_corr = run_eof(sla_sector_preprocessed, n_modes=3, use_correlation=True)
steric_eof_corr = run_eof(steric_sector_preprocessed, n_modes=3, use_correlation=True)

sla_pcs_rm = running_mean_13(
    sla_eof["pcs_std"],
    sla_sector_preprocessed.index,
)
steric_pcs_rm = running_mean_13(
    steric_eof["pcs_std"],
    steric_sector_preprocessed.index,
)


# 11. PREPARE BATHYMETRY DISPLAY GRID
print("\nPreparing bathymetry display grid for map panels...")
bathy_dict = prepare_bathymetry_grid(GEBCO_FILE)

sla_map_fields = [
    broadcast_sector_loading_to_grid(
        sla_eof["eofs"][mode],
        bathy_dict["lon2d"],
        bathy_dict["lat2d"],
        bathy_dict["ocean_mask"],
    )
    for mode in range(3)
]

steric_map_fields = [
    broadcast_sector_loading_to_grid(
        steric_eof["eofs"][mode],
        bathy_dict["lon2d"],
        bathy_dict["lat2d"],
        bathy_dict["ocean_mask"],
    )
    for mode in range(3)
]


# 12. SAVE NUMERICAL OUTPUTS
loading_rows = []
for family_name, eof_result in [
    ("SLA_covariance", sla_eof),
    ("Steric_covariance", steric_eof),
    ("SLA_correlation", sla_eof_corr),
    ("Steric_correlation", steric_eof_corr),
]:
    for mode in range(3):
        normalized = normalize_loading_vector(eof_result["eofs"][mode])
        for sea_index, sea_code in enumerate(SEA_CODES):
            loading_rows.append(
                {
                    "family": family_name,
                    "mode": mode + 1,
                    "explained_variance_percent": eof_result["explained_variance"][mode],
                    "sea_code": sea_code,
                    "loading_raw": eof_result["eofs"][mode, sea_index],
                    "loading_normalized_for_plot": normalized[sea_index],
                }
            )
pd.DataFrame(loading_rows).to_csv(OUT_EOF_LOADINGS, index=False)

pc_rows = []
for family_name, time_index, eof_result, running_mean_df in [
    ("SLA_covariance", sla_sector_preprocessed.index, sla_eof, sla_pcs_rm),
    ("Steric_covariance", steric_sector_preprocessed.index, steric_eof, steric_pcs_rm),
]:
    for mode in range(3):
        for i, time_value in enumerate(time_index):
            pc_rows.append(
                {
                    "family": family_name,
                    "time": pd.to_datetime(time_value),
                    "mode": mode + 1,
                    "explained_variance_percent": eof_result["explained_variance"][mode],
                    "pc_standardized": eof_result["pcs_std"][i, mode],
                    "pc_13month_running_mean": running_mean_df.iloc[i, mode],
                }
            )
pd.DataFrame(pc_rows).to_csv(OUT_EOF_PCS, index=False)

supp_rows = []
for family_name, eof_result in [
    ("SLA_correlation", sla_eof_corr),
    ("Steric_correlation", steric_eof_corr),
]:
    for mode in range(3):
        supp_rows.append(
            {
                "family": family_name,
                "mode": mode + 1,
                "explained_variance_percent": eof_result["explained_variance"][mode],
            }
        )
pd.DataFrame(supp_rows).to_csv(OUT_SUPP_CORR, index=False)


# PART B — LATEST FULL-CORRECTED FIGURE 12 PLOTTING
# This is the later source revision:
# "PLOT-ONLY VERSION FROM PRECOMPUTED OUTPUTS
#  (FULL CORRECTED — NO SEA LABELS)".
#
# It reads the files produced in Part A and does not recompute the EOF analysis.

# FIGURE 12 — SECTOR-BASED EOF ANALYSIS
# PLOT-ONLY VERSION FROM PRECOMPUTED OUTPUTS (FULL CORRECTED — NO SEA LABELS)
#
# User-requested corrections:
#   1. Do NOT recompute sector matrices / EOF summaries.
#   2. Read precomputed outputs:
#        - Fig12_SLA_sector_matrix_216x13.csv
#        - Fig12_steric_sector_matrix_216x13.csv
#        - Fig12_EOF_sector_loadings_summary.csv
#        - Fig12_EOF_PC_timeseries_summary.csv
#        - Fig12_correlation_EOF_sensitivity_summary.csv
#   3. EOF maps:
#        - show Antarctica background
#        - do not include sea labels on EOF maps
#        - latitude labels: only 60°S and 90°S
#        - longitude spacing: 30°
#        - no 0° label
#        - labels outside subplot and no overlap
#   4. PC panels:
#        - clean layout
#        - x ticks rotated 45°
#   5. Bottom colorbar label must NOT overlap with footer box
#   6. Save figure at 1080 dpi


# 2. PATHS
PAPER2_DIR = Path("/content/drive/MyDrive/SAM_Thesis/paper2")
DATA_DIR = Path("/content/drive/MyDrive/SAM_Thesis/Data")

GEBCO_CANDIDATES = [
    DATA_DIR / "GEBCO_2024_CEC.nc",
    DATA_DIR / "GEBCO_2024_CF.nc",
]

SLA_MATRIX_CSV = PAPER2_DIR / "Fig12_SLA_sector_matrix_216x13.csv"
STERIC_MATRIX_CSV = PAPER2_DIR / "Fig12_steric_sector_matrix_216x13.csv"
EOF_LOADINGS_CSV = PAPER2_DIR / "Fig12_EOF_sector_loadings_summary.csv"
EOF_PC_CSV = PAPER2_DIR / "Fig12_EOF_PC_timeseries_summary.csv"
EOF_CORR_CSV = PAPER2_DIR / "Fig12_correlation_EOF_sensitivity_summary.csv"

OUT_PNG = PAPER2_DIR / "Figure12_sector_based_EOF_analysis_FROM_PRECOMPUTED_outputs_1080dpi.png"
OUT_PDF = PAPER2_DIR / "Figure12_sector_based_EOF_analysis_FROM_PRECOMPUTED_outputs.pdf"
OUT_LOG = PAPER2_DIR / "Figure12_sector_based_EOF_analysis_FROM_PRECOMPUTED_outputs_log.txt"

for required_path in [
    SLA_MATRIX_CSV, STERIC_MATRIX_CSV, EOF_LOADINGS_CSV, EOF_PC_CSV
]:
    if not required_path.exists():
        raise FileNotFoundError(f"Required file not found:\n{required_path}")


# 3. STYLING
SAVE_DPI = 1080
FIGSIZE = (16.5, 9.2)

LAT_MIN = -90.0
LAT_MAX = -60.0
LON_MIN = -180.0
LON_MAX = 180.0

FONT_TITLE = 13.0
FONT_AXIS = 10.0
FONT_TICK = 8.7
FONT_MAP_LABEL = 10
FONT_FOOT = 9.0
FONT_LEGEND = 9.0

EOF_CMAP = plt.get_cmap("coolwarm")
EOF_NORM = mcolors.TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
PC_COLORS = ["red", "blue", "green"]

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


# 4. ANTARCTIC SHELF-SEA SECTOR DEFINITIONS
# These longitude limits are required to broadcast each of the 13 EOF
# sector loadings onto the Antarctic map. They do not create text labels.
SEA_SECTORS = OrderedDict([
    ("WED", {"lon_min": -60.0,  "lon_max": -20.0}),
    ("KHV", {"lon_min": -20.0,  "lon_max":  10.0}),
    ("RLS", {"lon_min":  10.0,  "lon_max":  35.0}),
    ("LAZ", {"lon_min":  35.0,  "lon_max":  60.0}),
    ("COS", {"lon_min":  60.0,  "lon_max":  90.0}),
    ("COO", {"lon_min":  90.0,  "lon_max": 115.0}),
    ("DAV", {"lon_min": 115.0,  "lon_max": 130.0}),
    ("MAW", {"lon_min": 130.0,  "lon_max": 150.0}),
    ("DUR", {"lon_min": 150.0,  "lon_max": 170.0}),
    ("SOM", {"lon_min": 170.0,  "lon_max": -160.0}),
    ("ROS", {"lon_min": -160.0, "lon_max": -130.0}),
    ("AMU", {"lon_min": -130.0, "lon_max": -100.0}),
    ("BEL", {"lon_min": -100.0, "lon_max": -60.0}),
])

SEA_CODES = list(SEA_SECTORS.keys())


# 5. HELPERS
def find_gebco_file():
    for path in GEBCO_CANDIDATES:
        if path.exists():
            return path
    raise FileNotFoundError(
        "No GEBCO file found. Checked:\n" + "\n".join([str(p) for p in GEBCO_CANDIDATES])
    )


def open_dataset_safely(path):
    attempts = []
    for engine in [None, "netcdf4", "h5netcdf", "scipy"]:
        try:
            kwargs = dict(decode_times=True, mask_and_scale=True, cache=False)
            if engine is not None:
                kwargs["engine"] = engine
            ds = xr.open_dataset(path, **kwargs)
            print(f"Opened {Path(path).name} with engine={engine or 'xarray-default'}")
            return ds
        except Exception as e:
            attempts.append(f"{engine}: {e}")
    raise RuntimeError("Could not open dataset:\n" + str(path) + "\n" + "\n".join(attempts))


def detect_coordinate(dataset, kind):
    aliases = {
        "lat": ["lat", "latitude", "y", "nav_lat"],
        "lon": ["lon", "longitude", "x", "nav_lon"],
    }
    for name in list(dataset.coords) + list(dataset.variables):
        lower = name.lower()
        da = dataset[name]
        std = str(da.attrs.get("standard_name", "")).lower()
        units = str(da.attrs.get("units", "")).lower()
        axis = str(da.attrs.get("axis", "")).upper()

        if kind == "lat":
            if lower in aliases["lat"] or std == "latitude" or axis == "Y" or "degrees_north" in units:
                return name
        elif kind == "lon":
            if lower in aliases["lon"] or std == "longitude" or axis == "X" or "degrees_east" in units:
                return name
    raise KeyError(f"Could not detect coordinate {kind}")


def standardize_dataset(dataset):
    rename_map = {}
    lat_name = detect_coordinate(dataset, "lat")
    lon_name = detect_coordinate(dataset, "lon")
    if lat_name != "lat":
        rename_map[lat_name] = "lat"
    if lon_name != "lon":
        rename_map[lon_name] = "lon"
    if rename_map:
        dataset = dataset.rename(rename_map)

    lon = (((dataset["lon"].astype(float) + 180.0) % 360.0) - 180.0)
    dataset = dataset.assign_coords(lon=lon)

    lon_values = np.asarray(dataset["lon"].values)
    _, idx = np.unique(lon_values, return_index=True)
    dataset = dataset.isel(lon=np.sort(idx))

    dataset = dataset.sortby("lon").sortby("lat")
    return dataset


def detect_main_variable(dataset, preferred_names):
    for name in preferred_names:
        if name in dataset.data_vars:
            return name
    candidates = []
    for name in dataset.data_vars:
        if {"lat", "lon"}.issubset(set(dataset[name].dims)):
            candidates.append(name)
    if len(candidates) == 0:
        raise KeyError("Could not detect main variable.")
    print("Candidate bathymetry/data variables:", candidates)
    return candidates[0]


def normalize_loading_vector(vector):
    vector = np.asarray(vector, dtype=np.float64)
    max_abs = np.nanmax(np.abs(vector))
    if not np.isfinite(max_abs) or max_abs == 0:
        max_abs = 1.0
    return vector / max_abs


def sector_mask(lon2d, lat2d, lon_min, lon_max):
    lat_cond = (lat2d >= LAT_MIN) & (lat2d <= LAT_MAX)
    if lon_max < lon_min:
        lon_cond = (lon2d >= lon_min) | (lon2d < lon_max)
    else:
        lon_cond = (lon2d >= lon_min) & (lon2d < lon_max)
    return lat_cond & lon_cond


def circular_boundary_path():
    theta = np.linspace(0, 2*np.pi, 361)
    center = np.array([0.5, 0.5])
    radius = 0.5
    vertices = np.vstack([np.sin(theta), np.cos(theta)]).T * radius + center
    return MplPath(vertices)


def prepare_bathymetry_grid(gebco_file):
    ds = standardize_dataset(open_dataset_safely(gebco_file))
    var = detect_main_variable(ds, ["elevation", "z", "depth", "Band1"])
    bathy = ds[var].sel(lat=slice(LAT_MIN, LAT_MAX))

    n_lat = bathy.sizes["lat"]
    n_lon = bathy.sizes["lon"]
    lat_step = max(1, int(np.ceil(n_lat / 420)))
    lon_step = max(1, int(np.ceil(n_lon / 720)))

    bathy = bathy.isel(lat=slice(None, None, lat_step), lon=slice(None, None, lon_step))
    lon2d, lat2d = np.meshgrid(
        bathy["lon"].values.astype(np.float64),
        bathy["lat"].values.astype(np.float64),
    )

    elevation = bathy.values.astype(np.float32)
    ocean_mask = np.isfinite(elevation) & (elevation < 0.0)
    land_mask = np.isfinite(elevation) & (elevation >= 0.0)

    return {
        "dataset": ds,
        "var": var,
        "bathy": bathy,
        "lon2d": lon2d,
        "lat2d": lat2d,
        "elevation": elevation,
        "ocean_mask": ocean_mask,
        "land_mask": land_mask,
    }


def broadcast_sector_loading_to_grid(loading_vector, lon2d, lat2d, ocean_mask):
    field = np.full(lon2d.shape, np.nan, dtype=np.float32)
    norm_load = normalize_loading_vector(loading_vector)

    for i, sea_code in enumerate(SEA_CODES):
        info = SEA_SECTORS[sea_code]
        mask = sector_mask(lon2d, lat2d, info["lon_min"], info["lon_max"]) & ocean_mask
        field[mask] = norm_load[i]

    return field


def sector_mid_lon(lon_min, lon_max):
    if lon_max < lon_min:
        span = (lon_max + 360.0) - lon_min
        mid = lon_min + span / 2.0
        if mid > 180.0:
            mid -= 360.0
        return mid
    return 0.5 * (lon_min + lon_max)


def style_pc_axis(ax, title):
    ax.set_title(title, fontsize=FONT_TITLE, fontweight="bold", pad=8)
    ax.set_ylabel("Amplitude (std. dev.)", fontsize=FONT_AXIS, fontweight="bold")
    ax.set_xlabel("Year", fontsize=FONT_AXIS, fontweight="bold")
    ax.tick_params(axis="both", labelsize=FONT_TICK)

    ax.xaxis.set_major_locator(mdates.YearLocator(3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_ha("right")
        tick.set_fontweight("bold")
    for tick in ax.get_yticklabels():
        tick.set_fontweight("bold")

    ax.grid(True, axis="x", linestyle=":", color="0.72", linewidth=0.8)
    ax.axhline(0.0, color="0.35", linewidth=0.9, linestyle="--")


def add_polar_grid_and_labels(ax):
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=ccrs.PlateCarree())
    ax.set_boundary(circular_boundary_path(), transform=ax.transAxes)

    # Grid lines: lon every 30°, lat at 60 and 75 (for geometry), but labels only 60 and 90
    ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        xlocs=np.arange(-180, 181, 30),
        ylocs=[-60, -75],
        linewidth=0.45,
        color="0.60",
        linestyle="--",
        alpha=0.65,
        zorder=2,
    )

    # Longitude labels outside each subplot, every 30°, no 0°
    lon_labels = [
        (-180, "180°"),
        (-150, "150°W"),
        (-120, "120°W"),
        (-90, "90°W"),
        (-60, "60°W"),
        (-30, "30°W"),
        (30, "30°E"),
        (60, "60°E"),
        (120, "120°"),
        (150, "150°E"),
        (180, "180°"),
    ]
    for lon, txt in lon_labels:
        ax.text(
            lon, -58.4, txt,
            transform=ccrs.PlateCarree(),
            fontsize=FONT_MAP_LABEL,
            fontweight="bold",
            ha="center", va="center",
            clip_on=False, zorder=30,
        )

    # Latitude labels: only 60°S and 90°S
    ax.text(
        0, -60.7, "60°S",
        transform=ccrs.PlateCarree(),
        fontsize=FONT_MAP_LABEL,
        fontweight="bold",
        ha="center", va="center",
        zorder=30,
        bbox=dict(boxstyle="round,pad=0.05", fc="white", ec="none", alpha=0.50),
    )
    ax.text(
        0.5, 0.5, "90°S",
        transform=ax.transAxes,
        fontsize=FONT_MAP_LABEL + 0.2,
        fontweight="bold",
        ha="center", va="center",
        zorder=31,
    )





def plot_eof_map(ax, field, bathy_dict, title):
    add_polar_grid_and_labels(ax)

    mesh = ax.pcolormesh(
        bathy_dict["lon2d"],
        bathy_dict["lat2d"],
        field,
        transform=ccrs.PlateCarree(),
        cmap=EOF_CMAP,
        norm=EOF_NORM,
        shading="auto",
        zorder=4,
    )

    # Antarctica land
    land_binary = np.where(bathy_dict["land_mask"], 1.0, np.nan)
    ax.contourf(
        bathy_dict["lon2d"], bathy_dict["lat2d"], land_binary,
        levels=[0.5, 1.5],
        colors=["0.85"],
        transform=ccrs.PlateCarree(),
        zorder=8,
    )

    # Coastline
    ax.contour(
        bathy_dict["lon2d"], bathy_dict["lat2d"], bathy_dict["elevation"],
        levels=[0.0],
        colors="0.35",
        linewidths=0.7,
        transform=ccrs.PlateCarree(),
        zorder=9,
    )

    # Sea-name labels are intentionally omitted from all EOF maps.
    ax.set_title(title, fontsize=FONT_TITLE, fontweight="bold", pad=10)
    return mesh


# 6. READ PRECOMPUTED TABLES
print("Reading precomputed outputs...")

sla_matrix = pd.read_csv(SLA_MATRIX_CSV, parse_dates=["time"]).set_index("time")
steric_matrix = pd.read_csv(STERIC_MATRIX_CSV, parse_dates=["time"]).set_index("time")

loadings_df = pd.read_csv(EOF_LOADINGS_CSV)
pcs_df = pd.read_csv(EOF_PC_CSV, parse_dates=["time"])

corr_df = None
if EOF_CORR_CSV.exists():
    try:
        corr_df = pd.read_csv(EOF_CORR_CSV)
    except Exception:
        corr_df = None

gebco_file = find_gebco_file()
bathy_dict = prepare_bathymetry_grid(gebco_file)

print("Using GEBCO:", gebco_file)


# 7. EXTRACT COVARIANCE EOF LOADINGS AND PCS
def get_mode_loading(family_name, mode_number):
    sub = loadings_df[
        (loadings_df["family"] == family_name) &
        (loadings_df["mode"] == mode_number)
    ].copy()

    if sub.empty:
        raise ValueError(f"No loadings found for {family_name}, mode {mode_number}")

    sub["sea_code"] = pd.Categorical(sub["sea_code"], categories=SEA_CODES, ordered=True)
    sub = sub.sort_values("sea_code")

    explained = float(sub["explained_variance_percent"].iloc[0])

    if "loading_normalized_for_plot" in sub.columns:
        load = sub["loading_normalized_for_plot"].to_numpy(dtype=float)
    else:
        load = sub["loading_raw"].to_numpy(dtype=float)
        load = normalize_loading_vector(load)

    return explained, load


def get_pc_running_mean(family_name):
    sub = pcs_df[pcs_df["family"] == family_name].copy()
    if sub.empty:
        raise ValueError(f"No PC data found for {family_name}")

    modes = sorted(sub["mode"].unique().tolist())
    pivot = sub.pivot(index="time", columns="mode", values="pc_13month_running_mean").sort_index()
    pivot = pivot[[1, 2, 3]]
    explained = []
    for m in [1, 2, 3]:
        explained.append(float(sub.loc[sub["mode"] == m, "explained_variance_percent"].iloc[0]))
    return pivot, explained


sla_expl = []
sla_fields = []
for mode in [1, 2, 3]:
    ev, load = get_mode_loading("SLA_covariance", mode)
    sla_expl.append(ev)
    sla_fields.append(
        broadcast_sector_loading_to_grid(
            load,
            bathy_dict["lon2d"],
            bathy_dict["lat2d"],
            bathy_dict["ocean_mask"],
        )
    )

steric_expl = []
steric_fields = []
for mode in [1, 2, 3]:
    ev, load = get_mode_loading("Steric_covariance", mode)
    steric_expl.append(ev)
    steric_fields.append(
        broadcast_sector_loading_to_grid(
            load,
            bathy_dict["lon2d"],
            bathy_dict["lat2d"],
            bathy_dict["ocean_mask"],
        )
    )

sla_pcs_rm, sla_pc_expl = get_pc_running_mean("SLA_covariance")
steric_pcs_rm, steric_pc_expl = get_pc_running_mean("Steric_covariance")


# 8. PLOT FIGURE
fig = plt.figure(figsize=FIGSIZE)

# tighter and cleaner layout; more bottom room so cbar label doesn't overlap footer
gs = fig.add_gridspec(
    nrows=2,
    ncols=4,
    left=0.040,
    right=0.985,
    top=0.965,
    bottom=0.170,
    wspace=0.14,
    hspace=0.26,
)

proj = ccrs.SouthPolarStereo()

# --- top row ---
ax_a = fig.add_subplot(gs[0, 0], projection=proj)
plot_eof_map(ax_a, sla_fields[0], bathy_dict, f"(a) SLA EOF1 ({sla_expl[0]:.1f}%)")

ax_b = fig.add_subplot(gs[0, 1], projection=proj)
plot_eof_map(ax_b, sla_fields[1], bathy_dict, f"(b) SLA EOF2 ({sla_expl[1]:.1f}%)")

ax_c = fig.add_subplot(gs[0, 2], projection=proj)
plot_eof_map(ax_c, sla_fields[2], bathy_dict, f"(c) SLA EOF3 ({sla_expl[2]:.1f}%)")

ax_d = fig.add_subplot(gs[0, 3])
for i, color in enumerate(PC_COLORS, start=1):
    ax_d.plot(
        sla_pcs_rm.index,
        sla_pcs_rm[i].values,
        color=color,
        linewidth=1.25,
        label=f"PC{i} ({sla_pc_expl[i-1]:.1f}%)"
    )
style_pc_axis(ax_d, "(d) SLA PCs (13-month running mean)")
ax_d.legend(loc="upper right", fontsize=FONT_LEGEND, frameon=False)

# --- bottom row ---
ax_e = fig.add_subplot(gs[1, 0], projection=proj)
plot_eof_map(ax_e, steric_fields[0], bathy_dict, f"(e) Steric EOF1 ({steric_expl[0]:.1f}%)")

ax_f = fig.add_subplot(gs[1, 1], projection=proj)
plot_eof_map(ax_f, steric_fields[1], bathy_dict, f"(f) Steric EOF2 ({steric_expl[1]:.1f}%)")

ax_g = fig.add_subplot(gs[1, 2], projection=proj)
plot_eof_map(ax_g, steric_fields[2], bathy_dict, f"(g) Steric EOF3 ({steric_expl[2]:.1f}%)")

ax_h = fig.add_subplot(gs[1, 3])
for i, color in enumerate(PC_COLORS, start=1):
    ax_h.plot(
        steric_pcs_rm.index,
        steric_pcs_rm[i].values,
        color=color,
        linewidth=1.25,
        label=f"PC{i} ({steric_pc_expl[i-1]:.1f}%)"
    )
style_pc_axis(ax_h, "(h) Steric PCs (13-month running mean)")
ax_h.legend(loc="upper right", fontsize=FONT_LEGEND, frameon=False)

# Shared colorbar for EOF maps: move it well above footer box
cax = fig.add_axes([0.075, 0.103, 0.53, 0.017])
cb = fig.colorbar(
    plt.cm.ScalarMappable(norm=EOF_NORM, cmap=EOF_CMAP),
    cax=cax,
    orientation="horizontal",
    extend="both",
)
cb.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
cb.set_label("EOF loading", fontsize=FONT_AXIS, fontweight="bold", labelpad=2)
cb.ax.tick_params(labelsize=FONT_TICK, length=2.5, pad=1.5)
for tick in cb.ax.get_xticklabels():
    tick.set_fontweight("bold")

# Footer box: place lower so it does not overlap cbar label
footer = (
    "Input: 216 months × 13 sectors    |    "
    "Covariance EOFs (main analysis)    |    "
    "Explained variance in parentheses (%)    |    "
    "PCs are 13-month running means"
)
fig.text(
    0.515,
    0.047,
    footer,
    ha="center",
    va="center",
    fontsize=FONT_FOOT,
    fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.34", fc="white", ec="0.55", lw=0.8),
)


# 9. SAVE
fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white", edgecolor="none")
fig.savefig(OUT_PNG, dpi=SAVE_DPI, bbox_inches="tight", facecolor="white", edgecolor="none")
plt.close(fig)

log_lines = [
    "=" * 88,
    "FIGURE 12 — SECTOR-BASED EOF ANALYSIS (PLOT-ONLY FROM PRECOMPUTED OUTPUTS)",
    "=" * 88,
    "",
    "Used precomputed files:",
    f"SLA sector matrix       : {SLA_MATRIX_CSV}",
    f"Steric sector matrix    : {STERIC_MATRIX_CSV}",
    f"EOF loadings summary    : {EOF_LOADINGS_CSV}",
    f"PC time-series summary  : {EOF_PC_CSV}",
    f"Correlation sensitivity : {EOF_CORR_CSV}",
    "",
    "Corrections applied:",
    "• No recomputation of sector matrices / EOF statistics.",
    "• Antarctica shown from GEBCO bathymetry.",
    "• Sea-name labels omitted from all EOF maps.",
    "• Latitude labels only: 60°S and 90°S.",
    "• Longitude labels every 30° outside each subplot, with no 0° label.",
    "• PC x ticks rotated 45°.",
    "• Colorbar and footer separated to avoid overlap.",
    "• Output saved at 1080 dpi.",
    "",
    f"Output PNG: {OUT_PNG}",
    f"Output PDF: {OUT_PDF}",
]
OUT_LOG.write_text("\n".join(log_lines), encoding="utf-8")

print("\n".join(log_lines))
print("\nFigure 12 plot-only corrected version completed successfully.")
