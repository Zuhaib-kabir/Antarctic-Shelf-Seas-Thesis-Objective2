# Seasonal Vertical Temperature–Salinity Profiles
# Second Group of Antarctic Marginal Seas
#
# Seas:
# MAW, DUR, SOM, ROS, AMU, BEL
#
# FINAL LAYOUT
# Temperature tick labels : first row only
# Salinity tick labels    : last row only
# Season labels           : 22 pt
# Sea labels              : 22 pt
# Temperature label       : 18 pt
# Salinity label          : 18 pt
# Depth label             : 18 pt
# X-axis tick labels      : 17 pt
# Y-axis tick labels      : 18 pt
# Legend/text-box words   : 16 pt
#
# Additional improvements:
#   - Compact vertical spacing
#   - Season and temperature labels placed closer together
#   - Dynamic legend placement
#   - Legend remains fully inside each panel
#   - Opaque legend prevents curves passing through text
#   - Ross Sea dateline-crossing sector handled correctly

#Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')

# 0) INSTALL REQUIRED PACKAGES — GOOGLE COLAB
!pip -q install netCDF4 numpy pandas matplotlib


# 1) IMPORTS
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from datetime import datetime, timedelta
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# 2) INPUT AND OUTPUT PATHS
ARGO_NC = (
    "/content/drive/MyDrive/SAM_Thesis/Data/"
    "argo_SO_profiles_2001_2025_cleaned_gridded.nc"
)

out_dir = "/content/drive/MyDrive/SAM_Thesis/Fig"
os.makedirs(out_dir, exist_ok=True)

out_file = os.path.join(
    out_dir,
    "Obj1Fig12_SO_Seasonal_TS_Profiles_Final_Compact.png"
)


# 3) ANALYSIS SETTINGS
YEAR_START = 2008
YEAR_END = 2025

LAT_MIN = -90.0
LAT_MAX = -60.0

PLOT_P_MIN = 0.0
PLOT_P_MAX = 2000.0

# Six rows instead of seven, so slightly reduced height.
FIGSIZE = (27, 25.5)

DPI_SAVE = 1080


# 4) FONT AND LINE SETTINGS
FONT_SEASON_TITLE = 22
FONT_SEA_LABEL = 22

FONT_TEMP_LABEL = 18
FONT_SAL_LABEL = 18
FONT_DEPTH_LABEL = 18

FONT_X_TICK = 17
FONT_Y_TICK = 18

FONT_LEGEND = 16
FONT_NO_DATA = 16

PROFILE_LINEWIDTH = 2.8
SPINE_LINEWIDTH = 1.2
GRID_LINEWIDTH = 0.65


# 5) GLOBAL MATPLOTLIB SETTINGS
plt.rcParams.update({
    "font.family": "DejaVu Sans",

    "axes.labelweight": "bold",
    "axes.linewidth": SPINE_LINEWIDTH,

    "xtick.labelsize": FONT_X_TICK,
    "ytick.labelsize": FONT_Y_TICK,

    "legend.fontsize": FONT_LEGEND,

    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",

    "mathtext.default": "regular",
})


# 6) ANTARCTIC MARGINAL-SEA SECTORS
sea_info = [
    ("MAW", "Mawson Sea",          "90°E–130°E",   90,  130),
    ("DUR", "D'Urville Sea",       "130°E–150°E", 130,  150),
    ("SOM", "Somov Sea",           "150°E–170°E", 150,  170),
    ("ROS", "Ross Sea",            "170°E–130°W", 170, -130),
    ("AMU", "Amundsen Sea",        "130°W–100°W", -130, -100),
    ("BEL", "Bellingshausen Sea",  "100°W–60°W",  -100,  -60),
]

sea_codes = [
    item[0]
    for item in sea_info
]

sea_names = {
    item[0]: item[1]
    for item in sea_info
}

sea_ranges = {
    item[0]: (item[3], item[4])
    for item in sea_info
}


# 7) HELPER FUNCTIONS
def norm_lon(lon):
    """
    Convert longitude values to the range -180° to 180°.
    """
    return ((lon + 180.0) % 360.0) - 180.0


def lon_in_range(lon_arr, lon_min, lon_max):
    """
    Select longitude values within a specified sea sector.

    Supports:
      1. Standard longitude sectors
      2. Dateline-crossing sectors

    Example:
      Ross Sea = 170°E to 130°W
               = 170 to -130
    """
    if lon_min <= lon_max:
        return (
            (lon_arr >= lon_min) &
            (lon_arr <= lon_max)
        )

    return (
        (lon_arr >= lon_min) |
        (lon_arr <= lon_max)
    )


# 8) SOUTHERN HEMISPHERE SEASONS
SEASON_ORDER = [
    "SON",
    "DJF",
    "MAM",
    "JJA",
]

SEASON_TITLE = {
    "SON": "Spring",
    "DJF": "Summer",
    "MAM": "Autumn",
    "JJA": "Winter",
}

SEASON_MONTHS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
}


def month_to_season(month):
    """
    Convert calendar month to Southern Hemisphere season.
    """
    for season, season_months in SEASON_MONTHS.items():

        if month in season_months:
            return season

    return None



# 9) READ ARGO DATA
print("Reading Argo dataset...")

with Dataset(ARGO_NC, mode="r") as ds:

    P = np.asarray(
        ds.variables["PRES_GRID"][:],
        dtype=np.float64
    )

    LAT = np.asarray(
        ds.variables["LATITUDE"][:],
        dtype=np.float64
    )

    LON = norm_lon(
        np.asarray(
            ds.variables["LONGITUDE"][:],
            dtype=np.float64
        )
    )

    JULD = np.asarray(
        ds.variables["JULD"][:],
        dtype=np.float64
    )

    TEMP = np.asarray(
        ds.variables["TEMP"][:],
        dtype=np.float64
    )

    PSAL = np.asarray(
        ds.variables["PSAL"][:],
        dtype=np.float64
    )

print("Argo dataset loaded successfully.")


# 10) CONVERT ARGO JULIAN DAYS TO DATETIME
argo_base_date = datetime(1950, 1, 1)

time = np.array([
    np.datetime64(
        argo_base_date + timedelta(days=float(day))
    )
    if np.isfinite(day)
    else np.datetime64("NaT")
    for day in JULD
])

time_pd = pd.DatetimeIndex(
    pd.to_datetime(time)
)

years = time_pd.year.to_numpy()
months = time_pd.month.to_numpy()

seasons = np.array(
    [
        month_to_season(month)
        if pd.notna(month)
        else None
        for month in months
    ],
    dtype=object
)


# 11) BASIC PROFILE AND PRESSURE MASKS
mask_base = (
    np.isfinite(LAT) &
    np.isfinite(LON) &
    (LAT >= LAT_MIN) &
    (LAT <= LAT_MAX) &
    (years >= YEAR_START) &
    (years <= YEAR_END) &
    pd.notna(seasons)
)

pressure_mask = (
    np.isfinite(P) &
    (P >= PLOT_P_MIN) &
    (P <= PLOT_P_MAX)
)

pressure_indices = np.where(
    pressure_mask
)[0]

P_selected = P[
    pressure_indices
]

if pressure_indices.size == 0:
    raise ValueError(
        "No pressure levels were found between "
        f"{PLOT_P_MIN:.0f} and {PLOT_P_MAX:.0f} m."
    )


# 12) COLLECT SEASONAL MEAN PROFILES
def collect_mean_profiles(sea_code, season_key):
    """
    Calculate seasonal mean temperature and salinity profiles
    for a selected Antarctic marginal sea.

    Returns
    -------
    mean_salinity : numpy.ndarray
        Mean salinity profile.

    mean_temperature : numpy.ndarray
        Mean temperature profile.

    depth : numpy.ndarray
        Selected pressure/depth levels.

    n_profiles : int
        Number of profiles included.
    """

    lon_min, lon_max = sea_ranges[sea_code]

    profile_mask = (
        mask_base &
        lon_in_range(
            LON,
            lon_min,
            lon_max
        ) &
        (seasons == season_key)
    )

    profile_indices = np.where(
        profile_mask
    )[0]

    if profile_indices.size == 0:
        return None

    temperature_subset = TEMP[
        profile_indices
    ][:, pressure_indices]

    salinity_subset = PSAL[
        profile_indices
    ][:, pressure_indices]

    # Basic physical quality-control ranges
    temperature_subset = np.where(
        np.isfinite(temperature_subset) &
        (temperature_subset > -3.0) &
        (temperature_subset < 20.0),
        temperature_subset,
        np.nan
    )

    salinity_subset = np.where(
        np.isfinite(salinity_subset) &
        (salinity_subset > 0.0) &
        (salinity_subset < 42.0),
        salinity_subset,
        np.nan
    )

    # Some depth levels may contain only missing values.
    with warnings.catch_warnings():

        warnings.simplefilter(
            "ignore",
            category=RuntimeWarning
        )

        mean_temperature = np.nanmean(
            temperature_subset,
            axis=0
        )

        mean_salinity = np.nanmean(
            salinity_subset,
            axis=0
        )

    valid_temperature_levels = np.sum(
        np.isfinite(mean_temperature)
    )

    valid_salinity_levels = np.sum(
        np.isfinite(mean_salinity)
    )

    if (
        valid_temperature_levels < 10 or
        valid_salinity_levels < 10
    ):
        return None

    return (
        mean_salinity,
        mean_temperature,
        P_selected,
        profile_indices.size
    )


# 13) AXIS LIMITS AND TICK LOCATIONS
# Temperature axis
T_XMIN = -2.0
T_XMAX = 2.5

T_XTICKS = np.arange(
    T_XMIN,
    T_XMAX + 0.001,
    0.5
)

# Salinity axis
S_XMIN = 33.75
S_XMAX = 34.75

S_XTICKS = np.arange(
    S_XMIN,
    S_XMAX + 0.001,
    0.25
)

# Depth axis
DEPTH_TICKS = np.arange(
    0,
    2000 + 1,
    250
)


# 14) AUTOMATIC LEGEND POSITION
def choose_legend_position(
    mean_salinity,
    mean_temperature,
    depth
):
    """
    Compare lower-left and lower-right legend positions.

    The position containing fewer temperature and salinity
    profile points is selected. This reduces overlap between
    the legend and plotted curves.
    """

    # Convert depth to normalized axis coordinates.
    # Surface is near 1.0 and 2000 m is near 0.0.
    y_norm = 1.0 - (
        (depth - PLOT_P_MIN) /
        (PLOT_P_MAX - PLOT_P_MIN)
    )

    # Normalize salinity to the main-axis coordinate range.
    salinity_norm = (
        (mean_salinity - S_XMIN) /
        (S_XMAX - S_XMIN)
    )

    # Normalize temperature to the secondary-axis range.
    temperature_norm = (
        (mean_temperature - T_XMIN) /
        (T_XMAX - T_XMIN)
    )

    # Approximate legend dimensions including a safety margin.
    legend_width = 0.46
    legend_height = 0.23

    lower_y_min = 0.015
    lower_y_max = lower_y_min + legend_height

    # Lower-left candidate
    left_x_min = 0.015
    left_x_max = left_x_min + legend_width

    left_salinity_overlap = np.sum(
        np.isfinite(salinity_norm) &
        np.isfinite(y_norm) &
        (salinity_norm >= left_x_min) &
        (salinity_norm <= left_x_max) &
        (y_norm >= lower_y_min) &
        (y_norm <= lower_y_max)
    )

    left_temperature_overlap = np.sum(
        np.isfinite(temperature_norm) &
        np.isfinite(y_norm) &
        (temperature_norm >= left_x_min) &
        (temperature_norm <= left_x_max) &
        (y_norm >= lower_y_min) &
        (y_norm <= lower_y_max)
    )

    left_score = (
        left_salinity_overlap +
        left_temperature_overlap
    )

    # Lower-right candidate
    right_x_max = 0.985
    right_x_min = right_x_max - legend_width

    right_salinity_overlap = np.sum(
        np.isfinite(salinity_norm) &
        np.isfinite(y_norm) &
        (salinity_norm >= right_x_min) &
        (salinity_norm <= right_x_max) &
        (y_norm >= lower_y_min) &
        (y_norm <= lower_y_max)
    )

    right_temperature_overlap = np.sum(
        np.isfinite(temperature_norm) &
        np.isfinite(y_norm) &
        (temperature_norm >= right_x_min) &
        (temperature_norm <= right_x_max) &
        (y_norm >= lower_y_min) &
        (y_norm <= lower_y_max)
    )

    right_score = (
        right_salinity_overlap +
        right_temperature_overlap
    )

    if left_score <= right_score:
        return "lower left", (0.018, 0.018)

    return "lower right", (0.982, 0.018)


# 15) CREATE FIGURE

n_rows = len(sea_codes)
n_columns = len(SEASON_ORDER)

fig, axes = plt.subplots(
    nrows=n_rows,
    ncols=n_columns,
    figsize=FIGSIZE,
    sharey=True,
    squeeze=False
)

# Compact spacing:
#   - Temperature ticks only appear on the first row.
#   - Salinity ticks only appear on the last row.
#   - Season titles are positioned close to temperature labels.
fig.subplots_adjust(
    left=0.120,
    right=0.975,
    top=0.925,
    bottom=0.075,

    wspace=0.30,
    hspace=0.14
)


# 16) PLOT ALL SEA–SEASON PANELS
for row_index, sea_code in enumerate(sea_codes):

    for column_index, season_key in enumerate(SEASON_ORDER):

        ax = axes[
            row_index,
            column_index
        ]

        
        # MAIN SALINITY AXIS
        ax.set_xlim(
            S_XMIN,
            S_XMAX
        )

        ax.set_ylim(
            PLOT_P_MAX,
            PLOT_P_MIN
        )

        ax.set_xticks(
            S_XTICKS
        )

        ax.set_yticks(
            DEPTH_TICKS
        )

        ax.tick_params(
            axis="y",
            which="major",
            colors="black",
            labelsize=FONT_Y_TICK,
            length=6,
            width=1.2,
            direction="out",
            pad=6
        )

        # SALINITY TICKS — LAST ROW ONLY
        if row_index == n_rows - 1:

            ax.tick_params(
                axis="x",
                which="major",

                bottom=True,
                top=False,

                labelbottom=True,
                labeltop=False,

                colors="blue",
                labelsize=FONT_X_TICK,

                length=6,
                width=1.2,

                direction="out",
                pad=5
            )

        else:
            ax.tick_params(
                axis="x",
                which="both",

                bottom=False,
                top=False,

                labelbottom=False,
                labeltop=False
            )

        # GRID
        
        ax.grid(
            True,
            which="major",
            linewidth=GRID_LINEWIDTH,
            linestyle="-",
            color="0.75",
            alpha=0.30,
            zorder=0
        )

        ax.set_axisbelow(True)

        # MAIN-AXIS BORDERS
        for spine in ax.spines.values():

            spine.set_linewidth(
                SPINE_LINEWIDTH
            )

            spine.set_color(
                "black"
            )

        # SECONDARY TEMPERATURE AXIS
        
        ax_top = ax.twiny()

        ax_top.set_xlim(
            T_XMIN,
            T_XMAX
        )

        ax_top.set_xticks(
            T_XTICKS
        )

        # Keep background transparent.
        ax_top.patch.set_visible(False)

        ax_top.spines["top"].set_linewidth(
            SPINE_LINEWIDTH
        )

        ax_top.spines["top"].set_color(
            "black"
        )

        
        # TEMPERATURE TICKS — FIRST ROW ONLY
        if row_index == 0:

            ax_top.tick_params(
                axis="x",
                which="major",

                top=True,
                bottom=False,

                labeltop=True,
                labelbottom=False,

                colors="red",
                labelsize=FONT_X_TICK,

                length=6,
                width=1.2,

                direction="out",
                pad=5
            )

        else:
            ax_top.tick_params(
                axis="x",
                which="both",

                top=False,
                bottom=False,

                labeltop=False,
                labelbottom=False
            )

            # Remove unnecessary top spine in inner rows.
            ax_top.spines["top"].set_visible(False)

        # RETRIEVE SEASONAL MEAN PROFILE
        profile_output = collect_mean_profiles(
            sea_code,
            season_key
        )

        if profile_output is not None:

            (
                mean_salinity,
                mean_temperature,
                depth,
                n_profiles
            ) = profile_output

            # SALINITY PROFILE
            ax.plot(
                mean_salinity,
                depth,

                color="blue",
                linewidth=PROFILE_LINEWIDTH,

                solid_capstyle="round",
                solid_joinstyle="round",

                zorder=4
            )

            # TEMPERATURE PROFILE
            ax_top.plot(
                mean_temperature,
                depth,

                color="red",
                linewidth=PROFILE_LINEWIDTH,

                solid_capstyle="round",
                solid_joinstyle="round",

                zorder=5
            )

            # AUTOMATIC LEGEND LOCATION
            legend_location, legend_anchor = (
                choose_legend_position(
                    mean_salinity,
                    mean_temperature,
                    depth
                )
            )

            # LEGEND HANDLES
            legend_handles = [
                Patch(
                    facecolor="blue",
                    edgecolor="blue",
                    label="Salinity"
                ),

                Patch(
                    facecolor="red",
                    edgecolor="red",
                    label="Temperature"
                ),

                Line2D(
                    [],
                    [],
                    linestyle="none",
                    marker=None,
                    color="none",
                    label=f"Nprof = {n_profiles:,}"
                )
            ]

            
            # LEGEND / TEXT BOX
            legend = ax.legend(
                handles=legend_handles,

                loc=legend_location,
                bbox_to_anchor=legend_anchor,
                bbox_transform=ax.transAxes,

                fontsize=FONT_LEGEND,

                frameon=True,
                fancybox=False,
                shadow=False,

                facecolor="white",
                edgecolor="0.30",
                framealpha=1.0,

                borderpad=0.58,
                labelspacing=0.32,

                handlelength=1.25,
                handleheight=0.90,
                handletextpad=0.52,

                borderaxespad=0.0,
                columnspacing=0.40,

                ncol=1
            )

            # Keep legend above temperature and salinity lines.
            legend.set_zorder(100)

            legend_frame = legend.get_frame()

            legend_frame.set_linewidth(
                1.0
            )

            legend_frame.set_facecolor(
                "white"
            )

            legend_frame.set_alpha(
                1.0
            )

            for legend_text in legend.get_texts():

                legend_text.set_fontsize(
                    FONT_LEGEND
                )

                legend_text.set_fontweight(
                    "normal"
                )

                legend_text.set_color(
                    "black"
                )

                legend_text.set_clip_on(
                    False
                )

        else:
            # NO-DATA TEXT BOX
            ax.text(
                0.50,
                0.50,
                "No data",

                transform=ax.transAxes,

                ha="center",
                va="center",

                fontsize=FONT_NO_DATA,
                fontweight="bold",
                color="dimgray",

                bbox={
                    "boxstyle": "round,pad=0.42",
                    "facecolor": "white",
                    "edgecolor": "0.35",
                    "linewidth": 1.0,
                    "alpha": 1.0
                },

                zorder=100,
                clip_on=True
            )

        
        # FIRST-ROW SEASON AND TEMPERATURE LABELS
        if row_index == 0:

            ax.set_title(
                SEASON_TITLE[season_key],

                fontsize=FONT_SEASON_TITLE,
                fontweight="bold",
                color="black",

                # Smaller pad reduces the gap between the
                # season title and temperature-axis label.
                pad=30
            )

            ax_top.set_xlabel(
                "Temperature (°C)",

                color="red",
                fontsize=FONT_TEMP_LABEL,
                fontweight="bold",

                # Small label padding keeps it near the ticks.
                labelpad=5
            )

        # LAST-ROW SALINITY LABELS
        if row_index == n_rows - 1:

            ax.set_xlabel(
                "Salinity (psu)",

                color="blue",
                fontsize=FONT_SAL_LABEL,
                fontweight="bold",

                labelpad=10
            )

        # FIRST-COLUMN DEPTH LABELS
        if column_index == 0:

            ax.set_ylabel(
                "Depth (m)",

                fontsize=FONT_DEPTH_LABEL,
                fontweight="bold",
                color="black",

                labelpad=12
            )

            ax.tick_params(
                axis="y",
                labelleft=True
            )

        else:
            ax.tick_params(
                axis="y",
                labelleft=False
            )

    # 17) LEFT-SIDE SEA ABBREVIATION
    axes[row_index, 0].text(
        -0.300,
        0.50,
        sea_code,

        transform=axes[row_index, 0].transAxes,

        rotation=90,
        rotation_mode="anchor",

        ha="center",
        va="center",

        fontsize=FONT_SEA_LABEL,
        fontweight="bold",
        color="black",

        clip_on=False,
        zorder=200
    )


# 18) SAVE FINAL FIGURE
fig.savefig(
    out_file,

    dpi=DPI_SAVE,
    bbox_inches="tight",

    facecolor="white",
    edgecolor="none",

    pad_inches=0.18
)

plt.show()
plt.close(fig)

print("\nFigure saved successfully:")
print(out_file)
