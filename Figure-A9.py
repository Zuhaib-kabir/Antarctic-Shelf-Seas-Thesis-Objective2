# Obj1Fig13_SO
# Seasonal T–S Diagrams with Water-Mass Classification
# First Seven Antarctic Marginal Seas
#
# FINAL COMPACT LAYOUT
# Salinity tick labels       : last row only
# Temperature tick labels    : first column only
# Season labels              : 22 pt
# Sea labels                 : 22 pt
# Axis labels                : 18 pt
# X-axis tick labels         : 17 pt
# Y-axis tick labels         : 18 pt
# Nprof text                 : 16 pt
# Top legend text            : 19 pt
#
# Improvements:
#   - Reduced gap between top legend and season labels
#   - Larger legend text
#   - Nprof box moved downward
#   - Upper density labels remain more visible
#   - Compact horizontal and vertical panel spacing

#Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')


# 0) INSTALL REQUIRED PACKAGES — GOOGLE COLAB

!pip -q install netCDF4 gsw numpy pandas matplotlib



# 1) IMPORTS

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gsw

from netCDF4 import Dataset
from datetime import datetime, timedelta
from matplotlib.lines import Line2D



# 2) INPUT AND OUTPUT PATHS

ARGO_NC = (
    "/content/drive/MyDrive/SAM_Thesis/Data/"
    "argo_SO_profiles_2001_2025_cleaned_gridded.nc"
)

out_dir = "/content/drive/MyDrive/SAM_Thesis/Fig"
os.makedirs(out_dir, exist_ok=True)

OUT_PNG = os.path.join(
    out_dir,
    "Obj1Fig13_SO_Seasonal_TS_WaterMass_Final_Compact.png"
)



# 3) ANALYSIS AND PLOTTING SETTINGS

YEAR_START = 2008
YEAR_END = 2025

LAT_MIN = -90.0
LAT_MAX = -60.0

PLOT_P_MIN = 0.0
PLOT_P_MAX = 2000.0

MAX_POINTS_PER_PANEL = 70000

# Grid used for sigma-theta contours
S_GRID = np.linspace(
    33.0,
    35.2,
    170
)

T_GRID = np.linspace(
    -2.2,
    8.0,
    170
)

DPI_SAVE = 1080

# Large enough for 7 × 4 panels while retaining compact spacing
FIGSIZE = (27, 31)



# 4) FONT AND GRAPHICAL SETTINGS
FONT_SEASON_TITLE = 22
FONT_SEA_LABEL = 22

FONT_X_LABEL = 18
FONT_Y_LABEL = 18

FONT_X_TICK = 17
FONT_Y_TICK = 18

FONT_NPROF = 16
FONT_NO_DATA = 16

# Increased legend font size
FONT_LEGEND = 19

FONT_CONTOUR = 10

SCATTER_SIZE = 3.0
SCATTER_ALPHA = 0.65

SPINE_LINEWIDTH = 1.2
GRID_LINEWIDTH = 0.65
CONTOUR_LINEWIDTH = 0.75



# 5) GLOBAL MATPLOTLIB SETTINGS

plt.rcParams.update({
    "font.family": "DejaVu Sans",

    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
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
    ("WED", "Weddell Sea",          "60°W–20°W", -60, -20),
    ("KHV", "King Haakon VII Sea",  "20°W–0°",   -20,   0),
    ("RLS", "Riiser-Larsen Sea",    "0°–10°E",     0,  10),
    ("LAZ", "Lazarev Sea",          "10°E–30°E",   10,  30),
    ("COS", "Cosmonauts Sea",       "30°E–50°E",   30,  50),
    ("COO", "Cooperation Sea",      "50°E–70°E",   50,  70),
    ("DAV", "Davis Sea",            "70°E–90°E",   70,  90),
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
    lon = np.asarray(
        lon,
        dtype=float
    )

    return (
        (lon + 180.0) % 360.0
    ) - 180.0


def lon_in_range(
    lon_arr,
    lon_min,
    lon_max
):
    """
    Select longitude values inside a specified sector.

    Supports:
      1. Standard longitude ranges
      2. Dateline-crossing longitude ranges
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



# 9) WATER-MASS CATEGORIES

WM_ORDER = [
    "AASW",
    "WW",
    "mCDW",
    "CDW",
    "AABW",
]

WM_FULL = {
    "AASW": "Antarctic Surface Water (AASW)",
    "WW": "Winter Water (WW)",
    "mCDW": "Modified Circumpolar Deep Water (mCDW)",
    "CDW": "Circumpolar Deep Water (CDW)",
    "AABW": "Antarctic Bottom Water (AABW)",
}

WM_COLOR = {
    "AASW": "cyan",
    "WW": "blue",
    "mCDW": "red",
    "CDW": "purple",
    "AABW": "green",
}


# 10) WATER-MASS CLASSIFICATION FUNCTION

def classify_watermass(
    salinity,
    temperature,
    pressure
):
    """
    Simplified water-mass classification for Antarctic
    shelf-sea temperature–salinity diagrams.
    """

    # Antarctic Bottom Water
    if (
        pressure >= 1500 and
        temperature <= 0.5
    ):
        return "AABW"

    # Surface waters
    if pressure < 200:

        if temperature <= -0.5:
            return "WW"

        return "AASW"

    # Subsurface and deep waters
    if pressure >= 200:

        if temperature >= 1.0:
            return "CDW"

        return "mCDW"

    return "mCDW"



# 11) READ ARGO DATA

print("Reading Argo dataset:")
print(ARGO_NC)

with Dataset(
    ARGO_NC,
    mode="r"
) as ds:

    P = np.asarray(
        ds.variables["PRES_GRID"][:],
        dtype=np.float64
    )

    LAT = np.asarray(
        ds.variables["LATITUDE"][:],
        dtype=np.float64
    )

    LON = np.asarray(
        ds.variables["LONGITUDE"][:],
        dtype=np.float64
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

LON = norm_lon(LON)

print("Argo dataset loaded successfully.")



# 12) CONVERT ARGO JULIAN DAYS TO DATETIME
argo_base_date = datetime(
    1950,
    1,
    1
)

time = np.full(
    JULD.shape,
    np.datetime64("NaT"),
    dtype="datetime64[ns]"
)

valid_time_mask = np.isfinite(
    JULD
)

time[valid_time_mask] = np.array([
    np.datetime64(
        argo_base_date +
        timedelta(days=float(day))
    )
    for day in JULD[valid_time_mask]
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


# 13) BASIC PROFILE AND PRESSURE MASKS
mask_base = (
    np.isfinite(LAT) &
    np.isfinite(LON) &
    (LAT >= LAT_MIN) &
    (LAT <= LAT_MAX) &
    (years >= YEAR_START) &
    (years <= YEAR_END) &
    pd.notna(seasons)
)

print(
    "Profiles in selected latitude band and period:",
    int(mask_base.sum()),
    "/",
    len(mask_base)
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


# 14) COLLECT AND CLASSIFY PANEL DATA
def collect_and_classify(
    sea_code,
    season_key
):
    """
    Collect valid temperature, salinity and pressure values
    for one sea-season panel and classify each point into a
    water-mass category.
    """

    lon_min, lon_max = sea_ranges[
        sea_code
    ]

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

    pressure_2d = np.broadcast_to(
        P_selected[None, :],
        temperature_subset.shape
    )

    temperature_flat = temperature_subset.reshape(
        -1
    )

    salinity_flat = salinity_subset.reshape(
        -1
    )

    pressure_flat = pressure_2d.reshape(
        -1
    )

    valid_data = (
        np.isfinite(temperature_flat) &
        np.isfinite(salinity_flat) &
        np.isfinite(pressure_flat)
    )

    # Basic physical-range quality control
    valid_data &= (
        (temperature_flat > -3.0) &
        (temperature_flat < 20.0)
    )

    valid_data &= (
        (salinity_flat > 0.0) &
        (salinity_flat < 42.0)
    )

    temperature_flat = temperature_flat[
        valid_data
    ]

    salinity_flat = salinity_flat[
        valid_data
    ]

    pressure_flat = pressure_flat[
        valid_data
    ]

    if temperature_flat.size < 200:
        return None

    # Reproducible random subsampling
    if temperature_flat.size > MAX_POINTS_PER_PANEL:

        random_seed = (
            1000 +
            sea_codes.index(sea_code) * 10 +
            SEASON_ORDER.index(season_key)
        )

        rng = np.random.default_rng(
            random_seed
        )

        selected_indices = rng.choice(
            temperature_flat.size,
            size=MAX_POINTS_PER_PANEL,
            replace=False
        )

        temperature_flat = temperature_flat[
            selected_indices
        ]

        salinity_flat = salinity_flat[
            selected_indices
        ]

        pressure_flat = pressure_flat[
            selected_indices
        ]

    watermass_labels = np.array(
        [
            classify_watermass(
                salinity,
                temperature,
                pressure
            )
            for salinity, temperature, pressure in zip(
                salinity_flat,
                temperature_flat,
                pressure_flat
            )
        ],
        dtype=object
    )

    return (
        salinity_flat,
        temperature_flat,
        pressure_flat,
        watermass_labels,
        profile_indices.size
    )


# 15) SIGMA-THETA CONTOURS
SS, TT = np.meshgrid(
    S_GRID,
    T_GRID
)

reference_pressure = np.zeros_like(
    SS
)

absolute_salinity = gsw.SA_from_SP(
    SS,
    reference_pressure,
    0.0,
    -70.0
)

conservative_temperature = gsw.CT_from_t(
    absolute_salinity,
    TT,
    reference_pressure
)

SIGMA0 = gsw.sigma0(
    absolute_salinity,
    conservative_temperature
)

SIGMA_LEVELS = np.arange(
    26.0,
    28.6,
    0.2
)


# 16) AXIS SETTINGS
XMIN = 33.0
XMAX = 35.2

YMIN = -2.2
YMAX = 8.0

XTICKS = np.arange(
    33.0,
    35.201,
    0.25
)

YTICKS = np.arange(
    -2.0,
    8.1,
    2.0
)


# 17) CREATE FIGURE
n_rows = len(
    sea_codes
)

n_columns = len(
    SEASON_ORDER
)

fig, axes = plt.subplots(
    nrows=n_rows,
    ncols=n_columns,
    figsize=FIGSIZE,
    sharex=True,
    sharey=True,
    squeeze=False
)

# First-row panels moved upward to reduce the gap between
# the top legend and the seasonal headings.
fig.subplots_adjust(
    left=0.115,
    right=0.980,

    # Increased from 0.865 to 0.895
    top=0.895,

    bottom=0.070,

    wspace=0.14,
    hspace=0.13
)


# 18) DRAW ALL SEA-SEASON PANELS
for row_index, sea_code in enumerate(
    sea_codes
):

    for column_index, season_key in enumerate(
        SEASON_ORDER
    ):

        ax = axes[
            row_index,
            column_index
        ]

        # AXIS LIMITS AND TICKS
        ax.set_xlim(
            XMIN,
            XMAX
        )

        ax.set_ylim(
            YMIN,
            YMAX
        )

        ax.set_xticks(
            XTICKS
        )

        ax.set_yticks(
            YTICKS
        )

        # DENSITY CONTOURS
        contour_set = ax.contour(
            SS,
            TT,
            SIGMA0,

            levels=SIGMA_LEVELS,

            linewidths=CONTOUR_LINEWIDTH,
            alpha=0.85,

            zorder=1
        )

        ax.clabel(
            contour_set,

            inline=True,
            inline_spacing=3,

            fontsize=FONT_CONTOUR,

            fmt="%.1f"
        )

        # RETRIEVE AND CLASSIFY PANEL DATA
        panel_output = collect_and_classify(
            sea_code,
            season_key
        )

        if panel_output is None:

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
                    "boxstyle": "round,pad=0.40",
                    "facecolor": "white",
                    "edgecolor": "0.40",
                    "linewidth": 1.0,
                    "alpha": 1.0
                },

                zorder=100,
                clip_on=True
            )

        else:

            (
                salinity_flat,
                temperature_flat,
                pressure_flat,
                watermass_labels,
                n_profiles
            ) = panel_output

            
            # WATER-MASS SCATTER POINTS
            for watermass in WM_ORDER:

                watermass_mask = (
                    watermass_labels == watermass
                )

                if np.any(
                    watermass_mask
                ):

                    ax.scatter(
                        salinity_flat[watermass_mask],
                        temperature_flat[watermass_mask],

                        s=SCATTER_SIZE,
                        alpha=SCATTER_ALPHA,

                        c=WM_COLOR[watermass],

                        edgecolors="none",
                        linewidths=0,

                        rasterized=True,
                        zorder=3
                    )

            # NPROF TEXT INSIDE PANEL
            #
            # Moved downward from y=0.965 to y=0.925 so the
            # upper sigma-density contour labels remain visible.
            
            ax.text(
                0.025,
                0.885,
                f"Nprof = {n_profiles:,}",

                transform=ax.transAxes,

                ha="left",
                va="top",

                fontsize=FONT_NPROF,
                fontweight="bold",
                color="black",

                bbox={
                    "boxstyle": "square,pad=0.20",
                    "facecolor": "white",
                    "edgecolor": "0.65",
                    "linewidth": 0.7,
                    "alpha": 0.92
                },

                zorder=50,
                clip_on=True
            )

        
        # GRID
        ax.grid(
            True,
            which="major",

            linewidth=GRID_LINEWIDTH,
            linestyle="-",

            color="0.75",
            alpha=0.22,

            zorder=0
        )

        ax.set_axisbelow(
            True
        )

        # PANEL BORDERS
        for spine in ax.spines.values():

            spine.set_linewidth(
                SPINE_LINEWIDTH
            )

            spine.set_color(
                "black"
            )

        # X TICKS AND SALINITY LABEL — LAST ROW ONLY
        if row_index == n_rows - 1:

            ax.tick_params(
                axis="x",
                which="major",

                bottom=True,
                top=False,

                labelbottom=True,
                labeltop=False,

                labelsize=FONT_X_TICK,

                length=6,
                width=1.2,

                direction="out",
                pad=5
            )

            for tick_label in ax.get_xticklabels():

                tick_label.set_rotation(
                    45
                )

                tick_label.set_ha(
                    "right"
                )

                tick_label.set_rotation_mode(
                    "anchor"
                )

                tick_label.set_fontweight(
                    "bold"
                )

            ax.set_xlabel(
                "Salinity",

                fontsize=FONT_X_LABEL,
                fontweight="bold",
                color="black",

                labelpad=11
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

        
        # Y TICKS AND TEMPERATURE LABEL — FIRST COLUMN ONLY
        if column_index == 0:

            ax.tick_params(
                axis="y",
                which="major",

                left=True,
                right=False,

                labelleft=True,
                labelright=False,

                labelsize=FONT_Y_TICK,

                length=6,
                width=1.2,

                direction="out",
                pad=6
            )

            for tick_label in ax.get_yticklabels():

                tick_label.set_fontweight(
                    "bold"
                )

            ax.set_ylabel(
                "Temperature (°C)",

                fontsize=FONT_Y_LABEL,
                fontweight="bold",
                color="black",

                labelpad=12
            )

        else:

            ax.tick_params(
                axis="y",
                which="both",

                left=False,
                right=False,

                labelleft=False,
                labelright=False
            )

        
        # SEASON TITLES — FIRST ROW ONLY
          if row_index == 0:

            ax.set_title(
                SEASON_TITLE[season_key],

                fontsize=FONT_SEASON_TITLE,
                fontweight="bold",
                color="black",

                pad=12
            )

    
    # 19) LEFT-SIDE SEA ABBREVIATION
    axes[row_index, 0].text(
        -0.305,
        0.50,
        sea_code,

        transform=axes[
            row_index,
            0
        ].transAxes,

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



# 20) TOP WATER-MASS LEGEND
legend_handles = [
    Line2D(
        [0],
        [0],

        marker="o",
        linestyle="None",

        markersize=12,

        markerfacecolor=WM_COLOR[watermass],
        markeredgecolor=WM_COLOR[watermass],

        label=WM_FULL[watermass]
    )
    for watermass in WM_ORDER
]

top_legend = fig.legend(
    handles=legend_handles,

    loc="upper center",

    # Reduced from 0.972 to 0.955 to move legend downward
    bbox_to_anchor=(0.50, 0.955),

    ncol=3,

    fontsize=FONT_LEGEND,

    frameon=True,
    fancybox=False,
    shadow=False,

    facecolor="white",
    edgecolor="0.45",
    framealpha=1.0,

    borderpad=0.65,
    labelspacing=0.45,

    handlelength=1.20,
    handletextpad=0.50,

    columnspacing=1.25
)

top_legend.set_zorder(
    300
)

legend_frame = top_legend.get_frame()

legend_frame.set_linewidth(
    1.0
)

legend_frame.set_facecolor(
    "white"
)

legend_frame.set_alpha(
    1.0
)

for legend_text in top_legend.get_texts():

    legend_text.set_fontsize(
        FONT_LEGEND
    )

    legend_text.set_fontweight(
        "bold"
    )

    legend_text.set_color(
        "black"
    )


# 21) SAVE FINAL FIGURE
fig.savefig(
    OUT_PNG,

    dpi=DPI_SAVE,
    bbox_inches="tight",

    facecolor="white",
    edgecolor="none",

    pad_inches=0.18
)

plt.show()
plt.close(fig)

print("\nFigure saved successfully:")
print(OUT_PNG)
