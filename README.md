# Antarctic Shelf Seas Thesis — Objective 2

This repository contains the Python workflows used for **Objective 2 of the Antarctic shelf-seas thesis**, focusing on **Southern Ocean hydrography, steric sea-level variability, observed sea-level anomaly (SLA), freshwater and salinity structure, upper-ocean stability, trends, EOF variability, local forcing, and AAO/SAM attribution**.

The code is organized as a set of main-figure and appendix/supplementary-figure workflows. Most analyses cover **2008–2025** and are designed primarily for execution in **Google Colab**, with input and output files stored in Google Drive.

---

## Objective 2 analysis scope

The repository combines Argo, satellite, reanalysis, bathymetric, and EN4-derived products to examine:

- seasonal temperature–salinity structure;
- water-mass characteristics across Antarctic marginal seas;
- total, thermosteric, and halosteric sea-level components;
- satellite SLA versus EN4 total steric height;
- steric decomposition and closure;
- freshwater content and surface salinity;
- upper-ocean stratification and mixed-layer depth;
- sea-wise monthly climatological cycles;
- study-period sea-level trends;
- sector-based EOF modes and principal components;
- atmospheric and oceanic forcing;
- AAO/SAM relationships;
- Argo observational coverage;
- Argo–EN4 steric validation;
- 0–1000 m versus 0–2000 m integration-depth sensitivity.

Most workflows use the same **13 Antarctic shelf-sea sectors**:

1. Weddell Sea (WED)
2. King Haakon VII Sea (KHV)
3. Riiser-Larsen Sea (RLS)
4. Lazarev Sea (LAZ)
5. Cosmonauts Sea (COS)
6. Cooperation Sea (COO)
7. Davis Sea (DAV)
8. Mawson Sea (MAW)
9. D'Urville Sea (DUR)
10. Somov Sea (SOM)
11. Ross Sea (ROS)
12. Amundsen Sea (AMU)
13. Bellingshausen Sea (BEL)

The Ross Sea is treated as a date-line-crossing sector where required.

---

## Repository contents

### Main figure scripts

| Script | Main purpose |
|---|---|
| [`Figure-10.py`](./Figure-10.py) | Seasonal Argo temperature–salinity profiles for the first seven Antarctic marginal seas |
| [`Figure-11.py`](./Figure-11.py) | Seasonal Argo temperature–salinity profiles for the remaining six Antarctic marginal seas |
| [`Figure-12.py`](./Figure-12.py) | Quantitative integration of T–S-derived hydrographic properties and steric variability across 13 shelf seas |
| [`Figure-13.py`](./Figure-13.py) | Observed SLA versus EN4 total steric-height seasonal anomalies |
| [`Figure-14.py`](./Figure-14.py) | EN4 seasonal total/thermosteric/halosteric decomposition and closure |
| [`Figure-15.py`](./Figure-15.py) | EN4 freshwater storage, surface-salinity anomaly, stratification, MLD, and SIC |
| [`Figure-16.py`](./Figure-16.py) | Sea-wise monthly climatological cycles of SLA and steric components |
| [`Figure-17.py`](./Figure-17.py) | 2008–2025 SLA and steric study-period trends |
| [`Figure-18.py`](./Figure-18.py) | Sector-based EOF analysis of SLA and total steric variability |
| [`Figure-19.py`](./Figure-19.py) | Local forcing and AAO attribution across the 13 Antarctic shelf seas |

### Appendix / supplementary scripts

| Script | Main purpose |
|---|---|
| [`Figure-A5.py`](./Figure-A5.py) | Argo observational coverage and profile-depth support |
| [`Figure-A6.py`](./Figure-A6.py) | Seasonal steric components using local objective mapping |
| [`Figure-A7.py`](./Figure-A7.py) | Argo–EN4 0–1000 m steric validation |
| [`Figure-A8.py`](./Figure-A8.py) | EN4 0–1000 m versus 0–2000 m integration-depth sensitivity |
| [`Figure-A9.py`](./Figure-A9.py) | Seasonal T–S diagrams and water-mass classification for the first seven seas |
| [`Figure-A10.py`](./Figure-A10.py) | Seasonal T–S diagrams and water-mass classification for the remaining six seas |

---

# Main figure workflows

## `Figure-10.py` — Seasonal vertical T–S profiles: first seven seas

This script analyzes seasonal vertical **temperature and salinity profiles** from the cleaned/gridded Argo dataset for:

```text
WED, KHV, RLS, LAZ, COS, COO, DAV
```

### Main settings

```text
Period:       2008–2025
Latitude:     90°S–60°S
Pressure:     0–2000 dbar
Seasons:      SON, DJF, MAM, JJA
```

The script:

- reads `PRES_GRID`, latitude, longitude, temperature, salinity, and Argo time information;
- normalizes longitude to `[-180°, 180°)`;
- assigns profiles to Southern Hemisphere seasons;
- assigns profiles to the appropriate Antarctic sea sector;
- calculates seasonal mean vertical temperature and salinity structure;
- produces a compact multi-row seasonal profile layout.

Main Argo input:

```text
argo_SO_profiles_2001_2025_cleaned_gridded.nc
```

---

## `Figure-11.py` — Seasonal vertical T–S profiles: second sea group

This workflow applies the same seasonal vertical-profile approach to:

```text
MAW, DUR, SOM, ROS, AMU, BEL
```

The Ross Sea date-line crossing is explicitly handled.

The analysis uses:

```text
Period:       2008–2025
Latitude:     90°S–60°S
Pressure:     0–2000 dbar
```

and the same Southern Hemisphere seasonal definitions used in `Figure-10.py`.

---

## `Figure-12.py` — T–S-derived hydrography and steric variability

This is a complete quantitative integration workflow across all **13 Antarctic shelf seas**.

It combines Argo hydrographic information with processed steric and freshwater products.

### Main inputs

```text
argo_SO_profiles_2001_2025_cleaned_gridded.nc
thermosteric_height_monthly_2008_2025_SO.nc
halosteric_height_monthly_2008_2025_SO.nc
freshwater_content_monthly_2008_2025_SO.nc
```

### Main processing

The script:

- reads and quality-controls Argo profiles;
- applies the source-defined water-mass criteria;
- calculates sea × season hydrographic metrics;
- calculates seasonal thermosteric, halosteric, and freshwater metrics;
- merges the 13 seas × 4 seasons into **52 analysis records**;
- calculates regression/correlation statistics;
- saves analysis-ready CSV and Excel tables;
- renders the latest corrected seasonal-trajectory figure.

### Main analysis outputs

```text
Fig06_TS_steric_metrics_13seas_seasonal.csv
Fig06_TS_steric_metrics_13seas_seasonal.xlsx
Fig06_regression_statistics.csv
```

---

## `Figure-13.py` — Observed SLA versus EN4 total steric height

This workflow compares satellite-observed **sea-level anomaly (SLA)** with **EN4 total steric-height anomaly integrated over 0–1000 m**.

Both fields are represented on the native EN4 1° grid.

### Scientific design

- study period: **2008–2025**;
- Southern Ocean domain south of 60°S;
- fixed common SLA–EN4 mask;
- OSTIA ocean mask;
- GEBCO minimum-water-depth criterion of **1000 m**;
- SLA minimum valid fraction of **0.70**;
- EN4 valid-layer support required through the full record;
- no Argo objective mapping in this workflow.

Seasonal anomaly is defined as:

```text
seasonal climatology − all-month mean
```

For each season, the script reports:

- spatial Pearson correlation;
- RMSD;
- mean bias;
- number of common cells;
- spatial coverage.

Bias is defined as:

```text
EN4 total steric − observed SLA
```

### Main outputs

```text
Fig07_SLA_EN4_total_steric_statistics.csv
Fig07_fixed_common_mask_EN4_grid.nc
Fig07_SLA_EN4_seasonal_anomalies_common_mask.nc
```

---

## `Figure-14.py` — EN4 seasonal steric decomposition

This workflow decomposes EN4 steric-height variability into:

- total steric height;
- thermosteric height;
- halosteric height.

The seasonal anomaly is:

```text
seasonal climatology − 2008–2025 all-month mean
```

### Steric closure

The script evaluates:

```text
residual = total − (thermosteric + halosteric)
```

using common finite cells across the component fields.

The source workflow uses a closure tolerance of:

```text
±1 cm
```

and displays the GEBCO **1000 m** contour.

### Main outputs

```text
Fig08_EN4_steric_decomposition_closure_metrics.csv
Fig08_EN4_steric_decomposition_seasonal_anomalies.nc
Fig08_EN4_steric_decomposition_summary.txt
```

---

## `Figure-15.py` — Freshwater storage and upper-ocean stability

This is a **low-RAM / Colab-safe** workflow that reads the monthly source products sequentially rather than loading the full 216-month MLD and SIC arrays simultaneously.

### Figure variables

1. freshwater-content seasonal climatology;
2. surface-salinity seasonal departure from the 2008–2025 all-month mean;
3. stratification seasonal climatology.

The stratification panels also include:

- solid mixed-layer-depth contours;
- dashed **15% SIC** contour.

### Main precomputed product

```text
Fig09_EN4_freshwater_stratification_seasonal_fields.nc
```

The final plotting section reads this seasonal NetCDF rather than recomputing the monthly processing.

---

## `Figure-16.py` — Sea-wise monthly climatological cycles

This workflow calculates Jan–Dec climatological cycles for:

- SLA;
- EN4 total steric height;
- EN4 thermosteric height;
- EN4 halosteric height.

The output is calculated for the **13 Antarctic shelf seas**.

### Method

- interpolate SLA to the EN4 grid;
- build one fixed SLA–EN4 common mask;
- calculate area-weighted monthly sea-wise time series;
- calculate Jan–Dec climatologies;
- express each month as a departure from that sea's annual mean;
- track SLA coverage;
- identify low-coverage cells.

### Main settings

```text
Study period:                         2008–2025
SLA core-validity threshold:          0.70
Low-coverage SLA threshold:           0.70
Minimum valid years per climatology:  5
Output units:                         cm
```

### Main analysis output

```text
Fig10_sea_wise_monthly_climatological_cycles.csv
```

---

## `Figure-17.py` — Study-period trends

This script calculates 2008–2025 trends in:

- SLA;
- total steric height;
- thermosteric height;
- halosteric height.

### Trend method

1. interpolate SLA to the EN4 grid;
2. apply one fixed SLA–EN4 common mask;
3. remove each grid cell's calendar-month climatology;
4. calculate annual deseasonalized anomalies;
5. calculate **Sen slopes**;
6. test significance with the **Hamed–Rao modified Mann–Kendall test**;
7. apply **Benjamini–Hochberg FDR correction**;
8. calculate sea-wise trends;
9. calculate **95% confidence intervals** with a 3-year circular moving-block bootstrap.

Trend units are:

```text
cm decade⁻¹
```

### Main analysis outputs

```text
Fig11_gridcell_trend_statistics_2008_2025.nc
Fig11_13sea_Sen_slopes_bootstrap_CI.csv
Fig11_deseasonalized_annual_anomalies_common_grid.nc
```

---

## `Figure-18.py` — Sector-based EOF analysis

This workflow analyzes monthly SLA and total steric variability using **13-sector EOF analysis**.

### Input matrices

```text
216 months × 13 Antarctic sectors
```

### Preprocessing

- remove each sector's monthly climatology;
- detrend each sector.

### EOF design

- **covariance EOFs** are the main analysis;
- **correlation EOFs** are calculated as a sensitivity analysis;
- the first **three EOF modes** are retained;
- principal components are standardized;
- 13-month running-mean PCs are also calculated.

### Numerical outputs

```text
Fig12_SLA_sector_matrix_216x13.csv
Fig12_steric_sector_matrix_216x13.csv
Fig12_EOF_sector_loadings_summary.csv
Fig12_EOF_PC_timeseries_summary.csv
Fig12_correlation_EOF_sensitivity_summary.csv
```

The final plotting stage uses those precomputed outputs rather than recalculating the EOF solution.

---

## `Figure-19.py` — Local forcing and AAO attribution

This is the repository's large **low-RAM / checkpointed attribution workflow**.

It reduces monthly gridded datasets to 13 area-weighted sea-sector time series and examines the relationships among steric sea level, local forcing, and AAO/SAM variability.

### Low-RAM strategy

The script:

- processes one NetCDF variable at a time;
- loads one monthly 2-D field at a time;
- immediately reduces fields to 13 sector means;
- closes each large dataset before moving to the next;
- saves small checkpoint CSV files;
- reuses cached sector/composite products when possible.

### Surface heat-flux proxy

The source-defined proxy is:

```text
Qnet* = SSRD + STRD + signed(SLHF) + signed(SSHF)
```

The script explicitly labels this as a **surface heat-flux proxy**, because the available radiation fields are downward SSRD and STRD rather than complete net shortwave and net longwave terms.

### Additional processing

The workflow includes:

- wind-stress curl;
- SST;
- OHC;
- MLD;
- SIC;
- freshwater content;
- surface-salinity anomaly;
- snowfall;
- stratification;
- AAO/SAM;
- deseasonalized correlations;
- regression coefficients;
- lag correlations;
- low-RAM AAO–SLA composites;
- online bootstrap sign counts.

### Main processed outputs

```text
Fig13_sector_driver_timeseries_LOW_RAM.csv
Fig13a_SLA_steric_correlations.csv
Fig13b_thermosteric_driver_correlations.csv
Fig13c_halosteric_driver_correlations.csv
Fig13d_thermosteric_regression_coefficients.csv
Fig13e_halosteric_regression_coefficients.csv
Fig13g_AAO_SLA_lag_correlations.csv
Fig13h_AAO_SLA_composite_LOW_RAM.nc
```

---

# Appendix / supplementary workflows

## `Figure-A5.py` — Argo coverage

This script characterizes Southern Ocean Argo sampling during **2001–2025**.

### Panels

- Spring profile-count map;
- Summer profile-count map;
- Autumn profile-count map;
- Winter profile-count map;
- sea × season profile-count heatmap;
- percentage of profiles reaching 1000 and 2000 dbar.

### Corrected depth-support calculation

For each profile, maximum pressure is determined from the deepest `PRES_GRID` level where **both TEMP and PSAL are finite**.

The workflow classifies profiles by whether they reach:

```text
1000 dbar
2000 dbar
```

and saves seasonal profile-count maps plus sea-wise support tables.

---

## `Figure-A6.py` — Seasonal steric components with objective mapping

This workflow displays seasonal climatologies of:

- total steric height;
- thermosteric height;
- halosteric height.

It applies **local objective mapping / optimal interpolation** to the Argo-derived steric products.

The code uses a Gaussian distance-based mapping framework with explicit settings for:

- length scale;
- search radius;
- minimum/maximum neighbors;
- noise-to-signal ratio.

Small black dots identify original valid Argo grid cells.

---

## `Figure-A7.py` — Argo–EN4 0–1000 m validation

This script performs direct profile-to-grid validation of Argo-derived and EN4 steric quantities.

### Validation panels

- Argo vs EN4 total steric height;
- Argo vs EN4 thermosteric height;
- Argo vs EN4 halosteric height;
- sea-wise mean bias;
- seasonal Pearson correlation;
- number of matched Argo profiles.

### Low-RAM / resumable design

- large NetCDF inputs can be staged to `/content`;
- `h5netcdf` is attempted before `netCDF4`;
- Argo T/S arrays are processed in contiguous chunks;
- fixed EN4 reference profiles are cached;
- completed profile chunks are cached as CSV files;
- the completed matched-profile table can be reused.

### Key matching settings

```text
Period:                         2008–2025
Layer:                          0–1000 m
Argo profile chunk:             1500 profiles
Maximum interpolation gap:      100 m
Maximum EN4-cell distance:      120 km
Maximum accepted top depth:     15 m
Required bottom support:        1000 m
Temperature guard:              −3.5 to 15 °C
Salinity guard:                 20 to 40
```

The steric calculations use the source workflow's TEOS-10 decomposition.

---

## `Figure-A8.py` — EN4 integration-depth sensitivity

This workflow compares steric variability integrated over:

```text
0–1000 m
0–2000 m
```

for:

- total steric height;
- thermosteric height;
- halosteric height.

### Fixed common spatial rule

Every comparison uses the same spatial mask requiring:

- GEBCO depth ≥2000 m;
- valid 0–1000 m total/thermo/halo data;
- valid 0–2000 m total/thermo/halo data;
- the configured valid-data fraction through the full study period.

The source setting is:

```text
COMMON_VALID_FRACTION = 1.00
```

so all 216 months must be valid for inclusion.

The difference is defined as:

```text
0–2000 m − 0–1000 m
```

The default diagnostic is:

```text
Winter climatology − Summer climatology
```

If the 0–2000 m products do not already exist, the script can create them month by month from raw EN4 temperature and salinity using the same fixed 2008–2025 grid-cell reference-state method.

---

## `Figure-A9.py` — Seasonal T–S diagrams: first seven seas

This workflow creates seasonal temperature–salinity diagrams with water-mass classification for:

```text
WED, KHV, RLS, LAZ, COS, COO, DAV
```

### Main features

- Argo profiles from 2008–2025;
- data to 2000 dbar;
- Southern Hemisphere seasons;
- TEOS-10 calculations using `gsw`;
- sigma-theta density contours;
- water-mass classification;
- profile-count annotation;
- compact multi-panel layout.

---

## `Figure-A10.py` — Seasonal T–S diagrams: second sea group

This script applies the same T–S/water-mass analysis to:

```text
MAW, DUR, SOM, ROS, AMU, BEL
```

The Ross Sea date-line crossing is explicitly handled.

---

# Principal datasets referenced by the repository

The workflows use a combination of observational, gridded, and reanalysis products, including:

- cleaned/gridded Argo hydrographic profiles;
- satellite-observed SLA;
- EN4 temperature and salinity;
- EN4 total steric height;
- EN4 thermosteric height;
- EN4 halosteric height;
- EN4 freshwater content;
- EN4 surface salinity and salinity anomaly;
- EN4 stratification;
- EN4 ocean heat content;
- OSTIA sea-ice fraction;
- OSTIA sea-surface temperature;
- mixed-layer depth;
- GEBCO bathymetry;
- ERA5 radiation and turbulent heat-flux variables;
- ERA5 wind;
- ERA5 snowfall;
- NOAA CPC SAM/AAO index.

Large observational and gridded data files are referenced through local/Google Drive paths and are not stored in this repository.

---

## Google Drive directory structure

The scripts primarily reference:

```text
/content/drive/MyDrive/SAM_Thesis/Data/
/content/drive/MyDrive/SAM_Thesis/Processed/
/content/drive/MyDrive/SAM_Thesis/Processed/EN4_NetCDF_inventory/
/content/drive/MyDrive/SAM_Thesis/Fig/
/content/drive/MyDrive/SAM_Thesis/paper2/
```

Before running a workflow, confirm that the required input files exist at the configured paths.

---

## Python environment

The repository is designed primarily for **Google Colab**.

Packages used across the workflows include:

```text
numpy
pandas
xarray
dask
netCDF4
h5netcdf
cftime
scipy
matplotlib
cartopy
gsw
statsmodels
openpyxl
```

Many of the later organized scripts automatically check for missing packages and install them with `pip`. Some earlier scripts use Colab/Jupyter `!pip` commands directly.

---

## Typical workflow

1. Open the required script in Google Colab.
2. Mount Google Drive.
3. Confirm the source datasets are available at the configured paths.
4. Allow installation of missing Python packages.
5. Run the script from top to bottom.
6. Review the generated NetCDF, CSV, Excel, PNG, PDF, and text outputs in the configured output directory.

Some computationally intensive workflows support:

- cached intermediate products;
- checkpoint CSV files;
- reuse switches;
- local staging under `/content`;
- low-RAM month-by-month processing.

These options are intended to make long Colab analyses more resilient and reproducible.

---

## Reproducibility notes

The code preserves important methodological settings directly in the workflows, including:

- the **2008–2025** main analysis period;
- **2001–2025** Argo coverage where applicable;
- the same 13 Antarctic shelf-sea sectors;
- Southern Hemisphere season definitions;
- explicit DJF handling;
- longitude normalization to `[-180°, 180°)`;
- common-mask requirements;
- bathymetric depth constraints;
- TEOS-10 calculations;
- fixed grid-cell reference profiles in corrected EN4 steric calculations;
- Sen-slope trend estimation;
- Hamed–Rao modified Mann–Kendall testing;
- Benjamini–Hochberg FDR correction;
- moving-block bootstrap confidence intervals;
- low-RAM and checkpointed processing where required.

Users reproducing the analyses should review the paths and settings in each script before execution.

---

## Repository structure

```text
.
├── Figure-10.py
├── Figure-11.py
├── Figure-12.py
├── Figure-13.py
├── Figure-14.py
├── Figure-15.py
├── Figure-16.py
├── Figure-17.py
├── Figure-18.py
├── Figure-19.py
├── Figure-A5.py
├── Figure-A6.py
├── Figure-A7.py
├── Figure-A8.py
├── Figure-A9.py
├── Figure-A10.py
├── .gitignore
├── LICENSE
└── README.md
```

---

## License

This repository is distributed under the license included in [`LICENSE`](./LICENSE).

---

## Repository scope

This repository provides the formatted and reproducible Python workflows for **Objective 2**, documenting Antarctic shelf-sea hydrography, steric sea-level variability, observed SLA, freshwater and salinity structure, upper-ocean stability, trends, EOF modes, atmospheric/oceanic forcing, AAO attribution, Argo validation, and depth-sensitivity analyses across the Southern Ocean.
