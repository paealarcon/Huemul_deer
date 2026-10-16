# Huemul_deer — Reproducibility package

Scripts, input data, and processed outputs required to reproduce the full
analysis reported in Alarcón et al., *Movement behaviour and habitat selection
to guide conservation of an endangered Neotropical deer* (JMAMM-2026-037).

## Folder layout

    Huemul_deer/
    ├── Huemul_deer.Rproj
    ├── data/
    │   ├── collarData/      raw & cleaned GPS collar data
    │   └── layers/          spatial covariates (habitat shapefiles,
    │                        MCP polygons, river/stream zonification,
    │                        static raster stack)
    ├── lsm/                 Latent-State Model (Fig. 1, SD4, SD5)
    ├── akde/                AKDE home-range estimation (SD6)
    ├── hrsf/                Landscape-scale habitat selection
    │                        (Fig. 2, Fig. 3, SD7, SD8)
    └── issf/                Patch-scale habitat selection
                             (Fig. 4, SD9, SD10)

## Running the pipeline

> **IMPORTANT.** Always start by double-clicking **`Huemul_deer.Rproj`** at
> the root of this folder. RStudio will open with the working directory
> already set to the project root, which is what every script expects. Do
> NOT open individual `.R` files from Windows Explorer or your file
> manager — the working directory will be wrong and relative paths will fail.

Once the project is open, run the scripts in this order:

### 1. Latent-State Model
    lsm/scripts/01a_lsmEstimation.R
    lsm/scripts/01b_lsmClassification.R

### 2. Home ranges
    akde/scripts/02_homeRangeEstimation.R

### 3. Landscape-scale habitat selection
    hrsf/scripts/00_prepare_used_available.R
    hrsf/scripts/03a_hrsfEstimation_forestIntercept.R
    hrsf/scripts/03b_hrsfEstimation_steppeIntercept.R
    hrsf/scripts/03c_BoyceIndex.R
    hrsf/scripts/03d_suitabilityMaps.R
    hrsf/scripts/03h_export_availability.R

### 4. Patch-scale habitat selection
    issf/scripts/04_prep_huemul_steps.R
    issf/scripts/04a_hissfEstimation.R
    issf/scripts/04b_hissfSimulation.R

### Figures and tables
    lsm/scripts/01c_Figure1.R                      -> Fig. 1
    lsm/scripts/01d_FigureSD5_NSD_timeseries.R     -> SD5
    hrsf/scripts/03e_Figure2_main.R                -> Fig. 2
    hrsf/scripts/03e_FigureSD7_supp.R              -> SD7 figure
    hrsf/scripts/03f_Figure3_main.R                -> Fig. 3
    hrsf/scripts/03f_Figure3_supp.R                -> SD8 figure
    hrsf/scripts/03g_Table2_hRSF.R                 -> SD7 table
    issf/scripts/04c_Figure4.R                     -> Fig. 4
    issf/scripts/04d_Table3_iSSF.R                 -> SD9 table

## Pre-computed outputs

Each module's `outputs/` folder already contains the `.rda` / `.rds`
fits and `.csv` tables, so figures and tables can be reproduced without
re-running the Bayesian models.

## Software

- R 4.4.1
- R-INLA (stable branch, 2024)
- JAGS 4.3.1 + jagsUI
- Packages: ctmm, terra, sf, ggplot2, patchwork, cowplot, elevatr

