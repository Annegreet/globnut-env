# globnut-env

This repo contains the code for the article "Nutrient stoichiometry mediates nitrogen deposition effects on plant diversity" (in prep).

## Environmental data

`/env_data` contains the scripts to extract and process environmental data for analysis. The repo `/env_data/outputs` contains the data produced by these scripts.

-   `Extract_geology.R` - for computing the rock acidity variable

-   `Extract_glacial_extent.R`  - for computing the soil age variable

-   `Extract_N_deposition.R` - to calculate the N deposition variable and fig. 1b

## Preparing data

-   `01_Gridding_globnut.R` - to create the grid for random effects and sub-sampling

-    `01_Species_diversity_indices.R` - to calculate species diversity indices

-   `02_Plot_selection.R` - compiles data set for analysis and performs sub-sampling

## Main analysis

Note: figures in manuscripts are post-processed in Adobe Illustrator, so figures produced by the script will differ in appearance.

-   `03_Geography_limitation.R`  - produces input for fig 1a and 1c

-   `03_SEM_GlobNut.R` - produces SEM, displayed in fig 2 and supplementary information fig. 1.3

-   `03_Linear_mixed_models.Rmd` - linear mixed modelling, fig 3 and supplementary information fig. S2.1 and tables 2.1-2.4

-   `03_Zi_Beta_regression.Rmd`  - calculates distance metrics, runs zero-inflated beta regression, produces fig. 4

Data to run these scripts is available at: <https://osf.io/pgxrt/>

## R packages

Users can use the `renv` R package to make sure the same R packages and versions are used. Run `renv::restore()` to install missing packages.
