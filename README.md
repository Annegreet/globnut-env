# globnut-env

This repo contains the code for the manuscript "Nutrient stoichiometry mediates nitrogen deposition effects on plant diversity" (in prep).
Note that figures in the manuscript are post-processed in Adobe Illustrator and appear different in R than presented in the manuscript. 

## Environmental data

`/env_data` contains the scripts to extract and process environmental data for analysis. The repo `/env_data/outputs` contains the data produced by these scripts.

-   `Extract_N_deposition.R` - to calculate the N deposition variable and fig. 1b

-   `Extract_DEM.R` - to extract elevation from DEM
  
-   `Extract_species_pool_Cai - to extract species pool values as described in the manuscript

## Preparing data

-   `01_Gridding_globnut.R` - to create the grid for random effects and sub-sampling

-   `01_Species_diversity_indices.R` - to calculate species diversity indices

-   `02_Plot_selection.R` - compiles data set for analysis and performs sub-sampling

## Main analysis

Note: figures in manuscripts are post-processed in Adobe Illustrator, so figures produced by the script will differ in appearance.

-   `03_Geography_limitation.R`  - produces input for fig 1a and 1c

-   `03_SEM.R` - produces SEM, displayed in fig 2 and supplementary information figures and tables

-   `03_Alpha_diversity_analysis.Rmd` - linear mixed modeling, fig 3, and supplementary information figures and tables

-   `03_Beta_diversity_analysis.Rmd`  - iNEXT beta diversity estimation, linear mixed modeling,  and produces fig. 4 and supplementary information tables

Data to run these scripts is available at: <https://osf.io/pgxrt/>

## R packages

Users can use the `renv` R package to make sure the same R packages and versions are used. Run `renv::restore()` to install missing packages.
