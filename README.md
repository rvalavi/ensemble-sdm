# Ensemble Species Distribution Modelling under Climate Change

An R workflow for fitting and evaluating ensemble species distribution models (SDMs) and mapping habitat suitability under historical and projected climates.

## Species of interest

- Leopard cat (*Prionailurus bengalensis*)
- Cycad (*Zamia prasina*, also recorded as *Zamia polymorpha*)

Both examples use bundled GBIF records with `1` for a presence and `0` for a background location. The modelling scripts remove duplicate records within each raster cell and class, then extract environmental predictors at the remaining locations.

## Climate and environmental covariates

The workflow fits a principal component analysis (PCA) to 19 CHELSA v2.1 bioclimatic variables for 1981–2010. It uses the same PCA to transform GFDL-ESM4 projections under SSP5-8.5 for 2041–2070 and 2071–2100, so the climate predictors are comparable across periods. Each species uses its supplied topographic position index raster in `data/Topo/` for all three periods.

## How to run the scripts

Open [ensemble-sdm.Rproj](ensemble-sdm.Rproj) in RStudio or start R from the repository root. The climate download uses `curl`, and the default MaxEnt fit requires Java. Install Java before installing `rJava`.

```r
install.packages(c(
    "magrittr", "terra", "geodata", "tidyverse", "sf", "blockCV", "janitor",
    "dismo", "gbm", "mgcv", "glmnet", "precrec", "ecospat", "randomForest",
    "rJava", "predicts", "remotes"
))
remotes::install_github("rvalavi/myspatial")
```

Check that `dismo::maxent()` returns `TRUE` before fitting. [dismo manual](https://cran.r-project.org/web/packages/dismo/refman/dismo.html)

Run the scripts in order in the same R session.

1. Run [`R/climate_data.R`](R/climate_data.R) to download the climate rasters and crop them to the study regions.
2. Create `data/CHELSA_data/PCA/`, then run [`R/climate_pca.R`](R/climate_pca.R) to save the PCA fits and climate predictor rasters.
3. Use the bundled species data. To refresh the GBIF records, install `spocc` and `rgbif`, then run [`R/species_data.R`](R/species_data.R).
4. Create `outputs/plant/` or `outputs/cat/`, then run [`R/modelling_plant.R`](R/modelling_plant.R) or [`R/modelling_cat.R`](R/modelling_cat.R) to fit and project the corresponding ensemble. The cat example may take longer to run.

## Usage

Use `model_data` and `covar_rast` from either modelling script, or supply a training data frame and a `SpatRaster` of predictors. The training data should contain an `occ` column coded as 0 or 1 and one column per predictor. Match the raster layer names to the predictor columns.

```r
source("R/models.R")

model <- ensemble(x = model_data, y = "occ")
suitability <- terra::predict(covar_rast, model, na.rm = TRUE)
terra::plot(suitability)
```

## Evaluation and outputs

Both examples use five-fold spatial block cross-validation and report ROC AUC, area under the precision–recall curve, and the Boyce index. They then refit the ensemble using all available records. The ensemble scales each model's predictions to 0–1 using its training prediction range, then averages the scaled scores. It reuses those ranges for all three periods.

| GeoTIFF | Period |
| --- | --- |
| `pred_current.tif` | 1981–2010 baseline |
| `pred_f1.tif` | 2041–2070 |
| `pred_f2.tif` | 2071–2100 |

The scripts also plot response curves. The plant example produces a multivariate environmental similarity surface (MESS) plot to flag late-century conditions outside the training ranges.

## Choosing ensemble models

The default ensemble combines a lasso-regularised GLM, a GAM, boosted regression trees (`GBM`), random forest (`RF`), and MaxEnt. Use `models` to select a subset. For presence–absence data, omit `Maxent` from the model list.

Random Forest balances classes by default. Set `balance_classes = FALSE` to disable it, as shown in this presence–absence example.

```r
model <- ensemble(
    x = model_data,
    y = "occ",
    models = c("GLM", "GAM", "GBM", "RF"),
    balance_classes = FALSE
)
```

## References

Valavi, R., Guillera-Arroita, G., Lahoz-Monfort, J. J., & Elith, J. (2022). Predictive performance of presence-only species distribution models: a benchmark study with reproducible code. *Ecological Monographs*, 92(1), e01486. [doi:10.1002/ecm.1486](https://doi.org/10.1002/ecm.1486)
