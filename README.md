# Technical usage notes for Highly Pathogenic Avian Influenza (HPAI) SDM

The code in this repository was used to conduct the analysis outlined in the following bioRxiv preprint:
* Sarah Hayes, Joe Hilton, Joaquin Mould-Quevedo, Christl Donnelly, Matthew Baylis, Liam Brierley (2024) &quot;Ecology and environment predict spatially stratified risk of H5 highly pathogenic avian influenza clade 2.3.4.4b in wild birds across Europe.&quot; <i>bioRxiv</i> doi:10.1101/2024.07.17.603912

The covariate rasters calculated in the data preparation section of our analysis and used to train our BART models are available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15344899.svg)](https://doi.org/10.5281/zenodo.15344899).

Below we outline the steps involved in using the code to conduct the analysis.

## Data preparation

### Avian influenza (AI) data 

Avian influenza detection data is compiled into a .csv file using publicly available data sources on avian influenza in wild birds from FAO (EMPRES-i data product) and WOAH (WAHIS data product) in *avian_flu_scripts/create_flu_csv_no_bvbrc.R*. Only high pathogenicity avian influenza (HPAI) subtypes are included. 

### Construction of training and testing data sets 

*resample_training_data.R* defines scope of training and testing periods. Cleaned avian influenza data is read in before filtering to H5 subtypes only (assumed to be clade 2.3.4.4b based on date/location), before assigning training/test divisions:

Training set A: 10/8/2016 to 9/8/2020 (primarily 2016/2017 H5N8 outbreak and sporadic 2018 H5N6 infections)

Test set A: 10/8/2020 to 9/8/2021 (discrete 2020/21 H5N8 outbreak)
 
Training set B: 10/8/2021 to 28/2/2023 (start of 2021 H5N1 outbreak)

Test set B: 28/2/2023 to 29/2/2024 (most recent complete calendar year of H5N1 outbreak at time of data extraction)

Raw records are plotted in time and space as `data_fig.png`, before stratifying by bird behavioural season and assigning cells with known influenza reports on a 10km^2 grid as positives (regardless of number of cases/reports).

### Pseudoabsence generation

Both pseudoabsence generation and data thinning are handled by *resample_training_data.R*

Pseudoabsences are selected for each season of each dataset at an approximately 1:1 ratio relative to the number of positives from same 10km^2 grid, in a sample weighted by log-eBird surveillance records (see *Weighting layer* below) and with further restriction that they must be greater than 25km from any positive detection.

### Thinning of avian influenza incidence data 

Positives and pseudoabsences on the 10km^2 grid are then thinned independently by environmental filtering (sometimes called “occfilt env”/”occfilter env”): this divides the environmental space into m x n strata (where m = number of environment layers and n = number of bins) and takes a stratified sample from each (see [Varela et al. 2014](https://doi.org/10.1111/j.1600-0587.2013.00441.x)).

We thin based on 10 purely environmental layers each divided into 6 strata (distance to coast, distance to inland water, max elevation, mean zero-degree isotherm, mean precipitation, mean relative humidity, mean monthly temperature, diurnal temperature range, seasonal temperature variation, mean normalised difference vegetation index (NDVI)). Only one data point is retained per combination of strata in this 10-dimensional environmental space.

Thinned positives and thinned pseudoabsences are finally combined to create single training and test set files for each season of each dataset, which are saved in folder `training_sets`.



### Variable manipulation

#### Weighting layer

*prepping_ebird_obs_layer.R* reads in the [eBird Basic Dataset (EBD)](https://science.ebird.org/en/use-ebird-data/download-ebird-data-products) and filters to reported citizen science bird sightings in Europe between 10/8/2016 and 29/2/2024 based on spatial extent. These are then mapped to a 10km^2 grid and summed as counts of unique date/user/latitude-longitude combinations (regardless of species sighted or not sighted) to represent bird surveillance accessibility. The outputted raster is used in pseudoabsence selection, but is **not** used as input to BART models.

#### Environmental

Scripts for processing of the environmental covariates are in the folder *variable_manipulation*.

#### Ecological 

A preliminary analysis to identify species-level factors associated with avian flu host status is performed in *host_species_analysis_all_ET.R* (host_species_analysis.R looks at only those species present in Europe according to eBird and is thus less comprehensive than *host_species_analysis_all_ET.R*, which covers a maximal set of species found in eBird, EltonTraits, and the IUCN Red List). Based on the findings around variable importance obtained in this script, we choose specific species-level factors to include in our geospatial model, outlined below. The rasters of ecological covariates are constructed in *make_eco_phylo_rasters.R.* 

- Species richness: multiply eBird spatial percentage prevalence by Callaghan et al global population estimates, and count how many species have estimated prevalence >1 in each cell 
- Predictor variables for taxon x: estimated absolute prevalence of Ardeidae, Arenaria and Calidris, Laridae, Anserinae, Aythyini, and Anatinae; we take a conservative view of what genera are in Anatinae 
- Phylogenetic distance to host species: estimated absolute prevalence across species which are within 50 phylogenetic units (i.e. MYA since divergence) of a known host species 
- Behavioural risk factors: absolute prevalence of species that IUCN labels as 
  - Congregative 
  - Migratory 
- Dietary risk factors: eBird spatial % prevalence multiplied by Callaghan et al. global population estimate multiplied by % of diet/foraging behaviour contributed by risky behaviour/foodsource 
  - Foraging around water surface 
  - Foraging below water surface 
  - Plants in diet 
  - Scavenging in diet 
  - Endotherms in diet 


## BART model 

Model construction, fitting, and analysis is done across these scripts: 

1\. *assemble_covariates.R* loads in all of the covariate layers and generates a single multilayer raster, all_q_covs_10k.tif, containing all the covariates at 10km2 resolution. We also generate four other multilayer rasters, q*_covs_10k.tif which contain the non-seasonally-variant covariates (elevation, etc.) as well as seasonally-variant covariates for the given season. These latter files are not used anywhere but are saved in case they are needed in alternative iterations of the model.

2\. *assemble_datasets.R* loads in the multilayer covariate raster generated in *assemble_covariates.R* along with all of the training/test coordinates. For each coordinate dataset (14 overall: one training set and one test set for each of the non-breeding, pre-breeding migration, and post-breeding migration seasons for period A, and one training set and one test set for each of the four behavioural seasons for period B) it creates a CSV file containing the complete machine learning model input for that season. The .csv file contains the positive/pseudoabsence status of each coordinate datapoint (column y), the country that coordinate lies in (column ri), and the values of the covariates at that coordinate. 

3a. *fit_avian_flu_model_with_cv.R* loads in each dataset generated in assemble_covariates.R and performs the BART fitting for each of the datasets. The embarcadero package does not natively allow for custom BART parameters in the variable selection function bart.step. We would like to use this functionality since we can then improve model performance by optimising these parameters through crossvalidation and using them in the model. We thus define versions of a few embarcadero functions, rewritten to take BART parameters as arguments. We also define functions based on operations within the embarcadero::summary function which calculate test performance metrics (true skills statistic (TSS), area under the curve (AUC), sensitivity, specificity) for a fitted model. For each dataset (i.e. each behavioural season for periods A and B) we first perform 5-fold cross-validation using a for loop sweeping over a range of values for the BART parameters k, base, and power, and choose the values that give the highest mean AUC across the five folds. We then fit a BART model with these parameters and save it as `cv_model_\*\_Q\*.rds`. We then use our customised version of the bart.step function to perform variable selection, leaving us with an optimally parsimonious model. We save this model as `cv_model_with_vs_\*_Q\*.rds`. Here "Q" stands for "Quarter", i.e. behavioural season, with Q1 being shorthand for the nonbreeding season, Q2 for the pre-breeding migration season, Q3 for the breeding season, and Q4 for the post-breeding migration season.

3b. *fit_avian_flu_model_with_ri.R* does the same tasks as *fit_avian_flu_model_with_cv.R*, but adds a random intercept on country to the model. This requires some light “exploitation” of the embarcadero and dbarts packages – dbarts throws an error if custom values of the BART parameter k are passed to a random intercept model as function parameters in a nested function, but does not if this parameter is passed as a global variable. We can thus make the random intercept + custom parameters structure work by defining and redefining a global variable K_OPT, the optimal value for k obtained through crossvalidation. The other two BART parameters (power and base) do not have any such problems. Analogously to in *fit_avian_flu_model_with_cv.R* , the fitted models from this script are saved as `ri_model_\*\_Q\*.rds` and `ri_model_with_vs_\*_Q\*.rds`.

4\. *predictions_from_fitted_models.R* generates Europe-level geospatial avian flu presence probability projections with 95% credible intervals (i.e. 2.5% and 97.5% percentile layers) for each dataset/season combination based on the model fits obtained in *fit_avian_flu_model_with_cv.R*. These are saved as `predictions_\*\_Q\*.rds` and plots produced via *plot_projections.R*. Test metrics are calculated for each model and saved as `metrics_\*_Q\*.rds`. For the random intercept models, this script just calculates test performance metrics and does not attempt to generate predictions.

5\. We provide bash shell scripts with names in the format *sbatch_fit_\*.sh* specifying the different combinations of random intercepts and cross-seasonal covariates which can be used to run the *fit_avian_flu_model_\*.R* scripts on Unix-based high-performance computing environments. For the models without random intercepts, these shell scripts run each chain of the MCMC in serial; running *fit_avian_flu_model_with_cv.R* will run the chains in parallel.

 
### Model interpretation

##### Variable importance

Both variable importance and partial dependence are handled by *varimp_pd_plot.R*.

Models fits for each dataset/season combination are loaded and variable importance calculated as relative frequency of retention of each individual covariate across component trees within posterior draws of the BART model, rescaled by the maximum frequency count. Means and standard deviations are then calculated in relative variable importance for each covariate and plotted as `variable_importance_quarterly_\*.png`

##### Partial dependence

Partial dependence is then calculated for model fits for each dataset/season combination, for either the top six covariates per model by average variable importance (if `ALL_VARS_FOR_PD = FALSE`) or all covariates per model (if `ALL_VARS_FOR_PD = TRUE`). Partial dependence is calculated for continuous covariates by generating model predicted probabilities for shifting values of the focal covariate whilst averaging over all observed values for other covariates (marginal predictions), and for binary covariates by generating predictions for `x = 0` and `x = 1` in the same way. Medians, 2.5th, and 97.5th percentiles in predicted probabilities are then calculated for each covariate and plotted as either `partial_dependence_top_mean_*.png` (if `ALL_VARS_FOR_PD = FALSE`) or `partial_dependence_\*_Q\*.png`(if `ALL_VARS_FOR_PD = TRUE`).

Plots of variable importance and partial dependence combining all dataset/season combinations are also generated as `variable_importance_quarterly_combi.png` and `partial_dependence_top_mean_combi.png`


