# avian_flu_sdm

## Data preparation

*Avian influenza (AI) data* 

Compiling the outcome data using publicly available data on avian influenza in wild birds from FAO and WOAH.

Main script for this is *create_flu_csv_no_bvbrc.R*


*Pseudoabsence generation*  

Select pseudoabsences at ~1:1 ratio from same 10km^2 grid, weighted sample by eBird surveillance records (10km^2 gridded counts of unique date/user/lat-long combinations recorded between 7/10/2005 and 30/6/2023, regardless of species sighted (or not sighted)) 

Pseudoabsences are also sampled with restriction they must be > 25km from positives 

 

*Thinning of AI data* 

Thin records and pseudoabsences on the 10km^2 grid independently by environmental filtering (sometimes caled “occfilt env”/”occfilter env”): divide environmental space into m x n strata (where m = number of environment layers and n = number of bins) and stratified sample from each 

We thin based on 9 purely environmental layers taking just one if a feature represented by multiple (distance to coast, distance to inland water, max elevation, diurnal temp range, precipitation, humidity, mean monthy temp, temp seasonality, ndvi) and use 6 bins 

 

*Test and train data sets* 

Training set A: all flu records before 2020/21 H5N8 outbreak (therefore includes 2017 H5N8 outbreak), plus pseudoabsences and thinned  
Test set A: 2020/21 H5N8 outbreak 
 
Training set B: all flu records post-2020/21 H5N8 outbreak (vast majority is 2021- H5N1 outbreak) 

 

*Variables for inclusion:* 

**Environmental**

Scripts for processing of the environmental variables are in the folder *variable_manipulation*

**Ecological** 

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

1. *assemble_covariates.R* loads in all of the covariate layers and generates a single multilayer raster, all_q_covs_10k.tif, containing all the covariates at 10km2 resolution. We also generate four other multilayer rasters, q*_covs_10k.tif which only contain the atemporal covariates (elevation etc) and the seasonally-specific ones for one season. These latter files are not used anywhere but are saved in case they are needed in alternative iterations of the model. 
2. *assemble_datasets.R* loads in the multilayer covariate raster generated in *assemble_covariates.R* along with all of the training/test coordinates. For each coordinate dataset (12 overall, A/B, Q1 to Q4, test/train for A, training only for B) it creates a CSV file containing the complete machine learning model output for that season. The CSV file contains the pos/neg status of each coordinate datapoint (column y), the country that coordinate lies in (column ri), and the values of the covariates at that coordinate. 
3. *fit_avian_flu_model_with_cv.R* loads in each dataset generated in assemble_covariates.R and performs the BART fitting for each of the datasets. The embarcadero package does not natively allow for custom BART parameters in the variable selection function bart.step. We would like to use this functionality since we can then improve model performance by optimising these parameters through crossvalidation and using them in the model. We thus define versions of a few embarcadero functions, rewritten to take BART parameters as arguments. We also define functions based on operations within the embarcadero::summary function which calculate test performance metrics (TSS, AUC, sensitivity, specificity) for a fitted model. For each dataset (A/B, Q1 to Q4) we first perform 5-fold crossvalidation using a for loop sweeping over a range of values for the BART parameters k, base, and power, and choose the values that give the highest mean AUC across the five folds. We then fit a BART model with these parameters and save it as cv_model_*_Q*.rds. We then use our customised version of the bart.step function to perform variable selection, leaving us with an optimally parsimonious model. We save this model as cv_model_with_vs_*_Q*.rds. 
4. *predictions_from_fitted_models.R* generates Europe-level geospatial avian flu presence probability projections with 95% confidence intervals (i.e. 2.5% and 97.5% layers) for each dataset/quarter combination based on the model fits obtained in *fit_avian_flu_model_with_cv.R*.  These are saved as predictions_*_Q*.rds. For dataset A, test metrics are calculated for each model and saved as metrics_A_Q*.rds. There is a global option B_TEST_DATA_AVAILABLE which, when set to TRUE, will cause metrics also to be calculated and saved for dataset B. This is currently set to FALSE as no such test sets exist as of 27/3/2024. An alternative would be to calculate training performance for the dataset B models – this is not currently done in the script. 

 

For step 3, there is an alternative: 

*fit_avian_flu_model_with_ri.rds* does the same tasks as *fit_avian_flu_model_with_cv.rds*, but adds a random intercept on country to the model. This requires some light “exploitation” of the embarcadero and dbarts packages – dbarts throws an error if custom values of the BART parameter k are passed to a random intercept model as function parameters in a nested function, but does not if this parameter is passed as a global variable. We can thus make the random intercept + custom parameters structure work by defining and redefining a global variable K_OPT, the optimal value for k obtained through crossvalidation. The other two BART parameters (power and base) do not have any such problems. Analogously to in *fit_avian_flu_model_with_cv.rds* , the fitted models from this script are saved as ri_model_*_Q*.rds and ri_model_with_vs_*_Q*.rds 
