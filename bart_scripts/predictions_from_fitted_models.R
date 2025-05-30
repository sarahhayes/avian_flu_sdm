# In this script we load the fitted models, calculate performance metrics, and
# perform Europe-wide predictions.

SKIP_AQ3 <- TRUE # Set to TRUE to skip small period A quarter 3 dataset

# Optional command line arguments, must be passed as strings:
args <- commandArgs(trailingOnly = T)
if (length(args)<4){
  INCLUDE_CROSSTERMS <- "with-crossterms-multichain" # Set to "no-crossterms" to do model without crossterms or "with-crossterms" to do model with crossterms
}else{
  INCLUDE_CROSSTERMS <- args[4]
}
if (length(args)<3){
  CV_OR_RI <- "cv" # Set to "cv" to do crossvalidated model or "ri" to do crossvalidation + random intercept model
}else{
  CV_OR_RI <- args[3]
}
if (length(args)<2){
  # Set path to folder containing data, and where output will be stored
  PATH_TO_COVS <- "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/"
}else{
  PATH_TO_COVS <- args[2]
}
if (length(args)<1){
  # Set path to folder containing data, and where output will be stored
  PATH_TO_OUTPUTS <- ""
}else{
  PATH_TO_OUTPUTS <- args[1]
}

# Decide whether to use models with ecological season boundaries
# If you want to distinguish between models based on bird behavioural seasons vs
# calendar quarters you can do it by storing the fitted models in the directory
# specified here.
ECO_SEASONS <- TRUE
if (ECO_SEASONS){
  PATH_TO_MODELS <- paste(PATH_TO_OUTPUTS,
                          "",
                          sep = "")
}else{
  PATH_TO_MODELS <- PATH_TO_OUTPUTS
}

SAVE_OUTPUTS <- TRUE

B_TEST_DATA_AVAILABLE <- TRUE # Set to TRUE if we have test data for dataset B, otherwise this should be FALSE

library(dplyr)
library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

#### Functions for calculating performance metrics ####

get_threshold <- function(object){
  if(class(object)=='rbart') {
    fitobj <- object$fit[[1]]
  }
  if(class(object)=='bart') {
    fitobj <- object$fit
  }
  
  true.vector <- fitobj$data@y 
  
  pred <- prediction(colMeans(pnorm(object$yhat.train)), true.vector)
  
  perf.tss <- performance(pred,"sens","spec")
  tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
  tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
  thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))])
  return(thresh)
}

get_sens_and_spec <- function(sdm, xtest, ytest, ri, cutoff){
  pred <- xtest[which(complete.cases(xtest)), ] %>%
    stats::predict(object=sdm, type = "bart", group.by=ri) %>%
    pnorm() %>%
    colMeans() %>%
    prediction(labels = ytest[which(complete.cases(xtest))])
  perf <-  performance(pred, measure = "sens", x.measure = "spec")
  tss_list <- (perf@x.values[[1]] + perf@y.values[[1]] - 1)
  tss_df <- data.frame(alpha=perf@alpha.values[[1]],tss=tss_list)
  # cutoff <- min(tss_df$alpha[which(tss_df$tss==max(tss_df$tss))])
  sens <- perf@x.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  spec <- perf@y.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  tss <- tss_df[which.min(abs(tss_df$alpha-cutoff)),'tss']
  auc <- performance(pred,"auc")@y.values[[1]]
  return(pairlist("pred"=pred,
                  "perf"=perf,
                  "sens"=sens,
                  "spec"=spec,
                  "tss"=tss,
                  "auc"=auc))
}

get_sens_and_spec_for_single_model <- function(pred){
  perf <-  performance(pred, measure = "sens", x.measure = "spec")
  tss_list <- (perf@x.values[[1]] + perf@y.values[[1]] - 1)
  tss_df <- data.frame(alpha=perf@alpha.values[[1]],tss=tss_list)
  # cutoff <- min(tss_df$alpha[which(tss_df$tss==max(tss_df$tss))])
  sens <- perf@x.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  spec <- perf@y.values[[1]][which.min(abs(perf@alpha.values[[1]]-cutoff))]
  tss <- tss_df[which.min(abs(tss_df$alpha-cutoff)),'tss']
  auc <- performance(pred,"auc")@y.values[[1]]
  return(pairlist("sens"=sens,
                  "spec"=spec,
                  "tss"=tss,
                  "auc"=auc))
}

get_sens_and_spec_ci <- function(sdm, xtest, ytest, ri, cutoff){
  predmat <- xtest[which(complete.cases(xtest)), ] %>%
    stats::predict(object=sdm, type = "bart", group.by=ri) %>%
    pnorm()
  pred_by_model <- sapply(1:nrow(predmat),
                          FUN=function(i){predmat[i,] %>%
                              prediction(labels = ytest[which(complete.cases(xtest))])})
  metrics_by_model <- lapply(1:length(pred_by_model),
                             FUN=function(i){
                               get_sens_and_spec_for_single_model(pred_by_model[[i]])})
  sens_by_model <- sapply(1:length(metrics_by_model),
                         FUN=function(i){metrics_by_model[[i]]$sens})
  spec_by_model <- sapply(1:length(metrics_by_model),
                         FUN=function(i){metrics_by_model[[i]]$spec})
  tss_by_model <- sapply(1:length(metrics_by_model),
                         FUN=function(i){metrics_by_model[[i]]$tss})
  auc_by_model <- sapply(1:length(metrics_by_model),
                         FUN=function(i){metrics_by_model[[i]]$auc})
  return(pairlist("sens"=quantile(sens_by_model, c(.025, .975), names=FALSE),
                  "spec"=quantile(spec_by_model, c(.025, .975), names=FALSE),
                  "tss"=quantile(tss_by_model, c(.025, .975), names=FALSE),
                  "auc"=quantile(auc_by_model, c(.025, .975), names=FALSE)))
}

#### Bring in year-round covariates ####

covstack <- raster::stack(paste(PATH_TO_COVS, "quarterly_covariates_eco_seasons/all_q_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

#### Dataset A Q1 ####

test_data <- read.csv("training_sets/test_data_eco_seasons_A_Q1.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q1.rds", sep = ""))
sdm_A_Q1 <- sdm
cutoff <- get_threshold(sdm)
if (CV_OR_RI=="ri"){
  metrics_A_Q1 <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
  metric_cis_A_Q1 <- get_sens_and_spec_ci(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metrics_A_Q1 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
  metric_cis_A_Q1 <- get_sens_and_spec_ci(sdm, xtest, ytest, NULL, cutoff)
}

if (SAVE_OUTPUTS){
  save(metrics_A_Q1,
       metric_cis_A_Q1,
       file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_A_Q1.rds", sep = ""))
}

if (CV_OR_RI=="cv"){
  # Generate risk map with percentiles
  pred_layers_A_Q1 <- predict(object = sdm,
                        x.layers = covstack,
                        quantiles = c(0.025, 0.975),
                        splitby = 20
  )
  names(pred_layers_A_Q1) <- c("Mean",
                               "Lower 95 percent confidence bound",
                               "Upper 95 percent confidence bound")
  if (SAVE_OUTPUTS){
    save(pred_layers_A_Q1,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_A_Q1.rds", sep = ""))
  }
}

#### Dataset A Q2 ####

test_data <- read.csv("training_sets/test_data_eco_seasons_A_Q2.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q2.rds", sep = ""))
sdm_A_Q2 <- sdm
cutoff <- get_threshold(sdm)
if (CV_OR_RI=="ri"){
  metrics_A_Q2 <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
  metric_cis_A_Q2 <- get_sens_and_spec_ci(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metrics_A_Q2 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
  metric_cis_A_Q2 <- get_sens_and_spec_ci(sdm, xtest, ytest, NULL, cutoff)
}

if (SAVE_OUTPUTS){
  save(metrics_A_Q2,
       metric_cis_A_Q2,
       file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_A_Q2.rds", sep = ""))
}

if (CV_OR_RI=="cv"){
  # Generate risk map with percentiles
  pred_layers_A_Q2 <- predict(object = sdm,
                              x.layers = covstack,
                              quantiles = c(0.025, 0.975),
                              splitby = 20
  )
  names(pred_layers_A_Q2) <- c("Mean",
                               "Lower 95 percent confidence bound",
                               "Upper 95 percent confidence bound")
  if (SAVE_OUTPUTS){
    save(pred_layers_A_Q2,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_A_Q2.rds", sep = ""))
  }
}


#### Dataset A Q3 ####

if (!SKIP_AQ3){
  test_data <- read.csv("training_sets/test_data_eco_seasons_A_Q3.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  
  load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q3.rds", sep = ""))
  sdm_A_Q3 <- sdm
  cutoff <- get_threshold(sdm)
  if (CV_OR_RI=="ri"){
    metrics_A_Q3 <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
    metric_cis_A_Q3 <- get_sens_and_spec_ci(sdm, xtest, ytest, countrytest, cutoff)
  }else{
    metrics_A_Q3 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
    metric_cis_A_Q3 <- get_sens_and_spec_ci(sdm, xtest, ytest, NULL, cutoff)
  }
  
  if (SAVE_OUTPUTS){
    save(metrics_A_Q3,
         metric_cis_A_Q3,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_A_Q3.rds", sep = ""))
  }
  
  if (CV_OR_RI=="cv"){
    # Generate risk map with percentiles
    pred_layers_A_Q3 <- predict(object = sdm,
                                x.layers = covstack,
                                quantiles = c(0.025, 0.975),
                                splitby = 20
    )
    names(pred_layers_A_Q3) <- c("Mean",
                                 "Lower 95 percent confidence bound",
                                 "Upper 95 percent confidence bound")
    if (SAVE_OUTPUTS){
      save(pred_layers_A_Q3,
           file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_A_Q3.rds", sep = ""))
    }
  }

}
#### Dataset A Q4 ####

test_data <- read.csv("training_sets/test_data_eco_seasons_A_Q4.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q4.rds", sep = ""))
sdm_A_Q4 <- sdm
cutoff <- get_threshold(sdm)
if (CV_OR_RI=="ri"){
  metrics_A_Q4 <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
  metric_cis_A_Q4 <- get_sens_and_spec_ci(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metrics_A_Q4 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
  metric_cis_A_Q4 <- get_sens_and_spec_ci(sdm, xtest, ytest, NULL, cutoff)
}
{
cat("Dataset A Q1, test metrics:\n sens =",
    signif(metrics_A_Q1$sens, digits=2),
    " (",
    signif(metric_cis_A_Q1$sens, digits=2),
    "),\n spec =",
    signif(metrics_A_Q1$spec, digits=2),
    " (",
    signif(metric_cis_A_Q1$spec, digits=2),
    "),\n AUC =",
    signif(metrics_A_Q1$auc, digits=2),
    " (",
    signif(metric_cis_A_Q1$auc, digits=2),
    ")\n")
cat("Dataset A Q2, test metrics:\n sens =",
    signif(metrics_A_Q2$sens, digits=2),
    " (",
    signif(metric_cis_A_Q2$sens, digits=2),
    "),\n spec =",
    signif(metrics_A_Q2$spec, digits=2),
    " (",
    signif(metric_cis_A_Q2$spec, digits=2),
    "),\n AUC =",
    signif(metrics_A_Q2$auc, digits=2),
    " (",
    signif(metric_cis_A_Q2$auc, digits=2),
    ")\n")
if (!SKIP_AQ3){
  cat("Dataset A Q3, test metrics:\n sens =",
      signif(metrics_A_Q3$sens, digits=2),
      " (",
      signif(metric_cis_A_Q3$sens, digits=2),
      "),\n spec =",
      signif(metrics_A_Q3$spec, digits=2),
      " (",
      signif(metric_cis_A_Q3$spec, digits=2),
      "),\n AUC =",
      signif(metrics_A_Q3$auc, digits=2),
      " (",
      signif(metric_cis_A_Q3$auc, digits=2),
      ")\n")
}
cat("Dataset A Q4, test metrics:\n sens =",
    signif(metrics_A_Q4$sens, digits=2),
    " (",
    signif(metric_cis_A_Q4$sens, digits=2),
    "),\n spec =",
    signif(metrics_A_Q4$spec, digits=2),
    " (",
    signif(metric_cis_A_Q4$spec, digits=2),
    "),\n AUC =",
    signif(metrics_A_Q4$auc, digits=2),
    " (",
    signif(metric_cis_A_Q4$auc, digits=2),
    ")\n")
}

if (SAVE_OUTPUTS){
  save(metrics_A_Q4,
       metric_cis_A_Q4,
       file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_A_Q4.rds", sep = ""))
}

if (CV_OR_RI=="cv"){
  # Generate risk map with percentiles
  pred_layers_A_Q4 <- predict(object = sdm,
                              x.layers = covstack,
                              quantiles = c(0.025, 0.975),
                              splitby = 20
  )
  names(pred_layers_A_Q4) <- c("Mean",
                               "Lower 95 percent confidence bound",
                               "Upper 95 percent confidence bound")
  if (SAVE_OUTPUTS){
    save(pred_layers_A_Q4,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_A_Q4.rds", sep = ""))
  }
}


#### Dataset B Q1 ####

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q1.rds", sep = ""))
sdm_B_Q1 <- sdm

if (B_TEST_DATA_AVAILABLE){
  test_data <- read.csv("training_sets/test_data_eco_seasons_B_Q1.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  cutoff <- get_threshold(sdm)
  if (CV_OR_RI=="ri"){
    metrics_B_Q1 <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
    metric_cis_B_Q1 <- get_sens_and_spec_ci(sdm, xtest, ytest, countrytest, cutoff)
  }else{
    metrics_B_Q1 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
    metric_cis_B_Q1 <- get_sens_and_spec_ci(sdm, xtest, ytest, NULL, cutoff)
  }
  
  if (SAVE_OUTPUTS){
    save(metrics_B_Q1,
         metric_cis_B_Q1,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_B_Q1.rds", sep = ""))
  }
}

if (CV_OR_RI=="cv"){
  # Generate risk map with percentiles
  pred_layers_B_Q1 <- predict(object = sdm,
                              x.layers = covstack,
                              quantiles = c(0.025, 0.975),
                              splitby = 20
  )
  names(pred_layers_B_Q1) <- c("Mean",
                               "Lower 95 percent confidence bound",
                               "Upper 95 percent confidence bound")
  if (SAVE_OUTPUTS){
    save(pred_layers_B_Q1,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_B_Q1.rds", sep = ""))
  }
}




#### Dataset B Q2 ####

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q2.rds", sep = ""))
sdm_B_Q2 <- sdm

if (B_TEST_DATA_AVAILABLE){
  test_data <- read.csv("training_sets/test_data_eco_seasons_B_Q2.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  cutoff <- get_threshold(sdm)
  if (CV_OR_RI=="ri"){
    metrics_B_Q2 <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
    metric_cis_B_Q2 <- get_sens_and_spec_ci(sdm, xtest, ytest, countrytest, cutoff)
  }else{
    metrics_B_Q2 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
    metric_cis_B_Q2 <- get_sens_and_spec_ci(sdm, xtest, ytest, NULL, cutoff)
  }
  
  if (SAVE_OUTPUTS){
    save(metrics_B_Q2,
         metric_cis_B_Q2,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_B_Q2.rds", sep = ""))
  }
}

if (CV_OR_RI=="cv"){
  # Generate risk map with percentiles
  pred_layers_B_Q2 <- predict(object = sdm,
                              x.layers = covstack,
                              quantiles = c(0.025, 0.975),
                              splitby = 20
  )
  names(pred_layers_B_Q2) <- c("Mean",
                               "Lower 95 percent confidence bound",
                               "Upper 95 percent confidence bound")
  if (SAVE_OUTPUTS){
    save(pred_layers_B_Q2,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_B_Q2.rds", sep = ""))
  }
}





#### Dataset B Q3 ####

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q3.rds", sep = ""))
sdm_B_Q3 <- sdm

if (B_TEST_DATA_AVAILABLE){
  test_data <- read.csv("training_sets/test_data_eco_seasons_B_Q3.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  cutoff <- get_threshold(sdm)
  if (CV_OR_RI=="ri"){
    metrics_B_Q3 <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
    metric_cis_B_Q3 <- get_sens_and_spec_ci(sdm, xtest, ytest, countrytest, cutoff)
  }else{
    metrics_B_Q3 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
    metric_cis_B_Q3 <- get_sens_and_spec_ci(sdm, xtest, ytest, NULL, cutoff)
  }
  
  if (SAVE_OUTPUTS){
    save(metrics_B_Q3,
         metric_cis_B_Q3,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_B_Q3.rds", sep = ""))
  }
}

if (CV_OR_RI=="cv"){
  # Generate risk map with percentiles
  pred_layers_B_Q3 <- predict(object = sdm,
                              x.layers = covstack,
                              quantiles = c(0.025, 0.975),
                              splitby = 20
  )
  names(pred_layers_B_Q3) <- c("Mean",
                               "Lower 95 percent confidence bound",
                               "Upper 95 percent confidence bound")
  if (SAVE_OUTPUTS){
    save(pred_layers_B_Q3,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_B_Q3.rds", sep = ""))
  }
}





#### Dataset B Q4 ####

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q4.rds", sep = ""))
sdm_B_Q4 <- sdm

if (B_TEST_DATA_AVAILABLE){
  test_data <- read.csv("training_sets/test_data_eco_seasons_B_Q4.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  cutoff <- get_threshold(sdm)
  if (CV_OR_RI=="ri"){
    metrics_B_Q4 <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
    metric_cis_B_Q4 <- get_sens_and_spec_ci(sdm, xtest, ytest, countrytest, cutoff)
  }else{
    metrics_B_Q4 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
    metric_cis_B_Q4 <- get_sens_and_spec_ci(sdm, xtest, ytest, NULL, cutoff)
  }
  cat("Dataset B Q1, test metrics:\n sens =",
      signif(metrics_B_Q1$sens, digits=2),
      " (",
      signif(metric_cis_B_Q1$sens, digits=2),
      "),\n spec =",
      signif(metrics_B_Q1$spec, digits=2),
      " (",
      signif(metric_cis_B_Q1$spec, digits=2),
      "),\n AUC =",
      signif(metrics_B_Q1$auc, digits=2),
      " (",
      signif(metric_cis_B_Q1$auc, digits=2),
      ")\n")
  cat("Dataset B Q2, test metrics:\n sens =",
      signif(metrics_B_Q2$sens, digits=2),
      " (",
      signif(metric_cis_B_Q2$sens, digits=2),
      "),\n spec =",
      signif(metrics_B_Q2$spec, digits=2),
      " (",
      signif(metric_cis_B_Q2$spec, digits=2),
      "),\n AUC =",
      signif(metrics_B_Q2$auc, digits=2),
      " (",
      signif(metric_cis_B_Q2$auc, digits=2),
      ")\n")
  cat("Dataset B Q3, test metrics:\n sens =",
      signif(metrics_B_Q3$sens, digits=2),
      " (",
      signif(metric_cis_B_Q3$sens, digits=2),
      "),\n spec =",
      signif(metrics_B_Q3$spec, digits=2),
      " (",
      signif(metric_cis_B_Q3$spec, digits=2),
      "),\n AUC =",
      signif(metrics_B_Q3$auc, digits=2),
      " (",
      signif(metric_cis_B_Q3$auc, digits=2),
      ")\n")
  cat("Dataset B Q4, test metrics:\n sens =",
      signif(metrics_B_Q4$sens, digits=2),
      " (",
      signif(metric_cis_B_Q4$sens, digits=2),
      "),\n spec =",
      signif(metrics_B_Q4$spec, digits=2),
      " (",
      signif(metric_cis_B_Q4$spec, digits=2),
      "),\n AUC =",
      signif(metrics_B_Q4$auc, digits=2),
      " (",
      signif(metric_cis_B_Q4$auc, digits=2),
      ")\n")
  
  if (SAVE_OUTPUTS){
    save(metrics_B_Q4,
         metric_cis_B_Q4,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_B_Q4.rds", sep = ""))
  }
}

if (CV_OR_RI=="cv"){
  # Generate risk map with percentiles
  pred_layers_B_Q4 <- predict(object = sdm,
                              x.layers = covstack,
                              quantiles = c(0.025, 0.975),
                              splitby = 20
  )
  names(pred_layers_B_Q4) <- c("Mean",
                               "Lower 95 percent confidence bound",
                               "Upper 95 percent confidence bound")
  if (SAVE_OUTPUTS){
    save(pred_layers_B_Q4,
         file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_B_Q4.rds", sep = ""))
  }
}


