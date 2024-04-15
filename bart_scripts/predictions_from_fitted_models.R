# In this script we load the fitted models, calculate performance metrics, and
# perform Europe-wide predictions.

# Optional command line arguments, must be passed as strings:
args <- commandArgs(trailingOnly = T)
if (length(args)<4){
  INCLUDE_CROSSTERMS <- "no-crossterms" # Set to "no-crossterms" to do model without crossterms or "with-crossterms" to do model with crossterms
}else{
  INCLUDE_CROSSTERMS <- as.logical(args[4])
}
if (length(args)<3){
  CV_OR_RI <- "cv" # Set to "cv" to do crossvalidated model or "ri" to do crossvalidation + random intercept model
}else{
  CV_OR_RI <- as.logical(args[3])
}
if (length(args)<2){
  # Set path to folder containing data, and where output will be stored
  PATH_TO_COVS <- "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/"
}else{
  PATH_TO_COVS <- args[2]
}
if (length(args)<1){
  # Set path to folder containing data, and where output will be stored
  PATH_TO_OUTPUTS <- "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/"
}else{
  PATH_TO_OUTPUTS <- args[1]
}

SAVE_OUTPUTS <- TRUE

B_TEST_DATA_AVAILABLE <- FALSE # Set to TRUE if we have test data for dataset B, otherwise this should be FALSE

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


#### Bring in year-round covariates ####

covstack <- raster::stack(paste(PATH_TO_COVS, "quarterly_covariates/all_q_covs_10k.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")

#### Dataset A Q1 ####

test_data <- read.csv("training_sets/test_data_A_Q1.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q1.rds", sep = ""))
sdm_A_Q1 <- sdm
cutoff <- get_threshold(sdm)
metrics_A_Q1 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
cat("Dataset A Q1, test metrics: sens =",
    metrics_A_Q1$sens,
    ", spec =",
    metrics_A_Q1$spec,
    ", AUC =",
    metrics_A_Q1$auc, "\n")

if (SAVE_OUTPUTS){
  save(metrics_A_Q1,
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_A_Q1.rds", sep = ""))
}

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
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_A_Q1.rds", sep = ""))
}


#### Dataset A Q2 ####

test_data <- read.csv("training_sets/test_data_A_Q2.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q2.rds", sep = ""))
sdm_A_Q2 <- sdm
cutoff <- get_threshold(sdm)
metrics_A_Q2 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
cat("Dataset A Q2, test metrics: sens =",
    metrics_A_Q2$sens,
    ", spec =",
    metrics_A_Q2$spec,
    ", AUC =",
    metrics_A_Q2$auc, "\n")

if (SAVE_OUTPUTS){
  save(metrics_A_Q2,
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_A_Q2.rds", sep = ""))
}

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
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_A_Q2.rds", sep = ""))
}


#### Dataset A Q3 ####

test_data <- read.csv("training_sets/test_data_A_Q3.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q3.rds", sep = ""))
sdm_A_Q3 <- sdm
cutoff <- get_threshold(sdm)
metrics_A_Q3 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
cat("Dataset A Q3, test metrics: sens =",
    metrics_A_Q3$sens,
    ", spec =",
    metrics_A_Q3$spec,
    ", AUC =",
    metrics_A_Q3$auc, "\n")

if (SAVE_OUTPUTS){
  save(metrics_A_Q3,
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_A_Q3.rds", sep = ""))
}

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
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_A_Q3.rds", sep = ""))
}


#### Dataset A Q4 ####

test_data <- read.csv("training_sets/test_data_A_Q4.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q4.rds", sep = ""))
sdm_A_Q4 <- sdm
cutoff <- get_threshold(sdm)
metrics_A_Q4 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
cat("Dataset A Q4, test metrics: sens =",
    metrics_A_Q4$sens,
    ", spec =",
    metrics_A_Q4$spec,
    ", AUC =",
    metrics_A_Q4$auc, "\n")

if (SAVE_OUTPUTS){
  save(metrics_A_Q4,
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_A_Q4.rds", sep = ""))
}

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
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_A_Q4.rds", sep = ""))
}


#### Dataset B Q1 ####

load(file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q1.rds", sep = ""))
sdm_B_Q1 <- sdm

if (B_TEST_DATA_AVAILABLE){
  test_data <- read.csv("training_sets/test_data_B_Q1.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  cutoff <- get_threshold(sdm)
  metrics_B_Q1 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
  cat("Dataset B Q1, test metrics: sens =",
      metrics_B_Q1$sens,
      ", spec =",
      metrics_B_Q1$spec,
      ", AUC =",
      metrics_B_Q1$auc, "\n")
  
  if (SAVE_OUTPUTS){
    save(metrics_B_Q1,
         file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_B_Q1.rds", sep = ""))
  }
}

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
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_B_Q1.rds", sep = ""))
}




#### Dataset B Q2 ####

load(file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q2.rds", sep = ""))
sdm_B_Q2 <- sdm

if (B_TEST_DATA_AVAILABLE){
  test_data <- read.csv("training_sets/test_data_B_Q2.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  cutoff <- get_threshold(sdm)
  metrics_B_Q2 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
  cat("Dataset B Q2, test metrics: sens =",
      metrics_B_Q2$sens,
      ", spec =",
      metrics_B_Q2$spec,
      ", AUC =",
      metrics_B_Q2$auc, "\n")
  
  if (SAVE_OUTPUTS){
    save(metrics_B_Q2,
         file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_B_Q2.rds", sep = ""))
  }
}

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
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_B_Q2.rds", sep = ""))
}





#### Dataset B Q3 ####

load(file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q3.rds", sep = ""))
sdm_B_Q3 <- sdm

if (B_TEST_DATA_AVAILABLE){
  test_data <- read.csv("training_sets/test_data_B_Q3.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  cutoff <- get_threshold(sdm)
  metrics_B_Q3 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
  cat("Dataset B Q3, test metrics: sens =",
      metrics_B_Q3$sens,
      ", spec =",
      metrics_B_Q3$spec,
      ", AUC =",
      metrics_B_Q3$auc, "\n")
  
  if (SAVE_OUTPUTS){
    save(metrics_B_Q3,
         file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_B_Q3.rds", sep = ""))
  }
}

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
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_B_Q3.rds", sep = ""))
}





#### Dataset B Q4 ####

load(file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q4.rds", sep = ""))
sdm_B_Q4 <- sdm

if (B_TEST_DATA_AVAILABLE){
  test_data <- read.csv("training_sets/test_data_B_Q4.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  cutoff <- get_threshold(sdm)
  metrics_B_Q4 <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
  cat("Dataset B Q4, test metrics: sens =",
      metrics_B_Q4$sens,
      ", spec =",
      metrics_B_Q4$spec,
      ", AUC =",
      metrics_B_Q4$auc, "\n")
  
  if (SAVE_OUTPUTS){
    save(metrics_B_Q4,
         file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_metrics_B_Q4.rds", sep = ""))
  }
}

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
       file = paste(PATH_TO_OUTPUTS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_predictions_B_Q4.rds", sep = ""))
}


