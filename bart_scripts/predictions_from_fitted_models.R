# In this script we load the fitted models, calculate performance metrics, and
# perform Europe-wide predictions.

# Set path to folder containing data, and where output will be stored
PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

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

#### Quarter 1 without random intercept ####

test_data <- read.csv("training_sets/test_data_A_Q1.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_A_Q1.rds", sep = ""))
q1_no_ri_sdm <- sdm
cutoff <- get_threshold(sdm)
metrics_Q1_no_ri <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
cat("Q1, no RI: sens=",
    metrics_Q1_no_ri$sens,
    "spec=",
    metrics_Q1_no_ri$spec,
    "AUC=",
    metrics_Q1_no_ri$auc)

#### Quarter 1 with random intercept ####

load(file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q1.rds", sep = ""))
q1_ri_sdm <- sdm
cutoff <- get_threshold(sdm)
metrics_Q1_ri <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
cat("Q1, with RI: sens=",
    metrics_Q1_ri$sens,
    "spec=",
    metrics_Q1_ri$spec,
    "AUC=",
    metrics_Q1_ri$auc)

#### Quarter 2 without random intercept ####

test_data <- read.csv("training_sets/test_data_A_Q2.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_model_with_vs_A_Q2.rds", sep = ""))
q2_no_ri_sdm <- sdm
cutoff <- get_threshold(sdm)
metrics_Q2_no_ri <- get_sens_and_spec(sdm, xtest, ytest, NULL, cutoff)
cat("Q2, no RI: sens=",
    metrics_Q2_no_ri$sens,
    "spec=",
    metrics_Q2_no_ri$spec,
    "AUC=",
    metrics_Q2_no_ri$auc)

#### Quarter 2 with random intercept ####

load(file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q2.rds", sep = ""))
q2_ri_sdm <- sdm
cutoff <- get_threshold(sdm)
metrics_Q2_ri <- get_sens_and_spec(sdm, xtest, ytest, countrytest, cutoff)
cat("Q2, with RI: sens=",
    metrics_Q2_ri$sens,
    "spec=",
    metrics_Q2_ri$spec,
    "AUC=",
    metrics_Q2_ri$auc)
