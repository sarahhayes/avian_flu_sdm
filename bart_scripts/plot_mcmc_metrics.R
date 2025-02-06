# In this script we plot the convergence behaviour of the MCMC used to fit BART

# In this script we load the fitted models, calculate performance metrics, and
# perform Europe-wide predictions.

SKIP_AQ3 <- TRUE # Set to TRUE to skip small period A quarter 3 dataset

# Optional command line arguments, must be passed as strings:
args <- commandArgs(trailingOnly = T)
if (length(args)<4){
  INCLUDE_CROSSTERMS <- "with-crossterms" # Set to "no-crossterms" to do model without crossterms or "with-crossterms" to do model with crossterms
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
  PATH_TO_OUTPUTS <- "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/"
}else{
  PATH_TO_OUTPUTS <- args[1]
}

# Decide whether to use models with ecological season boundaries
ECO_SEASONS <- TRUE
if (ECO_SEASONS){
  PATH_TO_MODELS <- paste(PATH_TO_OUTPUTS,
                          "fitted-BART-models-eco-seasons/",
                          sep = "")
}else{
  PATH_TO_MODELS <- PATH_TO_OUTPUTS
}

SAVE_OUTPUTS <- TRUE

B_TEST_DATA_AVAILABLE <- TRUE # Set to TRUE if we have test data for dataset B, otherwise this should be FALSE

library(dplyr)
library(embarcadero)
library(ggplot2)
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

# Get metrics for every particle in the posterior distribution
get_sens_and_spec_post <- function(sdm, xtest, ytest, ri, cutoff){
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
  return(data.frame("sample"=1:length(auc_by_model),
                  "sens"=sens_by_model,
                  "spec"=spec_by_model,
                  "tss"=tss_by_model,
                  "auc"=auc_by_model))
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
  metric_post_A_Q1 <- get_sens_and_spec_post(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metric_post_A_Q1 <- get_sens_and_spec_post(sdm, xtest, ytest, NULL, cutoff)
}

p_A_Q1 <- ggplot(metric_post_A_Q1, aes(x=sample)) +
  geom_line(aes(y=auc), colour = 'red') +
  geom_line(aes(y=sens), colour = 'blue') +
  geom_line(aes(y=spec), colour = 'green')

#### Dataset A Q2 ####

test_data <- read.csv("training_sets/test_data_eco_seasons_A_Q2.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q2.rds", sep = ""))
sdm_A_Q2 <- sdm
cutoff <- get_threshold(sdm)
if (CV_OR_RI=="ri"){
  metric_post_A_Q2 <- get_sens_and_spec_post(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metric_post_A_Q2 <- get_sens_and_spec_post(sdm, xtest, ytest, NULL, cutoff)
}

p_A_Q2 <- ggplot(metric_post_A_Q2, aes(x=sample)) +
  geom_line(aes(y=auc), colour = 'red') +
  geom_line(aes(y=sens), colour = 'blue') +
  geom_line(aes(y=spec), colour = 'green')

if (!SKIP_AQ3){
  #### Dataset A Q3 ####
  
  test_data <- read.csv("training_sets/test_data_eco_seasons_A_Q3.csv")
  xtest <- test_data %>% dplyr::select(!("y"|"ri"))
  ytest <- test_data$y
  countrytest <- test_data$ri
  
  load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_A_Q3.rds", sep = ""))
  sdm_A_Q3 <- sdm
  cutoff <- get_threshold(sdm)
  if (CV_OR_RI=="ri"){
    metric_post_A_Q3 <- get_sens_and_spec_post(sdm, xtest, ytest, countrytest, cutoff)
  }else{
    metric_post_A_Q3 <- get_sens_and_spec_post(sdm, xtest, ytest, NULL, cutoff)
  }
  
  p_A_Q3 <- ggplot(metric_post_A_Q3, aes(x=sample)) +
    geom_line(aes(y=auc), colour = 'red') +
    geom_line(aes(y=sens), colour = 'blue') +
    geom_line(aes(y=spec), colour = 'green')
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
  metric_post_A_Q4 <- get_sens_and_spec_post(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metric_post_A_Q4 <- get_sens_and_spec_post(sdm, xtest, ytest, NULL, cutoff)
}

p_A_Q4 <- ggplot(metric_post_A_Q4, aes(x=sample)) +
  geom_line(aes(y=auc), colour = 'red') +
  geom_line(aes(y=sens), colour = 'blue') +
  geom_line(aes(y=spec), colour = 'green')

#### Dataset B Q1 ####

test_data <- read.csv("training_sets/test_data_eco_seasons_B_Q1.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q1.rds", sep = ""))
sdm_B_Q1 <- sdm
cutoff <- get_threshold(sdm)
if (CV_OR_RI=="ri"){
  metric_post_B_Q1 <- get_sens_and_spec_post(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metric_post_B_Q1 <- get_sens_and_spec_post(sdm, xtest, ytest, NULL, cutoff)
}

p_B_Q1 <- ggplot(metric_post_B_Q1, aes(x=sample)) +
  geom_line(aes(y=auc), colour = 'red') +
  geom_line(aes(y=sens), colour = 'blue') +
  geom_line(aes(y=spec), colour = 'green')

#### Dataset B Q2 ####

test_data <- read.csv("training_sets/test_data_eco_seasons_B_Q2.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q2.rds", sep = ""))
sdm_B_Q2 <- sdm
cutoff <- get_threshold(sdm)
if (CV_OR_RI=="ri"){
  metric_post_B_Q2 <- get_sens_and_spec_post(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metric_post_B_Q2 <- get_sens_and_spec_post(sdm, xtest, ytest, NULL, cutoff)
}

p_B_Q2 <- ggplot(metric_post_B_Q2, aes(x=sample)) +
  geom_line(aes(y=auc), colour = 'red') +
  geom_line(aes(y=sens), colour = 'blue') +
  geom_line(aes(y=spec), colour = 'green')


#### Dataset B Q3 ####

test_data <- read.csv("training_sets/test_data_eco_seasons_B_Q3.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q3.rds", sep = ""))
sdm_B_Q3 <- sdm
cutoff <- get_threshold(sdm)
if (CV_OR_RI=="ri"){
  metric_post_B_Q3 <- get_sens_and_spec_post(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metric_post_B_Q3 <- get_sens_and_spec_post(sdm, xtest, ytest, NULL, cutoff)
}

p_B_Q3 <- ggplot(metric_post_B_Q3, aes(x=sample)) +
  geom_line(aes(y=auc), colour = 'red') +
  geom_line(aes(y=sens), colour = 'blue') +
  geom_line(aes(y=spec), colour = 'green')


#### Dataset B Q4 ####

test_data <- read.csv("training_sets/test_data_eco_seasons_B_Q4.csv")
xtest <- test_data %>% dplyr::select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

load(file = paste(PATH_TO_MODELS, "fitted-BART-models-", INCLUDE_CROSSTERMS,"/", CV_OR_RI, "_model_with_vs_B_Q4.rds", sep = ""))
sdm_B_Q4 <- sdm
cutoff <- get_threshold(sdm)
if (CV_OR_RI=="ri"){
  metric_post_B_Q4 <- get_sens_and_spec_post(sdm, xtest, ytest, countrytest, cutoff)
}else{
  metric_post_B_Q4 <- get_sens_and_spec_post(sdm, xtest, ytest, NULL, cutoff)
}

p_B_Q4 <- ggplot(metric_post_B_Q4, aes(x=sample)) +
  geom_line(aes(y=auc), colour = 'red') +
  geom_line(aes(y=sens), colour = 'blue') +
  geom_line(aes(y=spec), colour = 'green')

# Gather metrics for plotting
all_metrics <- cbind(metric_post_A_Q1 %>% rename_with(~paste0(., "_A_Q1"), -sample),
                     metric_post_A_Q2 %>% dplyr::select(-sample) %>% rename_with(~paste0(., "_A_Q2")),
                     metric_post_A_Q4 %>% dplyr::select(-sample) %>% rename_with(~paste0(., "_A_Q4")),
                     metric_post_B_Q1 %>% dplyr::select(-sample) %>% rename_with(~paste0(., "_B_Q1")),
                     metric_post_B_Q2 %>% dplyr::select(-sample) %>% rename_with(~paste0(., "_B_Q2")),
                     metric_post_B_Q3 %>% dplyr::select(-sample) %>% rename_with(~paste0(., "_B_Q3")),
                     metric_post_B_Q4 %>% dplyr::select(-sample) %>% rename_with(~paste0(., "_B_Q4"))) %>%
              pivot_longer(-sample)
auc_to_plot <- all_metrics %>% filter(grepl("auc", name))
p_auc <- ggplot(auc_to_plot, aes(x=sample, y=value, colour=name)) +
  ggtitle("AUC by MCMC sample") +
  geom_line()

sens_to_plot <- all_metrics %>% filter(grepl("sens", name))
p_sens <- ggplot(sens_to_plot, aes(x=sample, y=value, colour=name)) +
  ggtitle("Sensitivity by MCMC sample") +
  geom_line()

spec_to_plot <- all_metrics %>% filter(grepl("spec", name))
p_spec <- ggplot(spec_to_plot, aes(x=sample, y=value, colour=name)) +
  ggtitle("Specificity by MCMC sample") +
  geom_line()
