 # In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors.
# This script includes cross-validation of BART parameters.

SAVE_FITS <- FALSE
SAVE_PLOTS <- FALSE

# Set path to folder containing data, and where output will be stored
PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(dbarts)
library(dplyr)
library(rasterVis)
library(embarcadero)
library(raster)
library(sp)
library(terra)
library(rworldmap)
set.seed(12345)

crossvalidate_from_known_folds <- function(data,
                                           folds,
                                           k_vals,
                                           power_vals,
                                           base_vals){
  cv_results <- data.frame(k=numeric(),
                           power=numeric(),
                           base=numeric(),
                           k1=numeric(),
                           k2=numeric(),
                           k3=numeric(),
                           k4=numeric(),
                           k5=numeric())
  cv_results[1:length(k_vals)*length(power_vals)*length(base_vals), ] <- 0
  for (i in 1:length(k_vals)){
    for (j in 1:length(power_vals)){
      for (k in 1:length(base_vals)){
      }
        for (fold in folds){
          model <- rbart_vi(y ~ . - ri,
                                  data[-fold,], # THIS NEEDS TO BE CHANGED!
                                  group.by = ri,
                                  group.by.test = ri,
                                  test = data[fold,], # AND THIS!
                                  k = k_vals[i],
                                  power = power_vals[j],
                                  base = base_vals[k],
                                  n.chains = 1L,
                                  n.threads = 1L,
                                  keepTrees = TRUE)
          mean_err <- mean(model$residuals)
          
        }
      }
    }
}

#### Period A Q1 ####

training_data <- read.csv("training_sets/training_data_A_Q1.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_A_Q1.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                    training_data,
                    group.by = ri,
                    group.by.test = ri,
                    test = test_data,
                    n.chains = 1L,
                    n.threads = 1L,
                    keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_A_Q1.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q1.rds", sep = ""))
}


#### Period A Q2 ####

training_data <- read.csv("training_sets/training_data_A_Q2.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_A_Q2.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_A_Q2.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q2.rds", sep = ""))
}

#### Period A Q3 ####

training_data <- read.csv("training_sets/training_data_A_Q3.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_A_Q3.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_A_Q3.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q3.rds", sep = ""))
}

#### Period A Q4 ####

training_data <- read.csv("training_sets/training_data_A_Q4.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_A_Q4.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_A_Q4.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_A_Q4.rds", sep = ""))
}

#### Period B Q1 ####

training_data <- read.csv("training_sets/training_data_B_Q1.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_B_Q1.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_B_Q1.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_B_Q1.rds", sep = ""))
}



#### Period B Q2 ####

training_data <- read.csv("training_sets/training_data_B_Q2.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_B_Q2.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_B_Q2.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_B_Q2.rds", sep = ""))
}


#### Period B Q3 ####

training_data <- read.csv("training_sets/training_data_B_Q3.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_B_Q3.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_B_Q3.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_B_Q3.rds", sep = ""))
}


#### Period B Q4 ####

training_data <- read.csv("training_sets/training_data_B_Q4.csv")
xtrain <- training_data %>% select(!("y"|"ri"))
ytrain <- training_data$y
countrytrain <- training_data$ri

test_data <- read.csv("training_sets/test_data_B_Q4.csv")
xtest <- test_data %>% select(!("y"|"ri"))
ytest <- test_data$y
countrytest <- test_data$ri

# Initialise model
basic_model <- rbart_vi(y ~ . - ri,
                        training_data,
                        group.by = ri,
                        group.by.test = ri,
                        test = test_data,
                        n.chains = 1L,
                        n.threads = 1L,
                        keepTrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)


if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_B_Q4.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 ri.data = training_data$ri,
                 # k = k_retune,
                 # power = power_retune,
                 # base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_with_vs_B_Q4.rds", sep = ""))
}