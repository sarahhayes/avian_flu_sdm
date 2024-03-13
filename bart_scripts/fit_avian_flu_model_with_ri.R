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

# Eyeball RI terms:
ranef_means <- data.frame(name = names(basic_model$ranef.mean),
                          rem = unname(basic_model$ranef.mean))
ggplot(ranef_means) + geom_col(aes(rem, name)) + theme(axis.title=element_blank())

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_Q1.rds", sep = ""))
}

# retuned_model <- retune(basic_model)
# 
# k_retune = retuned_model$fit$model@node.hyperprior@k
# power_retune = retuned_model$fit$model@tree.prior@power
# base_retune = retuned_model$fit$model@tree.prior@base

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
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_sdm_Q1.rds", sep = ""))
}

# Eyeball RI terms:
vs_ranef_means <- data.frame(name = names(sdm$ranef.mean),
                          rem = unname(sdm$ranef.mean))
p <- ggplot(vs_ranef_means) + geom_col(aes(rem, name)) + theme(axis.title=element_blank())
p

covstack <- aggregate(covstack, fact = 10)
gc()

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack,
                      ri.data = country_rast,
                      ri.name = "country_rast",
                      quantiles = c(0.025, 0.975),
                      splitby = 20
                      )
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_prediction_Q1.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q1_map.png")
}
# Plot map
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q1')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q1_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q1')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q1_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q1')
if (SAVE_PLOTS){
  dev.off()
}

################################################################################
# Do the Q2 analysis

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")


# Load training data
training_coords <- readRDS("training_sets/training_coords_Q2.RDS")
cov_df <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]

if (SAVE_FITS){
  save(xtrain, ytrain, xtest, ytest, file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/q2_train_test_data.rds", sep = ""))
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_Q2.rds", sep = ""))
}

retuned_model <- retune(basic_model)

k_retune = retuned_model$fit$model@node.hyperprior@k
power_retune = retuned_model$fit$model@tree.prior@power
base_retune = retuned_model$fit$model@tree.prior@base

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = k_retune,
                 power = power_retune,
                 base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_sdm_Q2.rds", sep = ""))
}

covstack_lores <- aggregate(covstack, fact = 10)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_prediction_Q2.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q2_map.png")
}
# Plot map
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q2')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q2_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q2')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q2_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q2')
if (SAVE_PLOTS){
  dev.off()
}

################################################################################
# Do the Q3 analysis

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")


# Load training data
training_coords <- readRDS("training_sets/training_coords_Q3.RDS")
cov_df <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]

if (SAVE_FITS){
  save(xtrain, ytrain, xtest, ytest, file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/q3_train_test_data.rds", sep = ""))
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_Q3.rds", sep = ""))
}

retuned_model <- retune(basic_model)

k_retune = retuned_model$fit$model@node.hyperprior@k
power_retune = retuned_model$fit$model@tree.prior@power
base_retune = retuned_model$fit$model@tree.prior@base

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = k_retune,
                 power = power_retune,
                 base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_sdm_Q3.rds", sep = ""))
}

covstack_lores <- aggregate(covstack, fact = 10)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_prediction_Q3.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q3_map.png")
}
# Plot map
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q3')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q3_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q3')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q3_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q3')
if (SAVE_PLOTS){
  dev.off()
}

################################################################################
# Do the Q4 analysis

covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs.tif", sep = ""))

# Drop unclassified land layer
covstack <- dropLayer(covstack, "lc_17")


# Load training data
training_coords <- readRDS("training_sets/training_coords_Q4.RDS")
cov_df <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]

if (SAVE_FITS){
  save(xtrain, ytrain, xtest, ytest, file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/q4_train_test_data.rds", sep = ""))
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/ri_model_Q4.rds", sep = ""))
}

retuned_model <- retune(basic_model)

k_retune = retuned_model$fit$model@node.hyperprior@k
power_retune = retuned_model$fit$model@tree.prior@power
base_retune = retuned_model$fit$model@tree.prior@base

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 k = k_retune,
                 power = power_retune,
                 base = base_retune,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_sdm_Q4.rds", sep = ""))
}

covstack_lores <- aggregate(covstack, fact = 10)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/cv_prediction_Q4.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q4_map.png")
}
# Plot map
plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Mean prediction, Q4')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q4_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q4')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/cv_q4_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q4')
if (SAVE_PLOTS){
  dev.off()
}