 # In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors

SAVE_FITS <- FALSE
BUILD_COVS <- TRUE

library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

# Function for directly getting optimal cutoff
get_threshold <- function(object){
  fitobj <- object$fit
  
  true.vector <- fitobj$data@y 
  
  pred <- prediction(colMeans(pnorm(object$yhat.train)), true.vector)
  
  perf.tss <- performance(pred,"sens","spec")
  tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
  tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
  thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))])
  return(thresh)
}

if (BUILD_COVS){
  # set the crs we want to use
  crs <- "epsg:3035"
  
  euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later
  
  eco_paths <- list.files("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Eco-Rasters",
                              pattern = "*.tif",
                              full.names = TRUE)
  eco_lyrnames <- eco_paths %>%
                  sub(pattern = "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Eco-Rasters/",
                      replacement = "") %>%
                  sub(pattern = "_rast.tif",
                      replacement = "")
  
  eco_layers <- rast(lapply(eco_paths, rast))
  
  # Separate by quarters:
  q1_eco_layers <- eco_layers[[seq(1, nlyr(eco_layers), 4)]]
  q2_eco_layers <- eco_layers[[seq(2, nlyr(eco_layers), 4)]]
  q3_eco_layers <- eco_layers[[seq(3, nlyr(eco_layers), 4)]]
  q4_eco_layers <- eco_layers[[seq(4, nlyr(eco_layers), 4)]]
  set.names(q1_eco_layers, eco_lyrnames)
  set.names(q2_eco_layers, eco_lyrnames)
  set.names(q3_eco_layers, eco_lyrnames)
  set.names(q4_eco_layers, eco_lyrnames)
  
  # Now do environmental:
  env_paths <- list.files("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental rasters/",
                          pattern = "*.tif",
                          full.names = TRUE)
  env_paths <- env_paths[-grep("landcover_output_full", env_paths)]
  env_paths <- env_paths[-grep("duck_density_2015", env_paths)]
  env_lyrnames <- env_paths %>%
    sub(pattern = "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental rasters//",
        replacement = "") %>%
    sub(pattern = ".tif",
        replacement = "")
  
  env_layers <- rast(lapply(env_paths, rast))
  names(env_layers) <- env_lyrnames
  
  all_excludes <- grep("quart", env_lyrnames, value = TRUE)
  q1_excludes <- grep("first", all_excludes, value = TRUE, invert = TRUE)
  q2_excludes <- grep("second", all_excludes, value = TRUE, invert = TRUE)
  q3_excludes <- grep("third", all_excludes, value = TRUE, invert = TRUE)
  q4_excludes <- grep("fourth", all_excludes, value = TRUE, invert = TRUE)
  
  q1_env_layers <- subset(env_layers, setdiff(env_lyrnames, q1_excludes))
  q2_env_layers <- subset(env_layers, setdiff(env_lyrnames, q2_excludes))
  q3_env_layers <- subset(env_layers, setdiff(env_lyrnames, q3_excludes))
  q4_env_layers <- subset(env_layers, setdiff(env_lyrnames, q4_excludes))
  
  cov_coast <- read.csv("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental variable csvs/dist_to_coast_output.csv") %>%
    rename(x = X, y = Y) %>%
    select(-X.1) %>% 
    relocate(x,y) %>%
    rast(type = "xyz", crs = "epsg:3035")
  cov_water <- read.csv("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental variable csvs/dist_to_water_output.csv") %>%
    select(-ID) %>% 
    relocate(x,y) %>%
    rast(type = "xyz", crs = "epsg:3035")
  cov_alt <- read.csv("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental variable csvs/elevation_outputs.csv") %>%
    select(-ID, -X) %>% 
    relocate(x,y) %>%
    rast(type = "xyz", crs = "epsg:3035")
  
  landcover_rast <- rast("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental rasters/landcover_output_full.tif")
  n_landtypes <- 17
  landcover_layers <- lapply(1:17, FUN = function(i){landcover_rast == i}) %>%
    rast()
  names(landcover_layers) <- lapply(1:17,
                                    FUN = function(i){paste("lc_", i, sep = "")})
  
  q1_covs <- c(resample(q1_eco_layers, env_layers),
                         q1_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water,
                         cov_alt)
  q2_covs <- c(resample(q2_eco_layers, env_layers),
                         q2_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water,
                         cov_alt)
  q3_covs <- c(resample(q3_eco_layers, env_layers),
                        q3_env_layers,
                        landcover_layers,
                        cov_coast,
                        cov_water,
                        cov_alt)
  q4_covs <- c(resample(q4_eco_layers, env_layers),
                         q4_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water,
                         cov_alt)
  
  writeRaster(q1_covs, "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q1_covs.tif", overwrite = TRUE)
  writeRaster(q2_covs, "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q2_covs.tif", overwrite = TRUE)
  writeRaster(q3_covs, "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q2_covs.tif", overwrite = TRUE)
  writeRaster(q4_covs, "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q2_covs.tif", overwrite = TRUE)
}

################################################################################
# Do the Q1 analysis

covstack <- raster::stack("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q1_covs.tif")

# Load training data
training_coords <- readRDS("training_sets/training_coords_Q1.RDS")
cov_df <- data.frame(raster::extract(covstack, training_coords[, 1:2]))

n_pts <- nrow(training_coords)
n_training <- round(.75 * n_pts)

training <- sample(1:nrow(cov_df), n_training)
xtrain <- cov_df[training, ]
ytrain <- training_coords$pos[training]
xtest <- cov_df[-training, ]
ytest <- training_coords$pos[-training]

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = "../../fitted-BART-models/basic_model_Q1.rds")
}

cutoff <- get_threshold(basic_model)

# Check performance manually:
ytrain_pos <- ytrain[which(complete.cases(xtrain))] == 1
yhat.train <- plogis(basic_model$yhat.train)
yhat.maj <- colSums(yhat.train) >= cutoff*nrow(yhat.train)
false_neg_rate <- length(which((!yhat.maj)&ytrain_pos))/length(which(ytrain_pos))
false_pos_rate <- length(which(yhat.maj&!ytrain_pos))/length(which(!ytrain_pos))
misclass_rate <- (length(which((!yhat.maj)&ytrain_pos)) +
                    length(which(yhat.maj&!ytrain_pos))) / length(ytrain_pos)
cat("Train false negative rate is",
    false_neg_rate,
    ".\n Train false positive rate is",
    false_pos_rate,
    ".\n Train misclassification rate is",
    misclass_rate,
    ".\n")

# And for test data:
ytest_pos <- ytest[which(complete.cases(xtest))] == 1
yhat.test <- plogis(basic_model$yhat.test)
yhat.maj <- colSums(yhat.test) >= cutoff*nrow(yhat.test)

# Calculate error rates
false_neg_rate <- length(which((!yhat.maj)&ytest_pos))/length(which(ytest_pos))
false_pos_rate <- length(which(yhat.maj&!ytest_pos))/length(which(!ytest_pos))
misclass_rate <- (length(which((!yhat.maj)&ytest_pos)) +
                    length(which(yhat.maj&!ytest_pos))) / length(ytest_pos)
cat("Test false negative rate is",
    false_neg_rate,
    ".\n Test false positive rate is",
    false_pos_rate,
    ".\n Test misclassification rate is",
    misclass_rate,
    ".\n")

covstack_lores <- aggregate(covstack, fact = 10)


sdm <- bart.step(x.data = cov_df,
                 y.data = training_coords$pos,
                 full = TRUE,
                 quiet = TRUE)
if (SAVE_FITS){
  save(sdm,
       file = "../../fitted-BART-models/sdm_Q1.rds")
}

cutoff <- get_threshold(basic_model)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
                      )
if (SAVE_FITS){
  save(pred_layer,
       file = "../../fitted-BART-models/prediction_Q1.rds")
}

# Plot map
plot(pred_layer[[1]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk, Q1')
plot(pred_layer[[2]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q1')
plot(pred_layer[[3]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q1')

plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk, Q1')

# Extract names of variables in optimised model:
sdm_vars <- dimnames(sdm$fit$data@x)[[2]]

layers_for_spartial <- grep("lc_",
                            sdm_vars,
                            value = TRUE,
                            invert = TRUE)

pd_maps <- spartial(sdm,
               covstack_lores,
               layers_for_spartial)
if (SAVE_FITS){
  save(pred_layer,
       file = "../../fitted-BART-models/spartial_Q1.rds")
}
plot(pd_maps)

################################################################################
# Do the Q2 analysis

covstack <- raster::stack("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q2_covs.tif")

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

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = "../../fitted-BART-models/basic_model_Q2.rds")
}

cutoff <- get_threshold(basic_model)

# Check performance manually:
ytrain_pos <- ytrain[which(complete.cases(xtrain))] == 1
yhat.train <- plogis(basic_model$yhat.train)
yhat.maj <- colSums(yhat.train) >= cutoff*nrow(yhat.train)
false_neg_rate <- length(which((!yhat.maj)&ytrain_pos))/length(which(ytrain_pos))
false_pos_rate <- length(which(yhat.maj&!ytrain_pos))/length(which(!ytrain_pos))
misclass_rate <- (length(which((!yhat.maj)&ytrain_pos)) +
                    length(which(yhat.maj&!ytrain_pos))) / length(ytrain_pos)
cat("Train false negative rate is",
    false_neg_rate,
    ".\n Train false positive rate is",
    false_pos_rate,
    ".\n Train misclassification rate is",
    misclass_rate,
    ".\n")

# And for test data:
ytest_pos <- ytest[which(complete.cases(xtest))] == 1
yhat.test <- plogis(basic_model$yhat.test)
yhat.maj <- colSums(yhat.test) >= cutoff*nrow(yhat.test)

# Calculate error rates
false_neg_rate <- length(which((!yhat.maj)&ytest_pos))/length(which(ytest_pos))
false_pos_rate <- length(which(yhat.maj&!ytest_pos))/length(which(!ytest_pos))
misclass_rate <- (length(which((!yhat.maj)&ytest_pos)) +
                    length(which(yhat.maj&!ytest_pos))) / length(ytest_pos)
cat("Test false negative rate is",
    false_neg_rate,
    ".\n Test false positive rate is",
    false_pos_rate,
    ".\n Test misclassification rate is",
    misclass_rate,
    ".\n")

covstack_lores <- aggregate(covstack, fact = 10)

sdm <- bart.step(x.data = cov_df,
                 y.data = training_coords$pos,
                 full = TRUE,
                 quiet = TRUE)
if (SAVE_FITS){
  save(sdm,
       file = "../../fitted-BART-models/sdm_Q2.rds")
}

cutoff <- get_threshold(basic_model)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = "../../fitted-BART-models/prediction_Q2.rds")
}

# Plot map
plot(pred_layer[[1]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk, Q2')
plot(pred_layer[[2]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q2')
plot(pred_layer[[3]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q2')

plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk, Q2')

# Extract names of variables in optimised model:
sdm_vars <- dimnames(sdm$fit$data@x)[[2]]

layers_for_spartial <- grep("lc_",
                            sdm_vars,
                            value = TRUE,
                            invert = TRUE)

pd_maps <- spartial(sdm,
                    covstack_lores,
                    layers_for_spartial)
if (SAVE_FITS){
  save(pred_layer,
       file = "../../fitted-BART-models/spartial_Q2.rds")
}
plot(pd_maps)


################################################################################
# Do the Q3 analysis

covstack <- raster::stack("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q3_covs.tif")

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

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = "../../fitted-BART-models/basic_model_Q3.rds")
}

cutoff <- get_threshold(basic_model)

# Check performance manually:
ytrain_pos <- ytrain[which(complete.cases(xtrain))] == 1
yhat.train <- plogis(basic_model$yhat.train)
yhat.maj <- colSums(yhat.train) >= cutoff*nrow(yhat.train)
false_neg_rate <- length(which((!yhat.maj)&ytrain_pos))/length(which(ytrain_pos))
false_pos_rate <- length(which(yhat.maj&!ytrain_pos))/length(which(!ytrain_pos))
misclass_rate <- (length(which((!yhat.maj)&ytrain_pos)) +
                    length(which(yhat.maj&!ytrain_pos))) / length(ytrain_pos)
cat("Train false negative rate is",
    false_neg_rate,
    ".\n Train false positive rate is",
    false_pos_rate,
    ".\n Train misclassification rate is",
    misclass_rate,
    ".\n")

# And for test data:
ytest_pos <- ytest[which(complete.cases(xtest))] == 1
yhat.test <- plogis(basic_model$yhat.test)
yhat.maj <- colSums(yhat.test) >= cutoff*nrow(yhat.test)

# Calculate error rates
false_neg_rate <- length(which((!yhat.maj)&ytest_pos))/length(which(ytest_pos))
false_pos_rate <- length(which(yhat.maj&!ytest_pos))/length(which(!ytest_pos))
misclass_rate <- (length(which((!yhat.maj)&ytest_pos)) +
                    length(which(yhat.maj&!ytest_pos))) / length(ytest_pos)
cat("Test false negative rate is",
    false_neg_rate,
    ".\n Test false positive rate is",
    false_pos_rate,
    ".\n Test misclassification rate is",
    misclass_rate,
    ".\n")

covstack_lores <- aggregate(covstack, fact = 10)


sdm <- bart.step(x.data = cov_df,
                 y.data = training_coords$pos,
                 full = TRUE,
                 quiet = TRUE)
if (SAVE_FITS){
  save(sdm,
       file = "../../fitted-BART-models/sdm_Q3.rds")
}

cutoff <- get_threshold(basic_model)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = "../../fitted-BART-models/prediction_Q3.rds")
}

# Plot map
plot(pred_layer[[1]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk, Q3')
plot(pred_layer[[2]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q3')
plot(pred_layer[[3]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q3')

plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk, Q3')

# Extract names of variables in optimised model:
sdm_vars <- dimnames(sdm$fit$data@x)[[2]]

layers_for_spartial <- grep("lc_",
                            sdm_vars,
                            value = TRUE,
                            invert = TRUE)

pd_maps <- spartial(sdm,
                    covstack_lores,
                    layers_for_spartial)
if (SAVE_FITS){
  save(pred_layer,
       file = "../../fitted-BART-models/spartial_Q3.rds")
}
plot(pd_maps)


################################################################################
# Do the Q4 analysis

covstack <- raster::stack("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/q4_covs.tif")

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

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = "../../fitted-BART-models/basic_model_Q4.rds")
}

cutoff <- get_threshold(basic_model)

# Check performance manually:
ytrain_pos <- ytrain[which(complete.cases(xtrain))] == 1
yhat.train <- plogis(basic_model$yhat.train)
yhat.maj <- colSums(yhat.train) >= cutoff*nrow(yhat.train)
false_neg_rate <- length(which((!yhat.maj)&ytrain_pos))/length(which(ytrain_pos))
false_pos_rate <- length(which(yhat.maj&!ytrain_pos))/length(which(!ytrain_pos))
misclass_rate <- (length(which((!yhat.maj)&ytrain_pos)) +
                    length(which(yhat.maj&!ytrain_pos))) / length(ytrain_pos)
cat("Train false negative rate is",
    false_neg_rate,
    ".\n Train false positive rate is",
    false_pos_rate,
    ".\n Train misclassification rate is",
    misclass_rate,
    ".\n")

# And for test data:
ytest_pos <- ytest[which(complete.cases(xtest))] == 1
yhat.test <- plogis(basic_model$yhat.test)
yhat.maj <- colSums(yhat.test) >= cutoff*nrow(yhat.test)

# Calculate error rates
false_neg_rate <- length(which((!yhat.maj)&ytest_pos))/length(which(ytest_pos))
false_pos_rate <- length(which(yhat.maj&!ytest_pos))/length(which(!ytest_pos))
misclass_rate <- (length(which((!yhat.maj)&ytest_pos)) +
                    length(which(yhat.maj&!ytest_pos))) / length(ytest_pos)
cat("Test false negative rate is",
    false_neg_rate,
    ".\n Test false positive rate is",
    false_pos_rate,
    ".\n Test misclassification rate is",
    misclass_rate,
    ".\n")

covstack_lores <- aggregate(covstack, fact = 10)

sdm <- bart.step(x.data = cov_df,
                 y.data = training_coords$pos,
                 full = TRUE,
                 quiet = TRUE)
if (SAVE_FITS){
  save(sdm,
       file = "../../fitted-BART-models/sdm_Q4.rds")
}

cutoff <- get_threshold(basic_model)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
)
if (SAVE_FITS){
  save(pred_layer,
       file = "../../fitted-BART-models/prediction_Q4.rds")
}

# Plot map
plot(pred_layer[[1]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk, Q4')
plot(pred_layer[[2]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q4')
plot(pred_layer[[3]]>cutoff,
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q4')

plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk, Q4')

# Extract names of variables in optimised model:
sdm_vars <- dimnames(sdm$fit$data@x)[[2]]

layers_for_spartial <- grep("lc_",
                            sdm_vars,
                            value = TRUE,
                            invert = TRUE)

pd_maps <- spartial(sdm,
                    covstack_lores,
                    layers_for_spartial)
if (SAVE_FITS){
  save(pred_layer,
       file = "../../fitted-BART-models/spartial_Q4.rds")
}
plot(pd_maps)
