 # In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors

SAVE_FITS <- TRUE
SAVE_PLOTS <- FALSE
BUILD_COVS <- FALSE
BUILD_DATA <- FALSE

# Set path to folder containing data, and where output will be stored
PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

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
  
  eco_paths <- list.files(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Eco-Rasters", sep = ""),
                          pattern = "*.tif",
                          full.names = TRUE)
  eco_lyrnames <- eco_paths %>%
    sub(pattern = paste(PATH_TO_DATA, "AI_S2_SDM_storage/Eco-Rasters/", sep = ""),
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
  env_paths <- list.files(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental rasters/10k_res", sep = ""),
                          pattern = "*.tif",
                          full.names = TRUE)
  env_paths <- env_paths[-grep("landcover_output_full", env_paths)]
  env_paths <- env_paths[-grep("chicken_density_2015", env_paths)]
  env_paths <- env_paths[-grep("duck_density_2015", env_paths)]
  env_paths <- env_paths[-grep("mean_tmax", env_paths)]
  env_paths <- env_paths[-grep("mean_tmin", env_paths)]
  env_paths <- env_paths[-grep("isotherm_midnight", env_paths)]
  env_paths <- env_paths[-grep("isotherm_min", env_paths)]
  env_lyrnames <- env_paths %>%
    sub(pattern = paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental rasters/10k_res/", sep = ""),
        replacement = "") %>%
    sub(pattern = "_10kres",
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
  
  cov_coast <- read.csv(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental rasters/10k_res/dist_to_coast_10kres.csv", sep = "")) %>%
    rename(x = X, y = Y) %>%
    select(-X.1) %>% 
    relocate(x,y) %>%
    rast(type = "xyz", crs = "epsg:3035")
  cov_water <- read.csv(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental rasters/10k_res/dist_to_water_output_10kres.csv", sep = "")) %>%
    select(-ID) %>% 
    relocate(x,y) %>%
    rast(type = "xyz", crs = "epsg:3035")
  
  landcover_rast <- rast(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Environmental rasters/10k_res/landcover_output_full_2022_10kres.tif", sep = ""))
  n_landtypes <- 17
  landcover_layers <- lapply(1:n_landtypes,
                             FUN = function(i){landcover_rast == i}) %>%
    rast()
  names(landcover_layers) <- lapply(1:n_landtypes,
                                    FUN = function(i){paste("lc_", i, sep = "")})
  lc_lookup <- pairlist("lc_1" = "Water_bodies",
                        "lc_2" = "Evergreen_Needleleaf_Forests",
                        "lc_3" = "Evergreen_Broadleaf_Forests",
                        "lc_4" = "Deciduous_Needleleaf_Forests",
                        "lc_5" = "Deciduous_Broadleaf_Forests",
                        "lc_6" = "Mixed_Forests",
                        "lc_7" = "Closed_Shrublands",
                        "lc_8" = "Open_Shrublands",
                        "lc_9" = "Woody_Savannas",
                        "lc_10" = "Savannas",
                        "lc_11" = "Grasslands",
                        "lc_12" = "Permanent_Wetlands",
                        "lc_13" = "Croplands",
                        "lc_14" = "Urban_and_Built-up_Lands",
                        "lc_15" = "Cropland/Natural_Vegetation_Mosaics",
                        "lc_16" = "Non-Vegetated_Lands",
                        "lc_17" = "Unclassified")
  names(landcover_layers) <- lapply(names(landcover_layers),
                                    FUN = function(n){lc_lookup[[n]]})
  # Remove layers that are either all True or all False:
  hom_layers <- sapply(1:n_landtypes,
                       FUN = function(i){diff(minmax(landcover_layers[[i]]))==0}) %>%
    which()
  if (length(hom_layers)>0){
    landcover_layers <- landcover_layers[[-hom_layers]]
  }
  
  q1_covs <- c(resample(q1_eco_layers, env_layers, method = "near"),
                         q1_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water)  
  writeRaster(q1_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""), overwrite = TRUE)
  rm(q1_covs)
  gc()
  q2_covs <- c(resample(q2_eco_layers, env_layers, method = "near"),
                         q2_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water)
  writeRaster(q2_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs_10k.tif", sep = ""), overwrite = TRUE)
  rm(q2_covs)
  gc()
  q3_covs <- c(resample(q3_eco_layers, env_layers, method = "near"),
                        q3_env_layers,
                        landcover_layers,
                        cov_coast,
                        cov_water)
  writeRaster(q3_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs_10k.tif", sep = ""), overwrite = TRUE)
  rm(q3_covs)
  gc()
  q4_covs <- c(resample(q4_eco_layers, env_layers, method = "near"),
                         q4_env_layers,
                         landcover_layers,
                         cov_coast,
                         cov_water)
  writeRaster(q4_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs_10k.tif", sep = ""), overwrite = TRUE)
  rm(q4_covs)
  gc()

}

################################################################################
# Do the A Q1 analysis

if (BUILD_DATA){
  covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""))
  
  # Drop unclassified land layer
  covstack <- dropLayer(covstack, "lc_17")
  
  # Load training data
  training_coords <- readRDS("training_sets/training_coords_A_Q1.RDS")
  ytrain <- training_coords$pos
  xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtrain)))
  xtrain <- xtrain[-bad_rows, ]
  ytrain <- ytrain[-bad_rows]
  
  test_coords <- readRDS("training_sets/test_coords_A_Q1.RDS")
  ytest <- test_coords$const
  xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtest)))
  xtest <- xtest[-bad_rows, ]
  ytest <- ytest[-bad_rows]
}else{
  training_data <- read.csv("training_sets/training_data_A_Q1.csv")
  xtrain <- training_data[, names(training_data) != "y"]
  ytrain <- training_data$y
  
  test_data <- read.csv("training_sets/test_data_A_Q1.csv")
  xtest <- test_data[, names(test_data) != "y"]
  ytest <- test_data$y
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

ypred <- colMeans(pnorm(basic_model$yhat.test))
basic_model_prec <- length(which(ypred>get_threshold(basic_model))) / length(ypred)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_A_Q1.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/sdm_A_Q1.rds", sep = ""))
}

# Only do projections here if we have built the dataset and have covstack in workspace
if (BUILD_DATA){

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
                      )
ypred <- raster::extract(pred_layer[[2]], test_coords[-bad_rows, 1:2])
sdm_prec <- length(which(ypred>get_threshold(sdm))) / length(ypred)
if (SAVE_FITS){
  save(pred_layer,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/prediction_Q1.rds", sep = ""))
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/q1_map.png")
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
  png(filename = "../../bartfit-plots/q1_lbound.png")
}
plot(pred_layer[[2]],
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound, Q1')
if (SAVE_PLOTS){
  dev.off()
}

if (SAVE_PLOTS){
  png(filename = "../../bartfit-plots/q1_ubound.png")
}
plot(pred_layer[[3]],
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound, Q1')
if (SAVE_PLOTS){
  dev.off()
}}

################################################################################
# Do the A Q2 analysis

if (BUILD_DATA){
  covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs_10k.tif", sep = ""))
  
  # Drop unclassified land layer
  covstack <- dropLayer(covstack, "lc_17")
  
  # Load training data
  training_coords <- readRDS("training_sets/training_coords_A_Q2.RDS")
  ytrain <- training_coords$pos
  xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtrain)))
  xtrain <- xtrain[-bad_rows, ]
  ytrain <- ytrain[-bad_rows]
  
  test_coords <- readRDS("training_sets/test_coords_A_Q2.RDS")
  ytest <- test_coords$const
  xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtest)))
  xtest <- xtest[-bad_rows, ]
  ytest <- ytest[-bad_rows]
}else{
  training_data <- read.csv("training_sets/training_data_A_Q2.csv")
  xtrain <- training_data[, names(training_data) != "y"]
  ytrain <- training_data$y
  
  test_data <- read.csv("training_sets/test_data_A_Q2.csv")
  xtest <- test_data[, names(test_data) != "y"]
  ytest <- test_data$y
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

ypred <- colMeans(pnorm(basic_model$yhat.test))
basic_model_prec <- length(which(ypred>get_threshold(basic_model))) / length(ypred)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_A_Q2.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/sdm_A_Q2.rds", sep = ""))
}

# Only do projections here if we have built the dataset and have covstack in workspace
if (BUILD_DATA){
  
  # Generate risk map
  pred_layer <- predict(object = sdm,
                        x.layers = covstack,
                        quantiles = c(0.025, 0.975),
                        splitby = 20
  )
  ypred <- raster::extract(pred_layer[[2]], test_coords[-bad_rows, 1:2])
  sdm_prec <- length(which(ypred>get_threshold(sdm))) / length(ypred)
  if (SAVE_FITS){
    save(pred_layer,
         file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/prediction_Q2.rds", sep = ""))
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q2_map.png")
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
    png(filename = "../../bartfit-plots/q2_lbound.png")
  }
  plot(pred_layer[[2]],
       box = FALSE,
       axes = FALSE,
       main = '2.5% posterior bound, Q2')
  if (SAVE_PLOTS){
    dev.off()
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q2_ubound.png")
  }
  plot(pred_layer[[3]],
       box = FALSE,
       axes = FALSE,
       main = '97.5% posterior bound, Q2')
  if (SAVE_PLOTS){
    dev.off()
  }}

################################################################################
# Do the A Q3 analysis

if (BUILD_DATA){
  covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs_10k.tif", sep = ""))
  
  # Drop unclassified land layer
  covstack <- dropLayer(covstack, "lc_17")
  
  # Load training data
  training_coords <- readRDS("training_sets/training_coords_A_Q3.RDS")
  ytrain <- training_coords$pos
  xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtrain)))
  xtrain <- xtrain[-bad_rows, ]
  ytrain <- ytrain[-bad_rows]
  
  test_coords <- readRDS("training_sets/test_coords_A_Q3.RDS")
  ytest <- test_coords$const
  xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtest)))
  xtest <- xtest[-bad_rows, ]
  ytest <- ytest[-bad_rows]
}else{
  training_data <- read.csv("training_sets/training_data_A_Q3.csv")
  xtrain <- training_data[, names(training_data) != "y"]
  ytrain <- training_data$y
  
  test_data <- read.csv("training_sets/test_data_A_Q3.csv")
  xtest <- test_data[, names(test_data) != "y"]
  ytest <- test_data$y
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

ypred <- colMeans(pnorm(basic_model$yhat.test))
basic_model_prec <- length(which(ypred>get_threshold(basic_model))) / length(ypred)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_A_Q3.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/sdm_A_Q3.rds", sep = ""))
}

# Only do projections here if we have built the dataset and have covstack in workspace
if (BUILD_DATA){
  
  # Generate risk map
  pred_layer <- predict(object = sdm,
                        x.layers = covstack,
                        quantiles = c(0.025, 0.975),
                        splitby = 20
  )
  ypred <- raster::extract(pred_layer[[2]], test_coords[-bad_rows, 1:2])
  sdm_prec <- length(which(ypred>get_threshold(sdm))) / length(ypred)
  if (SAVE_FITS){
    save(pred_layer,
         file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/prediction_Q3.rds", sep = ""))
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q3_map.png")
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
    png(filename = "../../bartfit-plots/q3_lbound.png")
  }
  plot(pred_layer[[2]],
       box = FALSE,
       axes = FALSE,
       main = '2.5% posterior bound, Q3')
  if (SAVE_PLOTS){
    dev.off()
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q3_ubound.png")
  }
  plot(pred_layer[[3]],
       box = FALSE,
       axes = FALSE,
       main = '97.5% posterior bound, Q3')
  if (SAVE_PLOTS){
    dev.off()
  }}

################################################################################
# Do the A Q4 analysis

if (BUILD_DATA){
  covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs_10k.tif", sep = ""))
  
  # Drop unclassified land layer
  covstack <- dropLayer(covstack, "lc_17")
  
  # Load training data
  training_coords <- readRDS("training_sets/training_coords_A_Q4.RDS")
  ytrain <- training_coords$pos
  xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtrain)))
  xtrain <- xtrain[-bad_rows, ]
  ytrain <- ytrain[-bad_rows]
  
  test_coords <- readRDS("training_sets/test_coords_A_Q4.RDS")
  ytest <- test_coords$const
  xtest <- data.frame(raster::extract(covstack, test_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtest)))
  xtest <- xtest[-bad_rows, ]
  ytest <- ytest[-bad_rows]
}else{
  training_data <- read.csv("training_sets/training_data_A_Q4.csv")
  xtrain <- training_data[, names(training_data) != "y"]
  ytrain <- training_data$y
  
  test_data <- read.csv("training_sets/test_data_A_Q4.csv")
  xtest <- test_data[, names(test_data) != "y"]
  ytest <- test_data$y
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    x.test = xtest,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

ypred <- colMeans(pnorm(basic_model$yhat.test))
basic_model_prec <- length(which(ypred>get_threshold(basic_model))) / length(ypred)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_A_Q4.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/sdm_A_Q4.rds", sep = ""))
}

# Only do projections here if we have built the dataset and have covstack in workspace
if (BUILD_DATA){
  
  # Generate risk map
  pred_layer <- predict(object = sdm,
                        x.layers = covstack,
                        quantiles = c(0.025, 0.975),
                        splitby = 20
  )
  ypred <- raster::extract(pred_layer[[2]], test_coords[-bad_rows, 1:2])
  sdm_prec <- length(which(ypred>get_threshold(sdm))) / length(ypred)
  if (SAVE_FITS){
    save(pred_layer,
         file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/prediction_Q4.rds", sep = ""))
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q4_map.png")
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
    png(filename = "../../bartfit-plots/q4_lbound.png")
  }
  plot(pred_layer[[2]],
       box = FALSE,
       axes = FALSE,
       main = '2.5% posterior bound, Q4')
  if (SAVE_PLOTS){
    dev.off()
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q4_ubound.png")
  }
  plot(pred_layer[[3]],
       box = FALSE,
       axes = FALSE,
       main = '97.5% posterior bound, Q4')
  if (SAVE_PLOTS){
    dev.off()
  }}

################################################################################
# Do the B Q1 analysis

if (BUILD_DATA){
  covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q1_covs_10k.tif", sep = ""))
  
  # Drop unclassified land layer
  covstack <- dropLayer(covstack, "lc_17")
  
  # Load training data
  training_coords <- readRDS("training_sets/training_coords_B_Q1.RDS")
  ytrain <- training_coords$pos
  xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtrain)))
  xtrain <- xtrain[-bad_rows, ]
  ytrain <- ytrain[-bad_rows]
}else{
  training_data <- read.csv("training_sets/training_data_B_Q1.csv")
  xtrain <- training_data[, names(training_data) != "y"]
  ytrain <- training_data$y
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_B_Q1.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/sdm_B_Q1.rds", sep = ""))
}

# Only do projections here if we have built the dataset and have covstack in workspace
if (BUILD_DATA){
  
  # Generate risk map
  pred_layer <- predict(object = sdm,
                        x.layers = covstack,
                        quantiles = c(0.025, 0.975),
                        splitby = 20
  )
  ypred <- raster::extract(pred_layer[[2]], test_coords[-bad_rows, 1:2])
  sdm_prec <- length(which(ypred>get_threshold(sdm))) / length(ypred)
  if (SAVE_FITS){
    save(pred_layer,
         file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/prediction_Q1.rds", sep = ""))
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q1_map.png")
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
    png(filename = "../../bartfit-plots/q1_lbound.png")
  }
  plot(pred_layer[[2]],
       box = FALSE,
       axes = FALSE,
       main = '2.5% posterior bound, Q1')
  if (SAVE_PLOTS){
    dev.off()
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q1_ubound.png")
  }
  plot(pred_layer[[3]],
       box = FALSE,
       axes = FALSE,
       main = '97.5% posterior bound, Q1')
  if (SAVE_PLOTS){
    dev.off()
  }}


################################################################################
# Do the B Q2 analysis

if (BUILD_DATA){
  covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q2_covs_10k.tif", sep = ""))
  
  # Drop unclassified land layer
  covstack <- dropLayer(covstack, "lc_17")
  
  # Load training data
  training_coords <- readRDS("training_sets/training_coords_B_Q2.RDS")
  ytrain <- training_coords$pos
  xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtrain)))
  xtrain <- xtrain[-bad_rows, ]
  ytrain <- ytrain[-bad_rows]
}else{
  training_data <- read.csv("training_sets/training_data_B_Q2.csv")
  xtrain <- training_data[, names(training_data) != "y"]
  ytrain <- training_data$y
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_B_Q2.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/sdm_B_Q2.rds", sep = ""))
}

# Only do projections here if we have built the dataset and have covstack in workspace
if (BUILD_DATA){
  
  # Generate risk map
  pred_layer <- predict(object = sdm,
                        x.layers = covstack,
                        quantiles = c(0.025, 0.975),
                        splitby = 20
  )
  ypred <- raster::extract(pred_layer[[2]], test_coords[-bad_rows, 1:2])
  sdm_prec <- length(which(ypred>get_threshold(sdm))) / length(ypred)
  if (SAVE_FITS){
    save(pred_layer,
         file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/prediction_Q2.rds", sep = ""))
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q2_map.png")
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
    png(filename = "../../bartfit-plots/q2_lbound.png")
  }
  plot(pred_layer[[2]],
       box = FALSE,
       axes = FALSE,
       main = '2.5% posterior bound, Q2')
  if (SAVE_PLOTS){
    dev.off()
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q2_ubound.png")
  }
  plot(pred_layer[[3]],
       box = FALSE,
       axes = FALSE,
       main = '97.5% posterior bound, Q2')
  if (SAVE_PLOTS){
    dev.off()
  }}


################################################################################
# Do the B Q3 analysis

if (BUILD_DATA){
  covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q3_covs_10k.tif", sep = ""))
  
  # Drop unclassified land layer
  covstack <- dropLayer(covstack, "lc_17")
  
  # Load training data
  training_coords <- readRDS("training_sets/training_coords_B_Q3.RDS")
  ytrain <- training_coords$pos
  xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtrain)))
  xtrain <- xtrain[-bad_rows, ]
  ytrain <- ytrain[-bad_rows]
}else{
  training_data <- read.csv("training_sets/training_data_B_Q3.csv")
  xtrain <- training_data[, names(training_data) != "y"]
  ytrain <- training_data$y
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_B_Q3.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/sdm_B_Q3.rds", sep = ""))
}

# Only do projections here if we have built the dataset and have covstack in workspace
if (BUILD_DATA){
  
  # Generate risk map
  pred_layer <- predict(object = sdm,
                        x.layers = covstack,
                        quantiles = c(0.025, 0.975),
                        splitby = 20
  )
  ypred <- raster::extract(pred_layer[[2]], test_coords[-bad_rows, 1:2])
  sdm_prec <- length(which(ypred>get_threshold(sdm))) / length(ypred)
  if (SAVE_FITS){
    save(pred_layer,
         file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/prediction_Q3.rds", sep = ""))
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q3_map.png")
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
    png(filename = "../../bartfit-plots/q3_lbound.png")
  }
  plot(pred_layer[[2]],
       box = FALSE,
       axes = FALSE,
       main = '2.5% posterior bound, Q3')
  if (SAVE_PLOTS){
    dev.off()
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q3_ubound.png")
  }
  plot(pred_layer[[3]],
       box = FALSE,
       axes = FALSE,
       main = '97.5% posterior bound, Q3')
  if (SAVE_PLOTS){
    dev.off()
  }}


################################################################################
# Do the B Q4 analysis

if (BUILD_DATA){
  covstack <- raster::stack(paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates/q4_covs_10k.tif", sep = ""))
  
  # Drop unclassified land layer
  covstack <- dropLayer(covstack, "lc_17")
  
  # Load training data
  training_coords <- readRDS("training_sets/training_coords_B_Q4.RDS")
  ytrain <- training_coords$pos
  xtrain <- data.frame(raster::extract(covstack, training_coords[, 1:2]))
  bad_rows <- which(is.na(rowSums(xtrain)))
  xtrain <- xtrain[-bad_rows, ]
  ytrain <- ytrain[-bad_rows]
}else{
  training_data <- read.csv("training_sets/training_data_B_Q4.csv")
  xtrain <- training_data[, names(training_data) != "y"]
  ytrain <- training_data$y
}

# Initialise model
basic_model <- bart(xtrain,
                    ytrain,
                    keeptrees = TRUE)
invisible(basic_model$fit$state)
summary(basic_model)

if (SAVE_FITS){
  save(basic_model,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/basic_model_B_Q4.rds", sep = ""))
}

sdm <- bart.step(x.data = xtrain,
                 y.data = ytrain,
                 full = TRUE,
                 quiet = TRUE)
invisible(sdm$fit$state)
summary(sdm)
if (SAVE_FITS){
  save(sdm,
       file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/sdm_B_Q4.rds", sep = ""))
}

# Only do projections here if we have built the dataset and have covstack in workspace
if (BUILD_DATA){
  
  # Generate risk map
  pred_layer <- predict(object = sdm,
                        x.layers = covstack,
                        quantiles = c(0.025, 0.975),
                        splitby = 20
  )
  ypred <- raster::extract(pred_layer[[2]], test_coords[-bad_rows, 1:2])
  sdm_prec <- length(which(ypred>get_threshold(sdm))) / length(ypred)
  if (SAVE_FITS){
    save(pred_layer,
         file = paste(PATH_TO_DATA, "AI_S2_SDM_storage/fitted-BART-models/prediction_Q4.rds", sep = ""))
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q4_map.png")
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
    png(filename = "../../bartfit-plots/q4_lbound.png")
  }
  plot(pred_layer[[2]],
       box = FALSE,
       axes = FALSE,
       main = '2.5% posterior bound, Q4')
  if (SAVE_PLOTS){
    dev.off()
  }
  
  if (SAVE_PLOTS){
    png(filename = "../../bartfit-plots/q4_ubound.png")
  }
  plot(pred_layer[[3]],
       box = FALSE,
       axes = FALSE,
       main = '97.5% posterior bound, Q4')
  if (SAVE_PLOTS){
    dev.off()
  }}
