 # In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors

library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

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

# Restrict to quarter 2:
eco_layers <- eco_layers[[seq(2, nlyr(eco_layers), 4)]]
set.names(eco_layers, eco_lyrnames)

# Now do environmental:
env_paths <- list.files("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental rasters/",
                        pattern = "*.tif",
                        full.names = TRUE)
env_paths <- env_paths[-grep("landcover_type1_full_raster", env_paths)]
env_paths <- env_paths[-grep("first_quart", env_paths)]
env_paths <- env_paths[-grep("third_quart", env_paths)]
env_paths <- env_paths[-grep("fourth_quart", env_paths)]
env_lyrnames <- env_paths %>%
  sub(pattern = "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental rasters//",
      replacement = "") %>%
  sub(pattern = ".tif",
      replacement = "")

env_layers <- rast(lapply(env_paths, rast))
names(env_layers) <- gsub("6_Ch_2015_Aw", "chicken_density", names(env_layers))
names(env_layers) <- gsub("6_Dk_2015_Aw", "duck_density", names(env_layers))
names(env_layers) <- gsub("^mean$", "ndvi_second_quart", names(env_layers))

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

# landcover_rast <- rast("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental rasters/landcover_type1_full_raster.tif")
# crs(landcover_rast) <- crs

second_quart_covs <- c(resample(eco_layers, env_layers),
                       # landcover_rast,
                       env_layers,
                       cov_coast,
                       cov_water,
                       cov_alt)
rm(eco_layers,
   # landcover_rast,
   env_layers,
   cov_coast,
   cov_water,
   cov_alt)

writeRaster(second_quart_covs, "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/second_quart_covs.tif", overwrite = TRUE)
covstack <- raster::stack("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/quarterly_covariates/second_quart_covs.tif")

# Load training data
training_coords <- readRDS("training_sets/training_coords_Q2.RDS")
cov_df <- raster::extract(second_quart_covs, training_coords[, 1:2])
cov_df <- cov_df[, -1]

training <- sample(1:nrow(cov_df), 1600)
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

# Check performance manually:
ytrain_pos <- ytrain[which(complete.cases(xtrain))] == 1
yhat.train <- plogis(basic_model$yhat.train)
yhat.maj <- colSums(yhat.train) >= .36*nrow(yhat.train)
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
yhat.maj <- colSums(yhat.test) >= .4*nrow(yhat.test)

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

sdm <- bart.step(x.data = cov_df,
                 y.data = training_coords$pos,
                 full = TRUE,
                 quiet = TRUE)

euro_10k <- rast("output/euro_rast_10k.tif")
covstack_lores <- aggregate(covstack, fact = 10)

# Generate risk map
pred_layer <- predict(object = sdm,
                      x.layers = covstack_lores,
                      quantiles = c(0.025, 0.975),
                      splitby = 20
                      )

# Plot map
plot(pred_layer[[1]]>0.4,
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk')
plot(pred_layer[[2]]>0.4,
     box = FALSE,
     axes = FALSE,
     main = '2.5% posterior bound')
plot(pred_layer[[3]]>0.4,
     box = FALSE,
     axes = FALSE,
     main = '97.5% posterior bound')

plot(pred_layer[[1]],
     box = FALSE,
     axes = FALSE,
     main = 'Avian flu risk')
