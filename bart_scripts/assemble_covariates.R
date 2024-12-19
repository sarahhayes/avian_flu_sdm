# In this script we assemble covariate layers for avian flu model

# Set path to folder containing data, and where output will be stored
PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(embarcadero)
library(raster)
library(terra)
set.seed(12345)

# set the crs we want to use
crs <- "epsg:3035"

euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later

eco_paths <- list.files(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Eco-Rasters-Eco-Seasons", sep = ""),
                        pattern = "*.tif",
                        full.names = TRUE)
eco_lyrnames <- eco_paths %>%
  sub(pattern = paste(PATH_TO_DATA, "AI_S2_SDM_storage/Eco-Rasters-Eco-Seasons/", sep = ""),
      replacement = "") %>%
  sub(pattern = "_rast.tif",
      replacement = "") %>%
  lapply(FUN = function(str){c(paste(str, "_q1", sep=""),
                               paste(str, "_q2", sep=""),
                               paste(str, "_q3", sep=""),
                               paste(str, "_q4", sep=""))}) %>%
  unlist()

eco_layers <- rast(lapply(eco_paths, rast))
set.names(eco_layers, eco_lyrnames)

# Separate by quarters:
q1_eco_layers <- eco_layers[[seq(1, nlyr(eco_layers), 4)]]
q2_eco_layers <- eco_layers[[seq(2, nlyr(eco_layers), 4)]]
q3_eco_layers <- eco_layers[[seq(3, nlyr(eco_layers), 4)]]
q4_eco_layers <- eco_layers[[seq(4, nlyr(eco_layers), 4)]]

# Now do environmental:
env_paths <- list.files(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Env-Rasters-Eco-Seasons", sep = ""),
                        pattern = "*.tif",
                        full.names = TRUE)
env_paths <- env_paths[-grep("landcover_output_full", env_paths)]
env_paths <- env_paths[-grep("chicken_density_2015", env_paths)]
env_paths <- env_paths[-grep("duck_density_2015", env_paths)]
# env_paths <- env_paths[-grep("mean_tmax", env_paths)]
# env_paths <- env_paths[-grep("mean_tmin", env_paths)]
# env_paths <- env_paths[-grep("isotherm_midnight", env_paths)]
env_paths <- env_paths[-grep("isotherm_min", env_paths)]
env_lyrnames <- env_paths %>%
  sub(pattern = paste(PATH_TO_DATA, "AI_S2_SDM_storage/Env-Rasters-Eco-Seasons/", sep = ""),
      replacement = "") %>%
  sub(pattern = "_10kres",
      replacement = "") %>%
  sub(pattern = ".tif",
      replacement = "")

env_layers <- rast(lapply(env_paths, rast))
names(env_layers) <- env_lyrnames

all_excludes <- grep("quart|_q", env_lyrnames, value = TRUE)
q1_excludes <- grep("first|q1", all_excludes, value = TRUE, invert = TRUE)
q2_excludes <- grep("second|q2", all_excludes, value = TRUE, invert = TRUE)
q3_excludes <- grep("third|q3", all_excludes, value = TRUE, invert = TRUE)
q4_excludes <- grep("fourth|q4", all_excludes, value = TRUE, invert = TRUE)

q1_env_layers <- subset(env_layers, setdiff(env_lyrnames, q1_excludes))
q2_env_layers <- subset(env_layers, setdiff(env_lyrnames, q2_excludes))
q3_env_layers <- subset(env_layers, setdiff(env_lyrnames, q3_excludes))
q4_env_layers <- subset(env_layers, setdiff(env_lyrnames, q4_excludes))

cov_coast <- read.csv(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Env-Rasters-Eco-Seasons/dist_to_coast_10kres.csv", sep = "")) %>%
  rename(x = X, y = Y) %>%
  dplyr::select(-X.1) %>% 
  relocate(x,y) %>%
  rast(type = "xyz", crs = "epsg:3035")
cov_water <- read.csv(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Env-Rasters-Eco-Seasons/dist_to_water_output_10kres.csv", sep = "")) %>%
  dplyr::select(-ID) %>% 
  relocate(x,y) %>%
  rast(type = "xyz", crs = "epsg:3035")

landcover_rast <- rast(paste(PATH_TO_DATA, "AI_S2_SDM_storage/Env-Rasters-Eco-Seasons/landcover_output_full_2022_10kres.tif", sep = ""))
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

all_covs <- c(resample(eco_layers, env_layers, method = "near"),
              env_layers,
              landcover_layers,
              cov_coast,
              cov_water)  
writeRaster(all_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates_eco_seasons/all_q_covs_10k.tif", sep = ""), overwrite = TRUE)
rm(all_covs)
gc()

q1_covs <- c(resample(q1_eco_layers, env_layers, method = "near"),
             q1_env_layers,
             landcover_layers,
             cov_coast,
             cov_water)  
writeRaster(q1_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates_eco_seasons/q1_covs_10k.tif", sep = ""), overwrite = TRUE)
rm(q1_covs)
gc()
q2_covs <- c(resample(q2_eco_layers, env_layers, method = "near"),
             q2_env_layers,
             landcover_layers,
             cov_coast,
             cov_water)
writeRaster(q2_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates_eco_seasons/q2_covs_10k.tif", sep = ""), overwrite = TRUE)
rm(q2_covs)
gc()
q3_covs <- c(resample(q3_eco_layers, env_layers, method = "near"),
             q3_env_layers,
             landcover_layers,
             cov_coast,
             cov_water)
writeRaster(q3_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates_eco_seasons/q3_covs_10k.tif", sep = ""), overwrite = TRUE)
rm(q3_covs)
gc()
q4_covs <- c(resample(q4_eco_layers, env_layers, method = "near"),
             q4_env_layers,
             landcover_layers,
             cov_coast,
             cov_water)
writeRaster(q4_covs, paste(PATH_TO_DATA, "AI_S2_SDM_storage/quarterly_covariates_eco_seasons/q4_covs_10k.tif", sep = ""), overwrite = TRUE)
rm(q4_covs)
gc()