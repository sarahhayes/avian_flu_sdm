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

# Restrict to quarter 1:
eco_layers <- eco_layers[[seq(1, nlyr(eco_layers), 4)]]
set.names(eco_layers, eco_lyrnames)

# Now do environmental:
env_paths <- list.files("../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental rasters/",
                        pattern = "*.tif",
                        full.names = TRUE)
env_paths <- env_paths[-grep("landcover_type1_full_raster", env_paths)]
env_paths <- env_paths[-grep("second_quart", env_paths)]
env_paths <- env_paths[-grep("third_quart", env_paths)]
env_paths <- env_paths[-grep("fourth_quart", env_paths)]
env_lyrnames <- env_paths %>%
  sub(pattern = "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/Environmental rasters//",
      replacement = "") %>%
  sub(pattern = ".tif",
      replacement = "")

env_layers <- rast(lapply(env_paths, rast))

first_quart_rast <- c(resample(eco_layers, env_layers), env_layers)
