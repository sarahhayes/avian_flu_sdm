## In this script we are taking the raster layers and extracting the 
## values from each of the cells into a csv for feeding into the model

rm(list = ls())

library(tidyverse)
library(terra)

# start with one file

prec_q1 <- terra::rast("data/variables/climate/climate_prepped/mean_prec_first_quart.tif")
prec_q1 # using this to check extent, resolution etc.

plot(prec_q1)

## make the blank raster
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later
blank_3035 <- rast(crs=crs, extent=euro_ext, res = 1000)
blank_3035

# now we want to make a points object using the centre of each pixel.
points_3035 <- terra::as.points(blank_3035)
points_3035

tictoc::tic()
rast_res <- terra::extract(prec_q1, points_3035, method = "simple", xy = T)
tictoc::toc()

# Took 1 min to do

# I believe you can make a raster stack and do the process in a single step. 
# make lists of files so can import simultaneously

climate_files <- c("data/variables/climate/climate_prepped/mean_prec_first_quart.tif",
                   "data/variables/climate/climate_prepped/mean_prec_second_quart.tif",
                   "data/variables/climate/climate_prepped/mean_prec_third_quart.tif",
                   "data/variables/climate/climate_prepped/mean_prec_fourth_quart.tif",
                   "data/variables/climate/climate_prepped/mean_diff_first_quart.tif",
                   "data/variables/climate/climate_prepped/mean_diff_second_quart.tif",
                   "data/variables/climate/climate_prepped/mean_diff_third_quart.tif",
                   "data/variables/climate/climate_prepped/mean_diff_fourth_quart.tif",
                   "data/variables/climate/climate_prepped/mean_tmin_first_quart.tif",
                   "data/variables/climate/climate_prepped/mean_tmin_second_quart.tif",
                   "data/variables/climate/climate_prepped/mean_tmin_third_quart.tif",
                   "data/variables/climate/climate_prepped/mean_tmin_fourth_quart.tif",
                   "data/variables/climate/climate_prepped/mean_tmax_first_quart.tif",
                   "data/variables/climate/climate_prepped/mean_tmax_second_quart.tif",
                   "data/variables/climate/climate_prepped/mean_tmax_third_quart.tif",
                   "data/variables/climate/climate_prepped/mean_tmax_fourth_quart.tif"
                   )



multi_layer_climate <- terra::rast(climate_files)

tictoc::tic()
climate_res <- terra::extract(multi_layer_climate, points_3035, 
                              method = "simple", xy = T)
tictoc::toc()

head(climate_res)

write.csv(climate_res, "variable_manipulation/variable_outputs/climate_output.csv")

# # do with a small one first
# small_ext <- terra::ext(2000000, 2005000, 1000000, 1005000) # swap to base raster later
# blank_small <- rast(crs=crs, extent=small_ext, res = 1000)
# values(blank_small) <- seq(1:9)
# points_small <- terra::as.points(blank_small) 
# 
# plot(blank_small)
# plot(points_small, add = T)
# # this appears to be the points in the middle of the raster. 
# 
# small_res <- terra::extract(blank_small, points_small, method = "simple", xy = T)
