## Prepping human density data

rm(list = ls())
library(tidyverse)
library(terra)

dens <- terra::rast("data/variables/human_density/gpw_v4_population_density_rev11_2020_2pt5_min.tif")

dens
#plot(dens)
plot(app(dens, function(x) log(x+1))) # log transform for better viewing

# read in the reference raster
blank_3035 <- terra::rast("output/euro_rast_10k.tif")
blank_3035
plot(blank_3035)

dens_euro <- terra::project(x = dens, y = blank_3035, method = "bilinear")
dens_euro
#plot(dens_euro)
plot(app(dens_euro, function(x) log(x+1))) # log transform for better viewing

# # make a points object using the centre of each pixel from the blank raster
# points_3035 <- terra::as.points(blank_3035)
# points_3035
# #plot(points_3035)
# 
# tictoc::tic()
# dens_res <- terra::extract(dens_euro, points_3035, method = "simple", xy = T)
# tictoc::toc()
# 
# ## rename the columns
# dens_res <- rename(dens_res, "human_density_20" = "gpw_v4_population_density_rev11_2020_2pt5_min")
# 
# # save as a large combo csv rather than separate ones
# 
# write.csv(dens_res,
#           "variable_manipulation/variable_outputs/human_density_output.csv",
#           row.names = F)

# also save the rasters
terra::writeRaster(dens_euro, "variable_manipulation/variable_outputs/human_density.tif", overwrite = TRUE)
