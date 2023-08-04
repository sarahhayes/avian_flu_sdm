# 02/08/2023
# Creating the raster which will be used to 'mask' all the other layers
# this will be of the extent that we have decided
# but will only be the land

rm(list = ls())

library(tidyverse)
library(terra)

zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
plot(zipmap)

# change projection and crop 
crs <- "epsg:3035"
euro_ext <- terra::ext(2600000, 7000000, 1500000, 6400000) 
large_ext <- terra::ext(2000000, 9000000, 1000000, 7000000)

# using a blank raster doesn't seem to work with the vector so do step by step 
euro_map_crop <- terra::project(x = zipmap, y = crs) 
euro_map_crop_prj <- terra::crop(euro_map_crop, euro_ext) 
euro_map_crop
euro_map_crop_prj

plot(euro_map_crop)
plot(euro_map_crop_prj)

# writeVector(euro_map_crop_prj, "output/euro_map.shp", overwrite=TRUE)

# create the blank raster
blank_3035 <- terra::rast(crs=crs, extent=euro_ext, res = 1000)
blank_3035

eurorast <- rasterize(euro_map_crop_prj, blank_3035)
plot(eurorast)

# writeRaster(eurorast, "output/euro_rast.tif", overwrite = TRUE)

## Also create a 10km raster in case want to look at it
blank_10k <- terra::rast(crs=crs, extent=euro_ext, res = 10000)

eurorast_10k <- rasterize(euro_map_crop_prj, blank_10k)
plot(eurorast_10k)
# writeRaster(eurorast_10k, "output/euro_rast_10k.tif", overwrite = T)

# larger shapefile for use in distance to couast
large_ext <- terra::ext(1800000, 9000000, 1000000, 7200000)

# using a blank raster doesn't seem to work with the vector so do step by step 
large_map_crop <- terra::project(x = zipmap, y = crs) 
large_map_crop_prj <- terra::crop(large_map_crop, large_ext) 
plot(large_map_crop_prj)
plot(euro_map_crop_prj, col = "red", add = T, axes = F)

# writeVector(large_map_crop_prj, "output/large_map.shp", overwrite = T)
