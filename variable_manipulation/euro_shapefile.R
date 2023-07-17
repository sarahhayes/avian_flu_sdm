# make a map of Europe polygon for use in visualising our data/results.
# 09/06/2023
# Using terra package to avoid issues with rgdal which is no longer being maintained

rm(list = ls())

library(terra)

zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                  layer = "CNTR_RG_03M_2020_4326")
plot(zipmap)

# change projection and crop 
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) 

# using the blank raster doesn't seem to work with the vector so do step by step 
euro_map_crop <- terra::project(x = zipmap, y = crs) 
euro_map_crop2 <- terra::crop(euro_map_crop, euro_ext) 
euro_map_crop
euro_map_crop2

plot(euro_map_crop)
plot(euro_map_crop2)

#writeVector(euro_map_crop2, "output/euro_map.shp", overwrite=TRUE)


