## 21/09/2023
## Trying to get some humidity data for use within the SDMs.
## A couple of different sources considered

# install.packages("ncdf4")

rm(list = ls())

library(terra)
library(raster)

## Looking at data from Copernicus project.

files_list <- list.files("data/variables/climate/humidity_raw", pattern="\\.nc$")
files_list_full <- paste("data/variables/climate/humidity_raw/",files_list, sep = "")

hum_stack <- terra::rast(files_list_full)
hum_stack

plot(hum_stack[[1]])
hum_stack[[60]]

# Split into the different quarters
first_quart_stack <- hum_stack[[1:90]]
first_quart_stack
second_quart_stack <- hum_stack[[91:181]]
second_quart_stack
third_quart_stack <- hum_stack[[182:273]]
third_quart_stack
fourth_quart_stack <- hum_stack[[274:365]]
fourth_quart_stack

# find the mean of the humidity in each quarter
mean_humidity_q1 <- mean(first_quart_stack)
mean_humidity_q1
plot(mean_humidity_q1)

mean_humidity_q2 <- mean(second_quart_stack)
mean_humidity_q2
plot(mean_humidity_q2)

mean_humidity_q3 <- mean(third_quart_stack)
mean_humidity_q3
plot(mean_humidity_q3)

mean_humidity_q4 <- mean(fourth_quart_stack)
mean_humidity_q4
plot(mean_humidity_q4)

# Still global extent , WGS 84 projection and 0.1 degree resolution so need to change to the resolution we want. 

# Transforming the data
blank_3035 <- terra::rast("output/euro_rast.tif")

mean_hum_q1_prj <- terra::project(mean_humidity_q1, blank_3035, method = "near")
# using nearest neghbour as we are going from larger resolution to smaller. 
mean_hum_q1_prj
plot(mean_hum_q1_prj)

mean_hum_q2_prj <- terra::project(mean_humidity_q2, blank_3035, method = "near")
mean_hum_q2_prj
plot(mean_hum_q2_prj)

mean_hum_q3_prj <- terra::project(mean_humidity_q3, blank_3035, method = "near")
mean_hum_q3_prj
plot(mean_hum_q3_prj)

mean_hum_q4_prj <- terra::project(mean_humidity_q4, blank_3035, method = "near")
mean_hum_q4_prj
plot(mean_hum_q4_prj)

#terra::writeRaster(mean_hum_q1_prj, "variable_manipulation/variable_outputs/mean_relative_humidity_q1.tif")
#terra::writeRaster(mean_hum_q2_prj, "variable_manipulation/variable_outputs/mean_relative_humidity_q2.tif")
#terra::writeRaster(mean_hum_q3_prj, "variable_manipulation/variable_outputs/mean_relative_humidity_q3.tif")
#terra::writeRaster(mean_hum_q4_prj, "variable_manipulation/variable_outputs/mean_relative_humidity_q4.tif")



# This part is using data from the Met Office. But resolution is not great and we have missing data for some of our study area. 
# https://www.metoffice.gov.uk/hadobs/hadisdh/

humid_terra <- terra::rast("data/variables/climate/climate_raw/HadISDH.landq.4.5.1.2022f_FLATgridHOM5by5_anoms9120.nc")

humid_raster <- raster::brick("data/variables/climate/climate_raw/HadISDH.landq.4.5.1.2022f_FLATgridHOM5by5_anoms9120.nc")

humid_terra
humid_raster

# The terra version has only one layer, whilst the raster one has 600. 

humid_layer1 <- humid_raster[[1]]
humid_layer1
plot(humid_layer1)

humid_layer599 <- humid_raster[[599]]
humid_layer599
plot(humid_layer599)

humid_layer600 <- humid_raster[[600]]
humid_layer600
plot(humid_layer600)

# Looks like the data are monthly with the 12 most recent layers being the data from 2022. 

blank_3035 <- terra::rast("output/euro_rast.tif")
blank_3035
plot(blank_3035)

raster::writeRaster(humid_layer600, "data/variables/climate/climate_raw/humidity_dec_2022.tif")

humid_dec_2022 <- terra::rast("data/variables/climate/climate_raw/humidity_dec_2022.tif")

humid_12_2022 <- terra::project(x = humid_dec_2022, y = blank_3035, method = "near")
humid_12_2022
plot(humid_12_2022)
map_boundary_terra <- terra::vect("output/euro_map.shp")
plot(map_boundary_terra, add= T)

####
raster::writeRaster(humid_layer599, "data/variables/climate/climate_raw/humidity_nov_2022.tif")
humid_nov_2022 <- terra::rast("data/variables/climate/climate_raw/humidity_nov_2022.tif")

humid_11_2022 <- terra::project(x = humid_nov_2022, y = blank_3035, method = "near")
humid_11_2022
plot(humid_11_2022)
plot(map_boundary_terra, add= T)
