## 24/05/2023

## The purpose of this script is to take the raw data for the climate variables and transform 
## into the data needed for the model. 
## Seasons will be calculated quarterly, starting with January 
## Jan- Mar, Apr - Jun, Jul - Sept, Oct - Dec

rm(list = ls())

library(tidyverse)
library(terra)


files_list <- c("Data/IA_cdl_2015.tif", "Data/IA_cdl_2016.tif")

# make lists of files so can import simultaneously

files_tmin_first_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-01.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-02.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-03.tif")



multi_layer_tmin_first_quart <- terra::rast(files_tmin_first_quart)

# look at the details of the raster
multi_layer_tmin_first_quart

# calculate the mean over the 3 month period
mean_tmin_first_quart <- terra::app(multi_layer_tmin_first_quart, mean)
mean_tmin_first_quart

# now change projection and crop 
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later

mean_tmin_first_quart_crs <- terra::project(x = mean_tmin_first_quart, y = crs, method = "near") 
mean_tmin_first_quart_crs_crop <- crop(x = mean_tmin_first_quart_crs, y = euro_ext )

plot(mean_tmin_first_quart_crs_crop)

terra::writeRaster(mean_tmin_first_quart_crs_crop, 
                   "data/variables/climate/climate_prepped/mean_tmin_first_quart.tif",
                   overwrite = T) 

check_rast <- terra::rast("data/variables/climate/climate_prepped/mean_tmin_first_quart.tif")

## check that match pre- and post-save
mean_tmin_first_quart_crs_crop
check_rast
## make a function to do this for efficiency

climate_mean_fun <- function(data_list, set_crs, set_extent, name_to_save){
  multi_layer_file <- terra::rast(data_list)
  mean_multi_layer <- terra::app(multi_layer_file, mean)
  mean_multi_layer_crs <- terra::project(x = mean_multi_layer, y = set_crs, method = "near") 
  mean_multi_layer_crs_crop <- crop(x = mean_multi_layer_crs, y = set_extent )
  terra::writeRaster(mean_multi_layer_crs_crop, paste("data/variables/climate/climate_prepped/",
                                                      name_to_save, ".tif", sep = ""), overwrite = T) 
 
}


files_tmin_second_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-04.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-05.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-06.tif")

files_tmin_third_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-07.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-08.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-09.tif")

files_tmin_fourth_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-10.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-11.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-12.tif")


mean_tmin_first_quart <- climate_mean_fun(files_tmin_first_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_tmin_first_quart")
mean_tmin_second_quart <- climate_mean_fun(files_tmin_second_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_tmin_second_quart")
mean_tmin_third_quart <- climate_mean_fun(files_tmin_third_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_tmin_third_quart")
mean_tmin_fourth_quart <- climate_mean_fun(files_tmin_fourth_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_tmin_fourth_quart")

# min_first <- terra::rast("data/variables/climate/climate_prepped/mean_tmin_first_quart.tif")
# plot(min_first)
# min_second <- terra::rast("data/variables/climate/climate_prepped/mean_tmin_second_quart.tif")
# plot(min_second)
# min_third <- terra::rast("data/variables/climate/climate_prepped/mean_tmin_third_quart.tif")
# plot(min_third)
# min_fourth <- terra::rast("data/variables/climate/climate_prepped/mean_tmin_fourth_quart.tif")
# plot(min_fourth)


## Repeat for maximum temperature 

files_tmax_first_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-01.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-02.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-03.tif")

files_tmax_second_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-04.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-05.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-06.tif")

files_tmax_third_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-07.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-08.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-09.tif")

files_tmax_fourth_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-10.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-11.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-12.tif")


mean_tmax_first_quart <- climate_mean_fun(files_tmax_first_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_tmax_first_quart")
mean_tmax_second_quart <- climate_mean_fun(files_tmax_second_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_tmax_second_quart")
mean_tmax_third_quart <- climate_mean_fun(files_tmax_third_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_tmax_third_quart")
mean_tmax_fourth_quart <- climate_mean_fun(files_tmax_fourth_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_tmax_fourth_quart")


## And for precipitation 

files_prec_first_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-01.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-02.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-03.tif")

files_prec_second_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-04.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-05.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-06.tif")

files_prec_third_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-07.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-08.tif",
                            "data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-09.tif")

files_prec_fourth_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-10.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-11.tif",
                             "data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-12.tif")


mean_prec_first_quart <- climate_mean_fun(files_prec_first_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_prec_first_quart")
mean_prec_second_quart <- climate_mean_fun(files_prec_second_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_prec_second_quart")
mean_prec_third_quart <- climate_mean_fun(files_prec_third_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_prec_third_quart")
mean_prec_fourth_quart <- climate_mean_fun(files_prec_fourth_quart, set_crs = crs, set_extent = euro_ext, name_to_save = "mean_prec_fourth_quart")

