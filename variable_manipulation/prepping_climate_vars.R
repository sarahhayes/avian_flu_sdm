## 24/05/2023

## The purpose of this script is to take the raw data for the climate variables and transform 
## into the data needed for the model. 
## Seasons will be calculated quarterly, starting with January 
## Jan- Mar, Apr - Jun, Jul - Sept, Oct - Dec

rm(list = ls())

library(tidyverse)
library(terra)

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

# Create a blank raster with appropriate projection and extent
blank_raster <- rast(crs=crs, extent=euro_ext)

mean_tmin_first_quart_crs <- terra::project(x = mean_tmin_first_quart, y = crs, method = "near") 
mean_tmin_first_quart_crs_crop <- crop(x = mean_tmin_first_quart_crs, y = euro_ext )

# crop and project simultaneously using the blank raster
mean_tmin_first_quart_raster_crop <- terra::project(x = mean_tmin_first_quart, y = blank_3035,
                                                    method = "near") 

plot(mean_tmin_first_quart_crs_crop)
plot(mean_tmin_first_quart_raster_crop)

mean_tmin_first_quart_crs_crop
mean_tmin_first_quart_raster_crop

#terra::writeRaster(mean_tmin_first_quart_crs_crop, 
#                   "data/variables/climate/climate_prepped/mean_tmin_first_quart.tif",
#                   overwrite = T) 

#check_rast <- terra::rast("data/variables/climate/climate_prepped/mean_tmin_first_quart.tif")

## check that match pre- and post-save
# mean_tmin_first_quart_crs_crop
# check_rast
## make a function to do this for efficiency

# climate_mean_fun <- function(data_list, set_crs, set_extent, name_to_save){
#   multi_layer_file <- terra::rast(data_list)
#   mean_multi_layer <- terra::app(multi_layer_file, mean)
#   mean_multi_layer_crs <- terra::project(x = mean_multi_layer, y = set_crs, method = "near")
#   mean_multi_layer_crs_crop <- crop(x = mean_multi_layer_crs, y = set_extent )
#   terra::writeRaster(mean_multi_layer_crs_crop, paste("data/variables/climate/climate_prepped/",
#                                                       name_to_save, ".tif", sep = ""), overwrite = T)
# 
# }

climate_mean_fun <- function(data_list, blank_raster, name_to_save){
  multi_layer_file <- terra::rast(data_list)
  mean_multi_layer <- terra::app(multi_layer_file, mean)
  mean_multi_layer_crs <- terra::project(x = mean_multi_layer, y = blank_raster, method = "near")
  terra::writeRaster(mean_multi_layer_crs, paste("data/variables/climate/climate_prepped/",
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


mean_tmin_first_quart <- climate_mean_fun(files_tmin_first_quart, blank_raster = blank_3035, name_to_save = "mean_tmin_first_quart")
mean_tmin_second_quart <- climate_mean_fun(files_tmin_second_quart, blank_raster = blank_3035, name_to_save = "mean_tmin_second_quart")
mean_tmin_third_quart <- climate_mean_fun(files_tmin_third_quart, blank_raster = blank_3035, name_to_save = "mean_tmin_third_quart")
mean_tmin_fourth_quart <- climate_mean_fun(files_tmin_fourth_quart, blank_raster = blank_3035, name_to_save = "mean_tmin_fourth_quart")


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


mean_tmax_first_quart <- climate_mean_fun(files_tmax_first_quart, blank_raster = blank_3035, name_to_save = "mean_tmax_first_quart")
mean_tmax_second_quart <- climate_mean_fun(files_tmax_second_quart, blank_raster = blank_3035, name_to_save = "mean_tmax_second_quart")
mean_tmax_third_quart <- climate_mean_fun(files_tmax_third_quart, blank_raster = blank_3035, name_to_save = "mean_tmax_third_quart")
mean_tmax_fourth_quart <- climate_mean_fun(files_tmax_fourth_quart, blank_raster = blank_3035, name_to_save = "mean_tmax_fourth_quart")


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


mean_prec_first_quart <- climate_mean_fun(files_prec_first_quart, blank_raster = blank_3035, name_to_save = "mean_prec_first_quart")
mean_prec_second_quart <- climate_mean_fun(files_prec_second_quart, blank_raster = blank_3035, name_to_save = "mean_prec_second_quart")
mean_prec_third_quart <- climate_mean_fun(files_prec_third_quart, blank_raster = blank_3035, name_to_save = "mean_prec_third_quart")
mean_prec_fourth_quart <- climate_mean_fun(files_prec_fourth_quart, blank_raster = blank_3035, name_to_save = "mean_prec_fourth_quart")

## Next step is to calculate the temperature variability. 
## This will be done by extracting min from max for each month and then combining to find the average

month_1_temps <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-01.tif", 
                   "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-01.tif") 
  
multi_layer_temps_1 <- terra::rast(month_1_temps)
multi_layer_temps_1

# get the difference between the max and min for each cell 
diff_temps_1 <- terra::diff(multi_layer_temps_1)
diff_temps_1
plot(diff_temps_1)

# check it is doing what you think
max_1 <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-01.tif")
min_1 <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-01.tif")

# Raster data can be thought of as a matrix, but in a SpatRaster it is more commonly treated 
# as a vector. Cells are numbered from the upper left cell to the upper right cell
# and then continuing on the left side of the next row, and so on until 
# the last cell at the lower-right side of the raster

# lots of the first rows and columns have NaN as no value so pick vals likely to be 'on land'
max_vals <- values(max_1)[9000000:9000050]
min_vals <- values(min_1)[9000000:9000050]

max_vals
min_vals

diff_vals <- abs(max_vals - min_vals)
diff_vals

values(diff_temps_1)[9000000:9000050]

diff_vals == values(diff_temps_1)[9000000:9000050]
## seems to be doing what we want

## now crop and project
diff_temps_1_prj <- terra::project(x = diff_temps_1, y = blank_raster, method = "near")
diff_temps_1_prj

## Pop this into a function to make the 12 diff rasters

diff_fun <- function(min_rast, max_rast, name_to_save){
  month_temps <- c(min_rast, max_rast)
  multi_layer_temps <- terra::rast(month_temps)
  diff_temps <- terra::diff(multi_layer_temps)
  diff_temps_prj <- terra::project(x = diff_temps, y = blank_raster, method = "near")
  terra::writeRaster(diff_temps_prj, paste("data/variables/climate/climate_prepped/",
                                                 name_to_save, ".tif", sep = ""), overwrite = T)
  
  }

jan_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-01.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-01.tif",
                      "jan_diffs")
# quick check
jan_diffs
diff_temps_1_prj
values(jan_diffs)[1070:1080]
values(diff_temps_1_prj)[1070:1080]

# Now generate the diff rasters for all the months 
feb_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-02.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-02.tif",
                      "feb_diffs")
                      
mar_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-03.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-03.tif",
                      "mar_diffs")

apr_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-04.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-04.tif",
                      "apr_diffs")

may_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-05.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-05.tif",
                      "may_diffs")

jun_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-06.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-06.tif",
                      "jun_diffs")

jul_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-07.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-07.tif",
                      "jul_diffs")

aug_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-08.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-08.tif",
                      "aug_diffs")

sep_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-09.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-09.tif",
                      "sep_diffs")

oct_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-10.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-10.tif",
                      "oct_diffs")

nov_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-11.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-11.tif",
                      "nov_diffs")

dec_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-12.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-12.tif",
                      "dec_diffs")


## and now get the means as previous
## And for precipitation 

files_diff_first_quart <- c("data/variables/climate/climate_prepped/jan_diffs.tif",
                            "data/variables/climate/climate_prepped/feb_diffs.tif",
                            "data/variables/climate/climate_prepped/mar_diffs.tif")

files_diff_second_quart <- c("data/variables/climate/climate_prepped/apr_diffs.tif",
                             "data/variables/climate/climate_prepped/may_diffs.tif",
                             "data/variables/climate/climate_prepped/jun_diffs.tif")

files_diff_third_quart <- c("data/variables/climate/climate_prepped/jul_diffs.tif",
                            "data/variables/climate/climate_prepped/aug_diffs.tif",
                            "data/variables/climate/climate_prepped/sep_diffs.tif")

files_diff_fourth_quart <- c("data/variables/climate/climate_prepped/oct_diffs.tif",
                             "data/variables/climate/climate_prepped/nov_diffs.tif",
                             "data/variables/climate/climate_prepped/dec_diffs.tif")


mean_diff_first_quart <- climate_mean_fun(files_diff_first_quart, blank_raster = blank_3035, name_to_save = "mean_diff_first_quart")
mean_diff_second_quart <- climate_mean_fun(files_diff_second_quart, blank_raster = blank_3035, name_to_save = "mean_diff_second_quart")
mean_diff_third_quart <- climate_mean_fun(files_diff_third_quart, blank_raster = blank_3035, name_to_save = "mean_diff_third_quart")
mean_diff_fourth_quart <- climate_mean_fun(files_diff_fourth_quart, blank_raster = blank_3035, name_to_save = "mean_diff_fourth_quart")

