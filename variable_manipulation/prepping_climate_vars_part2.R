## 12/10/2023

## The purpose of this script is to take the raw data for the climate variables and transform 
## into data needed for the model. This is a second script calculating diurnal temp range and mean temperature
## as per the information given in this paper https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/joc.3711?saml_referrer
## Seasons will be calculated quarterly, starting with January 
## Jan- Mar, Apr - Jun, Jul - Sept, Oct - Dec

rm(list = ls())

library(tidyverse)
library(terra)

# make lists of files so can import simultaneously
# 
# files_tmin_first_quart <- c("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-01.tif",
#                             "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-02.tif",
#                             "data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-03.tif")
# 
# multi_layer_tmin_first_quart <- terra::rast(files_tmin_first_quart)


blank_3035 <- terra::rast("output/euro_rast.tif")
blank_3035
#terra::xyFromCell(blank_3035, 1) # coordinates of the centre of the first cell

plot(blank_3035)


## Use a function to make the 12 diff rasters with the difference between max and min temps 

diff_fun <- function(min_rast, max_rast, blank_raster, name_to_save){
  month_temps <- c(min_rast, max_rast)
  multi_layer_temps <- terra::rast(month_temps)
  diff_temps <- terra::diff(multi_layer_temps)
  diff_temps_prj <- terra::project(x = diff_temps, y = blank_raster, method = "near")
  terra::writeRaster(diff_temps_prj, paste("data/variables/climate/climate_prepped/",
                                           name_to_save, ".tif", sep = ""), overwrite = T)
  
}

jan_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-01.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-01.tif",
                      blank_raster =  blank_3035,
                      "jan_diffs")

# Now generate the diff rasters for all the months 
feb_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-02.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-02.tif",
                      blank_raster =  blank_3035,
                      "feb_diffs")

mar_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-03.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-03.tif",
                      blank_raster =  blank_3035,
                      "mar_diffs")

apr_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-04.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-04.tif",
                      blank_raster =  blank_3035,
                      "apr_diffs")

may_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-05.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-05.tif",
                      blank_raster =  blank_3035,
                      "may_diffs")

jun_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-06.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-06.tif",
                      blank_raster =  blank_3035,
                      "jun_diffs")

jul_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-07.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-07.tif",
                      blank_raster =  blank_3035,
                      "jul_diffs")

aug_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-08.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-08.tif",
                      blank_raster =  blank_3035,
                      "aug_diffs")

sep_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-09.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-09.tif",
                      blank_raster =  blank_3035,
                      "sep_diffs")

oct_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-10.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-10.tif",
                      blank_raster =  blank_3035,
                      "oct_diffs")

nov_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-11.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-11.tif",
                      blank_raster =  blank_3035,
                      "nov_diffs")

dec_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-12.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-12.tif",
                      blank_raster =  blank_3035,
                      "dec_diffs")

## Using the formula from the paper, to get the mean temperature, we can either do:
## Max temp - diurnal range/2
## Or
## Min temp + diurnal range/2 

jan_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-01.tif")
feb_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-02.tif")
mar_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-03.tif")
apr_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-04.tif")
may_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-05.tif")
jun_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-06.tif")
jul_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-07.tif")
aug_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-08.tif")
sep_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-09.tif")
oct_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-10.tif")
nov_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-11.tif")
dec_min <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-12.tif")


jan_min <- terra::project(jan_min, y = blank_3035, method = "near")
jan_mean <- jan_min + (jan_diffs/2)
plot(jan_mean)
plot(jan_min)
plot(jan_diffs)

feb_min <- terra::project(feb_min, y = blank_3035, method = "near")
feb_mean <- feb_min + (feb_diffs/2)

mar_min <- terra::project(mar_min, y = blank_3035, method = "near")
mar_mean <- mar_min + (mar_diffs/2)

apr_min <- terra::project(apr_min, y = blank_3035, method = "near")
apr_mean <- apr_min + (apr_diffs/2)

may_min <- terra::project(may_min, y = blank_3035, method = "near")
may_mean <- may_min + (may_diffs/2)

jun_min <- terra::project(jun_min, y = blank_3035, method = "near")
jun_mean <- jun_min + (jun_diffs/2)

jul_min <- terra::project(jul_min, y = blank_3035, method = "near")
jul_mean <- jul_min + (jul_diffs/2)

aug_min <- terra::project(aug_min, y = blank_3035, method = "near")
aug_mean <- aug_min + (aug_diffs/2)

sep_min <- terra::project(sep_min, y = blank_3035, method = "near")
sep_mean <- sep_min + (sep_diffs/2)

oct_min <- terra::project(oct_min, y = blank_3035, method = "near")
oct_mean <- oct_min + (oct_diffs/2)

nov_min <- terra::project(nov_min, y = blank_3035, method = "near")
nov_mean <- nov_min + (nov_diffs/2)

dec_min <- terra::project(dec_min, y = blank_3035, method = "near")
dec_mean <- dec_min + (dec_diffs/2)

plot(aug_mean)

# now get the mean for the quarters
files_mean_first_quart <- c(jan_mean, feb_mean, mar_mean)
files_mean_second_quart <- c(apr_mean, may_mean, jun_mean)
files_mean_third_quart <- c(jul_mean, aug_mean, sep_mean)
files_mean_fourth_quart <- c(oct_mean, nov_mean, dec_mean)


# calculate the mean over the 3 month period

mean_mean_first_quart <- terra::app(files_mean_first_quart, mean)
mean_mean_first_quart
plot(mean_mean_first_quart)

mean_mean_second_quart <- terra::app(files_mean_second_quart, mean)
mean_mean_third_quart <- terra::app(files_mean_third_quart, mean)
mean_mean_fourth_quart <- terra::app(files_mean_fourth_quart, mean)

par(mfrow = c(2,2))
plot(mean_mean_first_quart)
plot(mean_mean_second_quart)
plot(mean_mean_third_quart)
plot(mean_mean_fourth_quart)

terra::writeRaster(mean_mean_first_quart, "data/variables/climate/climate_prepped/mean_temp_q1.tif")
terra::writeRaster(mean_mean_second_quart, "data/variables/climate/climate_prepped/mean_temp_q2.tif")
terra::writeRaster(mean_mean_third_quart, "data/variables/climate/climate_prepped/mean_temp_q3.tif")
terra::writeRaster(mean_mean_fourth_quart, "data/variables/climate/climate_prepped/mean_temp_q4.tif")
