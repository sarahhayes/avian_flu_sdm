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


#blank_3035 <- terra::rast("output/euro_rast.tif")
blank_3035 <- terra::rast("output/euro_rast_10k.tif")
blank_3035
#terra::xyFromCell(blank_3035, 1) # coordinates of the centre of the first cell

plot(blank_3035)

## Use a function to make the 12 diff rasters with the difference between max and min temps 

diff_fun <- function(min_rast, max_rast, blank_raster, name_to_save){
  month_temps <- c(min_rast, max_rast)
  multi_layer_temps <- terra::rast(month_temps)
  diff_temps <- terra::diff(multi_layer_temps)
#  diff_temps_prj <- terra::project(x = diff_temps, y = blank_raster, method = "near") # used when producing 1k res as making smaller
  diff_temps_prj <- terra::project(x = diff_temps, y = blank_raster, method = "bilinear") # used for 10k res as making larger
  terra::writeRaster(diff_temps_prj, paste("data/variables/climate/climate_prepped/",
                                           name_to_save, ".tif", sep = ""), overwrite = T)
  
}

jan_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-01.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-01.tif",
                      blank_raster =  blank_3035,
                   #   "jan_diffs")
                      "jan_diffs_10kres")

# Now generate the diff rasters for all the months 
feb_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-02.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-02.tif",
                      blank_raster =  blank_3035,
                    #  "feb_diffs")
                      "feb_diffs_10kres")

mar_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-03.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-03.tif",
                      blank_raster =  blank_3035,
                    #  "mar_diffs")
                      "mar_diffs_10kres")

apr_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-04.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-04.tif",
                      blank_raster =  blank_3035,
                    #  "apr_diffs")
                      "apr_diffs_10kres")

may_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-05.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-05.tif",
                      blank_raster =  blank_3035,
                     # "may_diffs")
                      "may_diffs_10kres")

jun_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-06.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-06.tif",
                      blank_raster =  blank_3035,
                     # "jun_diffs")
                      "jun_diffs_10kres")

jul_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-07.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-07.tif",
                      blank_raster =  blank_3035,
                     # "jul_diffs")
                      "jul_diffs_10kres")

aug_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-08.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-08.tif",
                      blank_raster =  blank_3035,
                     # "aug_diffs")
                      "aug_diffs_10kres")

sep_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-09.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-09.tif",
                      blank_raster =  blank_3035,
                    #  "sep_diffs")
                      "sep_diffs_10kres")

oct_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-10.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-10.tif",
                      blank_raster =  blank_3035,
                    #  "oct_diffs")
                      "oct_diffs_10kres")

nov_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-11.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-11.tif",
                      blank_raster =  blank_3035,
                    #  "nov_diffs")
                      "nov_diffs_10kres")

dec_diffs <- diff_fun("data/variables/climate/climate_raw/wc2.1_2.5m_tmin_2018-12.tif",
                      "data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-12.tif",
                      blank_raster =  blank_3035,
                    #  "dec_diffs")
                      "dec_diffs_10kres")

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


jan_min <- terra::project(jan_min, y = blank_3035, method = "bilinear")
jan_mean <- jan_min + (jan_diffs/2)
plot(jan_mean)
plot(jan_min)
plot(jan_diffs)

feb_min <- terra::project(feb_min, y = blank_3035, method = "bilinear")
feb_mean <- feb_min + (feb_diffs/2)

mar_min <- terra::project(mar_min, y = blank_3035, method = "bilinear")
mar_mean <- mar_min + (mar_diffs/2)

apr_min <- terra::project(apr_min, y = blank_3035, method = "bilinear")
apr_mean <- apr_min + (apr_diffs/2)

may_min <- terra::project(may_min, y = blank_3035, method = "bilinear")
may_mean <- may_min + (may_diffs/2)

jun_min <- terra::project(jun_min, y = blank_3035, method = "bilinear")
jun_mean <- jun_min + (jun_diffs/2)

jul_min <- terra::project(jul_min, y = blank_3035, method = "bilinear")
jul_mean <- jul_min + (jul_diffs/2)

aug_min <- terra::project(aug_min, y = blank_3035, method = "bilinear")
aug_mean <- aug_min + (aug_diffs/2)

sep_min <- terra::project(sep_min, y = blank_3035, method = "bilinear")
sep_mean <- sep_min + (sep_diffs/2)

oct_min <- terra::project(oct_min, y = blank_3035, method = "bilinear")
oct_mean <- oct_min + (oct_diffs/2)

nov_min <- terra::project(nov_min, y = blank_3035, method = "bilinear")
nov_mean <- nov_min + (nov_diffs/2)

dec_min <- terra::project(dec_min, y = blank_3035, method = "bilinear")
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

# terra::writeRaster(mean_mean_first_quart, "data/variables/climate/climate_prepped/mean_temp_q1.tif")
# terra::writeRaster(mean_mean_second_quart, "data/variables/climate/climate_prepped/mean_temp_q2.tif")
# terra::writeRaster(mean_mean_third_quart, "data/variables/climate/climate_prepped/mean_temp_q3.tif")
# terra::writeRaster(mean_mean_fourth_quart, "data/variables/climate/climate_prepped/mean_temp_q4.tif")

# terra::writeRaster(mean_mean_first_quart, "data/variables/climate/climate_prepped/mean_temp_q1_10kres.tif")
# terra::writeRaster(mean_mean_second_quart, "data/variables/climate/climate_prepped/mean_temp_q2_10kres.tif")
# terra::writeRaster(mean_mean_third_quart, "data/variables/climate/climate_prepped/mean_temp_q3_10kres.tif")
# terra::writeRaster(mean_mean_fourth_quart, "data/variables/climate/climate_prepped/mean_temp_q4_10kres.tif")



## labels in the means are still as for min. Plot the march mins next to the mean just to do a visual check
## that isn't the same. 

plot(files_mean_first_quart)
plot(mar_min, add = T)

files_mean_first_quart

q1_min_means <- terra::app(files_mean_first_quart, min)
q1_min_means # seems to have selected the minimums of the files
q1_max_means <- terra::app(files_mean_first_quart, max)
q1_max_min  <- c(q1_min_means, q1_max_means)
q1_max_min
q1_mean_range <- terra::diff(q1_max_min)
q1_mean_range
plot(q1_mean_range)

# run some checks

test_dat <- seq(1,50000, 500)

q1_test_min <- q1_min_means[test_dat]
q1_test_max <- q1_max_means[test_dat]
q1_test_range <- q1_mean_range[test_dat]

test_df <- cbind(q1_test_min, q1_test_max, q1_test_range)
colnames(test_df) <- c("min", "max", "diff")

test_df$check_range <- test_df$max  - test_df$min

nrow(test_df[which(test_df$check_range != test_df$diff),])

## Seems to be OK so run for the other quarters. 
q2_min_means <- terra::app(files_mean_second_quart, min)
q2_max_means <- terra::app(files_mean_second_quart, max)
q2_max_min  <- c(q2_min_means, q2_max_means)
q2_mean_range <- terra::diff(q2_max_min)
plot(q2_mean_range)

## Q3
q3_min_means <- terra::app(files_mean_third_quart, min)
files_mean_third_quart
q3_min_means
q3_max_means <- terra::app(files_mean_third_quart, max)
q3_max_min  <- c(q3_min_means, q3_max_means)
q3_mean_range <- terra::diff(q3_max_min)
plot(q3_mean_range)


## Q4

q4_min_means <- terra::app(files_mean_fourth_quart, min)
q4_max_means <- terra::app(files_mean_fourth_quart, max)
files_mean_fourth_quart
q4_max_means
q4_max_min  <- c(q4_min_means, q4_max_means)
q4_mean_range <- terra::diff(q4_max_min)
plot(q4_mean_range)


# terra::writeRaster(q1_mean_range, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q1.tif")
# terra::writeRaster(q2_mean_range, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q2.tif")
# terra::writeRaster(q3_mean_range, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q3.tif")
# terra::writeRaster(q4_mean_range, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q4.tif")

terra::writeRaster(q1_mean_range, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q1_10kres.tif")
terra::writeRaster(q2_mean_range, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q2_10kres.tif")
terra::writeRaster(q3_mean_range, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q3_10kres.tif")
terra::writeRaster(q4_mean_range, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q4_10kres.tif")


