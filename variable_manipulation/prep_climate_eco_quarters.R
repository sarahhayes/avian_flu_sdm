## 07/12/2024

## The purpose of this script is to take the raw data for the climate variables and transform 
## into the data needed for the model. 
## This script is using the ecological quarters instead of the seasonal quarters

### Data are worldclim historical monthly data
### https://www.worldclim.org/data/monthlywth.html
### https://www.worldclim.org/data/downscaling.html
### https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.06/
### https://www.nature.com/articles/s41597-020-0453-3 to check understand variables correctly 

### TMN is calculated as TMP − 0.5*DTR, and TMX as TMP + 0.5*DTR.

rm(list = ls())

library(tidyverse)
library(terra)

# we need a reference raster for later. 
blank_3035 <- terra::rast("output/euro_rast_10k.tif") # for 10km res
blank_3035
terra::xyFromCell(blank_3035, 1) # coordinates of the centre of the first cell

plot(blank_3035)

## Read in precipitation data 

jan_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-01.tif")
feb_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-02.tif")
mar_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-03.tif")
apr_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-04.tif")
may_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-05.tif")
jun_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-06.tif")
jul_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-07.tif")
aug_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-08.tif")
sep_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-09.tif")
oct_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-10.tif")
nov_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-11.tif")
dec_prec <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_prec_2018-12.tif")

### For ecological variables, have different quarters. 
### Use number of days in each month to make the means.
### Weights below

q1_wts <- c(1,31,31,28) # c(nov, dec, jan, feb)
q2_wts <- c(31,30,31,6) # (march, april, may, june)
q3_wts <- c(24,31,9) # (june, july, aug)
q4_wts <- c(22,30,31,29) # (aug, sep, oct, nov)


q1_prec_files <- c(nov_prec, dec_prec, jan_prec, feb_prec)
q2_prec_files <- c(mar_prec, apr_prec, may_prec, jun_prec)
q3_prec_files <- c(jun_prec, jul_prec, aug_prec)
q4_prec_files <- c(aug_prec, sep_prec, oct_prec, nov_prec)

q1_prec_wt_mean <- terra::weighted.mean(q1_prec_files, q1_wts)
q1_prec_wt_mean_prj <- terra::project(x = q1_prec_wt_mean, y = blank_3035, method = "bilinear") # re-project using the extent of blank raster and bilinear interpolation. When going to 10km res
q1_prec_wt_mean_prj
plot(q1_prec_wt_mean_prj)
#terra::writeRaster(q1_prec_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_prec_first_quart_eco_rasts.tif", overwrite = T)

# Q2
q2_prec_wt_mean <- terra::weighted.mean(q2_prec_files, q2_wts)
q2_prec_wt_mean_prj <- terra::project(x = q2_prec_wt_mean, y = blank_3035, method = "bilinear") 
q2_prec_wt_mean_prj
plot(q2_prec_wt_mean_prj)
#terra::writeRaster(q2_prec_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_prec_second_quart_eco_rasts.tif", overwrite = T)

# Q3
q3_prec_wt_mean <- terra::weighted.mean(q3_prec_files, q3_wts)
q3_prec_wt_mean_prj <- terra::project(x = q3_prec_wt_mean, y = blank_3035, method = "bilinear") 
q3_prec_wt_mean_prj
plot(q3_prec_wt_mean_prj)
#terra::writeRaster(q3_prec_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_prec_third_quart_eco_rasts.tif", overwrite = T)

# Q4
q4_prec_wt_mean <- terra::weighted.mean(q4_prec_files, q4_wts)
q4_prec_wt_mean_prj <- terra::project(x = q4_prec_wt_mean, y = blank_3035, method = "bilinear") 
q4_prec_wt_mean_prj
plot(q4_prec_wt_mean_prj)
#terra::writeRaster(q4_prec_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_prec_fourth_quart_eco_rasts.tif", overwrite = T)

####################################################################

### TEMPERATURE
## Next step is to calculate the temperature variability. 
## This will be done by extracting min from max for each month and then combining to find the average

# minimums
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

# maximums
jan_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-01.tif")
feb_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-02.tif")
mar_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-03.tif")
apr_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-04.tif")
may_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-05.tif")
jun_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-06.tif")
jul_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-07.tif")
aug_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-08.tif")
sep_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-09.tif")
oct_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-10.tif")
nov_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-11.tif")
dec_max <- terra::rast("data/variables/climate/climate_raw/wc2.1_2.5m_tmax_2018-12.tif")

jan_diff <- abs(jan_max - jan_min)
feb_diff <- abs(feb_max - feb_min)
mar_diff <- abs(mar_max - mar_min)
apr_diff <- abs(apr_max - apr_min)
may_diff <- abs(may_max - may_min)
jun_diff <- abs(jun_max - jun_min)
jul_diff <- abs(jul_max - jul_min)
aug_diff <- abs(aug_max - aug_min)
sep_diff <- abs(sep_max - sep_min)
oct_diff <- abs(oct_max - oct_min)
nov_diff <- abs(nov_max - nov_min)
dec_diff <- abs(dec_max - dec_min)


q1_diff_files <- c(nov_diff, dec_diff, jan_diff, feb_diff)
q2_diff_files <- c(mar_diff, apr_diff, may_diff, jun_diff)
q3_diff_files <- c(jun_diff, jul_diff, aug_diff)
q4_diff_files <- c(aug_diff, sep_diff, oct_diff, nov_diff)

#Q1
# Note - these temperatures are 2m temperature
q1_diff_wt_mean <- terra::weighted.mean(q1_diff_files, q1_wts)
q1_diff_wt_mean_prj <- terra::project(x = q1_diff_wt_mean, y = blank_3035, method = "bilinear") 
q1_diff_wt_mean_prj
plot(q1_diff_wt_mean_prj)
#terra::writeRaster(q1_diff_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_diff_first_quart_eco_rasts.tif", overwrite = T)


#Q2
q2_diff_wt_mean <- terra::weighted.mean(q2_diff_files, q2_wts)
q2_diff_wt_mean_prj <- terra::project(x = q2_diff_wt_mean, y = blank_3035, method = "bilinear") 
q2_diff_wt_mean_prj
plot(q2_diff_wt_mean_prj)
#terra::writeRaster(q2_diff_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_diff_second_quart_eco_rasts.tif", overwrite = T)


#Q3
q3_diff_wt_mean <- terra::weighted.mean(q3_diff_files, q3_wts)
q3_diff_wt_mean_prj <- terra::project(x = q3_diff_wt_mean, y = blank_3035, method = "bilinear") 
q3_diff_wt_mean_prj
plot(q3_diff_wt_mean_prj)
#terra::writeRaster(q3_diff_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_diff_third_quart_eco_rasts.tif", overwrite = T)

#Q4
q4_diff_wt_mean <- terra::weighted.mean(q4_diff_files, q4_wts)
q4_diff_wt_mean_prj <- terra::project(x = q4_diff_wt_mean, y = blank_3035, method = "bilinear") 
q4_diff_wt_mean_prj
plot(q4_diff_wt_mean_prj)
#terra::writeRaster(q4_diff_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_diff_fourth_quart_eco_rasts.tif", overwrite = T)


###############################################################################
# Mean temp

## Estimate the mean temperature using the formula from the paper
## To get the mean temperature, we can either do:
## Max temp - diurnal range/2
## Or
## Min temp + diurnal range/2 

jan_mean <- jan_min + (jan_diff/2)
feb_mean <- feb_min + (feb_diff/2)
mar_mean <- mar_min + (mar_diff/2)
apr_mean <- apr_min + (apr_diff/2)
may_mean <- may_min + (may_diff/2)
jun_mean <- jun_min + (jun_diff/2)
jul_mean <- jul_min + (jul_diff/2)
aug_mean <- aug_min + (aug_diff/2)
sep_mean <- sep_min + (sep_diff/2)
oct_mean <- oct_min + (oct_diff/2)
nov_mean <- nov_min + (nov_diff/2)
dec_mean <- dec_min + (dec_diff/2)

q1_mean_files <- c(nov_mean, dec_mean, jan_mean, feb_mean)
q2_mean_files <- c(mar_mean, apr_mean, may_mean, jun_mean)
q3_mean_files <- c(jun_mean, jul_mean, aug_mean)
q4_mean_files <- c(aug_mean, sep_mean, oct_mean, nov_mean)

#Q1
q1_mean_wt_mean <- terra::weighted.mean(q1_mean_files, q1_wts)
q1_mean_wt_mean_prj <- terra::project(x = q1_mean_wt_mean, y = blank_3035, method = "bilinear") 
q1_mean_wt_mean_prj
plot(q1_mean_wt_mean_prj)
#terra::writeRaster(q1_mean_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_mean_first_quart_eco_rasts.tif", overwrite = T)

#Q2
q2_mean_wt_mean <- terra::weighted.mean(q2_mean_files, q2_wts)
q2_mean_wt_mean_prj <- terra::project(x = q2_mean_wt_mean, y = blank_3035, method = "bilinear") 
q2_mean_wt_mean_prj
plot(q2_mean_wt_mean_prj)
#terra::writeRaster(q2_mean_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_mean_second_quart_eco_rasts.tif", overwrite = T)

#Q3
q3_mean_wt_mean <- terra::weighted.mean(q3_mean_files, q3_wts)
q3_mean_wt_mean_prj <- terra::project(x = q3_mean_wt_mean, y = blank_3035, method = "bilinear") 
q3_mean_wt_mean_prj
plot(q3_mean_wt_mean_prj)
#terra::writeRaster(q3_mean_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_mean_third_quart_eco_rasts.tif", overwrite = T)

#Q4
q4_mean_wt_mean <- terra::weighted.mean(q4_mean_files, q4_wts)
q4_mean_wt_mean_prj <- terra::project(x = q4_mean_wt_mean, y = blank_3035, method = "bilinear") 
q4_mean_wt_mean_prj
plot(q4_mean_wt_mean_prj)
#terra::writeRaster(q4_mean_wt_mean_prj, 
#                   "data/variables/climate/climate_prepped/mean_mean_fourth_quart_eco_rasts.tif", overwrite = T)


###########################################################################
## Quarterly variation in mean

# for this we are looking at any month that's mostly in the season,
# i.e. difference between the maximum mean of any month with at least 50% of its days in the season and
#  the minimum mean of any month with at least 50% of its days in the season

# Reminder of the cutoffs
#q1_wts <- c(1,31,31,28) # c(nov, dec, jan, feb)
#q2_wts <- c(31,30,31,6) # (march, april, may, june)
#q3_wts <- c(24,31,9) # (june, july, aug)
#q4_wts <- c(22,30,31,29) # (aug, sep, oct, nov)

q1_quart_var_files <- c(dec_mean, jan_mean, feb_mean)
q2_quart_var_files <- c(mar_mean, apr_mean, may_mean)
q3_quart_var_files <- c(jun_mean, jul_mean)
q4_quart_var_files <- c(aug_mean, sep_mean, oct_mean, nov_mean)

#Q1
q1_min_means <- terra::app(q1_quart_var_files, min)
q1_min_means # seems to have selected the minimums of the files
q1_max_means <- terra::app(q1_quart_var_files, max)
q1_max_min  <- c(q1_min_means, q1_max_means)
q1_max_min
q1_mean_range <- terra::diff(q1_max_min)
q1_mean_range
plot(q1_mean_range)
q1_mean_range_prj <- terra::project(x = q1_mean_range, y = blank_3035, method = "bilinear") 
plot(q1_mean_range_prj)
# terra::writeRaster(q1_mean_range_prj, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q1_eco_rasts.tif")

#Q2
q2_min_means <- terra::app(q2_quart_var_files, min)
q2_min_means 
q2_max_means <- terra::app(q2_quart_var_files, max)
q2_max_min  <- c(q2_min_means, q2_max_means)
q2_max_min
q2_mean_range <- terra::diff(q2_max_min)
q2_mean_range
plot(q2_mean_range)
q2_mean_range_prj <- terra::project(x = q2_mean_range, y = blank_3035, method = "bilinear") 
plot(q2_mean_range_prj)
#terra::writeRaster(q2_mean_range_prj, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q2_eco_rasts.tif")

#Q3
q3_min_means <- terra::app(q3_quart_var_files, min)
q3_max_means <- terra::app(q3_quart_var_files, max)
q3_max_min  <- c(q3_min_means, q3_max_means)
q3_mean_range <- terra::diff(q3_max_min)
plot(q3_mean_range)
q3_mean_range_prj <- terra::project(x = q3_mean_range, y = blank_3035, method = "bilinear") 
plot(q3_mean_range_prj)
#terra::writeRaster(q3_mean_range_prj, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q3_eco_rasts.tif")

#Q4
q4_min_means <- terra::app(q4_quart_var_files, min)
q4_max_means <- terra::app(q4_quart_var_files, max)
q4_max_min  <- c(q4_min_means, q4_max_means)
q4_mean_range <- terra::diff(q4_max_min)
plot(q4_mean_range)
q4_mean_range_prj <- terra::project(x = q4_mean_range, y = blank_3035, method = "bilinear") 
plot(q4_mean_range_prj)
#terra::writeRaster(q4_mean_range_prj, "data/variables/climate/climate_prepped/variation_in_quarterly_mean_temp_q4_eco_rasts.tif")
