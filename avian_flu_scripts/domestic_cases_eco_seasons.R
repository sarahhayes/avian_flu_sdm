## 22/03/2025
## Creating the avian flu database of domestic cases
## Plan is to overlay this on our risk map 

rm(list = ls())

library(readxl)
library(tidyverse)
library(terra)
library(RColorBrewer)

## bring in the different data sets and combine to one large data set. 
## processing has been done in domestic_cases.R when originally did this analysis for the quarters. 

hpai <- read.csv("data/flu_data/prepped_data/hpai_pos_domestic_birds.csv")

hpai <- hpai %>% dplyr::select(!X.1)

# makes sure the observation date is in date format
hpai$observation.date <- as.Date(hpai$observation.date, "%Y-%m-%d")

# Add a year and month and year_month to the data

hpai$year <- lubridate::year(hpai$observation.date)
hpai$month <- lubridate::month(hpai$observation.date)
hpai$month_year <-  format(as.Date(hpai$observation.date), "%Y-%m")
hpai$week <- format(hpai$observation.date, "%Y Week %W")
hpai$week_num <- lubridate::isoweek(hpai$observation.date)
hpai$day_yr <- lubridate::yday(hpai$observation.date)
hpai$day_mo <- lubridate::day(hpai$observation.date)

range(hpai$observation.date)
head(hpai)
# Now to plot the time series by serotype 

serotype_data <- as.data.frame(table(hpai$serotype_HN))

## To make sure we include all weeks in a plot, even those without any data, we need a dataframe with 
## all the weeks and months,
## # Now to plot the time series by serotype 

## start on 3rd Jan 2005 as that's a Monday
week_calendar <- data.frame(date=seq(as.Date("2005-01-03"), as.Date("2024-04-28"), by="day")) %>% 
  mutate(week_num=isoweek(date),year=year(date)) %>%
  group_by(year,week_num) %>% 
  summarise(weekdate=min(date)) 
week_calendar$week_of_study <- seq(1, nrow(week_calendar), by = 1)

# combine with the data
dates <- merge(hpai, week_calendar)

######

### Now want to create a raster at 10km resolution where we label a pixel as positive if it has had a domestic oputbreak
### Then will want to do with counts of numbers of outbreaks 
### and then perhaps look at it in the context of the poultry density? 
### For this latter option, might need to 'cut out' our wild bird risk map based on presence of poultry and overlay? 
### Will need some thought

euro_rast <- terra::rast("output/euro_rast_10k.tif")
plot(euro_rast)

# extract points
p <- cbind(hpai$X, hpai$Y)

plot(euro_rast)
points(p, pch = 18, cex = 0.5)

# rasterize points as a matrix
points_rast <- rasterize(p, euro_rast, value = 1) # this one just marks as present or not. Doesn't give an idea of number. 

plot(points_rast)
points_rast

# We want these split into the eco_seasons and also the periods over which we split the data. 

# data on cases in wild birds spanned August 2016 to Feb Feb 2024
# For A) models were trained on cases data from 10/8/2016 to 9/8/2020 
# before testing on an outbreak from 10/8/2020 to 9/8/2021 
# For B) models were trained on H5N1 cases data from 10/8/2021 to 28/2/2023
# before testing on H5N1 cases data from 1/3/2023 to 29/2/2024 

# So the full periods are

# period A - 10/08/16 to 09/08/2021
# preiod B - 10/08/2021 to 29/02/2024 

dom_a <- hpai %>% filter(observation.date <= "2021-08-09" & observation.date >= "2016-08-10")
dom_b <- hpai %>% filter(observation.date >= "2021-08-10" & observation.date <= "2024-02-29")

# Ecological seasons are:
# Q1 30th November to 28th (or 29th) Feb
# Q2 1st March to 6th June
# Q3 7th June to 9th August
# Q4 10th August to 29th November

dom_a_q1 <- dom_a %>% dplyr::filter(month == 11 & day_mo == 30 | month == 12 | month == 1 | 
                               month == 2 & day_mo <=29)
range(dom_a_q1$observation.date)

dom_a_q2 <- dom_a %>% dplyr::filter(month == 3 | month == 4 | month == 5 | 
                                      month == 6 & day_mo <= 6 )
range(dom_a_q2$observation.date)

dom_a_q3 <- dom_a %>% dplyr::filter(month ==6 & day_mo >=7 | month == 7 |
                                      month == 8 & day_mo <=9)
range(dom_q_q3$observation.date)

dom_a_q4 <- dom_a %>% dplyr::filter(month == 8 & day_mo >= 10 | month == 9 | month == 10 |
                                      month == 11 & day_mo <= 29)
range(dom_a_q4$observation.date)

# check the numbers add up.  
nrow(dom_a_q1) + nrow(dom_a_q2) + nrow(dom_a_q3) + nrow(dom_a_q4) == nrow(dom_a)

## Now repeat for period b

dom_b_q1 <- dom_b %>% dplyr::filter(month == 11 & day_mo == 30 | month == 12 | month == 1 | 
                                      month == 2 & day_mo <=29)
range(dom_b_q1$observation.date)

dom_b_q2 <- dom_b %>% dplyr::filter(month == 3 | month == 4 | month == 5 | 
                                      month == 6 & day_mo <= 6 )
range(dom_b_q2$observation.date)

dom_b_q3 <- dom_b %>% dplyr::filter(month ==6 & day_mo >=7 | month == 7 |
                                      month == 8 & day_mo <=9)
range(dom_b_q3$observation.date)

dom_b_q4 <- dom_b %>% dplyr::filter(month == 8 & day_mo >= 10 | month == 9 | month == 10 |
                                      month == 11 & day_mo <= 29)
range(dom_b_q4$observation.date)

# check the numbers add up.  
nrow(dom_b_q1) + nrow(dom_b_q2) + nrow(dom_b_q3) + nrow(dom_b_q4) == nrow(dom_b)


# function to make the raster
presence_rast  <- function(pts, ref_rast){
  p <- cbind(pts$X, pts$Y)
  #  points(p, pch = 18, cex = 0.5) # this just plots the points 
  points_rast <- rasterize(p, ref_rast, value = 1)
}

dom_a_q1_presence_rast <- presence_rast(dom_a_q1, euro_rast)
plot(dom_a_q1_presence_rast, col = "orange")

dom_a_q2_presence_rast <- presence_rast(dom_a_q2, euro_rast)
dom_a_q3_presence_rast <- presence_rast(dom_a_q3, euro_rast)
dom_a_q4_presence_rast <- presence_rast(dom_a_q4, euro_rast)
# repeat for b

dom_b_q1_presence_rast <- presence_rast(dom_b_q1, euro_rast)
dom_b_q2_presence_rast <- presence_rast(dom_b_q2, euro_rast)
dom_b_q3_presence_rast <- presence_rast(dom_b_q3, euro_rast)
dom_b_q4_presence_rast <- presence_rast(dom_b_q4, euro_rast)

 # writeRaster(dom_a_q1_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q1_eco_pres_abs.tif")
 # writeRaster(dom_a_q2_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q2_eco_pres_abs.tif")
 # writeRaster(dom_a_q3_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q3_eco_pres_abs.tif")
 # writeRaster(dom_a_q4_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q4_eco_pres_abs.tif")
 # 
 # writeRaster(dom_b_q1_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q1_eco_pres_abs.tif")
 # writeRaster(dom_b_q2_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q2_eco_pres_abs.tif")
 # writeRaster(dom_b_q3_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q3_eco_pres_abs.tif")
 # writeRaster(dom_b_q4_presence_rast, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q4_eco_pres_abs.tif")

## One issue with this is that you may well just get rasters that overlay the wild bird cases.
## Probably useful to compare the two rasters and see if there are many cells that are independent

## Next part is to count how many unique occurrences there are in each cell.
## We have removed duplicates based on date and location, but possible that we might want to only count an 
## occurrence in a cell once in a day? 

## first will just count them all.

## issue with the code below is that using raster makes the crs tricky. 
# total_count_function <- function(points_obj){
#   p <- cbind(points_obj$X, points_obj$Y)  # make the points object
#   empty_rast <- raster::raster(euro_rast) # create an empty raster
#   empty_rast[] = 0 # set value as 0
#   counts = table(raster::cellFromXY(empty_rast,p)) # make a table of the counts for the cells
#   empty_rast[as.numeric(names(counts))] = counts #add the counts to the raster
#   return(empty_rast)
# }

total_count_function <- function(points_obj){
  p <- cbind(points_obj$X, points_obj$Y)  # make the points object
  empty_rast <- terra::rast(euro_rast ) # create an empty raster
  empty_rast[] = 0 # set value as 0
  counts = table(terra::cellFromXY(empty_rast,p)) # make a table of the counts for the cells
  cells_with_cases <- as.numeric(names(counts))
  cell_vals <- values(empty_rast)
  cell_vals[cells_with_cases] <-counts[] #  counts_vect
  values(empty_rast) <- cell_vals #add the counts to the raster
  return(empty_rast)
}

points_obj <- dom_a_q1

dom_a_q1_rast_counts <- total_count_function(dom_a_q1)
dom_a_q1_rast_counts
plot(dom_a_q1_rast_counts)
table(values(dom_a_q1_rast_counts))

breakpoints <- c(0, 1, 5, 10, 20, 50, 100, 150)
colors <- c("white", rev(RColorBrewer::brewer.pal(7, "Reds")))

euro_shp <- terra::vect("output/euro_map.shp")

plot(dom_a_q1_rast_counts, breaks = breakpoints, col = colors, box = F, axes = F)
plot(euro_shp, add = T)

## create the other rasters
dom_a_q2_rast_counts <- total_count_function(dom_a_q2)
dom_a_q3_rast_counts <- total_count_function(dom_a_q3)
dom_a_q4_rast_counts <- total_count_function(dom_a_q4)

dom_b_q1_rast_counts <- total_count_function(dom_b_q1)
dom_b_q1_rast_counts
dom_b_q2_rast_counts <- total_count_function(dom_b_q2)
dom_b_q3_rast_counts <- total_count_function(dom_b_q3)
dom_b_q4_rast_counts <- total_count_function(dom_b_q4)

# # save the rasters
# writeRaster(dom_a_q1_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q1_eco_counts.tif")
# writeRaster(dom_a_q2_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q2_eco_counts.tif")
# writeRaster(dom_a_q3_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q3_eco_counts.tif")
# writeRaster(dom_a_q4_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_a_q4_eco_counts.tif")
# 
# writeRaster(dom_b_q1_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q1_eco_counts.tif")
# writeRaster(dom_b_q2_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q2_eco_counts.tif")
# writeRaster(dom_b_q3_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q3_eco_counts.tif")
# writeRaster(dom_b_q4_rast_counts, "data/flu_data/prepped_data/domestic_cases_rasters/dom_b_q4_eco_counts.tif")




