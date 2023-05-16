### Climate data 21/04/2023


# Using the historical data from WorldClim https://www.worldclim.org/data/monthlywth.html 
# Notes on .tif files from https://inbo.github.io/tutorials/tutorials/spatial_standards_raster/

rm(list = ls())

library(tidyverse)
library(raster)
library(rnaturalearth)

test <- raster("data/variables/climate/wc2.1_2.5m_tmin_2018-12.tif")
test
spplot(test)

# below are the coordinates that are specified when I Google the boundaries of Europe.
eur_extent <- c(-25, 45, 34, 72 )

test_europe <- crop(test, eur_extent)
test_europe

spplot(test_europe)

# We might want to settle on a shape file and use tht? 
