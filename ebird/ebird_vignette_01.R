## Preparing and looking at some eBird data

# Instructions taken from 
# https://cornelllabofornithology.github.io/ebird-best-practices/index.html

rm(list = ls())

#install.packages("remotes")
#remotes::install_github("mstrimas/ebppackages")

#auk::auk_set_ebd_path("ebird/ebird_data")

library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)

map_proj <- "WGS84"# st_crs(102003) # I can't get this crs to work so I have used a different one for now
ne_land <- read_sf("ebird/ebird_data/book_example_data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf("ebird/ebird_data/book_example_data/gis-data.gpkg", "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("ebird/ebird_data/book_example_data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("ebird/ebird_data/book_example_data/gis-data.gpkg", "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

# read in the ebird data from the example data
ebird <- read.csv("ebird/ebird_data/book_example_data/ebd_woothr_june_bcr27_zf.csv")

# prepare ebird data for mapping
ebird_sf <- ebird %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = map_proj) %>% 
  select(species_observed)

# map
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
plot(st_geometry(ebird_sf), col = NA)
# contextual gis data
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
plot(bcr, col = "#cccccc", border = NA, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
# ebird observations
# not observed
plot(st_geometry(ebird_sf),
     pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
     add = TRUE)
# observed
plot(filter(ebird_sf, species_observed) %>% st_geometry(),
     pch = 19, cex = 0.3, col = alpha("#4daf4a", 1),
     add = TRUE)
# legend
legend("bottomright", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("eBird checklists", "Wood Thrush sightings"),
       pch = 19)
box()
par(new = TRUE, mar = c(0, 0, 3, 0))
title("Wood Thrush eBird Observations\nJune 2010-2019, BCR 27")
