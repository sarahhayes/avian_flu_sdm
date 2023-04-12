## running through the vignette from: https://ebird.github.io/ebirdst/articles/rasters.html

rm(list = ls())

library(ebirdst)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(geodata)
library(ggplot2)
library(fields)
extract <- terra::extract

# download the example yellow-bellied sapsucker data
# this simplified dataset doesn't require an access key
path <- ebirdst_download("example_data")
# alternatively, you can use the following to get the data path
# if the data package has already been downloaded
path <- get_species_path("example_data")

# load seasonal mean relative abundance at low res
abd_seasonal <- load_raster(path, 
                            product = "abundance", 
                            period = "seasonal",
                            metric = "mean",
                            resolution = "lr")

# get the seasons corresponding to each layer
names(abd_seasonal)


# extract just the breeding season relative abundance
abd_breeding <- abd_seasonal[["breeding"]]
# We can get the dates and quality scores associated with each of these seasons
# by filtering the ebirdst_runs data frame.

ebirdst_runs %>% 
  # note that the example data are for yellow-bellied sapsucker
  filter(common_name == "Yellow-bellied Sapsucker") %>% 
  glimpse()


# When plotting could try just using plot.... 
plot(abd_breeding, axes = FALSE)
#This doesn't work - vignette explains why!

# Part of the reason is that will do a global map unless we crop
# boundary polygon for michigan
mi <- ne_states(iso_a2 = "US", returnclass = "sf") %>% 
  filter(postal == "MI") %>% 
  # project to same coordinate reference system as the raster data
  st_transform(st_crs(abd_seasonal))

# crop data to michigan
abd_breeding_mi <- crop(abd_breeding, mi)

# map the cropped data
plot(abd_breeding_mi, axes = FALSE)

## Next section in vignette talks about the projection of the data

# Projection
# The raster data are all provided in the same equal area sinusoidal projection as
# NASA MODIS data. While this projection is suitable for analysis, 
# it is not ideal for mapping since it introduces significant distortion. 
# Instead, as part of the Status and Trends workflow, custom species-specific 
# projections are provided that are optimized for the region that the species 
# occurs within. We can access the projection for Yellow-bellied Sapsucker with
# load_fac_map_parameters(), then transform the raster data to this 
# custom projection.

# load the mapping parameters
fac_parameters <- load_fac_map_parameters(path)
crs <- fac_parameters$custom_projection

# transform to the custom projection using nearest neighbor resampling
abd_projected <- project(abd_breeding_mi, crs, method = "near")

# map the cropped and projected data
plot(abd_projected, axes = FALSE)

# change the bins for the abundance

# quantiles of non-zero values
v <- values(abd_projected)
v <- v[!is.na(v) & v > 0]
bins <- quantile(v, seq(0, 1, by = 0.1))
# add a bin for 0
bins <- c(0, bins)

# status and trends palette
pal <- abundance_palette(length(bins) - 2)
# add a color for zero
pal <- c("#e6e6e6", pal)

# map using the quantile bins
plot(abd_projected, breaks = bins, col = pal, axes = FALSE)


# Finally, we’ll add state and country boundaries to provide some context. 
# The R package rnaturalearth is an excellent source of attribution 
# free contextual GIS data.

# natural earth boundaries
countries <- ne_countries(returnclass = "sf") %>% 
  st_geometry() %>% 
  st_transform(crs)
states <- ne_states(iso_a2 = "US", returnclass = "sf") %>% 
  st_geometry() %>% 
  st_transform(crs)

# define the map extent with the michigan polygon
mi_ext <- mi %>% 
  st_geometry() %>% 
  st_transform(crs)
plot(mi_ext)
# add basemap
plot(countries, col = "#cfcfcf", border = "#888888", add = TRUE)
# add data
plot(abd_projected, 
     breaks = bins, col = pal, 
     axes = FALSE, legend = FALSE, add = TRUE)
# add boundaries
plot(countries, col = NA, border = "#888888", lwd = 3, add = TRUE)
plot(states, col = NA, border = "#888888", add = TRUE)

# add legend using the fields package
# label the bottom, middle, and top
labels <- quantile(bins, c(0, 0.5, 1))
label_breaks <- seq(0, 1, length.out = length(bins))
image.plot(zlim = c(0, 1), breaks = label_breaks, col = pal,
           smallplot = c(0.90, 0.93, 0.15, 0.85),
           legend.only = TRUE,
           axis.args = list(at = c(0, 0.5, 1), 
                            labels = round(labels, 2),
                            col.axis = "black", fg = NA,
                            cex.axis = 0.9, lwd.ticks = 0,
                            line = -0.5))
