# Map of abundance for a selected species for use in poster

rm(list = ls())

library(ebirdst)
library(terra)
library(sf)
library(fields)
library(rnaturalearth)
library(raster)
#remotes::install_github("ropensci/rnaturalearthhires")

## look at the species table provided

sp_df <- ebirdst_runs

#################################################################################
# First do an example with one species only to work out the steps. 

# download data - we will look at the herring gull 
# path <- ebirdst_download(species = "hergul")
path <- ebirdst_download(species = "mallar3")

# load relative abundance raster stack with 52 layers, one for each week
abd <- load_raster(path = path, resolution = "lr") # currently set as lr = low resolution for speed
#abd <- load_raster(path = path, resolution = "hr")

# set the crs we want to use
crs <- "epsg:3035"

#change the projection of the bird data
abd_prj <- terra::project(x = abd, y = crs, method = "near") # ideally (according to help file) would use our template 
# Spatraster of the area we want as y. 

# Now we want to crop
e <- terra::ext(2000000, 6000000, 1000000, 5500000)
crop_abd <- crop(x = abd_prj, y = e )
plot(crop_abd[[1]], axes = F)

# change the colours
# quantiles of non-zero values
v <- values(crop_abd)
v <- v[!is.na(v) & v > 0]
range(v)
 bins <- quantile(v, seq(0, 1, by = 0.1)) # splits into the 10% quantiles
 bins <- c(0, bins) # add a bin for 0
bins <- c(0,0.1,1,2,3,4,5,10,50,100)

# status and trends palette
pal <- abundance_palette(length(bins) - 2)
# add a color for zero
pal <- c("#e6e6e6", pal)
pal <- c("grey", pal)

#alternative colour scheme
pal <- abundance_palette(length(bins) - 1)

# map using the quantile bins
plot(crop_abd[[1]], breaks = bins, col = pal, axes = FALSE)


crop_abd
winter <- mean(crop_abd[[1:13]])
autumn <- mean(crop_abd[[26:39]])
plot(winter, breaks = bins, col = pal, axes = FALSE)
plot(autumn, breaks = bins, col = pal, axes = F)
