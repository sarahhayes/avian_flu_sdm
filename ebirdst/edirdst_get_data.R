# We want a script that will bring in the ebird weekly abundance layers
# for all the species that are found in Europe. 

# We want to then save these layers 

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
path <- ebirdst_download(species = "hergul")
# path <- ebirdst_download(species = "mallar3")

# load relative abundance raster stack with 52 layers, one for each week
abd <- load_raster(path = path, resolution = "lr") # currently set as lr = low resolution for speed
# abd <- load_raster(path = path, resolution = "hr")
abd

# select one layer of the raster
abd1 <-abd[[1]]
plot(abd1)

# set the crs we want to use
crs <- "epsg:3035"

# Define an area to crop to
e <- terra::ext(2000000, 6000000, 1000000, 5500000)

# Create a blank raster with appropriate projection and extent
blank_3035 <- rast(crs=crs, extent=e, res=9042.959)

# get reference data from the rnaturalearth package for Europe
wh_europe <- ne_countries(continent = "europe",
                          returnclass = "sf") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

# Test to see if projecting to cropped version of 3035 projection is quicker.
# This needs to all be run as one block to work.

# First try projecting to full global-level coordinates:
#change the projection of the bird data
start.time <- Sys.time()
abd_prj_global <- project(x = abd, y = crs, method = "near") # ideally (according to help file) would use our template 
end.time <- Sys.time()
cat("Projecting to basic 3035 coords took ",
    as.numeric(difftime(end.time, start.time, units="secs")),
    " seconds.")

# Now try projecting to cropped version:
start.time <- Sys.time()
abd_prj_cropped <- project(x = abd, y = blank_3035, method = "near") # ideally (according to help file) would use our template 
end.time <- Sys.time()
cat("Projecting to cropped 3035 coords took ",
    as.numeric(difftime(end.time, start.time, units="secs")),
    " seconds.")

plot(abd_prj_global[[1]]) # This is now nicely centred on Europe
# Now we want to crop
crop_abd <- crop(x = abd_prj_global, y = e )
plot(crop_abd[[1]])

# Compare plot we got from cropping first - should look the same up to resolution.
plot(abd_prj_cropped[[1]])

# Have a look at what we have created
crop_abd
crop1 <- crop_abd[[1]]
crop1


# terra::writeRaster(crop_abd, "ebirdst/output_layers/her_gul.tif", overwrite = T) 
# read back in to check 
# her_st <- terra::rast("ebirdst/output_layers/her_gul.tif")
# plot(her_st)  
  
###################################################################################

## We want to loop this process through all of the species codes that are stored in 
## "ebird/codes_for_europe_clean.R"

##** NB We want to confirm extent and projection and possibly produce the raster that
##** we will be using as a template for all the layers first.  

euro_bird_codes <- read.csv("ebird/codes_for_europe_clean.csv")

#start with small subset and use lr for speed. Change to hr once happy running OK
first_ten <- euro_bird_codes$code[1:10]

# set the crs we want to use
crs <- "epsg:3035"

euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later

# Create a blank raster with appropriate projection and extent
blank_3035 <- rast(crs=crs, extent=euro_ext, res=9042.959)

species_sel <- first_ten[1]
path <- ebirdst_download(species = species_sel)
bird_rast <- load_raster(path = path, resolution = "lr")
bird_rast <- project(x = bird_rast, y = blank_3035, method = "near") 


# Loop over other 9 species:
loop.start <- Sys.time()
for (i in 2:length(first_ten)) {
  species_sel <- first_ten[i]
  path <- ebirdst_download(species = species_sel)
  this_rast <- load_raster(path = path, resolution = "lr")
  this_rast <- project(x = this_rast, y = blank_3035, method = "near") 
  bird_rast <-  c(bird_rast, this_rast)
  time.now <- Sys.time()
  time_remaining <- (length(first_ten) - i) * 
    as.numeric(difftime(time.now, loop.start, units="mins"))/(i-1)
  cat(time.now - loop.start,
      " seconds elapsed since start, estimated ",
      time_remaining,
      " remaining.\n")
}

# Add directory to store output layers if one does not exist
dir.create(file.path("ebirdst", "output_layers"), showWarnings = FALSE)
terra::writeRaster(bird_rast, paste("ebirdst/output_layers/first_ten.tif", sep = ""), overwrite = T) 

par(mfrow = c(5, 2))
# Do some plots of abundance at quarterly intervals
for (week in c(1,14,27,40)){
  for (spec_idx in 1:length(first_ten)) {
    plot(bird_rast[[52*(spec_idx-1) + week]])
  }
}


# Projection
# The raster data are all provided in the same equal area sinusoidal projection as
# NASA MODIS data. While this projection is suitable for analysis, 
# it is not ideal for mapping since it introduces significant distortion. 
# Instead, as part of the Status and Trends workflow, custom species-specific 
# projections are provided that are optimized for the region that the species 
# occurs within. We can access the projection for the species with
# load_fac_map_parameters(), then transform the raster data to this 
# custom projection.

# However, I don't think we will want to use this projection as will need a uniform projection for all species, 
# rather than a species-specific one. 

# load the mapping parameters
# fac_parameters <- load_fac_map_parameters(path)
# crs <- fac_parameters$custom_projection

# transform to the custom projection using nearest neighbour resampling
# abd1_projected <- project(abd1, crs, method = "near")

# map the cropped and projected data
# plot(abd1_projected, axes = FALSE)
# abd1_projected
