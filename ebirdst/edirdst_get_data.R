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

# get reference data from the rnaturalearth package for Europe
wh_europe <- ne_countries(continent = "europe",
                          returnclass = "sf") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

#change the projection of the bird data
abd_prj <- project(x = abd, y = crs, method = "near") # ideally (according to help file) would use our template 
# Spatraster of the area we want as y. 

plot(abd_prj[[1]]) # This is now nicely centred on Europe
# Now we want to crop

e <- terra::ext(2000000, 6000000, 1000000, 5500000)
crop_abd <- crop(x = abd_prj, y = e )
plot(crop_abd[[1]])

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

species_sel <- first_ten[1]
path <- ebirdst_download(species = species_sel)
bird_rast <-  load_raster(path = path, resolution = "lr")

par(mfrow = c(5, 2))
plot(bird_rast[[26]])

# Add directory to store output layers if one does not exist
dir.create(file.path("ebirdst", "output_layers"), showWarnings = FALSE)
for (i in 2:length(first_ten)) {
  species_sel <- first_ten[i]
  path <- ebirdst_download(species = species_sel)
  bird_rast <-  c(bird_rast, load_raster(path = path, resolution = "lr"))
  plot(bird_rast[[52*(i-1) + 26]])
}


proj_abd <- project(x = bird_rast, y = crs, method = "near") 
for (i in 1:length(first_ten)) {
  plot(bird_rast[[52*(i-1) + 1]])
}


crop_abd <- crop(x = proj_abd, y = euro_ext )
for (i in 1:length(first_ten)) {
  plot(crop_abd[[52*(i-1) + 26]])
}

terra::writeRaster(crop_abd, paste("ebirdst/output_layers/first_ten.tif", sep = ""), overwrite = T) 



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
