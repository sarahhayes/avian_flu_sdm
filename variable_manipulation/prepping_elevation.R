# Previously downloaded files from USGS Earth Explorer for elevation using SRTM void filled
# However, these don't extend above 60 latitude so do not include Scandinavia when re-assesmble
# (Script to reassemble at bottom of this file)
# As such looking at alternatives. 

#install.packages("geodata")

rm(list =ls())
library(tidyverse)
library(geodata)
library(terra)

global_elev <- elevation_global(res = 0.5, path = "data/variables")

# make blank raster
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later

blank_raster <- rast(crs=crs, extent=euro_ext, res = 1000)

elev_euro <- terra::project(x = global_elev, y = blank_raster, method = "near")
elev_euro

plot(elev_euro)

###------------------------------------------------------------------
### Alternative using elevatr package

#install.packages("elevatr")
library(elevatr)

# make the template raster 
crs <- "epsg:3035"

blank_rast <- raster::raster(xmn = 2000000, xmx = 6000000, 
                             ymn = 1000000, ymx = 5500000, 
                             crs = crs, res = 1000)

tictoc::tic()
try <- get_elev_raster(locations = blank_rast, z = 7)
tictoc::toc()
# need to have a look at all the arguments to see what is going on 
plot(try)
try

# see if takes longer to get more zoom (400-500m approx). Get warning if try to do zoom =10
tictoc::tic()
try_2 <- get_elev_raster(locations = blank_rast, z = 8)
tictoc::toc()
plot(try_2)
try_2

# Works, but need to look at details to work out best way to get data and to re-project etc

# Will need to work out how to deal with cells that contain sea.
######################################################################################

## Re-assembling the separate elevation tif files into one file

rm(list = ls())

library(tidyverse)
library(terra)

# using the code on this page to join the tifs together
# https://stackoverflow.com/questions/50234139/using-mosaic-in-r-to-merge-multiple-geotiff

# first build a virtual raster file 

vrt(
  x = list.files(path = "data/variables/elevation/raw_3s_35N72N25W65E_secondtry",
                 pattern = "*.tif$", full.names = TRUE), 
  filename = "dem.vrt", overwrite = T
)

# afterwards read it as if it was a normal raster:
dem <- rast("dem.vrt")
dem
plot(dem)

# Not sure why it doesn't appear to include Scandinavia as that was definitely inside 
# the box I selected from USGS

# set projection and extent 
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later

# Create a blank raster with appropriate projection and extent
blank_3035 <- rast(crs=crs, extent=euro_ext, res = 1000)
blank_3035

####-----------------------------------------------------------
# just trialing this initially. Actually may not want to change the resolution of this
# as if go to 1000 we will lose detail? Explore/consider other options
# At present just trying to find the missing sections!! 

elev_euro <- terra::project(x = dem, y = blank_3035, method = "near")
elev_euro
plot(elev_euro)

