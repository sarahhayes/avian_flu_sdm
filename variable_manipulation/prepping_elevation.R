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
