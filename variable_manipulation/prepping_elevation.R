# Previously downloaded files from USGS Earth Explorer for elevation using SRTM void filled
# However, these don't extend above 60 latitude so do not include Scandinavia when re-assesmble
# (Script to reassemble at bottom of this file)
# As such looking at alternatives. 

#install.packages("geodata")

rm(list =ls())
library(tidyverse)
library(geodata)
library(terra)

global_elev_5 <- elevation_global(res = 5, path = "data/variables/elevation/geodata_global_res5")
global_elev_5
# global raster in standard long/lat. Resolution is 0.083. This is in degrees. 
# This is 5 arc minutes which is approx 10km

# read in reference raster
euro_rast_10 <- terra::rast("output/euro_rast_10k.tif") # first look on larger scale

elev_euro_5_10 <- terra::project(x = global_elev_5, y = euro_rast_10, method = "mode") #
# for now using mode whilst working it out. 
elev_euro_5_10
plot(elev_euro_5_10)

### That seems reasonable and quick. 
### Now let's look at a smaller res for both parts

#euro_rast <- terra::rast("output/euro_rast.tif") # for 1k res
euro_rast <- terra::rast("output/euro_rast_10k.tif") # for 10k res

global_elev_0.5 <- elevation_global(res = 0.5, path = "data/variables/elevation/geodata_global_res0.5")
global_elev_0.5 # 0.00833 degrees (0r 30 seconds).  from this webpage: https://longlovemyu.com/metrics_gis/
## this is about 928m at the equator. 

head(global_elev_0.5)
plot(global_elev_0.5)

##  first_looking at max.

elev_euro_0.5_max <- terra::project(x = global_elev_0.5, y = euro_rast, method = "max") #
elev_euro_0.5_max
head(elev_euro_0.5_max)
plot(elev_euro_0.5_max)

# then minimum. 
elev_euro_0.5_min <- terra::project(x = global_elev_0.5, y = euro_rast, method = "min") #
elev_euro_0.5_min
plot(elev_euro_0.5_min)

# next thing is to look at the difference between min and max values in a cell. 
elev_euro_0.5_diff <- elev_euro_0.5_max - elev_euro_0.5_min
elev_euro_0.5_diff # diff of 2 km in places in 1k res and up to 3.5 in 10k res 
plot(elev_euro_0.5_diff)

# wondering about also incorporating some kind of average to get an idea of the predominat elevation at that area.
# I think mode might work best for this? 
elev_euro_0.5_mode <- terra::project(x = global_elev_0.5, y = euro_rast, method = "mode")
elev_euro_0.5_mode # diff of 2 km in places in 1k res and up to 3.5 in 10k res 
plot(elev_euro_0.5_mode)

dev.off()
par(mfrow = c(2,2))
plot(elev_euro_0.5_min, main = "min")
plot(elev_euro_0.5_max, main = "max")
plot(elev_euro_0.5_diff, main = "diff")
plot(elev_euro_0.5_mode, main = "mode")



#save these rasters
# terra::writeRaster(elev_euro_0.5_min, 
#                    "variable_manipulation/variable_outputs/elevation_min.tif")
# terra::writeRaster(elev_euro_0.5_max, 
#                    "variable_manipulation/variable_outputs/elevation_max.tif")
# terra::writeRaster(elev_euro_0.5_diff, 
#                    "variable_manipulation/variable_outputs/elevation_diff.tif")

# terra::writeRaster(elev_euro_0.5_min,
#                    "variable_manipulation/variable_outputs/elevation_min_10kres.tif")
# terra::writeRaster(elev_euro_0.5_max,
#                    "variable_manipulation/variable_outputs/elevation_max_10kres.tif")
# terra::writeRaster(elev_euro_0.5_diff,
#                    "variable_manipulation/variable_outputs/elevation_diff_10kres.tif")
# terra::writeRaster(elev_euro_0.5_mode,
#                    "variable_manipulation/variable_outputs/elevation_mode_10kres.tif")



### extract the data we want at the centre points of our grid
# make a points object using the centre of each pixel from the blank raster
points_3035 <- terra::as.points(euro_rast)
points_3035

tictoc::tic()
elev_res_max <- terra::extract(elev_euro_0.5_max, points_3035, method = "simple", xy = T)
tictoc::toc()
head(elev_res_max)
elev_res_max <- rename(elev_res_max, "elev_max" = "wc2.1_30s_elev")

elev_res_min <- terra::extract(elev_euro_0.5_min, points_3035, method = "simple", xy = T)
head(elev_res_min)
elev_res_min <- rename(elev_res_min, "elev_min" = "wc2.1_30s_elev")


elev_res_diff <- terra::extract(elev_euro_0.5_diff, points_3035, method = "simple", xy = T)
head(elev_res_diff)
elev_res_diff <- rename(elev_res_diff, "elev_diff" = "wc2.1_30s_elev")

elev_res_all <- full_join(elev_res_min, elev_res_max,  by = c("ID", "x","y" )) %>%
  full_join(elev_res_diff)
head(elev_res_all)

elev_res_all <- dplyr::select(elev_res_all, all_of(c("ID", "x", "y", "elev_min",
                                                      "elev_max", "elev_diff")))


# write.csv(elev_res_all, "variable_manipulation/variable_outputs/elevation_outputs.csv")

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

