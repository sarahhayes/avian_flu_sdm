## Prepping livestock data
## There are two version of the data. These are DA and AW. Details are 
## in the metadata file saved to the One Drive alongside the data. 

rm(list = ls())
library(tidyverse)
library(terra)

chucks <- terra::rast("data/variables/livestock/chickens/6_Ch_2015_Aw.tif")
chucks
plot(chucks)

# set projection and extent 
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later

# Create a blank raster with appropriate projection and extent
blank_3035 <- rast(crs=crs, extent=euro_ext, res = 1000)
blank_3035

chucks_euro <- terra::project(x = chucks, y = blank_3035, method = "near")
chucks_euro

# make a points object using the centre of each pixel from the blank raster
points_3035 <- terra::as.points(blank_3035)
points_3035

tictoc::tic()
chuck_res <- terra::extract(chucks_euro, points_3035, method = "simple", xy = T)
tictoc::toc()

### Do the same for the ducks 

ducks <- terra::rast("data/variables/livestock/ducks/6_Dk_2015_Aw.tif")
ducks
plot(ducks)

ducks_euro <- terra::project(x = ducks, y = blank_3035, method = "near")
ducks_euro # maximum value is smaller - as would expect at the different resolution

tictoc::tic()
duck_res <- terra::extract(ducks_euro, points_3035, method = "simple", xy = T)
tictoc::toc()


#write.csv(chuck_res, "variable_manipulation/variable_outputs/chicken_output.csv")
#write.csv(duck_res, "variable_manipulation/variable_outputs/duck_output.csv")
