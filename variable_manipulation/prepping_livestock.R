## Prepping livestock data
## There are two version of the data. These are DA and AW. Details are 
## in the metadata file saved to the One Drive alongside the data. 

rm(list = ls())
library(tidyverse)
library(terra)

chucks <- terra::rast("data/variables/livestock/chickens/6_Ch_2015_Aw.tif")
chucks_2010 <- terra::rast("data/variables/livestock/chickens/6_Ch_2010_Aw.tif")

chucks # 0.083 degrees = 5 arc mins which is approx 9.28km at the equator
plot(chucks)
chucks_2010
plot(chucks_2010)

# # set projection and extent 
# crs <- "epsg:3035"
# euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later
# 
# # Create a blank raster with appropriate projection and extent
# blank_3035 <- rast(crs=crs, extent=euro_ext, res = 1000)
# blank_3035

# read in the reference raster
#blank_3035 <- terra::rast("output/euro_rast.tif") # for 1k res
blank_3035 <- terra::rast("output/euro_rast_10k.tif") # for 10k res
blank_3035
plot(blank_3035)

chucks_euro <- terra::project(x = chucks, y = blank_3035, method = "bilinear")
chucks_euro
plot(chucks_euro)

chucks_euro_10 <- terra::project(x = chucks_2010, y = blank_3035, method = "bilinear")
chucks_euro_10
plot(chucks_euro_10)


# make a points object using the centre of each pixel from the blank raster
points_3035 <- terra::as.points(blank_3035)
points_3035
#plot(points_3035)

tictoc::tic()
chuck_res <- terra::extract(chucks_euro, points_3035, method = "simple", xy = T)
tictoc::toc()

chuck_res_10 <- terra::extract(chucks_euro_10, points_3035, method = "simple", xy = T)

head(chuck_res)
head(chuck_res_10)
table(chuck_res$`6_Ch_2015_Aw`)
range(chuck_res$`6_Ch_2015_Aw`, na.rm = T)

### Do the same for the ducks 

ducks <- terra::rast("data/variables/livestock/ducks/6_Dk_2015_Aw.tif")
ducks
plot(ducks)

ducks_euro_15 <- terra::project(x = ducks, y = blank_3035, method = "bilinear")
ducks_euro_15 # maximum value is smaller - as would expect at the different resolution

#make points object
tictoc::tic()
duck_res_15 <- terra::extract(ducks_euro_15, points_3035, method = "simple", xy = T)
tictoc::toc()

## rename the columns
chuck_res <- rename(chuck_res, "chicken_density_2015" = "6_Ch_2015_Aw")
range(chuck_res$chicken_density_2015, na.rm = T)
head(chuck_res_10)
chuck_res_10 <- rename(chuck_res_10, "chicken_density_2010" = "6_Ch_2010_Aw")

head(duck_res_15)
duck_res_15 <- rename(duck_res_15, "duck_density_2015" = "6_Dk_2015_Aw")
range(duck_res_15$duck_density_2015, na.rm = T)

duck_chuck_res <- full_join(chuck_res, duck_res_15) %>%
  full_join(chuck_res_10)

head(duck_chuck_res)

## Now bring in the 2010 duck data

ducks_10 <- terra::rast("data/variables/livestock/ducks/6_Dk_2010_Aw.tif")
ducks_10
plot(ducks_10)

ducks_euro_10 <- terra::project(x = ducks_10, y = blank_3035, method = "bilinear")
ducks_euro_10 
plot(ducks_euro_10)

tictoc::tic()
duck_res_10 <- terra::extract(ducks_euro_10, points_3035, method = "simple", xy = T)
tictoc::toc()

## rename the columns

head(duck_res_10)
duck_res_10 <- rename(duck_res_10, "duck_density_2010" = "6_Dk_2010_Aw")
range(duck_res_10$duck_density_2010, na.rm = T)

duck_chuck_res <- full_join(duck_chuck_res, duck_res_10)

head(duck_chuck_res)

duck_chuck_res <- dplyr::select(duck_chuck_res, 
                                all_of(c("ID", "x" , "y","chicken_density_2010", "chicken_density_2015",
                                         "duck_density_2015", "duck_density_2010")))

head(duck_chuck_res)

# save as a large combo csv rather than separate ones

write.csv(duck_chuck_res,
         "variable_manipulation/variable_outputs/chicken_duck_output.csv",
         row.names = F)

#write.csv(chuck_res, "variable_manipulation/variable_outputs/chicken_output.csv")
#write.csv(chuck_res_10, "variable_manipulation/variable_outputs/chicken_output_2010.csv")
#write.csv(duck_res_15, "variable_manipulation/variable_outputs/duck_output_2015.csv")
#write.csv(duck_res_10, "variable_manipulation/variable_outputs/duck_output_2010.csv")


# also save the rasters
# terra::writeRaster(chucks_euro,
#                   "variable_manipulation/variable_outputs/chicken_density.tif")
# terra::writeRaster(chucks_euro_10,
#                   "variable_manipulation/variable_outputs/chicken_density_2010.tif")

# terra::writeRaster(ducks_euro_15,
#                   "variable_manipulation/variable_outputs/duck_density_2015.tif")
# terra::writeRaster(ducks_euro_10,
#                    "variable_manipulation/variable_outputs/duck_density_2010.tif")


# also save the rasters
 terra::writeRaster(chucks_euro,
                   "variable_manipulation/variable_outputs/chicken_density_10kres.tif")
 terra::writeRaster(chucks_euro_10,
                   "variable_manipulation/variable_outputs/chicken_density_2010_10kres.tif")
 terra::writeRaster(ducks_euro_15,
                   "variable_manipulation/variable_outputs/duck_density_2015_10kres.tif")
 terra::writeRaster(ducks_euro_10,
                    "variable_manipulation/variable_outputs/duck_density_2010_10kres.tif")

