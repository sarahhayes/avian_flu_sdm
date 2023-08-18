### 09/08/2023
### Looking at ways to get the land cover data
## thi script is based on this webpage https://rspatial.org/modis/2-download.html

#install.packages("remotes")
#remotes::install_github("rspatial/luna")

rm(list = ls())

library(terra)
library(luna)
library(sf)


# lists all products that are currently searchable
prod <- getProducts()
head(prod)

modis <- getProducts("^MOD|^MYD|^MCD")
head(modis)

product <- "MCD12Q1"
# product <- "MOD09A1"
# To learn more about a specific product you can launch a webpage
# productInfo(product)

# Once we finalize the product we want to use, 
# we define some parameters for the data we want: 
# product name, start and end date, and area of interest.

start <- "2020-01-01"
end <- "2020-12-31"

# To define the area of interest, we can define a spatial extent,
# or use an object that has an extent. 

# map_boundary <- st_read("output/euro_map.shp")
# map_boundary <- terra::vect("output/euro_map.shp")
# eurounion <- st_union(map_boundary)
# plot(eurounion)
# eurounion
# 
# map_boundary_euro <- eurounion
# plot(map_boundary_euro)

# we need to extent of the area we are interested in in long/lat. 
map_rast <- terra::rast("output/euro_rast.tif")
map_rast_lonlat <- terra::project(map_rast, "epsg:4326")
map_rast_lonlat

# Let’s now find out what MODIS data is available for this area. 
# We can search the data available from a NASA server

mf <- luna::getModis(product, start, end, aoi=map_rast_lonlat, 
                     download = FALSE, version = "061")
mf


 #To download the tiles, usually you would download them to a folder 
# where you save the data for your project. 

datadir <- file.path("data/landcover_data/euro")
dir.create(datadir, showWarnings=FALSE)

lc_data <- luna::getModis(product, start, end, aoi=map_rast_lonlat, download=TRUE,
                     path=datadir, version = "061",
                     username="sarahhayes", password="NASATigtogs43!")


# Now that we have downloaded some MODIS data, we can explore and visualize it.
# First create a SpatRaster object from the file created on the previous page.

# this is looking at just one tile
lc_path <- file.path(datadir, "MCD12Q1.A2020001.h15v03.061.2022171190052.hdf")
#library(terra)
r <- terra::rast(lc_path[1])
r
plot(r)
plot(r$LC_Type1)
r$LC_Type1

# I want all the tiles and just the first layer from each of the files in the hdf folder
files <- dir("data/landcover_data/euro", pattern = ".hdf")
files_vect <- paste("data/landcover_data/euro/", files, sep = "")

for (i in 1:length(files_vect)) {
  rr <- terra::rast(files_vect[[i]], lyrs = "LC_Type1")
  terra::writeRaster(rr, paste("data/landcover_data/landcover_LC_Type1/",i, ".tif", sep ="" ),
                     overwrite = T)
}

## combine them into a single raster
vrt(
  x = list.files(path = "data/landcover_data/landcover_LC_Type1",
                 pattern = "*.tif$", full.names = TRUE), 
  filename = "dem.vrt", overwrite = T
)

# afterwards read it as if it was a normal raster:
dem <- rast("dem.vrt")
dem
plot(dem)

## save this raster
# terra::writeRaster(dem, "data/landcover_data/landcover_type1_full_raster.tif")
dem <- terra::rast("data/landcover_data/landcover_type1_full_raster.tif")
plot(dem)
## Now we need to change projection and crop

lc1_rast <- terra::project(dem, map_rast)
lc1_rast
plot(lc1_rast)

lc1_factor <- terra::as.factor(lc1_rast)
plot(lc1_factor)
lc1_factor

# save this raster 
#writeRaster(lc1_factor, "variable_manipulation/variable_outputs/landcover_output_full.tif")

# make a points object using the centre of each pixel from the blank raster
points_3035 <- terra::as.points(map_rast)
points_3035

tictoc::tic()
lc_res_full_fact <- terra::extract(lc1_factor, points_3035, method = "simple", xy = T)
tictoc::toc()
summary(lc_res_full_fact[is.na(lc_res_full_fact)])

# write.csv(lc_res_full_fact,
#           "variable_manipulation/variable_outputs/landcover_output_full_factor.csv",
#           row.names = F)
#  


####
masked <- terra::mask(lc1_rast, map_rast)
plot(masked)
table(values(masked))

masked_int <- terra::as.int(masked)
masked_int
plot(masked_int)
head(masked_int)
table(values(masked_int))

names(masked_int)
values(masked_int)

masked_factor <- terra::as.factor(masked)
masked_factor
plot(masked_factor)

names(masked_factor) <- "landcover" 
names(masked_factor)

# make a points object using the centre of each pixel from the blank raster
points_3035 <- terra::as.points(map_rast)
points_3035

tictoc::tic()
lc_res <- terra::extract(masked_factor, points_3035, method = "simple", xy = T)
tictoc::toc()

 # write.csv(lc_res,
 #          "variable_manipulation/variable_outputs/landcover_output.csv",
 #          row.names = F)
 # 