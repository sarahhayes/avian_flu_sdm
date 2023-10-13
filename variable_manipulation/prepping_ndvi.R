### Prepping NDVI data
### 8th August 2023

## Plan to use the MODIS data 

rm(list =ls())
library(tidyverse)
library(sf)
library(MODIStsp)
library(terra)
library(ggplot2)
library(viridis)
#library(rgdal)

### using MODIS
### info found here: https://docs.ropensci.org/MODIStsp/articles/analyze.html
### and here: https://rspatialdata.github.io/vegetation.html

MODIStsp_get_prodlayers("(M*D13A2)")

# Downloading the boundary 
#map_boundary_euro <- rgeoboundaries::geoboundaries("Mongolia")
map_boundary_terra <- terra::vect("output/euro_map.shp")
ext(map_boundary_terra)

rast <- terra::rast("output/euro_rast.tif")
rast
map_lonlat <- terra::project(rast, "epsg:4326")
map_lonlat

map_boundary <- st_read("output/euro_map.shp")
plot(map_boundary)

#map_boundary_lonlat <- st_transform(map_boundary, crs = crs(map_lonlat))

eurounion <- st_union(map_boundary)
#eurounion <- st_union(map_boundary_lonlat)
plot(eurounion)
eurounion

map_boundary_euro <- eurounion
plot(map_boundary_euro)

# Defining filepath to save downloaded spatial file
spatial_filepath <- "data/vegetation_data/euro.shp"
#spatial_filepath <- "data/vegetation_data/euro_lonlat.shp"

# Saving downloaded spatial file on to our computer
st_write(map_boundary_euro, paste0(spatial_filepath))

MODIStsp_get_prodnames()

## NB - as of 13/10/2023 the code below does not work with the updated version of R!!! 
## Can use bbox as the bounding box but hard to work with the MODIS projection. 
## As such this was run on my old laptop which has an old version of R on it.

tictoc::tic()
MODIStsp(gui             = FALSE,
         out_folder      = "data/vegetation_data/vegetation_data_2022",
         out_folder_mod  = "data/vegetation_data/vegetation_data_2022",
         selprod = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
         bandsel = c("NDVI"),
         user            = "sarahhayes" ,
         prod_version    = "061",  
         password        = "NASATigtogs43!",
         start_date      = "2022.01.01", 
         end_date        = "2022.12.31", 
         verbose         = FALSE,
         spatmeth        = "file",
         spafile         = spatial_filepath,
         output_proj = "epsg:3035",
         out_format      = "GTiff")
tictoc::toc()


# Reading in the downloaded NDVI raster data
# NDVI_raster <- raster(here::here("VegetationData/mongolia/VI_16Days_1Km_v6/NDVI/MYD13A2_NDVI_2020_153.tif"))

# A google search suggests the MOD is Terra and the MYD is Aqua so assume we want Terra? 

NDVI_raster <- terra::rast("data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_129.tif")
NDVI_raster
plot(NDVI_raster)

# MYD_raster <- terra::rast("data/vegetation_data/euro/VI_16Days_1Km_v61/NDVI/MYD13A2_NDVI_2019_009.tif")
# MYD_raster
# plot(MYD_raster)

# Try and import them all in a stack
# rastlist <- list.files(path = "data/vegetation_data/euro/VI_16Days_1Km_v61/NDVI", 
#                        pattern='.tif$', all.files= T, full.names= T)

# ndvi_stack <- terra::rast(rastlist)
# this works but we probably want them in the separate seasons so do as 4 sections. 
# these files are generated q 16 days so won't fit exactly into quarters. 

first_quart <- c("data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_001.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_017.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_033.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_049.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_065.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_081.tif")


second_quart <- c("data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_097.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_113.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_129.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_145.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_161.tif",
                 "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_177.tif")

third_quart <- c("data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_193.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_209.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_225.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_241.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_257.tif")

fourth_quart <- c("data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_273.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_289.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_305.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_321.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_337.tif",
                  "data/vegetation_data/vegetation_data_2022/VI_16Days_1Km_v61/NDVI/battling_with_ndvi/MOD13A2_NDVI_2022_353.tif")

first_stack <- terra::rast(first_quart)
second_stack <- terra::rast(second_quart)
third_stack <- terra::rast(third_quart)
fourth_stack <- terra::rast(fourth_quart)

## calculate the mean
# for the very early and very late parts of the year, some of the time there
## are missing data for the far north. Poss related to snow?? 
## So calculate the mean but remove any NAs. 
first_stack_mean <- mean(first_stack, na.rm = T)
plot(first_stack_mean)

second_stack_mean <- mean(second_stack)
plot(second_stack_mean)

third_stack_mean <- mean(third_stack)
plot(third_stack_mean)

fourth_stack_mean <- mean(fourth_stack, na.rm = T)
plot(fourth_stack_mean)

# Transforming the data
blank_3035 <- terra::rast("output/euro_rast.tif")

first_stack_prj <- terra::project(first_stack_mean, blank_3035)
# as we haven't specified a method, the default will be bilinear interpolation as 
# the value is numeric
first_stack_prj
plot(first_stack_prj)

second_stack_prj <- terra::project(second_stack_mean, blank_3035)
third_stack_prj <- terra::project(third_stack_mean, blank_3035)
fourth_stack_prj <- terra::project(fourth_stack_mean, blank_3035)

second_stack_mean
plot(second_stack_prj)

# Cropping the data
#NDVI_raster <- raster::mask(NDVI_raster, as_Spatial(map_boundary))

# Dividing values by 10000 to have NDVI values between -1 and 1
first_stack_prj <- first_stack_prj/10000
first_stack_prj
plot(first_stack_prj)

second_stack_prj <- second_stack_prj/10000
third_stack_prj <- third_stack_prj/10000
fourth_stack_prj <- fourth_stack_prj/10000

plot(fourth_stack_prj)

# save the rasters
terra::writeRaster(first_stack_prj, 
                    "variable_manipulation/variable_outputs/ndvi_first_quart_2022.tif")
terra::writeRaster(second_stack_prj, 
                    "variable_manipulation/variable_outputs/ndvi_second_quart_2022.tif")
terra::writeRaster(third_stack_prj, 
                    "variable_manipulation/variable_outputs/ndvi_third_quart_2022.tif")
terra::writeRaster(fourth_stack_prj, 
                    "variable_manipulation/variable_outputs/ndvi_fourth_quart_2022.tif")


## As don't use the points - below not updated for 2022. CSV is thus still 2019

## extract the data 
points_3035 <- terra::as.points(blank_3035)
points_3035

tictoc::tic()
ndvi_res_first <- terra::extract(first_stack_prj, points_3035, method = "simple", xy = T)
tictoc::toc()

ndvi_res_second <- terra::extract(second_stack_prj, points_3035, method = "simple", xy = T)
ndvi_res_third <- terra::extract(third_stack_prj, points_3035, method = "simple", xy = T)
ndvi_res_fourth <- terra::extract(fourth_stack_prj, points_3035, method = "simple", xy = T)

ndvi_res_first <- rename(ndvi_res_first, mean_ndvi_first_quarter = mean)
ndvi_res_second <- rename(ndvi_res_second, mean_ndvi_second_quarter = mean)
ndvi_res_third <- rename(ndvi_res_third, mean_ndvi_third_quarter = mean)
ndvi_res_fourth <- rename(ndvi_res_fourth, mean_ndvi_fourth_quarter = mean)


full_ndvi_res <- inner_join(ndvi_res_first, ndvi_res_second) %>%
  inner_join(ndvi_res_third) %>%
  inner_join(ndvi_res_fourth)
 
head(full_ndvi_res)

full_ndvi_res <- dplyr::select(full_ndvi_res, all_of(c("ID", "x", "y", 
                                                        "mean_ndvi_first_quarter",
                                                        "mean_ndvi_second_quarter",
                                                        "mean_ndvi_third_quarter",
                                                        "mean_ndvi_fourth_quarter")))

head(full_ndvi_res)

#write.csv(full_ndvi_res, "variable_manipulation/variable_outputs/ndvi_output.csv")






# # Converting the raster object into a dataframe
# NDVI_df <- as.data.frame(NDVI_raster, xy = TRUE, na.rm = TRUE)
# rownames(NDVI_df) <- c()
# 
# # Visualising using ggplot2
# ggplot() +
#   geom_raster(
#     data = NDVI_df,
#     aes(x = x, y = y, fill = MYD13A2_NDVI_2020_153)
#   ) +
#   geom_sf(data = map_boundary, inherit.aes = FALSE, fill = NA) +
#   scale_fill_viridis(name = "NDVI") +
#   labs(
#     title = "NDVI (Normalized Difference Vegetation Index) in Mongolia",
#     subtitle = "01-06-2020",
#     x = "Longitude",
#     y = "Latitude"
#   ) +
#   theme_minimal()
# 
# 
# 
# 

# library(MODIStsp)
# MODIStsp(gui             = FALSE,
#          out_folder      = "data/temp",
#          out_folder_mod  = "data/temp",
#          selprod         = "LandCover_Type_Yearly_500m (MCD12Q1)",
#          bandsel         = "LC1", 
#          user            = "sarahhayes" ,
#          password        = "NASATigtogs43!",
#          start_date      = "2019.01.01", 
#          end_date        = "2019.12.31", 
#          verbose         = FALSE,
#          spatmeth        = "file",
#          spafile         = spatial_filepath,
#          out_format      = "GTiff")
# 


# remotes::install_github("wmgeolab/rgeoboundaries")
library(raster)
library(here)
library(ggplot2)
library(viridis)
library(dplyr)

# Downloading the boundary of Zimbabwe
#map_boundary <- geoboundaries("Zimbabwe")

# Reading in the downloaded landcover raster data
IGBP_raster <- raster(here::here("data/landcover_data/euro_map/LandCover_Type_Yearly_500m_v6/LC1/MCD12Q1_LC1_2019_001.tif"))

# Transforming data
IGBP_raster <- projectRaster(IGBP_raster, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Cropping data
IGBP_raster <- raster::mask(IGBP_raster, as_Spatial(map_boundary))

# Converting the raster object into a dataframe and converting the IGBP classification into a factor
IGBP_df <- as.data.frame(IGBP_raster, xy = TRUE, na.rm = TRUE) %>%
  mutate(MCD12Q1_LC1_2019_001 = as.factor(round(MCD12Q1_LC1_2019_001)))
rownames(IGBP_df) <- c()
# Renaming IGBP classification levels
levels(IGBP_df$MCD12Q1_LC1_2019_001) <- c( "Evergreen needleleaf forests",
                                           "Evergreen broadleaf forests",
                                           "Deciduous needleleaf forests",
                                           "Deciduous broadleaf forests",
                                           "Mixed forests",
                                           "Closed shrublands",
                                           "Open shrublands",
                                           "Woody savannas",
                                           "Savannas",
                                           "Grasslands",
                                           "Permanent wetlands",
                                           "Croplands",
                                           "Urban and built-up lands",
                                           "Cropland/natural vegetation mosaics",
                                           "Snow and ice",
                                           "Barren",
                                           "Water bodies")
# Visualising using ggplot2
ggplot() + 
  geom_raster(data = IGBP_df,
              aes(x = x, y = y, fill = MCD12Q1_LC1_2019_001)) +
  geom_sf(data = map_boundary, inherit.aes = FALSE, fill = NA) +
  scale_fill_viridis(name = "Land Cover Type", discrete = TRUE) +
  labs(title = "Land Cover classification in Zimbabwe",
       subtitle = "01-01-2019 - 31-12-2019",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()



# make blank raster
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later

blank_raster <- rast(crs=crs, extent=euro_ext, res = 1000)

elev_euro <- terra::project(x = global_elev, y = blank_raster, method = "near")
elev_euro

plot(elev_euro)