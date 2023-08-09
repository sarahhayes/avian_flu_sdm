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
map_boundary <- st_read("output/euro_map.shp")
plot(map_boundary)

eurounion <- st_union(map_boundary)
plot(eurounion)
eurounion

map_boundary_euro <- eurounion


# Defining filepath to save downloaded spatial file
spatial_filepath <- "data/vegetation_data/euro.shp"
# Saving downloaded spatial file on to our computer
st_write(map_boundary_euro, paste0(spatial_filepath), append = F)

plot(map_boundary_euro)

MODIStsp_get_prodnames()

tictoc::tic()
MODIStsp(gui             = FALSE,
         out_folder      = "data/vegetation_data",
         out_folder_mod  = "data/vegetation_data",
         selprod = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
         bandsel = c("NDVI"),
         user            = "sarahhayes" ,
         prod_version    = "061",  
         password        = "NASATigtogs43!",
         start_date      = "2019.01.01", 
         end_date        = "2019.12.31", 
         verbose         = FALSE,
         spatmeth        = "file",
         spafile         = spatial_filepath,
         out_format      = "GTiff")
tictoc::toc()


# Reading in the downloaded NDVI raster data
# NDVI_raster <- raster(here::here("VegetationData/mongolia/VI_16Days_1Km_v6/NDVI/MYD13A2_NDVI_2020_153.tif"))

NDVI_raster <- terra::rast("data/vegetation_data/euro/VI_16Days_1Km_v61/NDVI/MOD13A2_NDVI_2019_001.tif")
NDVI_raster
plot(NDVI_raster)

map_boundary_terra

# Transforming the data

blank_3035 <- terra::rast("output/euro_rast.tif")

#NDVI_raster <- raster::projectRaster(NDVI_raster, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
NDVI_raster2 <- terra::project(NDVI_raster, map_boundary_terra)
NDVI_raster3 <- terra::project(NDVI_raster, blank_3035)

plot(NDVI_raster2)
plot(NDVI_raster3)

# Cropping the data
NDVI_raster <- raster::mask(NDVI_raster, as_Spatial(map_boundary))

# Dividing values by 10000 to have NDVI values between -1 and 1
gain(NDVI_raster) <- 0.0001

# Converting the raster object into a dataframe
NDVI_df <- as.data.frame(NDVI_raster, xy = TRUE, na.rm = TRUE)
rownames(NDVI_df) <- c()

# Visualising using ggplot2
ggplot() +
  geom_raster(
    data = NDVI_df,
    aes(x = x, y = y, fill = MYD13A2_NDVI_2020_153)
  ) +
  geom_sf(data = map_boundary, inherit.aes = FALSE, fill = NA) +
  scale_fill_viridis(name = "NDVI") +
  labs(
    title = "NDVI (Normalized Difference Vegetation Index) in Mongolia",
    subtitle = "01-06-2020",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()





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