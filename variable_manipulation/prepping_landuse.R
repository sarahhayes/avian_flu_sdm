### Prepping landuse data
### 27th July 2023

### There are 2 options for this:
### 1) Use the USGS Earth explorer to download the MODIS dat and process manually
### 2) Use geoData package to download landuse. This has fewer categories 

rm(list =ls())
# library(tidyverse)
# library(geodata)
# library(terra)
# 
# land_cover <- landcover(var = "trees", path = "data/variables")
# plot(land_cover)
# land_cover
## for this you have to call in each of the landcover types separately. 
## It then tells you what proportion of that landcover type is in each cell
## this is an option, but look at an alternative too


### using MODIS. F
### info found here: https://rspatialdata.github.io/land_cover.html

library(MODIStsp)
MODIStsp_get_prodlayers("MCD12Q1")

# remotes::install_github("wmgeolab/rgeoboundaries")
# install.packages("sf")
# library(rgeoboundaries)
library(sf)

# Downloading the country boundary of France
#map_boundary <- geoboundaries("France")

#map_boundary <- terra::vect("output/euro_map.shp")
map_boundary <- st_read("output/euro_map.shp")
plot(map_boundary)

eurounion <- st_union(map_boundary)
plot(eurounion)
eurounion

map_boundary_euro <- eurounion


# Defining filepath to save downloaded spatial file
spatial_filepath <- "LandCoverData/euro.shp"
# Saving downloaded spatial file on to our computer
st_write(map_boundary_euro, paste0(spatial_filepath))


library(MODIStsp)
tictoc::tic()
MODIStsp(gui             = FALSE,
         out_folder      = "data/landcover_data",
         out_folder_mod  = "data/landcover_data",
         selprod         = "LandCover_Type_Yearly_500m (MCD12Q1)",
         bandsel         = "LC1", 
         user            = "sarahhayes" ,
         password        = "NASATigtogs43!",
         start_date      = "2019.01.01", 
         end_date        = "2019.12.31", 
         verbose         = FALSE,
         spatmeth        = "file",
         spafile         = spatial_filepath,
         out_format      = "GTiff")
tictoc::toc()

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