## Alternative way to get landcover

# first bit from https://flograttarola.com/post/modis-downloads/

rm(list = ls())

library(MODIStsp)
library(terra)
library(sf)
library(tidyverse)


MODIStsp(gui             = FALSE,
         out_folder      = 'data/landcover_data',
         out_folder_mod  = 'data/landcover_data',
         selprod         = 'LandCover_Type_Yearly_500m (MCD12Q1)',
         bandsel         = 'LC1',
         sensor          = 'Terra',
         user            = 'sarahhayes',
         password        = 'NASATigtogs43!',  
         start_date      = '2020.01.01',
         end_date        = '2020.12.31',
         verbose         = TRUE,
         #bbox            =  c(2000000, 1000000, 6000000, 5500000), # xmin, ymin, xmax, ymax
         bbox            =  c(5328683, 3250589, 5330000, 3260000), # xmin, ymin, xmax, ymax
         spatmeth        = 'bbox',
         out_format      = 'GTiff',
         compress        = 'LZW',
         #out_projsel     = 'User Defined',
         #output_proj     = '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',
         delete_hdf      = TRUE,
         parallel        = TRUE
)

crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 6000000, 1000000, 5500000) # swap to base raster later
blank_raster <- rast(crs=crs, extent=euro_ext, res = 1000)




# top line is if have multiple files
# so could look at multiple years and get mode rather than picking a year
landcover_pos_files <- list.files('data/landcover_data/LandCover_Type_Yearly_500m_v6/LC1', '201[4-9]|202[0-9]', full.names = T)
landcover_pos_c <- rast(landcover_pos_files)
land_pos <- modal(landcover_pos_c)
names(land_pos) <- 'land_pos'
rm(landcover_pos_c)

land_pos



plot(land_pos, main='Land cover 2020', axes=T)




# bit below here is from https://rspatialdata.github.io/land_cover.html

# Reading in the downloaded landcover raster data
IGBP_raster <- raster(here::here("data/temp/france/LandCover_Type_YEarly_500m_v6/LC1/MCD12Q1_LC1_2019_001.tif"))

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

