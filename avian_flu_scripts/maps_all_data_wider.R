### Looking at the distribution of cases outside of Europe but 
### within bordering countries

## First bring in all the data sets

rm(list = ls())

library(tidyverse)
library(readxl)

# WAHIS
wahis <- read_excel("data/flu_data/raw_data/WOAH_july_2023.xlsx", sheet = 2)

# contains all diseases and species so filter these down to the ones we want

inf_data <- wahis[which(wahis$disease_eng %in% 
                          c("High pathogenicity avian influenza viruses (poultry) (Inf. with)",
                            "Influenza A virus (Inf. with)",
                            "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")),]

inf_a_hpai <- dplyr::filter(inf_data, disease_eng == "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)" )

#FAO 
# use the full data as want to include asia as well 
fao_data  <- read.csv("data/flu_data/raw_data/epidemiology-raw-data_202308011345.csv")

table(fao_data$Region)
## filter to Europe and Asia

fao_data <- dplyr::filter(fao_data, Region %in% c("Europe", "Asia"))
colnames(fao_data)

fao_data <-  dplyr::rename(fao_data, observation.date = "Observation.date..dd.mm.yyyy.")
fao_data <-  dplyr::rename(fao_data, report.date = "Report.date..dd.mm.yyyy.")

# There are rows with no entries for observation date so remove these
fao_data <- fao_data[which(fao_data$observation.date != ""),]

# we only want to include wild birds
table(fao_data$Species)


# BVBRC
bvbrc_data <- read.csv("data/flu_data/raw_data/BVBRC_surveillance.csv")

columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                    "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                    "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health", 
                    "Subtype", "Strain")

bvbrc_data_slim <- dplyr::select(bvbrc_data, all_of(columns_needed))

bvbrc_data_slim <- filter(bvbrc_data_slim, Host.Natural.State == "Wild")

# remove the ones without a year
bvbrc_data_slim <- drop_na(bvbrc_data_slim, Collection.Year)

# There are samples labelled "env" which I assume are environmental so remove these
bvbrc_data_slim <- bvbrc_data_slim[which(bvbrc_data_slim$Host.Species != "Env"),]

# separate positive and negative
pos_data <- filter(bvbrc_data_slim, Pathogen.Test.Result == "Positive")
neg_data <- filter(bvbrc_data_slim, Pathogen.Test.Result == "Negative")


# Bring in the maps

library(terra)
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
plot(zipmap)
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 10000000, 1000000, 10000000) 

# change projection and extent. 
# using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs) 
euro_map_crop <- terra::crop(euro_map, euro_ext)
plot(euro_map)
plot(euro_map_crop)

# Transform the data to spatial points
pts_woah <- terra::vect(inf_a_hpai, geom=c("Longitude", "Latitude"),
                       crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
pts_woah <- terra::project(pts_woah,  "epsg:3035")

##
pts_fao <- terra::vect(fao_data, geom=c("Longitude", "Latitude"),
                        crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
pts_fao <- terra::project(pts_fao,  "epsg:3035")

## 
pts_bvbrc_pos <- terra::vect(pos_data, geom=c("Collection.Longitude", "Collection.Latitude"),
                        crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
pts_bvbrc_pos <- terra::project(pts_bvbrc_pos,  "epsg:3035")

pts_bvbrc_neg <- terra::vect(neg_data, geom=c("Collection.Longitude", "Collection.Latitude"),
                             crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
pts_bvbrc_neg <- terra::project(pts_bvbrc_pos,  "epsg:3035")


#pdf("plots/alldata_pos_map.pdf", width = 5, height = 6)
plot(euro_map_crop, col = "white", background = "azure2", main = "All positive")
plot(pts_bvbrc_pos, add = T, col = "red", pch = 3)
plot(pts_fao, add= T, col = 'blue', pch = 3)
plot(pts_woah, add = T, col = "orange", pch = 3)
abline( v = 7000000, h = 6400000)
abline(h = 1550000, v = 2600000)
dev.off()
