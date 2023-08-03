## Creating the avian flu data csv for use within the BART model

rm(list = ls())

library(readxl)
library(tidyverse)

## bring in the different data sets and combine to one large data set. 

## first do fao as it's the biggest
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

unwanted_sp <- c("Unspecified Mammal", "Stone Marten", 
                 "South American Coati (Nasua Nasua):Procyonidae-Carnivora",
                 "Red Fox", "Polecat", "Polar Fox", "Nyctereutes Viverrinus (Japanese Racoon Dog)",
                 "Mink", "Harbor Seal (Phoca Vitulina):Phocidae-Carnivora", 
                 "Gray Seal (Halichoerus Grypus):Phocidae-Carnivora", "Fox",
                 "European Pine Marten", "Civet", "Caspian Seal (Pusa Caspica)")

fao_data <- fao_data[-which(fao_data$Species %in% unwanted_sp),]

fao_data_trim <- dplyr::select(fao_data, all_of(c("Latitude", "Longitude", "observation.date")))
fao_data_trim$source <- "fao"
fao_data_trim$observation.date <- as.Date(fao_data_trim$observation.date, "%d/%m/%Y")

# WAHIS
wahis <- read_excel("data/flu_data/raw_data/WOAH_july_2023.xlsx", sheet = 2)

# contains all diseases and species so filter these down to the ones we want

inf_data <- wahis[which(wahis$disease_eng %in% 
                          c("High pathogenicity avian influenza viruses (poultry) (Inf. with)",
                            "Influenza A virus (Inf. with)",
                            "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")),]

inf_a_hpai <- dplyr::filter(inf_data, disease_eng == "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)" )
table(inf_a_hpai$region) # all regions represented

inf_a_hpai_trim <- dplyr::select(inf_a_hpai, all_of(c("Latitude", "Longitude", "Outbreak_start_date")))
inf_a_hpai_trim$source <- "woah"


# BVBRC
bvbrc_data <- read.csv("data/flu_data/raw_data/BVBRC_surveillance.csv")

columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                    "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                    "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health", 
                    "Subtype", "Strain")

bvbrc_data_slim <- dplyr::select(bvbrc_data, all_of(columns_needed))
table(bvbrc_data_slim$Collection.Country) # still contains global data as required

bvbrc_data_slim <- filter(bvbrc_data_slim, Host.Natural.State == "Wild")

# remove the ones without a year
bvbrc_data_slim <- drop_na(bvbrc_data_slim, Collection.Year)

# There are samples labelled "env" which I assume are environmental so remove these
bvbrc_data_slim <- bvbrc_data_slim[which(bvbrc_data_slim$Host.Species != "Env"),]

# separate positive and negative
pos_data <- filter(bvbrc_data_slim, Pathogen.Test.Result == "Positive")
neg_data <- filter(bvbrc_data_slim, Pathogen.Test.Result == "Negative")

bv_pos_data_trim <- dplyr::select(pos_data, all_of(c("Collection.Date", 
                                                     "Collection.Latitude",
                                                     "Collection.Longitude")))

bv_pos_data_trim$source <- "bvbrc"
bv_pos_data_trim$Collection.Date <- as.Date(bv_pos_data_trim$Collection.Date,
                                            "%Y-%m-%d")


bv_neg_data_trim <- dplyr::select(neg_data, all_of(c("Collection.Date", 
                                                     "Collection.Latitude",
                                                     "Collection.Longitude")))
bv_neg_data_trim$source <- "bvbrc"
bv_neg_data_trim$Collection.Date <- as.Date(bv_neg_data_trim$Collection.Date,
                                            "%Y-%m-%d")

# We now have the trimmed data sets for all the sources. 
# combine them and add a section for pos or neg

# first need all the names to match so we can rbind

colnames(fao_data_trim)
colnames(inf_a_hpai_trim)
colnames(bv_pos_data_trim)
inf_a_hpai_trim <- rename(inf_a_hpai_trim, observation.date = Outbreak_start_date)
bv_pos_data_trim <- rename(bv_pos_data_trim, Latitude = Collection.Latitude,
                           Longitude = Collection.Longitude,
                           observation.date = Collection.Date)

bv_pos_data_trim <- dplyr::select(bv_pos_data_trim, c(Latitude, Longitude, observation.date, source))
colnames(bv_pos_data_trim)

all_ai_data <- rbind(fao_data_trim, inf_a_hpai_trim, bv_pos_data_trim)
all_ai_data$flu <- 1 # denotes a positive

## deal with the negatives

bv_neg_data_trim <- rename(bv_neg_data_trim, Latitude = Collection.Latitude,
                           Longitude = Collection.Longitude,
                           observation.date = Collection.Date)

bv_neg_data_trim <- dplyr::select(bv_neg_data_trim, c(Latitude, Longitude, observation.date, source))
colnames(bv_neg_data_trim)
bv_neg_data_trim$flu <- 0 # denotes negative test

all_ai_data <- rbind(all_ai_data, bv_neg_data_trim)

## check for any NA values
sum(is.na(all_ai_data$Latitude))
sum(is.na(all_ai_data$Longitude))

# remove the rows that are NA for location 
all_ai_data <- drop_na(all_ai_data, Latitude)
sum(is.na(all_ai_data$Latitude))

###################################################################################
## make the points into a raster. 
## Using info from here: 
## # https://gis.stackexchange.com/questions/458522/use-r-to-create-raster-from-data-frame-with-non-uniformly-gridded-points-and-cat

library(sf)
# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

point_data <- st_as_sf(x = all_ai_data, 
                        coords = c("Longitude", "Latitude"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
class(point_data)
# change the projections
point_data <- sf::st_transform(point_data, 3035)
#plot(point_data)
point_data
# Now a spatial points object in the right projection

points_sp <- sf::as_Spatial(point_data)

points_vect <- terra::vect(points_sp)

#read in the raster we are using for the project
euro_rast <- terra::rast("output/euro_rast.tif")
euro_rast

points_rast <- terra::rasterize(points_vect, euro_rast,"flu")
plot(points_rast)
points_rast
head(points_rast)

#read in the map
euro_map <- terra::vect("output/euro_map.shp")

plot(euro_map)
plot(points_rast, add = TRUE)

table(values(points_rast))

## Hard to see. Look at with larger raster

#read in the raster we are using for the project
euro_rast_10k <- terra::rast("output/euro_rast_10k.tif")
euro_rast_10k

points_rast_10k <- terra::rasterize(points_vect, euro_rast_10k,"flu")
plot(points_rast_10k)
points_rast_10k
head(points_rast_10k)
table(values(points_rast_10k))
# Easier to visualise that it probably is working OK. 
# I think the 1km are just too small to see. 

#Perhaps just look at UK? 
GB_ext <- terra::ext(3300000, 3800000, 3100000, 4100000)
GB_crop <- terra::crop(euro_map, GB_ext)
plot(GB_crop)

GB_rast <- terra::crop(points_rast, GB_ext)
plot(GB_crop)
plot(GB_rast, add = T, axes = F)

###
# Create the csv file for avian flu cases
# make a points object using the centre of each pixel from the ref raster
points_euro_rast <- terra::as.points(euro_rast)
points_euro_rast

tictoc::tic()
flu_res <- terra::extract(points_rast, points_euro_rast, method = "simple", xy = T)
tictoc::toc()

table(flu_res$last)
table(values(points_euro_rast))
table(values(points_rast))
euro_rast
nrow(flu_res)
points_rast

# why have we got fewer values in the results table rather than the raster? 