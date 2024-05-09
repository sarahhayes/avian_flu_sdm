## 07/05/2024
## Creating the avian flu database of domestic animal cases - focus on mammals in this script
## Plan is to overlay this on our risk map 

rm(list = ls())

library(readxl)
library(tidyverse)
library(terra)
library(RColorBrewer)

## bring in the different data sets and combine to one large data set. 

## first do fao 
## 
## FAO

# when downloading these data, go to the epidemiology page, select the virus, 
# select the date (01/01/9170 - 30/06/2023), select domestic and then download raw data from tab in upper right corner
# use the full data as want to include asia as well 
# I have selected those just listed as domestic and not included 'captive'. I did look at captive on the fao page but these were more
# rare/zoo species. First will focus on just domestic outbreaks. 

#fao_data  <- read.csv("data/flu_data/raw_data/fao_domestic_cases_raw_data.csv")
fao_data  <- read.csv("data/flu_data/raw_data/epidemiology-raw-data_domestic_apr_2024.csv") # later download from april 2024

table(fao_data$Region)

colnames(fao_data)
table(fao_data$Diagnosis.status)
table(fao_data$Animal.type)

# the rogue entry is misaligned and is in cats in Korea so remove
fao_data <- fao_data %>% dplyr::filter(Animal.type != "Confirmed")

# rename the date columns
fao_data <-  dplyr::rename(fao_data, any_of(c(observation.date = "Observation.date..dd.mm.yyyy.",
                                              report.date = "Report.date..dd.mm.yyyy.")))

# There are 1684 rows with no entries for observation date 
nrow(fao_data[which(fao_data$observation.date == ""),])

no_obs_data <- fao_data[which(fao_data$observation.date == ""),]
table(no_obs_data$Region)
# Very few of these are in europe so will remove them

fao_data <- fao_data[which(fao_data$observation.date != ""),]

# we only want to include domestic species this time
table(fao_data$Species) # selected just those listed as 'domestic' from website. So these should all be domestic. 

# Select the entries that are mammals 

mammal_sp_fao <- c("Canine", "Cats", "Cattle", "Dogs", "Donkey", "Ferret", "Goats", "House Crow", 
                     "Swine")

fao_data <- fao_data %>% filter(Species %in% mammal_sp_fao)
table(fao_data$Species)

# Rename some of them with very long names

fao_data_trim <- dplyr::select(fao_data, all_of(c("Latitude", "Longitude", "observation.date", "Species", "Country", "Serotype")))
fao_data_trim$Epi_unit <- NA # for later merge with WAHIS
fao_data_trim$source <- "fao"
fao_data_trim$observation.date <- as.Date(fao_data_trim$observation.date, "%d/%m/%Y")


####### Now do WAHIS

# WAHIS
wahis <- read.csv("data/flu_data/raw_data/WOAH_april_2024.csv") # later download from april 2024

# contains all diseases and species so filter these down to the ones we want

inf_data <- wahis[which(wahis$disease_eng %in% 
                          c("High pathogenicity avian influenza viruses (poultry) (Inf. with)",
                            "Influenza A virus (Inf. with)",
                            "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")),]


## first select those labelled as wild = F

inf_a_hpai <- dplyr::filter(inf_data, is_wild == F)
table(inf_a_hpai$disease_eng)

poss_cases <- dplyr::filter(inf_a_hpai, disease_eng == "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")
table(poss_cases$Epi_unit)
table(poss_cases$Species)

wahis_mamms_vect <- c("Bovine", "Cats", "Dogs", "Goats", "Mustelidae (dom)")

# select the mammals
inf_a_hpai <- dplyr::filter(inf_a_hpai, Species %in% wahis_mamms_vect)
table(inf_a_hpai$region) # all regions represented
table(inf_a_hpai$disease_eng)

inf_a_hpai_trim <- dplyr::select(inf_a_hpai, all_of(c("Latitude", "Longitude", "Outbreak_start_date", "Species", 
                                                      "country", "sero_sub_genotype_eng", "Epi_unit")))
inf_a_hpai_trim$source <- "woah"

# We now have the trimmed data sets for all the sources. 
# first need all the names to match so we can rbind

colnames(fao_data_trim)
colnames(inf_a_hpai_trim)
inf_a_hpai_trim <- rename(inf_a_hpai_trim, observation.date = Outbreak_start_date, 
                          Country = country, Serotype = sero_sub_genotype_eng)



all_ai_data <- rbind(fao_data_trim, inf_a_hpai_trim)

## check for any NA values
sum(is.na(all_ai_data$Latitude))
sum(is.na(all_ai_data$Longitude))

# None without location data

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
proj_crs <- 3035
point_data <- sf::st_transform(point_data, proj_crs)
#plot(point_data)
point_data # Now a spatial points object in the right projection
class(point_data)
head(point_data)

point_data_df <- point_data %>%
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2])
head(point_data_df)
st_geometry(point_data)
range(point_data_df$X)
range(point_data_df$Y) # checking that we have labelled the columns correctly

## Now select only those rows which are in the area that we want

euro_shp <- terra::vect("output/euro_map.shp")
ext(euro_shp)

ai_data_prj_area <- point_data_df %>%
  dplyr::filter(X >= 2600000 & X <= 7000000) %>%
  dplyr::filter(Y >= 1500000 & Y <= 6400000)


table(ai_data_prj_area$Species)

# There are only 15 entries for domestic mammals. 
# Not really enough to follow this through 