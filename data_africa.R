## 21/03/2024
## Exploring data availability across Africa

rm(list = ls())

library(readxl)
library(tidyverse)
library(terra)
library(RColorBrewer)

## bring in the different data sets and combine to one large data set. 

## first do fao as it's the biggest
#FAO 
# use the full data as want to include asia as well 
fao_data  <- read.csv("data/flu_data/raw_data/epidemiology-raw-data_202308011345.csv")

table(fao_data$Region)

colnames(fao_data)
# rename the date columns
fao_data <-  dplyr::rename(fao_data, observation.date = "Observation.date..dd.mm.yyyy.")
fao_data <-  dplyr::rename(fao_data, report.date = "Report.date..dd.mm.yyyy.")

# There are rows with no entries for observation date so remove these
fao_data <- fao_data[which(fao_data$observation.date != ""),]

# we only want to include wild birds
table(fao_data$Species)


fao_data_trim <- dplyr::select(fao_data, all_of(c("Latitude", "Longitude", "observation.date", "Species", "Country", "Serotype", "Region")))
fao_data_trim$source <- "fao"
fao_data_trim$observation.date <- as.Date(fao_data_trim$observation.date, "%d/%m/%Y")

fao_data_trim_africa <- fao_data_trim %>% filter(Region %in% c("Africa"))

# WAHIS
wahis <- read_excel("data/flu_data/raw_data/WOAH_july_2023.xlsx", sheet = 2)

# contains all diseases and species so filter these down to the ones we want

inf_data <- wahis[which(wahis$disease_eng %in% 
                          c("High pathogenicity avian influenza viruses (poultry) (Inf. with)",
                            "Influenza A virus (Inf. with)",
                            "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")),]

inf_a_hpai <- dplyr::filter(inf_data, disease_eng == "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)" )
table(inf_a_hpai$region) # all regions represented
colnames(inf_a_hpai)
table(inf_a_hpai$is_wild)
table(inf_a_hpai$Epi_unit)
table(inf_a_hpai$Epi_unit, inf_a_hpai$is_wild)

# I think we want to filter by is.wild = T 
# This will leave in some animals that are listed as wild in a zoo - but these could have been wild animals 
# that were found within zoo grounds? Or they could be captive wild... Looking at the numbers, the latter seems more likely

zoo_cases <- inf_a_hpai[which(inf_a_hpai$Epi_unit == "Zoo"),]

table(zoo_cases$wild_type) # this doesn't show the many NAs

# so perhaps we should remove the zoo captives? So select all those that are wild and then remove those that
# are captive 
inf_a_hpai_wild <- inf_a_hpai[which(inf_a_hpai$is_wild == T),]
inf_a_hpai_wild <- inf_a_hpai_wild[which(inf_a_hpai_wild$wild_type != "captive"),]
table(inf_a_hpai_wild$Epi_unit, inf_a_hpai_wild$is_wild)

table(inf_a_hpai_wild$wild_type)

colnames(inf_a_hpai_wild)
inf_a_hpai_trim <- dplyr::select(inf_a_hpai_wild, all_of(c("Latitude", "Longitude", "Outbreak_start_date", "Species", 
                                                           "country", "sero_sub_genotype_eng", "region")))
inf_a_hpai_trim$source <- "woah"

inf_a_hpai_africa <- inf_a_hpai_trim %>% filter(region %in% c("Africa"))

# BVBRC
bvbrc_data <- read.csv("data/flu_data/raw_data/BVBRC_surveillance.csv")

table(bvbrc_data$Collection.Country)

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

bv_pos_data_trim <- dplyr::select(pos_data, all_of(c("Collection.Date", 
                                                     "Collection.Latitude",
                                                     "Collection.Longitude",
                                                     "Host.Species",
                                                     "Collection.Country",
                                                     "Subtype")))

bv_pos_data_trim$source <- "bvbrc"

bv_pos_data_trim$Collection.Date <- as.Date(bv_pos_data_trim$Collection.Date,
                                            "%Y-%m-%d")

africa_countries <- c("Cameroon", "Central African Republic", "Egypt", "Gabon", "Ghana", "Republic of the Congo", 
                      "South Africa")


bv_africa <- bv_pos_data_trim %>% filter(Collection.Country %in% africa_countries)

## we only want HPAI 
table(bv_africa$Subtype)
# filter those that contain H5 or H7 - although these are not necessarily HPAI
poss_hpai <- bv_africa %>% filter(str_detect(Subtype,"H5") | str_detect(Subtype, "H7"))
poss_hpai$Region <- "Africa"

# We now have the trimmed data sets for all the sources. 
# combine them and add a section for pos or neg

# first need all the names to match so we can rbind

colnames(fao_data_trim_africa)
colnames(inf_a_hpai_africa)
colnames(poss_hpai)
inf_a_hpai_africa <- rename(inf_a_hpai_africa, observation.date = Outbreak_start_date, 
                                   Country = country, Serotype = sero_sub_genotype_eng, Region = region)
poss_hpai <- rename(poss_hpai, Latitude = Collection.Latitude,
                    Longitude = Collection.Longitude,
                    observation.date = Collection.Date,
                    Species  = Host.Species,
                    Country = Collection.Country,
                    Serotype = Subtype)

poss_hpai <- dplyr::select(poss_hpai, c(Latitude, Longitude, observation.date, Species,
                                        Country, Serotype, source, Region))

africa_data <- rbind(fao_data_trim_africa, inf_a_hpai_africa, poss_hpai)


## check for any NA values
sum(is.na(africa_data$Latitude))
sum(is.na(africa_data$Longitude))

# remove the rows that are NA for location 
africa_data <- drop_na(africa_data, Latitude)
sum(is.na(africa_data$Latitude))

## Below is from the Europe code with a few obvious entries added - will need to be re-done if we do look at Asia and the Americas

unwanted_sp <- c("Unspecified Mammal", "Stone Marten", 
                 "South American Coati (Nasua Nasua):Procyonidae-Carnivora",
                 "Red Fox", "Polecat", "Polar Fox", "Nyctereutes Viverrinus (Japanese Racoon Dog)",
                 "Mink", "Harbor Seal (Phoca Vitulina):Phocidae-Carnivora", 
                 "Gray Seal (Halichoerus Grypus):Phocidae-Carnivora", "Fox",
                 "European Pine Marten", "Civet", "Caspian Seal (Pusa Caspica)",
                 "Halichoerus grypus", "Phoca vitulina", "Mustela putorius", 
                 "Ursus americanus", "Ursus arctos", "Virginia Opossum (Didelphis Virginiana)",
                 "Vulpes vulpes", "Wild Cat", "Wild Fox", "Procyon lotor", 
                 "Chilean Dolphin (Cephalorhynchus Eutropia)", "Al indeterminatum fau")

## Remove from the data so just birds

#ai_pos_birds <- 
#  africa_data[-which(africa_data$Species %in% unwanted_sp),]

ai_pos_birds <- africa_data %>% dplyr::filter(!Species %in% unwanted_sp)
  

## Ideally we want only the HPAI

table(ai_pos_birds$Serotype)
# remove those which state LPAI

ai_pos_birds <- ai_pos_birds %>% filter(!str_detect(Serotype, "LPAI"))

library(sf)
# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

point_data <- st_as_sf(x = ai_pos_birds, 
                       coords = c("Longitude", "Latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
class(point_data)

point_data_df <- point_data %>%
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2])
head(point_data_df)

# Add a year and month and year_month to the data

point_data_df$year <- lubridate::year(point_data_df$observation.date)
point_data_df$month <- lubridate::month(point_data_df$observation.date)
point_data_df$month_year <-  format(as.Date(point_data_df$observation.date), "%Y-%m")
point_data_df$week <- format(point_data_df$observation.date, "%Y Week %W")
point_data_df$week_num <- lubridate::isoweek(point_data_df$observation.date)

range(point_data_df$observation.date)

## make plots to show distribution
library(rnaturalearth)
spdf_world <- ne_download(scale = 110, type = "countries")
plot(spdf_world)

africa_map <- spdf_world %>% subset(., REGION_UN == "Africa")
plot(africa_map)

africa_ai_birds <- point_data

#png("plots/africa_pos_data.png", width = 480, height = 480)
par(mar = c(0,0,0,0))
plot(africa_map)
plot(africa_ai_birds[which(africa_ai_birds$source == "fao"),], add = T, pch = 18, col = "red", cex = 0.6)
plot(africa_ai_birds[which(africa_ai_birds$source == "woah"),], add = T, pch = 18, col = "blue", cex = 0.6)
plot(africa_ai_birds[which(africa_ai_birds$source == "bvbrc"),], add = T, pch = 18, col = "green", cex = 0.6)
legend(-40,20, legend=c("FAO", "WOAH", "BV-BRC"),  
       fill = c("red", "blue", "green"),
       cex = 0.8
)
#dev.off()



############################################################
# If we don't filter to remove the LPAI



## look at subtype. 
table(bv_africa$Subtype)



# We now have the trimmed data sets for all the sources. 
# combine them and add a section for pos or neg

# first need all the names to match so we can rbind

colnames(fao_data_trim_africa)
colnames(inf_a_hpai_africa)
colnames(bv_africa)

all_bv_africa  <- rename(bv_africa, Latitude = Collection.Latitude,
                    Longitude = Collection.Longitude,
                    observation.date = Collection.Date,
                    Species  = Host.Species,
                    Country = Collection.Country,
                    Serotype = Subtype)

all_bv_africa$Region <- "Africa"

all_bv_africa <- dplyr::select(all_bv_africa, c(Latitude, Longitude, observation.date, Species,
                                        Country, Serotype, source, Region))

africa_data_all <- rbind(fao_data_trim_africa, inf_a_hpai_africa, all_bv_africa)


## check for any NA values
sum(is.na(africa_data_all$Latitude))
sum(is.na(africa_data_all$Longitude))


## Below is from the Europe code with a few obvious entries added - will need to be re-done if we do look at Asia and the Americas

unwanted_sp <- c("Unspecified Mammal", "Stone Marten", 
                 "South American Coati (Nasua Nasua):Procyonidae-Carnivora",
                 "Red Fox", "Polecat", "Polar Fox", "Nyctereutes Viverrinus (Japanese Racoon Dog)",
                 "Mink", "Harbor Seal (Phoca Vitulina):Phocidae-Carnivora", 
                 "Gray Seal (Halichoerus Grypus):Phocidae-Carnivora", "Fox",
                 "European Pine Marten", "Civet", "Caspian Seal (Pusa Caspica)",
                 "Halichoerus grypus", "Phoca vitulina", "Mustela putorius", 
                 "Ursus americanus", "Ursus arctos", "Virginia Opossum (Didelphis Virginiana)",
                 "Vulpes vulpes", "Wild Cat", "Wild Fox", "Procyon lotor", 
                 "Chilean Dolphin (Cephalorhynchus Eutropia)", "Al indeterminatum fau")

## Remove from the data so just birds

ai_birds <- africa_data_all %>% dplyr::filter(!Species %in% unwanted_sp)


library(sf)
# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

all_point_data <- st_as_sf(x = ai_birds, 
                       coords = c("Longitude", "Latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


all_point_data_df <- all_point_data %>%
  dplyr::mutate(X = sf::st_coordinates(.)[,1],
                Y = sf::st_coordinates(.)[,2])


plot(africa_map)

#png("plots/africa_all_data.png", width = 480, height = 480)
par(mar = c(0,0,0,0))
plot(africa_map)
plot(all_point_data[which(all_point_data$source == "fao"),], add = T, pch = 18, col = "red", cex = 0.6)
plot(all_point_data[which(all_point_data$source == "woah"),], add = T, pch = 18, col = "blue", cex = 0.6)
plot(all_point_data[which(all_point_data$source == "bvbrc"),], add = T, pch = 18, col = "green", cex = 0.6)
legend(-40,20, legend=c("FAO", "WOAH", "BV-BRC"),  
       fill = c("red", "blue", "green"),
       cex = 0.8
)
#dev.off()


