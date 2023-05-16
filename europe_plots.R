# Script for plots of cases on map of Europe. 
# First draft - Jan 2023. 
# The areas/countries that are visualised can be tweaked if we want to include others 
# This is likely to be determined by the final data we use

rm(list = ls())
library(tidyverse)
library(zoo)
library(utils)
library(rgdal)
library(ggplot2)
library(raster)

# Read in the data we want to plot
# first the bvbrc data

bvbrc_pos_data <- read.csv("data/flu_data/BVBRC_wild_captivewild_positive.csv")
bvbrc_neg_data <- read.csv("data/flu_data/BVBRC_wild_captivewild_negative.csv")

columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                    "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                    "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health")

pos_data <- dplyr::select(bvbrc_pos_data, all_of(columns_needed))
neg_data <- dplyr::select(bvbrc_neg_data, all_of(columns_needed))

# select just the wild species
pos_data <- dplyr::filter(pos_data, Host.Natural.State == "Wild" )
neg_data <- dplyr::filter(neg_data, Host.Natural.State == "Wild")

euro_countries <- c("Bulgaria", "Georgia", "Germany", "Greece", "Greenland", "Hungary", "Iceland", 
                    "Latvia", "Lithuania","Netherlands", "Romania", "Russia", "Sweden", "Turkey")

pos_data_europe <-filter(pos_data, Collection.Country %in% euro_countries)
neg_data_europe <- filter(neg_data, Collection.Country %in% euro_countries)

pos_data_europe %>%
  group_by(Collection.Country) %>%
  count()

## there are two wrongly entered collection dates in the negative data.
## Will initially change to what I think but might be better to remove? 

neg_data_europe[which(neg_data_europe$Collection.Date == "0002-10-11T00:00:00Z"),
                "Collection.Date"] <- "2002-10-11T00:00:00Z"

neg_data_europe[which(neg_data_europe$Collection.Date == "0005-12-09T00:00:00Z"),
                "Collection.Date"] <- "2005-12-09T00:00:00Z"


### And now read in the FAO data
### Just using the wild birds for now - not captive
empres_europe_pos <- read.csv("data/flu_data/empresi_wild_confirmed_europe.csv")

# we don't need to humans affected or human deaths columns. 
empres_europe_pos <- dplyr::select(empres_europe_pos, !c(Humans.affected, Human.deaths))

# convert the dates of observed and reported dates into date format and extract year
empres_europe_pos$obs_date <- strptime(empres_europe_pos$Observation.date..dd.mm.yyyy.,
                                       format = "%d/%m/%Y")
empres_europe_pos$obs_year <- format(empres_europe_pos$obs_date, format = "%Y")

empres_europe_pos$rep_date <- strptime(empres_europe_pos$Report.date..dd.mm.yyyy.,
                                       format = "%d/%m/%Y")
empres_europe_pos$rep_year <- format(empres_europe_pos$rep_date, format = "%Y")

# how many NA values in each column of the data
colSums(is.na(empres_europe_pos))

# No missing data in any of the columns except for the observation date. 
# Let's see how much difference there is between the dates

empres_europe_pos$obs_rep_time_diff <- abs(empres_europe_pos$obs_date - empres_europe_pos$rep_date)
units(empres_europe_pos$obs_rep_time_diff) <- "days"
empres_europe_pos$obs_rep_time_diff <- round(empres_europe_pos$obs_rep_time_diff)

### Prepare the maps

### Numerous options of maps to plot. Have looked at others in plotting_europe.R

zipmap <- readOGR(dsn = "C:/Users/hayes/OneDrive - Nexus365/Documents/GitHub/avian_flu_sdm/data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                  layer = "CNTR_RG_03M_2020_4326")


euro_countries <- c("AL", "AD", "AM", "AT", "AZ", "BY", "BE", "BA", "BG", 
                    "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "GE", "DE", 
                    "GR", "HU", "IC", "IE", "IR", "IT", "KZ", "XK", "LV", "LI", 
                    "LT", "LU", "MT", "MD", "MC", "ME", "NL", "MK", "NO", 
                    "PL", "PT", "RO", "RU", "SM", "RS", "SK", "SI", "ES", 
                    "SE", "CH", "UA", "UK", "VA", "TR", "IS") # adding iceland for now

euro_map <- zipmap[zipmap$CNTR_ID %in% euro_countries, ]
plot(euro_map)

# Transform the data to spatial points
bvbrc_xy_neg <- neg_data_europe[,c("Collection.Longitude", "Collection.Latitude")]
bvbrc_spdf_neg <- SpatialPointsDataFrame(coords = bvbrc_xy_neg, data = neg_data_europe, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

bvbrc_xy_pos <- pos_data_europe[,c("Collection.Longitude", "Collection.Latitude")]
bvbrc_spdf_pos <- SpatialPointsDataFrame(coords = bvbrc_xy_pos, data = pos_data_europe, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

empres_xy_pos <- empres_europe_pos[,c("Longitude", "Latitude")]
empres_spdf_pos <- SpatialPointsDataFrame(coords = empres_xy_pos, data = empres_europe_pos,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))


### 
crop_euro_map <- crop(euro_map, extent(-25, 50, 34, 72))
plot(crop_euro_map, col = "white", bg = "azure2")
plot(bvbrc_spdf_neg, add = T, pch = 20)
plot(bvbrc_spdf_pos, add = T, col = "red", pch = 18)

## There is a gap in the east that is coloured as water just because i've removed these
## these countries, but it isn't actually water.
## This can be changed once we have confirmed area

## too many if we also overlay FAO so plot separately
plot(crop_euro_map, col = "white", bg = "azure2")
plot(empres_spdf_pos, add = T, col = "blue", pch = 22)

par(mar = c(1,1,1,1))
# plot all positive points from both data sets
plot(crop_euro_map, col = "white", border = "dark grey", bg = "azure2", )
plot(empres_spdf_pos, add = T, col = "red", pch = 20, cex = 0.5)
plot(bvbrc_spdf_pos, add = T, col = "red", pch = 20, cex = 0.5)
# plot(bvbrc_spdf_neg, add = T, col = "orange", pch = 18, cex = 0.5)

pdf("plots/allpos.pdf", height=10, width=10)
plot(crop_euro_map, col = "white", border = "dark grey", bg = "azure2", )
plot(empres_spdf_pos, add = T, col = "red", pch = 20, cex = 0.5)
plot(bvbrc_spdf_pos, add = T, col = "red", pch = 20, cex = 0.5)
dev.off()
