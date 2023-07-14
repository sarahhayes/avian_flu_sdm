### Look at the bv-brc data.
### These are filtered to be avian hosts with influenza a
### All years and all areas. 
### Downloaded positive and negative data so will need to separate. 
### didn't include the 'not tested' samples

rm(list = ls())
library(tidyverse)
library(zoo)

# This is wild and captive wild
bvbrc_data <- read.csv("data/flu_data/raw_data/BVBRC_surveillance.csv")

columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                    "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                    "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health", 
                    "Subtype", "Strain")

bvbrc_data_slim <- dplyr::select(bvbrc_data, all_of(columns_needed))

table(bvbrc_data_slim$Collection.Country)
## identify the countries that are in Europe

euro_countries <- c("Bulgaria", "Georgia", "Germany", "Greece", "Greenland", "Hungary", "Iceland", 
                    "Latvia", "Lithuania","Netherlands", "Romania", "Russia", "Sweden", "Turkey")

data_europe <-filter(bvbrc_data_slim, Collection.Country %in% euro_countries)

table(data_europe$Host.Natural.State) # select wild only
data_europe <- filter(data_europe, Host.Natural.State == "Wild")

table(data_europe$Pathogen.Test.Result) # separate positive and negative
pos_data_europe <- filter(data_europe, Pathogen.Test.Result == "Positive")
neg_data_europe <- filter(data_europe, Pathogen.Test.Result == "Negative")

pos_data_europe %>%
  group_by(Collection.Country) %>%
  count()

neg_data_europe %>%
  group_by(Collection.Country) %>%
  count()


## there are two wrongly entered collection dates in the negative data.
## Opt to remove 
## Doing step-by-step to ensure removing just one at a time 

nd <- neg_data_europe[which(neg_data_europe$Collection.Date != "0002-10-11T00:00:00Z"),] 
nd <- nd[which(nd$Collection.Date != "0005-12-09T00:00:00Z"),]


# Save these files 

#write.csv(pos_data_europe, "data/flu_data/prepped_data/bvbrc_pos_europe.csv")
#write.csv(nd, "data/flu_data/prepped_data/bvbrc_neg_europe.csv")



