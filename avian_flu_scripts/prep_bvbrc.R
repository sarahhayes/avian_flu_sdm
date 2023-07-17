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

# check if have dates for all the samples

sum(is.na(data_europe$Collection.Year))
# remove the ones without a year
data_europe <- drop_na(data_europe, Collection.Year)

table(data_europe$Host.Species)

# There are samples labelled "env" which I assume are environmental so remove these
data_europe <- data_europe[which(data_europe$Host.Species != "Env"),]

table(data_europe$Pathogen.Test.Result) # separate positive and negative
pos_data_europe <- filter(data_europe, Pathogen.Test.Result == "Positive")
neg_data_europe <- filter(data_europe, Pathogen.Test.Result == "Negative")

pos_data_europe %>%
  group_by(Collection.Country) %>%
  count()

neg_data_europe %>%
  group_by(Collection.Country) %>%
  count()


# Save these files 

#write.csv(pos_data_europe, "data/flu_data/prepped_data/bvbrc_pos_europe.csv", row.names = F)
#write.csv(neg_data_europe, "data/flu_data/prepped_data/bvbrc_neg_europe.csv", row.names = F)



