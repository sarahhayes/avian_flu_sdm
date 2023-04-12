## 19/01/2023

## GISAID

## Can we use the GISAID data to look at location? 
## First looking at small sample of data that was filtered to be 
## Europe only 01/01/2022 - 19/01/2023 and to be avian type a only 

rm(list = ls())

library(tidyverse)
library(readxl)

# all_gis_dat <- read_excel("data/flu_data/gisaid_epiflu_isolates_2022_2023.xls", sheet = 1)
all_gis_dat <- read_excel("data/flu_data/gisaid_epiflu_isolates_2010_2023.xls", sheet = 1)


# Filter the columns we want

# colnames(all_gis_dat)

gis_dat <- dplyr::select(all_gis_dat, all_of(c("Isolate_Name", "Subtype", "Location", 
                                            "Host", "Collection_Date")))

## We now want to split the location into all of its component parts - separated by / currently

## As don't know the maximum number of splits, make enough new columns so that it
## reports all rows filled with NA in last columns

g2 <- gis_dat %>%
  separate(Location, c("Loc_A","Loc_B", "Loc_C", "Loc_D", "Loc_E", "Loc_F", "Loc_G"), sep = "/")


## I think this would still require manual examination of the results to see which
## contained GPS. 
## Alternative might be to search for entries with GPS location. 

gis_gps <- gis_dat %>%
  mutate(Location = str_to_lower(Location)) %>% 
  filter(str_detect(Location, ("gps|zip"))) %>%
  separate(Location, c("Loc_A","Loc_B", "Loc_C", "Loc_D", "Loc_E", "Loc_F", "Loc_G"), sep = "/")

## alternative would be to select rows that contain numbers in location

gis_num <- gis_dat %>%
  mutate(Location = str_to_lower(Location)) %>% 
  filter(str_detect(Location, "[0-9]")) %>%
  separate(Location, c("Loc_A","Loc_B", "Loc_C", "Loc_D", "Loc_E", "Loc_F", "Loc_G"), sep = "/")

  
