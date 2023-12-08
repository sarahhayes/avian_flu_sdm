### Look at the fao empres-i data.
### These are filtered to be influenza-avian in wild hosts
### Set date as 01/01/1970 - 30/06/2023
### Need to look for mammals in data. 

rm(list = ls())
library(tidyverse)
library(zoo)

# This is wild only on empres-i website 
# the first file ("fao_world_wild_full...." comes from using the overview tab on the website
# and downloading separate csv files for certain years to avoid the row limit and then combining
# The third set "epidemiology_raw_data..." come from the epidemiology page and the download raw data
# tab in the upper right corner. Not 100% sure how I got the second set!!!)

fao_data1<- read.csv("data/flu_data/raw_data/fao_world_wild_full_1970_2023.csv")
fao_data2 <- read.csv("data/flu_data/raw_data/overview-raw-data_202307141702.csv")
fao_data3  <- read.csv("data/flu_data/raw_data/epidemiology-raw-data_202308011345.csv")

table(fao_data1$Region)
table(fao_data2$Region)
table(fao_data3$Region)

## filter to Europe

data_europe1 <- fao_data1[which(fao_data1$Region == "Europe"),]
data_europe2 <- fao_data2
data_europe3 <- fao_data3[which(fao_data3$Region == "Europe"),]
colnames(data_europe3)

data_europe2 <-  dplyr::rename(data_europe2, observation.date = "Observation.date..dd.mm.yyyy.")
data_europe2 <-  dplyr::rename(data_europe2, report.date = "Report.date..dd.mm.yyyy.")

data_europe3 <-  dplyr::rename(data_europe3, observation.date = "Observation.date..dd.mm.yyyy.")
data_europe3 <-  dplyr::rename(data_europe3, report.date = "Report.date..dd.mm.yyyy.")

# There are rows with no entries for observation date so remove these

# remove the ones without a date
data_europe1 <- data_europe1[which(data_europe1$observation.date != ""),]
data_europe2 <- data_europe2[which(data_europe2$observation.date != ""),]
data_europe3 <- data_europe3[which(data_europe3$observation.date != ""),]

data_europe1 <- data_europe1 %>%
  dplyr::filter(grepl("wild|Wild", Species))

# when filter for wild species, fewer data in the first data and these are labelled very unhelpfully
# in terms of species. 

data_europe2$date <- as.Date(data_europe2$observation.date, "%d/%m/%Y")
head(data_europe2)

data_europe3$date <- as.Date(data_europe3$observation.date, "%d/%m/%Y")
head(data_europe3)

data_europe2a <- data_europe2[which(data_europe2$date < "2023-07-01"),]
data_europe3a <- data_europe3[which(data_europe3$date < "2023-07-01"),]

geo1 <- data_europe2a[,c("Latitude", "Longitude", "observation.date", "Event.ID")]
geo2 <- data_europe3a[,c("Latitude", "Longitude", "observation.date", "Event.ID")]


diff_df <- setdiff(geo2, geo1)
diff_ids <- diff_df$Event.ID

examine_diffs <- data_europe3a[which(data_europe3a$Event.ID %in% diff_ids),]

# I think it might be the report dates that differ

data_europe2a$rep_date <- as.Date(data_europe2a$report.date, "%d/%m/%Y")
data_europe3a$rep_date <- as.Date(data_europe3a$report.date, "%d/%m/%Y")

range(data_europe2a$rep_date)
range(data_europe3a$rep_date)

diff_reps <- data_europe3a[which(data_europe3a$rep_date > "2023-07-13"),]

# still two missing! 

two_missing <- diff_reps$Event.ID
two_ids <- setdiff(diff_ids, two_missing)

ID2 <- data_europe3a[which(data_europe3a$Event.ID %in% two_ids),]

# Plan to use the data_europe3a dataset

# Save this file 

#write.csv(data_europe3a, "data/flu_data/prepped_data/fao_europe.csv", row.names = F)



