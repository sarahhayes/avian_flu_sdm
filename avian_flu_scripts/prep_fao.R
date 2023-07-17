### Look at the fao empres-i data.
### These are filtered to be influenza-avian in wild hosts
### Set date as 01/01/1970 - 30/06/2023
### Need to look for mammals in data. 

rm(list = ls())
library(tidyverse)
library(zoo)

# This is wild and captive wild

#fao_data <- read.csv("data/flu_data/raw_data/FAO Latest Reported Events.csv")
fao_data <- read.csv("data/flu_data/raw_data/overview-raw-data_202307141702.csv")

colnames(fao_data)
table(fao_data$Region)

## filter to Europe

#data_europe <- fao_data[which(fao_data$Region == "Europe"),]
data_europe <- fao_data
colnames(data_europe)

data_europe <-  dplyr::rename(data_europe, observation.date = "Observation.date..dd.mm.yyyy.")
data_europe <-  dplyr::rename(data_europe, report.date = "Report.date..dd.mm.yyyy.")

# There are rows with no entries for observation date so remove these

sum(is.na(data_europe$observation.date))
nrow(data_europe[which(data_europe$observation.date == ""),])

# remove the ones without a date
data_europe <- data_europe[which(data_europe$observation.date != ""),]

# Save this file 

#write.csv(data_europe, "data/flu_data/prepped_data/fao_europe.csv", row.names = F)



