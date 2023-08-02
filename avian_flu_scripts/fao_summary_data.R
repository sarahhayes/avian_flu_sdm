### Look at the prepped fao data to produce some summaries

rm(list = ls())
library(tidyverse)
library(zoo)

# This is wild species only
fao_data <- read.csv("data/flu_data/prepped_data/fao_europe.csv")


# These data contain the source which is listed for some as OIE. Might be useful when comoaring
# these data to WAHIS to have this variable. 

# look at a few things to check as expected
table(fao_data$Country)
table(fao_data$Animal.type)
species_raw <- as.data.frame(table(fao_data$Species))

## Lots and lots of species and not entered in a very uniform way. 
### Can also see that there are some mammals in here. 

# From a manual look through this table we have:
# Fox; Gray seal; Harbor Seal; Mink; Red fox; South American Coati; Stone Marten; Unspecified mammal
# there are also lots listed as 'Domestic:...'

# i think we want to remove these

fao_data <- fao_data %>%
  dplyr::filter(grepl("wild|Wild", Species))


table(fao_data$Disease)
table(fao_data$Diagnosis.source)



# check if there are any NA values in the observation date. 

## Subtypes from the positive cases

fao_subtype <- as.data.frame(table(fao_data$Serotype))

fao_subtype_country <- as.data.frame(table(fao_data$Country, 
                                             fao_data$Serotype)) %>%
  filter(Freq != 0) 

# Something else we could look at is splitting these things by LPAI and HPAI

#write.csv(fao_subtype, "output/fao_summaries/fao_subtypes.csv", row.names = F)
#write.csv(fao_subtype_country, "output/fao_summaries/fao_subtypes_countries.csv", row.names = F)

# counts per country

# make a year column. 
fao_data$obs_year <- as.Date(fao_data$observation.date, "%d/%m/%Y")
head(fao_data$obs_year)
fao_data$obs_year <- as.numeric(format(fao_data$obs_year, "%Y"))
head(fao_data$obs_year)

fao_data_counts <- fao_data %>%
  dplyr::count(as.factor(Country),as.factor(obs_year),.drop=FALSE) %>%
  rename(Country = "as.factor(Country)", Year = "as.factor(obs_year)")

ggplot(fao_data_counts, aes(x = Year, y = n, 
                            fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "FAO samples from wild birds in Europe")


## bit too busy so split just by country and not by year
fao_data %>%
  dplyr::count(as.factor(Country), .drop=FALSE) %>%
  rename(Country = "as.factor(Country)") %>%
  ggplot(., aes(x = Country, y = n, fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", 
       title = "Positive samples from wild birds in Europe 
       (FAO)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
# ggsave("plots/by_country_pos_fao.png" )#, width = 5, height = 5)


## plots of time series of the data
## plot a time series of the cases. 

fao_data$date <- as.Date(fao_data$observation.date, "%d/%m/%Y")
head(fao_data$date)
fao_data$month_year <- format(fao_data$date, format = "%m/%Y")
fao_data$month <- format(fao_data$date, format = "%m")

fao_counts <- fao_data %>% group_by(month_year) %>% count()
head(fao_counts)
fao_counts$date <- as.Date(zoo::as.yearmon(fao_counts$month_year, "%m/%Y"))
range(fao_counts$date)

range(fao_counts$month_year)
ggplot(fao_counts, aes(x = date, y = n)) + 
  geom_line( col = "red", linewidth = 1) + 
  labs(title = "Positive cases in wild birds in Europe 
       (FAO)", 
       x = "Date", y = "Count") + 
  theme_bw() + 
  scale_x_date(date_breaks = "2 years",
               date_minor_breaks = "1 year", date_labels = "%Y") # + 

#ggsave("plots/timeline_fao_pos.png")

## look at entries per month
month_counts <- as.data.frame(table(fao_data$month))
month_counts$Var1 <- as.character(month_counts$Var1)
ggplot(data = month_counts, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", col = "orange", fill = "orange") + 
  labs(x = "Month", title = "Monthly positive counts fao data") 



# Next part is to plot the data
# Prepare the maps
# (Numerous options of maps to plot. Have looked at others in plotting_europe.R)


library(terra)
zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
plot(zipmap)
crs <- "epsg:3035"
euro_ext <- terra::ext(2000000, 9000000, 1000000, 9000000) 

# change projection and extent. 
# using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs) 
euro_map_crop <- terra::crop(euro_map, euro_ext)
plot(euro_map)
plot(euro_map_crop)

# Transform the data to spatial points
pts_pos <- terra::vect(fao_data, geom=c("Longitude", "Latitude"),
                       crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
pts_pos <- terra::project(pts_pos,  "epsg:3035")


#pdf("plots/fao_pos_map.pdf", width = 5, height = 6)
plot(euro_map_crop, col = "white", background = "azure2", main = "FAO positive")
plot(pts_pos, add = T, col = "red", pch = 18)
abline(v = 8000000)
#dev.off()

### Need to look at the species that are included in the data 

# Will need to remove the mammals
# and try and find a way to identify the species