### Look at the prepped fao data to produce some summaries

rm(list = ls())
library(tidyverse)
library(zoo)

# This is wild species only
woah_data <- read.csv("data/flu_data/prepped_data/woah_europe.csv")

# look at a few things to check as expected
table(woah_data$country)

# Doesn't contain any species data

## Subtypes from the positive cases

woah_subtype <- as.data.frame(table(woah_data$sero_sub_genotype_eng))

sum(woah_subtype$Freq) # This tells us that we don't have strain for all of them 
nrow(woah_data) - sum(woah_subtype$Freq) # 25 that don't have strain/subtype.

woah_subtype_country <- as.data.frame(table(woah_data$country, 
                                             woah_data$sero_sub_genotype_eng)) %>%
  filter(Freq != 0) 

#write.csv(woah_subtype, "output/woah_summaries/fao_subtypes.csv", row.names = F)
#write.csv(woah_subtype_country, "output/woah_summaries/fao_subtypes_countries.csv", row.names = F)

# counts per country

# make a year column. 
head(woah_data$Outbreak_start_date)
woah_data$obs_year <- as.Date(woah_data$Outbreak_start_date, "%Y-%m-%d")
head(woah_data$obs_year)
woah_data$obs_year <- as.numeric(format(woah_data$obs_year, "%Y"))
head(woah_data$obs_year)

woah_data_counts <- woah_data %>%
  dplyr::count(as.factor(country),as.factor(obs_year),.drop=FALSE) %>%
  rename(Country = "as.factor(country)", Year = "as.factor(obs_year)")

ggplot(woah_data_counts, aes(x = Year, y = n, 
                            fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "WOAH samples from Europe")


## bit too busy so split just by country and not by year
woah_data %>%
  dplyr::count(as.factor(country), .drop=FALSE) %>%
  rename(Country = "as.factor(country)") %>%
  ggplot(., aes(x = Country, y = n, fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", 
       title = "Positive samples from wild birds in Europe 
       (FAO)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
#ggsave("plots/by_country_woah.png" )#, width = 5, height = 5)


## plots of time series of the data
## plot a time series of the cases. 
woah_data$date <- strptime(woah_data$Outbreak_start_date,
                          format = "%Y-%m-%d")
head(woah_data$date)
range(woah_data$date)
class(woah_data$date)

#woah_data$month_year <- format(woah_data$date, format = "%m-%Y")
woah_data$month_year <- as.Date(zoo::as.yearmon(woah_data$date, "%Y-%m"))
head(woah_data$month_year)
range(woah_data$month_year)

woah_counts <- woah_data %>% group_by(month_year) %>% count()
class(woah_counts$month_year)
range(woah_counts$month_year)

ggplot(woah_counts, aes(x = month_year, y = n)) + 
  geom_line( col = "red", linewidth = 1) + 
  labs(title = "Positive cases in Europe 
       (WOAH)", 
       x = "Date", y = "Count") + 
  theme_bw() + 
  scale_x_date(date_breaks = "2 years",
               date_minor_breaks = "1 year", date_labels = "%Y") # + 

#ggsave("plots/timeline_woah_pos.png")

class(woah_data$date)
## look at entries per month
woah_data$month <- lubridate::month(woah_data$date)
head(woah_data$month)
range(woah_data$month)

month_counts <- as.data.frame(table(woah_data$month))
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
pts_pos <- terra::vect(woah_data, geom=c("Longitude", "Latitude"),
                       crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")
pts_pos <- terra::project(pts_pos,  "epsg:3035")


#pdf("plots/woah_pos_map.pdf", width = 5, height = 6)
plot(euro_map_crop, col = "white", background = "azure2", main = "WOAH positive")
plot(pts_pos, add = T, col = "red", pch = 18)
abline( v = 8000000)
#dev.off()

### Need to look at the species that are included in the data 

# Will need to remove the mammals
# and try and find a way to identify the species