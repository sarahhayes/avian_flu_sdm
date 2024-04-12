## 011/04/2024
## Preparing data for the Americas 

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

# use this if filtering by specified region 
fao_data_trim_americas <- fao_data_trim %>% filter(Region == "Americas")

# WAHIS
wahis <- read_excel("data/flu_data/raw_data/WOAH_july_2023.xlsx", sheet = 2)

# contains all diseases and species so filter these down to the ones we want

inf_data <- wahis[which(wahis$disease_eng %in% 
                          c("High pathogenicity avian influenza viruses (poultry) (Inf. with)",
                            "Influenza A virus (Inf. with)",
                            "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")),]


inf_a_hpai <- dplyr::filter(inf_data, disease_eng %in% c("Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)",
                                                         "High pathogenicity avian influenza viruses (poultry) (Inf. with)"))

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
#inf_a_hpai_wild <- inf_a_hpai_wild[which(inf_a_hpai_wild$wild_type != "captive"),] #10,328 observations as removes NAs
inf_a_hpai_wild <- dplyr::filter(inf_a_hpai_wild, (wild_type != "captive"|is.na(wild_type))) #14,449 observations

table(inf_a_hpai_wild$Epi_unit, inf_a_hpai_wild$is_wild)
zoo_cases_wild <- inf_a_hpai_wild[which(inf_a_hpai_wild$Epi_unit == "Zoo"),]

table(inf_a_hpai_wild$wild_type)

colnames(inf_a_hpai_wild)
inf_a_hpai_trim <- dplyr::select(inf_a_hpai_wild, all_of(c("Latitude", "Longitude", "Outbreak_start_date", "Species", 
                                                           "country", "sero_sub_genotype_eng", "region")))
inf_a_hpai_trim$source <- "woah"

inf_a_hpai_americas <- inf_a_hpai_trim %>% filter(region == "Americas")

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
# bvbrc_data_slim <- bvbrc_data_slim[which(bvbrc_data_slim$Host.Species != "Env"),]
bvbrc_data_slim <- dplyr::filter(bvbrc_data_slim, Host.Species != "Env")

# separate positive and negative
pos_data <- filter(bvbrc_data_slim, Pathogen.Test.Result == "Positive")

bv_pos_data_trim <- dplyr::select(pos_data, all_of(c("Collection.Date", 
                                                     "Collection.Latitude",
                                                     "Collection.Longitude",
                                                     "Host.Species",
                                                     "Collection.Country",
                                                     "Subtype")))

bv_pos_data_trim$source <- "zbvbrc"
bv_pos_data_trim$Collection.Date <- as.Date(bv_pos_data_trim$Collection.Date,
                                            "%Y-%m-%d")

americas_countries <- c("Argentina", "Brazil", "Canada", "Chile", "Colombia", "Guatemala", "Greenland",
                        "United States Of America", "USA" )


bv_americas <- bv_pos_data_trim %>% filter(Collection.Country %in% c(americas_countries)) 
bv_americas$Region <- "Americas"

table(bv_americas$Region, bv_americas$Collection.Country)


## we only want HPAI 
table(bv_americas$Subtype)
# filter those that contain H5 or H7 - although these are not necessarily HPAI
poss_hpai <- bv_americas %>% filter(str_detect(Subtype,"H5") | str_detect(Subtype, "H7"))

table(poss_hpai$Subtype)

# We now have the trimmed data sets for all the sources. 
# first need all the names to match so we can rbind

colnames(fao_data_trim_americas)
colnames(inf_a_hpai_americas)
colnames(poss_hpai)
inf_a_hpai_americas <- rename(inf_a_hpai_americas, observation.date = Outbreak_start_date, 
                                   Country = country, Serotype = sero_sub_genotype_eng, Region = region)
poss_hpai <- rename(poss_hpai, Latitude = Collection.Latitude,
                    Longitude = Collection.Longitude,
                    observation.date = Collection.Date,
                    Species  = Host.Species,
                    Country = Collection.Country,
                    Serotype = Subtype)


poss_hpai <- dplyr::select(poss_hpai, c(Latitude, Longitude, observation.date, Species,
                                        Country, Serotype, source, Region))

americas_data <- rbind(fao_data_trim_americas, inf_a_hpai_americas, poss_hpai)


## check for any NA values
sum(is.na(americas_data$Latitude))
sum(is.na(americas_data$Longitude))

# remove the rows that are NA for location 
americas_data <- drop_na(americas_data, Latitude)
sum(is.na(americas_data$Latitude))

## Want to remove mammals

table(americas_data$Species)

species_list <- read.csv("avian_flu_scripts/detecting_mammals_asia_americas.csv")

unwanted_sp <- species_list[which(species_list$class == "Mammalia"), "species"]

## manually extract those not automatically done. 
unwanted_sp_manual_ext <- c("Al indeterminatum fau",
                            "Bottlenose Dolphin (Tursiops Truncatus)", 
                            "Chilean Dolphin (Cephalorhynchus Eutropia)",
                            "Civet", 
                            "Coati (Gen. Nasua)",
                            "Common Raccoon (Procyon Lotor)",
                            "Cougar (Puma Concolor)",
                            "Environmental",
                            "Fox",
                            "HieraaÃ«tus fasciatus",
                            "Kodiak Grizzly Bear (Ursus Arctos Horribilis)",
                            "Nyctereutes Viverrinus (Japanese Racoon Dog)",
                            "Skunk (Mephitis Mephitis)",
                            "Unspecified Mammal")


unwanted_sp <- c(unwanted_sp, unwanted_sp_manual_ext)

## Remove from the data so just birds

ai_pos_birds <- 
  americas_data[-which(americas_data$Species %in% unwanted_sp),]

## Ideally we want only the HPAI

table(ai_pos_birds$Serotype)
# remove those which state LPAI
ai_pos_birds <- ai_pos_birds %>% filter(!str_detect(Serotype, "LPAI"))
table(ai_pos_birds$Serotype)
table(ai_pos_birds$source)
table(ai_pos_birds$Serotype, ai_pos_birds$source)

# Produce a new column 
ai_pos_birds$serotype_HN <- ai_pos_birds$Serotype
# Now remove extra info
ai_pos_birds$serotype_HN <-  gsub("HPAI","",as.character(ai_pos_birds$serotype_HN))
ai_pos_birds$serotype_HN<-toupper(ai_pos_birds$serotype_HN) # so they match
ai_pos_birds$serotype_HN <- trimws(ai_pos_birds$serotype_HN) # trim any white space in the entries to aid matching
table(ai_pos_birds$serotype_HN)


library(sf)
# change into sf object so can change projection of point. 
# Assuming FAO are entered in standard lon/lat format

point_data <- st_as_sf(x = ai_pos_birds, 
                       coords = c("Longitude", "Latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
class(point_data)

# change the projections
proj_crs <- 8858
point_data <- sf::st_transform(point_data, proj_crs)
point_data # Now a spatial points object in the right projection
class(point_data)
head(point_data)

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

library(data.table)
dt_pos <- data.table(point_data_df)

table(dt_pos$Country) # a couple of different spellings for a number of the countries
dt_pos[which(dt_pos$Country == "United States of America"), "Country"] <- "USA"
dt_pos[which(dt_pos$Country == "United States Of America"), "Country"] <- "USA"
table(dt_pos$Country)


no_duplicates_date_and_loc <- unique(dt_pos, by = c("X", "Y", "observation.date"))
no_duplicates_all <- unique(dt_pos, by = c("X", "Y", "observation.date", "source", "Species", "Country"))
no_dups_all_except_source <- unique(dt_pos, by = c("X", "Y", "observation.date", "Species", "Country"))
no_dups_date_loc_serotype <- unique(dt_pos, by = c("X", "Y", "observation.date", "serotype_HN"))

# The middle two are the same, suggesting that there are none that are the same for all details except source. 
# The final one (no_dups_date_loc_serotype) contains 89 more than the one filtered just on date and location.
# On inspection, these are indeed listed as different serotypes in the same location and same date. could be an error,
# but no way to check
## No duplicates_all will still have some that are the same place and time but different species. However, when looking
## at these, there are instances where we have fao and woah reporting a same location and date but one has the common name
## for the bird and one has the scientific name. I think we basically want to know if a location has AI in wild birds
## so we don't mind which species. 
## As such, no_dups_date_loc_serotype is probably the one to use. 

## ideally, if we are removing duplicates, we want to remove the bv-brc ones as these are the ones we are less confident in. 
## As we have labelled bvbrc with a z, these should be removed as they will be the last row. 


no_dups_date_loc_serotype <- sf::st_as_sf(no_dups_date_loc_serotype)
# split by the two regions

americas_ai_birds <- no_dups_date_loc_serotype

table(americas_ai_birds$source)

## make plots to show distribution

americas_map <- terra::vect("output/am_vect_poly.shp")
plot(americas_map)

#png("plots/americas_pos_data_hpai.png", width = 480, height = 480)
par(mar = c(0,0,0,0))
plot(americas_map)
plot(americas_ai_birds[which(americas_ai_birds$source == "fao"),], add = T, pch = 18, col = "red", cex = 0.6)
plot(americas_ai_birds[which(americas_ai_birds$source == "woah"),], add = T, pch = 18, col = "blue", cex = 0.6)
plot(americas_ai_birds[which(americas_ai_birds$source == "zbvbrc"),], add = T, pch = 18, col = "green", cex = 0.6)
legend(-4000000,20, legend=c("FAO", "WOAH", "BV-BRC"),  
       fill = c("red", "blue", "green"),
       cex = 0.8
)
#dev.off()


## Plot the timeseries
## Need to ensure that have all the weeks represented as currently misses out those with no cases 

## To make sure we include all weeks in a plot, even those without any data, we need a dataframe with 
## all the weeks and months,

range(americas_ai_birds$observation.date)

## start on 5th Jan 1987 as that's a Monday
week_calendar_americas <- data.frame(date=seq(as.Date("1987-01-05"), as.Date("2023-07-27"), by="day")) %>% 
  mutate(week_num=isoweek(date),year=year(date)) %>%
  group_by(year,week_num) %>% 
  summarise(weekdate=min(date)) 
week_calendar_americas$week_of_study <- seq(1, nrow(week_calendar_americas), by = 1)

# combine with the data
dates_americas <- merge(americas_ai_birds, week_calendar_americas)

# plot directly
ggplot(dates_americas, aes(x=weekdate)) + geom_histogram(bins = nrow(week_calendar_americas), col = "blue") + 
  labs(y = "Number of weekly cases", x = "Year", size = 18) +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey90"))
#ggsave("plots/hist_by_week_hpai_only_americas.png")


# could also make a dataframe of the counts of the number in each week
count_dates_americas <- dates_americas %>% 
  count(weekdate)

weekly_counts_americas <- left_join(week_calendar_americas, count_dates_americas)
weekly_counts_americas[which(is.na(weekly_counts_americas$n)), "n"] <- 0


## Also then can have as a line plot
ggplot(data = weekly_counts_americas, aes(x = weekdate, y = n)) +
  geom_line(lwd = 0.9, col = "blue")  + 
  labs(y = "Number of weekly cases", x = "Year", size = 18) +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "grey90"))
#ggsave("plots/lineplot_by_week_hpai_only_americas.png")

# Now to plot the time series by serotype 
serotype_data_americas <- as.data.frame(table(americas_ai_birds$serotype_HN))

# For quite a number of the subtypes there are only a few entries. 
# For the line plot by subtype, only use those that have >10 entries
serotype_data_ten_or_more_americas <- serotype_data_americas[which(serotype_data_americas$Freq > 10),]
serotype_data_ten_or_more_americas <- serotype_data_ten_or_more_americas[which(serotype_data_ten_or_more_americas$Var1 != ""),]

serotype_over_ten_americas <- serotype_data_ten_or_more_americas$Var1 %>% droplevels()

subtype_plot_data_americas <- dates_americas[which(dates_americas$serotype_HN %in% serotype_over_ten_americas),]
table(subtype_plot_data_americas$serotype_HN)

range(subtype_plot_data_americas$year)


yearlabs_americas <- c( "1988", "1989", "1990", "1991", "1992", "1993", "1994", "1995",
                        "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004",
                        "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", 
                        "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")

weekbreaks_americas <- week_calendar_americas[which(week_calendar_americas$week_num == 1), "week_of_study"]
weekbreaks_americas <- weekbreaks_americas[["week_of_study"]]

ggplot(subtype_plot_data_americas, aes(x = week_of_study, fill = serotype_HN)) + 
  geom_bar(position = "stack") + 
  labs(y = "Number of weekly cases", x = "Year", size = 18, fill = "Serotype") +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        panel.background = element_rect(fill = "grey90", colour = "grey90"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "white")) +
  scale_x_continuous(breaks = weekbreaks_americas, labels = yearlabs_americas)
# ggsave("plots/subtypes/serotype_by_week_hpai_only_americas.png")

# This might need looking at again depending on what we wish to focus on. Perhaps just look at certain part of plot? 
# Below is a small adaptation just to look at 2005 onwards 

ggplot(subtype_plot_data_americas[which(subtype_plot_data_americas$year >=2005),], aes(x = week_of_study, fill = serotype_HN)) + 
  geom_bar(position = "stack") + 
  labs(y = "Number of weekly cases", x = "Year", size = 18, fill = "Serotype") +
  theme(axis.text.x = element_text(angle=90, margin = margin(t = 0.1, r = 0.2, b = 0.2, l = 0.3, unit = "cm"), 
                                   face = "bold", size = 10, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 10),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 12),
        panel.background = element_rect(fill = "grey90", colour = "grey90"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',
                                        colour = "white")) +
  scale_x_continuous(breaks = weekbreaks_americas, labels = yearlabs_americas)
ggsave("plots/subtypes/serotype_by_week_hpai_only_americas_2005onwards.png")





