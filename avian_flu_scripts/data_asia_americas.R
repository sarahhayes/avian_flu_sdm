## 03/01/2024
## Exploring data availability across Asia and the Americas 

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

fao_data_trim_asia_americas <- fao_data_trim %>% filter(Region %in% c("Americas", "Asia"))

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

inf_a_hpai_americas_asia <- inf_a_hpai_trim %>% filter(region %in% c("Asia", "Americas"))

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

americas_countries <- c("Argentina", "Brazil", "Canada", "Chile", "Colombia", "Guatemala", 
                        "United States Of America", "USA" )

asia_countries <- c("Bangladesh", "Bhutan", "Cambodia", "China", "Georgia", "Japan", "Lebanon", "Mongolia", 
                    "Oman", "Russia", "Sri Lanka", "Taiwan", "Thailand", "Turkey", "Viet Nam")

bv_americas_asia <- bv_pos_data_trim %>% filter(Collection.Country %in% c(americas_countries, asia_countries)) 
bv_americas_asia$Region <- "Asia"
bv_americas_asia[which(bv_americas_asia$Collection.Country %in% americas_countries),"Region"] <- "Americas"

table(bv_americas_asia$Region)

## we only want HPAI 
table(bv_americas_asia$Serotype)
# filter those that contain H5 or H7 - although these are not necessarily HPAI
poss_hpai <- bv_americas_asia %>% filter(str_detect(Subtype,"H5") | str_detect(Subtype, "H7"))

# We now have the trimmed data sets for all the sources. 
# combine them and add a section for pos or neg

# first need all the names to match so we can rbind

colnames(fao_data_trim_asia_americas)
colnames(inf_a_hpai_americas_asia)
colnames(poss_hpai)
inf_a_hpai_americas_asia <- rename(inf_a_hpai_americas_asia, observation.date = Outbreak_start_date, 
                          Country = country, Serotype = sero_sub_genotype_eng, Region = region)
poss_hpai <- rename(poss_hpai, Latitude = Collection.Latitude,
                           Longitude = Collection.Longitude,
                           observation.date = Collection.Date,
                           Species  = Host.Species,
                           Country = Collection.Country,
                           Serotype = Subtype)

poss_hpai <- dplyr::select(poss_hpai, c(Latitude, Longitude, observation.date, Species,
                                                      Country, Serotype, source, Region))

americas_asia_data <- rbind(fao_data_trim_asia_americas, inf_a_hpai_americas_asia, poss_hpai)


## check for any NA values
sum(is.na(americas_asia_data$Latitude))
sum(is.na(americas_asia_data$Longitude))

# remove the rows that are NA for location 
americas_asia_data <- drop_na(americas_asia_data, Latitude)
sum(is.na(americas_asia_data$Latitude))

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
                 "Chilean Dolphin (Cephalorhynchus Eutropia)")

## Remove from the data so just birds

ai_pos_birds <- 
  americas_asia_data[-which(americas_asia_data$Species %in% unwanted_sp),]

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


# split by the two regions

asia_ai_birds <- point_data_df %>% filter(Region == "Asia")
americas_ai_birds <- point_data_df %>% filter(Region == "Americas")

## make plots to show distribution
library(rnaturalearth)
spdf_world <- ne_download(scale = 110, type = "countries")
plot(spdf_world)

americas_map <- spdf_world %>% subset(., REGION_UN == "Americas")
plot(americas_map)
asia_map <- spdf_world %>% subset(., REGION_UN == "Asia")
plot(asia_map)

#png("plots/americas_pos_data.png", width = 480, height = 480)
par(mar = c(0,0,0,0))
plot(americas_map)
plot(americas_ai_birds[which(americas_ai_birds$source == "fao"),], add = T, pch = 18, col = "red", cex = 0.6)
plot(americas_ai_birds[which(americas_ai_birds$source == "woah"),], add = T, pch = 18, col = "blue", cex = 0.6)
plot(americas_ai_birds[which(americas_ai_birds$source == "bvbrc"),], add = T, pch = 18, col = "green", cex = 0.6)
legend(-40,20, legend=c("FAO", "WOAH", "BV-BRC"),  
       fill = c("red", "blue", "green"),
       cex = 0.8
)
#dev.off()

#png("plots/asia_pos_data.png", width = 480, height = 480)
par(mar = c(0,0,0,0))
plot(asia_map)
plot(asia_ai_birds[which(asia_ai_birds$source == "fao"),], add = T, pch = 18, col = "red", cex = 0.6)
plot(asia_ai_birds[which(asia_ai_birds$source == "woah"),], add = T, pch = 18, col = "blue", cex = 0.6)
plot(asia_ai_birds[which(asia_ai_birds$source == "bvbrc"),], add = T, pch = 18, col = "green", cex = 0.6)
legend(40,0, legend=c("FAO", "WOAH", "BV-BRC"),  
       fill = c("red", "blue", "green"),
       cex = 0.8
)
#dev.off()

## Plot the timeseries
## Need to ensure that have all the weeks represented as currently misses out those with no cases 

## To make sure we include all weeks in a plot, even those without any data, we need a dataframe with 
## all the weeks and months,

range(americas_ai_birds$observation.date)

#Start on 5th Jan 1987 as that's a Monday and no data before middle of the year
start.plot_americas <- c("1987-01-05", "2023-12-31")
week_test_americas <- as.character(seq(as.Date(start.plot_americas[1]), as.Date(start.plot_americas[2]), by="weeks"))

all_dates_americas <- as.data.frame(matrix(nrow = length(week_test_americas), ncol = 1))
all_dates_americas[,1]  <- week_test_americas
colnames(all_dates_americas) <- "week_from_start"
all_dates_americas$date <- as.Date(all_dates_americas$week_from_start, "%Y-%m-%d")

all_dates_americas$year <- lubridate::year(all_dates_americas$date)
all_dates_americas$month <- lubridate::month(all_dates_americas$date)
all_dates_americas$month_year <- format(all_dates_americas$date, "%Y-%m")
all_dates_americas$week <- format(all_dates_americas$date, "%Y Week %W")
all_dates_americas$week_num <- lubridate::isoweek(all_dates_americas$date)

# by week 
subtype_counts_weekyear_americas <- americas_ai_birds %>%
  group_by(week) %>% count()

# merge with the weekly counts
all_dates_americas_data <- left_join(all_dates_americas, subtype_counts_weekyear_americas)
all_dates_americas_data$year_of_study <- all_dates_americas_data$year - 1987
all_dates_americas_data$week_seq <- all_dates_americas_data$week_num + (52*all_dates_americas_data$year_of_study)
all_dates_americas_data[which(is.na(all_dates_americas_data$n)),"n"] <- 0

weekbreaks_americas <- all_dates_americas_data %>% filter(grepl('Week 01', week))
weekbreaks_americas <- unique(weekbreaks_americas[,"week_seq"])
weekbreaks_americas


yearlabs_americas <- c("1987", "1988", "1989", "1990", "1991", "1992", "1993", "1994", "1995",
              "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004",
              "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", 
              "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")

par(mar = c(3,3,3,3))
ggplot(data = all_dates_americas_data, 
       aes(x = week_seq, y = n)) +
  geom_line(lwd = 0.5, col = "orange") + 
  scale_color_brewer(palette="Paired", na.translate = F) + 
  theme(axis.text.x = element_text(angle=90, margin = margin(t = .2, unit = "cm"), 
                                   face = "bold", size = 12, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 12),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 14),
        plot.background = element_rect(fill = "white"))+
  #  panel.background = element_rect(fill = "white", colour = "grey50")) + 
  labs(y = "Number of weekly cases", x = "Year", size = 18, title = "Americas") +
  scale_x_continuous(breaks = weekbreaks_americas, labels = yearlabs_americas)
ggsave("plots/time_series_americas.png")


# now repeat for asia

range(asia_ai_birds$observation.date)

#Start on 5th Jan 2004 as that's a Monday and no data before middle of the year
start.plot_asia <- c("2004-01-05", "2023-12-31")
week_test_asia <- as.character(seq(as.Date(start.plot_asia[1]), as.Date(start.plot_asia[2]), by="weeks"))

all_dates_asia <- as.data.frame(matrix(nrow = length(week_test_asia), ncol = 1))
all_dates_asia[,1]  <- week_test_asia
colnames(all_dates_asia) <- "week_from_start"
all_dates_asia$date <- as.Date(all_dates_asia$week_from_start, "%Y-%m-%d")

all_dates_asia$year <- lubridate::year(all_dates_asia$date)
all_dates_asia$month <- lubridate::month(all_dates_asia$date)
all_dates_asia$month_year <- format(all_dates_asia$date, "%Y-%m")
all_dates_asia$week <- format(all_dates_asia$date, "%Y Week %W")
all_dates_asia$week_num <- lubridate::isoweek(all_dates_asia$date)

# by week 
subtype_counts_weekyear_asia <- asia_ai_birds %>%
  group_by(week) %>% count()

# merge with the weekly counts
all_dates_asia_data <- left_join(all_dates_asia, subtype_counts_weekyear_asia)
all_dates_asia_data$year_of_study <- all_dates_asia_data$year - 2004
all_dates_asia_data$week_seq <- all_dates_asia_data$week_num + (52*all_dates_asia_data$year_of_study)
all_dates_asia_data[which(is.na(all_dates_asia_data$n)),"n"] <- 0

weekbreaks_asia <- all_dates_asia_data %>% filter(grepl('Week 01', week))
weekbreaks_asia <- unique(weekbreaks_asia[,"week_seq"])
weekbreaks_asia


yearlabs_asia <- c("2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", 
                       "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")

par(mar = c(3,3,3,3))
ggplot(data = all_dates_asia_data, 
       aes(x = week_seq, y = n)) +
  geom_line(lwd = 0.5, col = "turquoise") + 
  scale_color_brewer(palette="Paired", na.translate = F) + 
  theme(axis.text.x = element_text(angle=90, margin = margin(t = .2, unit = "cm"), 
                                   face = "bold", size = 12, vjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 12),
        plot.margin = margin(1,1,1,1, "cm"), 
        axis.title =  element_text(size = 14),
        plot.background = element_rect(fill = "white"))+
  #  panel.background = element_rect(fill = "white", colour = "grey50")) + 
  labs(y = "Number of weekly cases", x = "Year", size = 18, title = "Asia") +
  scale_x_continuous(breaks = weekbreaks_asia, labels = yearlabs_asia)
ggsave("plots/time_series_asia.png")

