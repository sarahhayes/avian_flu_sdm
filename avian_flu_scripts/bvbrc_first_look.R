### Quick look at the bv-brc data.
### These are filtered to include wild and captive-wild birds only
### All years and all areas. 
### Downloaded positive and negative data separately. 
### didn't include the 'not tested' samples

rm(list = ls())
library(tidyverse)
library(zoo)

# This is wild and captive wild
pos_data <- read.csv("data/flu_data/older_downloads/BVBRC_wild_captivewild_positive.csv")
neg_data <- read.csv("data/flu_data/older_downloads/BVBRC_wild_captivewild_negative.csv")
# if we just want to use wild data. 
pos_data <- read.csv("data/flu_data/older_downloads/BVBRC_wild_captivewild_positive.csv") %>%
  filter(Host.Natural.State == "Wild")
neg_data <- read.csv("data/flu_data/older_downloads/BVBRC_wild_captivewild_negative.csv")%>%
  filter(Host.Natural.State == "Wild")


columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                     "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                     "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health", 
                    "Subtype", "Strain")

pos_data <- dplyr::select(pos_data, all_of(columns_needed))
neg_data <- dplyr::select(neg_data, all_of(columns_needed))

table(pos_data$Collection.Country)
table(neg_data$Collection.Country)

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

bvbrc_subtype <- as.data.frame(table(pos_data_europe$Subtype))

bvbrc_subtype_country <- as.data.frame(table(pos_data_europe$Collection.Country, 
                                             pos_data_europe$Subtype)) %>%
  filter(Freq != 0) %>%
  filter(Var2 != "") %>%
  filter(!Var2 %in% c("HxNx", "HxNy", "HXNX"))

#write.csv(bvbrc_subtype, "output/bvbrc_subtypes.csv", row.names = F)
#write.csv(bvbrc_subtype_country, "output/bvbrc_subtypes_countries.csv", row.names = F)

# counts per country
pos_data_europe_counts <- pos_data_europe %>%
  dplyr::count(as.factor(Collection.Country),as.factor(Collection.Year),.drop=FALSE) %>%
  rename(Country = "as.factor(Collection.Country)", Year = "as.factor(Collection.Year)")

neg_data_europe_counts <- neg_data_europe %>%
  count(as.factor(Collection.Country),as.factor(Collection.Year),.drop=FALSE) %>%
  rename(Country = "as.factor(Collection.Country)", Year = "as.factor(Collection.Year)")


ggplot(pos_data_europe_counts, aes(x = Year, y = n, 
       fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "Positive samples from wild birds in Europe")


ggplot(neg_data_europe_counts, aes(x = Year, y = n, 
                                   fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "Negative samples from wild birds in Europe")


## just by country and not by year
pos_data_europe %>%
  dplyr::count(as.factor(Collection.Country), .drop=FALSE) %>%
  rename(Country = "as.factor(Collection.Country)") %>%
  ggplot(., aes(x = Country, y = n, fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", 
       title = "Positive samples from wild birds in Europe 
       (BV-BRC)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
#ggsave("plots/by_country_pos_bvbrc.png" )#, width = 5, height = 5)

neg_data_europe %>%
  dplyr::count(as.factor(Collection.Country), .drop=FALSE) %>%
  rename(Country = "as.factor(Collection.Country)") %>%
  ggplot(., aes(x = Country, y = n, fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "Negative samples from wild  
       birds in Europe (BV-BRC)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
# ggsave("plots/by_country_neg_bvbrc.png")


## plots of time series of the data
## plot a time series of the cases. 
pos_data_europe$date <- strptime(pos_data_europe$Collection.Date,
                                       format = "%Y-%m-%d")

pos_data_europe$month_year <- format(pos_data_europe$date, format = "%m/%Y")
pos_data_europe$month <- format(pos_data_europe$date, format = "%m")

pos_counts <- pos_data_europe %>% group_by(month_year) %>% count()
pos_counts$date <- as.Date(zoo::as.yearmon(pos_counts$month_year, "%m/%Y"))

ggplot(pos_counts, aes(x = date, y = n)) + 
  geom_line( col = "red", size = 1) + 
  labs(title = "Positive cases in wild birds in Europe 
       (BV-BRC)", 
       x = "Date", y = "Count") + 
  theme_bw() + 
  xlim(c(as.Date("2003-01-01"), as.Date("2023-01-01")))

#ggsave("plots/timeline_bvbrc_pos.png")

## look at entries per month
month_counts <- as.data.frame(table(pos_data_europe$month))
month_counts$Var1 <- as.character(month_counts$Var1)
ggplot(data = month_counts, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", col = "orange", fill = "orange") + 
  labs(x = "Month", title = "Monthly positive counts bvbrc data") 


# now negative
neg_data_europe$date <- strptime(neg_data_europe$Collection.Date,
                                 format = "%Y-%m-%d")

neg_data_europe$month_year <- format(neg_data_europe$date, format = "%m/%Y")
neg_data_europe$month <- format(neg_data_europe$date, format = "%m")


neg_counts <- neg_data_europe %>% group_by(month_year) %>% count()
neg_counts$date <- as.Date(zoo::as.yearmon(neg_counts$month_year, "%m/%Y"))

# There are two dates that are wrong
str(neg_data_europe)

ggplot(neg_counts, aes(x = date, y = n)) + 
  geom_line( col = "red", size = 1) + 
  labs(title = "Negative cases in wild birds in Europe 
       (BV-BRC)", 
       x = "Date", y = "Count") + 
  theme_bw() + 
  xlim(c(as.Date("2003-01-01"), as.Date("2023-01-01")))

#ggsave("plots/timeline_bvbrc_neg.png")

# plot of monthly counts of negative cases
month_counts_neg <- as.data.frame(table(neg_data_europe$month))
month_counts_neg$Var1 <- as.character(month_counts_neg$Var1)
ggplot(data = month_counts_neg, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", col = "orangered3", fill = "orangered3", ) + 
  labs(x = "Month", title = "Monthly negative counts bvbrc data")


### Look at the species that are included in the data 

species_table <- as.data.frame(table(pos_data_europe$Host.Species))
#write.csv(species_table, "output/species_bvbrc_pos.csv")


 
