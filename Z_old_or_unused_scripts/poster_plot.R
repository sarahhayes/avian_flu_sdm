# 06/12/2022

### Looking at the data from the FAO empres-i data set.
### this was filtered to be confirmed cases of avian influenza in wild birds only. 

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(zoo)

empres_europe_pos <- read.csv("data/flu_data/empresi_wild_confirmed_europe.csv")
# empres_europe_pos_capt <- read.csv("data/flu_data/empresi_captive_confirmed_europe.csv")

# empres_europe_pos <- rbind(empres_europe_pos, empres_europe_pos_capt)

# we don't need to humans affected of human deaths columns. 
empres_europe_pos <- dplyr::select(empres_europe_pos, !c(Humans.affected, Human.deaths))

# convert the dates of observed and reported dates into date format and extract year
empres_europe_pos$obs_date <- strptime(empres_europe_pos$Observation.date..dd.mm.yyyy.,
                                       format = "%d/%m/%Y")
empres_europe_pos$obs_year <- format(empres_europe_pos$obs_date, format = "%Y")
empres_europe_pos$obs_month <- format(empres_europe_pos$obs_date, format = "%m")

empres_europe_pos$rep_date <- strptime(empres_europe_pos$Report.date..dd.mm.yyyy.,
                                       format = "%d/%m/%Y")
empres_europe_pos$rep_year <- format(empres_europe_pos$rep_date, format = "%Y")
empres_europe_pos$rep_month <- format(empres_europe_pos$rep_date, format = "%m")

empres_europe_pos$obs_month_year <- format(empres_europe_pos$obs_date, format = "%m/%Y")

obs_counts <- empres_europe_pos %>% group_by(obs_month_year) %>% count()
obs_counts$date <- as.Date(zoo::as.yearmon(obs_counts$obs_month_year, "%m/%Y"))



# Bring in bv-brc data 

pos_data <- read.csv("data/flu_data/BVBRC_wild_captivewild_positive.csv") %>%
  filter(Host.Natural.State == "Wild")

columns_needed <- c("Collection.Date", "Collection.Year", "Collection.Country", 
                    "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                    "Host.Species", "Host.Group", "Host.Natural.State", "Host.Health", 
                    "Subtype", "Strain")

pos_data <- dplyr::select(pos_data, all_of(columns_needed))

euro_countries <- c("Bulgaria", "Georgia", "Germany", "Greece", "Greenland", "Hungary", "Iceland", 
                    "Latvia", "Lithuania","Netherlands", "Romania", "Russia", "Sweden", "Turkey")

pos_data_europe <-filter(pos_data, Collection.Country %in% euro_countries)

pos_data_europe$date <- strptime(pos_data_europe$Collection.Date,
                                 format = "%Y-%m-%d")

pos_data_europe$month_year <- format(pos_data_europe$date, format = "%m/%Y")
pos_data_europe$month <- format(pos_data_europe$date, format = "%m")

pos_counts <- pos_data_europe %>% group_by(month_year) %>% count()
pos_counts$date <- as.Date(zoo::as.yearmon(pos_counts$month_year, "%m/%Y"))


# combine the two data sets

emp <- obs_counts[,c("n", "date")]
bv <- pos_counts[,c("n", "date")]


full_dat <- rbind(emp, bv)

sums <- full_dat %>% group_by(date) %>%
  summarise(total = sum(n))
  

ggplot(sums, aes(x = date, y = total)) + 
  geom_line( col = "red", size = 1) + 
  labs(title = "", x = "Date", y = "Count" ) + 
  theme_bw(base_size = 20, base_family = "bold") + 
  xlim(c(as.Date("2003-01-01"), as.Date("2023-01-01")))

#ggsave("plots/timeline_empres_fao.png")


png("plots/joint_time_series.png", width = 560, height = 490)
ggplot(sums, aes(x = date, y = total)) + 
  geom_line( col = "red", size = 1.0) + 
  labs(title = "", x = "Date", y = "Number of cases" ) + 
  theme_bw(base_size = 20, base_family = "bold") + 
  xlim(c(as.Date("2005-01-01"), as.Date("2023-01-01")))
dev.off()
