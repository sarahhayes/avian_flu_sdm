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

# how many NA values in each column of the data
colSums(is.na(empres_europe_pos))

# look at the data available for the subtypes that are available 

empres_euro_subtypes <- as.data.frame(table(empres_europe_pos$Serotype)) %>%
  separate(Var1, c("subtype", "pathogenicity"))

empres_euro_subtypes[which(is.na(empres_euro_subtypes$pathogenicity)), "pathogenicity"] <-
  empres_euro_subtypes[which(is.na(empres_euro_subtypes$pathogenicity)), "subtype"]
  
# write.csv(empres_euro_subtypes, "output/empres_subtypes_summary.csv")

# No missing data in any of the columns except for the observation date. 
# Let's see how much difference there is between the dates

empres_europe_pos$obs_rep_time_diff <- abs(empres_europe_pos$obs_date - empres_europe_pos$rep_date)
units(empres_europe_pos$obs_rep_time_diff) <- "days"
empres_europe_pos$obs_rep_time_diff <- round(empres_europe_pos$obs_rep_time_diff)


## plot a time series of the cases. 
empres_europe_pos$obs_month_year <- format(empres_europe_pos$obs_date, format = "%m/%Y")

obs_counts <- empres_europe_pos %>% group_by(obs_month_year) %>% count()
obs_counts$date <- as.Date(zoo::as.yearmon(obs_counts$obs_month_year, "%m/%Y"))

ggplot(obs_counts, aes(x = date, y = n)) + 
  geom_line( col = "red", size = 1) + 
  labs(title = "Positive cases in wild birds in Europe 
       (Empres-i)", x = "Date", y = "Count") + 
  theme_bw() + 
  xlim(c(as.Date("2003-01-01"), as.Date("2023-01-01")))

#ggsave("plots/timeline_empres.png")


table(empres_europe_pos$Country)
table(empres_europe_pos$Diagnosis.status) # just checking they are all confirmed! 
sort(unique(empres_europe_pos$Country))

table(empres_europe_pos$Species)
## will need to filter for anything that isn't a bird as think there are a few mammals

pos_data_counts <- empres_europe_pos %>%
  dplyr::count(as.factor(Country),as.factor(obs_year),.drop=FALSE) %>%
  rename(Country = "as.factor(Country)", Year = "as.factor(obs_year)")

ggplot(pos_data_counts, aes(x = Year, y = n, 
                                   fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "Positive samples from wild and captive birds in Europe (empres-i)")

# Bit too busy to be of much use. 


# could group by subregion instead.
empres_europe_pos %>%
  dplyr::count(as.factor(Subregion),as.factor(obs_year),.drop=FALSE) %>%
  rename(Subregion = "as.factor(Subregion)", Year = "as.factor(obs_year)") %>%
ggplot(., aes(x = Year, y = n, 
                            fill = Subregion, colour = Subregion)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Year", y = "Frequency", title = "Positive samples from wild and captive birds in Europe (empres-i)")

# or just by country over whole period.
country_data <- empres_europe_pos %>%
  dplyr::count(as.factor(Country),.drop=FALSE) %>%
  rename(Country = "as.factor(Country)") %>%
  ggplot(., aes(x = Country, y = n, 
               fill = Country, colour = Country)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Country", y = "Frequency", title = "Positive samples from wild birds in Europe (empres-i)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
#ggsave("plots/by_country_empres.png")

country_data

## table of species
species_table <- as.data.frame(table(empres_europe_pos$Species))
# write.csv(species_table, "output/species_empresi_pos.csv")

# Look at how many entries per month

month_counts <- as.data.frame(table(empres_europe_pos$obs_month))
month_counts$Var1 <- as.character(month_counts$Var1)
ggplot(data = month_counts, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", col = "turquoise", fill = "turquoise") + 
  labs(x = "Month", title = "Monthly positive counts EMPRES-i data") 
