## Looking at the WAHIS data. 
## also filtering on whether overlap with FAO

rm(list = ls())

library(tidyverse)
library(readxl)

wahis <- read_excel("data/flu_data/infur_20230609.xlsx", sheet = 2)

# contains all diseases and species so filter these down to the ones we want

colnames(wahis)
unique(wahis$disease_eng)
disease_inc <- as.data.frame(table(wahis$disease_eng))

disease_inc[63,"Var1"]
disease_inc[71,"Var1"]
disease_inc[72,"Var1"]

# want to check out all three of these categories

inf_data <- wahis[which(wahis$disease_eng %in% 
                          c("High pathogenicity avian influenza viruses (poultry) (Inf. with)",
                           "Influenza A virus (Inf. with)",
                           "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)")),]

table(inf_data$region)
table(is.na(inf_data$region)) # no missing data for this variable
# filter to Europe only

inf_data_eur <- filter(inf_data, region == "Europe")

table(inf_data_eur$Terra_Aqua) # all terrestrial. The only options for this are terrestrial or aqua 
table(inf_data_eur$Epi_unit)
dz_unit_df <- as.data.frame(table(inf_data_eur$disease_eng, inf_data_eur$Epi_unit))
