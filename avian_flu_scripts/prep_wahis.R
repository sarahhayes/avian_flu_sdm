## Looking at the WAHIS data. 
## also filtering on whether overlap with FAO

rm(list = ls())

library(tidyverse)
library(readxl)

wahis <- read_excel("data/flu_data/raw_data/WOAH_july_2023.xlsx", sheet = 2)

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

# I think the data we want are the ones labelled
# "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)"


inf_a_hpai <- dplyr::filter(inf_data_eur, disease_eng == "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)" )

# not sure what the "Influenza A virus (infection with)" category refers to. 
# Have a little look at it
inf_a <- dplyr::filter(inf_data_eur, disease_eng == "Influenza A virus (Inf. with)")
table(inf_a$country)
table(inf_a$Epi_unit)
# They are farm cases. 
range(inf_a$`event_start date`)
# And all occur 2009 - 2010. 

table(inf_a_hpai$country)
table(inf_a_hpai$Epi_unit)
nrow(is.na(inf_a_hpai$`event_start date`))
# there are some zoo locations here. Remove these? 

inf_a_hpai_slim <- dplyr::select(inf_a_hpai, -c(strain_fr, strain_esp, sero_sub_genotype_fr,
                                            sero_sub_genotype_esp, disease_fr, disease_esp, 
                                            Report_Nat_Ref, Outbreak_Nat_Ref, level3_area_id,
                                            level3_unique_code, level3_name, level2_area_id, 
                                            level2_unique_code, level2_name, level1_area_id, 
                                            level1_unique_code, level1_name))

                
# write.csv(inf_a_hpai_slim, "data/flu_data/prepped_data/woah_europe.csv", row.names = F)             
