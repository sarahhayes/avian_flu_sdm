# Taking the codes that are left for birds once filter for Europe on ebird webpage
# This is 708 species
# https://science.ebird.org/en/status-and-trends/species?regionCode=eu

rm(list = ls())

library(readxl)

sp_codes <- read.csv("ebird/ebird_species_europe_copy.csv", header = F)

remove_rows <- sort(c(seq(from = 1, to = nrow(sp_codes), by = 4), seq(from = 3, to = nrow(sp_codes), by = 4)))

keep_rows <- sort(c(seq(from = 2, to = nrow(sp_codes), by = 4), seq(from = 4, to = nrow(sp_codes), by = 4)))

sp_codes_keep <- as.data.frame(sp_codes[keep_rows,])
colnames(sp_codes_keep) <- "vars"

name_rows <- seq(from = 1, to = nrow(sp_codes_keep), by = 2)

library(dplyr)
library(tidyr)

df <- sp_codes_keep %>%
  mutate(ind = rep(c(1, 2),length.out = n())) %>%
  group_by(ind) %>%
  mutate(id = row_number()) %>%
  pivot_wider(names_from = ind, values_from = vars) %>%
  dplyr::select(-id)

colnames(df) <- c("species", "code")

df2 <- df %>% 
  mutate(across(where(is.character), stringr::str_trim))

df3 <- df2
df3$code <- df3$code %>% 
  stringr::str_replace( "<https://science.ebird.org/en/status-and-trends/species/", "") %>%
  stringr::str_replace(">", "")

#write.csv(df3, "ebird/codes_for_europe_clean.csv")
