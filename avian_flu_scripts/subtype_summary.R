# Subtype compile

library(tidyverse)
library(magrittr)

bdf <- read.csv("data_offline\\Avian flu data\\Subtypes\\bvbrc_subtypes.csv") %>% 
  set_colnames(c("Subtype", "Freq")) %>% 
  mutate(Class = case_when(grepl("H7|H5", Subtype) ~ "Potentially HPAI",
                           Subtype == "" ~ "Missing",
                           TRUE ~ "LPAI"))

fdf <- read.csv("data_offline\\Avian flu data\\Subtypes\\fao_subtypes.csv") %>% 
  separate(Var1, into = c("Subtype", "Class")) %>% 
  set_colnames(c("Subtype", "Class", "Freq"))
fdf[17,1:2] <- c(NA, "HPAI")
fdf[18,1:2] <- c(NA, "LPAI")

wdf <- read.csv("data_offline\\Avian flu data\\Subtypes\\woah_subtypes.csv") %>% 
  set_colnames(c("Subtype", "Freq")) %>% 
  mutate(Class = "HPAI")

bdf %>% group_by(Class) %>% summarise(n = sum(Freq)) %>% mutate(p = n*100/sum(n))
fdf %>% group_by(Class) %>% summarise(n = sum(Freq)) %>% mutate(p = n*100/sum(n))
wdf %>% group_by(Class) %>% summarise(n = sum(Freq)) %>% mutate(p = n*100/sum(n))

bind_rows(bdf, fdf, wdf) %>% select(Subtype) %>% distinct %>% filter(!(grepl("x|X|y|Y|n", Subtype))) %>% filter(!is.na(Subtype))
