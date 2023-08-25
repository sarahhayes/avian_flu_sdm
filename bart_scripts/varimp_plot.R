# Varimp plot

rm(list = ls())

library(embarcadero)

varimp_summ <- list()

for(i in 1:4){
  
  load(paste0("output/fitted-BART-models/sdm_Q",i,".rds"))
  
  varimp_raw <- sdm$varcount %>% 
    as.data.frame %>%
    mutate(./rowSums(.)) %>%
    mutate(./max(.))
  
  varimp_summ[[i]] <-  bind_rows(varimp_raw %>% summarise(across(everything(), mean)),
                                 varimp_raw %>% summarise(across(everything(), sd))) %>%
    t() %>%
    data.frame(var = row.names(.), Q = paste0("Q",i)) %>%
    relocate(var) %>%
    rename("mean" = "X1", "sd" = "X2")
  
}  


pal <- c("#2271B2",
         "#F748A5",
         "#359B73",
         "#e69f00")

df <- varimp_summ %>% 
  bind_rows %>%
  mutate(var = gsub("_first_quart$|_second_quart$|_third_quart$|_fourth_quart$", "", var),  # determine plot labels
         var = case_when(var == "around_surf" ~ "abundance: surface-feeders", 
                         var == "below_surf" ~ "abundance: sub-surface feeders", 
                         var == "host_dist" ~ "avg. phylo dist to host", 
                         var == "migr" ~ "abundance: migratory", 
                         var == "chicken_density" ~ "chicken density", 
                         var == "duck_density_2010" ~ "duck density", 
                         var == "mean_diff" ~ "temperature range", 
                         var == "mean_tmax" ~ "max temperature", 
                         var == "mean_tmin" ~ "min temperature", 
                         var == "mean_prec" ~ "total rainfall", 
                         var == "dist_to_coast_km" ~ "dist. to coast", 
                         var == "dist_to_water" ~ "dist. to inland water", 
                         var == "elev_min" ~ "min altitude", 
                         var == "elev_max" ~ "max altitude",
                         var == "elev_diff" ~ "altitude range",
                         var == "ndvi" ~ "vegetation index", 
                         var == "lc_1" ~ "water bodies",
                         var == "lc_2" ~ "evergreen needleleaf forests",
                         var == "lc_3" ~ "evergreen broadleaf forests",
                         var == "lc_4" ~ "deciduous needleleaf forests",
                         var == "lc_5" ~ "deciduous broadleaf forests",
                         var == "lc_6" ~ "mixed forests",
                         var == "lc_7" ~ "closed shrublands",
                         var == "lc_8" ~ "open_shrublands",
                         var == "lc_9" ~ "woody savannas",
                         var == "lc_10" ~ "savannas",
                         var == "lc_11" ~ "grasslands",
                         var == "lc_12" ~ "permanent wetlands",
                         var == "lc_13" ~ "croplands",
                         var == "lc_14" ~ "urband and built-up lands",
                         var == "lc_15" ~ "cropland/natural vegetation mosaics",
                         var == "lc_16" ~ "non-vegetated lands",
                         var == "lc_17" ~ "unclassified land"
         ),
         var = fct_relevel(var, c("abundance: surface-feeders",         # determine plot order
                                  "abundance: sub-surface feeders", 
                                  "abundance: migratory",
                                  "avg. phylo dist to host", 
                                  "chicken density", 
                                  "duck density", 
                                  "temperature range", 
                                  "max temperature", 
                                  "min temperature", 
                                  "total rainfall", 
                                  "dist. to coast", 
                                  "dist. to inland water", 
                                  "min altitude", 
                                  "max altitude",
                                  "altitude range",
                                  "vegetation index", 
                                  "water bodies",
                                  "evergreen needleleaf forests",
                                  "evergreen broadleaf forests",
                                  "deciduous needleleaf forests",
                                  "deciduous broadleaf forests",
                                  "mixed forests",
                                  "closed shrublands",
                                  "open_shrublands",
                                  "woody savannas",
                                  "savannas",
                                  "grasslands",
                                  "permanent wetlands",
                                  "croplands",
                                  "urband and built-up lands",
                                  "cropland/natural vegetation mosaics",
                                  "non-vegetated lands",
                                  "unclassified land"
         )),
         upper = mean + sd, lower = mean - sd)
  g <- ggplot(df, aes(x = var, y = mean, ymin = lower, ymax = upper, color = Q)) + 
  geom_errorbar(width=0, position=position_dodge(0.6)) + 
  geom_point(position=position_dodge(0.6)) +
    geom_vline(xintercept=seq(1.5, nrow(df)-0.5, 1), 
               lwd=.5, colour="grey75") + 
  scale_colour_manual("",
                      breaks = c("Q1", "Q2", "Q3", "Q4"),
                      values = pal) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle = 35, hjust = 1), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key.size=unit(0.8, "lines"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.background = element_rect(colour = "grey50"),
        axis.title.x=element_blank(),
        plot.margin = margin(10, 3, 3, 80),
        panel.grid.major.x = element_blank()) +
  ylab("Relative variable importance")

ggsave("plots/variable_importance_quarterly.png", plot = g, width = 10, height = 4.5)

