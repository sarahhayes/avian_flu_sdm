# Varimp plot

rm(list = ls())

PATH_TO_DATA <- "../../../OneDrive - The University of Liverpool/"

library(tidyverse)
library(dbarts)

varimp_summ <- list()

for(i in 1:4){
  
  load(file = paste(PATH_TO_DATA,
                    "AI_S2_SDM_storage/fitted-BART-models/sdm_Q",
                    i,
                    ".rds",
                    sep = ""))
  
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
                         var == "cong" ~ "abundance: congregative", 
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
         ))
ordered_var <- unique(as.character(df$var[order(df$mean, decreasing = TRUE)]))
df <- df %>% mutate(var = fct_relevel(var, ordered_var),
         upper = mean + sd, lower = mean - sd)

fig_varimp <- ggplot(df, aes(x = var, y = mean, ymin = lower, ymax = upper, color = Q)) + 
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

ggsave("plots/variable_importance_quarterly_reorder.png", plot = fig_varimp, width = 10, height = 4.5)


# Rank variables by variable importance (sum; those important in more models will rank higher)
varimp_rank <- varimp_summ %>%
  bind_rows %>%
  mutate(var = gsub("_first_quart$|_second_quart$|_third_quart$|_fourth_quart$", "", var)) %>%
  group_by(var) %>% 
  summarise(rank = sum(mean)) %>%
  arrange(-rank)

# # Calc and bind pd across all quarters for four most important variables by varimp
# # TO DO: line up variable names with _first_quart etc as in the below
#
# pd_summ <- replicate(4, vector("list", 4), simplify = FALSE) # initialise empty lists
# 
# for(i in 1:4){
#   
#   load(paste0("output/fitted-BART-models/sdm_Q",i,".rds"))
#   
#   for(j in 1:4){
#     
#     # Adapted from embarcadero::partial
#     raw <- sdm$fit$data@x[, varimp_rank$var[j]]
#     lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
#     pd <- pdbart(sdm, xind = varimp_rank$var[j], levs = lev, pl = FALSE)
#     
#     pd_summ[[i]][[j]] <- data.frame(var = varimp_rank$var[j], 
#                                     x = unlist(lev), 
#                                     y = pd$fd[[1]] %>% apply(., 2, median),
#                                     upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025),
#                                     lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975),
#                                     Q = paste0("Q",i))
#   }  
# }
# 
# 
# fig_pd_best <- pd_summ %>% 
#   bind_rows %>%
#   mutate(var = gsub("_first_quart$|_second_quart$|_third_quart$|_fourth_quart$", "", var)) +
#   ggplot(aes(x = x, y = y, ymin = lower, ymax = upper, fill = Q, color = Q)) +
#   geom_ribbon(alpha = 0.15, colour = NA) +
#   geom_line(lwd = 1.5, alpha = 0.6) +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_colour_manual("",
#                       breaks = c("Q1", "Q2", "Q3", "Q4"),
#                       values = pal) +
#   scale_fill_manual("",
#                     breaks = c("Q1", "Q2", "Q3", "Q4"),
#                     values = pal) +
#   ylab("Probability") +
#   theme_bw() + 
#   theme(axis.title.x=element_blank()) +
#   facet_wrap(~ var, scales = "free_x")
# 
# 
# ggsave("plots/partial_dependence_quarterly.png", plot = fig_pd_best, width = 10, height = 4.5)

# Calc and bind pd across all quarters for given variables by name
# TO DO: catch cases where variables are binary (land cover)

vars <- c("elev_min","below_surf")

pd_summ <- replicate(4, vector("list", length(vars)), simplify = FALSE) # initialise empty lists

for(i in 1:4){
  
  load(file = paste(PATH_TO_DATA,
                    "AI_S2_SDM_storage/fitted-BART-models/sdm_Q",
                    i,
                    ".rds",
                    sep = ""))
  
  varimp_raw <- sdm$varcount %>% 
    as.data.frame %>%
    mutate(./rowSums(.)) %>%
    mutate(./max(.))
  
  for(j in 1:length(vars)){
    
    # Adapted from embarcadero::partial
    fullvarname <- sdm$fit$data@x %>% as.data.frame %>% select(matches(vars[j])) %>% names
    if (length(fullvarname) == 1) {
    raw <- sdm$fit$data@x[, fullvarname]
    lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
    pd <- pdbart(sdm, xind = fullvarname, levs = lev, pl = FALSE)
    
    pd_summ[[i]][[j]] <- data.frame(var = vars[j], 
                                    x = unlist(lev), 
                                    y = pd$fd[[1]] %>% apply(., 2, median) %>% pnorm,
                                    upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025) %>% pnorm,
                                    lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975) %>% pnorm,
                                    Q = paste0("Q",i))
    }
  }  
}

fig_pd_chosen <- pd_summ %>%
  bind_rows %>%
  mutate(var = gsub("_first_quart$|_second_quart$|_third_quart$|_fourth_quart$", "", var),
         var = case_when(var == "below_surf" ~ "abundance: sub-surface feeders",
                         var == "elev_min" ~ "min altitude (m)"),
  ) %>%
  ggplot(aes(x = x, y = y, ymin = lower, ymax = upper, fill = Q, color = Q)) +
  geom_ribbon(alpha = 0.08, colour = NA) +
  geom_line(lwd = 0.8, alpha = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_colour_manual("",
                      breaks = c("Q1", "Q2", "Q3", "Q4"),
                      values = pal) +
  scale_fill_manual("",
                    breaks = c("Q1", "Q2", "Q3", "Q4"),
                    values = pal) +
  ylab("Probability") +
  theme_bw() +
  theme(axis.title.x=element_blank()) +
  facet_wrap(~ var, scales = "free_x")


ggsave("plots/partial_dependence_poster.png", plot = fig_pd_chosen, width = 6, height = 3)
