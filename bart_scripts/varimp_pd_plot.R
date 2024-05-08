# Varimp plot

rm(list = ls())

ALL_VARS_FOR_PD = FALSE # If set to true we do partial dependence on everything

# Optional command line arguments, must be passed as strings:
args <- commandArgs(trailingOnly = T)
if (length(args)<4){
  INCLUDE_CROSSTERMS <- "with-crossterms" # Set to "no-crossterms" to do model without crossterms or "with-crossterms" to do model with crossterms
}else{
  INCLUDE_CROSSTERMS <- args[4]
}
if (length(args)<3){
  CV_OR_RI <- "cv" # Set to "cv" to do crossvalidated model or "ri" to do crossvalidation + random intercept model
}else{
  CV_OR_RI <- args[3]
}
if (length(args)<1){
  # Set path to folder containing data, and where output will be stored
  PATH_TO_OUTPUTS <- "../../../OneDrive - The University of Liverpool/AI_S2_SDM_storage/"
}else{
  PATH_TO_OUTPUTS <- args[1]
}

library(tidyverse)
library(dbarts)
library(patchwork)

#### First do dataset A ####

varimp_summ <- list()

for(idx in 1:4){
  
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_model_with_vs_A_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
  varimp_raw <- sdm$varcount %>% 
    as.data.frame %>%
    mutate(./rowSums(.)) %>%
    mutate(./max(.))
  
  varimp_summ[[idx]] <-  bind_rows(varimp_raw %>% summarise(across(everything(), mean)),
                                   varimp_raw %>% summarise(across(everything(), sd))) %>% # IQR, 95% CI?
    t() %>%
    data.frame(Q = paste0("A Q",idx)) %>%
    rownames_to_column(var = "var") %>%
    rename("mean" = "X1", "sd" = "X2") %>%
    mutate(var = gsub("first_quart", "q1", var),
           var = gsub("second_quart", "q2", var),
           var = gsub("third_quart", "q3", var),
           var = gsub("fourth_quart", "q4", var)) %>%
    mutate(varq = str_extract(var, "_q\\d") %>% str_sub(-1) %>% as.numeric) %>%
    rowwise %>%
    mutate(var = gsub("_q[1-4]", paste0("_lag", idx-varq), var)) %>%
    mutate(var = gsub("-3", "1", var),
           var = gsub("-2", "2", var),
           var = gsub("-1", "3", var)) %>%
    select(-varq)
  
}  

if (INCLUDE_CROSSTERMS=="with-crossterms"){
  # Calculate a dataframe outlining amount of cross-season interaction
  crossterm_df <- lapply(1:4,
                         FUN = function(i){
                           c(length(grep("lag0",varimp_summ[[i]]$var)),
                             length(grep("lag1",varimp_summ[[i]]$var)),
                             length(grep("lag2",varimp_summ[[i]]$var)),
                             length(grep("lag3",varimp_summ[[i]]$var))
                           )}) %>%
    as.data.frame(row.names = c("Q1 model",
                                "Q2 model",
                                "Q3 model",
                                "Q4 model"),
                  col.names = c("lag 0 covs",
                                "lag 1 covs",
                                "lag 2 covs",
                                "lag 3 covs"))
}

pal <- c("#2271B2",
         "#F748A5",
         "#359B73",
         "#e69f00")

df <- varimp_summ %>% 
  bind_rows %>%
  # Set plotting labels
  mutate(var = gsub("anatinae" , "abundance: anatinae", var),
         var = gsub("anserinae" , "abundance: anserinae", var),
         var = gsub("ardeidae" , "abundance: ardeidae", var),
         var = gsub("arenaria_calidris" , "combined abundance: arenaria/calidris", var),
         var = gsub("aythyini" , "abundance: aythyini", var),
         var = gsub("laridae" , "abundance: laridae", var),
         var = gsub("around_surf" , "abundance: surface-feeders", var),
         var = gsub("below_surf" , "abundance: sub-surface feeders", var),
         var = gsub("plant" , "abundance: plant diet", var),
         var = gsub("scav" , "abundance: scavengers", var),
         var = gsub("vend" , "abundance: endotherm diet", var),
         var = gsub("host_dist" , "avg. phylo dist to host", var),
         var = gsub("migr" , "abundance: migratory", var),
         var = gsub("cong" , "abundance: congregative", var),
         var = gsub("species_richness" , "species richness", var),
         var = gsub("chicken_density_2010" , "chicken density", var),
         var = gsub("duck_density_2010" , "duck density", var),
         var = gsub("mean_relative_humidity" , "mean relative humidity", var),
         var = gsub("mean_diff" , "temperature range", var),
         var = gsub("variation_in_quarterly_mean_temp" , "temperature variation", var),
         var = gsub("mean_temp" , "mean temperature", var),
         var = gsub("mean_tmax" , "max temperature", var),
         var = gsub("mean_tmin" , "min temperature", var),
         var = gsub("isotherm_mean" , "mean isotherm", var),
         var = gsub("isotherm_midday_days_below1" , "frequency of days w/ isotherm <1", var),
         var = gsub("mean_prec" , "total rainfall", var),
         var = gsub("dist_to_coast_km" , "dist. to coast", var),
         var = gsub("dist_to_water" , "dist. to inland water", var),
         var = gsub("elevation_min" , "min altitude", var),
         var = gsub("elevation_max" , "max altitude", var),
         var = gsub("elevation_mode" , "modal altitude", var),
         var = gsub("elevation_diff" , "altitude range", var),
         var = gsub("ndvi" , "vegetation index", var),
         var = gsub("Water_bodies" , "water bodies", var),
         var = gsub("Evergreen_Needleleaf_Forests" , "evergreen needleleaf forests", var),
         var = gsub("Evergreen_Broadleaf_Forests" , "evergreen broadleaf forests", var),
         var = gsub("Deciduous_Needleleaf_Forests" , "deciduous needleleaf forests", var),
         var = gsub("Deciduous_Broadleaf_Forests" , "deciduous broadleaf forests", var),
         var = gsub("Mixed_Forests" , "mixed forests", var),
         var = gsub("Closed_Shrublands" , "closed shrublands", var),
         var = gsub("Open_Shrublands" , "open_shrublands", var),
         var = gsub("Woody_Savannas" , "woody savannas", var),
         var = gsub("Savannas" , "savannas", var),
         var = gsub("Grasslands" , "grasslands", var),
         var = gsub("Permanent_Wetlands" , "permanent wetlands", var),
         var = gsub("Croplands" , "croplands", var),
         var = gsub("Urban_and_Built-up_Lands" , "urband and built-up lands", var),
         var = gsub("Cropland/Natural_Vegetation_Mosaics" , "cropland/natural vegetation mosaics", var),
         var = gsub("Non-Vegetated_Lands" , "non-vegetated lands", var),
         var = gsub("Unclassified" , "unclassified land", var),
         var = gsub("_2022", "", var),
         var = gsub("_lag0", "", var),
         var = gsub("_lag", ", q-", var)
  )

ordered_var <- unique(as.character(df$var[order(df$mean, decreasing = TRUE)]))
df <- df %>% mutate(var = fct_relevel(var, ordered_var),
                    upper = mean + sd, lower = mean - sd) %>% mutate(var = fct_relevel(var, ordered_var)) %>% arrange(var)

fig_varimp <- ggplot(df, aes(x = var, y = mean, ymin = lower, ymax = upper, color = Q)) + 
  geom_errorbar(width=0, position=position_dodge(0.6)) + 
  geom_point(position=position_dodge(0.6)) +
  geom_vline(xintercept=seq(1.5, nrow(df)-0.5, 1), 
             lwd=.5, colour="grey75") + 
  scale_colour_manual("",
                      breaks = c("A Q1", "A Q2", "A Q3", "A Q4"),
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

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_variable_importance_quarterly_A.png", sep=""),
       plot = fig_varimp,
       width = 10,
       height = 4.5)


# Calc and bind pd across all quarters for given variables by name. All you need to do is set names of variables of interest (or don't do anything and let it plot all of them) :)

if (ALL_VARS_FOR_PD){
  vars <- gsub("\\..*",
               "",
               x=rownames(df)) %>%
    unique()
}else{
  vars = c("variation_in_quarterly_mean_temp_lag0", "variation_in_quarterly_mean_temp_lag2")
}

pd_summ <- replicate(4, vector("list", length(vars)), simplify = FALSE) # initialise empty lists

for(idx in 1:4){
  
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_model_with_vs_A_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
  colnames(sdm$fit$data@x) <- data.frame("var" = colnames(sdm$fit$data@x)) %>%
    mutate(var = gsub("first_quart", "q1", var),
           var = gsub("second_quart", "q2", var),
           var = gsub("third_quart", "q3", var),
           var = gsub("fourth_quart", "q4", var)) %>%
    mutate(varq = str_extract(var, "_q\\d") %>% str_sub(-1) %>% as.numeric) %>%
    rowwise %>%
    mutate(var = gsub("_q[1-4]", paste0("_lag", idx-varq), var)) %>%
    mutate(var = gsub("-3", "1", var),
           var = gsub("-2", "2", var),
           var = gsub("-1", "3", var)) %>%
    select(-varq) %>% pull(var)
  
  for(j in 1:length(vars)){
    
    # Adapted from embarcadero::partial
    fullvarname <- sdm$fit$data@x %>% as.data.frame %>% dplyr::select(matches(vars[j])) %>% names
    if (length(fullvarname) == 1) {
      
      # Check if the predictor variable is binary
      if (all(sdm$fit$data@x[, fullvarname] %in% 0:1)){
        raw <- sdm$fit$data@x[, fullvarname]
        lev <- list(c(0,1))
        pd <- pdbart(sdm, xind = fullvarname, levs = lev, pl = FALSE)
        
      } else {
        raw <- sdm$fit$data@x[, fullvarname]
        lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
        pd <- pdbart(sdm, xind = fullvarname, levs = lev, pl = FALSE)
        
      }
      
      pd_summ[[idx]][[j]] <- data.frame(var = vars[j], 
                                        x = unlist(lev), 
                                        y = pd$fd[[1]] %>% apply(., 2, median) %>% pnorm,
                                        upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025) %>% pnorm,
                                        lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975) %>% pnorm,
                                        Q = paste0("A Q",idx))
    }
  }  
}

# Make individual panels of a partial dependence plot, automatically assigning line or point depending on whether variable is continuous or binary

panel_list <- vector("list", length(vars)) # initialise empty lists

for(j in 1:length(vars)){
  
  df <- pd_summ %>%
    bind_rows %>%
    filter(var == vars[j]) %>%
    # Set plotting labels
    mutate(var = gsub("anatinae" , "abundance: anatinae", var),
           var = gsub("anserinae" , "abundance: anserinae", var),
           var = gsub("ardeidae" , "abundance: ardeidae", var),
           var = gsub("arenaria_calidris" , "combined abundance: arenaria/calidris", var),
           var = gsub("aythyini" , "abundance: aythyini", var),
           var = gsub("laridae" , "abundance: laridae", var),
           var = gsub("around_surf" , "abundance: surface-feeders", var),
           var = gsub("below_surf" , "abundance: sub-surface feeders", var),
           var = gsub("plant" , "abundance: plant diet", var),
           var = gsub("scav" , "abundance: scavengers", var),
           var = gsub("vend" , "abundance: endotherm diet", var),
           var = gsub("host_dist" , "avg. phylo dist to host", var),
           var = gsub("migr" , "abundance: migratory", var),
           var = gsub("cong" , "abundance: congregative", var),
           var = gsub("species_richness" , "species richness", var),
           var = gsub("chicken_density_2010" , "chicken density", var),
           var = gsub("duck_density_2010" , "duck density", var),
           var = gsub("mean_relative_humidity" , "mean relative humidity", var),
           var = gsub("mean_diff" , "temperature range", var),
           var = gsub("variation_in_quarterly_mean_temp" , "temperature variation", var),
           var = gsub("mean_temp" , "mean temperature", var),
           var = gsub("mean_tmax" , "max temperature", var),
           var = gsub("mean_tmin" , "min temperature", var),
           var = gsub("isotherm_mean" , "mean isotherm", var),
           var = gsub("isotherm_midday_days_below1" , "frequency of days w/ isotherm <1", var),
           var = gsub("mean_prec" , "total rainfall", var),
           var = gsub("dist_to_coast_km" , "dist. to coast", var),
           var = gsub("dist_to_water" , "dist. to inland water", var),
           var = gsub("elevation_min" , "min altitude", var),
           var = gsub("elevation_max" , "max altitude", var),
           var = gsub("elevation_mode" , "modal altitude", var),
           var = gsub("elevation_diff" , "altitude range", var),
           var = gsub("ndvi" , "vegetation index", var),
           var = gsub("Water_bodies" , "water bodies", var),
           var = gsub("Evergreen_Needleleaf_Forests" , "evergreen needleleaf forests", var),
           var = gsub("Evergreen_Broadleaf_Forests" , "evergreen broadleaf forests", var),
           var = gsub("Deciduous_Needleleaf_Forests" , "deciduous needleleaf forests", var),
           var = gsub("Deciduous_Broadleaf_Forests" , "deciduous broadleaf forests", var),
           var = gsub("Mixed_Forests" , "mixed forests", var),
           var = gsub("Closed_Shrublands" , "closed shrublands", var),
           var = gsub("Open_Shrublands" , "open_shrublands", var),
           var = gsub("Woody_Savannas" , "woody savannas", var),
           var = gsub("Savannas" , "savannas", var),
           var = gsub("Grasslands" , "grasslands", var),
           var = gsub("Permanent_Wetlands" , "permanent wetlands", var),
           var = gsub("Croplands" , "croplands", var),
           var = gsub("Urban_and_Built-up_Lands" , "urband and built-up lands", var),
           var = gsub("Cropland/Natural_Vegetation_Mosaics" , "cropland/natural vegetation mosaics", var),
           var = gsub("Non-Vegetated_Lands" , "non-vegetated lands", var),
           var = gsub("Unclassified" , "unclassified land", var),
           var = gsub("_2022", "", var),
           var = gsub("_lag0", "", var),
           var = gsub("_lag", ", q-", var)
    )
  
  xlab <- unique(df$var)
  
  panel_list[[j]] <- df %>% 
    ggplot(aes(x = x, y = y, ymin = lower, ymax = upper, fill = Q, color = Q)) +
    {if(!(all(df$x %in% 0:1)))
      list(
        geom_ribbon(alpha = 0.08, colour = NA),
        geom_line(lwd = 0.8, alpha = 0.4),
        scale_x_continuous(expand = c(0, 0))
      )} +
    {if(all(df$x %in% 0:1)) 
      list(
        geom_pointrange(aes(x = factor(x)), size = 0.7, position = position_dodge(width = 0.5), alpha = 0.8)
      )} +
    scale_colour_manual("",
                        breaks = c("A Q1", "A Q2", "A Q3", "A Q4"),
                        values = pal) +
    scale_fill_manual("",
                      breaks = c("A Q1", "A Q2", "A Q3", "A Q4"),
                      values = pal) +
    xlab(xlab) +
    ylab("Probability") +
    theme_bw() +
    {if(all(df$x %in% 0:1)) 
      list(
        guides(colour = "none", fill = "none") 
      )}
  
}

fig_pd_chosen <- wrap_plots(panel_list, ncol = 2) +
  plot_layout(guides = "collect", axis_titles = "collect")

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_partial_dependence_A.png", sep=""),
       plot = fig_pd_chosen,
       width = 9,
       height = 3)


#### Now do dataset B ####

varimp_summ <- list()

for(idx in 1:4){
  
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_model_with_vs_B_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
  varimp_raw <- sdm$varcount %>% 
    as.data.frame %>%
    mutate(./rowSums(.)) %>%
    mutate(./max(.))
  
  varimp_summ[[idx]] <-  bind_rows(varimp_raw %>% summarise(across(everything(), mean)),
                                   varimp_raw %>% summarise(across(everything(), sd))) %>% # IQR, 95% CI?
    t() %>%
    data.frame(Q = paste0("B Q",idx)) %>%
    rownames_to_column(var = "var") %>%
    rename("mean" = "X1", "sd" = "X2") %>%
    mutate(var = gsub("first_quart", "q1", var),
           var = gsub("second_quart", "q2", var),
           var = gsub("third_quart", "q3", var),
           var = gsub("fourth_quart", "q4", var)) %>%
    mutate(varq = str_extract(var, "_q\\d") %>% str_sub(-1) %>% as.numeric) %>%
    rowwise %>%
    mutate(var = gsub("_q[1-4]", paste0("_lag", idx-varq), var)) %>%
    mutate(var = gsub("-3", "1", var),
           var = gsub("-2", "2", var),
           var = gsub("-1", "3", var)) %>%
    select(-varq)
  
}  

if (INCLUDE_CROSSTERMS=="with-crossterms"){
  # Calculate a dataframe outlining amount of cross-season interaction
  crossterm_df <- lapply(1:4,
                         FUN = function(i){
                           c(length(grep("lag0",varimp_summ[[i]]$var)),
                             length(grep("lag1",varimp_summ[[i]]$var)),
                             length(grep("lag2",varimp_summ[[i]]$var)),
                             length(grep("lag3",varimp_summ[[i]]$var))
                           )}) %>%
    as.data.frame(row.names = c("Q1 model",
                                "Q2 model",
                                "Q3 model",
                                "Q4 model"),
                  col.names = c("lag 0 covs",
                                "lag 1 covs",
                                "lag 2 covs",
                                "lag 3 covs"))
}

pal <- c("#2271B2",
         "#F748A5",
         "#359B73",
         "#e69f00")

df <- varimp_summ %>% 
  bind_rows %>%
  # Set plotting labels
  mutate(var = gsub("anatinae" , "abundance: anatinae", var),
         var = gsub("anserinae" , "abundance: anserinae", var),
         var = gsub("ardeidae" , "abundance: ardeidae", var),
         var = gsub("arenaria_calidris" , "combined abundance: arenaria/calidris", var),
         var = gsub("aythyini" , "abundance: aythyini", var),
         var = gsub("laridae" , "abundance: laridae", var),
         var = gsub("around_surf" , "abundance: surface-feeders", var),
         var = gsub("below_surf" , "abundance: sub-surface feeders", var),
         var = gsub("plant" , "abundance: plant diet", var),
         var = gsub("scav" , "abundance: scavengers", var),
         var = gsub("vend" , "abundance: endotherm diet", var),
         var = gsub("host_dist" , "avg. phylo dist to host", var),
         var = gsub("migr" , "abundance: migratory", var),
         var = gsub("cong" , "abundance: congregative", var),
         var = gsub("species_richness" , "species richness", var),
         var = gsub("chicken_density_2010" , "chicken density", var),
         var = gsub("duck_density_2010" , "duck density", var),
         var = gsub("mean_relative_humidity" , "mean relative humidity", var),
         var = gsub("mean_diff" , "temperature range", var),
         var = gsub("variation_in_quarterly_mean_temp" , "temperature variation", var),
         var = gsub("mean_temp" , "mean temperature", var),
         var = gsub("mean_tmax" , "max temperature", var),
         var = gsub("mean_tmin" , "min temperature", var),
         var = gsub("isotherm_mean" , "mean isotherm", var),
         var = gsub("isotherm_midday_days_below1" , "frequency of days w/ isotherm <1", var),
         var = gsub("mean_prec" , "total rainfall", var),
         var = gsub("dist_to_coast_km" , "dist. to coast", var),
         var = gsub("dist_to_water" , "dist. to inland water", var),
         var = gsub("elevation_min" , "min altitude", var),
         var = gsub("elevation_max" , "max altitude", var),
         var = gsub("elevation_mode" , "modal altitude", var),
         var = gsub("elevation_diff" , "altitude range", var),
         var = gsub("ndvi" , "vegetation index", var),
         var = gsub("Water_bodies" , "water bodies", var),
         var = gsub("Evergreen_Needleleaf_Forests" , "evergreen needleleaf forests", var),
         var = gsub("Evergreen_Broadleaf_Forests" , "evergreen broadleaf forests", var),
         var = gsub("Deciduous_Needleleaf_Forests" , "deciduous needleleaf forests", var),
         var = gsub("Deciduous_Broadleaf_Forests" , "deciduous broadleaf forests", var),
         var = gsub("Mixed_Forests" , "mixed forests", var),
         var = gsub("Closed_Shrublands" , "closed shrublands", var),
         var = gsub("Open_Shrublands" , "open_shrublands", var),
         var = gsub("Woody_Savannas" , "woody savannas", var),
         var = gsub("Savannas" , "savannas", var),
         var = gsub("Grasslands" , "grasslands", var),
         var = gsub("Permanent_Wetlands" , "permanent wetlands", var),
         var = gsub("Croplands" , "croplands", var),
         var = gsub("Urban_and_Built-up_Lands" , "urband and built-up lands", var),
         var = gsub("Cropland/Natural_Vegetation_Mosaics" , "cropland/natural vegetation mosaics", var),
         var = gsub("Non-Vegetated_Lands" , "non-vegetated lands", var),
         var = gsub("Unclassified" , "unclassified land", var),
         var = gsub("_2022", "", var),
         var = gsub("_lag0", "", var),
         var = gsub("_lag", ", q-", var)
  )

ordered_var <- unique(as.character(df$var[order(df$mean, decreasing = TRUE)]))
df <- df %>% mutate(var = fct_relevel(var, ordered_var),
                    upper = mean + sd, lower = mean - sd) %>% mutate(var = fct_relevel(var, ordered_var)) %>% arrange(var)

fig_varimp <- ggplot(df, aes(x = var, y = mean, ymin = lower, ymax = upper, color = Q)) + 
  geom_errorbar(width=0, position=position_dodge(0.6)) + 
  geom_point(position=position_dodge(0.6)) +
  geom_vline(xintercept=seq(1.5, nrow(df)-0.5, 1), 
             lwd=.5, colour="grey75") + 
  scale_colour_manual("",
                      breaks = c("B Q1", "B Q2", "B Q3", "B Q4"),
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

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_variable_importance_quarterly_B.png", sep=""),
       plot = fig_varimp,
       width = 10,
       height = 4.5)


# Calc and bind pd across all quarters for given variables by name. All you need to do is set names of variables of interest (or don't do anything and let it plot all of them) :)

if (ALL_VARS_FOR_PD){
  vars <- gsub("\\..*",
               "",
               x=rownames(df)) %>%
    unique()
}else{
  vars = c("dist_to_coast_km", "plant_lag0")
}

pd_summ <- replicate(4, vector("list", length(vars)), simplify = FALSE) # initialise empty lists

for(idx in 1:4){
  
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_model_with_vs_B_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
  colnames(sdm$fit$data@x) <- data.frame("var" = colnames(sdm$fit$data@x)) %>%
    mutate(var = gsub("first_quart", "q1", var),
           var = gsub("second_quart", "q2", var),
           var = gsub("third_quart", "q3", var),
           var = gsub("fourth_quart", "q4", var)) %>%
    mutate(varq = str_extract(var, "_q\\d") %>% str_sub(-1) %>% as.numeric) %>%
    rowwise %>%
    mutate(var = gsub("_q[1-4]", paste0("_lag", idx-varq), var)) %>%
    mutate(var = gsub("-3", "1", var),
           var = gsub("-2", "2", var),
           var = gsub("-1", "3", var)) %>%
    select(-varq) %>% pull(var)
  
  for(j in 1:length(vars)){
    
    # Adapted from embarcadero::partial
    fullvarname <- sdm$fit$data@x %>% as.data.frame %>% dplyr::select(matches(vars[j])) %>% names
    if (length(fullvarname) == 1) {
      
      # Check if the predictor variable is binary
      if (all(sdm$fit$data@x[, fullvarname] %in% 0:1)){
        raw <- sdm$fit$data@x[, fullvarname]
        lev <- list(c(0,1))
        pd <- pdbart(sdm, xind = fullvarname, levs = lev, pl = FALSE)
        
      } else {
        raw <- sdm$fit$data@x[, fullvarname]
        lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
        pd <- pdbart(sdm, xind = fullvarname, levs = lev, pl = FALSE)
        
      }
      
      pd_summ[[idx]][[j]] <- data.frame(var = vars[j], 
                                        x = unlist(lev), 
                                        y = pd$fd[[1]] %>% apply(., 2, median) %>% pnorm,
                                        upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025) %>% pnorm,
                                        lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975) %>% pnorm,
                                        Q = paste0("B Q",idx))
    }
  }  
}

# Make individual panels of a partial dependence plot, automatically assigning line or point depending on whether variable is continuous or binary

panel_list <- vector("list", length(vars)) # initialise empty lists

for(j in 1:length(vars)){
  
  df <- pd_summ %>%
    bind_rows %>%
    filter(var == vars[j]) %>%
    # Set plotting labels
    mutate(var = gsub("anatinae" , "abundance: anatinae", var),
           var = gsub("anserinae" , "abundance: anserinae", var),
           var = gsub("ardeidae" , "abundance: ardeidae", var),
           var = gsub("arenaria_calidris" , "combined abundance: arenaria/calidris", var),
           var = gsub("aythyini" , "abundance: aythyini", var),
           var = gsub("laridae" , "abundance: laridae", var),
           var = gsub("around_surf" , "abundance: surface-feeders", var),
           var = gsub("below_surf" , "abundance: sub-surface feeders", var),
           var = gsub("plant" , "abundance: plant diet", var),
           var = gsub("scav" , "abundance: scavengers", var),
           var = gsub("vend" , "abundance: endotherm diet", var),
           var = gsub("host_dist" , "avg. phylo dist to host", var),
           var = gsub("migr" , "abundance: migratory", var),
           var = gsub("cong" , "abundance: congregative", var),
           var = gsub("species_richness" , "species richness", var),
           var = gsub("chicken_density_2010" , "chicken density", var),
           var = gsub("duck_density_2010" , "duck density", var),
           var = gsub("mean_relative_humidity" , "mean relative humidity", var),
           var = gsub("mean_diff" , "temperature range", var),
           var = gsub("variation_in_quarterly_mean_temp" , "temperature variation", var),
           var = gsub("mean_temp" , "mean temperature", var),
           var = gsub("mean_tmax" , "max temperature", var),
           var = gsub("mean_tmin" , "min temperature", var),
           var = gsub("isotherm_mean" , "mean isotherm", var),
           var = gsub("isotherm_midday_days_below1" , "frequency of days w/ isotherm <1", var),
           var = gsub("mean_prec" , "total rainfall", var),
           var = gsub("dist_to_coast_km" , "dist. to coast", var),
           var = gsub("dist_to_water" , "dist. to inland water", var),
           var = gsub("elevation_min" , "min altitude", var),
           var = gsub("elevation_max" , "max altitude", var),
           var = gsub("elevation_mode" , "modal altitude", var),
           var = gsub("elevation_diff" , "altitude range", var),
           var = gsub("ndvi" , "vegetation index", var),
           var = gsub("Water_bodies" , "water bodies", var),
           var = gsub("Evergreen_Needleleaf_Forests" , "evergreen needleleaf forests", var),
           var = gsub("Evergreen_Broadleaf_Forests" , "evergreen broadleaf forests", var),
           var = gsub("Deciduous_Needleleaf_Forests" , "deciduous needleleaf forests", var),
           var = gsub("Deciduous_Broadleaf_Forests" , "deciduous broadleaf forests", var),
           var = gsub("Mixed_Forests" , "mixed forests", var),
           var = gsub("Closed_Shrublands" , "closed shrublands", var),
           var = gsub("Open_Shrublands" , "open_shrublands", var),
           var = gsub("Woody_Savannas" , "woody savannas", var),
           var = gsub("Savannas" , "savannas", var),
           var = gsub("Grasslands" , "grasslands", var),
           var = gsub("Permanent_Wetlands" , "permanent wetlands", var),
           var = gsub("Croplands" , "croplands", var),
           var = gsub("Urban_and_Built-up_Lands" , "urband and built-up lands", var),
           var = gsub("Cropland/Natural_Vegetation_Mosaics" , "cropland/natural vegetation mosaics", var),
           var = gsub("Non-Vegetated_Lands" , "non-vegetated lands", var),
           var = gsub("Unclassified" , "unclassified land", var),
           var = gsub("_2022", "", var),
           var = gsub("_lag0", "", var),
           var = gsub("_lag", ", q-", var)
    )
  
  xlab <- unique(df$var)
  
  panel_list[[j]] <- df %>% 
    ggplot(aes(x = x, y = y, ymin = lower, ymax = upper, fill = Q, color = Q)) +
    {if(!(all(df$x %in% 0:1)))
      list(
        geom_ribbon(alpha = 0.08, colour = NA),
        geom_line(lwd = 0.8, alpha = 0.4),
        scale_x_continuous(expand = c(0, 0))
      )} +
    {if(all(df$x %in% 0:1)) 
      list(
        geom_pointrange(aes(x = factor(x)), size = 0.7, position = position_dodge(width = 0.5), alpha = 0.8)
      )} +
    scale_colour_manual("",
                        breaks = c("B Q1", "B Q2", "B Q3", "B Q4"),
                        values = pal) +
    scale_fill_manual("",
                      breaks = c("B Q1", "B Q2", "B Q3", "B Q4"),
                      values = pal) +
    xlab(xlab) +
    ylab("Probability") +
    theme_bw() +
    {if(all(df$x %in% 0:1)) 
      list(
        guides(colour = "none", fill = "none") 
      )}
  
}

fig_pd_chosen <- wrap_plots(panel_list, ncol = 2) +
  plot_layout(guides = "collect", axis_titles = "collect")

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_partial_dependence_B.png", sep=""),
       plot = fig_pd_chosen,
       width = 9,
       height = 3)
