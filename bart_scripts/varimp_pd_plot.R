# Varimp plot

rm(list = ls())

ALL_VARS_FOR_PD = TRUE # If set to true we do partial dependence on everything

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

#PATH_TO_OUTPUTS <- "H:/Working/avian_flu_sdm/output/"

library(tidyverse)
library(dbarts)
library(patchwork)
library(ggnewscale)

### Set plot variable order (should match tables)

ordered_var <- c(
  "mean relative humidity",
  "temperature range",
  "mean temperature",
  "temperature variation",
  "max temperature",
  "min temperature",
  "mean isotherm",
  "frequency of days w/ isotherm <1m",
  "total rainfall",
  "min altitude",
  "max altitude",
  "modal altitude",
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
  "urban and built-up lands",
  "cropland/natural vegetation mosaics",
  "non-vegetated lands",
  "unclassified land",
  "dist. to coast",
  "dist. to inland water",
  "chicken density",
  "dom. duck density", 
  "abundance: anatinae", 
  "abundance: anserinae", 
  "abundance: ardeidae",
  "combined abundance: arenaria/calidris",
  "abundance: aythyini",
  "abundance: laridae", 
  "abundance: surface-feeders",
  "abundance: sub-surface feeders",
  "abundance: plant diet",
  "abundance: scavengers",
  "abundance: endotherm diet",
  "abundance: congregative",
  "abundance: migratory",
  "species richness",
  "avg. phylo dist to host"
)

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
    mutate(lag = idx-varq) %>%
    mutate(lag = case_when(
      lag == -3 ~ 1,
      lag == -2 ~ 2,
      lag == -1 ~ 3,
      TRUE ~ lag
    )) %>%
    mutate(var = gsub("_q[1-4]", paste0("_lag", lag), var)) %>%
    mutate(lag = as.character(lag)) %>%
    mutate(lag = replace_na(lag, "static")) %>%
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
         var = gsub("duck_density_2010" , "dom. duck density", var),
         var = gsub("mean_relative_humidity" , "mean relative humidity", var),
         var = gsub("mean_diff" , "temperature range", var),
         var = gsub("variation_in_quarterly_mean_temp" , "temperature variation", var),
         var = gsub("mean_temp" , "mean temperature", var),
         var = gsub("mean_tmax" , "max temperature", var),
         var = gsub("mean_tmin" , "min temperature", var),
         var = gsub("isotherm_mean" , "mean isotherm", var),
         var = gsub("isotherm_midday_days_below1" , "frequency of days w/ isotherm <1m", var),
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
         var = gsub("Urban_and_Built-up_Lands" , "urban and built-up lands", var),
         var = gsub("Cropland/Natural_Vegetation_Mosaics" , "cropland/natural vegetation mosaics", var),
         var = gsub("Non-Vegetated_Lands" , "non-vegetated lands", var),
         var = gsub("Unclassified" , "unclassified land", var),
         var = gsub("_2022", "", var),
         var = gsub("_lag.$", "", var),
  )

#ordered_var <- unique(as.character(df$var[order(df$mean, decreasing = TRUE)])) # Order variables by mean varimp
df <- df %>% mutate(var = fct_relevel(var, ordered_var),
                    upper = mean + sd, lower = mean - sd) %>% mutate(var = fct_relevel(var, ordered_var)) %>% arrange(var)

fig_varimp_A <- ggplot(df, aes(x = mean, y = var, xmin = lower, xmax = upper, color = lag)) + 
  geom_errorbar(width=0, position=position_dodge(0.6)) + 
  geom_point(position=position_dodge(0.6)) +
  geom_hline(yintercept=seq(1.5, nrow(df)-0.5, 1), lwd=.5, colour="grey75") +  # add custom gridlines
  geom_hline(yintercept=c(4.5,5.5,10.5), lwd=1, colour="grey75") + # separate bioclimatic, topographic, livestock, host ecological variables - must set manually each recalculation!
  scale_x_continuous(position="top") +
  scale_y_discrete(limits=rev) +
  scale_colour_manual(name = element_blank(),
                      breaks = c("0", "1", "2", "3"),
                      values = c(pal, "black")) +
  theme_bw(base_size = 16) +
  theme(legend.title=element_blank(), 
        legend.position=c(-0.325, 0.975),
        legend.key.size=unit(0.8, "lines"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.background = element_rect(colour = "grey50"),
        axis.title.y=element_blank(),
        plot.margin = margin(3, 3, 3, 3),
        panel.grid.major.y = element_blank()) +
  facet_wrap(~ Q, nrow=1, strip.position="bottom") +
  xlab("Relative variable importance")

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_variable_importance_quarterly_A.png", sep=""),
       plot = fig_varimp_A,
       width = 10,
       height = 4.5)

# Calc and bind pd across all quarters for given variables by name. All you need to do is set names of variables of interest (or don't do anything and let it plot all of them) :)

for(idx in 1:4){
  
  if (ALL_VARS_FOR_PD){
    
    load(file = paste(PATH_TO_OUTPUTS,
                      "fitted-BART-models-",
                      INCLUDE_CROSSTERMS,
                      "/",
                      CV_OR_RI,
                      "_model_with_vs_A_Q",
                      idx,
                      ".rds",
                      sep = ""))
    
    vars <- sdm$varcount %>% as.data.frame %>% names
    
  }else{
    vars <- c("variation_in_quarterly_mean_temp_q1",  
              "variation_in_quarterly_mean_temp_q2",
              "variation_in_quarterly_mean_temp_q3",
              "variation_in_quarterly_mean_temp_q4")
  }
  
  pd_summ <- list() # initialise empty list
  
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_model_with_vs_A_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
  for(j in 1:length(vars)){
    
    # Adapted from embarcadero::partial
    if (vars[j] %in% (sdm$fit$data@x %>% attr("term.labels"))) {
      
      # Check if the predictor variable is binary
      if (all(sdm$fit$data@x[, vars[j]] %in% 0:1)){
        raw <- sdm$fit$data@x[, vars[j]]
        lev <- list(c(0,1))
        pd <- pdbart(sdm, xind = vars[j], levs = lev, pl = FALSE)
        
      } else {
        raw <- sdm$fit$data@x[, vars[j]]
        lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
        pd <- pdbart(sdm, xind = vars[j], levs = lev, pl = FALSE)
        
      }
      
      pd_summ[[j]] <- data.frame(var = vars[j], 
                                 x = unlist(lev), 
                                 y = pd$fd[[1]] %>% apply(., 2, median) %>% pnorm,
                                 upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025) %>% pnorm,
                                 lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975) %>% pnorm,
                                 Q = paste0("A Q",idx)) %>%
        mutate(var = gsub("first_quart", "q1", var),
               var = gsub("second_quart", "q2", var),
               var = gsub("third_quart", "q3", var),
               var = gsub("fourth_quart", "q4", var)) %>%
        mutate(varq = str_extract(var, "_q\\d") %>% str_sub(-1) %>% as.numeric) %>%
        rowwise %>%
        mutate(lag = idx-varq) %>%
        mutate(lag = case_when(
          lag == -3 ~ 1,
          lag == -2 ~ 2,
          lag == -1 ~ 3,
          TRUE ~ lag
        )) %>%
        mutate(var = gsub("_q[1-4]", paste0("_lag", lag), var)) %>%
        mutate(lag = as.character(lag)) %>%
        mutate(lag = replace_na(lag, "static")) %>%
        select(-varq)
    }
  }  
  
  # Make individual panels of a partial dependence plot, automatically assigning line or point depending on whether variable is continuous or binary
  
  panel_list <- list() # initialise empty lists
  plot_list <- list()
  
  pd_summ <- pd_summ %>%
    bind_rows %>%
    as.data.frame %>%
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
           var = gsub("duck_density_2010" , "dom. duck density", var),
           var = gsub("mean_relative_humidity" , "mean relative humidity", var),
           var = gsub("mean_diff" , "temperature range", var),
           var = gsub("variation_in_quarterly_mean_temp" , "temperature variation", var),
           var = gsub("mean_temp" , "mean temperature", var),
           var = gsub("mean_tmax" , "max temperature", var),
           var = gsub("mean_tmin" , "min temperature", var),
           var = gsub("isotherm_mean" , "mean isotherm", var),
           var = gsub("isotherm_midday_days_below1" , "frequency of days w/ isotherm <1m", var),
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
           var = gsub("Urban_and_Built-up_Lands" , "urban and built-up lands", var),
           var = gsub("Cropland/Natural_Vegetation_Mosaics" , "cropland/natural vegetation mosaics", var),
           var = gsub("Non-Vegetated_Lands" , "non-vegetated lands", var),
           var = gsub("Unclassified" , "unclassified land", var),
           var = gsub("_2022", "", var),
           var = gsub("_lag.$", "", var),
    ) %>%
    mutate(var = fct_relevel(var, ordered_var)) %>% 
    arrange(var)
  
  for(j in 1:length(unique(pd_summ$var))){
    
    df <- pd_summ %>%
      filter(var == unique(pd_summ$var)[j])
    
    xlab <- unique(df$var)
    
    panel_list[[j]] <- df %>% 
      ggplot(aes(x = x, y = y, ymin = lower, ymax = upper, fill = lag, color = lag)) +
      {if(!(all(df$x %in% 0:1)))
        list(
          geom_ribbon(alpha = 0.08, colour = NA),
          geom_line(lwd = 0.8, alpha = 0.4),
          scale_x_continuous(expand = c(0, 0)),
          scale_fill_manual(name = element_blank(),
                            breaks = c("0", "1", "2", "3"),
                            values = setNames(c(pal, "black"), c("0","1","2","3","static"))),
          scale_colour_manual(name = element_blank(),
                              breaks = c("0", "1", "2", "3"),
                              values = setNames(c(pal, "black"), c("0","1","2","3","static")))
        )} +
      {if(all(df$x %in% 0:1)) 
        list(
          geom_pointrange(aes(x = factor(x)), size = 0.7, position = position_dodge(width = 0.5), alpha = 0.8, color = "black", fill = "black")
        )} +
      xlab(xlab) +
      ylab("Probability") +
      theme_bw() +
      list(guides(colour = "none", fill = "none"))
    
  }
  
  dummy <- data.frame(x = rep(c(0),4), 
                      y = rep(c(0),4),
                      lower = rep(c(0),4),
                      upper = rep(c(0),4),
                      lag = c("0","1","2","3"))
  
  plot_list[[idx]] <- wrap_plots(panel_list, ncol = 3) +
    plot_layout(guides = "collect", axis_titles = "collect") +
    # Add legend manually
    new_scale_colour() +
    new_scale_fill() +
    geom_point(aes(colour = lag), data = dummy, alpha = 0) +
    scale_colour_manual(name = element_blank(),
                        breaks = c("0", "1", "2", "3"),
                        values = setNames(pal, c("0","1","2","3")))
  
  ggsave(paste("plots/",
               INCLUDE_CROSSTERMS,
               "_",
               CV_OR_RI,
               "_partial_dependence_A_Q", idx, ".png", sep=""),
         plot = plot_list[[idx]],
         width = 12,
         height = 6)
  
}

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
    mutate(lag = idx-varq) %>%
    mutate(lag = case_when(
      lag == -3 ~ 1,
      lag == -2 ~ 2,
      lag == -1 ~ 3,
      TRUE ~ lag
    )) %>%
    mutate(var = gsub("_q[1-4]", paste0("_lag", lag), var)) %>%
    mutate(lag = as.character(lag)) %>%
    mutate(lag = replace_na(lag, "static")) %>%
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
         var = gsub("duck_density_2010" , "dom. duck density", var),
         var = gsub("mean_relative_humidity" , "mean relative humidity", var),
         var = gsub("mean_diff" , "temperature range", var),
         var = gsub("variation_in_quarterly_mean_temp" , "temperature variation", var),
         var = gsub("mean_temp" , "mean temperature", var),
         var = gsub("mean_tmax" , "max temperature", var),
         var = gsub("mean_tmin" , "min temperature", var),
         var = gsub("isotherm_mean" , "mean isotherm", var),
         var = gsub("isotherm_midday_days_below1" , "frequency of days w/ isotherm <1m", var),
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
         var = gsub("Urban_and_Built-up_Lands" , "urban and built-up lands", var),
         var = gsub("Cropland/Natural_Vegetation_Mosaics" , "cropland/natural vegetation mosaics", var),
         var = gsub("Non-Vegetated_Lands" , "non-vegetated lands", var),
         var = gsub("Unclassified" , "unclassified land", var),
         var = gsub("_2022", "", var),
         var = gsub("_lag.$", "", var),
  )

#ordered_var <- unique(as.character(df$var[order(df$mean, decreasing = TRUE)])) # Order variables by mean varimp
df <- df %>% mutate(var = fct_relevel(var, ordered_var),
                    upper = mean + sd, lower = mean - sd) %>% mutate(var = fct_relevel(var, ordered_var)) %>% arrange(var)

fig_varimp_B <- ggplot(df, aes(x = mean, y = var, xmin = lower, xmax = upper, color = lag)) + 
  geom_errorbar(width=0, position=position_dodge(0.6)) + 
  geom_point(position=position_dodge(0.6)) +
  geom_hline(yintercept=seq(1.5, nrow(df)-0.5, 1), lwd=.5, colour="grey75") +  # add custom gridlines
  geom_hline(yintercept=c(9.5,15.5), lwd=1, colour="grey75") + # separate bioclimatic, topographic, livestock, host ecological variables - must set manually each recalculation!
  scale_x_continuous(position="top") +
  scale_y_discrete(limits=rev) +
  scale_colour_manual(name = element_blank(),
                      breaks = c("0", "1", "2", "3"),
                      values = c(pal, "black")) +
  theme_bw(base_size = 16) +
  theme(legend.title=element_blank(), 
        legend.position=c(-0.325, 0.975),
        legend.key.size=unit(0.8, "lines"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.background = element_rect(colour = "grey50"),
        axis.title.y=element_blank(),
        plot.margin = margin(3, 3, 3, 3),
        panel.grid.major.y = element_blank()) +
  facet_wrap(~ Q, nrow=1, strip.position="bottom") +
  xlab("Relative variable importance")

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_variable_importance_quarterly_B.png", sep=""),
       plot = fig_varimp_B,
       width = 10,
       height = 5.5)


# Calc and bind pd across all quarters for given variables by name. All you need to do is set names of variables of interest (or don't do anything and let it plot all of them) :)

for(idx in 1:4){
  
  if (ALL_VARS_FOR_PD){
    
    load(file = paste(PATH_TO_OUTPUTS,
                      "fitted-BART-models-",
                      INCLUDE_CROSSTERMS,
                      "/",
                      CV_OR_RI,
                      "_model_with_vs_B_Q",
                      idx,
                      ".rds",
                      sep = ""))
    
    vars <- sdm$varcount %>% as.data.frame %>% names
    
  }else{
    vars <- c("variation_in_quarterly_mean_temp_q1",  
              "variation_in_quarterly_mean_temp_q2",
              "variation_in_quarterly_mean_temp_q3",
              "variation_in_quarterly_mean_temp_q4")
  }
  
  pd_summ <- list() # initialise empty list
  
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_model_with_vs_B_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
  for(j in 1:length(vars)){
    
    # Adapted from embarcadero::partial
    if (vars[j] %in% (sdm$fit$data@x %>% attr("term.labels"))) {
      
      # Check if the predictor variable is binary
      if (all(sdm$fit$data@x[, vars[j]] %in% 0:1)){
        raw <- sdm$fit$data@x[, vars[j]]
        lev <- list(c(0,1))
        pd <- pdbart(sdm, xind = vars[j], levs = lev, pl = FALSE)
        
      } else {
        raw <- sdm$fit$data@x[, vars[j]]
        lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
        pd <- pdbart(sdm, xind = vars[j], levs = lev, pl = FALSE)
        
      }
      
      pd_summ[[j]] <- data.frame(var = vars[j], 
                                 x = unlist(lev), 
                                 y = pd$fd[[1]] %>% apply(., 2, median) %>% pnorm,
                                 upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025) %>% pnorm,
                                 lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975) %>% pnorm,
                                 Q = paste0("B Q",idx)) %>%
        mutate(var = gsub("first_quart", "q1", var),
               var = gsub("second_quart", "q2", var),
               var = gsub("third_quart", "q3", var),
               var = gsub("fourth_quart", "q4", var)) %>%
        mutate(varq = str_extract(var, "_q\\d") %>% str_sub(-1) %>% as.numeric) %>%
        rowwise %>%
        mutate(lag = idx-varq) %>%
        mutate(lag = case_when(
          lag == -3 ~ 1,
          lag == -2 ~ 2,
          lag == -1 ~ 3,
          TRUE ~ lag
        )) %>%
        mutate(var = gsub("_q[1-4]", paste0("_lag", lag), var)) %>%
        mutate(lag = as.character(lag)) %>%
        mutate(lag = replace_na(lag, "static")) %>%
        select(-varq)
    }
  }  
  
  # Make individual panels of a partial dependence plot, automatically assigning line or point depending on whether variable is continuous or binary
  
  panel_list <- list() # initialise empty lists
  plot_list <- list()
  
  pd_summ <- pd_summ %>%
    bind_rows %>%
    as.data.frame %>%
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
           var = gsub("duck_density_2010" , "dom. duck density", var),
           var = gsub("mean_relative_humidity" , "mean relative humidity", var),
           var = gsub("mean_diff" , "temperature range", var),
           var = gsub("variation_in_quarterly_mean_temp" , "temperature variation", var),
           var = gsub("mean_temp" , "mean temperature", var),
           var = gsub("mean_tmax" , "max temperature", var),
           var = gsub("mean_tmin" , "min temperature", var),
           var = gsub("isotherm_mean" , "mean isotherm", var),
           var = gsub("isotherm_midday_days_below1" , "frequency of days w/ isotherm <1m", var),
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
           var = gsub("Urban_and_Built-up_Lands" , "urban and built-up lands", var),
           var = gsub("Cropland/Natural_Vegetation_Mosaics" , "cropland/natural vegetation mosaics", var),
           var = gsub("Non-Vegetated_Lands" , "non-vegetated lands", var),
           var = gsub("Unclassified" , "unclassified land", var),
           var = gsub("_2022", "", var),
           var = gsub("_lag.$", "", var),
    ) %>%
    mutate(var = fct_relevel(var, ordered_var)) %>% 
    arrange(var)
  
  for(j in 1:length(unique(pd_summ$var))){
    
    df <- pd_summ %>%
      filter(var == unique(pd_summ$var)[j])
    
    xlab <- unique(df$var)
    
    panel_list[[j]] <- df %>% 
      ggplot(aes(x = x, y = y, ymin = lower, ymax = upper, fill = lag, color = lag)) +
      {if(!(all(df$x %in% 0:1)))
        list(
          geom_ribbon(alpha = 0.08, colour = NA),
          geom_line(lwd = 0.8, alpha = 0.4),
          scale_x_continuous(expand = c(0, 0)),
          scale_fill_manual(name = element_blank(),
                            breaks = c("0", "1", "2", "3"),
                            values = setNames(c(pal, "black"), c("0","1","2","3","static"))),
          scale_colour_manual(name = element_blank(),
                              breaks = c("0", "1", "2", "3"),
                              values = setNames(c(pal, "black"), c("0","1","2","3","static")))
        )} +
      {if(all(df$x %in% 0:1)) 
        list(
          geom_pointrange(aes(x = factor(x)), size = 0.7, position = position_dodge(width = 0.5), alpha = 0.8, color = "black", fill = "black")
        )} +
      xlab(xlab) +
      ylab("Probability") +
      theme_bw() +
      list(guides(colour = "none", fill = "none"))
    
  }
  
  dummy <- data.frame(x = rep(c(0),4), 
                      y = rep(c(0),4),
                      lower = rep(c(0),4),
                      upper = rep(c(0),4),
                      lag = c("0","1","2","3"))
  
  plot_list[[idx]] <- wrap_plots(panel_list, ncol = 3) +
    plot_layout(guides = "collect", axis_titles = "collect") +
    # Add legend manually
    new_scale_colour() +
    new_scale_fill() +
    geom_point(aes(colour = lag), data = dummy, alpha = 0) +
    scale_colour_manual(name = element_blank(),
                        breaks = c("0", "1", "2", "3"),
                        values = setNames(pal, c("0","1","2","3")))
  
  ggsave(paste("plots/",
               INCLUDE_CROSSTERMS,
               "_",
               CV_OR_RI,
               "_partial_dependence_B_Q", idx, ".png", sep=""),
         plot = plot_list[[idx]],
         width = 12,
         height = 6)
  
}

