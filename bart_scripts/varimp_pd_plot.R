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

#PATH_TO_OUTPUTS <- "H:/Working/avian_flu_sdm/output/"

library(tidyverse)
library(dbarts)
library(patchwork)

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
    data.frame(Q = paste0("Q",idx)) %>%
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
    mutate(lag = replace_na(lag, "ns")) %>%
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
  mutate(orig_var = var,
         var = gsub("anatinae" , "abundance: anatinae", var),
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
         orig_var = gsub("_2022", "", orig_var),
         orig_var = gsub("_lag.$", "", orig_var),
         var = gsub("_2022", "", var),
         var = gsub("_lag.$", "", var),
         Q = case_when(Q == "Q1" ~ "Q1 (Jan - Mar)",
                       Q == "Q2" ~ "Q2 (Apr - Jun)",
                       Q == "Q3" ~ "Q3 (Jul - Sep)",
                       Q == "Q4" ~ "Q4 (Oct - Dec)"
         )
  )

top_var <- df %>% group_by(var, orig_var) %>% summarise(g_mean = mean(mean), count = n()) # Data frame to order variables by mean varimp or n times selected
df <- df %>% mutate(var = fct_relevel(var, ordered_var),
                    upper = mean + sd, lower = mean - sd) %>% mutate(var = fct_relevel(var, ordered_var)) %>% arrange(var)

fig_varimp_A <- ggplot(df, aes(x = mean, y = var, xmin = lower, xmax = upper, color = lag)) + 
  geom_errorbar(width=0, position=position_dodge(0.6)) + 
  geom_point(position=position_dodge(0.6)) +
  geom_hline(yintercept=seq(1.5, nrow(df)-0.5, 1), lwd=.5, colour="grey75") +  # add custom gridlines
  geom_hline(yintercept=c(6.5,11.5,12.5), lwd=1.5, colour="grey75") + # separate bioclimatic, topographic, livestock, host ecological variables - must set manually each recalculation!
  scale_x_continuous() +
  scale_y_discrete() +
  scale_colour_manual(name = element_blank(),
                      breaks = c("ns", "0", "1", "2", "3"),
                      values = c("black", pal)) +
  coord_flip() +
  facet_wrap(~ Q, nrow=4, strip.position="right") +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(size = 11, angle = 30, hjust = 1), 
        legend.title=element_blank(), 
        legend.position=c(-0.325, 0.975),
        legend.key.size=unit(0.8, "lines"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.background = element_rect(colour = "grey50"),
        axis.title.x=element_blank(),
        plot.margin = margin(3, 3, 3, 20),
        panel.grid.major.y = element_blank(),
        strip.text.y = element_text(size = 10)) +
  xlab("Relative variable importance")

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_variable_importance_quarterly_A.png", sep=""),
       plot = fig_varimp_A,
       width = 10,
       height = 6)

# Calc and bind pd across all quarters for given variables by name. All you need to do is set names of variables of interest (or don't do anything and let it plot all of them) :)

# Define variables of interest
varorigs <- top_var %>% arrange(-g_mean) %>% pull(orig_var) %>% .[1:6]
varlabs <- top_var %>% arrange(-g_mean) %>% pull(var) %>% .[1:6]

pd_summ <- lapply(1:4, function(x) list()) # initialise empty list of lists

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
  
  if (ALL_VARS_FOR_PD){
    vars <- sdm$varcount %>% as.data.frame %>% names
  }else{
    vars <- varorigs
  }
  
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_model_with_vs_A_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
  for(j in 1:length(vars)){ # for variables for interest..
    
    for(k in grep(vars[j], (sdm$fit$data@x %>% attr("term.labels")))){ # capture any matching (i.e., across all lags)
      
      # Check if the predictor variable is binary
      if (all(sdm$fit$data@x[, k] %in% 0:1)){
        raw <- sdm$fit$data@x[, k]
        lev <- list(c(0,1))
        pd <- pdbart(sdm, xind = k, levs = lev, pl = FALSE)
        
      } else {
        raw <- sdm$fit$data@x[, k]
        lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
        pd <- pdbart(sdm, xind = k, levs = lev, pl = FALSE)
        
      }
      
      pd_var <- data.frame(var = sdm$fit$data@x %>% attr("term.labels") %>% .[k], 
                           x = unlist(lev), 
                           y = pd$fd[[1]] %>% apply(., 2, median) %>% pnorm,
                           upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025) %>% pnorm,
                           lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975) %>% pnorm,
                           Q = paste0("Q",idx)) %>%
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
        mutate(lag = replace_na(lag, "ns")) %>%
        select(-varq)
      
      pd_summ[[idx]] <- append(pd_summ[[idx]], list(pd_var))
      
    }
  }
}  

if (ALL_VARS_FOR_PD){
  
  for(idx in 1:4){
    
    # Make individual panels of a partial dependence plot, automatically assigning line or point depending on whether variable is continuous or binary
    panel_list <- list() # initialise empty lists
    plot_list <- list()
    
    pd_df <- pd_summ[[idx]] %>%
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
      mutate(Q = as.factor(Q),
             lag = as.factor(lag)) %>%
      arrange(var)
    
    for(j in 1:length(unique(pd_df$var))){
      
      df <- pd_df %>%
        filter(var == unique(pd_df$var)[j])
      
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
                              values = setNames(c(pal, "black"), c("0","1","2","3","ns")),
                              drop = FALSE),
            scale_colour_manual(name = element_blank(),
                                breaks = c("0", "1", "2", "3"),
                                values = setNames(c(pal, "black"), c("0","1","2","3","ns")),
                                drop = FALSE)
          )} +
        {if(all(df$x %in% 0:1)) 
          list(
            geom_pointrange(aes(x = factor(x)), size = 0.7, position = position_dodge(width = 0.5), alpha = 0.8, color = "black", fill = "black")
          )} +
        xlab(xlab) +
        ylab("Probability") +
        theme_bw()
      
    }
    
    plot_list[[idx]] <- wrap_plots(panel_list, ncol = 3) +
      plot_layout(guides = "collect", axis_titles = "collect")
    
    ggsave(paste("plots/",
                 INCLUDE_CROSSTERMS,
                 "_",
                 CV_OR_RI,
                 "_partial_dependence_A_Q", idx, ".png", sep=""),
           plot = plot_list[[idx]],
           width = 12,
           height = 6)
    
  }
  
}else{
  
  pd_df <- pd_summ %>%
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
    mutate(Q = as.factor(Q),
           lag = as.factor(lag)) %>%
    arrange(var)
  
  panel_list <- list()
  
  for(j in 1:length(unique(pd_df$var))){
    
    xlab <- unique(pd_df$var)[j]
    
    panel_list[[j]] <- pd_df %>% 
      filter(var == unique(pd_df$var)[j]) %>%
      ggplot(aes(x = x, y = y, ymin = lower, ymax = upper, fill = Q, color = Q, lty = lag)) +
      {if(!(all(pd_df$x %in% 0:1)))
        list(
          geom_line(lwd = 0.8, alpha = 0.7),
          scale_x_continuous(expand = c(0, 0))
        )} +
      {if(all(pd_df$x %in% 0:1)) 
        list(
          geom_pointrange(aes(x = factor(x)), size = 0.7, position = position_dodge(width = 0.5))
        )} +
      scale_colour_manual("",
                          breaks = c("Q1", "Q2", "Q3", "Q4"),
                          values = pal,
                          drop = FALSE) +
      scale_fill_manual("",
                        breaks = c("Q1", "Q2", "Q3", "Q4"),
                        values = pal,
                        drop = FALSE) +
      scale_linetype_manual(values = c("solid", "longdash", "twodash", "dotted"),
                            drop = FALSE) +
      xlab(xlab) +
      ylab("Probability") +
      theme_bw()
  }
  
  fig_pd_chosen <- wrap_plots(panel_list, ncol = 2) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = 'bottom')
    
  ggsave(paste("plots/",
               INCLUDE_CROSSTERMS,
               "_",
               CV_OR_RI,
               "_partial_dependence_top_mean_A.png", sep=""),
         plot = fig_pd_chosen,
         width = 9,
         height = 7)
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
    data.frame(Q = paste0("Q",idx)) %>%
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
    mutate(lag = replace_na(lag, "ns")) %>%
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
  mutate(orig_var = var,
         var = gsub("anatinae" , "abundance: anatinae", var),
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
         orig_var = gsub("_2022", "", orig_var),
         orig_var = gsub("_lag.$", "", orig_var),
         var = gsub("_2022", "", var),
         var = gsub("_lag.$", "", var),
         Q = case_when(Q == "Q1" ~ "Q1 (Jan - Mar)",
                       Q == "Q2" ~ "Q2 (Apr - Jun)",
                       Q == "Q3" ~ "Q3 (Jul - Sep)",
                       Q == "Q4" ~ "Q4 (Oct - Dec)"
         )
  )

top_var <- df %>% group_by(var, orig_var) %>% summarise(g_mean = mean(mean), count = n()) # Data frame to order variables by mean varimp or n times selected
df <- df %>% mutate(var = fct_relevel(var, ordered_var),
                    upper = mean + sd, lower = mean - sd) %>% mutate(var = fct_relevel(var, ordered_var)) %>% arrange(var)

fig_varimp_B <- ggplot(df, aes(x = mean, y = var, xmin = lower, xmax = upper, color = lag)) + 
  geom_errorbar(width=0, position=position_dodge(0.6)) + 
  geom_point(position=position_dodge(0.6)) +
  geom_hline(yintercept=seq(1.5, nrow(df)-0.5, 1), lwd=.5, colour="grey75") +  # add custom gridlines
  geom_hline(yintercept=c(7.5,13.5), lwd=1.5, colour="grey75") + # separate bioclimatic, topographic, livestock, host ecological variables - must set manually each recalculation!
  scale_x_continuous() +
  scale_y_discrete() +
  scale_colour_manual(name = element_blank(),
                      breaks = c("ns", "0", "1", "2", "3"),
                      values = c("black", pal)) +
  coord_flip() +
  facet_wrap(~ Q, nrow=4, strip.position="right") +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(size = 11, angle = 30, hjust = 1), 
        legend.title=element_blank(), 
        legend.position=c(-0.325, 0.975),
        legend.key.size=unit(0.8, "lines"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.background = element_rect(colour = "grey50"),
        axis.title.x=element_blank(),
        plot.margin = margin(3, 3, 3, 20),
        panel.grid.major.y = element_blank(),
        strip.text.y = element_text(size = 10)) +
  xlab("Relative variable importance")

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_variable_importance_quarterly_B.png", sep=""),
       plot = fig_varimp_B,
       width = 10,
       height = 6)


# Calc and bind pd across all quarters for given variables by name. All you need to do is set names of variables of interest (or don't do anything and let it plot all of them) :)

# Define variables of interest
varorigs <- top_var %>% arrange(-g_mean) %>% pull(orig_var) %>% .[1:6]
varlabs <- top_var %>% arrange(-g_mean) %>% pull(var) %>% .[1:6]

pd_summ <- lapply(1:4, function(x) list()) # initialise empty list of lists

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
  
  if (ALL_VARS_FOR_PD){
    vars <- sdm$varcount %>% as.data.frame %>% names
  }else{
    vars <- varorigs
  }
  
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_model_with_vs_B_Q",
                    idx,
                    ".rds",
                    sep = ""))
  
  for(j in 1:length(vars)){ # for variables for interest..
    
    for(k in grep(vars[j], (sdm$fit$data@x %>% attr("term.labels")))){ # capture any matching (i.e., across all lags)
      
      # Check if the predictor variable is binary
      if (all(sdm$fit$data@x[, k] %in% 0:1)){
        raw <- sdm$fit$data@x[, k]
        lev <- list(c(0,1))
        pd <- pdbart(sdm, xind = k, levs = lev, pl = FALSE)
        
      } else {
        raw <- sdm$fit$data@x[, k]
        lev <- list(seq(min(raw), max(raw), ((max(raw) - min(raw))/15)))
        pd <- pdbart(sdm, xind = k, levs = lev, pl = FALSE)
        
      }
      
      pd_var <- data.frame(var = sdm$fit$data@x %>% attr("term.labels") %>% .[k], 
                           x = unlist(lev), 
                           y = pd$fd[[1]] %>% apply(., 2, median) %>% pnorm,
                           upper = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.025) %>% pnorm,
                           lower = pd$fd[[1]] %>% apply(.,  2, quantile, probs = 0.975) %>% pnorm,
                           Q = paste0("Q",idx)) %>%
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
        mutate(lag = replace_na(lag, "ns")) %>%
        select(-varq)
      
      pd_summ[[idx]] <- append(pd_summ[[idx]], list(pd_var))
      
    }
  }
}  

if (ALL_VARS_FOR_PD){
  
  for(idx in 1:4){
    
    # Make individual panels of a partial dependence plot, automatically assigning line or point depending on whether variable is continuous or binary
    panel_list <- list() # initialise empty lists
    plot_list <- list()
    
    pd_df <- pd_summ[[idx]] %>%
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
      mutate(Q = as.factor(Q),
             lag = as.factor(lag)) %>%
      arrange(var)
    
    for(j in 1:length(unique(pd_df$var))){
      
      df <- pd_df %>%
        filter(var == unique(pd_df$var)[j])
      
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
                              values = setNames(c(pal, "black"), c("0","1","2","3","ns")),
                              drop = FALSE),
            scale_colour_manual(name = element_blank(),
                                breaks = c("0", "1", "2", "3"),
                                values = setNames(c(pal, "black"), c("0","1","2","3","ns")),
                                drop = FALSE)
          )} +
        {if(all(df$x %in% 0:1)) 
          list(
            geom_pointrange(aes(x = factor(x)), size = 0.7, position = position_dodge(width = 0.5), alpha = 0.8, color = "black", fill = "black")
          )} +
        xlab(xlab) +
        ylab("Probability") +
        theme_bw()
      
    }
    
    plot_list[[idx]] <- wrap_plots(panel_list, ncol = 3) +
      plot_layout(guides = "collect", axis_titles = "collect")
    
    ggsave(paste("plots/",
                 INCLUDE_CROSSTERMS,
                 "_",
                 CV_OR_RI,
                 "_partial_dependence_B_Q", idx, ".png", sep=""),
           plot = plot_list[[idx]],
           width = 12,
           height = 6)
    
  }
  
}else{
  
  pd_df <- pd_summ %>%
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
    mutate(Q = as.factor(Q),
           lag = as.factor(lag)) %>%
    arrange(var)
  
  panel_list <- list()
  
  for(j in 1:length(unique(pd_df$var))){
    
    xlab <- unique(pd_df$var)[j]
    
    panel_list[[j]] <- pd_df %>% 
      filter(var == unique(pd_df$var)[j]) %>%
      ggplot(aes(x = x, y = y, ymin = lower, ymax = upper, fill = Q, color = Q, lty = lag)) +
      {if(!(all(pd_df$x %in% 0:1)))
        list(
          geom_line(lwd = 0.8, alpha = 0.7),
          scale_x_continuous(expand = c(0, 0))
        )} +
      {if(all(pd_df$x %in% 0:1)) 
        list(
          geom_pointrange(aes(x = factor(x)), size = 0.7, position = position_dodge(width = 0.5))
        )} +
      scale_colour_manual("",
                          breaks = c("Q1", "Q2", "Q3", "Q4"),
                          values = pal,
                          drop = FALSE) +
      scale_fill_manual("",
                        breaks = c("Q1", "Q2", "Q3", "Q4"),
                        values = pal,
                        drop = FALSE) +
      scale_linetype_manual(values = c("solid", "longdash", "twodash", "dotted"),
                            drop = FALSE) +
      xlab(xlab) +
      ylab("Probability") +
      theme_bw()
  }
  
  fig_pd_chosen <- wrap_plots(panel_list, ncol = 2) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = 'bottom')
  
  ggsave(paste("plots/",
               INCLUDE_CROSSTERMS,
               "_",
               CV_OR_RI,
               "_partial_dependence_top_mean_B.png", sep=""),
         plot = fig_pd_chosen,
         width = 9,
         height = 7)
}



# Combination plots

fig_varimp_combi <- wrap_plots(list(fig_varimp_A, fig_varimp_B), ncol = 2) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = 'bottom',
        plot.tag = element_text(size = 14))

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_variable_importance_quarterly_combi.png", sep=""),
       plot = fig_varimp_combi,
       width = 16,
       height = 8)
