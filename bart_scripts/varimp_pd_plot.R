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

### Set plot variable order (should match tables) and names

plot_labels <- data.frame(
  var = c(
    "mean_relative_humidity",
    "mean_diff",
    "mean_temp",
    "variation_in_quarterly_mean_temp",
    "mean_tmax",
    "mean_tmin",
    "isotherm_mean",
    "isotherm_midday_days_below1",
    "mean_prec",
    "elevation_min",
    "elevation_max",
    "elevation_mode",
    "elevation_diff",
    "ndvi",
    "Water_bodies",
    "Evergreen_Needleleaf_Forests",
    "Evergreen_Broadleaf_Forests", 
    "Deciduous_Needleleaf_Forests",
    "Deciduous_Broadleaf_Forests",
    "Mixed_Forests",
    "Closed_Shrublands",
    "Open_Shrublands",
    "Woody_Savannas",
    "Savannas",
    "Grasslands",
    "Permanent_Wetlands",
    "Croplands",
    "Urban_and_Built-up_Lands",
    "Cropland/Natural_Vegetation_Mosaics",
    "Non-Vegetated_Lands",
    "Unclassified",
    "dist_to_coast_km",
    "dist_to_water",
    "chicken_density_2010",
    "duck_density_2010",
    "anatinae",
    "anserinae",
    "ardeidae",
    "arenaria_calidris",
    "aythyini",
    "laridae",
    "around_surf",
    "below_surf",
    "plant",
    "scav",
    "vend",
    "cong",
    "migr",
    "species_richness",
    "host_dist"
  ),
  label = c(
    "mean relative humidity",
    "monthly temp. range",
    "mean temp.",
    "quarterly temp. variation",
    "max. temp.",
    "min. temp.",
    "mean 0° isotherm",
    "n days 0° isotherm <1m",
    "total precipitation",
    "min. elevation",
    "max. elevation",
    "modal elevation",
    "elevation diff.",
    "vegetation index",
    "water body",
    "evergreen needleleaf forest",
    "evergreen broadleaf forest",
    "deciduous needleleaf forest",
    "deciduous broadleaf forest",
    "mixed forest",
    "closed shrubland",
    "open shrubland",
    "woody savanna",
    "savanna",
    "grassland",
    "permanent wetland",
    "cropland",
    "urban and built-up land",
    "cropland/natural vegetation mosaic",
    "non-vegetated land",
    "unclassified land",
    "dist. to coast",
    "dist. to inland water",
    "chicken density",
    "dom. duck density", 
    "Anatinae", 
    "Anserinae", 
    "Ardeidae",
    "Arenaria/Calidris",
    "Aythyini",
    "Laridae", 
    "water surface feeding",
    "below water surface feeding",
    "plant diet",
    "scavenging diet",
    "endotherm vert. diet",
    "congregative birds",
    "migratory birds",
    "species richness",
    "below phylo threshold"
  ),
  label_units = c(
    "mean relative humidity (%)",
    "monthly temp. range (°C)",
    "mean temp. (°C)",
    "quarterly temp. variation (°C)",
    "max. temp. (°C)",
    "min. temp. (°C)",
    "mean 0° isotherm (m)",
    "n days 0° isotherm <1m",
    "total precipitation (mm)",
    "min. elevation (m)",
    "max. elevation (m)",
    "modal elevation (m)",
    "elevation diff. (m)",
    "vegetation index",
    "water body",
    "evergreen needleleaf forest",
    "evergreen broadleaf forest",
    "deciduous needleleaf forest",
    "deciduous broadleaf forest",
    "mixed forest",
    "closed shrubland",
    "open shrubland",
    "woody savanna",
    "savanna",
    "grassland",
    "permanent wetland",
    "cropland",
    "urban and built-up land",
    "cropland/natural vegetation mosaic",
    "non-vegetated land",
    "unclassified land",
    "dist. to coast (km)",
    "dist. to inland water (m)",
    "chicken density (headcount)",
    "dom. duck density (headcount)", 
    "Anatinae (log(abundance index))", 
    "Anserinae (log(abundance index))", 
    "Ardeidae (log(abundance index))",
    "Arenaria/Calidris (log(abundance index))",
    "Aythyini (log(abundance index))",
    "Laridae (log(abundance index))", 
    "water surface feeding (log(abundance index))",
    "below water surface feeding (log(abundance index))",
    "plant diet (log(abundance index))",
    "scavenging diet (log(abundance index))",
    "endotherm vert. diet (log(abundance index))",
    "congregative birds (log(abundance index))",
    "migratory birds (log(abundance index))",
    "species richness",
    "below phylo threshold (log(abundance index))"
  )
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
  mutate(var = gsub("_2022", "", var),
         var = gsub("_lag.$", "", var),
         Q = case_when(Q == "Q1" ~ "Q1 (Jan - Mar)",
                       Q == "Q2" ~ "Q2 (Apr - Jun)",
                       Q == "Q3" ~ "Q3 (Jul - Sep)",
                       Q == "Q4" ~ "Q4 (Oct - Dec)"
         )
  ) %>% 
  left_join(plot_labels)

top_var <- df %>% group_by(var, label) %>% summarise(g_mean = mean(mean), count = n()) # Data frame to order variables by mean varimp or n times selected
df <- df %>% mutate(label = fct_relevel(label, plot_labels$label) %>% suppressWarnings(),
                    upper = mean + sd, lower = mean - sd) %>% mutate(label = fct_relevel(label, plot_labels$label) %>% suppressWarnings()) %>% arrange(label)

fig_varimp_A <- ggplot(df, aes(x = mean, y = label, xmin = lower, xmax = upper, color = lag)) + 
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
selected_vars <- top_var %>% arrange(-g_mean) %>% pull(var) %>% .[1:6]

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
    vars <- selected_vars
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
      mutate(var = gsub("_2022", "", var),
             var = gsub("_lag.$", "", var)
      ) %>% 
      left_join(plot_labels) %>%
      mutate(label_units = as.factor(label_units),
             label_units = fct_relevel(label_units, plot_labels$label_units) %>% suppressWarnings()) %>% 
      mutate(Q = as.factor(Q),
             lag = as.factor(lag)) %>%
      arrange(label_units)
    
    for(j in 1:length(unique(pd_df$var))){
      
      df <- pd_df %>%
        filter(var == unique(pd_df$var)[j])
      
      xlab <- unique(df$label_units)
      
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
        {if(all(df$var %in% c("anatinae", "anserinae", "ardeidae", "arenaria_calidris", "aythyini", 
                              "laridae", "around_surf", "below_surf", "plant", "scav", "vend", 
                              "cong", "migr", "host_dist"))) 
          list(
            scale_x_continuous(trans = "log10", expand = c(0,0))
          )} +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        xlab(xlab) +
        ylab("Probability") +
        theme_bw()
      
    }
    
    plot_list[[idx]] <- wrap_plots(panel_list, ncol = 3) +
      plot_annotation(tag_levels = 'i') +
      plot_layout(guides = "collect")
    
    ggsave(paste("plots/",
                 INCLUDE_CROSSTERMS,
                 "_",
                 CV_OR_RI,
                 "_partial_dependence_A_Q", idx, ".png", sep=""),
           plot = plot_list[[idx]],
           width = 12,
           height = length(plot_list[[idx]])*2/3)
    
  }
  
}else{
  
  pd_df <- pd_summ %>%
    bind_rows %>%
    mutate(var = gsub("_2022", "", var),
           var = gsub("_lag.$", "", var)
    ) %>% 
    left_join(plot_labels) %>%
    mutate(label_units = as.factor(label_units),
           label_units = fct_relevel(label_units, plot_labels$label_units) %>% suppressWarnings()) %>%
    mutate(Q = as.factor(Q),
           lag = as.factor(lag)) %>%
    arrange(label_units)
  
  panel_list <- list()
  
  for(j in 1:length(unique(pd_df$var))){
    
    xlab <- unique(pd_df$label_units)[j]
    
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
      {if(all(df$x %in% c("anatinae", "anserinae", "ardeidae", "arenaria_calidris", "aythyini", 
                          "laridae", "around_surf", "below_surf", "plant", "scav", "vend", 
                          "cong", "migr", "host_dist"))) 
        list(
          scale_x_continuous(trans = "log10", expand = c(0,0))
        )} +
      scale_colour_manual("",
                          breaks = c("Q1", "Q2", "Q3", "Q4", "ns"),
                          values = c(pal, "black"),
                          drop = FALSE) +
      scale_fill_manual("",
                        breaks = c("Q1", "Q2", "Q3", "Q4", "ns"),
                        values = c(pal, "black"),
                        drop = FALSE) +
      scale_linetype_manual("",
                            breaks = c("ns",0,1,2,3),
                            values = c("solid", "solid", "longdash", "twodash", "dotted"),
                            drop = FALSE) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      xlab(xlab) +
      ylab("Probability") +
      theme_bw()
  }
  
  fig_pd_chosen_A <- wrap_plots(panel_list, ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom', legend.key.size = unit(1.05, "cm"),
          plot.tag = element_text(size = 12))
  
  ggsave(paste("plots/",
               INCLUDE_CROSSTERMS,
               "_",
               CV_OR_RI,
               "_partial_dependence_top_mean_A.png", sep=""),
         plot = fig_pd_chosen_A,
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
  mutate(var = gsub("_2022", "", var),
         var = gsub("_lag.$", "", var),
         Q = case_when(Q == "Q1" ~ "Q1 (Jan - Mar)",
                       Q == "Q2" ~ "Q2 (Apr - Jun)",
                       Q == "Q3" ~ "Q3 (Jul - Sep)",
                       Q == "Q4" ~ "Q4 (Oct - Dec)"
         )
  ) %>% 
  left_join(plot_labels)

top_var <- df %>% group_by(var, label) %>% summarise(g_mean = mean(mean), count = n()) # Data frame to order variables by mean varimp or n times selected
df <- df %>% mutate(label = fct_relevel(label, plot_labels$label) %>% suppressWarnings(),
                    upper = mean + sd, lower = mean - sd) %>% mutate(label = fct_relevel(label, plot_labels$label) %>% suppressWarnings()) %>% arrange(label)

fig_varimp_B <- ggplot(df, aes(x = mean, y = label, xmin = lower, xmax = upper, color = lag)) + 
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
selected_vars <- top_var %>% arrange(-g_mean) %>% pull(var) %>% .[1:6]

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
    vars <- selected_vars
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
      mutate(var = gsub("_2022", "", var),
             var = gsub("_lag.$", "", var)
      ) %>% 
      left_join(plot_labels) %>%
      mutate(label_units = as.factor(label_units),
             label_units = fct_relevel(label_units, plot_labels$label_units) %>% suppressWarnings()) %>% 
      mutate(Q = as.factor(Q),
             lag = as.factor(lag)) %>%
      arrange(label_units)
    
    for(j in 1:length(unique(pd_df$var))){
      
      df <- pd_df %>%
        filter(var == unique(pd_df$var)[j])
      
      xlab <- unique(df$label_units)
      
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
        {if(all(df$var %in% c("anatinae", "anserinae", "ardeidae", "arenaria_calidris", "aythyini", 
                              "laridae", "around_surf", "below_surf", "plant", "scav", "vend", 
                              "cong", "migr", "host_dist"))) 
          list(
            scale_x_continuous(trans = "log10", expand = c(0,0))
          )} +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        xlab(xlab) +
        ylab("Probability") +
        theme_bw()
      
    }
    
    plot_list[[idx]] <- wrap_plots(panel_list, ncol = 3) +
      plot_annotation(tag_levels = 'i') +
      plot_layout(guides = "collect")
    
    ggsave(paste("plots/",
                 INCLUDE_CROSSTERMS,
                 "_",
                 CV_OR_RI,
                 "_partial_dependence_B_Q", idx, ".png", sep=""),
           plot = plot_list[[idx]],
           width = 12,
           height = length(plot_list[[idx]])*2/3)
    
  }
  
}else{
  
  pd_df <- pd_summ %>%
    bind_rows %>%
    mutate(var = gsub("_2022", "", var),
           var = gsub("_lag.$", "", var)
    ) %>% 
    left_join(plot_labels) %>%
    mutate(label_units = as.factor(label_units),
           label_units = fct_relevel(label_units, plot_labels$label_units) %>% suppressWarnings()) %>%
    mutate(Q = as.factor(Q),
           lag = as.factor(lag)) %>%
    arrange(label_units)
  
  panel_list <- list()
  
  for(j in 1:length(unique(pd_df$var))){
    
    xlab <- unique(pd_df$label_units)[j]
    
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
      {if(all(df$x %in% c("anatinae", "anserinae", "ardeidae", "arenaria_calidris", "aythyini", 
                          "laridae", "around_surf", "below_surf", "plant", "scav", "vend", 
                          "cong", "migr", "host_dist"))) 
        list(
          scale_x_continuous(trans = "log10", expand = c(0,0))
        )} +
      scale_colour_manual("",
                          breaks = c("Q1", "Q2", "Q3", "Q4", "ns"),
                          values = c(pal, "black"),
                          drop = FALSE) +
      scale_fill_manual("",
                        breaks = c("Q1", "Q2", "Q3", "Q4", "ns"),
                        values = c(pal, "black"),
                        drop = FALSE) +
      scale_linetype_manual("",
                            breaks = c("ns",0,1,2,3),
                            values = c("solid", "solid", "longdash", "twodash", "dotted"),
                            drop = FALSE) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      xlab(xlab) +
      ylab("Probability") +
      theme_bw()
  }
  
  fig_pd_chosen_B <- wrap_plots(panel_list, ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom', legend.key.size = unit(1.05, "cm"),
          plot.tag = element_text(size = 12))
  
  ggsave(paste("plots/",
               INCLUDE_CROSSTERMS,
               "_",
               CV_OR_RI,
               "_partial_dependence_top_mean_B.png", sep=""),
         plot = fig_pd_chosen_B,
         width = 9,
         height = 7)
}



# Combination plots

fig_varimp_combi <- wrap_plots(list(fig_varimp_A, fig_varimp_B), ncol = 2) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = 'bottom',
        plot.tag = element_text(size = 14))

# Ensure same y ranges for all plots; adapted from https://stackoverflow.com/a/60018428
flipped_y_range <- c(ggplot_build(fig_varimp_combi[[1]])$layout$panel_scales_x[[1]]$range$range,
                     ggplot_build(fig_varimp_combi[[2]])$layout$panel_scales_x[[1]]$range$range)

fig_varimp_combi <- fig_varimp_combi & xlim(min(flipped_y_range), max(flipped_y_range))

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_variable_importance_quarterly_combi.png", sep=""),
       plot = fig_varimp_combi,
       width = 16,
       height = 8)

fig_pd_combi <- (fig_pd_chosen_A)/(fig_pd_chosen_B) +
  plot_annotation(tag_levels = list(c("A.i", "ii", "iii", "iv", "v", "vi", "B.i", "ii", "iii", "iv", "v", "vi"))) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = 'bottom',
        plot.tag = element_text(size = 12))

ggsave(paste("plots/",
             INCLUDE_CROSSTERMS,
             "_",
             CV_OR_RI,
             "_partial_dependence_top_mean_combi.png", sep=""),
       plot = fig_pd_combi,
       width = 9,
       height = 14.5)
