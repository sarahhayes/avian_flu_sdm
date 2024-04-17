# In this script we plot risk projections using the outputs from
# fit_avian_flu_model.R.

rm(list = ls())

# Optional command line arguments, must be passed as strings:
args <- commandArgs(trailingOnly = T)
if (length(args)<4){
  INCLUDE_CROSSTERMS <- "no-crossterms" # Set to "no-crossterms" to do model without crossterms or "with-crossterms" to do model with crossterms
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

library(embarcadero)
library(raster)
library(terra)
library(viridis)
set.seed(12345)

preds <- list()


for (idx in 1:4){
  # Slightly hacky way to load in predictions and give it an arbitrary name.
  # This is needed because readRDS appears to be deprecated and the load
  # function brings in the original variable name. If you do x<-load(x.rds) then
  # x just gives you the variable name; when we pipe with get() we can get the
  # value of the thing with that variable name, provided it's been loaded in.
  load(file = paste(PATH_TO_OUTPUTS,
                    "fitted-BART-models-",
                    INCLUDE_CROSSTERMS,
                    "/",
                    CV_OR_RI,
                    "_predictions_A_Q",
                    idx,
                    ".rds",
                    sep = ""))
  pred_layer <- load(file = paste(PATH_TO_OUTPUTS,
                                     "fitted-BART-models-",
                                     INCLUDE_CROSSTERMS,
                                     "/",
                                     CV_OR_RI,
                                     "_predictions_A_Q",
                                     idx,
                                     ".rds",
                                     sep = "")) %>% get
  preds[[idx]] <- pred_layer
  names(preds[[idx]]) <- c(paste0("Q",idx),
                           paste0("Q",idx,"_2.5th_percentile"),
                           paste0("Q",idx,"_97.5th_percentile"))
  
}

png("plots/projections.png", width = 7, height = 7,
    units = "in", res = 330)
spplot(stack(preds[[1]][[1]],
             preds[[2]][[1]], 
             preds[[3]][[1]], 
             preds[[4]][[1]]),
       col.regions = viridis_pal()(100),
       at = seq(0,1,0.01),
       cex = 0.8)
grid::grid.text("Probability", x=grid::unit(0.98, "npc"), y=grid::unit(0.50, "npc"), rot=-90)
dev.off()


png("plots/q1_uncertainty.png", width = 21, height = 7,
    units = "in", res = 330)
spplot(stack(preds[[1]][[2]],
             preds[[1]][[1]],
             preds[[1]][[3]]),
       col.regions = viridis_pal()(100),
       at = seq(0,1,0.01),
       cex = 0.8)
grid::grid.text("Probability", x=grid::unit(0.98, "npc"), y=grid::unit(0.50, "npc"), rot=-90)
dev.off()


png("plots/q2_uncertainty.png", width = 21, height = 7,
    units = "in", res = 330)
spplot(stack(preds[[2]][[2]],
             preds[[2]][[1]],
             preds[[2]][[3]]),
       col.regions = viridis_pal()(100),
       at = seq(0,1,0.01),
       cex = 0.8)
grid::grid.text("Probability", x=grid::unit(0.98, "npc"), y=grid::unit(0.50, "npc"), rot=-90)
dev.off()


png("plots/q3_uncertainty.png", width = 21, height = 7,
    units = "in", res = 330)
spplot(stack(preds[[3]][[2]],
             preds[[3]][[1]],
             preds[[3]][[3]]),
       col.regions = viridis_pal()(100),
       at = seq(0,1,0.01),
       cex = 0.8)
grid::grid.text("Probability", x=grid::unit(0.98, "npc"), y=grid::unit(0.50, "npc"), rot=-90)
dev.off()


png("plots/q4_uncertainty.png", width = 21, height = 7,
    units = "in", res = 330)
spplot(stack(preds[[4]][[2]],
             preds[[4]][[1]],
             preds[[4]][[3]]),
       col.regions = viridis_pal()(100),
       at = seq(0,1,0.01),
       cex = 0.8)
grid::grid.text("Probability", x=grid::unit(0.98, "npc"), y=grid::unit(0.50, "npc"), rot=-90)
dev.off()


pdf("plots/all_uncertainty.pdf", paper="a4", width = 8, height = 11.3)
spplot(stack(preds[[1]][[2]],
             preds[[1]][[1]],
             preds[[1]][[3]],
             preds[[2]][[2]], 
             preds[[2]][[1]], 
             preds[[2]][[3]], 
             preds[[3]][[2]], 
             preds[[3]][[1]], 
             preds[[3]][[3]], 
             preds[[4]][[2]],
             preds[[4]][[1]],
             preds[[4]][[3]]),
       col.regions = viridis_pal()(100),
       at = seq(0,1,0.01),
       cex = 0.8)
grid::grid.text("Probability", x=grid::unit(0.98, "npc"), y=grid::unit(0.50, "npc"), rot=-90)
dev.off()

################################################################################
# Plot case data

pos_data <- read.csv(paste(PATH_TO_OUTPUTS,
                           "AI_S2_SDM_storage/Avian flu data/pos_points_proj_area_all_sources_duplicates_removed.csv",
                           sep = ""))

zipmap <- terra::vect(x = "data/gis_europe/CNTR_RG_03M_2020_4326.shp.zip",
                      layer = "CNTR_RG_03M_2020_4326")
crs <- "epsg:3035"
euro_ext <- extent(preds[[1]])

# change projection and extent. 
# using quite a generous extent whilst plotting as looking at where to set the boundaries
euro_map <- terra::project(x = zipmap, y = crs) 
euro_map_crop <- terra::crop(euro_map, euro_ext)

pts_pos <- terra::vect(pos_data, geom=c("X", "Y"),
                       crs =  "+proj=longlat +ellps=WGS84 +datum=WGS84")

w <- 7
h <- w * (euro_ext@ymax - euro_ext@ymin)/(euro_ext@xmax - euro_ext@xmin)
png("plots/poster_positives.png", width = w, height = h,
    units = "in", res = 330)
plot(euro_map_crop,
     col = "white",
     background = "azure2",
     axes = FALSE,
     buffer = FALSE,
     xmin = euro_ext@xmin,
     mar = c(0, 0, 0, 0))
plot(pts_pos, add = T, col = "red", pch = 16, cex = .3)
dev.off()