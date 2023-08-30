# In this script we plot risk projections using the outputs from
# fit_avian_flu_model.R.

rm(list = ls())

library(embarcadero)
library(raster)
library(terra)
library(viridis)
set.seed(12345)

preds <- list()


for (idx in 1:4){
  load(file = paste("output/fitted-BART-models/prediction_Q",
                    idx,
                    ".rds",
                    sep = ""))
  preds[[idx]] <- pred_layer
  names(preds[[idx]]) <- c(paste0("Q",idx),
                           paste0("Q",idx,"_2.5th_percentile"),
                           paste0("Q",idx,"_97.5th_percentile"))
  
}

png("plots/all_poster.png", width = 7, height = 7,
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



