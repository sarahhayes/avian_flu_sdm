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
  preds[[idx]] <- pred_layer[[1]]
  names(preds[[idx]]) <- paste0("Q",idx)
  
}

png("plots/all_poster.png", width = 7, height = 7,
    units = "in", res = 330)
spplot(stack(preds),
       col.regions = viridis_pal()(100),
       at = seq(0,1,0.01),
       cex = 0.8)
grid::grid.text("Probability", x=grid::unit(0.98, "npc"), y=grid::unit(0.50, "npc"), rot=-90)
dev.off()




