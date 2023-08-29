# In this script we plot risk projections using the outputs from
# fit_avian_flu_model.R.

library(embarcadero)
library(raster)
library(terra)
library(viridis)
set.seed(12345)

for (idx in 1:4){
  load(file = paste("output/fitted-BART-models/prediction_Q",
                    idx,
                    ".rds",
                    sep = ""))
       png(filename = paste("plots/q", idx, "_poster.png", sep = ""))
  plot(pred_layer[[1]],
       box = FALSE,
       axes = FALSE,
       col = viridis_pal()(100),
       main = paste("Q", idx, sep = ""),
       legend.args = list(text = "Probability", side = 2))
  dev.off()
  
}
