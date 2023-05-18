# change the colours
# quantiles of non-zero values
v <- values(crop_abd)
v <- v[!is.na(v) & v > 0]
range(v)
bins <- quantile(v, seq(0, 1, by = 0.1)) # splits into the 10% quantiles
bins <- c(0, bins) # add a bin for 0
bins <- c(1e-10, 1e-6, 1e-3, 1e-2, 1.5e-2, 1e-1, 1.5e-1, 1,100)

# status and trends palette
pal <- abundance_palette(length(bins) - 2)
# add a color for zero
pal <- c("#e6e6e6", pal)
pal <- c("grey", pal)

#alternative colour scheme
pal <- abundance_palette(length(bins) - 1)
pal <- abundance_palette(length(bins))
# map using the quantile bins
plot(crop_abd[[1]], breaks = bins, col = pal, axes = FALSE)


crop_abd
winter <- mean(crop_abd[[1:13]])
autumn <- mean(crop_abd[[27:39]])
plot(crop_abd[[26]], axes = F, col = "#DCDCDC", legend = F)
plot(winter, breaks = bins, col = pal, axes = FALSE, add = T)

plot(crop_abd[[26]], axes = F, col = "grey", legend = F)
plot(autumn, breaks = bins, col = pal, axes = F, add = T)

v_wint <- values(winter)
v_wint <- v_wint[!is.na(v_wint) & v_wint > 0]
range(v_wint)

v_aut <- values(autumn)
v_aut <- v_aut[!is.na(v_aut) & v_aut > 0]
range(v_aut)


bins_wint <- quantile(v_wint, seq(0, 1, by = 0.1)) # splits into the 10% quantiles
bins_wint <- c(0, bins_wint) # add a bin for 0
#bins <- c(0,0.1,1,2,3,4,5,10,50,100)

# status and trends palette
pal_wint <- abundance_palette(length(bins_wint) - 2)

# add a color for zero
pal_wint <- c("#e6e6e6", pal_wint)
#pal <- c("grey", pal)

# alternative colour scheme
# pal <- abundance_palette(length(bins) - 1)

# map using the quantile bins
plot(winter, breaks = bins_wint, col = pal_wint, axes = FALSE)
plot(autumn, breaks = bins_wint, col = pal_wint, axes = FALSE)


set_bins <-  c(0,0.1,0.5,1,2,3,4,5,10,50,100)
pal_set <- abundance_palette(length(set_bins) - 2)
# add a color for zero
pal_set <- c("#e6e6e6", pal_set)

# map using the quantile bins
save
plot(winter, breaks = set_bins, col = pal_set, axes = FALSE)
plot(autumn, breaks = set_bins, col = pal_set, axes = FALSE)


png("plots/mall_wint.png", height=480, width=480)
plot(winter, breaks = set_bins, col = pal_set, axes = FALSE)
dev.off()

png("plots/mall_aut.png", height=480, width=480)
plot(autumn, breaks = set_bins, col = pal_set, axes = FALSE)
dev.off()
