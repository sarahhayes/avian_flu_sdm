# In this script we investigate patterns in the Callaghan et al. population size
# estimates to see how using median might differ from using upper 95% boundary

pop_size_df <- read.csv("callaghan_pop_estimates.csv")
no_pops <- nrow(pop_size_df)

med_order <- order(pop_size_df$Abundance.estimate)
med_vect <- pop_size_df$Abundance.estimate[med_order]
u95_vect <- pop_size_df$X95..Upper.CI[med_order]
l95_vect <- pop_size_df$X95..Lower.CI[med_order]

u95_diffs <- u95_vect[2:no_pops] - u95_vect[1:no_pops-1]
median_diffs <- med_vect[2:no_pops] - med_vect[1:no_pops-1]

diff_scale <- u95_diffs[which(median_diffs>0)] / 
  median_diffs[which(median_diffs>0)]
log_diff_scale <- sign(diff_scale) * log(abs(diff_scale))

plot(sort(log_diff_scale),
     ylab = "log(D(U95)/D(median))")
plot(med_vect[2:no_pops][which(median_diffs>0)],
     log_diff_scale,
     log = "x",
     xlab = "Median abundance estimate",
     ylab = "log(D(U95)/D(median))")

plot(sort(u95_vect / med_vect),
     log = "y",
     ylab = "U95/median")

plot(med_vect,
     u95_vect / med_vect,
     log = "xy",
     xlab = "Median abundance estimate",
     ylab = "D(U95)/D(median)")

plot(med_vect,
     u95_vect,
     log = "xy",
     xlab = "Median abundance estimate",
     ylab = "U95")

int_width <- u95_vect - l95_vect
plot(med_vect,
     int_width / med_vect,
     log = "xy",
     xlab = "Median abundance estimate",
     ylab = "Width of 95% CI")
plot(med_vect,
     int_width,
     log = "xy",
     xlab = "Median abundance estimate",
     ylab = "Width of 95% CI")

idx <- which(med_vect>0)
ratio <- int_width[idx] / med_vect[idx]
hist(ratio,
     breaks = as.integer(max(int_width[idx] / med_vect[idx])),
     xlab = "Width of CI / median",
     main = paste("Mean=",
                  mean(ratio),
                  ", std=",
                  sd(ratio)))
