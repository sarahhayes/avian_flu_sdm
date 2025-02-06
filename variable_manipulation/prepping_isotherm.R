## Preparing data on the zero degree isotherm for inclusion within the model. 
## Data are from copernicus
## https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
## The data are from midday on each day

rm(list = ls())

library(terra)

isotherm <- terra::rast("data/variables/isotherm/isocline_raw.grib")
isotherm

## Have a little look at some of the data
isotherm[[1]][1:10]
isotherm[[1]][500000:500010]
plot(isotherm[[1]])

iso1<- isotherm[[1]]
iso1
plot(iso1)
isotherm[[2]]

# Split into the different quarters as the raw data are daily 

## calendar quarters 
# first_quart_stack <- isotherm[[1:90]] 
# first_quart_stack
# second_quart_stack <- isotherm[[91:181]]
# second_quart_stack
# third_quart_stack <- isotherm[[182:273]]
# third_quart_stack
# fourth_quart_stack <- isotherm[[274:365]]
# fourth_quart_stack

### ecological quarters
first_quart_stack <- c(isotherm[[333:365]], isotherm[[1:59]])
first_quart_stack
second_quart_stack <- isotherm[[60:157]]
second_quart_stack
third_quart_stack <- isotherm[[158:221]]
third_quart_stack
fourth_quart_stack <- isotherm[[222:333]]
fourth_quart_stack


min_q1 <- min(first_quart_stack)
min_q1
first_quart_stack
first_quart_stack[[1]][500000:500010]
min_q1[500000:500010]
par(mfrow = c(1,2))
plot(first_quart_stack[[1]])
plot(min_q1)


### Which aspects of these data are we interested in? 
### Paper on isocline talks about AI outbreaks around the time the isocline occurs. 
### Minimum value in the data is 0. So if we select the minimum, we will know if the isocline occurred in
### that area and that quarter
### Could also try and look at the number of days that was zero in that pixel

## first look at minimum as most straightforward
min_q2 <- min(second_quart_stack)
min_q3 <- min(third_quart_stack)
min_q4 <- min(fourth_quart_stack)

# Still global extent , WGS 84 projection and 0.25 degree resolution so need to change to the resolution we want. 

# Transforming the data
# blank_3035 <- terra::rast("output/euro_rast.tif") # 1k res
blank_3035 <- terra::rast("output/euro_rast_10k.tif") # 10k res
blank_3035

# will use the minimum method so we know if we have the isotherm in that pixel
min_q1_prj <- terra::project(min_q1, blank_3035, method = "min")
min_q2_prj <- terra::project(min_q2, blank_3035, method = "min")
min_q3_prj <- terra::project(min_q3, blank_3035, method = "min")
min_q4_prj <- terra::project(min_q4, blank_3035, method = "min")

#pdf("plots/min_isotherm.pdf", height = 7, width = 7)
pdf("plots/min_isotherm_eco.pdf", height = 7, width = 7)
par(mfrow = c(2,2))
plot(min_q1_prj); plot(min_q2_prj); plot(min_q3_prj); plot(min_q4_prj)
mtext("Minimum isotherm", side = 3, line = - 2, outer = TRUE)
dev.off()


plot(min_q1_prj) # bit unexpected that most of the UK doesn't seem to have a zero in this quarter, but these data are at midday. 


#terra::writeRaster(min_q1_prj, "variable_manipulation/variable_outputs/isotherm_min_q1.tif")
#terra::writeRaster(min_q2_prj, "variable_manipulation/variable_outputs/isotherm_min_q2.tif")
#terra::writeRaster(min_q3_prj, "variable_manipulation/variable_outputs/isotherm_min_q3.tif")
#terra::writeRaster(min_q4_prj, "variable_manipulation/variable_outputs/isotherm_min_q4.tif")

# terra::writeRaster(min_q1_prj, "variable_manipulation/variable_outputs/isotherm_min_q1_eco_quarts.tif")
# terra::writeRaster(min_q2_prj, "variable_manipulation/variable_outputs/isotherm_min_q2_eco_quarts.tif")
# terra::writeRaster(min_q3_prj, "variable_manipulation/variable_outputs/isotherm_min_q3_eco_quarts.tif")
# terra::writeRaster(min_q4_prj, "variable_manipulation/variable_outputs/isotherm_min_q4_eco_quarts.tif")


# Next is to try to count number of times the value of 0 appears in the quarter. This will tell 
# us the number of days where the isocline was at ground level
# Some useful (but not perfect) info available here https://gis.stackexchange.com/questions/455349/counting-number-of-observations-per-pixel-in-a-rasterstack
# Alternative would be to do a for loop and add the counts, but suspect this would take quite a long time.... 

?terra::freq
freq(first_quart_stack)
freq(first_quart_stack, value = 0) # this just gives overall counts and isn't for each cell

zeros_q1 <- app(first_quart_stack, function(x) sum(x == 0))
zeros_q1
plot(zeros_q1)

# Or could use those within 1m of the ground? 
zero_one_q1 <- app(first_quart_stack, function(x) sum(x <=1))
zero_one_q1
plot(zero_one_q1)

par(mfrow = c(1,2))
plot(zeros_q1); plot(zero_one_q1)
# perhaps not looking obvious because we have the poles and that is skewing the same? 
# Let's re-project and review

zero_number_q1_prj <- terra::project(zeros_q1, blank_3035, method = "max")
zero_ones_q1_prj <- terra::project(zero_one_q1, blank_3035, method = "max")

plot(zero_number_q1_prj); plot(zero_ones_q1_prj)
# This seems feasible. It is picking out the mountain ranges etc.

# I think will use the less than 1m cutoff
zero_one_q2 <- app(second_quart_stack, function(x) sum(x <=1))
zero_ones_q2_prj <- terra::project(zero_one_q2, blank_3035, method = "max")
zero_ones_q2_prj

zero_one_q3 <- app(third_quart_stack, function(x) sum(x <=1))
zero_ones_q3_prj <- terra::project(zero_one_q3, blank_3035, method = "max")
zero_ones_q3_prj

zero_one_q4 <- app(fourth_quart_stack, function(x) sum(x <=1))
zero_ones_q4_prj <- terra::project(zero_one_q4, blank_3035, method = "max")
zero_ones_q4_prj

#pdf("plots/isotherm_number_days_below1_midday.pdf", height = 7, width = 7)
pdf("plots/isotherm_number_days_below1_midday_eco_quarts.pdf", height = 7, width = 7)
par(mfrow = c(2,2))
plot(zero_ones_q1_prj, main = "Q1")
plot(zero_ones_q2_prj, main = "Q2")
plot(zero_ones_q3_prj, main = "Q3")
plot(zero_ones_q4_prj, main = "Q4")
mtext("Number of days where isotherm is <1m", side = 3, outer = T, line = -1)
dev.off()


#terra::writeRaster(zero_ones_q1_prj, "variable_manipulation/variable_outputs/isotherm_midday_days_below1_q1.tif")
#terra::writeRaster(zero_ones_q2_prj, "variable_manipulation/variable_outputs/isotherm_midday_days_below1_q2.tif")
#terra::writeRaster(zero_ones_q3_prj, "variable_manipulation/variable_outputs/isotherm_midday_days_below1_q3.tif")
#terra::writeRaster(zero_ones_q4_prj, "variable_manipulation/variable_outputs/isotherm_midday_days_below1_q4.tif")

terra::writeRaster(zero_ones_q1_prj, "variable_manipulation/variable_outputs/isotherm_midday_days_below1_q1_eco_quarts.tif")
terra::writeRaster(zero_ones_q2_prj, "variable_manipulation/variable_outputs/isotherm_midday_days_below1_q2_eco_quarts.tif")
terra::writeRaster(zero_ones_q3_prj, "variable_manipulation/variable_outputs/isotherm_midday_days_below1_q3_eco_quarts.tif")
terra::writeRaster(zero_ones_q4_prj, "variable_manipulation/variable_outputs/isotherm_midday_days_below1_q4_eco_quarts.tif")

#### START HERE FOR EDITS

## Next we want to look at the mean monthly data

iso_mean <- terra::rast("data/variables/isotherm/isocline_2022_monthly_mean.grib")
iso_mean
## This one has 12 layers as monthly data

mean_first_quart_stack <- iso_mean[[1:3]]
mean_second_quart_stack <- iso_mean[[4:6]]
mean_third_quart_stack <- iso_mean[[7:9]]
mean_fourth_quart_stack <- iso_mean[[10:12]]

mean_iso_q1 <- mean(mean_first_quart_stack)
mean_iso_q1
mean_iso_q1_prj <- terra::project(mean_iso_q1, blank_3035, method = "bilinear")
mean_iso_q1_prj
plot(mean_iso_q1_prj)

mean_iso_q2 <- mean(mean_second_quart_stack)
mean_iso_q2_prj <- terra::project(mean_iso_q2, blank_3035, method = "bilinear")
mean_iso_q2_prj
plot(mean_iso_q2_prj)

mean_iso_q3 <- mean(mean_third_quart_stack)
mean_iso_q3_prj <- terra::project(mean_iso_q3, blank_3035, method = "bilinear")
mean_iso_q3_prj
plot(mean_iso_q3_prj)

mean_iso_q4 <- mean(mean_fourth_quart_stack)
mean_iso_q4_prj <- terra::project(mean_iso_q4, blank_3035, method = "bilinear")
mean_iso_q4_prj
plot(mean_iso_q4_prj)

# terra::writeRaster(mean_iso_q1_prj, "variable_manipulation/variable_outputs/isotherm_mean_q1.tif")
# terra::writeRaster(mean_iso_q2_prj, "variable_manipulation/variable_outputs/isotherm_mean_q2.tif")
# terra::writeRaster(mean_iso_q3_prj, "variable_manipulation/variable_outputs/isotherm_mean_q3.tif")
# terra::writeRaster(mean_iso_q4_prj, "variable_manipulation/variable_outputs/isotherm_mean_q4.tif")


pdf("plots/isotherm_mean.pdf", height = 7, width = 7)
par(mfrow = c(2,2))
plot(mean_iso_q1_prj, main = "Q1")
plot(mean_iso_q2_prj, main = "Q2")
plot(mean_iso_q3_prj, main = "Q3")
plot(mean_iso_q4_prj, main = "Q4")
mtext("Mean isotherm (m)", side = 3, outer = T, line = -1)
dev.off()

## Using the ecological quarters it is not as easy to use the monthly means. 
## However, given we have the daily data we could use the means from that

daily_mean_first_quart <- mean(first_quart_stack)
daily_mean_first_quart
daily_mean_iso_q1_prj <- terra::project(daily_mean_first_quart, blank_3035, method = "bilinear")
daily_mean_iso_q1_prj
plot(daily_mean_iso_q1_prj)

daily_mean_second_quart <- mean(second_quart_stack)
daily_mean_second_quart
daily_mean_iso_q2_prj <- terra::project(daily_mean_second_quart, blank_3035, method = "bilinear")
daily_mean_iso_q2_prj
plot(daily_mean_iso_q2_prj)

daily_mean_third_quart <- mean(third_quart_stack)
daily_mean_third_quart
daily_mean_iso_q3_prj <- terra::project(daily_mean_third_quart, blank_3035, method = "bilinear")
daily_mean_iso_q3_prj
plot(daily_mean_iso_q3_prj)

daily_mean_fourth_quart <- mean(fourth_quart_stack)
daily_mean_fourth_quart
daily_mean_iso_q4_prj <- terra::project(daily_mean_fourth_quart, blank_3035, method = "bilinear")
daily_mean_iso_q4_prj
plot(daily_mean_iso_q4_prj)

# terra::writeRaster(daily_mean_iso_q1_prj, "variable_manipulation/variable_outputs/isotherm_mean_q1_eco_quarts.tif")
# terra::writeRaster(daily_mean_iso_q2_prj, "variable_manipulation/variable_outputs/isotherm_mean_q2_eco_quarts.tif")
# terra::writeRaster(daily_mean_iso_q3_prj, "variable_manipulation/variable_outputs/isotherm_mean_q3_eco_quarts.tif")
# terra::writeRaster(daily_mean_iso_q4_prj, "variable_manipulation/variable_outputs/isotherm_mean_q4_eco_quarts.tif")
# 


###############
### Now repeat the daily process for data from midnight (c.f. midday) 

iso_midnight <- terra::rast("data/variables/isotherm/isotherm_daily_midnight.grib")
iso_midnight

# Split into the different quarters as the raw data are daily 
midnight_q1_stack <- iso_midnight[[1:90]]
midnight_q1_stack
midnight_q2_stack <- iso_midnight[[91:181]]
midnight_q2_stack
midnight_q3_stack <- iso_midnight[[182:273]]
midnight_q3_stack
midnight_q4_stack <- iso_midnight[[274:365]]
midnight_q4_stack


## Calculate the minimums
min_midnight_q1 <- min(midnight_q1_stack)
min_midnight_q1

## first look at minimum as most straightforward
min_midnight_q2 <- min(midnight_q2_stack)
min_midnight_q3 <- min(midnight_q3_stack)
min_midnight_q4 <- min(midnight_q4_stack)

# Still global extent , WGS 84 projection and 0.25 degree resolution so need to change to the resolution we want. 

# will use the minimum method so we know if we have the isotherm in that pixel
min_midnight_q1_prj <- terra::project(min_midnight_q1, blank_3035, method = "min")
min_midnight_q2_prj <- terra::project(min_midnight_q2, blank_3035, method = "min")
min_midnight_q3_prj <- terra::project(min_midnight_q3, blank_3035, method = "min")
min_midnight_q4_prj <- terra::project(min_midnight_q4, blank_3035, method = "min")

pdf("plots/min_isotherm_midnight.pdf", height = 7, width = 7)
par(mfrow = c(2,2))
plot(min_midnight_q1_prj); plot(min_midnight_q2_prj); plot(min_midnight_q3_prj); plot(min_midnight_q4_prj)
mtext("Minimum isotherm - midnight", side = 3, line = - 2, outer = TRUE)
dev.off()


#terra::writeRaster(min_midnight_q1_prj, "variable_manipulation/variable_outputs/isotherm_min_midnight_q1.tif")
#terra::writeRaster(min_midnight_q2_prj, "variable_manipulation/variable_outputs/isotherm_min_midnight_q2.tif")
#terra::writeRaster(min_midnight_q3_prj, "variable_manipulation/variable_outputs/isotherm_min_midnight_q3.tif")
#terra::writeRaster(min_midnight_q4_prj, "variable_manipulation/variable_outputs/isotherm_min_midnight_q4.tif")


# Next is to try to count number of times the value of 0 appears in the quarter. This will tell 
# us the number of days where the isocline was at ground level

# Using those within 1m of the ground? 
zero_one_q1_midnight <- app(midnight_q1_stack, function(x) sum(x <=1))
zero_one_q1_midnight
plot(zero_one_q1_midnight)

zero_ones_q1_prj_midnight <- terra::project(zero_one_q1_midnight, blank_3035, method = "max")

plot(zero_ones_q1_prj_midnight)
# This seems feasible. It is picking out the mountain ranges etc.

# I think will use the less than 1m cutoff
zero_one_q2_midnight <- app(midnight_q2_stack, function(x) sum(x <=1))
zero_ones_q2_prj_midnight <- terra::project(zero_one_q2_midnight, blank_3035, method = "max")
zero_ones_q2_prj_midnight

zero_one_q3_midnight <- app(midnight_q3_stack, function(x) sum(x <=1))
zero_ones_q3_prj_midnight <- terra::project(zero_one_q3_midnight, blank_3035, method = "max")
zero_ones_q3_prj_midnight

zero_one_q4_midnight <- app(midnight_q4_stack, function(x) sum(x <=1))
zero_ones_q4_prj_midnight <- terra::project(zero_one_q4_midnight, blank_3035, method = "max")
zero_ones_q4_prj_midnight

pdf("plots/isotherm_number_days_below1_midnight.pdf", height = 7, width = 7)
par(mfrow = c(2,2))
plot(zero_ones_q1_prj_midnight, main = "Q1")
plot(zero_ones_q2_prj_midnight, main = "Q2")
plot(zero_ones_q3_prj_midnight, main = "Q3")
plot(zero_ones_q4_prj_midnight, main = "Q4")
mtext("Number of days where isotherm is <1m", side = 3, outer = T, line = -1)
dev.off()


#terra::writeRaster(zero_ones_q1_prj_midnight, "variable_manipulation/variable_outputs/isotherm_midnight_days_below1_q1.tif")
#terra::writeRaster(zero_ones_q2_prj_midnight, "variable_manipulation/variable_outputs/isotherm_midnight_days_below1_q2.tif")
#terra::writeRaster(zero_ones_q3_prj_midnight, "variable_manipulation/variable_outputs/isotherm_midnight_days_below1_q3.tif")
#terra::writeRaster(zero_ones_q4_prj_midnight, "variable_manipulation/variable_outputs/isotherm_midnight_days_below1_q4.tif")


