# In this script we benchmark the performance of the project function as the
# number of layers increases.

rm(list = ls())

library(terra)
library(rnaturalearth)

increment_size = 20
no_increments = 25

# Load in the Europe projection
crs <- "epsg:3035"
wh_europe <- ne_countries(continent = "europe",
                          returnclass = "sf") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

# Createa single-layer test object to check it's working
test_rast <- rast(nrows=180, ncols=360, xmin=-180, xmax=180, ymin=-90, ymax=90)
values(test_rast) <- runif(ncell(test_rast))
test_prj <- project(test_rast, crs)
par(mfrow=c(1,2))
plot(test_rast, main="Original projection")
plot(test_prj, main=crs)

par(mfrow=c(1,1))
for (i in 1:1) {
  for (i in 1:(increment_size-1)) {
    new_lyr <- rast(nrows=180,
                    ncols=360,
                    xmin=-180,
                    xmax=180,
                    ymin=-90,
                    ymax=90)
    values(new_lyr) <- runif(ncell(new_lyr))
    test_rast <- c(test_rast, new_lyr)
  }
  start.time <- Sys.time()
  test_prj <- project(test_rast, crs)
  end.time <- Sys.time()
  benchmark_vals <- as.numeric(difftime(end.time, start.time, units="secs"))
}

loop.start <- Sys.time()
for (i in 2:no_increments) {
  for (j in 1:increment_size) {
    new_lyr <- rast(nrows=180,
                    ncols=360,
                    xmin=-180,
                    xmax=180,
                    ymin=-90,
                    ymax=90)
    values(new_lyr) <- runif(ncell(new_lyr))
    test_rast <- c(test_rast, new_lyr)
  }
  start.time <- Sys.time()
  test_prj <- project(test_rast, crs)
  end.time <- Sys.time()
  benchmark_vals <- c(benchmark_vals, as.numeric(difftime(end.time, start.time, units="secs")))
  cat("Iteration ", i, " of ", no_increments, "completed. \n")
  if (i>5){
    m_est = (benchmark_vals[i] - benchmark_vals[i-5])/5
  }
  else {
    m_est = (benchmark_vals[i] - benchmark_vals[i-1])
  }
  time_remaining = (no_increments - i) * benchmark_vals[i] + m_est * sum(1:(no_increments-i))
  cat(end.time - loop.start,
      " seconds elapsed since start, estimated ",
      time_remaining,
      " remaining.\n")
}

plot(increment_size*(1:no_increments),
     benchmark_vals,
     main="Benchmarking of project function",
     xlab="Number of layers",
     ylab="Execution time")
