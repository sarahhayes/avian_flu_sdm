### Looking at using CRU data in NetCDF format. 

## some tips taken from https://www.youtube.com/watch?v=jWy_jGZo2oc

rm(list = ls())

# install.packages("ncdf4")

library(ncdf4)
library(raster)
# library(tidyverse) this masks some of the raster commands so if want to use 
# tidyverse will specify inline. 

fn <- ("data/variables/cru_ts4.06.2011.2020.tmn.dat.nc") # set the file name

#nc_file_tmp <- nc_open("data/variables/cru_ts4.06.2011.2020.tmn.dat.nc") # can do with pathwa
nc_file_tmp <- nc_open(fn) # or with file name

print(nc_file_tmp)

# make in to a raster brick

rasbrick <- raster::brick(fn) # can open directly as a raster brick
rasbrick

#calculating some descriptive stats using the raster brick https://www.youtube.com/watch?v=ZNXRexFIQ5w

rasmean <- mean(rasbrick)

# to get standard dev we have to do a slightly different methods

rassd <- calc(rasbrick, fun = sd, na.rm = T)


library(viridis)
library(viridisLite)
library(ggplot2)
