# Trialling setting up the Modis data as per the instructions here:

# https://cornelllabofornithology.github.io/ebird-best-practices/covariates.html

rm(list = ls())

# installing and updating the necessary packages.
#install.packages("remotes")
#remotes::install_github("mstrimas/ebppackages")

library(cli)
# MODIS package has been archived so need to download from the archive as per instructions
# here https://stackoverflow.com/questions/24194409/how-do-i-install-a-package-that-has-been-archived-from-cran

# Download package tarball from CRAN archive

#url <- "http://cran.r-project.org/src/contrib/Archive/MODIS/MODIS_1.2.11.tar.gz"
#pkgFile <- "MODIS_1.2.11.tar.gz"
#download.file(url = url, destfile = pkgFile)

# Expand the zip file using whatever system functions are preferred

# look at the DESCRIPTION file in the expanded package directory

# Install dependencies list in the DESCRIPTION file

# install.packages(c("mapdata", "bitops", "mapedit", "ptw"))

# Install package
#install.packages(pkgs=pkgFile, type="source", repos=NULL)

# Delete package tarball
#unlink(pkgFile)


library(sf)
library(raster)
library(MODIS)
library(exactextractr)
library(viridis)
library(tidyverse)

MODIS::EarthdataLogin(usr = "sarahhayes", pwd = "EDLTigtogs43!")

MODIS:::checkTools("GDAL")

# getTile() # this allows you to select which tiles you want

eur_tile <- getTile() # temporary whilst working out how it works
eur_tile@tile # shows which ones you've picked

## Next bit shows how to get data over multiple years
# earliest year required
# begin_year <- format(min(ebird$observation_date), "%Y.01.01") # from vignette
begin_year <- "2017.01.01"
# end date for ebird data
# end_year <- format(max(ebird$observation_date), "%Y.12.31")
end_year <- "2018.12.31"

MODISoptions(check_earthdata_login = TRUE)

# download tiles and combine into a single raster for each year
tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01", 
                #extent = bcr %>% st_buffer(dist = 10000), 
                extent = eur_tile@extent,
                begin = begin_year, end = end_year, 
                outDirPath = "data", job = "modis",
                MODISserverOrder = "LPDAAC") %>% 
  pluck("MCD12Q1.006") %>% 
  unlist()

# rename tifs to have more descriptive names
new_names <- format(as.Date(names(tifs)), "%Y") %>% 
  sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
  file.path(dirname(tifs), .)
file.rename(tifs, new_names)