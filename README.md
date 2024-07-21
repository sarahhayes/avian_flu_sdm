# avian_flu_sdm

Scripts for avian influenza species distribution modelling project

For the species-level preliminary analysis and ecological covariate raster construction, users will need to have version 2.2021.3 of the ebirdst package installed. To install it, run the following R code:
```r
packageurl <- "https://cran.r-project.org/src/contrib/Archive/ebirdst/ebirdst_2.2021.3.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
```

**Folders:**

data - subfolders within this folder for storage of all the data needed for the model 
  
output - folder for storing useful outputs from scripts

plots - folder for storage of plot pngs/pdfs

**Files:**


elton_traits - reading in data from elton traits site
