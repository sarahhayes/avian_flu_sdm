# avian_flu_sdm

Scripts for avian influenza species distribution modelling project

For the species-level preliminary analysis and ecological covariate raster construction, users will need to have version 2.2021.3 of the ebirdst package installed. To install it, run the following R code:
```r
packageurl <- "https://cran.r-project.org/src/contrib/Archive/ebirdst/ebirdst_2.2021.3.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
```

**Folders:**

data - subfolders within this folder for storage of all the data needed for the model 
  
ebird - data and vignette explaining how to use EBird. May not be needed

ebirdst - vignettes for ebirdst which is likely to be the package that will be used. Scripts contain links to web pages with more details

output - folder for storing useful outputs from scripts

plots - folder for storage of plot pngs/pdfs

z_old_or_unused scripts - storage for my scripts which are no longer needed or are me trying things out 

**Files:**

bvbrc_first_look - initial exploration of bv-brc data

elton_traits - reading in data from elton traits site

empresi_first_look - initial exploration of EMPRES-i data

europe_plots - experimenting with plotting case data on map of Europe
