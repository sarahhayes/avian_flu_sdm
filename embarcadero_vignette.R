# embarcadero vignette. 
# 01/08/2023
# stored for me at file:///C:/Users/hayes/OneDrive%20-%20Nexus365/Documents/GitHub/embarcadero/vignettes/virtualbart.html

rm(list = ls())

# devtools::install_github('cjcarlson/embarcadero')
# install.packages("remotes")
# remotes::install_github("cran/RandomFieldsUtils")
 remotes::install_github("cran/RandomFields")
 remotes::install_github("ropensci/NLMR")

library(embarcadero, quietly = T)
library(dismo, quietly=T)
library(NLMR, quietly = T)
library(virtualspecies, quietly = T)
set.seed(42)

# Generate some imaginary variables
onelandscape <- function(x) {NLMR::nlm_gaussianfield(nrow = 150,
                                                     ncol = 150,
                                                     rescale = FALSE)}
climate <- stack(lapply(c(1:8), onelandscape))
xnames <- c('x1','x2','x3','x4','x5','x6','x7','x8')
names(climate) <- xnames

plot(climate[[1]],main='An imaginary variable')
# These variables are in a raster stack

# Generate the species' climatic niche

random.sp <- generateRandomSp(climate[[1:4]], 
                              # ^ These are the informative predictors
                              approach="pca",
                              relations='gaussian',
                              species.prevalence=0.5,
                              realistic.sp = TRUE,
                              PA.method='threshold')
class(random.sp)

# Generate some presences, and some absences, with imperfect detection
# sampleOccurrences is from the virtualspecies package 
sp.points <- sampleOccurrences(random.sp,
                               n=250,
                               type = 'presence-absence',
                               detection.probability = 0.9)

# Extract the associated climate values

occ <- SpatialPoints(sp.points$sample.points[,c('x','y')])
occ.df <- cbind(sp.points$sample.points,
                raster::extract(climate, occ))

# Finally, let's drop the long-lats and the "Real" ground truthed presence-absence 
# values, and just leave behind an "Observed" and the climate data

occ.df <- occ.df[,-c(1:3)]

# Building a basic BART model
sdm <- bart(y.train=occ.df[,'Observed'],
            x.train=occ.df[,xnames],
            keeptrees = TRUE) # It's very important this is set to TRUE
