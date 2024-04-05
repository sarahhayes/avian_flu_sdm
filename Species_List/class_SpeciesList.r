# 04/01/2024
# Script to process bird species in a more reproducible format. 

rm(list = ls())

library(tidyverse)
library(ebirdst)

# Class for reading, cleaning and storing species list
## CSV that is read in comes from: https://science.ebird.org/en/status-and-trends/species?regionCode=eu
## These are the birds for which ebird had status and trends data for region specified in filter 
## The list of species, after filtering to specified region, is copied and saved as csv file.  
setClass("SpeciesList",
	slots = list(
	# Data frame to store species list
	df = "data.frame"
	# Region for species list
	,region = "character"
	# Total number of species for given region
	,n = "integer"
))

# clean df so that species list has scientific and common names and codes
setGeneric("clean",function(obj) standardGeneric("clean"))
setMethod("clean","SpeciesList",
	function(obj){
		## Identify scientific family names (end with "dae"), family names in common english appear in the row after.
		n = grep("dae",obj@df$species)
		## Ensure no species names identified with family names (species names are binomial so have a space)
		m = grep(" ", obj@df$species[n])
		if(length(m)){
			n = n[-m]
		}
		## Removing family (scientific and common) names 
		## drop argument ensures obj@df remains a dataframe
		obj@df <- obj@df[-c(n,n+1),,drop=F]
		
		## Removing order name (entries ending with FORMES) and entries with word SPECIES in it 
		obj@df <- obj@df %>% 
			filter(!str_detect(species, "SPECIES")) %>% 
			filter(!str_detect(species, "FORMES")) 

		## The copy-paste from website creates duplicates which need to be removed
		obj@df <- unique(obj@df,by = c("species"))
			
		## Common (Scientific) name appears in odd (even) rows
		## creating new column name to identify name type
		obj@df$name_type <- rep(c("common","scientific"),nrow(obj@df)/2)
		## Providing unique id for each species name pair (1,1,2,2,3,3,...) to match scientific and common names
		obj@df$id <- sort(rep(1:(nrow(obj@df)/2),2))
		## Reshaping dataframe so scientific and common names are in separate columns
		obj@df <- spread(obj@df,key=name_type,val=species)
		## Removing id column
		obj@df <- obj@df[,-1,drop=F]
		
		## Get code species names from ebird 
		obj@df$common_code <- ebirdst::get_species(obj@df$common)
		obj@df$scientific_code <- ebirdst::get_species(obj@df$scientific)
		
		return(obj)
	}
)


# Checking if scientific and common codes match
# Checking if number of species match, i.e. correct number of species is obtained after raw data is processed
setGeneric("check",function(obj) standardGeneric("check"))
setMethod("check","SpeciesList",
	function(obj){
		cat(paste(
			"Species List for:"
			,obj@region
			,"\nNumber of species match:"
			,obj@n == nrow(obj@df)
			,"\nScientific and common codes match:"
			,0==length(which(obj@df$common_code != obj@df$scientific_code))
			,"\nNumber of species:"
			,obj@n
			,"\nNumber of species that match ebirdst_runs database:"
			,length(which(!is.na(obj@df$common_code)))
			,"\n\n"
		))
	}
)

# Save cleaned species list as csv
setGeneric("save",function(obj) standardGeneric("save"))
setMethod("save","SpeciesList",
	function(obj){
		## only saving those species which have a valid code in ebirdst_runs database
		temp <- subset(obj@df,!is.na(common_code))
		## Selecting columns required and renaming
		temp <- temp[,c("common_code","scientific","common")]
		colnames(temp) = c("code","scientific_name","common_name")
		## Setting scienfitic and common names to lower case
		temp$scientific_name <- tolower(temp$scientific_name)
		temp$common_name <- tolower(temp$common_name)
		## Saving as csv in Process_Species_List using name specified in region
		write.csv(temp,paste0("Processed_Species_List/",obj@region,".csv"),row.names=FALSE)
	}
)

# Constructor function
SpeciesList <- function(region){
	obj <- new("SpeciesList"
	## Assumes that data to be read is stored in folder "Raw_Species_List" and has name specified in region
	,df = read.csv(paste0("Raw_Species_List/",region,".txt"),sep="\t",header=F)
	,region = region
	)
	## df has one column so changing name to "species" for ease
	colnames(obj@df) <- "species"
	
	## Entries with SPECIES HIDE give the number of species for a family, so can use that to count total number of species
	n = grep("SPECIES HIDE",obj@df$species)
	obj@n <- sum(
		as.integer(sub(" SPECIES HIDE","",obj@df$species[n]))
	)
	
	obj <- clean(obj)
	
	return(obj)
}









                   