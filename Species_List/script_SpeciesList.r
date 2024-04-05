## Create species list using class_SpeciesList
## Raw species list obtained from https://science.ebird.org/en/status-and-trends/species and saved in folder Raw_Data

source("class_SpeciesList.r")

for(i in c(
	"north_america"
	,"central_america"
	,"south_america"
	,"asia"
	,"australia_and_territories"
	,"europe"
)){
	obj<-SpeciesList(i)
	check(obj)
	save(obj)
}

