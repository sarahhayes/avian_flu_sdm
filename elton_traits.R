
## Reading in the Elton Traits from the website

rm(list = ls())

library(data.table)
library(tidyverse)

myfile <- fread("https://www.esapubs.org/archive/ecol/E095/178/BirdFuncDat.txt")
