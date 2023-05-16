## Trialing a different package to get the MODIS data

rm(list = ls())

#install.packages("MODIStsp")
# install.packages("leaflet")
# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("shinyFiles")
# install.packages("shinyalert")
# install.packages("rappdirs")
# install.packages("shinyjs")
# install.packages("leafem")
# install.packages("mapedit")
# install.packages("magrittr")

library(leaflet)
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyalert)
library(rappdirs)
library(shinyjs)
library(leafem)
library(mapedit)
library(magrittr)
library(MODIStsp)

MODIStsp()
