# Example Using geoSIMEX 
# Using geoSIMEX with data from Uganda's third administrative division and dev. finance data from AidData
# AidData, REU
# Last Updated: 10/7/2016

# --- --- --- --- --- --- --- --- --- --- #
##### Loading Packages and Data #####
# --- --- --- --- --- --- --- --- --- --- #
library(raster)

source("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/R/geoSIMEX.R")
uga.adm3.df <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/merge_uga_adm3.csv")
uga.aiddata <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/UgandaAMP_GeocodedResearchRelease_Level1_v1.3/data/level_1a.csv")

# --- --- --- --- --- --- --- --- --- --- #
##### Prepping Data for Analysis #####
# --- --- --- --- --- --- --- --- --- --- #

# 1. get ADM data on AidData data to match uga.adm3.df data
# 2. do a check to make sure sure all names in aiddata data are in uga.adm3.df data

# --- --- --- --- --- --- --- --- --- --- #
##### Analysis #####
# --- --- --- --- --- --- --- --- --- --- #


