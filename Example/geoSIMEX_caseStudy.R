install.packages("geoSIMEX")
install.packages("sp")
install.packages("jsonlite")
install.packages("raster")
library(simexlibrary)
library(sp)
library(jsonlite)
library(raster)







source("https://github.com/itpir/geoSIMEX/blob/master/R/geoSIMEX.R")

##### Upload CSV and json #####
json.uganda <- fromJSON("~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX/Example/summary.json")
uga.adm <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/merge_uga_adm3.csv")

uga.aiddata <- subset.aiddata(json.uganda)


# Other
uga.aiddata <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/UgandaAMP_GeocodedResearchRelease_Level1_v1.3/data/level_1a.csv")
  

