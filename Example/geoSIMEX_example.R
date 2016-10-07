# Example Using geoSIMEX 
# Using geoSIMEX with data from Uganda's third administrative division and dev. finance data from AidData
# AidData, REU
# Last Updated: 10/7/2016

# --- --- --- --- --- --- --- --- --- --- #
##### Loading Packages and Data #####
# --- --- --- --- --- --- --- --- --- --- #
library(raster)
library(sp)

source("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/R/geoSIMEX.R")
uga.adm3.df <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/merge_uga_adm3.csv")
uga.aiddata <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/UgandaAMP_GeocodedResearchRelease_Level1_v1.3/data/level_1a.csv")

# --- --- --- --- --- --- --- --- --- --- #
##### Prepping Data for Analysis #####
# --- --- --- --- --- --- --- --- --- --- #

uga.aiddata$longitude[is.na(uga.aiddata$longitude)] <- uga.aiddata$longitude[!is.na(uga.aiddata$longitude)][1]
uga.aiddata$latitude[is.na(uga.aiddata$latitude)] <- uga.aiddata$latitude[!is.na(uga.aiddata$latitude)][1]

coordinates(uga.aiddata) <- ~longitude+latitude
proj4string(uga.aiddata) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

over(uga.aiddata, uga.1)

# Get Uganda GADM Data to Match with AidData Data
uga.0 <- getData('GADM', country='UGA', level=0)
uga.1 <- getData('GADM', country='UGA', level=1)
uga.2 <- getData('GADM', country='UGA', level=2)
uga.3 <- getData('GADM', country='UGA', level=3)









# 1. get ADM data on AidData data to match uga.adm3.df data
# 2. do a check to make sure sure all names in aiddata data are in uga.adm3.df data



# --- --- --- --- --- --- --- --- --- --- #
##### Analysis #####
# --- --- --- --- --- --- --- --- --- --- #

geoSIMEX_model <- geoSIMEX(model = naive_model, 
                           geoSIMEXvariable = "expectedAid", 
                           roiData = subcountyData, 
                           aidData = aidData, 
                           aid.amount = "aid",
                           iterations = simex_randProb_iterations, 
                           bins = bins, 
                           roi.area = "subcountyArea", 
                           roi.prob.aid = "prob.aid.roi.use", 
                           roi.pc1.name="subcounty", 
                           roi.pc2.name="county", 
                           roi.pc3.name="district", 
                           roi.pc4.name="region", 
                           roi.pc5.name="region", 
                           roi.pc6.name="country",  
                           aid.pc1.centroid.name="trueSubcounty", 
                           aid.precision.code="PC",
                           binary=FALSE,
                           sim_pc1=sim_pc1)


