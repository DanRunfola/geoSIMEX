# Example Using geoSIMEX 
# Using geoSIMEX with data from Uganda's third administrative division and dev. finance data from AidData
# AidData, REU
# Last Updated: 10/7/2016

# --- --- --- --- --- --- --- --- --- --- #
##### Loading Packages and Data #####
# --- --- --- --- --- --- --- --- --- --- #
library(raster)
library(sp)
library(rgdal)

source("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/R/geoSIMEX.R")
uga.adm3.df <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/merge_uga_adm3.csv")
uga.aiddata <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/UgandaAMP_GeocodedResearchRelease_Level1_v1.3/data/level_1a.csv")

# --- --- --- --- --- --- --- --- --- --- #
##### Prepping Data for Analysis #####
# --- --- --- --- --- --- --- --- --- --- #

##### Precision Code 6,8 projects don't have lat/long; give random lat/long
uga.aiddata$longitude[is.na(uga.aiddata$longitude)] <- uga.aiddata$longitude[!is.na(uga.aiddata$longitude)][1]
uga.aiddata$latitude[is.na(uga.aiddata$latitude)] <- uga.aiddata$latitude[!is.na(uga.aiddata$latitude)][1]

##### Projecting AidData Data
coordinates(uga.aiddata) <- ~longitude+latitude
proj4string(uga.aiddata) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

##### Putting ADM heirarchy into AidData data
  # NOTE: Can typically use getData to get ADM info; however, issue with server where this comes 
  # from... so might be resolved later? Consequently, I'm importing from my local drive. And
  # exporting the resulting dataset to github.
if(FALSE){
  setwd("~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX/Example/UGA_adm_shp")
  uga.3 <- readOGR(dsn = ".", layer = "UGA_adm3")
  
  aiddata.gadm.over <- over(uga.aiddata, uga.3)
  
  uga.aiddata$NAME_0 <- aiddata.gadm.over$NAME_0
  uga.aiddata$NAME_1 <- aiddata.gadm.over$NAME_1
  uga.aiddata$NAME_2 <- aiddata.gadm.over$NAME_2
  uga.aiddata$NAME_3 <- aiddata.gadm.over$NAME_3
  uga.aiddata <- uga.aiddata@data
  
  write.csv(uga.aiddata, "~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX/Example/uga_aiddata_gadm.csv")
} else{
  uga.aiddata <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/uga_aiddata_gadm.csv")
}


uga.3$area <- area(uga.3)

uga.adm3.df_a <- merge(uga.adm3.df, uga.3, by=c("NAME_0", "NAME_1", "NAME_2", "NAME_3"))

write.csv(uga.adm3.df, "~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX/Example/merge_uga_adm3.csv")

# --- --- --- --- --- --- --- --- --- --- #
##### Checks on Data #####
# --- --- --- --- --- --- --- --- --- --- #

# Making sure all ADM names in aiddata dataset are in ROI-level dataset
table(uga.aiddata$NAME_1 %in% uga.adm3.df$NAME_1)
uga.aiddata[uga.aiddata$NAME_1 %in% uga.adm3.df$NAME_1 == FALSE,]
uga.aiddata <- uga.aiddata[!is.na(uga.aiddata$NAME_1),]
table(uga.aiddata$NAME_1 %in% uga.adm3.df$NAME_1)

# --- --- --- --- --- --- --- --- --- --- #
##### Analysis #####
# --- --- --- --- --- --- --- --- --- --- #

##### Calculating Expected Aid in ROI
subcountyData$trueAid <- expected_aid_ROI(aidData=uga.aiddata, 
                                          roiData=uga.adm3.df, 
                                          probAidAssume=probAid_equal, 
                                          dollar_set=aidData$aid, 
                                          aid.precision.code="PC", 
                                          roi.pc1.name="subcounty", 
                                          roi.pc2.name="county", 
                                          roi.pc3.name="district", 
                                          roi.pc4.name="region", 
                                          roi.pc5.name="region", 
                                          roi.pc6.name="country", 
                                          aid.pc1.centroid.name="trueSubcounty")

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


