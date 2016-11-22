# Example Using geoSIMEX 
# Using geoSIMEX with data from Uganda's third administrative division and dev. finance data from AidData
# AidData, REU
# Last Updated: 10/7/2016

# I'm using R version 3.3.1 when running the below code

# --- --- --- --- --- --- --- --- --- --- #
##### Loading Packages and Data #####
# --- --- --- --- --- --- --- --- --- --- #
library(raster)
library(sp)
library(rgdal)
library(doBy)

source("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/R/geoSIMEX.R")
uga.adm3.df <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/merge_uga_adm3.csv")
uga.aiddata <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/UgandaAMP_GeocodedResearchRelease_Level1_v1.3/data/level_1a.csv")

#json.uganda <- fromJSON("~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX/Example/summary.json")
#uga.aiddata <- subset.aiddata(json.uganda)

# --- --- --- --- --- --- --- --- --- --- #
##### Prepping Data for Analysis #####
# --- --- --- --- --- --- --- --- --- --- #

# Subsetting
#uga.aiddata <- uga.aiddata[uga.aiddata$ad_sector_names == "Government and civil society, general",]
#uga.aiddata <- uga.aiddata[(uga.aiddata$transactions_end_year >= 2000) &
#                                     (uga.aiddata$transactions_end_year <= 2010),]

##### Precision Code 6,8 projects don't have lat/long; give random lat/long
uga.aiddata$longitude[is.na(uga.aiddata$longitude)] <- uga.aiddata$longitude[!is.na(uga.aiddata$longitude)][1]
uga.aiddata$latitude[is.na(uga.aiddata$latitude)] <- uga.aiddata$latitude[!is.na(uga.aiddata$latitude)][1]

uga.aiddata <- uga.aiddata[!is.na(uga.aiddata$transactions_end_year),]

nrow(uga.aiddata)

##### Projecting AidData Data
coordinates(uga.aiddata) <- ~longitude+latitude
proj4string(uga.aiddata) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

##### Putting ADM heirarchy into AidData data
  # NOTE: Can typically use getData to get ADM info; however, issue with server where this comes 
  # from... so might be resolved later? Consequently, I'm importing from my local drive. And
  # exporting the resulting dataset to github.
if(TRUE){
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

# --- --- --- --- --- --- --- --- --- --- #
##### Checks on Data #####
# --- --- --- --- --- --- --- --- --- --- #

##### ADM Names have to be unique (e.g., subcounties can't have same name)
uga.aiddata$NAME_0 <- as.character(uga.aiddata$NAME_0)
uga.aiddata$NAME_1 <- as.character(uga.aiddata$NAME_1)
uga.aiddata$NAME_2 <- as.character(uga.aiddata$NAME_2)
uga.aiddata$NAME_3 <- as.character(uga.aiddata$NAME_3)

uga.adm3.df$NAME_0 <- as.character(uga.adm3.df$NAME_0)
uga.adm3.df$NAME_1 <- as.character(uga.adm3.df$NAME_1)
uga.adm3.df$NAME_2 <- as.character(uga.adm3.df$NAME_2)
uga.adm3.df$NAME_3 <- as.character(uga.adm3.df$NAME_3)

# Some ADM 3s don't have unique names
table(as.numeric(table(uga.adm3.df$NAME_3)))

# ADM 2s have unique names
table(as.numeric(table(summaryBy(.~NAME_2, data=uga.adm3.df)$NAME_2)))

# ADM 1s have unique names
table(as.numeric(table(summaryBy(.~NAME_1, data=uga.adm3.df)$NAME_1)))

# Making unique NAME_3 name
uga.aiddata$NAME_3_unique <- paste(uga.aiddata$NAME_3, uga.aiddata$NAME_2) 
uga.adm3.df$NAME_3_unique <- paste(uga.adm3.df$NAME_3, uga.adm3.df$NAME_2) 

# Confirming NAME_3_unique have
table(as.numeric(table(uga.adm3.df$NAME_3_unique)))

##### Making sure all ADM names in aiddata dataset are in ROI-level dataset
table(uga.aiddata$NAME_1 %in% uga.adm3.df$NAME_1)
uga.aiddata[uga.aiddata$NAME_1 %in% uga.adm3.df$NAME_1 == FALSE,]
uga.aiddata <- uga.aiddata[!is.na(uga.aiddata$NAME_1),]
table(uga.aiddata$NAME_1 %in% uga.adm3.df$NAME_1)

##### Converting ADM NAMES to IDs [change to work with strings!]
  # Basically, make what's happening here happen interally early on in geoSIMEX
uga.adm3.df$NAME_3_unique.id <- as.numeric(as.factor(uga.adm3.df$NAME_3_unique))
uga.adm3.df$NAME_2.id <- as.numeric(as.factor(uga.adm3.df$NAME_2))
uga.adm3.df$NAME_1.id <- as.numeric(as.factor(uga.adm3.df$NAME_1))
uga.adm3.df$NAME_0.id <- as.numeric(as.factor(uga.adm3.df$NAME_0))

uga.adm3.df_IDS <- subset(uga.adm3.df, select=c(NAME_0, NAME_1, NAME_2, NAME_3_unique,
                                                NAME_0.id, NAME_1.id, NAME_2.id, NAME_3_unique.id))

uga.aiddata <- merge(uga.aiddata, uga.adm3.df_IDS, by=c("NAME_0", "NAME_1", "NAME_2", "NAME_3_unique"))
uga.aiddata.reduced <- uga.aiddata

# --- --- --- --- --- --- --- --- --- --- #
##### Analysis at ADM2 Level #####
# --- --- --- --- --- --- --- --- --- --- #

# Analysis at ADM3 levels takes a while to run. So, here doing analysis at ADM2 level
# See next section for an example of ADM3 analysis, though.

uga.adm2.df <- summaryBy(. ~ NAME_2, data=uga.adm3.df, keep.names=TRUE)

##### Look at subset of AidData

# Projects ended between 2000 and 2010



# General budget support aid
#uga.aiddata.reduced <- uga.aiddata.reduced[uga.aiddata.reduced$ad_sector_names %in% c("General budget support"),]

##### Calculating Spatial Uncertainty of Dataset
(lambda <- calc_lambda(uga.aiddata.reduced, 
                       uga.adm3.df, 
                       roi.area="area", 
                       aid.precision.code="precision_code", 
                       roi.pc1.name="NAME_3_unique.id", 
                       roi.pc2.name="NAME_2.id", 
                       roi.pc3.name="NAME_1.id", 
                       roi.pc4.name="NAME_1.id", 
                       roi.pc5.name="NAME_1.id", 
                       roi.pc6.name="NAME_0.id", 
                       aid.pc1.centroid.name="NAME_1.id"))

##### Calculating Expected Aid in ROI
uga.aiddata.reduced$total_commitments <- log(uga.aiddata.reduced$total_commitments+1)
uga.adm2.df$expected_aid <- expected_aid_ROI(aidData=uga.aiddata.reduced, 
                                             roiData=uga.adm2.df, 
                                             probAidAssume=uga.adm2.df$area, 
                                             dollar_set=uga.aiddata.reduced$total_commitments, 
                                             aid.precision.code="precision_code", 
                                             roi.pc1.name="NAME_2.id", 
                                             roi.pc2.name="NAME_2.id", 
                                             roi.pc3.name="NAME_1.id", 
                                             roi.pc4.name="NAME_1.id", 
                                             roi.pc5.name="NAME_1.id", 
                                             roi.pc6.name="NAME_0.id", 
                                             aid.pc1.centroid.name="NAME_2.id")

uga.adm <- uga.adm2.df
uga.aiddata <- uga.aiddata.reduced

#uga.adm2.df$NTL_change <- uga.adm2.df$ncc4_2010e - uga.adm2.df$ncc4_2000e
naive_model <- lm(ncc4_2010e ~ expected_aid + gpw3_2000e, data=uga.adm)
summary(naive_model)


geoSIMEX_model <- geoSIMEX(model = naive_model, 
                           geoSIMEXvariable = "expected_aid", 
                           roiData = uga.adm, 
                           aidData = uga.aiddata, 
                           aid.amount = "total_commitments",
                           iterations = 100, 
                           bins = 4, 
                           roi.area = "area", 
                           roi.prob.aid = "area", 
                           roi.pc1.name="NAME_2.id", 
                           roi.pc2.name="NAME_2.id", 
                           roi.pc3.name="NAME_1.id", 
                           roi.pc4.name="NAME_1.id", 
                           roi.pc5.name="NAME_1.id", 
                           roi.pc6.name="NAME_0.id", 
                           aid.pc1.centroid.name="NAME_2.id", 
                           aid.precision.code="precision_code",
                           binary=FALSE,
                           sim_pc1=TRUE)

geoSIMEX_model$lambda
summary(naive_model)
summary(geoSIMEX_model)
plot(geoSIMEX_model, variable="expected_aid")

# --- --- --- --- --- --- --- --- --- --- --- --- ---#
##### Analysis at ADM2 Level - Interaction Term #####
# --- --- --- --- --- --- --- --- --- --- --- --- ---#

naive_model <- lm(ncc4_2010e ~ expected_aid + gpw3_2000e + expected_aid:gpw3_2000e, data=uga.adm2.df)

geoSIMEX_model <- geoSIMEX(model = naive_model, 
                           geoSIMEXvariable = "expected_aid", 
                           roiData = uga.adm2.df, 
                           aidData = uga.aiddata.reduced, 
                           aid.amount = "total_commitments",
                           iterations = 100, 
                           bins = 4, 
                           roi.area = "area", 
                           roi.prob.aid = "area", 
                           roi.pc1.name="NAME_2.id", 
                           roi.pc2.name="NAME_2.id", 
                           roi.pc3.name="NAME_1.id", 
                           roi.pc4.name="NAME_1.id", 
                           roi.pc5.name="NAME_1.id", 
                           roi.pc6.name="NAME_0.id", 
                           aid.pc1.centroid.name="NAME_0.id", 
                           aid.precision.code="precision_code",
                           binary=FALSE,
                           sim_pc1=TRUE)

geoSIMEX_model$lambda
summary(naive_model)
summary(geoSIMEX_model)
plot(geoSIMEX_model, variable="expected_aid")
plot(geoSIMEX_model, variable="expected_aid:gpw3_2000e")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- #
##### Analysis at ADM2 Level - Aid as Binary Term #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# To treat aid as a binary variable (i.e., whether an ROI recievied aid or note),
# we use the probability that an ROI recieved aid instead of using expected
# dollars of aid.

##### Calculating Expected Aid in ROI
uga.adm2.df$prob_aid <- prob_aid_ROI(aidData=uga.aiddata.reduced, 
                                             roiData=uga.adm2.df, 
                                             probAidAssume=uga.adm2.df$area, 
                                             dollar_set=uga.aiddata.reduced$total_commitments, 
                                             aid.precision.code="precision_code", 
                                             roi.pc1.name="NAME_2.id", 
                                             roi.pc2.name="NAME_2.id", 
                                             roi.pc3.name="NAME_1.id", 
                                             roi.pc4.name="NAME_1.id", 
                                             roi.pc5.name="NAME_1.id", 
                                             roi.pc6.name="NAME_0.id", 
                                             aid.pc1.centroid.name="NAME_1.id")

naive_model <- lm(ncc4_2010e ~ prob_aid + gpw3_2000e, data=uga.adm2.df)

geoSIMEX_model <- geoSIMEX(model = naive_model, 
                           geoSIMEXvariable = "prob_aid", 
                           roiData = uga.adm2.df, 
                           aidData = uga.aiddata.reduced, 
                           aid.amount = "total_commitments",
                           iterations = 100, 
                           bins = 4, 
                           roi.area = "area", 
                           roi.prob.aid = "area", 
                           roi.pc1.name="NAME_2.id", 
                           roi.pc2.name="NAME_2.id", 
                           roi.pc3.name="NAME_1.id", 
                           roi.pc4.name="NAME_1.id", 
                           roi.pc5.name="NAME_1.id", 
                           roi.pc6.name="NAME_0.id", 
                           aid.pc1.centroid.name="NAME_0.id", 
                           aid.precision.code="precision_code",
                           binary=TRUE,
                           sim_pc1=TRUE)

geoSIMEX_model$lambda
summary(naive_model)
summary(geoSIMEX_model)
plot(geoSIMEX_model, variable="prob_aid")

# --- --- --- --- --- --- --- --- --- --- #
##### Analysis at ADM3 Level #####
# --- --- --- --- --- --- --- --- --- --- #

##### Calculating Spatial Uncertainty of Dataset
(lambda <- calc_lambda(uga.aiddata.reduced, 
            uga.adm3.df, 
            roi.area="area", 
            aid.precision.code="precision_code", 
            roi.pc1.name="NAME_3_unique.id", 
            roi.pc2.name="NAME_2.id", 
            roi.pc3.name="NAME_1.id", 
            roi.pc4.name="NAME_1.id", 
            roi.pc5.name="NAME_1.id", 
            roi.pc6.name="NAME_0.id", 
            aid.pc1.centroid.name="NAME_1.id"))


##### Calculating Expected Aid in ROI
uga.adm3.df$expected_aid <- expected_aid_ROI(aidData=uga.aiddata.reduced, 
                                          roiData=uga.adm3.df, 
                                          probAidAssume=uga.adm3.df$area, 
                                          dollar_set=uga.aiddata.reduced$total_commitments, 
                                          aid.precision.code="precision_code", 
                                          roi.pc1.name="NAME_3_unique.id", 
                                          roi.pc2.name="NAME_2.id", 
                                          roi.pc3.name="NAME_1.id", 
                                          roi.pc4.name="NAME_1.id", 
                                          roi.pc5.name="NAME_1.id", 
                                          roi.pc6.name="NAME_0.id", 
                                          aid.pc1.centroid.name="NAME_1.id")

naive_model <- lm(ncc4_2010e ~ expected_aid, data=uga.adm3.df)

geoSIMEX_model <- geoSIMEX(model = naive_model, 
                           geoSIMEXvariable = "expected_aid", 
                           roiData = uga.adm3.df, 
                           aidData = uga.aiddata.reduced, 
                           aid.amount = "total_commitments",
                           iterations = 100, 
                           bins = 3, 
                           roi.area = "area", 
                           roi.prob.aid = "area", 
                           roi.pc1.name="NAME_3_unique.id", 
                           roi.pc2.name="NAME_2.id", 
                           roi.pc3.name="NAME_1.id", 
                           roi.pc4.name="NAME_1.id", 
                           roi.pc5.name="NAME_1.id", 
                           roi.pc6.name="NAME_0.id", 
                           aid.pc1.centroid.name="NAME_0.id", 
                           aid.precision.code="precision_code",
                           binary=FALSE,
                           sim_pc1=TRUE)

geoSIMEX_model$lambda
summary(naive_model)
summary(geoSIMEX_model)
plot(geoSIMEX_model, variable="expected_aid")



