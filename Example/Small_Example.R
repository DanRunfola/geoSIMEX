# Installation
#source("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/R/geoSIMEX.R")
source("~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX/R/geoSIMEX_updating.R")

##### Load Data #####
# ADM Level Data
uga.adm <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/merge_uga_adm3.csv")

# Aid Project Level Data
uga.aiddata <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/UgandaAMP_GeocodedResearchRelease_Level1_v1.3/data/level_1a.csv")

##### Calculating Expected Amount of Aid in Each ADM #####
uga.adm$expected_aid <- expected_aid_ROI(aidData=uga.aiddata, 
                                         roiData=uga.adm, 
                                         roi.prob.aid=uga.adm$area, 
                                         dollar_set=uga.aiddata$total_commitments, 
                                         aid.precision.code=uga.aiddata$precision_code, 
                                         roi.pc1.name="NAME_2.id", 
                                         roi.pc2.name="NAME_2.id", 
                                         roi.pc3.name="NAME_1.id", 
                                         roi.pc4.name="NAME_1.id", 
                                         roi.pc5.name="NAME_1.id", 
                                         roi.pc6.name="NAME_0.id", 
                                         aid.pc1.centroid.name="NAME_2.id")

##### Run Naive Model #####
naive_model <- lm(ncc4_2010e ~ expected_aid + gpw3_2000e, data=uga.adm)

# View Results
summary(naive_model)

##### Run geoSIMEX Model #####
geoSIMEX_model <- geoSIMEX(model = naive_model, 
                           geoSIMEXvariable = "expected_aid", 
                           roiData = uga.adm, 
                           aidData = uga.aiddata, 
                           aid.amount = "total_commitments",
                           iterations = 500, 
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

# View Results
summary(geoSIMEX_model)
plot(geoSIMEX_model, variable = "expected_aid")