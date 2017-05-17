
# geoSIMEX

> geoSIMEX

A package designed to incorporate spatial imprecision into regression analysis. Method based on the Simulation and Extrapolation Method (SIMEX), which is one approach to incorporate measurement imprecision into regression analysis. geoSIMEX adapts SIMEX to handle spatial imprecision.   

## Installation
```r
library(devtools)
install_github("itpir/geoSIMEX")
```

## Usage

```r
##### Generating Country-Level Dataset #####
numSubcounties <- 120
numSubcountyInCounty <- 2 
numCountyInDistrict <- 3
numDistrictInRegion <- 2

N <- numSubcounties
subcounty <- 1:N
county <- rep(1:(N/numSubcountyInCounty), each=numSubcountyInCounty)
district <- rep(1:(N/(numSubcountyInCounty*numCountyInDistrict)), each=(numSubcountyInCounty*numCountyInDistrict))
region <- rep(1:(N/(numSubcountyInCounty*numCountyInDistrict*numDistrictInRegion)), each=(numSubcountyInCounty*numCountyInDistrict*numDistrictInRegion))
country <- 1

subcountyArea <- runif(N)
probAid_assumed <- runif(N)

subcountyData <- as.data.frame(cbind(subcounty,county,district,region,country,probAid_assumed,subcountyArea))

##### Creating Aid Dataset #####
numberProjects = 50
aidData <- as.data.frame(matrix(NA,nrow=numberProjects,ncol=3))
names(aidData) <- c("aid","trueSubcounty","PC")
aidData$aid <- runif(nrow(aidData)) * 100
probAid_true <- runif(N)
aidData$trueSubcounty <- sample(size=numberProjects,x=c(1:N), prob=probAid_true, replace=TRUE)
aidData$PC  <- sample(size=numberProjects, x=c(1,2,3,4,6), prob=runif(5), replace=TRUE)

# True Aid
aidData$PC.1s <- 1
subcountyData$trueAid <- expected_aid_ROI(aidData=aidData, 
                                         roiData=subcountyData, 
                                         roi.prob.aid="probAid_assumed", 
                                         aid.project.amount="aid", 
                                         aid.precision.code="PC.1s", 
                                         roi.pc1.name="subcounty", 
                                         roi.pc2.name="county", 
                                         roi.pc3.name="district", 
                                         roi.pc4.name="region", 
                                         roi.pc5.name="region", 
                                         roi.pc6.name="country", 
                                         aid.pc1.centroid.name="trueSubcounty")

# Wealth - 1 to 1 relation with aid
subcountyData$wealth <- subcountyData$trueAid + runif(nrow(subcountyData))

# Expected Value Aid 
subcountyData$expectedAid <- expected_aid_ROI(aidData=aidData, 
                                             aid.project.amount="aid", 
                                             aid.precision.code="PC", 
                                             aid.pc1.centroid.name="trueSubcounty",
                                             roiData=subcountyData, 
                                             roi.prob.aid="probAid_assumed", 
                                             roi.pc1.name="subcounty", 
                                             roi.pc2.name="county", 
                                             roi.pc3.name="district", 
                                             roi.pc4.name="region", 
                                             roi.pc5.name="region", 
                                             roi.pc6.name="country")

naive_model <- lm(wealth ~ expectedAid, data=subcountyData)

geoSIMEX_model <- geoSIMEX(model = naive_model, 
                          geoSIMEXvariable = "expectedAid", 
                          aidData = aidData, 
                          aid.project.amount = "aid",
                          aid.pc1.centroid.name="trueSubcounty", 
                          aid.precision.code="PC",
                          roiData = subcountyData, 
                          roi.area = "subcountyArea", 
                          roi.prob.aid = "probAid_assumed", 
                          roi.pc1.name="subcounty", 
                          roi.pc2.name="county", 
                          roi.pc3.name="district", 
                          roi.pc4.name="region", 
                          roi.pc5.name="region", 
                          roi.pc6.name="country")

##### Analyzing Results #####
summary(naive_model)
summary(geoSIMEX_model)
plot(geoSIMEX_model, variable="expectedAid")
```

## License

Copyright 2017 Robert Marty (ramarty@email.wm.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

