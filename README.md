
# geoSIMEX

> geoSIMEX

A package designed to incorporate spatial imprecision into regression analysis. Method based on the Simulation and Extrapolation Method (SIMEX), which is one approach to incorporate measurement imprecision into regression analysis. geoSIMEX adapts SIMEX to handle spatial imprecision.   

## Installation
```r
source("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/R/geoSIMEX.R")
```

## Usage

```r
##### Load Data #####
# ADM Level Data
uga.adm <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/merge_uga_adm3.csv")

# Aid Project Level Data
uga.aiddata <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/UgandaAMP_GeocodedResearchRelease_Level1_v1.3/data/level_1a.csv")


##### Calculating Expected Amount of Aid in Each ADM #####
uga.adm$expected_aid <- expected_aid_ROI(aidData=uga.aiddata, 
                                             roiData=uga.adm, 
                                             probAidAssume=uga.adm$area, 
                                             dollar_set=uga.aiddata$total_commitments, 
                                             aid.precision.code="precision_code", 
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
                           iterations = 25, 
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

```

## License

Copyright 2016 Robert Marty (ramarty@email.wm.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

