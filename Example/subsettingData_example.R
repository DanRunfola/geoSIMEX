library(jsonlite)

# When load geoSIMEX package, the following dataset will be loaded:
colombiaaims_geocodedresearchrelease_level1_v1_1_1 <- read.csv("~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX/Example/ColombiaAIMS_GeocodedResearchRelease_Level1_v1.1.1/data/level_1a.csv")

# Load JSON
json.colombia <- fromJSON("~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX/Example/summary.json")

subset.aiddata <- function(json){
  
  # Import dataset. Dataset names must be loaded from Rpackage and be the same that appears in JSON.
  geo.data <- eval(parse(text=json$dataset))
  
  # Number of filters to subset by
  num.filters <- length(json$request$options$filters)
  
  for(i in 1:num.filters){
    
    # Checking to make sure filter exists in dataset (if doesn't, skip). Filter name in JSON must match filter
    if(row.names(as.matrix(json$request$options$filters[i])) %in% names(geo.data)){ 
      geo.data <- geo.data[geo.data[row.names(as.matrix(json$request$options$filters[i]))] == as.character(json$request$options$filters[i]),]
    }
  }
  return(geo.data)
}

aid.data <- subset.aiddata(json.colombia)



