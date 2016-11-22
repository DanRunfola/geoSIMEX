library(jsonlite)
library(raster)


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

aiddata.fill <- function(geo.data){
  
  ##### Precision Code 6,8 projects don't have lat/long; give random lat/long
  geo.data$longitude[is.na(geo.data$longitude)] <- geo.data$longitude[!is.na(geo.data$longitude)][1]
  geo.data$latitude[is.na(geo.data$latitude)] <- geo.data$latitude[!is.na(geo.data$latitude)][1]
  
  # Spatial Merging GADM data with AidData Data
  coordinates(geo.data) <- ~longitude+latitude
  proj4string(geo.data) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  gadm.data <- getData('GADM', country=as.character(geo.data$recipients_iso3[1]), level=2)
  
  aiddata.gadm.over <- over(geo.data, gadm.data)
  
  geo.data$NAME_0 <- aiddata.gadm.over$NAME_0
  geo.data$NAME_1 <- aiddata.gadm.over$NAME_1
  geo.data$NAME_2 <- aiddata.gadm.over$NAME_2
  geo.data$NAME_3 <- aiddata.gadm.over$NAME_3
  geo.data <- geo.data@data
  
  return(geo.data)
}

aid.data <- subset.aiddata(json.colombia)
aid.data <- aiddata.fill(aid.data)






