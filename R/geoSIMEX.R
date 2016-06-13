##### SECTION 1: GEOSIMEX SIMULATION FUNCTION #####
setwd("~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX")
roxygen2::roxygenise()

#' @title Geographic SIMEX
#' 
#' @author AidData
#' 
#' @import parallel
#' 
#' @description
#' \code{geoSIMEX} Implementation of the geoSIMEX algorithm for models 
#' with spatial uncertainty. Package built to work with data from 
#' AidData's data extraction tool.
#' 
#' @usage
#' geoSIMEX(model, geoSIMEXvariable, roiData, aidData, aid.amount, 
#' iterations=500, bins=3, fitting.method = "quadratic", roi.area="area",  
#' roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", 
#' roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc5.name="pc6.id",  
#' aid.pc1.centroid.name="centroid.pc1.id", aid.precision.code="precision.code",
#' parallel=TRUE, mc.cores=2)
#' 
#' @param model the naive model
#' @param SIMEXvariable character containing the name of the variable with spatial uncertainty
#' @param roiData name of dataframe of ROI data 
#' @param aidData name of dataframe of aid project data
#' @param aid.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
#' @param iterations number of simulated error iterations
#' @param bins number of bins to group coefficients
#' @param fitting.method fitting method for the extrapolation. linear and quadratic are implemented.
#' @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param parallel parallelize iterations
#' @param mc.cores number of cores to use for parallelization
#' 
#' @details The values contained within roi.pc1.name and aid.pc1.centroid.name variables should be the same.
#' 
#' @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website].
#' 
#' @references Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. Journal of American Statistical Association, 89, 1314 – 1328
#' 
#' @examples 
#' set.seed(42)
#' 
#' # Generating Country Dataset
#' countryData <- as.data.frame(matrix(NA, nrow=120, ncol=0))
#' countryData$pc1.id <- 1:120
#' countryData$pc2.id <- rep(1:(120/3), each=3)
#' countryData$pc3.id <- rep(1:(120/(3*2)), each=(3*2))
#' countryData$pc4.id <- rep(1:(120/(3*2*4)), each=(3*2*4))
#' countryData$pc6.id <- 1
#' countryData$area <- rgamma(120, shape=2)
#' 
#' # Creating Aid Dataset Without Error
#' aidData <- as.data.frame(matrix(NA,nrow=100,ncol=0))
#' aidData$disbursement <- runif(100,0,1) * 10000000 
#' aidData$centroid.pc1.id <- sample(size=100,x=c(1:120), prob=rep(1/120, 120), replace=TRUE)
#' aidData$precision.code <- 1
#' 
#' # Adding True Aid to Country Dataset
#' countryData$aid <- expected_aid_ROI(aidData=aidData, roiData=countryData, probAidAssume=countryData$area, dollar_set=aidData$disbursement, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
#' 
#' # Defining True Relation Between Aid and Wealth
#' countryData$wealth <- countryData$aid + rnorm(120) * 0.1
#' 
#' # Creating Datasets with Uncertainty
#' countryData <- subset(countryData, select = -c(aid))
#' aidData$precision.code <- sample(size=100, x=c(1,2,3,4,6), prob=rep(1/5, 5), replace=TRUE)
#' 
#' # Calculating Expected Aid and Running Naive Model
#' countryData$Expected.Aid <- expected_aid_ROI(aidData=aidData, roiData=countryData, probAidAssume=countryData$area, dollar_set=aidData$disbursement, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
#' 
#' lm_naive <- lm(wealth ~ Expected.Aid, data=countryData)
#' 
#' # Implementing GeoSIMEX
#' lm_geoSIMEX <- geoSIMEX(model = lm_naive, 
#'                         geoSIMEXvariable = "Expected.Aid", 
#'                         roiData = countryData, 
#'                         aidData = aidData, 
#'                         aid.amount = "disbursement")
#' 
#' summarize(lm_geoSIMEX)
#' geoSIMEX_plot(lm_geoSIMEX, variable="Expected.Aid")
geoSIMEX <- function(model, 
                     geoSIMEXvariable, 
                     roiData, 
                     aidData, 
                     aid.amount,
                     iterations=500, 
                     bins=4, 
                     fitting.method = "quadratic", 
                     roi.area="area",  
                     roi.pc1.name="pc1.id", 
                     roi.pc2.name="pc2.id", 
                     roi.pc3.name="pc3.id", 
                     roi.pc4.name="pc4.id", 
                     roi.pc5.name="pc5.id", 
                     roi.pc6.name="pc6.id",  
                     aid.pc1.centroid.name="centroid.pc1.id", 
                     aid.precision.code="precision.code",
                     parallel=TRUE, 
                     mc.cores=2){
    
  ##### Creating Variables
  # Creating probability variable from area
  probAid_area <- as.matrix(roiData[roi.area] / sum(roiData[roi.area]))
  
  # Creating matrix of original precision codes
  precision.code.original <- as.matrix(aidData[aid.precision.code])
  
  ##### Calculating Maximum Lambda Denominator #####
  aidData[aid.precision.code] <- 6  
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAid_area, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)  
  param_set[param_set!=0] <- 1
  param_set <- param_set * as.matrix(roiData[roi.area])
  maxLambda_denom <- sum(colSums(param_set))
  
  ##### Calculate Naive Lambda #####
  aidData[aid.precision.code] <- precision.code.original 
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAid_area, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)  
  lambda_naive <- calcLambda(param_set, maxLambda_denom, as.matrix(roiData[roi.area]))
 
  ##### Simulating Data with Additional Error and Putting in One Dataframe #####
  prob.increase.list <- runif(iterations,0,1)
  
  geoSimulateError.results <- lapply(prob.increase.list, geoSimulateError, 
                                     aidData=aidData, 
                                     roiData=roiData, 
                                     probAidAssume=probAid_area, 
                                     PC_researcherSees=precision.code.original, 
                                     maxLambda_denom=maxLambda_denom,
                                     roi.area=roi.area,
                                     aid.precision.code=aid.precision.code, 
                                     roi.pc1.name=roi.pc1.name, 
                                     roi.pc2.name=roi.pc2.name, 
                                     roi.pc3.name=roi.pc3.name, 
                                     roi.pc4.name=roi.pc4.name, 
                                     roi.pc5.name=roi.pc5.name, 
                                     roi.pc6.name=roi.pc6.name, 
                                     aid.pc1.centroid.name=aid.pc1.centroid.name,
                                     aid.amount=aid.amount,
                                     model=model,
                                     geoSIMEXvariable=geoSIMEXvariable)
  
  geoSimulateError.results.df <- matrix(NA, ncol=ncol(geoSimulateError.results[[1]]),nrow=0)
  geoSimulateError.results.df <- as.data.frame(geoSimulateError.results.df)
  for(i in 1:length(geoSimulateError.results)){
    geoSimulateError.results.df <- rbind(geoSimulateError.results.df, geoSimulateError.results[[i]])
  }
  
    
  ##### Extrapolated Mean Coefficient #####  
  
  # Putting Into Bins
  extrapolatedMean.df <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(geoSimulateError.results.df)))
  
  minLambda <- min(geoSimulateError.results.df$lambda)
  maxLambda <- max(geoSimulateError.results.df$lambda)
  
  binSize <- (maxLambda - minLambda) / bins
  
  binSize_lb <- minLambda
  binSize_ub <- minLambda + binSize
  
  numIter <- (maxLambda - minLambda) / binSize
  for(i in 1:ceiling(numIter)){
    
    geoSimulateError.temp <- geoSimulateError.results.df[(geoSimulateError.results.df$lambda >= binSize_lb) & (geoSimulateError.results.df$lambda <= binSize_ub),]
    geoSimulateError.temp.mean <- colMeans(geoSimulateError.temp)
    
    extrapolatedMean.df <- rbind(extrapolatedMean.df, geoSimulateError.temp.mean)
    
    binSize_lb <- binSize_lb + binSize
    binSize_ub <- binSize_ub + binSize 
  }
  names(extrapolatedMean.df) <- names(geoSimulateError.results.df)
  
  # Regressions Using Bins
  extrapolatedMean.df$lambda_sq <- extrapolatedMean.df$lambda^2
  
  numVars <- ncol(extrapolatedMean.df) - 2
  coef.geoSIMEX <- matrix(NA, nrow=1, ncol=numVars)
  for(i in 1:numVars){
    coef.geoSIMEX[i] <- summary(lm(as.matrix(extrapolatedMean.df[i]) ~ lambda + lambda_sq, data = extrapolatedMean.df))$coefficients[1]
  }
  coef.geoSIMEX <- as.data.frame(coef.geoSIMEX)
  names(coef.geoSIMEX) <- head(names(extrapolatedMean.df), -2)
  
  ##### Bootstrap Standard Error #####
  numberBootIter <- nrow(geoSimulateError.results.df)
  
  geoSimulateError.results.df$lambda_sq <- geoSimulateError.results.df$lambda^2
  bootIter.list <- lapply(seq(1:numberBootIter), bootIter, 
                          geoSimulateError.results.df=geoSimulateError.results.df, 
                          bins=bins, 
                          numFromBin=1)
  
  boot.coefs.matrix <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(bootIter.list[[1]])))
  for(i in 1:length(bootIter.list)){
    boot.coefs.matrix <- rbind(boot.coefs.matrix, bootIter.list[[i]])
  }
    
  se.geoSIMEX <- apply(boot.coefs.matrix, 2, sd)
  
  ##### Collecting Results #####
  row.names(coef.geoSIMEX) <- "Coefficients"
  coef <- t(coef.geoSIMEX)
  se <- as.data.frame(se.geoSIMEX)
  names(se) <- "Std. Errors"
  coef.se <- cbind(coef,se)
  
  return(list(coefficients=coef.se,
              values=geoSimulateError.results.df,
              valuesMean=extrapolatedMean.df))
}

##### SECTION 2: FUNCTIONS TO BE INCLUDED IN PACKAGE #####
#' Calculates expected value of aid for each ROI
#' 
#' @title Expected Value Aid
#' 
#' @author AidData
#' 
#' @description
#' \code{expected_aid_ROI} Calculates expected value of aid.
#' 
#' @usage
#' expected_aid_ROI(aidData, roiData, probAidAssume, dollar_set, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
#'  
#' @param roiData name of dataframe of ROI data 
#' @param aidData name of dataframe of aid project data
#' @param aid.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
#' @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#'  
#' @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website].
expected_aid_ROI <- function(aidData, roiData, probAidAssume, dollar_set, aid.precision.code, roi.pc1.name, roi.pc2.name, roi.pc3.name, roi.pc4.name, roi.pc5.name, roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name){
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  aid <- dollar_expected_value(param_set=param_set, dollar_set=dollar_set)
  return(aid)
}



#' @title Spatial Uncertainty (Lambda)
#' 
#' @author AidData
#' 
#' @description
#' \code{calc_lambda} Calculates spatial uncertainty (lambda) of aid project dataset
#' 
#' @usage
#' calc_lambda(aidData, roiData, roi.area="area", aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
#'  
#' @param roiData name of dataframe of ROI data 
#' @param aidData name of dataframe of aid project data
#' @param aid.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
#' @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#'  
#' @example 
#' lambda.naive <- calc_lambda(aidData=aidData.sub@data, roiData=roiData)
#' 
#' @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website].

calc_lambda <- function(aidData, 
                        roiData, 
                        roi.area="area", 
                        aid.precision.code="precision.code", 
                        roi.pc1.name="pc1.id", 
                        roi.pc2.name="pc2.id", 
                        roi.pc3.name="pc3.id", 
                        roi.pc4.name="pc4.id", 
                        roi.pc5.name="pc5.id", 
                        roi.pc6.name="pc6.id", 
                        aid.pc1.centroid.name="centroid.pc1.id"){
  
  # Prepping Variables
  pc.original <- aidData[,aid.precision.code] 
  area <- roiData[,roi.area]
  probAidAssume <- area / sum(area)
  
  # Calculaitng Max Lambda
  aidData[,aid.precision.code] <- 6
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set[param_set!=0] <- 1
  param_set <- param_set*area
  maxLambda_denom <- sum(colSums(param_set))
  
  aidData[,aid.precision.code] <- pc.original
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set[param_set!=0] <- 1
  param_set <- param_set*area
  lambda <- sum(colSums(param_set)) / maxLambda_denom
  return(lambda)
}


##### SECTION 3: DEFINING FUNCTIONS THAT OTHER FUNCTIONS USE #####

calcLambda <- function(param_set, maxLambda_denom, area){
  param_set[param_set!=0] <- 1
  param_set <- param_set*area
  lambda <- sum(colSums(param_set)) / maxLambda_denom
  return(lambda)
}

# Creating Parameter Matrix
paramCol <- function(i, aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name){
  
  # Getting information about aid project 
  aidDataPrj_Temp <- aidData[i,]
  PC_temp <- as.numeric(aidDataPrj_Temp[aid.precision.code])
  
  # Getting information about subcounty & region that the aid project falls in
  subcounty_temp <- roiData[roiData[roi.pc1.name] == as.numeric(aidDataPrj_Temp[aid.pc1.centroid.name]),]
  
  # Making temporary dataset that will merge into
  paramCol <- roiData
  paramCol$prj <- 0
  
  if(PC_temp == 1){
    paramCol$prj[paramCol[roi.pc1.name] == as.numeric(subcounty_temp[roi.pc1.name])] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  if(PC_temp == 2){
    paramCol$prj[paramCol[roi.pc2.name] == as.numeric(subcounty_temp[roi.pc2.name])] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  if(PC_temp == 3){
    paramCol$prj[paramCol[roi.pc3.name] == as.numeric(subcounty_temp[roi.pc3.name])] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  if(PC_temp == 4){
    paramCol$prj[paramCol[roi.pc4.name] == as.numeric(subcounty_temp[roi.pc4.name])] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  if((PC_temp == 6) | (PC_temp == 8)){
    paramCol$prj[paramCol[roi.pc6.name] == as.numeric(subcounty_temp[roi.pc6.name])] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  row.names(paramCol) <- as.matrix(paramCol[roi.pc1.name])
  
  paramCol_return <- as.data.frame(paramCol$prj)
  row.names(paramCol_return) <- row.names(paramCol)
  names(paramCol_return) <- paste("prj",i,sep="")
  
  return(paramCol_return)
}

paramSet <- function(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name){
  param_set = as.data.frame(lapply(1:nrow(aidData), paramCol, aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name))
  return(param_set) 
}

dollar_expected_value <- function(param_set, dollar_set){    
  paramDollars <- lapply(1:length(dollar_set), function(i) param_set[,i] * dollar_set[i])
  paramDollars <- as.data.frame(paramDollars)
  dollars_expected <- rowSums(paramDollars)
  return(dollars_expected)
}

geoSimulateError <- function(probIncPC, aidData=aidData, roiData=roiData, probAidAssume=probAid_area, PC_researcherSees=precision.code.original, maxLambda_denom=maxLambda_denom, roi.area=roi.area, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name, aid.amount=aid.amount, model=model, geoSIMEXvariable=geoSIMEXvariable){
    
  # Initializing Results Matrix
  results <- matrix(NA, ncol=2,nrow=1)
  results <- as.data.frame(results)
  names(results) <- c("lambda","aid.expected.coef")
  
  # Simulating Error
  aidData[aid.precision.code] <- PC_researcherSees
  probStayPC <- 1 - probIncPC
    
  # ... Precision Code 1s
  aidData[aid.precision.code][aidData[aid.precision.code] == 1] <- sample(size=length(aidData[aid.precision.code][aidData[aid.precision.code] == 1]), 
                                        x=c(1,2,3,4,6), 
                                        prob=c(probStayPC,probIncPC/4,probIncPC/4,probIncPC/4,probIncPC/4), 
                                        replace=TRUE)
  
  # ... Precision Code 2s
  aidData[aid.precision.code][aidData[aid.precision.code] == 2] <- sample(size=length(aidData[aid.precision.code][aidData[aid.precision.code] == 2]), 
                                        x=c(2,3,4,6), 
                                        prob=c(probStayPC,probIncPC/3,probIncPC/3,probIncPC/3), 
                                        replace=TRUE)
  
  # ... Precision Code 3s
  aidData[aid.precision.code][aidData[aid.precision.code] == 3] <- sample(size=length(aidData[aid.precision.code][aidData[aid.precision.code] == 3]), 
                                        x=c(3,4,6), 
                                        prob=c(probStayPC,probIncPC/2,probIncPC/2), 
                                        replace=TRUE)
  
  # ... Precision Code 4s
  aidData[aid.precision.code][aidData[aid.precision.code] == 4] <- sample(size=length(aidData[aid.precision.code][aidData[aid.precision.code] == 4]), 
                                        x=c(4,6), 
                                        prob=c(probStayPC,probIncPC), 
                                        replace=TRUE)

  # Calculating paramSet
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  
  # Update Model With New Expected Aid
  df.temp <- model$model
  df.temp[geoSIMEXvariable] <- dollar_expected_value(param_set=param_set, dollar_set=as.matrix(aidData[aid.amount]))
  model.SIMEX <- update(model, data = df.temp)
  
  # Calculate Lambda
  lambda <- calcLambda(param_set, maxLambda_denom, as.matrix(roiData[roi.area]))
  
  # Collecting Results
  model.SIMEX.coefs <- as.data.frame(as.list(model.SIMEX$coefficients))
  names(model.SIMEX.coefs) <- names(model.SIMEX$coefficients)
  model.SIMEX.coefs$lambda <- lambda

  return(model.SIMEX.coefs)
}

bootIter <- function(i, geoSimulateError.results.df, bins, numFromBin){
  
  geoSimulateError.results.df$rand <- runif(nrow(geoSimulateError.results.df))
  geoSimulateError.results.df <- geoSimulateError.results.df[order(geoSimulateError.results.df$rand),] 
  
  results.boot <- matrix(NA,nrow=0,ncol=ncol(geoSimulateError.results.df))
  results.boot <- as.data.frame(results.boot)
  
  boot.coefs <- as.data.frame(matrix(NA, nrow=1, ncol=ncol(geoSimulateError.results.df)-3))
  names(boot.coefs) <- head(names(geoSimulateError.results.df), -3)
  
  # Randomly Pulling from Results Matrix to Get Subsample for Bootstrapping 
  minLambda <- min(geoSimulateError.results.df$lambda)
  maxLambda <- max(geoSimulateError.results.df$lambda)
  
  binSize <- (maxLambda - minLambda) / bins
    
  binSize_lb <- minLambda
  binSize_ub <- minLambda + binSize
  
  numIter <- (maxLambda - minLambda) / binSize
  for(i in 1:ceiling(numIter)){
    
    results.boot <- rbind(results.boot, geoSimulateError.results.df[(geoSimulateError.results.df$lambda >= binSize_lb) & (geoSimulateError.results.df$lambda <= binSize_ub),][1:numFromBin,])
    
    binSize_lb <- binSize_lb + binSize
    binSize_ub <- binSize_ub + binSize 
  }
  
  # Extrapolating 
  numVars <- ncol(results.boot) - 3
  coefs.geoSIMEX.boot <- matrix(NA, nrow=1, ncol=numVars)
  for(i in 1:numVars){
    coefs.geoSIMEX.boot[i] <- summary(lm(as.matrix(results.boot[i]) ~ lambda + lambda_sq, data = results.boot))$coefficients[1]
  }
  
  coefs.geoSIMEX.boot <- as.data.frame(coefs.geoSIMEX.boot)
  names(coefs.geoSIMEX.boot) <- head(names(results.boot), -3)
  
  return(coefs.geoSIMEX.boot)
}

##### SECTION 4: DEFINING FUNCTIONS THAT ACT ON geoSIMEX PACKAGE #####
summarize <- function(glm_geoSIMEX){
  return(glm_geoSIMEX$coefficients)
}

geoSIMEX_plot <- function(glm_geoSIMEX, variable){
  return(plot(glm_geoSIMEX$valuesMean$lambda, glm_geoSIMEX$valuesMean[,variable],
       xlab="lambda",
       ylab=variable,
       main=paste("geoSIMEX Plot of ",variable, sep=""),
       pch=20))
}



