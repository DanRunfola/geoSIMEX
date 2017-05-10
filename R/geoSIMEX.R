##### SECTION 1: GEOSIMEX SIMULATION FUNCTION #####
# setwd("~/Desktop/AidData/MeasureErrorsInEx/geoSIMEX/geoSIMEX")
#roxygen2::roxygenise()
#devtools::document()

# https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
# document()

geoSIMEX_est <- function(model, 
                         geoSIMEXvariable, 
                         roiData, 
                         aidData, 
                         aid.project.amount,
                         iterations, 
                         bins,
                         number.from.bin,
                         number.from.bin.average,
                         fitting.method, 
                         roi.area,  
                         roi.prob.aid,
                         roi.pc1.name, 
                         roi.pc2.name, 
                         roi.pc3.name, 
                         roi.pc4.name, 
                         roi.pc5.name, 
                         roi.pc6.name,  
                         aid.pc1.centroid.name, 
                         aid.precision.code,
                         binary,
                         sim_pc1,
                         parallel,
                         extrapolation,
                         mc.cores,
                         status.bar){
  
  ##### Converting roi.names into Characters #####
  
  roiData[,roi.pc1.name] <- as.character(roiData[,roi.pc1.name])
  roiData[,roi.pc2.name] <- as.character(roiData[,roi.pc2.name])
  roiData[,roi.pc3.name] <- as.character(roiData[,roi.pc3.name])
  roiData[,roi.pc4.name] <- as.character(roiData[,roi.pc4.name])
  roiData[,roi.pc5.name] <- as.character(roiData[,roi.pc5.name])
  roiData[,roi.pc6.name] <- as.character(roiData[,roi.pc6.name])
  
  aidData[,aid.pc1.centroid.name] <- as.character(aidData[,aid.pc1.centroid.name])
  
  ##### Subsetting ROI data
  for(var in names(model$model)){
    roiData <- roiData[!is.na(roiData[var]),]
  }
  
  ##### Creating Variables
  # Creating probability variable from area
  probAid_area <- as.matrix(roiData[roi.area] / sum(roiData[roi.area]))
  probAid <- as.matrix(roiData[roi.prob.aid] / sum(roiData[roi.prob.aid]))
  
  # Precision Code as list, not variable name
  aid.precision.code=aidData[,aid.precision.code]
      
  # Creating matrix of original precision codes
  precision.code.original <- aid.precision.code
  
  ##### Calculating Maximum Lambda Denominator #####
  aid.precision.code <- rep(6,length(aid.precision.code))  
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAid_area, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)  
  param_set[param_set!=0] <- 1
  param_set <- param_set * as.matrix(roiData[roi.area])
  maxLambda_denom <- sum(colSums(param_set))
  
  ##### Calculate Naive Model and Naive Lambda #####
  aid.precision.code <- precision.code.original 
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAid_area, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)  
  
  # Update aid variable
  if(binary){
    roiData[geoSIMEXvariable] <- prob_aid(param_set=param_set)
  } else{
    roiData[geoSIMEXvariable] <- dollar_expected_value(param_set=param_set, aid.project.amount=as.matrix(aidData[aid.project.amount]))
  }
  
  # Update Model
  df.temp <- model$model
  df.temp[,geoSIMEXvariable] <- roiData[,geoSIMEXvariable]
  
  # Update Model
  model_naive <- update(model, data = df.temp)
  
  # Calculate Lambda
  lambda_naive <- calcLambda(param_set, maxLambda_denom, as.matrix(roiData[roi.area]))
    
  ##### Simulating Data with Additional Error and Putting in One Dataframe #####
  prob.increase.list <- runif(iterations,0,1)
  
  status.bar = TRUE
  if(status.bar){
  geoSimulateError.results <- lapply_pb(prob.increase.list, geoSimulateError, 
                                       aidData=aidData, 
                                       roiData=roiData, 
                                       probAidAssume=probAid, 
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
                                       aid.project.amount=aid.project.amount,
                                       model=model,
                                       geoSIMEXvariable=geoSIMEXvariable,
                                       binary=binary,
                                       sim_pc1=sim_pc1)
  } else {
    geoSimulateError.results <- mclapply(prob.increase.list, geoSimulateError, 
                                         aidData=aidData, 
                                         roiData=roiData, 
                                         probAidAssume=probAid, 
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
                                         aid.project.amount=aid.project.amount,
                                         model=model,
                                         geoSIMEXvariable=geoSIMEXvariable,
                                         binary=binary,
                                         sim_pc1=sim_pc1,
                                         mc.cores=mc.cores)
  }
  
  geoSimulateError.results.df <- matrix(NA, ncol=ncol(geoSimulateError.results[[1]][["model.SIMEX.coefs"]]),nrow=0)
  geoSimulateError.results.df <- as.data.frame(geoSimulateError.results.df)
  
  geoSimulateError.results.df.se <- matrix(NA, ncol=ncol(geoSimulateError.results[[1]][["model.SIMEX.se"]]),nrow=0)
  geoSimulateError.results.df.se <- as.data.frame(geoSimulateError.results.df.se)
  
  for(i in 1:length(geoSimulateError.results)){
    geoSimulateError.results.df <- rbind(geoSimulateError.results.df, geoSimulateError.results[[i]]$model.SIMEX.coefs)
    geoSimulateError.results.df.se <- rbind(geoSimulateError.results.df.se, geoSimulateError.results[[i]]$model.SIMEX.se)
  }
  
  ##### Extrapolated Mean Coefficient #####  
  
  # Putting Into Bins
  extrapolatedMean.df <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(geoSimulateError.results.df)))
  
  minLambda <- min(geoSimulateError.results.df$lambda)
  maxLambda <- max(geoSimulateError.results.df$lambda)
  
  binSize <- (maxLambda - minLambda) / bins
  
  binSize_lb <- minLambda
  binSize_ub <- minLambda + binSize
  
  numIter <- bins
  for(i in 1:numIter){
    
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
    
    if(extrapolation=="linear"){
      coef.geoSIMEX[i] <- summary(lm(as.matrix(extrapolatedMean.df[i]) ~ lambda, data = extrapolatedMean.df))$coefficients[1]
    }
    
    if(extrapolation=="quadratic"){
      coef.geoSIMEX[i] <- summary(lm(as.matrix(extrapolatedMean.df[i]) ~ lambda + lambda_sq, data = extrapolatedMean.df))$coefficients[1]
    }
    
  }
  coef.geoSIMEX <- as.data.frame(coef.geoSIMEX)
  names(coef.geoSIMEX) <- head(names(extrapolatedMean.df), -2)
  
  ##### Bootstrap Standard Error #####
  numberBootIter <- nrow(geoSimulateError.results.df)
  
  geoSimulateError.results.df$lambda_sq <- geoSimulateError.results.df$lambda^2
  geoSimulateError.results.df.se$lambda_sq <- geoSimulateError.results.df.se$lambda^2
    
  bootIter.list <- mclapply(seq(1:numberBootIter), bootIter, 
                            geoSimulateError.results.df=geoSimulateError.results.df, 
                            geoSimulateError.results.df.se=geoSimulateError.results.df.se,
                            bins=bins, 
                            number.from.bin=number.from.bin,
                            mc.cores=mc.cores,
                            extrapolation=extrapolation,
                            number.from.bin.average=number.from.bin.average)
  
  boot.coefs.matrix <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(bootIter.list[[1]]$coefs.geoSIMEX.boot)))
  boot.se.matrix <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(bootIter.list[[1]]$se.geoSIMEX.boot)))
  
  for(i in 1:length(bootIter.list)){
    boot.coefs.matrix <- rbind(boot.coefs.matrix, bootIter.list[[i]]$coefs.geoSIMEX.boot)
    boot.se.matrix <- rbind(boot.se.matrix, bootIter.list[[i]]$se.geoSIMEX.boot)
  }
  
  calcSE <- function(i, boot.coefs.matrix, boot.se.matrix){
    se <- sqrt(sum((1/length(boot.coefs.matrix[,i]))*(boot.se.matrix[,i]^2 + (boot.coefs.matrix[,i] - mean(boot.coefs.matrix[,i]))^2)))
    return(se)
  }
  
  se.geoSIMEX <- lapply(1:ncol(boot.coefs.matrix), calcSE, boot.coefs.matrix=boot.coefs.matrix, boot.se.matrix=boot.se.matrix) 
  se.geoSIMEX <- as.data.frame(se.geoSIMEX)
  names(se.geoSIMEX) <- names(boot.coefs.matrix)
  
  ##### Standard Error by Parts #####
  model.var <- colMeans(boot.se.matrix^2)
  
  calcVar.imprecision <- function(i, boot.coefs.matrix){
    var <- sum((1/nrow(boot.coefs.matrix)) * (boot.coefs.matrix[,i] - mean(boot.coefs.matrix[,i]))^2)
    return(var)
  }
  
  var.imprecision <- lapply(1:ncol(boot.coefs.matrix), calcVar.imprecision, boot.coefs.matrix=boot.coefs.matrix)
  var.imprecision <- as.data.frame(var.imprecision)
  names(var.imprecision) <- names(boot.coefs.matrix)
  
  #sqrt(model.var + var.imprecision)
  #se.geoSIMEX
  
  ##### Collecting Results #####
  row.names(coef.geoSIMEX) <- "Coefficients"
  coef <- t(coef.geoSIMEX)
  df = model$df
  
  # Values so will work with stargazer
  # https://github.com/cran/stargazer/blob/master/R/stargazer-internal.R
  # NOTE: Fix so will reflect geoSIMEX coefficients and not naive coefficients?
  #residuals = model$residuals 
  residuals <- rep(NA, length(model$residuals))
  
  model.imprecision.variance <- rbind(model.var,var.imprecision)
  row.names(model.imprecision.variance) <- c("Model Variance", "Imprecision Variance")
  model.imprecision.variance <- as.data.frame(model.imprecision.variance)
  
  return(list(coefficients=coef.geoSIMEX,
              StdErr=se.geoSIMEX,
              df = df,
              variance.model.imprecision = model.imprecision.variance,
              residuals = residuals,
              naive.model = model_naive,
              lambda_naive = lambda_naive,
              simulations = iterations,
              geoSIMEXvariable = geoSIMEXvariable,
              values=geoSimulateError.results.df,
              values.se=geoSimulateError.results.df.se,
              valuesMean=extrapolatedMean.df))
}

#' @title Geographic SIMEX
#' 
#' @author AidData
#' 
#' @import parallel
#'
#' @export
#' 
#' @description
#' \code{geoSIMEX} Implementation of the geoSIMEX algorithm for models 
#' with spatial uncertainty. Package built to work with data from 
#' AidData's data extraction tool.
#' 
#' @usage
#' geoSIMEX(model, geoSIMEXvariable, roiData, aidData, aid.project.amount, 
#' iterations=500, bins=3, fitting.method = "quadratic", roi.area="area",  
#' roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", 
#' roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc5.name="pc6.id",  
#' aid.pc1.centroid.name="centroid.pc1.id", aid.precision.code="precision.code",
#' parallel=TRUE, mc.cores=2)
#' 
#' @param model the naive model
#' @param geoSIMEXvariable character containing the name of the variable with spatial uncertainty
#' @param aidData name of dataframe of aid project data
#' @param aid.project.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
#' @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roiData name of dataframe of ROI data 
#' @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.prob.aid character containing the name of the variable in the ROI dataset which contains the probability of an ROI receiving aid.
#' @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 2 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 3 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 4 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 6 and 8 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
#' @param iterations number of simulated error iterations
#' @param bins number of bins to group coefficients
#' @param fitting.method fitting method for the extrapolation. linear and quadratic are implemented.
#' @param mc.cores number of cores to use for parallelization
#' 
#' @details The values contained within roi.pc1.name and aid.pc1.centroid.name variables should be the same.
#' 
#' @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website]. 
#' 
#' Just like the lm() and glm() packages, geoSIMEX() is equipped to work with stargazer.
#' 
#' @references Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. Journal of American Statistical Association, 89, 1314 – 1328
#' 
#' @examples 
#' library(devtools)
#' install_github("itpir/geoSIMEX")
#' 
#' set.seed(500)
#' 
#' ##### Generating Country-Level Dataset #####
#' numSubcounties <- 120
#' numSubcountyInCounty <- 2 
#' numCountyInDistrict <- 3
#' numDistrictInRegion <- 2
#'
#' N <- numSubcounties
#' subcounty <- 1:N
#' county <- rep(1:(N/numSubcountyInCounty), each=numSubcountyInCounty)
#' district <- rep(1:(N/(numSubcountyInCounty*numCountyInDistrict)), each=(numSubcountyInCounty*numCountyInDistrict))
#' region <- rep(1:(N/(numSubcountyInCounty*numCountyInDistrict*numDistrictInRegion)), each=(numSubcountyInCounty*numCountyInDistrict*numDistrictInRegion))
#' country <- 1
#'
#' subcountyArea <- runif(N)
#' probAid_assumed <- runif(N)
#'
#' subcountyData <- as.data.frame(cbind(subcounty,county,district,region,country,probAid_assumed,subcountyArea))
#'
#' ##### Creating Aid Dataset #####
#' numberProjects = 50
#' aidData <- as.data.frame(matrix(NA,nrow=numberProjects,ncol=3))
#' names(aidData) <- c("aid","trueSubcounty","PC")
#' aidData$aid <- runif(nrow(aidData)) * 100
#' probAid_true <- runif(N)
#' aidData$trueSubcounty <- sample(size=numberProjects,x=c(1:N), prob=probAid_true, replace=TRUE)
#' aidData$PC  <- sample(size=numberProjects, x=c(1,2,3,4,6), prob=runif(5), replace=TRUE)
#'
#' # True Aid
#' aidData$PC.1s <- 1
#' subcountyData$trueAid <- expected_aid_ROI(aidData=aidData, 
#'                                          roiData=subcountyData, 
#'                                          roi.prob.aid="probAid_assumed", 
#'                                          aid.project.amount="aid", 
#'                                          aid.precision.code="PC.1s", 
#'                                          roi.pc1.name="subcounty", 
#'                                          roi.pc2.name="county", 
#'                                          roi.pc3.name="district", 
#'                                          roi.pc4.name="region", 
#'                                          roi.pc5.name="region", 
#'                                          roi.pc6.name="country", 
#'                                          aid.pc1.centroid.name="trueSubcounty")
#'
#' # Wealth - 1 to 1 relation with aid
#' subcountyData$wealth <- subcountyData$trueAid + runif(nrow(subcountyData))
#'
#' # Expected Value Aid 
#' subcountyData$expectedAid <- expected_aid_ROI(aidData=aidData, 
#'                                              aid.project.amount="aid", 
#'                                              aid.precision.code="PC", 
#'                                              aid.pc1.centroid.name="trueSubcounty",
#'                                              roiData=subcountyData, 
#'                                              roi.prob.aid="probAid_assumed", 
#'                                              roi.pc1.name="subcounty", 
#'                                              roi.pc2.name="county", 
#'                                              roi.pc3.name="district", 
#'                                              roi.pc4.name="region", 
#'                                              roi.pc5.name="region", 
#'                                              roi.pc6.name="country")
#' 
#' naive_model <- lm(wealth ~ expectedAid, data=subcountyData)
#' 
#' geoSIMEX_model <- geoSIMEX(model = naive_model, 
#'                           geoSIMEXvariable = "expectedAid", 
#'                           aidData = aidData, 
#'                           aid.project.amount = "aid",
#'                           aid.pc1.centroid.name="trueSubcounty", 
#'                           aid.precision.code="PC",
#'                           roiData = subcountyData, 
#'                           roi.area = "subcountyArea", 
#'                           roi.prob.aid = "probAid_assumed", 
#'                           roi.pc1.name="subcounty", 
#'                           roi.pc2.name="county", 
#'                           roi.pc3.name="district", 
#'                           roi.pc4.name="region", 
#'                           roi.pc5.name="region", 
#'                           roi.pc6.name="country")
#' 
#' ##### Analyzing Results #####
#' summary(naive_model)
#' summary(geoSIMEX_model)
#' plot(geoSIMEX_model, variable="expectedAid")
geoSIMEX <- function(x, ...) UseMethod("geoSIMEX")

geoSIMEX.default <- function(model, 
                             geoSIMEXvariable, 
                             roiData, 
                             aidData, 
                             aid.project.amount,
                             iterations=500, 
                             bins=4, 
                             number.from.bin=2,
                             number.from.bin.average=TRUE,
                             fitting.method = "quadratic", 
                             roi.area="Shape_Area",  
                             roi.prob.aid="Shape_Area",
                             roi.pc1.name="NAME_3", 
                             roi.pc2.name="NAME_2", 
                             roi.pc3.name="NAME_1", 
                             roi.pc4.name="NAME_1", 
                             roi.pc5.name="NAME_1", 
                             roi.pc6.name="NAME_0",  
                             aid.pc1.centroid.name="NAME_3", 
                             aid.precision.code="precision_code",
                             binary=FALSE,
                             sim_pc1=TRUE,
                             parallel=TRUE, 
                             extrapolation="quadratic",
                             mc.cores=1,
                             status.bar = TRUE){
  
  est <- geoSIMEX_est(model=model, 
                      geoSIMEXvariable=geoSIMEXvariable, 
                      roiData=roiData, 
                      aidData=aidData, 
                      aid.project.amount=aid.project.amount,
                      iterations=iterations, 
                      bins=bins, 
                      number.from.bin=number.from.bin,
                      number.from.bin.average=number.from.bin.average,
                      fitting.method = fitting.method, 
                      roi.area=roi.area, 
                      roi.prob.aid=roi.prob.aid,
                      roi.pc1.name=roi.pc1.name, 
                      roi.pc2.name=roi.pc2.name, 
                      roi.pc3.name=roi.pc3.name, 
                      roi.pc4.name=roi.pc4.name, 
                      roi.pc5.name=roi.pc5.name, 
                      roi.pc6.name=roi.pc6.name,  
                      aid.pc1.centroid.name=aid.pc1.centroid.name, 
                      aid.precision.code=aid.precision.code,
                      binary=binary,
                      sim_pc1=sim_pc1,
                      parallel=parallel, 
                      extrapolation=extrapolation,
                      mc.cores=mc.cores)
  
  est$call <- model$call
  
  class(est) <- "geoSIMEX"
  est
}

print.geoSIMEX <- function(x, ...){
  cat("Naive Model:\n")
  print(x$call)
  cat("\ngeoSIMEX-Variables: ", x$geoSIMEXvariable, sep="")
  cat("\nNumer of Simulations: ", x$simulations, sep= "")
  cat("\n\nCoefficients:\n")
  print(x$coefficients, row.names=FALSE)  
}

summary.geoSIMEX <- function(object, ...){
  
  Coef <- as.numeric(object$coefficients)
  StdErr <- as.numeric(object$StdErr)
  tval <- Coef / StdErr
  
  TAB <- cbind(Estimate = Coef,
               StdErr = StdErr,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  
  TAB <- as.data.frame(TAB)
  row.names(TAB) <- names(object$coefficients)
  names(TAB) <- c("Estimate","Std. Error","t.value","p.value")
  
  res <- list(call=object$call,
              coefficients=TAB)
  
  class(res) <- "summary.geoSIMEX"
  res
}

print.summary.geoSIMEX <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

plot.geoSIMEX <- function(x, variable, confInt = 95, allSimulations=TRUE, includeTitle=TRUE, name_variable="", ylim="default"){
  
  #name_variable = variable
  # NOTE: Can make things like "all simulations", "confidence bands", "% interval"
  # all parameters / options.
  
  if(name_variable == ""){
    variable_name = variable
  } else {
    variable_name = name_variable
  }
  
  # Calculating 95\% Confidence Bands
  t_critical <- qt((1-confInt/100), x$df)
  
  coef.geoSIMEX <- x$coefficients[,variable]
  se.geoSIMEX <- x$StdErr[,variable]
  CI.geoSIMEX <- as.data.frame(matrix(NA, nrow=2, ncol=2))
  names(CI.geoSIMEX) <- c("value", "lambda")
  CI.geoSIMEX$lambda <- 0
  CI.geoSIMEX$value[1] <- coef.geoSIMEX - t_critical*se.geoSIMEX
  CI.geoSIMEX$value[2] <- coef.geoSIMEX + t_critical*se.geoSIMEX
  
  coef.naive <- as.data.frame(t(x$naive.model$coefficients))[,variable]
  se.naive <- as.data.frame(summary(x$naive.model)$coefficients)[variable,2]
  CI.naive <- as.data.frame(matrix(NA, nrow=2, ncol=2))
  names(CI.naive) <- c("value", "lambda")
  CI.naive$lambda <- x$lambda_naive
  CI.naive$value[1] <- coef.naive - t_critical*se.naive
  CI.naive$value[2] <- coef.naive + t_critical*se.naive
  
  # All Variable Values to Define Min / Max ylim of Plot
  var.vales.all <- c(x$valuesMean[,variable], x$values[,variable], x$coefficients[,variable], CI.geoSIMEX$value[1], CI.geoSIMEX$value[2], CI.naive$value[1], CI.naive$value[2])
  
  title = ""
  
  if(includeTitle){
    
    title = paste("geoSIMEX Plot of ",variable_name, sep="")
    
  }
  

  
  # Plot Mean Lambda Values
  ylim==c(NA,NA)
  

  
  if(ylim[1] == "default"){
    
    plot(x$valuesMean$lambda, x$valuesMean[,variable],
         xlab = expression((lambda)),
         ylab = paste(variable_name, " coefficient", sep=""),
         main = title,
         xlim = c(0,1),
         ylim = c(min(var.vales.all), max(var.vales.all)),
         pch  = 16,
         col="red")
    
  } else{
    plot(x$valuesMean$lambda, x$valuesMean[,variable],
         xlab = expression((lambda)),
         ylab = paste(variable_name, " coefficient", sep=""),
         main = title,
         xlim = c(0,1),
         ylim = ylim,
         pch  = 16,
         col="red")
  }
  
  if(allSimulations){
    points(x$values$lambda, x$values[,variable],
           pch=19, cex = .15)
  }
  
  points(0,x$coefficients[,variable],
         pch=1)
  
  # Prepping Data for Lines
  var <- x$valuesMean[,variable]
  lambda <- x$valuesMean$lambda
  fit <- lm(var ~ lambda + I(lambda^2))
  
  newdat_sim = data.frame(lambda = seq(min(lambda), max(lambda), length.out = 100))
  newdat_sim$pred = predict(fit, newdata = newdat_sim)
  
  newdat_extrp = data.frame(lambda = seq(0, min(lambda), length.out = 100))
  newdat_extrp$pred = predict(fit, newdata = newdat_extrp)
  
  # Plotting Lines
  with(newdat_sim, lines(x = lambda, y = pred))
  with(newdat_extrp, lines(x = lambda, y = pred, lty = 2))
  
  # Plot Naive Model
  points(x$lambda_naive, coef.naive,
         pch="*")
  
  # Plotting Confidence Intervals
  lines(x=CI.geoSIMEX$lambda, y=CI.geoSIMEX$value, lty=5)
  lines(x=CI.naive$lambda, y=CI.naive$value, lty=1)
  
}

##### SECTION 1A: MODEL AVERAGING - RAND PROB #####

modelAverageRandProb_est <- function(iterations, 
                                     model, 
                                     geoSIMEXvariable, 
                                     roiData, 
                                     aidData, 
                                     aid.project.amount, 
                                     roi.pc1.name, 
                                     roi.pc2.name, 
                                     roi.pc3.name, 
                                     roi.pc4.name, 
                                     roi.pc5.name, 
                                     roi.pc6.name, 
                                     aid.pc1.centroid.name,
                                     aid.precision.code,
                                     binary, 
                                     parallel, 
                                     mc.cores){
  
  aid.precision.code <- aidData[[aid.precision.code]]
  
  param_set <- paramSet(aidData=aidData, roiData=roiData, probAidAssume=runif(nrow(roiData)), aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set.bin <- (param_set > 0)*1
  
  model.list <- mclapply(1:iterations, model_rand_prob, param_set.bin=param_set.bin, aidData=aidData, roiData=roiData, aid.project.amount=aid.project.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary)
  
  coef.df <- model.list[[1]]$model.SIMEX.coefs
  se.df <- model.list[[1]]$model.SIMEX.se
  for(i in 2:length(model.list)){
    coef.df <- rbind(coef.df, model.list[[i]]$model.SIMEX.coefs)
    se.df <- rbind(se.df, model.list[[i]]$model.SIMEX.se)
  }
  
  coef <- colSums(coef.df) / nrow(coef.df)
  
  var.df <- se.df^2
  
  se <- lapply(1:length(coef), function(i) sqrt(sum((1/iterations)*(var.df[,i] + (coef.df[,i] - coef[i])^2))))
  se <- as.data.frame(se)
  names(se) <- names(coef)
  
  coef <- as.data.frame(t(coef))
  row.names(coef) <- "Coefficients"
  row.names(se) <- "Std. Error"
  
  return(list(coefficients=coef,
              se=se,
              df=model$df,
              coef.all=coef.df,
              se.all=se.df,
              geoSIMEXvariable=geoSIMEXvariable,
              iterations=iterations,
              residuals = model$residuals))
}

modelAverageRandProb <- function(x, ...) UseMethod("modelAverageRandProb")

modelAverageRandProb.default <- function(iterations=500, 
                                         model,
                                         geoSIMEXvariable,
                                         roiData, 
                                         aidData,
                                         aid.project.amount, 
                                         roi.pc1.name="pc1.id", 
                                         roi.pc2.name="pc2.id", 
                                         roi.pc3.name="pc3.id", 
                                         roi.pc4.name="pc4.id", 
                                         roi.pc5.name="pc5.id", 
                                         roi.pc6.name="pc6.id", 
                                         aid.pc1.centroid.name="centroid.pc1.id",
                                         aid.precision.code="precision.code",
                                         binary=FALSE,
                                         parallel=TRUE, 
                                         mc.cores=1){
  
  est <- modelAverageRandProb_est(iterations=iterations, 
                                  model=model,
                                  geoSIMEXvariable=geoSIMEXvariable,
                                  roiData=roiData, 
                                  aidData=aidData,
                                  aid.project.amount=aid.project.amount, 
                                  roi.pc1.name=roi.pc1.name, 
                                  roi.pc2.name=roi.pc2.name, 
                                  roi.pc3.name=roi.pc3.name, 
                                  roi.pc4.name=roi.pc4.name, 
                                  roi.pc5.name=roi.pc5.name, 
                                  roi.pc6.name=roi.pc6.name, 
                                  aid.pc1.centroid.name=aid.pc1.centroid.name,
                                  aid.precision.code=aid.precision.code,
                                  binary=binary,
                                  parallel=parallel, 
                                  mc.cores=mc.cores)
  
  est$call <- model$call
  
  class(est) <- "modelAverageRandProb"
  est
}

print.modelAverageRandProb <- function(x, ...){
  cat("Naive Model:\n")
  print(x$call)
  cat("\ngeoSIMEX-Variables: ", x$geoSIMEXvariable, sep="")
  cat("\nNumer of Iterations: ", x$iterations, sep= "")
  cat("\n\nCoefficients:\n")
  print(x$coefficients, row.names=FALSE)  
}

summary.modelAverageRandProb<- function(object, ...){
  
  Coef <- as.numeric(object$coefficients)
  StdErr <- as.numeric(object$se)
  tval <- Coef / StdErr
  
  TAB <- cbind(Estimate = Coef,
               StdErr = StdErr,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  
  TAB <- as.data.frame(TAB)
  row.names(TAB) <- names(object$coefficients)
  names(TAB) <- c("Estimate","Std. Error","t.value","p.value")
  
  res <- list(call=object$call,
              coefficients=TAB)
  
  class(res) <- "summary.modelAverageRandProb"
  res
}

print.summary.modelAverageRandProb <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}


##### SECTION 2: FUNCTIONS TO BE INCLUDED IN PACKAGE #####
#' Calculates expected value of aid for each ROI
#' 
#' @title Expected Value Aid
#' 
#' @author AidData
#' 
#' @import parallel
#'
#' @export
#' 
#' @description
#' \code{expected_aid_ROI} Calculates expected value of aid.
#' 
#' @usage
#' expected_aid_ROI(aidData, roiData, probAidAssume, aid.project.amount, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
#'  
#' @param aidData name of dataframe of aid project data
#' @param aid.project.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
#' @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roiData name of dataframe of ROI data 
#' @param roi.prob.aid character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
#'  
#' @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website].
#' 
#' @examples 
#' library(devtools)
#' install_github("itpir/geoSIMEX")
#' 
#' set.seed(42)
#' 
#' ##### Generating Country-Level Dataset #####
#' numSubcounties <- 120
#' numSubcountyInCounty <- 2 
#' numCountyInDistrict <- 3
#' numDistrictInRegion <- 2
#' 
#' N <- numSubcounties
#' subcounty <- 1:N
#' county <- rep(1:(N/numSubcountyInCounty), each=numSubcountyInCounty)
#' district <- rep(1:(N/(numSubcountyInCounty*numCountyInDistrict)), each=(numSubcountyInCounty*numCountyInDistrict))
#' region <- rep(1:(N/(numSubcountyInCounty*numCountyInDistrict*numDistrictInRegion)), each=(numSubcountyInCounty*numCountyInDistrict*numDistrictInRegion))
#' country <- 1
#' 
#' subcountyArea <- runif(N)
#' probAid_assumed <- runif(N)
#' 
#' subcountyData <- as.data.frame(cbind(subcounty,county,district,region,country,probAid_assumed,subcountyArea))
#' 
#' ##### Creating Aid Dataset #####
#' numberProjects = 50
#' aidData <- as.data.frame(matrix(NA,nrow=numberProjects,ncol=3))
#' names(aidData) <- c("aid","trueSubcounty","PC")
#' aidData$aid <- runif(nrow(aidData)) * 100
#' probAid_true <- runif(N)
#' aidData$trueSubcounty <- sample(size=numberProjects,x=c(1:N), prob=probAid_true, replace=TRUE)
#' aidData$PC  <- sample(size=numberProjects, x=c(1,2,3,4,6), prob=runif(5), replace=TRUE)
#' 
#' ##### Calculating Expected Aid #####
#' subcountyData$expectedAid <- expected_aid_ROI(aidData=aidData, 
#'                                              aid.project.amount="aid", 
#'                                              aid.precision.code="PC", 
#'                                              aid.pc1.centroid.name="trueSubcounty",
#'                                              roiData=subcountyData, 
#'                                              roi.prob.aid="probAid_assumed", 
#'                                              roi.pc1.name="subcounty", 
#'                                              roi.pc2.name="county", 
#'                                              roi.pc3.name="district", 
#'                                              roi.pc4.name="region", 
#'                                              roi.pc5.name="region", 
#'                                              roi.pc6.name="country")
expected_aid_ROI <- function(aidData, roiData, roi.prob.aid, aid.project.amount, aid.precision.code, roi.pc1.name="ID_3", roi.pc2.name="ID_2", roi.pc3.name="ID_1", roi.pc4.name="ID_1", roi.pc5.name="ID_1", roi.pc6.name="ID_0", aid.pc1.centroid.name="ID_3"){

  aid.precision.code <- aidData[,aid.precision.code]   
  probAid_area <- as.matrix(roiData[,roi.prob.aid] / sum(roiData[,roi.prob.aid]))
  dollar_set_values <- aidData[,aid.project.amount]
  
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAid_area, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  aid <- dollar_expected_value(param_set=param_set, aid.project.amount=dollar_set_values)
  return(aid)
}

prob_aid_ROI <- function(aidData, roiData, probAidAssume, aid.project.amount, aid.precision.code, roi.pc1.name, roi.pc2.name, roi.pc3.name, roi.pc4.name, roi.pc5.name, roi.pc6.name, aid.pc1.centroid.name){
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set_probNoAid <- 1 - param_set
  prob.0.prjs <- apply(param_set_probNoAid, 1, function(x) prod(x))
  prob.atLeast1.prjs <- 1 - prob.0.prjs
  
  #prob.atLeast1.prjs[prob.atLeast1.prjs < 0.5] <- 0
  #prob.atLeast1.prjs[prob.atLeast1.prjs >= 0.5] <- 1
  return(prob.atLeast1.prjs)
}

#' @title Spatial Uncertainty (Lambda)
#' 
#' @author AidData
#' 
#' @import parallel
#'
#' @export
#' 
#' @description
#' \code{calc_lambda} Calculates spatial uncertainty (lambda) of aid project dataset
#' 
#' @usage
#' calc_lambda(aidData, roiData, roi.area="area", aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
#'  
#' @param roiData name of dataframe of ROI data 
#' @param aidData name of dataframe of aid project data
#' @param aid.project.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
#' @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
#' @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
#' @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#'  
#' @examples 
#' library(devtools)
#' install_github("itpir/geoSIMEX")
#' set.seed(42)
#' 
#' ##### Generating Country-Level Dataset #####
#' numSubcounties <- 120
#' numSubcountyInCounty <- 2 
#' numCountyInDistrict <- 3
#' numDistrictInRegion <- 2
#'
#' N <- numSubcounties
#' subcounty <- 1:N
#' county <- rep(1:(N/numSubcountyInCounty), each=numSubcountyInCounty)
#' district <- rep(1:(N/(numSubcountyInCounty*numCountyInDistrict)), each=(numSubcountyInCounty*numCountyInDistrict))
#' region <- rep(1:(N/(numSubcountyInCounty*numCountyInDistrict*numDistrictInRegion)), each=(numSubcountyInCounty*numCountyInDistrict*numDistrictInRegion))
#' country <- 1
#'
#' subcountyArea <- runif(N)
#' probAid_assumed <- runif(N)
#'
#' subcountyData <- as.data.frame(cbind(subcounty,county,district,region,country,probAid_assumed,subcountyArea))
#'
#' ##### Creating Aid Dataset #####
#' numberProjects = 50
#' aidData <- as.data.frame(matrix(NA,nrow=numberProjects,ncol=3))
#' names(aidData) <- c("aid","trueSubcounty","PC")
#' aidData$aid <- runif(nrow(aidData)) * 100
#' probAid_true <- runif(N)
#' aidData$trueSubcounty <- sample(size=numberProjects,x=c(1:N), prob=probAid_true, replace=TRUE)
#' aidData$PC  <- sample(size=numberProjects, x=c(1,2,3,4,6), prob=runif(5), replace=TRUE)
#'
#' lambda <- calc_lambda(aidData=aidData, 
#'                      roiData=subcountyData, 
#'                      roi.area="subcountyArea", 
#'                      aid.precision.code="PC", 
#'                      roi.pc1.name="subcounty", 
#'                      roi.pc2.name="county", 
#'                      roi.pc3.name="district", 
#'                      roi.pc4.name="region", 
#'                      roi.pc5.name="region", 
#'                      roi.pc6.name="country", 
#'                      aid.pc1.centroid.name="trueSubcounty")
#' lambda
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
  #aid.precision.code <- aidData[[aid.precision.code]]
  pc.original <- aidData[,aid.precision.code] 
  area <- roiData[,roi.area]
  probAidAssume <- area / sum(area)

  # Calculaitng Max Lambda
  pc.6.vector <- rep(6, nrow(aidData))
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=pc.6.vector, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set[param_set!=0] <- 1
  param_set <- param_set*area
  maxLambda_denom <- sum(colSums(param_set))
  
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=pc.original, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set[param_set!=0] <- 1
  param_set <- param_set*area
  lambda <- sum(colSums(param_set)) / maxLambda_denom
  return(lambda)
}

subset.aiddata <- function(json){
  
  #### Import Dataset
  # Eventually all datasets will be imported into R, so can just use this.
  #geo.data <- eval(parse(text=as.character(json$release_data[5])))
  # https://github.com/chrisalbon/code_r/blob/master/download-and-unzip-data.r
  # http://stackoverflow.com/questions/23899525/using-r-to-download-zipped-data-file-extract-and-import-csv  
  # http://stackoverflow.com/questions/32647779/r-download-and-unzip-file-to-store-in-data-frame
  # https://github.com/AidData-WM/public_datasets/tree/master/geocoded
  # "https://github.com/AidData-WM/public_datasets/geocoded/UgandaAIMS_GeocodedResearchRelease_Level1_v1.4.1.zip"
  
  geo.data.name <- json$release_data$dataset
  # Modifying geo.data.name (geocoded dataset name) to match the name in github:
  # 1. Changing to lowercase
  # 2. Replacing "_" with "."
  
  geo.data.name <- tolower(geo.data.name)
  geo.data.name <- gsub("_",".",geo.data.name, fixed=TRUE)
  
  # Downloading dataset from github
  github.datasets <- read_html("https://github.com/AidData-WM/public_datasets/tree/master/geocoded")
  github.datasets <- html_text(github.datasets)
  github.datasets <- strsplit(github.datasets, "\n")
  github.datasets <- github.datasets[[1]]
  github.datasets <- gsub(" ", "", github.datasets, fixed = TRUE)
  github.datasets.modified <- gsub("_", ".", github.datasets, fixed = TRUE)
  github.datasets.modified <- tolower(github.datasets.modified)
  github.dataset.name.zip <- github.datasets[github.datasets.modified %in% paste(geo.data.name, ".zip", sep="")]
  github.dataset.name <- gsub('.{4}$', '', github.dataset.name.zip)
  
  file.path <- paste("https://github.com/AidData-WM/public_datasets/blob/master/geocoded/",github.dataset.name.zip,"?raw=TRUE",sep="")
  
  download.file(file.path, destfile=github.dataset.name.zip)
  unzip(github.dataset.name.zip, exdir=".")
  
  geo.data <- read.csv(paste(github.dataset.name,"/data/level_1a.csv",sep=""))
  
  # Merging GADM names to data  
  geo.data$latitude[is.na(geo.data$latitude)] <- geo.data$latitude[!is.na(geo.data$latitude)][1]
  geo.data$longitude[is.na(geo.data$longitude)] <- geo.data$longitude[!is.na(geo.data$longitude)][1]
  
  coordinates(geo.data) <- ~longitude+latitude
  proj4string(geo.data) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  #### Extracting GADM Data, that matches AidData 
  gadm <- raster::getData('GADM', country=toupper(substr(json$boundary$name,1,3)), level=as.numeric(substr(json$boundary$name,8,8)))
  
  aiddata.gadm.over <- over(geo.data, gadm)
  
  #### Adding ADM Names to Aid Dataset 
  
  adm.level <- as.numeric(substr(json$boundary$name,8,8))
  
  for(i in 0:adm.level){
    
    geo.data[[paste("NAME_",i,sep="")]] <- aiddata.gadm.over[,paste("NAME_",i,sep="")]
    geo.data[[paste("ID_",i,sep="")]] <- aiddata.gadm.over[,paste("ID_",i,sep="")]
  }
  
  geo.data <- geo.data@data
  
  #### Subsetting by Filter
  
  # Number of filters to subset by
  num.filters <- length(json$release_data$filters)
  
  for(i in 1:num.filters){
    
    filter <- names(json$release_data$filters)[i]
    use.filter <- TRUE
    
    if( (filter == "donors") & (json$release_data$filters[i] == "All") ){
      use.filter <- FALSE
    }
    
    if(filter == "transaction_year"){
      
      # Loop through years, add "1" to keep.obs if meets condition to keep observation
      keep.obs <- rep(0,length(geo.data$transactions_start_year))
      for(year in matrix(unlist(json$release_data$filters[i]))[,1]){
        keep.obs <- keep.obs + as.numeric((year >= geo.data$transactions_start_year) & (year <= geo.data$transactions_end_year))
      }
      
      geo.data <- geo.data[keep.obs > 0,]
      
      use.filter <- FALSE
    }
    
    
    # Checking to ensure filter exists in dataset.
    # If doesn't exist, skip.
    # Filter name in JSON must match variable name in aid dataset
    
    if(use.filter){
      
      # Checking to make sure filter exists in dataset (if doesn't, skip). Filter name in JSON must match filter
      if(filter %in% names(geo.data)){ 
        geo.data <- geo.data[geo.data[filter][,1] %in% matrix(unlist(json$release_data$filters[i]))[,1],]
      }
      
    }
    
  }
  
  # Removing NA rows.
  geo.data <- geo.data[!is.na(geo.data$project_id),]
  
  return(geo.data)
}





##### SECTION 3: DEFINING FUNCTIONS THAT OTHER FUNCTIONS USE #####

# https://ryouready.wordpress.com/2010/01/11/progress-bars-in-r-part-ii-a-wrapper-for-apply-functions/
lapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}

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
  PC_temp <- aid.precision.code[i]
  
  # Getting information about subcounty & region that the aid project falls in
  subcounty_temp <- roiData[roiData[,roi.pc1.name] == aidDataPrj_Temp[,aid.pc1.centroid.name],]
  
  # Making temporary dataset that will merge into
  paramCol <- roiData
  paramCol$prj <- 0
  
  if(PC_temp == 1){
    paramCol$prj[paramCol[,roi.pc1.name] %in% subcounty_temp[,roi.pc1.name]] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  if(PC_temp == 2){
    paramCol$prj[paramCol[,roi.pc2.name] %in% subcounty_temp[,roi.pc2.name]] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  if(PC_temp == 3){
    paramCol$prj[paramCol[,roi.pc3.name] %in% subcounty_temp[,roi.pc3.name]] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  if(PC_temp == 4){
    paramCol$prj[paramCol[,roi.pc4.name] %in% subcounty_temp[,roi.pc4.name]] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  if((PC_temp == 6) | (PC_temp == 8)){
    paramCol$prj[paramCol[,roi.pc6.name] %in% subcounty_temp[,roi.pc6.name]] <- 1
    paramCol$prj <- paramCol$prj * probAidAssume
    paramCol$prj <- paramCol$prj / sum(paramCol$prj)  
  }
  
  #row.names(paramCol) <- paramCol[,roi.pc1.name]
  row.names(paramCol) <- paste(1:length(paramCol[,roi.pc1.name]), paramCol[,roi.pc1.name])
  
  paramCol_return <- as.data.frame(paramCol$prj)
  row.names(paramCol_return) <- row.names(paramCol)
  names(paramCol_return) <- paste("prj",i,sep="")
  
  return(paramCol_return)
}

paramSet <- function(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name){
  param_set = as.data.frame(lapply(1:nrow(aidData), paramCol, aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name))
  
  #param_set.a = lapply(1:nrow(aidData), paramCol, aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  
  return(param_set) 
}

dollar_expected_value <- function(param_set, aid.project.amount){    
  aid.project.amount[is.na(aid.project.amount)] <- 0
  paramDollars <- lapply(1:length(aid.project.amount), function(i) param_set[,i] * aid.project.amount[i])
  paramDollars <- as.data.frame(paramDollars)
  dollars_expected <- rowSums(paramDollars)
  return(dollars_expected)
}

prob_aid <- function(param_set){    
  param_set_probNoAid <- 1 - param_set
  prob.0.prjs <- apply(param_set_probNoAid, 1, function(x) prod(x))
  prob.atLeast1.prjs <- 1 - prob.0.prjs
  
  #prob.atLeast1.prjs[prob.atLeast1.prjs < 0.5] <- 0
  #prob.atLeast1.prjs[prob.atLeast1.prjs >= 0.5] <- 1
  return(prob.atLeast1.prjs)
}

geoSimulateError <- function(probIncPC, aidData=aidData, roiData=roiData, probAidAssume=probAid, PC_researcherSees=precision.code.original, maxLambda_denom=maxLambda_denom, roi.area=roi.area, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name, aid.project.amount=aid.project.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary, sim_pc1=sim_pc1){
  
  # Initializing Results Matrix
  results <- matrix(NA, ncol=2,nrow=1)
  results <- as.data.frame(results)
  names(results) <- c("lambda","aid.expected.coef")
  
  # Simulating Error
  aid.precision.code <- PC_researcherSees
  probStayPC <- 1 - probIncPC
  
  ### Equal chance of increasing precision codes
  uncertainty.level <- sample(1,x=c(0,1,2,3,4),prob=c(0.20,0.20,0.20,0.20,0.20))
  
  if(uncertainty.level == 0){
    # ... Precision Code 4s
    aid.precision.code[aid.precision.code==4] <- sample(size=length(aid.precision.code[aid.precision.code==4]), 
                                                                            x=c(4,6), 
                                                                            prob=c(probStayPC,probIncPC), 
                                                                            replace=TRUE)
    # ... Precision Code 3s
    aid.precision.code[aid.precision.code==3] <- sample(size=length(aid.precision.code[aid.precision.code==3]), 
                                                                            x=c(3,4,6), 
                                                                            prob=c(probStayPC,probIncPC/2,probIncPC/2), 
                                                                            replace=TRUE)
    
    # ... Precision Code 2s
    aid.precision.code[aid.precision.code==2] <- sample(size=length(aid.precision.code[aid.precision.code==2]), 
                                                                            x=c(2,3,4,6), 
                                                                            prob=c(probStayPC,probIncPC/3,probIncPC/3,probIncPC/3), 
                                                                            replace=TRUE)
    
    # ... Precision Code 1s
    if(sim_pc1){
      aid.precision.code[aid.precision.code==1] <- sample(size=length(aid.precision.code[aid.precision.code==1]), 
                                                                              x=c(1,2,3,4,6), 
                                                                              prob=c(probStayPC,probIncPC/4,probIncPC/4,probIncPC/4,probIncPC/4), 
                                                                              replace=TRUE)
    }
  }
  
  if(uncertainty.level == 1){
    # ... Precision Code 4s
    aid.precision.code[aid.precision.code==4] <- sample(size=length(aid.precision.code[aid.precision.code==4]), 
                                                        x=c(4,6), 
                                                        prob=c(probStayPC,probIncPC), 
                                                        replace=TRUE)
    # ... Precision Code 3s
    aid.precision.code[aid.precision.code==3] <- sample(size=length(aid.precision.code[aid.precision.code==3]), 
                                                        x=c(3,4,6), 
                                                        prob=c(probStayPC,probIncPC/3,probIncPC/2), 
                                                        replace=TRUE)
    
    # ... Precision Code 2s
    aid.precision.code[aid.precision.code==2] <- sample(size=length(aid.precision.code[aid.precision.code==2]), 
                                                        x=c(2,3,4,6), 
                                                        prob=c(probStayPC,probIncPC/5,probIncPC/4,probIncPC/3), 
                                                        replace=TRUE)
    
    # ... Precision Code 1s
    if(sim_pc1){
      aid.precision.code[aid.precision.code==1] <- sample(size=length(aid.precision.code[aid.precision.code==1]), 
                                                          x=c(1,2,3,4,6), 
                                                          prob=c(probStayPC,probIncPC/7,probIncPC/6,probIncPC/5,probIncPC/4), 
                                                          replace=TRUE)
    }
  }
  
  if(uncertainty.level == 2){
    # ... Precision Code 4s
    aid.precision.code[aid.precision.code==4] <- sample(size=length(aid.precision.code[aid.precision.code==4]), 
                                                        x=c(4,6), 
                                                        prob=c(probStayPC/3,probIncPC), 
                                                        replace=TRUE)
    # ... Precision Code 3s
    aid.precision.code[aid.precision.code==3] <- sample(size=length(aid.precision.code[aid.precision.code==3]), 
                                                        x=c(3,4,6), 
                                                        prob=c(probStayPC/3,probIncPC/2,probIncPC/1), 
                                                        replace=TRUE)
    
    # ... Precision Code 2s
    aid.precision.code[aid.precision.code==2] <- sample(size=length(aid.precision.code[aid.precision.code==2]), 
                                                        x=c(2,3,4,6), 
                                                        prob=c(probStayPC/3,probIncPC/5,probIncPC/2,probIncPC/1), 
                                                        replace=TRUE)
    
    # ... Precision Code 1s
    if(sim_pc1){
      aid.precision.code[aid.precision.code==1] <- sample(size=length(aid.precision.code[aid.precision.code==1]), 
                                                          x=c(1,2,3,4,6), 
                                                          prob=c(probStayPC/3,probIncPC/7,probIncPC/6,probIncPC/2,probIncPC/1), 
                                                          replace=TRUE)
    }
  }
  
  if(uncertainty.level == 3){
    # ... Precision Code 4s
    aid.precision.code[aid.precision.code==4] <- sample(size=length(aid.precision.code[aid.precision.code==4]), 
                                                        x=c(4,6), 
                                                        prob=c(probStayPC/15,probIncPC), 
                                                        replace=TRUE)
    # ... Precision Code 3s
    aid.precision.code[aid.precision.code==3] <- sample(size=length(aid.precision.code[aid.precision.code==3]), 
                                                        x=c(3,4,6), 
                                                        prob=c(probStayPC/15,probIncPC/3,probIncPC/1), 
                                                        replace=TRUE)
    
    # ... Precision Code 2s
    aid.precision.code[aid.precision.code==2] <- sample(size=length(aid.precision.code[aid.precision.code==2]), 
                                                        x=c(2,3,4,6), 
                                                        prob=c(probStayPC/15,probIncPC/11,probIncPC/3,probIncPC/1), 
                                                        replace=TRUE)
    
    # ... Precision Code 1s
    if(sim_pc1){
      aid.precision.code[aid.precision.code==1] <- sample(size=length(aid.precision.code[aid.precision.code==1]), 
                                                          x=c(1,2,3,4,6), 
                                                          prob=c(probStayPC/15,probIncPC/11,probIncPC/10,probIncPC/3,probIncPC/1), 
                                                          replace=TRUE)
    }
  }
  
  if(uncertainty.level == 4){
    # ... Precision Code 4s
    aid.precision.code[aid.precision.code==4] <- sample(size=length(aid.precision.code[aid.precision.code==4]), 
                                                        x=c(4,6), 
                                                        prob=c(probStayPC/50,probIncPC), 
                                                        replace=TRUE)
    # ... Precision Code 3s
    aid.precision.code[aid.precision.code==3] <- sample(size=length(aid.precision.code[aid.precision.code==3]), 
                                                        x=c(3,4,6), 
                                                        prob=c(probStayPC/50,probIncPC/5,probIncPC/1), 
                                                        replace=TRUE)
    
    # ... Precision Code 2s
    aid.precision.code[aid.precision.code==2] <- sample(size=length(aid.precision.code[aid.precision.code==2]), 
                                                        x=c(2,3,4,6), 
                                                        prob=c(probStayPC/50,probIncPC/20,probIncPC/5,probIncPC/1), 
                                                        replace=TRUE)
    
    # ... Precision Code 1s
    if(sim_pc1){
      aid.precision.code[aid.precision.code==1] <- sample(size=length(aid.precision.code[aid.precision.code==1]), 
                                                          x=c(1,2,3,4,6), 
                                                          prob=c(probStayPC/50,probIncPC/25,probIncPC/20,probIncPC/5,probIncPC/1), 
                                                          replace=TRUE)
    }
  }
  
  
  # Calculating paramSet
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  
  # Update aid variable
  if(binary){
    roiData[geoSIMEXvariable] <- prob_aid(param_set=param_set)
  } else{
    roiData[geoSIMEXvariable] <- dollar_expected_value(param_set=param_set, aid.project.amount=as.matrix(aidData[aid.project.amount]))
  }
  
  # Update Model
  df.temp <- model$model
  df.temp[geoSIMEXvariable] <- roiData[geoSIMEXvariable]
  
  # Update Model
  model.SIMEX <- update(model, data = df.temp)
  
  # Calculate Lambda
  lambda <- calcLambda(param_set, maxLambda_denom, as.matrix(roiData[roi.area]))
  
  # Collecting Results
  
  # Coefficients
  model.SIMEX.coefs <- as.data.frame(as.list(model.SIMEX$coefficients))
  names(model.SIMEX.coefs) <- names(model.SIMEX$coefficients)
  model.SIMEX.coefs$lambda <- lambda
  
  # Standard Error
  model.SIMEX.se <- as.data.frame(t(summary(model.SIMEX)$coefficients[,2]))
  names(model.SIMEX.se) <- names(model.SIMEX$coefficients)
  model.SIMEX.se$lambda <- lambda
  
  return(list(model.SIMEX.coefs=model.SIMEX.coefs,
              model.SIMEX.se=model.SIMEX.se))
}

geoSimulateError2 <- function(iter, aidData_1, aidData_2, roiData, probAidAssume_1, probAidAssume_2, PC_researcherSees_1, PC_researcherSees_2, maxLambda_denom_1, maxLambda_denom_2, roi.area_1, roi.area_2, aid.precision.code_1, aid.precision.code_2, roi.pc1.name, roi.pc2.name, roi.pc3.name, roi.pc4.name, roi.pc5.name, roi.pc6.name, aid.pc1.centroid.name_1, aid.pc1.centroid.name_2, aid.amount_1, aid.amount_2, model, geoSIMEXvariable_1, geoSIMEXvariable_2, binary, diagonal, sim_pc1){
  
  if(diagonal){
    probIncPC_1 <- runif(1)
    probIncPC_2 <- probIncPC_1
  } else {
    probIncPC_1 <- runif(1)
    probIncPC_2 <- runif(1)
  }
  
  # to make sure not don't get too correlated
  if(probIncPC_1 > .85){
    probIncPC_2 <- 1 - probIncPC_2
  }
  
  # Simulating Error
  aidData_1[aid.precision.code_1] <- PC_researcherSees_1
  aidData_2[aid.precision.code_2] <- PC_researcherSees_2
  
  probStayPC_1 <- 1 - probIncPC_1
  probStayPC_2 <- 1 - probIncPC_2
  
  # Randomly Increase Precision Codes
  
  # ... Precision Code 4s
  aidData_1[aid.precision.code_1][aidData_1[aid.precision.code_1] == 4] <- sample(size=length(aidData_1[aid.precision.code_1][aidData_1[aid.precision.code_1] == 4]), 
                                                                                  x=c(4,6), 
                                                                                  prob=c(probStayPC_1,probIncPC_1), 
                                                                                  replace=TRUE)
  
  aidData_2[aid.precision.code_2][aidData_2[aid.precision.code_2] == 4] <- sample(size=length(aidData_2[aid.precision.code_2][aidData_2[aid.precision.code_2] == 4]), 
                                                                                  x=c(4,6), 
                                                                                  prob=c(probStayPC_2,probIncPC_2), 
                                                                                  replace=TRUE)
  
  # ... Precision Code 3s
  aidData_1[aid.precision.code_1][aidData_1[aid.precision.code_1] == 3] <- sample(size=length(aidData_1[aid.precision.code_1][aidData_1[aid.precision.code_1] == 3]), 
                                                                                  x=c(3,4,6), 
                                                                                  prob=c(probStayPC_1,probIncPC_1/2,probIncPC_1/2), 
                                                                                  replace=TRUE)
  
  aidData_2[aid.precision.code_2][aidData_2[aid.precision.code_2] == 3] <- sample(size=length(aidData_2[aid.precision.code_2][aidData_2[aid.precision.code_2] == 3]), 
                                                                                  x=c(3,4,6), 
                                                                                  prob=c(probStayPC_2,probIncPC_2/2,probIncPC_2/2), 
                                                                                  replace=TRUE)
  
  # ... Precision Code 2s
  aidData_1[aid.precision.code_1][aidData_1[aid.precision.code_1] == 2] <- sample(size=length(aidData_1[aid.precision.code_1][aidData_1[aid.precision.code_1] == 2]), 
                                                                                  x=c(2,3,4,6), 
                                                                                  prob=c(probStayPC_1,probIncPC_1/3,probIncPC_1/3,probIncPC_1/3), 
                                                                                  replace=TRUE)
  
  aidData_2[aid.precision.code_2][aidData_2[aid.precision.code_2] == 2] <- sample(size=length(aidData_2[aid.precision.code_2][aidData_2[aid.precision.code_2] == 2]), 
                                                                                  x=c(2,3,4,6), 
                                                                                  prob=c(probStayPC_2,probIncPC_2/3,probIncPC_2/3,probIncPC_2/3), 
                                                                                  replace=TRUE)
  
  # ... Precision Code 1s
  if(sim_pc1){
    aidData_1[aid.precision.code_1][aidData_1[aid.precision.code_1] == 1] <- sample(size=length(aidData_1[aid.precision.code_1][aidData_1[aid.precision.code_1] == 1]), 
                                                                                    x=c(1,2,3,4,6), 
                                                                                    prob=c(probStayPC_1,probIncPC_1/4,probIncPC_1/4,probIncPC_1/4,probIncPC_1/4), 
                                                                                    replace=TRUE)
    
    aidData_2[aid.precision.code_2][aidData_2[aid.precision.code_2] == 1] <- sample(size=length(aidData_2[aid.precision.code_2][aidData_2[aid.precision.code_2] == 1]), 
                                                                                    x=c(1,2,3,4,6), 
                                                                                    prob=c(probStayPC_2,probIncPC_2/4,probIncPC_2/4,probIncPC_2/4,probIncPC_2/4), 
                                                                                    replace=TRUE)
    
  }
  
  # Calculating paramSet
  param_set_1 = paramSet(aidData=aidData_1, roiData=roiData, probAidAssume=probAidAssume_1, aid.precision.code=aid.precision.code_1, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_1)
  param_set_2 = paramSet(aidData=aidData_2, roiData=roiData, probAidAssume=probAidAssume_2, aid.precision.code=aid.precision.code_2, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_2)
  
  # Update aid variable
  if(binary){
    roiData[geoSIMEXvariable_1] <- prob_aid(param_set=param_set_1)
    roiData[geoSIMEXvariable_2] <- prob_aid(param_set=param_set_2)
  } else{
    roiData[geoSIMEXvariable_1] <- dollar_expected_value(param_set=param_set_1, dollar_set=as.matrix(aidData_1[aid.amount_1]))
    roiData[geoSIMEXvariable_2] <- dollar_expected_value(param_set=param_set_2, dollar_set=as.matrix(aidData_2[aid.amount_2]))
  }
  
  # Update Model
  df.temp <- model$model
  df.temp[geoSIMEXvariable_1] <- roiData[geoSIMEXvariable_1]
  df.temp[geoSIMEXvariable_2] <- roiData[geoSIMEXvariable_2]
  
  # Update Model
  model.SIMEX <- update(model, data = df.temp)
  
  # Calculate Lambda
  lambda_1 <- calcLambda(param_set_1, maxLambda_denom_1, as.matrix(roiData[roi.area_1]))
  lambda_2 <- calcLambda(param_set_2, maxLambda_denom_2, as.matrix(roiData[roi.area_2]))
  
  # Collecting Results
  
  # Coefficients
  model.SIMEX.coefs <- as.data.frame(as.list(model.SIMEX$coefficients))
  names(model.SIMEX.coefs) <- names(model.SIMEX$coefficients)
  model.SIMEX.coefs$lambda_1 <- lambda_1
  model.SIMEX.coefs$lambda_2 <- lambda_2
  
  # Standard Error
  model.SIMEX.se <- as.data.frame(t(summary(model.SIMEX)$coefficients[,2]))
  names(model.SIMEX.se) <- names(model.SIMEX$coefficients)
  model.SIMEX.se$lambda_1 <- lambda_1
  model.SIMEX.se$lambda_2 <- lambda_2
  
  return(list(model.SIMEX.coefs=model.SIMEX.coefs,
              model.SIMEX.se=model.SIMEX.se))
}


bootIter <- function(i, geoSimulateError.results.df, geoSimulateError.results.df.se, bins, number.from.bin, extrapolation, number.from.bin.average){
  
  geoSimulateError.results.df$rand <- runif(nrow(geoSimulateError.results.df))
  geoSimulateError.results.df.se$rand <- geoSimulateError.results.df$rand
  
  geoSimulateError.results.df <- geoSimulateError.results.df[order(geoSimulateError.results.df$rand),] 
  geoSimulateError.results.df.se <- geoSimulateError.results.df.se[order(geoSimulateError.results.df.se$rand),] 
  
  results.boot <- matrix(NA,nrow=0,ncol=ncol(geoSimulateError.results.df))
  results.boot <- as.data.frame(results.boot)
  
  results.boot.se <- matrix(NA,nrow=0,ncol=ncol(geoSimulateError.results.df.se))
  results.boot.se <- as.data.frame(results.boot.se)
  
  boot.coefs <- as.data.frame(matrix(NA, nrow=1, ncol=ncol(geoSimulateError.results.df)-3))
  names(boot.coefs) <- head(names(geoSimulateError.results.df), -3)
  
  boot.se <- as.data.frame(matrix(NA, nrow=1, ncol=ncol(geoSimulateError.results.df.se)-3))
  names(boot.se) <- head(names(geoSimulateError.results.df.se), -3)
  
  # Randomly Pulling from Results Matrix to Get Subsample for Bootstrapping 
  minLambda <- min(geoSimulateError.results.df$lambda)
  maxLambda <- max(geoSimulateError.results.df$lambda)
  
  binSize <- (maxLambda - minLambda) / bins
  
  binSize_lb <- minLambda
  binSize_ub <- minLambda + binSize
  
  numIter <- bins
  #numIter <- (maxLambda - minLambda) / binSize
  #for(i in 1:ceiling(numIter)){
  for(i in 1:numIter){
    
    results.boot <- rbind(results.boot, geoSimulateError.results.df[(geoSimulateError.results.df$lambda >= binSize_lb) & (geoSimulateError.results.df$lambda <= binSize_ub),][1:number.from.bin,])
    results.boot.se <- rbind(results.boot.se, geoSimulateError.results.df.se[(geoSimulateError.results.df.se$lambda >= binSize_lb) & (geoSimulateError.results.df.se$lambda <= binSize_ub),][1:number.from.bin,])
    
    binSize_lb <- binSize_lb + binSize
    binSize_ub <- binSize_ub + binSize
    
    if(i == (numIter - 1)){
      binSize_ub <-  binSize_ub + 1
    }
  }
    
  # Summarizing by Bins
  results.boot$bin <- rep(1:bins,each=number.from.bin)
  results.boot.se$bin <- rep(1:bins,each=number.from.bin)
  
  # Removing NAs
  results.boot <- results.boot[!is.na(results.boot$rand),]
  results.boot.se <- results.boot.se[!is.na(results.boot.se$rand),]
  
  # Take Means
  if(number.from.bin.average){
    results.boot <- aggregate(. ~ bin, results.boot, mean)
    results.boot.se <- aggregate(. ~ bin, results.boot.se, mean)
  }
  
  results.boot <- subset(results.boot, select = -c(bin))
  results.boot.se <- subset(results.boot.se, select = -c(bin))
  
  # Extrapolating 
  numVars <- ncol(results.boot) - 3
  coefs.geoSIMEX.boot <- matrix(NA, nrow=1, ncol=numVars)
  for(i in 1:numVars){
    
    if(extrapolation=="linear"){
      coefs.geoSIMEX.boot[i] <- summary(lm(as.matrix(results.boot[i]) ~ lambda, data = results.boot))$coefficients[1]
    }
    
    if(extrapolation=="quadratic"){
      coefs.geoSIMEX.boot[i] <- summary(lm(as.matrix(results.boot[i]) ~ lambda + lambda_sq, data = results.boot))$coefficients[1]
    }
    
  }
  coefs.geoSIMEX.boot <- as.data.frame(coefs.geoSIMEX.boot)
  names(coefs.geoSIMEX.boot) <- head(names(results.boot), -3)
  
  numVars <- ncol(results.boot.se) - 3
  se.geoSIMEX.boot <- matrix(NA, nrow=1, ncol=numVars)
  for(i in 1:numVars){
    
    if(extrapolation=="linear"){
      se.geoSIMEX.boot[i] <- summary(lm(as.matrix(results.boot.se[i]) ~ lambda, data = results.boot.se))$coefficients[1]
    }
    
    if(extrapolation=="quadratic"){
      se.geoSIMEX.boot[i] <- summary(lm(as.matrix(results.boot.se[i]) ~ lambda + lambda_sq, data = results.boot.se))$coefficients[1]
    }
    
  }
  se.geoSIMEX.boot <- as.data.frame(se.geoSIMEX.boot)
  names(se.geoSIMEX.boot) <- head(names(results.boot.se), -3)
  
  return(list(coefs.geoSIMEX.boot=coefs.geoSIMEX.boot,
              se.geoSIMEX.boot=se.geoSIMEX.boot))
}

# Model Averaging Change Probability 
model_rand_prob <- function(j, param_set.bin=param_set.bin, aidData=aidData, roiData=roiData, aid.project.amount=aid.project.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary){
  
  # Update aid variable
  if(binary){
    # WILL NEED THE RANDOMLY GENERATING PARAM SET HERE!
  } else{
    roiData[geoSIMEXvariable] <- dollar_expected_value_randProb(param_set.bin=param_set.bin, aid.project.amount=as.matrix(aidData[aid.project.amount]))
  }
  
  # Update Model
  df.temp <- model$model
  df.temp[geoSIMEXvariable] <- roiData[geoSIMEXvariable]
  
  # Update Model
  model.SIMEX <- update(model, data = df.temp)
  
  #### Collecting Results
  
  # Coefficients
  model.SIMEX.coefs <- as.data.frame(as.list(model.SIMEX$coefficients))
  names(model.SIMEX.coefs) <- names(model.SIMEX$coefficients)
  
  # Standard Error
  model.SIMEX.se <- as.data.frame(t(summary(model.SIMEX)$coefficients[,2]))
  names(model.SIMEX.se) <- names(model.SIMEX$coefficients)
  
  # can also return AIC here...
  return(list(model.SIMEX.coefs=model.SIMEX.coefs,
              model.SIMEX.se=model.SIMEX.se))
}

model_rand_prob2 <- function(j, 
                             param_set.bin_1=param_set.bin_1, 
                             param_set.bin_2=param_set.bin_2, 
                             aidData_1=aidData_1, 
                             aidData_2=aidData_2, 
                             roiData=roiData, 
                             aid.amount_1=aid.amount_1, 
                             aid.amount_2=aid.amount_2, 
                             model=model, 
                             geoSIMEXvariable_1=geoSIMEXvariable_1,
                             geoSIMEXvariable_2=geoSIMEXvariable_2, 
                             binary=binary){
  
  # Update aid variable
  if(binary){
    # WILL NEED THE RANDOMLY GENERATING PARAM SET HERE!
  } else{
    roiData[geoSIMEXvariable_1] <- dollar_expected_value_randProb(param_set.bin=param_set.bin_1, dollar_set=as.matrix(aidData_1[aid.amount_1]))
    roiData[geoSIMEXvariable_2] <- dollar_expected_value_randProb(param_set.bin=param_set.bin_2, dollar_set=as.matrix(aidData_2[aid.amount_2]))
  }
  
  #cor.vars <- cor(roiData[geoSIMEXvariable_1],roiData[geoSIMEXvariable_2])
  #while(cor.vars > 0.7){
  #  
  #}
  
  # Update Model
  df.temp <- model$model
  df.temp[geoSIMEXvariable_1] <- roiData[geoSIMEXvariable_1]
  df.temp[geoSIMEXvariable_2] <- roiData[geoSIMEXvariable_2]
  
  # Update Model
  model.SIMEX <- update(model, data = df.temp)
  
  #### Collecting Results
  
  # Coefficients
  model.SIMEX.coefs <- as.data.frame(as.list(model.SIMEX$coefficients))
  names(model.SIMEX.coefs) <- names(model.SIMEX$coefficients)
  
  # Standard Error
  model.SIMEX.se <- as.data.frame(t(summary(model.SIMEX)$coefficients[,2]))
  names(model.SIMEX.se) <- names(model.SIMEX$coefficients)
  
  # can also return AIC here...
  return(list(model.SIMEX.coefs=model.SIMEX.coefs,
              model.SIMEX.se=model.SIMEX.se))
}



dollar_expected_value_randProb <- function(j, param_set.bin, aid.project.amount){
  
  dist.type <- sample(1,x=c(1,2),prob=c(1,1))
  if(dist.type==1){
    probAidGuess.current <- runif(nrow(param_set.bin))
  }
  
  if(dist.type==2){
    probAidGuess.current <- rgamma(nrow(param_set.bin), shape=runif(1,1,20))
  }
  
  param_set.bin <- param_set.bin*probAidGuess.current
  param_set.bin.colsums <- colSums(param_set.bin)
  
  paramDollars <- lapply(1:length(aid.project.amount), function(i) param_set.bin[,i] / param_set.bin.colsums[i])
  paramDollars <- as.data.frame(paramDollars)
  
  paramDollars <- lapply(1:length(aid.project.amount), function(i) paramDollars[,i] * aid.project.amount[i])
  paramDollars <- as.data.frame(paramDollars)
  
  return(rowSums(paramDollars))
}


geoSimulate_realization <- function(j, param_set=param_set, roiData=roiData, aid.project.amount=aid.project.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary, aidData=aidData){
  
  # Update aid variable
  if(binary){
    temp.var <- realization_of_aid(param_set=param_set, aid.project.amount=as.matrix(aidData[aid.project.amount]))
    temp.var[temp.var > 0] <- 1
    roiData[geoSIMEXvariable] <- temp.var
  } else{
    roiData[geoSIMEXvariable] <- realization_of_aid(param_set=param_set, aid.project.amount=as.matrix(aidData[aid.project.amount]))
  }
  
  # Update Model
  df.temp <- model$model
  df.temp[geoSIMEXvariable] <- roiData[geoSIMEXvariable]
  
  # Update Model
  model.SIMEX <- update(model, data = df.temp)
  
  # Collecting Results
  
  # Coefficients
  model.SIMEX.coefs <- as.data.frame(as.list(model.SIMEX$coefficients))
  names(model.SIMEX.coefs) <- names(model.SIMEX$coefficients)
  
  # Standard Error
  model.SIMEX.se <- as.data.frame(t(summary(model.SIMEX)$coefficients[,2]))
  names(model.SIMEX.se) <- names(model.SIMEX$coefficients)
  
  return(list(model.SIMEX.coefs=model.SIMEX.coefs,
              model.SIMEX.se=model.SIMEX.se))
}

realization_of_aid <- function(param_set, aid.project.amount){
  
  param_set_1pjrPerROI <- function(i, param_set=param_set, aid.project.amount=aid.project.amount){
    id <- sample(size=1,x=1:length(param_set[,i]), prob=param_set[,i])
    col.temp <- as.data.frame(matrix(0,nrow=length(param_set[,i]), ncol=1))
    col.temp[id,1] <- 1
    return(col.temp)
  }
  
  param_set_realization <- as.data.frame(lapply(1:ncol(param_set), param_set_1pjrPerROI, param_set=param_set, aid.project.amount=aid.project.amount))
  
  paramDollars <- lapply(1:length(aid.project.amount), function(i) param_set_realization[,i] * aid.project.amount[i])
  paramDollars <- as.data.frame(paramDollars)
  dollars_realization <- rowSums(paramDollars)
  return(dollars_realization)
}

