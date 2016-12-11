##### SECTION 1: GEOSIMEX SIMULATION FUNCTION #####
#roxygen2::roxygenise()

library(parallel)

geoSIMEX_est <- function(model, 
                         geoSIMEXvariable, 
                         roiData, 
                         aidData, 
                         aid.amount,
                         iterations, 
                         bins, 
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
                         mc.cores){
  
  ##### Converting roi.names into IDs #####
  
  #roiData[roi.pc1.name] <- as.numeric(as.factor(roiData[roi.pc1.name]))
  #roi.pc1.name="NAME_2", 
  #roi.pc2.name="NAME_2", 
  #roi.pc3.name="NAME_2", 
  #roi.pc4.name="NAME_1", 
  #roi.pc5.name="NAME_1", 
  #roi.pc6.name="NAME_0", 
  
  ##### Creating Variables
  # Creating probability variable from area
  probAid_area <- as.matrix(roiData[roi.area] / sum(roiData[roi.area]))
  probAid <- as.matrix(roiData[roi.prob.aid] / sum(roiData[roi.prob.aid]))
    
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
  
  ##### FIND PROBABILITY THAT MAXIMIZES LIKELIHOOD HERE #####
    
  ##### Simulating Data with Additional Error and Putting in One Dataframe #####
  prob.increase.list <- runif(iterations,0,1)
  
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
                                       aid.amount=aid.amount,
                                       model=model,
                                       geoSIMEXvariable=geoSIMEXvariable,
                                       binary=binary,
                                       sim_pc1=sim_pc1,
                                       mc.cores=2)
  
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
  geoSimulateError.results.df.se$lambda_sq <- geoSimulateError.results.df.se$lambda^2
  
  bootIter.list <- mclapply(seq(1:numberBootIter), bootIter, 
                            geoSimulateError.results.df=geoSimulateError.results.df, 
                            geoSimulateError.results.df.se=geoSimulateError.results.df.se,
                            bins=bins, 
                            numFromBin=1,
                            mc.cores=2)
  
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
              naive.model = model,
              lambda_naive = lambda_naive,
              simulations = iterations,
              geoSIMEXvariable = geoSIMEXvariable,
              values=geoSimulateError.results.df,
              values.se=geoSimulateError.results.df.se,
              valuesMean=extrapolatedMean.df))
}

# @title Geographic SIMEX
# 
# @author AidData
# 
# @import parallel
# 
# @description
# \code{geoSIMEX} Implementation of the geoSIMEX algorithm for models 
# with spatial uncertainty. Package built to work with data from 
# AidData's data extraction tool.
# 
# @usage
# geoSIMEX(model, geoSIMEXvariable, roiData, aidData, aid.amount, 
# iterations=500, bins=3, fitting.method = "quadratic", roi.area="area",  
# roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", 
# roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc5.name="pc6.id",  
# aid.pc1.centroid.name="centroid.pc1.id", aid.precision.code="precision.code",
# parallel=TRUE, mc.cores=2)
# 
# @param model the naive model
# @param SIMEXvariable character containing the name of the variable with spatial uncertainty
# @param roiData name of dataframe of ROI data 
# @param aidData name of dataframe of aid project data
# @param aid.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
# @param iterations number of simulated error iterations
# @param bins number of bins to group coefficients
# @param fitting.method fitting method for the extrapolation. linear and quadratic are implemented.
# @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 2 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 3 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 4 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 6 and 8 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param parallel parallelize iterations
# @param mc.cores number of cores to use for parallelization
# 
# @details The values contained within roi.pc1.name and aid.pc1.centroid.name variables should be the same.
# 
# @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website]. 
# 
# Just like the lm() and glm() packages, geoSIMEX() is equipped to work with stargazer.
# 
# @references Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. Journal of American Statistical Association, 89, 1314 – 1328
# 
# @examples 
# set.seed(42)
# 
# # Generating Country Dataset
# countryData <- as.data.frame(matrix(NA, nrow=120, ncol=0))
# countryData$pc1.id <- 1:120
# countryData$pc2.id <- rep(1:(120/3), each=3)
# countryData$pc3.id <- rep(1:(120/(3*2)), each=(3*2))
# countryData$pc4.id <- rep(1:(120/(3*2*4)), each=(3*2*4))
# countryData$pc6.id <- 1
# countryData$area <- rgamma(120, shape=2)
# 
# # Creating Aid Dataset Without Error
# aidData <- as.data.frame(matrix(NA,nrow=100,ncol=0))
# aidData$disbursement <- runif(100,0,1) * 10000000 
# aidData$centroid.pc1.id <- sample(size=100,x=c(1:120), prob=rep(1/120, 120), replace=TRUE)
# aidData$precision.code <- 1
# 
# # Adding True Aid to Country Dataset
# countryData$aid <- expected_aid_ROI(aidData=aidData, roiData=countryData, probAidAssume=countryData$area, dollar_set=aidData$disbursement, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
# 
# # Defining True Relation Between Aid and Wealth
# countryData$wealth <- countryData$aid + rnorm(120) * 0.1
# 
# # Creating Datasets with Uncertainty
# countryData <- subset(countryData, select = -c(aid))
# aidData$precision.code <- sample(size=100, x=c(1,2,3,4,6), prob=rep(1/5, 5), replace=TRUE)
# 
# # Calculating Expected Aid and Running Naive Model
# countryData$Expected.Aid <- expected_aid_ROI(aidData=aidData, roiData=countryData, probAidAssume=countryData$area, dollar_set=aidData$disbursement, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
# 
# lm_naive <- lm(wealth ~ Expected.Aid, data=countryData)
# 
# # Implementing GeoSIMEX
# lm_geoSIMEX <- geoSIMEX(model = lm_naive, 
#                         geoSIMEXvariable = "Expected.Aid", 
#                         roiData = countryData, 
#                         aidData = aidData, 
#                         aid.amount = "disbursement")
# 
# summary(lm_geoSIMEX)
# plot(lm_geoSIMEX, variable="Expected.Aid") 

geoSIMEX <- function(x, ...) UseMethod("geoSIMEX")

geoSIMEX.default <- function(model, 
                             geoSIMEXvariable, 
                             roiData, 
                             aidData, 
                             aid.amount,
                             iterations=500, 
                             bins=4, 
                             fitting.method = "quadratic", 
                             roi.area="area",  
                             roi.prob.aid="prob",
                             roi.pc1.name="pc1.id", 
                             roi.pc2.name="pc2.id", 
                             roi.pc3.name="pc3.id", 
                             roi.pc4.name="pc4.id", 
                             roi.pc5.name="pc5.id", 
                             roi.pc6.name="pc6.id",  
                             aid.pc1.centroid.name="centroid.pc1.id", 
                             aid.precision.code="precision.code",
                             binary=FALSE,
                             sim_pc1=TRUE,
                             parallel=TRUE, 
                             mc.cores=2){
  
  est <- geoSIMEX_est(model=model, 
                      geoSIMEXvariable=geoSIMEXvariable, 
                      roiData=roiData, 
                      aidData=aidData, 
                      aid.amount=aid.amount,
                      iterations=iterations, 
                      bins=bins, 
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

plot.geoSIMEX <- function(x, variable, confInt = 95, allSimulations=FALSE, includeTitle=TRUE, name_variable="", ylim="default"){
  
  
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
       pch  = 16)
  
  } else{
    plot(x$valuesMean$lambda, x$valuesMean[,variable],
         xlab = expression((lambda)),
         ylab = paste(variable_name, " coefficient", sep=""),
         main = title,
         xlim = c(0,1),
         ylim = ylim,
         pch  = 16)
  }
  
  if(allSimulations){
    points(x$values$lambda, x$values[,variable],
           pch=".")
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

##### SECTION 1_BUFFER: GEOSIMEX WITH BUFFER #####

geoSIMEX_est <- function(model, 
                         geoSIMEXvariable, 
                         roiData, 
                         aidData, 
                         aid.amount,
                         iterations, 
                         bins, 
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
                         mc.cores){
  
  ##### Converting roi.names into IDs #####
  
  #roiData[roi.pc1.name] <- as.numeric(as.factor(roiData[roi.pc1.name]))
  #roi.pc1.name="NAME_2", 
  #roi.pc2.name="NAME_2", 
  #roi.pc3.name="NAME_2", 
  #roi.pc4.name="NAME_1", 
  #roi.pc5.name="NAME_1", 
  #roi.pc6.name="NAME_0", 
  
  ##### Creating Variables
  # Creating probability variable from area
  probAid_area <- as.matrix(roiData[roi.area] / sum(roiData[roi.area]))
  probAid <- as.matrix(roiData[roi.prob.aid] / sum(roiData[roi.prob.aid]))
  
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
  
  ##### FIND PROBABILITY THAT MAXIMIZES LIKELIHOOD HERE #####
  
  ##### Simulating Data with Additional Error and Putting in One Dataframe #####
  prob.increase.list <- runif(iterations,0,1)
  
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
                                       aid.amount=aid.amount,
                                       model=model,
                                       geoSIMEXvariable=geoSIMEXvariable,
                                       binary=binary,
                                       sim_pc1=sim_pc1,
                                       mc.cores=2)
  
  # ncol(...) is throwing the error. Figure out why? says subscript is out of bounds. ??
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
  geoSimulateError.results.df.se$lambda_sq <- geoSimulateError.results.df.se$lambda^2
  
  bootIter.list <- mclapply(seq(1:numberBootIter), bootIter, 
                            geoSimulateError.results.df=geoSimulateError.results.df, 
                            geoSimulateError.results.df.se=geoSimulateError.results.df.se,
                            bins=bins, 
                            numFromBin=1,
                            mc.cores=2)
  
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
              naive.model = model,
              lambda_naive = lambda_naive,
              simulations = iterations,
              geoSIMEXvariable = geoSIMEXvariable,
              values=geoSimulateError.results.df,
              values.se=geoSimulateError.results.df.se,
              valuesMean=extrapolatedMean.df))
}

# @title Geographic SIMEX
# 
# @author AidData
# 
# @import parallel
# 
# @description
# \code{geoSIMEX} Implementation of the geoSIMEX algorithm for models 
# with spatial uncertainty. Package built to work with data from 
# AidData's data extraction tool.
# 
# @usage
# geoSIMEX(model, geoSIMEXvariable, roiData, aidData, aid.amount, 
# iterations=500, bins=3, fitting.method = "quadratic", roi.area="area",  
# roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", 
# roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc5.name="pc6.id",  
# aid.pc1.centroid.name="centroid.pc1.id", aid.precision.code="precision.code",
# parallel=TRUE, mc.cores=2)
# 
# @param model the naive model
# @param SIMEXvariable character containing the name of the variable with spatial uncertainty
# @param roiData name of dataframe of ROI data 
# @param aidData name of dataframe of aid project data
# @param aid.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
# @param iterations number of simulated error iterations
# @param bins number of bins to group coefficients
# @param fitting.method fitting method for the extrapolation. linear and quadratic are implemented.
# @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 2 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 3 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 4 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 6 and 8 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param parallel parallelize iterations
# @param mc.cores number of cores to use for parallelization
# 
# @details The values contained within roi.pc1.name and aid.pc1.centroid.name variables should be the same.
# 
# @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website]. 
# 
# Just like the lm() and glm() packages, geoSIMEX() is equipped to work with stargazer.
# 
# @references Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. Journal of American Statistical Association, 89, 1314 – 1328
# 
# @examples 
# set.seed(42)
# 
# # Generating Country Dataset
# countryData <- as.data.frame(matrix(NA, nrow=120, ncol=0))
# countryData$pc1.id <- 1:120
# countryData$pc2.id <- rep(1:(120/3), each=3)
# countryData$pc3.id <- rep(1:(120/(3*2)), each=(3*2))
# countryData$pc4.id <- rep(1:(120/(3*2*4)), each=(3*2*4))
# countryData$pc6.id <- 1
# countryData$area <- rgamma(120, shape=2)
# 
# # Creating Aid Dataset Without Error
# aidData <- as.data.frame(matrix(NA,nrow=100,ncol=0))
# aidData$disbursement <- runif(100,0,1) * 10000000 
# aidData$centroid.pc1.id <- sample(size=100,x=c(1:120), prob=rep(1/120, 120), replace=TRUE)
# aidData$precision.code <- 1
# 
# # Adding True Aid to Country Dataset
# countryData$aid <- expected_aid_ROI(aidData=aidData, roiData=countryData, probAidAssume=countryData$area, dollar_set=aidData$disbursement, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
# 
# # Defining True Relation Between Aid and Wealth
# countryData$wealth <- countryData$aid + rnorm(120) * 0.1
# 
# # Creating Datasets with Uncertainty
# countryData <- subset(countryData, select = -c(aid))
# aidData$precision.code <- sample(size=100, x=c(1,2,3,4,6), prob=rep(1/5, 5), replace=TRUE)
# 
# # Calculating Expected Aid and Running Naive Model
# countryData$Expected.Aid <- expected_aid_ROI(aidData=aidData, roiData=countryData, probAidAssume=countryData$area, dollar_set=aidData$disbursement, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
# 
# lm_naive <- lm(wealth ~ Expected.Aid, data=countryData)
# 
# # Implementing GeoSIMEX
# lm_geoSIMEX <- geoSIMEX(model = lm_naive, 
#                         geoSIMEXvariable = "Expected.Aid", 
#                         roiData = countryData, 
#                         aidData = aidData, 
#                         aid.amount = "disbursement")
# 
# summary(lm_geoSIMEX)
# plot(lm_geoSIMEX, variable="Expected.Aid") 

geoSIMEX <- function(x, ...) UseMethod("geoSIMEX")

geoSIMEX.default <- function(model, 
                             geoSIMEXvariable, 
                             roiData, 
                             aidData, 
                             aid.amount,
                             iterations=500, 
                             bins=4, 
                             fitting.method = "quadratic", 
                             roi.area="area",  
                             roi.prob.aid="prob",
                             roi.pc1.name="pc1.id", 
                             roi.pc2.name="pc2.id", 
                             roi.pc3.name="pc3.id", 
                             roi.pc4.name="pc4.id", 
                             roi.pc5.name="pc5.id", 
                             roi.pc6.name="pc6.id",  
                             aid.pc1.centroid.name="centroid.pc1.id", 
                             aid.precision.code="precision.code",
                             binary=FALSE,
                             sim_pc1=TRUE,
                             parallel=TRUE, 
                             mc.cores=2){
  
  est <- geoSIMEX_est(model=model, 
                      geoSIMEXvariable=geoSIMEXvariable, 
                      roiData=roiData, 
                      aidData=aidData, 
                      aid.amount=aid.amount,
                      iterations=iterations, 
                      bins=bins, 
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

plot.geoSIMEX <- function(x, variable, confInt = 95, allSimulations=FALSE, includeTitle=TRUE, name_variable=""){
  
  
  #name_variable = variable
  # NOTE: Can make things like "all simulations", "confidence bands", "% interval"
  # all parameters / options.
  
  if(name_variable == ""){
    variable_name =  variable
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
  plot(x$valuesMean$lambda, x$valuesMean[,variable],
       xlab = expression((lambda)),
       ylab = paste(variable_name, " coefficient", sep=""),
       main = title,
       xlim = c(0,1),
       ylim = c(min(var.vales.all), max(var.vales.all)),
       pch  = 16)
  
  if(allSimulations){
    points(x$values$lambda, x$values[,variable],
           pch=".")
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

##### SECTION 1_1: GEOSIMEX SIMULATION FUNCTION, 2 VARS #####

geoSIMEX2_est <- function(model, 
                          geoSIMEXvariable_1, 
                          geoSIMEXvariable_2,
                          roiData, 
                          aidData_1,
                          aidData_2,
                          aid.amount_1,
                          aid.amount_2,
                          iterations, 
                          bins, 
                          fitting.method, 
                          roi.area_1,
                          roi.area_2, 
                          roi.pc1.name, 
                          roi.pc2.name, 
                          roi.pc3.name, 
                          roi.pc4.name, 
                          roi.pc5.name, 
                          roi.pc6.name,  
                          aid.pc1.centroid.name_1,
                          aid.pc1.centroid.name_2, 
                          aid.precision.code_1,
                          aid.precision.code_2,
                          binary,
                          sim_pc1,
                          average_lambda,
                          diagonal,
                          parallel, 
                          mc.cores){
  
  #average_lambda <- TRUE
  ##### Creating Variables
  # Creating probability variable from area
  probAid_area_1 <- as.matrix(roiData[roi.area_1] / sum(roiData[roi.area_1]))
  probAid_area_2 <- as.matrix(roiData[roi.area_2] / sum(roiData[roi.area_2]))
  
  # Creating matrix of original precision codes
  precision.code.original_1 <- as.matrix(aidData_1[aid.precision.code_1])
  precision.code.original_2 <- as.matrix(aidData_2[aid.precision.code_2])
  
  ##### Calculating Maximum Lambda Denominator #####
  aidData_1[aid.precision.code_1] <- 6  
  param_set = paramSet(aidData=aidData_1, roiData=roiData, probAidAssume=probAid_area_1, aid.precision.code=aid.precision.code_1, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_1)  
  param_set[param_set!=0] <- 1
  param_set <- param_set * as.matrix(roiData[roi.area_1])
  maxLambda_denom_1 <- sum(colSums(param_set))
  
  aidData_2[aid.precision.code_2] <- 6  
  param_set = paramSet(aidData=aidData_2, roiData=roiData, probAidAssume=probAid_area_2, aid.precision.code=aid.precision.code_2, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_2)  
  param_set[param_set!=0] <- 1
  param_set <- param_set * as.matrix(roiData[roi.area_2])
  maxLambda_denom_2 <- sum(colSums(param_set))
  
  ##### Calculate Naive Lambda #####
  aidData_1[aid.precision.code_1] <- precision.code.original_1 
  param_set = paramSet(aidData=aidData_1, roiData=roiData, probAidAssume=probAid_area_1, aid.precision.code=aid.precision.code_1, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_1)  
  lambda_naive_1 <- calcLambda(param_set, maxLambda_denom_1, as.matrix(roiData[roi.area_1]))
  
  aidData_2[aid.precision.code_2] <- precision.code.original_2 
  param_set = paramSet(aidData=aidData_2, roiData=roiData, probAidAssume=probAid_area_2, aid.precision.code=aid.precision.code_2, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_2)  
  lambda_naive_2 <- calcLambda(param_set, maxLambda_denom_2, as.matrix(roiData[roi.area_2]))
  
  ##### Simulating Data with Additional Error and Putting in One Dataframe #####
  #prob.increase.list <- runif(iterations,0,1)
  #prob.increase.list_2 <- runif(iterations,0,1)
  #prob.increase.list_2 <- prob.increase.list_1
  
  geoSimulateError.results <- mclapply(1:iterations, geoSimulateError2, 
                                       aidData_1=aidData_1, 
                                       aidData_2=aidData_2, 
                                       roiData=roiData, 
                                       probAidAssume_1=probAid_area_1, 
                                       probAidAssume_2=probAid_area_2,
                                       PC_researcherSees_1=precision.code.original_1, 
                                       PC_researcherSees_2=precision.code.original_2, 
                                       maxLambda_denom_1=maxLambda_denom_1,
                                       maxLambda_denom_2=maxLambda_denom_2,
                                       roi.area_1=roi.area_1,
                                       roi.area_2=roi.area_2,
                                       aid.precision.code_1=aid.precision.code_1, 
                                       aid.precision.code_2=aid.precision.code_2, 
                                       roi.pc1.name=roi.pc1.name, 
                                       roi.pc2.name=roi.pc2.name, 
                                       roi.pc3.name=roi.pc3.name, 
                                       roi.pc4.name=roi.pc4.name, 
                                       roi.pc5.name=roi.pc5.name, 
                                       roi.pc6.name=roi.pc6.name, 
                                       aid.pc1.centroid.name_1=aid.pc1.centroid.name_1,
                                       aid.pc1.centroid.name_2=aid.pc1.centroid.name_2,
                                       aid.amount_1=aid.amount_1,
                                       aid.amount_2=aid.amount_2,
                                       model=model,
                                       geoSIMEXvariable_1=geoSIMEXvariable_1,
                                       geoSIMEXvariable_2=geoSIMEXvariable_2,
                                       binary=binary,
                                       diagonal=diagonal,
                                       sim_pc1=sim_pc1,
                                       mc.cores=2)
  
  geoSimulateError.results.df <- matrix(NA, ncol=ncol(geoSimulateError.results[[1]][["model.SIMEX.coefs"]]),nrow=0)
  geoSimulateError.results.df <- as.data.frame(geoSimulateError.results.df)
  
  geoSimulateError.results.df.se <- matrix(NA, ncol=ncol(geoSimulateError.results[[1]][["model.SIMEX.se"]]),nrow=0)
  geoSimulateError.results.df.se <- as.data.frame(geoSimulateError.results.df.se)
  
  for(i in 1:iterations){
    geoSimulateError.results.df <- rbind(geoSimulateError.results.df, geoSimulateError.results[[i]][["model.SIMEX.coefs"]])
    geoSimulateError.results.df.se <- rbind(geoSimulateError.results.df.se, geoSimulateError.results[[i]][["model.SIMEX.se"]])
  }
  
  if(average_lambda){
    avg.lambda <- (geoSimulateError.results.df$lambda_1 + geoSimulateError.results.df$lambda_2)/2
    geoSimulateError.results.df$lambda_1 <- avg.lambda
    geoSimulateError.results.df$lambda_2  <- avg.lambda
  }
  
  # Adding Bin Number into Dataset
  bin_1_size <- (max(geoSimulateError.results.df$lambda_1) - min(geoSimulateError.results.df$lambda_1)) / bins
  for(b in 1:bins){
    geoSimulateError.results.df$bin_1[(geoSimulateError.results.df$lambda_1 >= (min(geoSimulateError.results.df$lambda_1) + bin_1_size*(b-1))) & 
                                        (geoSimulateError.results.df$lambda_1 <= (min(geoSimulateError.results.df$lambda_1) + bin_1_size*b))] <- b
  }
  
  bin_2_size <- (max(geoSimulateError.results.df$lambda_2) - min(geoSimulateError.results.df$lambda_2)) / bins
  for(b in 1:bins){
    geoSimulateError.results.df$bin_2[(geoSimulateError.results.df$lambda_2 >= (min(geoSimulateError.results.df$lambda_2) + bin_2_size*(b-1))) & 
                                        (geoSimulateError.results.df$lambda_2 <= (min(geoSimulateError.results.df$lambda_2) + bin_2_size*b))] <- b
  }
  
  geoSimulateError.results.df$bin_12 <- paste(as.character(geoSimulateError.results.df$bin_1),as.character(geoSimulateError.results.df$bin_2),sep="")
  geoSimulateError.results.df$bin_12 <- as.numeric(geoSimulateError.results.df$bin_12)
  
  geoSimulateError.results.df.se$bin_1 <- geoSimulateError.results.df$bin_1
  geoSimulateError.results.df.se$bin_2 <- geoSimulateError.results.df$bin_2
  geoSimulateError.results.df.se$bin_12 <- geoSimulateError.results.df$bin_12
  
  # Remove NAs
  geoSimulateError.results.df <- geoSimulateError.results.df[!is.na(geoSimulateError.results.df$bin_12),]
  geoSimulateError.results.df.se <- geoSimulateError.results.df.se[!is.na(geoSimulateError.results.df.se$bin_12),]
  
  # Only look at lambda along diagonal
  if(diagonal){
    geoSimulateError.results.df <- geoSimulateError.results.df[geoSimulateError.results.df$bin_1 == geoSimulateError.results.df$bin_2,]
    geoSimulateError.results.df.se <- geoSimulateError.results.df.se[geoSimulateError.results.df.se$bin_1 == geoSimulateError.results.df.se$bin_2,]
  }
  
  # Mean Bin Coefficient
  extrapolatedMean.df <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(geoSimulateError.results.df)))
  for(b in names(table(geoSimulateError.results.df$bin_12))){
    df.temp <- geoSimulateError.results.df[geoSimulateError.results.df$bin_12 == b,]
    df.mean.temp <- colMeans(df.temp)
    extrapolatedMean.df <- rbind(extrapolatedMean.df, df.mean.temp)
  }
  names(extrapolatedMean.df) <- names(df.temp)
  
  average_lambda
  # Mean Bin Coefficient Extrapolation
  extrapolatedMean.df$lambda_1_sq <- extrapolatedMean.df$lambda_1^2
  extrapolatedMean.df$lambda_2_sq <- extrapolatedMean.df$lambda_2^2
  
  if(average_lambda){
    numVars <- ncol(extrapolatedMean.df) - 7
    coef.geoSIMEX <- matrix(NA, nrow=1, ncol=numVars)
    for(i in 1:numVars){
      coef.geoSIMEX[i] <- summary(lm(as.matrix(extrapolatedMean.df[i]) ~ lambda_1 + lambda_1_sq, data = extrapolatedMean.df))$coefficients[1]
    }
    coef.geoSIMEX <- as.data.frame(coef.geoSIMEX)
    names(coef.geoSIMEX) <- head(names(extrapolatedMean.df), -7)
  } else{
    numVars <- ncol(extrapolatedMean.df) - 7
    coef.geoSIMEX <- matrix(NA, nrow=1, ncol=numVars)
    for(i in 1:numVars){
      coef.geoSIMEX[i] <- summary(lm(as.matrix(extrapolatedMean.df[i]) ~ lambda_1 + lambda_2 + lambda_1_sq + lambda_2_sq, data = extrapolatedMean.df))$coefficients[1]
    }
    coef.geoSIMEX <- as.data.frame(coef.geoSIMEX)
    names(coef.geoSIMEX) <- head(names(extrapolatedMean.df), -7)
  }
  
  ##### Bootstrap Standard Error #####
  numberBootIter <- nrow(geoSimulateError.results.df)
  
  geoSimulateError.results.df$lambda_1_sq <- geoSimulateError.results.df$lambda_1^2
  geoSimulateError.results.df$lambda_2_sq <- geoSimulateError.results.df$lambda_2^2
  geoSimulateError.results.df.se$lambda_1_sq <- geoSimulateError.results.df.se$lambda_1^2
  geoSimulateError.results.df.se$lambda_2_sq <- geoSimulateError.results.df.se$lambda_2^2
  
  bootIter.list <- mclapply(seq(1:numberBootIter), bootIter2, 
                            geoSimulateError.results.df=geoSimulateError.results.df, 
                            geoSimulateError.results.df.se=geoSimulateError.results.df.se,
                            average_lambda,
                            mc.cores=2)
  
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
  se <- se.geoSIMEX
  coef.se <- cbind(coef,t(se))
  df = model$df
  
  # Values so will work with stargazer
  # https://github.com/cran/stargazer/blob/master/R/stargazer-internal.R
  # NOTE: Fix so will reflect geoSIMEX coefficients and not naive coefficients?
  #residuals = model$residuals 
  residuals <- rep(NA, length(model$residuals))
  
  model.imprecision.variance <- rbind(model.var,var.imprecision)
  row.names(model.imprecision.variance) <- c("Model Variance", "Imprecision Variance")
  model.imprecision.variance <- as.data.frame(model.imprecision.variance)
  
  lambda_naive <- (lambda_naive_1 + lambda_naive_2) / 2
  
  return(list(coefficients=coef.geoSIMEX,
              StdErr=se.geoSIMEX,
              coef.se = coef.se,
              df = df,
              variance.model.imprecision = model.imprecision.variance,
              residuals = residuals,
              naive.model = model,
              lambda_1 = lambda_naive_1,
              lambda_2 = lambda_naive_2,
              lambda = lambda_naive,
              simulations = iterations,
              geoSIMEXvariable_1 = geoSIMEXvariable_1,
              geoSIMEXvariable_2 = geoSIMEXvariable_2,
              values=geoSimulateError.results.df,
              values.se=geoSimulateError.results.df.se,
              valuesMean=extrapolatedMean.df))
}

# @title Geographic SIMEX
# 
# @author AidData
# 
# @import parallel
# 
# @description
# \code{geoSIMEX} Implementation of the geoSIMEX algorithm for models 
# with spatial uncertainty. Package built to work with data from 
# AidData's data extraction tool.
# 
# @usage
# geoSIMEX(model, geoSIMEXvariable, roiData, aidData, aid.amount, 
# iterations=500, bins=3, fitting.method = "quadratic", roi.area="area",  
# roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", 
# roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc5.name="pc6.id",  
# aid.pc1.centroid.name="centroid.pc1.id", aid.precision.code="precision.code",
# parallel=TRUE, mc.cores=2)
# 
# @param model the naive model
# @param SIMEXvariable character containing the name of the variable with spatial uncertainty
# @param roiData name of dataframe of ROI data 
# @param aidData name of dataframe of aid project data
# @param aid.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
# @param iterations number of simulated error iterations
# @param bins number of bins to group coefficients
# @param fitting.method fitting method for the extrapolation. linear and quadratic are implemented.
# @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 2 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 3 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 4 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 6 and 8 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param parallel parallelize iterations
# @param mc.cores number of cores to use for parallelization
# 
# @details The values contained within roi.pc1.name and aid.pc1.centroid.name variables should be the same.
# 
# @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website]. 
# 
# Just like the lm() and glm() packages, geoSIMEX() is equipped to work with stargazer.
# 
# @references Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. Journal of American Statistical Association, 89, 1314 – 1328
# 
# @examples 
# set.seed(42)
# 
# # Generating Country Dataset
# countryData <- as.data.frame(matrix(NA, nrow=120, ncol=0))
# countryData$pc1.id <- 1:120
# countryData$pc2.id <- rep(1:(120/3), each=3)
# countryData$pc3.id <- rep(1:(120/(3*2)), each=(3*2))
# countryData$pc4.id <- rep(1:(120/(3*2*4)), each=(3*2*4))
# countryData$pc6.id <- 1
# countryData$area <- rgamma(120, shape=2)
# 
# # Creating Aid Dataset Without Error
# aidData <- as.data.frame(matrix(NA,nrow=100,ncol=0))
# aidData$disbursement <- runif(100,0,1) * 10000000 
# aidData$centroid.pc1.id <- sample(size=100,x=c(1:120), prob=rep(1/120, 120), replace=TRUE)
# aidData$precision.code <- 1
# 
# # Adding True Aid to Country Dataset
# countryData$aid <- expected_aid_ROI(aidData=aidData, roiData=countryData, probAidAssume=countryData$area, dollar_set=aidData$disbursement, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
# 
# # Defining True Relation Between Aid and Wealth
# countryData$wealth <- countryData$aid + rnorm(120) * 0.1
# 
# # Creating Datasets with Uncertainty
# countryData <- subset(countryData, select = -c(aid))
# aidData$precision.code <- sample(size=100, x=c(1,2,3,4,6), prob=rep(1/5, 5), replace=TRUE)
# 
# # Calculating Expected Aid and Running Naive Model
# countryData$Expected.Aid <- expected_aid_ROI(aidData=aidData, roiData=countryData, probAidAssume=countryData$area, dollar_set=aidData$disbursement, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
# 
# lm_naive <- lm(wealth ~ Expected.Aid, data=countryData)
# 
# # Implementing GeoSIMEX
# lm_geoSIMEX <- geoSIMEX(model = lm_naive, 
#                         geoSIMEXvariable = "Expected.Aid", 
#                         roiData = countryData, 
#                         aidData = aidData, 
#                         aid.amount = "disbursement")
# 
# summary(lm_geoSIMEX)
# plot(lm_geoSIMEX, variable="Expected.Aid") 

geoSIMEX2 <- function(x, ...) UseMethod("geoSIMEX2")

geoSIMEX2.default <- function(model, 
                              geoSIMEXvariable_1, 
                              geoSIMEXvariable_2, 
                              roiData, 
                              aidData_1, 
                              aidData_2,
                              aid.amount_1,
                              aid.amount_2,
                              iterations=500, 
                              bins=4, 
                              fitting.method = "quadratic", 
                              roi.area_1="area",  
                              roi.area_2="area",  
                              roi.pc1.name="pc1.id", 
                              roi.pc2.name="pc2.id", 
                              roi.pc3.name="pc3.id", 
                              roi.pc4.name="pc4.id", 
                              roi.pc5.name="pc5.id", 
                              roi.pc6.name="pc6.id",  
                              aid.pc1.centroid.name_1="centroid.pc1.id_1", 
                              aid.pc1.centroid.name_2="centroid.pc1.id_2", 
                              aid.precision.code_1="precision.code_1",
                              aid.precision.code_2="precision.code_2",
                              binary=FALSE,
                              sim_pc1=TRUE,
                              average_lambda=TRUE,
                              diagonal,
                              parallel=TRUE, 
                              mc.cores=2){
  
  est <- geoSIMEX2_est(model=model, 
                       geoSIMEXvariable_1=geoSIMEXvariable_1, 
                       geoSIMEXvariable_2=geoSIMEXvariable_2, 
                       roiData=roiData, 
                       aidData_1=aidData_1, 
                       aidData_2=aidData_2, 
                       aid.amount_1=aid.amount_1,
                       aid.amount_2=aid.amount_2,
                       iterations=iterations, 
                       bins=bins, 
                       fitting.method = fitting.method, 
                       roi.area_1=roi.area_1,  
                       roi.area_2=roi.area_2,  
                       roi.pc1.name=roi.pc1.name, 
                       roi.pc2.name=roi.pc2.name, 
                       roi.pc3.name=roi.pc3.name, 
                       roi.pc4.name=roi.pc4.name, 
                       roi.pc5.name=roi.pc5.name, 
                       roi.pc6.name=roi.pc6.name,  
                       aid.pc1.centroid.name_1=aid.pc1.centroid.name_1, 
                       aid.pc1.centroid.name_2=aid.pc1.centroid.name_2, 
                       aid.precision.code_1=aid.precision.code_1,
                       aid.precision.code_2=aid.precision.code_2,
                       binary=binary,
                       sim_pc1=sim_pc1,
                       average_lambda=average_lambda,
                       diagonal=diagonal,
                       parallel=parallel, 
                       mc.cores=mc.cores)
  
  est$call <- model$call
  
  class(est) <- "geoSIMEX2"
  est
}

print.geoSIMEX2 <- function(x, ...){
  cat("Naive Model:\n")
  print(x$call)
  cat("\ngeoSIMEX-Variables: ", x$geoSIMEXvariable_1, " ", x$geoSIMEXvariable_2, sep="")
  cat("\nNumer of Simulations: ", x$simulations, sep= "")
  cat("\n\nCoefficients:\n")
  print(x$coefficients, row.names=FALSE)  
}

summary.geoSIMEX2 <- function(object, ...){
  
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
  
  class(res) <- "summary.geoSIMEX2"
  res
}

print.summary.geoSIMEX2 <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

plot.geoSIMEX2 <- function(x, variable, confInt = 95, allSimulations=FALSE, includeTitle=TRUE, name_variable=""){
  
  variable_name = variable
  #name_variable = variable
  # NOTE: Can make things like "all simulations", "confidence bands", "% interval"
  # all parameters / options.
  
  if(name_variable == ""){
    variable_name = variable_name
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
  CI.naive$lambda <- x$lambda
  CI.naive$value[1] <- coef.naive - t_critical*se.naive
  CI.naive$value[2] <- coef.naive + t_critical*se.naive
  
  # All Variable Values to Define Min / Max ylim of Plot
  var.vales.all <- c(x$valuesMean[,variable], x$values[,variable], x$coefficients[,variable], CI.geoSIMEX$value[1], CI.geoSIMEX$value[2], CI.naive$value[1], CI.naive$value[2])
  
  title = ""
  
  if(includeTitle){
    
    title = paste("geoSIMEX Plot of ",variable_name, sep="")
    
  }
  
  x$valuesMean$lambda <- (x$valuesMean$lambda_1 + x$valuesMean$lambda_2) / 2
  # Plot Mean Lambda Values
  plot(x$valuesMean$lambda, x$valuesMean[,variable],
       xlab = expression((lambda)),
       ylab = paste(variable_name, " coefficient", sep=""),
       main = title,
       xlim = c(0,1),
       ylim = c(min(var.vales.all), max(var.vales.all)),
       pch  = 16)
  
  if(allSimulations){
    points(x$values$lambda, x$values[,variable],
           pch=".")
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
  points(x$lambda, coef.naive,
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
                                     aid.amount, 
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
  
  param_set <- paramSet(aidData=aidData, roiData=roiData, probAidAssume=runif(nrow(roiData)), aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set.bin <- (param_set > 0)*1
  
  model.list <- mclapply(1:iterations, model_rand_prob, param_set.bin=param_set.bin, aidData=aidData, roiData=roiData, aid.amount=aid.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary)
  
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
                                         aid.amount, 
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
                                         mc.cores=2){
  
  est <- modelAverageRandProb_est(iterations=iterations, 
                                  model=model,
                                  geoSIMEXvariable=geoSIMEXvariable,
                                  roiData=roiData, 
                                  aidData=aidData,
                                  aid.amount=aid.amount, 
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

##### SECTION 1A_1: MODEL AVERAGING - RAND PROB, 2 covs #####

modelAverageRandProb2_est <- function(iterations, 
                                      model, 
                                      geoSIMEXvariable_1, 
                                      geoSIMEXvariable_2, 
                                      roiData, 
                                      aidData_1, 
                                      aidData_2, 
                                      aid.amount_1, 
                                      aid.amount_2, 
                                      roi.pc1.name, 
                                      roi.pc2.name, 
                                      roi.pc3.name, 
                                      roi.pc4.name, 
                                      roi.pc5.name, 
                                      roi.pc6.name, 
                                      aid.pc1.centroid.name_1,
                                      aid.pc1.centroid.name_2,
                                      aid.precision.code_1,
                                      aid.precision.code_2,
                                      binary, 
                                      parallel, 
                                      mc.cores){
  
  param_set_1 <- paramSet(aidData=aidData_1, roiData=roiData, probAidAssume=runif(nrow(roiData)), aid.precision.code=aid.precision.code_1, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_1)
  param_set_2 <- paramSet(aidData=aidData_2, roiData=roiData, probAidAssume=runif(nrow(roiData)), aid.precision.code=aid.precision.code_2, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_2)
  param_set.bin_1 <- (param_set_1 > 0)*1
  param_set.bin_2 <- (param_set_2 > 0)*1
  
  # names(colSums(param_set_1)[colSums(param_set_1)==0])
  
  model.list <- mclapply(1:iterations, model_rand_prob2, 
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
                         binary=binary)
  
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
              geoSIMEXvariable_1=geoSIMEXvariable_1,
              geoSIMEXvariable_2=geoSIMEXvariable_2,
              iterations=iterations,
              residuals = model$residuals))
}

modelAverageRandProb2 <- function(x, ...) UseMethod("modelAverageRandProb2")

modelAverageRandProb2.default <- function(iterations=500, 
                                          model,
                                          geoSIMEXvariable_1,
                                          geoSIMEXvariable_2,
                                          roiData, 
                                          aidData_1,
                                          aidData_2,
                                          aid.amount_1, 
                                          aid.amount_2, 
                                          roi.pc1.name="pc1.id", 
                                          roi.pc2.name="pc2.id", 
                                          roi.pc3.name="pc3.id", 
                                          roi.pc4.name="pc4.id", 
                                          roi.pc5.name="pc5.id", 
                                          roi.pc6.name="pc6.id", 
                                          aid.pc1.centroid.name_1="centroid.pc1.id",
                                          aid.pc1.centroid.name_2="centroid.pc1.id",
                                          aid.precision.code_1="precision.code",
                                          aid.precision.code_2="precision.code",
                                          binary=FALSE,
                                          parallel=TRUE, 
                                          mc.cores=2){
  
  est <- modelAverageRandProb2_est(iterations=iterations, 
                                   model=model,
                                   geoSIMEXvariable_1=geoSIMEXvariable_1,
                                   geoSIMEXvariable_2=geoSIMEXvariable_2,
                                   roiData=roiData, 
                                   aidData_1=aidData_1,
                                   aidData_2=aidData_2,
                                   aid.amount_1=aid.amount_1, 
                                   aid.amount_2=aid.amount_2, 
                                   roi.pc1.name=roi.pc1.name, 
                                   roi.pc2.name=roi.pc2.name, 
                                   roi.pc3.name=roi.pc3.name, 
                                   roi.pc4.name=roi.pc4.name, 
                                   roi.pc5.name=roi.pc5.name, 
                                   roi.pc6.name=roi.pc6.name, 
                                   aid.pc1.centroid.name_1=aid.pc1.centroid.name_1,
                                   aid.pc1.centroid.name_2=aid.pc1.centroid.name_2,
                                   aid.precision.code_1=aid.precision.code_1,
                                   aid.precision.code_2=aid.precision.code_2,
                                   binary=binary,
                                   parallel=parallel, 
                                   mc.cores=mc.cores)
  
  est$call <- model$call
  
  class(est) <- "modelAverageRandProb2"
  est
}

print.modelAverageRandProb2 <- function(x, ...){
  cat("Naive Model:\n")
  print(x$call)
  cat("\ngeoSIMEX-Variables: ", x$geoSIMEXvariable_1,", ",x$geoSIMEXvariable_2, sep="")
  cat("\nNumer of Iterations: ", x$iterations, sep= "")
  cat("\n\nCoefficients:\n")
  print(x$coefficients, row.names=FALSE)  
}

summary.modelAverageRandProb2 <- function(object, ...){
  
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
  
  class(res) <- "summary.modelAverageRandProb2"
  res
}

print.summary.modelAverageRandProb2 <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}



##### SECTION 1B: MODEL AVERAGING BAYES - RAND PROB #####

modelAverageMaxLikeRandProb_est <- function(iterations, 
                                            model, 
                                            geoSIMEXvariable, 
                                            roiData, 
                                            aidData, 
                                            aid.amount, 
                                            roi.pc1.name, 
                                            roi.pc2.name, 
                                            roi.pc3.name, 
                                            roi.pc4.name, 
                                            roi.pc5.name, 
                                            roi.pc6.name, 
                                            aid.pc1.centroid.name,
                                            aid.precision.code,
                                            binary,
                                            threshold, 
                                            burn.in, 
                                            burn.in.check, 
                                            target.accept,
                                            proposal.width, 
                                            do.burn.in,
                                            parallel, 
                                            mc.cores){
  
  dollar_expected_value_fast <- function(param_set.bin, dollar_set, probAidGuess.current){
    
    param_set.bin <- param_set.bin*probAidGuess.current
    param_set.bin.colsums <- colSums(param_set.bin)
    
    paramDollars <- lapply(1:length(dollar_set), function(i) param_set.bin[,i] / param_set.bin.colsums[i])
    paramDollars <- as.data.frame(paramDollars)
    
    paramDollars <- lapply(1:length(dollar_set), function(i) paramDollars[,i] * dollar_set[i])
    paramDollars <- as.data.frame(paramDollars)
    
    return( rowSums(paramDollars) )
  }
  
  bayes_model_avg <- function(k, model.initial, probAidGuess.current, subcountyData, aidData, param_set.bin=param_set.bin, chain.length=chain.length, threshold=threshold, burn.in=burn.in, burn.in.check=burn.in.check, target.accept=target.accept, proposal.width=proposal.width, do.burn.in=do.burn.in){
    
    temp=k
    
    accept.vec <- matrix(NA, nrow=0, ncol=1)
    
    roiData[geoSIMEXvariable] <- dollar_expected_value_fast(param_set.bin=param_set.bin, dollar_set=aidData[,aid.amount], probAidGuess.current=probaid.current)
    df.temp <- model.initial$model
    df.temp[geoSIMEXvariable] <- roiData[geoSIMEXvariable]
    model <- update(model.initial, data = df.temp)
    
    likelihood.current <- logLik(model)
    
    probAidGuess.best <- probAidGuess.current
    aic.best <- extractAIC(model)[2]
    
    coefs <- summary(model)$coefficients[,1]
    se <- summary(model)$coefficients[,2]
    coefs <- as.data.frame(t(coefs))
    se <- as.data.frame(t(se))
    AIC <- as.data.frame(as.matrix(extractAIC(model)[2], nrow=1, ncol=1))
    names(AIC) <- "AIC"
    results.df.coef <- cbind(coefs, AIC)
    results.df.se <- cbind(se, AIC)
    
    for(i in 1:iterations){
      
      # Also to Try:
      # series of if statements for how update... start with large groups, then work down.
      
      # Proposed probabilities and model      
      probAidGuess.proposed <- probAidGuess.current + rnorm(length(probAidGuess.current),mean=0,sd=proposal.width)
      probAidGuess.proposed[probAidGuess.proposed < 0] <- runif(length(probAidGuess.proposed[probAidGuess.proposed < 0]))*proposal.width
      roiData[geoSIMEXvariable] <- dollar_expected_value_fast(param_set.bin=param_set.bin, dollar_set=aidData[,aid.amount], probAidGuess.current=probAidGuess.proposed)
      df.temp <- model.initial$model
      df.temp[geoSIMEXvariable] <- roiData[geoSIMEXvariable]
      model <- update(model.initial, data = df.temp)
      likelihood.proposed <- logLik(model)
      
      # Update 
      c <- likelihood.proposed - likelihood.current
      accept <- threshold <= exp(c)
      accept.vec <- rbind(accept.vec, accept)
      
      if(accept){
        
        likelihood.current <- likelihood.proposed
        probAidGuess.current <- probAidGuess.proposed
        coefs <- summary(model)$coefficients[,1]
        se <- summary(model)$coefficients[,2]
        coefs <- as.data.frame(t(coefs))
        se <- as.data.frame(t(se))
        AIC <- as.data.frame(as.matrix(extractAIC(model)[2], nrow=1, ncol=1))
        names(AIC) <- "AIC"
        results.df.coef.temp <- cbind(coefs, AIC)
        results.df.se.temp <- cbind(se, AIC)
        
        results.df.coef <- rbind(results.df.coef, results.df.coef.temp)
        results.df.se <- rbind(results.df.se, results.df.se.temp)
        
        # TAKING LAST VALUE OF CHAIN
        #probAidGuess.best <- probAidGuess.proposed
        
      }
      
      # UPDATE PROB AID GUESS
      # If has lower AIC than best, update aic.best and probAidGuess.best
      if(extractAIC(model)[2] < aic.best){
      #if(abs(extractAIC(model)[2]) < abs(aic.best)){  
        aic.best <- AIC(model)
        probAidGuess.best <- probAidGuess.proposed
      }
      
      # burn-in: adjust threshold for accepting
      if((i >= burn.in.check+1) & (i < burn.in) & (i%%burn.in.check==TRUE) & (do.burn.in)){
        acceptance.rate.current <- sum(accept.vec[(i-burn.in.check):i]) / length(accept.vec[(i-burn.in.check):i])
        
        acceptRatio = acceptance.rate.current / target.accept
        #proposal.width = proposal.width * acceptRatio
        #print(acceptance.rate.current)
        #print(threshold)
        #print("")
        
        threshold <- threshold * acceptRatio
        
        if(threshold == 0){
          threshold <- runif(1)
        }
      }
    }
    
    coef.df <- results.df.coef
    se.df <- results.df.se
    
    coef.df.lim <- coef.df
    coef.df.lim$AIC_change <- coef.df.lim$AIC - min(coef.df.lim$AIC)
    coef.df.lim <- coef.df.lim[coef.df.lim$AIC_change < 500,]
    coef.df.lim$model.likelihood <- exp(-(1/2)*coef.df.lim$AIC_change)
    coef.df.lim$w <- coef.df.lim$model.likelihood/sum(coef.df.lim$model.likelihood)
    
    coef.df.noAIC <- coef.df[1:(ncol(coef.df)-1)]
    coef.df.noAIC <- coef.df.noAIC * coef.df.lim$w
    coef <- colSums(coef.df.noAIC)    
    
    var.df <- se.df^2
    
    se <- lapply(1:length(coef), function(i) sqrt(sum(coef.df.lim$w*(var.df[,i] + (coef.df[,i] - coef[i])^2))))
    se <- as.data.frame(se)
    names(se) <- names(coef)
    
    coef <- as.data.frame(t(coef))
    row.names(coef) <- "Coefficients"
    row.names(se) <- "Std. Error"
    
    return(list(coef=coef,
                se=se,
                coef.chain=coef.df,
                se.chain=se.df,
                probAidGuess.current=probAidGuess.best,
                acceptance.rate.current=acceptance.rate.current,
                threshold=threshold))
  }
  
  param_set <- paramSet(aidData=aidData, roiData=roiData, probAidAssume=runif(nrow(roiData)), aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set.bin <- (param_set > 0)*1
  
  probaid.current <- runif(nrow(roiData))
  
  mod.avg.chain <- lapply(1, bayes_model_avg, 
                          model.initial=model, 
                          probAidGuess.current=probaid.current, 
                          subcountyData=roiData, 
                          aidData=aidData, 
                          param_set.bin=param_set.bin, 
                          chain.length=iterations, 
                          threshold=threshold, 
                          burn.in=burn.in, 
                          burn.in.check=burn.in.check, 
                          target.accept=target.accept, 
                          proposal.width=proposal.width, 
                          do.burn.in=do.burn.in)
  
  mod.avg.chain <- mod.avg.chain[[1]]
  
  return(list(coefficients=mod.avg.chain$coef,
              se=mod.avg.chain$se,
              df=model$df,
              coef.chain=mod.avg.chain$coef.chain,
              se.chain=mod.avg.chain$coef.chain,
              geoSIMEXvariable=geoSIMEXvariable,
              iterations=iterations,
              residuals = model$residuals,
              prob.aid.best.fit = mod.avg.chain$probAidGuess.current,
              acceptance.rate.current=mod.avg.chain$acceptance.rate.current,
              threshold=mod.avg.chain$threshold))
}

modelAverageMaxLikeRandProb <- function(x, ...) UseMethod("modelAverageMaxLikeRandProb")

modelAverageMaxLikeRandProb.default <- function(iterations=500, 
                                                model,
                                                geoSIMEXvariable,
                                                roiData, 
                                                aidData,
                                                aid.amount, 
                                                roi.pc1.name="pc1.id", 
                                                roi.pc2.name="pc2.id", 
                                                roi.pc3.name="pc3.id", 
                                                roi.pc4.name="pc4.id", 
                                                roi.pc5.name="pc5.id", 
                                                roi.pc6.name="pc6.id", 
                                                aid.pc1.centroid.name="centroid.pc1.id",
                                                aid.precision.code="precision.code",
                                                binary=FALSE,
                                                threshold=.1, 
                                                burn.in=400, 
                                                burn.in.check=100, 
                                                target.accept=.9,
                                                proposal.width=1, 
                                                do.burn.in=TRUE,
                                                parallel=TRUE, 
                                                mc.cores=2){
  
  est <- modelAverageMaxLikeRandProb_est(iterations=iterations, 
                                         model=model,
                                         geoSIMEXvariable=geoSIMEXvariable,
                                         roiData=roiData, 
                                         aidData=aidData,
                                         aid.amount=aid.amount, 
                                         roi.pc1.name=roi.pc1.name, 
                                         roi.pc2.name=roi.pc2.name, 
                                         roi.pc3.name=roi.pc3.name, 
                                         roi.pc4.name=roi.pc4.name, 
                                         roi.pc5.name=roi.pc5.name, 
                                         roi.pc6.name=roi.pc6.name, 
                                         aid.pc1.centroid.name=aid.pc1.centroid.name,
                                         aid.precision.code=aid.precision.code,
                                         binary=binary,
                                         threshold, 
                                         burn.in, 
                                         burn.in.check, 
                                         target.accept,
                                         proposal.width, 
                                         do.burn.in,
                                         parallel=parallel, 
                                         mc.cores=mc.cores)
  
  est$call <- model$call
  
  class(est) <- "modelAverageMaxLikeRandProb"
  est
}

print.modelAverageMaxLikeRandProb <- function(x, ...){
  cat("Naive Model:\n")
  print(x$call)
  cat("\ngeoSIMEX-Variables: ", x$geoSIMEXvariable, sep="")
  cat("\nNumer of Iterations: ", x$iterations, sep= "")
  cat("\n\nCoefficients:\n")
  print(x$coefficients, row.names=FALSE)  
}

summary.modelAverageMaxLikeRandProb <- function(object, ...){
  
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

print.summary.modelAverageMaxLikeRandProb <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

##### SECTION 1B_2: MODEL AVERAGING BAYES - RAND PROB, 2 Vars #####

modelAverageMaxLikeRandProb2_est <- function(iterations, 
                                             model, 
                                             geoSIMEXvariable_1, 
                                             geoSIMEXvariable_2,
                                             roiData, 
                                             aidData_1,
                                             aidData_2, 
                                             aid.amount_1, 
                                             aid.amount_2, 
                                             roi.pc1.name, 
                                             roi.pc2.name, 
                                             roi.pc3.name, 
                                             roi.pc4.name, 
                                             roi.pc5.name, 
                                             roi.pc6.name, 
                                             aid.pc1.centroid.name_1,
                                             aid.pc1.centroid.name_2,
                                             aid.precision.code_1,
                                             aid.precision.code_2,
                                             binary,
                                             threshold, 
                                             burn.in, 
                                             burn.in.check, 
                                             target.accept,
                                             proposal.width, 
                                             do.burn.in,
                                             parallel, 
                                             mc.cores){
  
  dollar_expected_value_fast <- function(param_set.bin, dollar_set, probAidGuess.current){
    
    param_set.bin <- param_set.bin*probAidGuess.current
    param_set.bin.colsums <- colSums(param_set.bin)
    
    paramDollars <- lapply(1:length(dollar_set), function(i) param_set.bin[,i] / param_set.bin.colsums[i])
    paramDollars <- as.data.frame(paramDollars)
    
    paramDollars <- lapply(1:length(dollar_set), function(i) paramDollars[,i] * dollar_set[i])
    paramDollars <- as.data.frame(paramDollars)
    
    return( rowSums(paramDollars) )
  }
  
  bayes_model_avg <- function(k, 
                              model.initial, 
                              probAidGuess.current_1, 
                              probAidGuess.current_2, 
                              subcountyData, 
                              aidData_1, 
                              aidData_2, 
                              param_set.bin_1=param_set.bin_1, 
                              param_set.bin_2=param_set.bin_2, 
                              chain.length=chain.length, 
                              threshold=threshold, 
                              burn.in=burn.in, 
                              burn.in.check=burn.in.check, 
                              target.accept=target.accept, 
                              proposal.width=proposal.width, 
                              do.burn.in=do.burn.in){
    
    temp=k
    
    accept.vec <- matrix(NA, nrow=0, ncol=1)
    
    roiData[geoSIMEXvariable_1] <- dollar_expected_value_fast(param_set.bin=param_set.bin_1, dollar_set=aidData_1[,aid.amount_1], probAidGuess.current=probaid.current_1)
    roiData[geoSIMEXvariable_2] <- dollar_expected_value_fast(param_set.bin=param_set.bin_2, dollar_set=aidData_2[,aid.amount_2], probAidGuess.current=probaid.current_2)
    df.temp <- model.initial$model
    df.temp[geoSIMEXvariable_1] <- roiData[geoSIMEXvariable_1]
    df.temp[geoSIMEXvariable_2] <- roiData[geoSIMEXvariable_2]
    model <- update(model.initial, data = df.temp)
    
    likelihood.current <- logLik(model)
    
    probAidGuess.best_1 <- probAidGuess.current_1
    probAidGuess.best_2 <- probAidGuess.current_2
    aic.best <- extractAIC(model)[2]
    
    coefs <- summary(model)$coefficients[,1]
    se <- summary(model)$coefficients[,2]
    coefs <- as.data.frame(t(coefs))
    se <- as.data.frame(t(se))
    AIC <- as.data.frame(as.matrix(extractAIC(model)[2], nrow=1, ncol=1))
    names(AIC) <- "AIC"
    results.df.coef <- cbind(coefs, AIC)
    results.df.se <- cbind(se, AIC)
    
    for(i in 1:iterations){
      
      # Also to Try:
      # series of if statements for how update... start with large groups, then work down.
      
      # Proposed probabilities and model      
      probAidGuess.proposed_1 <- probAidGuess.current_1 + rnorm(length(probAidGuess.current_1),mean=0,sd=proposal.width)
      probAidGuess.proposed_2 <- probAidGuess.current_2 + rnorm(length(probAidGuess.current_2),mean=0,sd=proposal.width)
      
      probAidGuess.proposed_1[probAidGuess.proposed_1 < 0] <- runif(length(probAidGuess.proposed_1[probAidGuess.proposed_1 < 0]))*proposal.width
      probAidGuess.proposed_2[probAidGuess.proposed_2 < 0] <- runif(length(probAidGuess.proposed_2[probAidGuess.proposed_2 < 0]))*proposal.width
      
      roiData[geoSIMEXvariable_1] <- dollar_expected_value_fast(param_set.bin=param_set.bin_1, dollar_set=aidData_1[,aid.amount_1], probAidGuess.current=probAidGuess.proposed_1)
      roiData[geoSIMEXvariable_2] <- dollar_expected_value_fast(param_set.bin=param_set.bin_2, dollar_set=aidData_2[,aid.amount_2], probAidGuess.current=probAidGuess.proposed_2)
      
      df.temp <- model.initial$model
      df.temp[geoSIMEXvariable_1] <- roiData[geoSIMEXvariable_1]
      df.temp[geoSIMEXvariable_2] <- roiData[geoSIMEXvariable_2]
      model <- update(model.initial, data = df.temp)
      likelihood.proposed <- logLik(model)
      
      # Update 
      c <- likelihood.proposed - likelihood.current
      accept <- threshold <= exp(c)
      accept.vec <- rbind(accept.vec, accept)
      
      if(accept){
        
        likelihood.current <- likelihood.proposed
        
        probAidGuess.current_1 <- probAidGuess.proposed_1
        probAidGuess.current_2 <- probAidGuess.proposed_2
        
        coefs <- summary(model)$coefficients[,1]
        se <- summary(model)$coefficients[,2]
        coefs <- as.data.frame(t(coefs))
        se <- as.data.frame(t(se))
        AIC <- as.data.frame(as.matrix(extractAIC(model)[2], nrow=1, ncol=1))
        names(AIC) <- "AIC"
        results.df.coef.temp <- cbind(coefs, AIC)
        results.df.se.temp <- cbind(se, AIC)
        
        results.df.coef <- rbind(results.df.coef, results.df.coef.temp)
        results.df.se <- rbind(results.df.se, results.df.se.temp)
        
        # TAKING LAST VALUE OF CHAIN
        #probAidGuess.best_1 <- probAidGuess.proposed_1
        #probAidGuess.best_2 <- probAidGuess.proposed_2
      }
      
      #UPDATE PROB AID GUESS
      #If has lower AIC than best, update aic.best and probAidGuess.best
      if(extractAIC(model)[2] < aic.best){
        #if(abs(extractAIC(model)[2]) < abs(aic.best)){  
        aic.best <- AIC(model)
        probAidGuess.best_1 <- probAidGuess.proposed_1
        probAidGuess.best_2 <- probAidGuess.proposed_2
      }
      acceptance.rate.current <- NA
      # burn-in: adjust threshold for accepting
      if((i >= burn.in.check+1) & (i < burn.in) & (i%%burn.in.check==TRUE) & (do.burn.in)){
        acceptance.rate.current <- sum(accept.vec[(i-burn.in.check):i]) / length(accept.vec[(i-burn.in.check):i])
        
        acceptRatio = acceptance.rate.current / target.accept
        #proposal.width = proposal.width * acceptRatio
        #print(acceptance.rate.current)
        #print(threshold)
        #print("")
        
        threshold <- threshold * acceptRatio
        
        if(threshold == 0){
          threshold <- runif(1)
        }
      }
    }
    
    coef.df <- results.df.coef
    se.df <- results.df.se
    
    coef.df.lim <- coef.df
    coef.df.lim$AIC_change <- coef.df.lim$AIC - min(coef.df.lim$AIC)
    coef.df.lim <- coef.df.lim[coef.df.lim$AIC_change < 500,]
    coef.df.lim$model.likelihood <- exp(-(1/2)*coef.df.lim$AIC_change)
    coef.df.lim$w <- coef.df.lim$model.likelihood/sum(coef.df.lim$model.likelihood)
    
    coef.df.noAIC <- coef.df[1:(ncol(coef.df)-1)]
    coef.df.noAIC <- coef.df.noAIC * coef.df.lim$w
    coef <- colSums(coef.df.noAIC)    
    
    var.df <- se.df^2
    
    se <- lapply(1:length(coef), function(i) sqrt(sum(coef.df.lim$w*(var.df[,i] + (coef.df[,i] - coef[i])^2))))
    se <- as.data.frame(se)
    names(se) <- names(coef)
    
    coef <- as.data.frame(t(coef))
    row.names(coef) <- "Coefficients"
    row.names(se) <- "Std. Error"
    
    return(list(coef=coef,
                se=se,
                coef.chain=coef.df,
                se.chain=se.df,
                probAidGuess.current_1=probAidGuess.best_1,
                probAidGuess.current_2=probAidGuess.best_2,
                acceptance.rate.current=acceptance.rate.current,
                threshold=threshold))
  }
  
  param_set_1 <- paramSet(aidData=aidData_1, roiData=roiData, probAidAssume=runif(nrow(roiData)), aid.precision.code=aid.precision.code_1, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_1)
  param_set_2 <- paramSet(aidData=aidData_2, roiData=roiData, probAidAssume=runif(nrow(roiData)), aid.precision.code=aid.precision.code_2, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name= roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name_2)
  param_set.bin_1 <- (param_set_1 > 0)*1
  param_set.bin_2 <- (param_set_2 > 0)*1
  
  probaid.current_1 <- runif(nrow(roiData))
  probaid.current_2 <- runif(nrow(roiData))
  
  mod.avg.chain <- lapply(1, bayes_model_avg, 
                          model.initial=model, 
                          probAidGuess.current_1=probaid.current_1, 
                          probAidGuess.current_2=probaid.current_2, 
                          subcountyData=roiData, 
                          aidData_1=aidData_1, 
                          aidData_2=aidData_2, 
                          param_set.bin_1=param_set.bin_1, 
                          param_set.bin_2=param_set.bin_2, 
                          chain.length=iterations, 
                          threshold=threshold, 
                          burn.in=burn.in, 
                          burn.in.check=burn.in.check, 
                          target.accept=target.accept, 
                          proposal.width=proposal.width, 
                          do.burn.in=do.burn.in)
  
  mod.avg.chain <- mod.avg.chain[[1]]
  
  return(list(coefficients=mod.avg.chain$coef,
              se=mod.avg.chain$se,
              df=model$df,
              coef.chain=mod.avg.chain$coef.chain,
              se.chain=mod.avg.chain$coef.chain,
              geoSIMEXvariable_1=geoSIMEXvariable_1,
              geoSIMEXvariable_2=geoSIMEXvariable_2,
              iterations=iterations,
              residuals = model$residuals,
              prob.aid.best.fit_1 = mod.avg.chain$probAidGuess.current_1,
              prob.aid.best.fit_2 = mod.avg.chain$probAidGuess.current_2,
              acceptance.rate.current=mod.avg.chain$acceptance.rate.current,
              threshold=mod.avg.chain$threshold))
}

modelAverageMaxLikeRandProb2 <- function(x, ...) UseMethod("modelAverageMaxLikeRandProb2")

modelAverageMaxLikeRandProb2.default <- function(iterations=500, 
                                                 model,
                                                 geoSIMEXvariable_1,
                                                 geoSIMEXvariable_2,
                                                 roiData, 
                                                 aidData_1,
                                                 aidData_2,
                                                 aid.amount_1, 
                                                 aid.amount_2, 
                                                 roi.pc1.name="pc1.id", 
                                                 roi.pc2.name="pc2.id", 
                                                 roi.pc3.name="pc3.id", 
                                                 roi.pc4.name="pc4.id", 
                                                 roi.pc5.name="pc5.id", 
                                                 roi.pc6.name="pc6.id", 
                                                 aid.pc1.centroid.name_1="centroid.pc1.id",
                                                 aid.pc1.centroid.name_2="centroid.pc1.id",
                                                 aid.precision.code_1="precision.code",
                                                 aid.precision.code_2="precision.code",
                                                 binary=FALSE,
                                                 threshold=.1, 
                                                 burn.in=400, 
                                                 burn.in.check=100, 
                                                 target.accept=.9,
                                                 proposal.width=1, 
                                                 do.burn.in=TRUE,
                                                 parallel=TRUE, 
                                                 mc.cores=2){
  
  est <- modelAverageMaxLikeRandProb2_est(iterations=iterations, 
                                          model=model,
                                          geoSIMEXvariable_1=geoSIMEXvariable_1,
                                          geoSIMEXvariable_2=geoSIMEXvariable_2,
                                          roiData=roiData, 
                                          aidData_1=aidData_1,
                                          aidData_2=aidData_2,
                                          aid.amount_1=aid.amount_1, 
                                          aid.amount_2=aid.amount_2, 
                                          roi.pc1.name=roi.pc1.name, 
                                          roi.pc2.name=roi.pc2.name, 
                                          roi.pc3.name=roi.pc3.name, 
                                          roi.pc4.name=roi.pc4.name, 
                                          roi.pc5.name=roi.pc5.name, 
                                          roi.pc6.name=roi.pc6.name, 
                                          aid.pc1.centroid.name_1=aid.pc1.centroid.name_1,
                                          aid.pc1.centroid.name_2=aid.pc1.centroid.name_2,
                                          aid.precision.code_1=aid.precision.code_1,
                                          aid.precision.code_2=aid.precision.code_2,
                                          binary=binary,
                                          threshold, 
                                          burn.in, 
                                          burn.in.check, 
                                          target.accept,
                                          proposal.width, 
                                          do.burn.in,
                                          parallel=parallel, 
                                          mc.cores=mc.cores)
  
  est$call <- model$call
  
  class(est) <- "modelAverageMaxLikeRandProb2"
  est
}

print.modelAverageMaxLikeRandProb2 <- function(x, ...){
  cat("Naive Model:\n")
  print(x$call)
  cat("\ngeoSIMEX-Variables: ", x$geoSIMEXvariable_1, " ",x$geoSIMEXvariable_1, sep="")
  cat("\nNumer of Iterations: ", x$iterations, sep= "")
  cat("\n\nCoefficients:\n")
  print(x$coefficients, row.names=FALSE)  
}

summary.modelAverageMaxLikeRandProb2 <- function(object, ...){
  
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
  
  class(res) <- "summary.modelAverageRandProb2"
  res
}

print.summary.modelAverageMaxLikeRandProb2 <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}


##### SECTION 1C: MODEL AVERAGING, RAND ROI #####

# Model Averaging Randomly Drop Project in ROI
modelAverageRandROI_est <- function(iterations, 
                                    aidData, 
                                    roiData, 
                                    roi.area, 
                                    aid.amount, 
                                    model, 
                                    geoSIMEXvariable, 
                                    binary, 
                                    aid.precision.code, 
                                    roi.pc1.name, 
                                    roi.pc2.name, 
                                    roi.pc3.name, 
                                    roi.pc4.name, 
                                    roi.pc5.name, 
                                    roi.pc6.name, 
                                    aid.pc1.centroid.name){
  
  # Calculating paramSet
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=roiData[,roi.area], aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  
  model.list <- lapply(1:iterations, geoSimulate_realization, param_set=param_set, roiData=roiData, aid.amount=aid.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary, aidData=aidData)
  
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
              residuals=model$residuals,
              geoSIMEXvariable=geoSIMEXvariable,
              iterations=iterations))
}

modelAverageRandROI <- function(x, ...) UseMethod("modelAverageRandROI")

modelAverageRandROI.default <- function(iterations=500, 
                                        aidData, 
                                        roiData, 
                                        roi.area, 
                                        aid.amount, 
                                        model, 
                                        geoSIMEXvariable, 
                                        binary=FALSE, 
                                        aid.precision.code="precision.code", 
                                        roi.pc1.name="pc1.id", 
                                        roi.pc2.name="pc2.id", 
                                        roi.pc3.name="pc3.id", 
                                        roi.pc4.name="pc4.id", 
                                        roi.pc5.name="pc5.id", 
                                        roi.pc6.name="pc6.id", 
                                        aid.pc1.centroid.name="centroid.pc1.id"){
  
  est <- modelAverageRandROI_est(iterations=iterations, 
                                 aidData=aidData, 
                                 roiData=roiData, 
                                 roi.area=roi.area, 
                                 aid.amount=aid.amount, 
                                 model=model, 
                                 geoSIMEXvariable=geoSIMEXvariable, 
                                 binary=binary, 
                                 aid.precision.code=aid.precision.code, 
                                 roi.pc1.name=roi.pc1.name, 
                                 roi.pc2.name=roi.pc2.name, 
                                 roi.pc3.name=roi.pc3.name, 
                                 roi.pc4.name=roi.pc4.name, 
                                 roi.pc5.name=roi.pc5.name, 
                                 roi.pc6.name=roi.pc6.name, 
                                 aid.pc1.centroid.name=aid.pc1.centroid.name)
  
  est$call <- model$call
  
  class(est) <- "modelAverageRandROI"
  est
}

print.modelAverageRandROI <- function(x, ...){
  cat("Naive Model:\n")
  print(x$call)
  cat("\ngeoSIMEX-Variables: ", x$geoSIMEXvariable, sep="")
  cat("\nNumer of Iterations: ", x$iterations, sep= "")
  cat("\n\nCoefficients:\n")
  print(x$coefficients, row.names=FALSE)  
}

summary.modelAverageRandROI <- function(object, ...){
  
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

print.summary.modelAverageRandROI <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

##### SECTION 2: FUNCTIONS TO BE INCLUDED IN PACKAGE #####
# Calculates expected value of aid for each ROI
# 
# @title Expected Value Aid
# 
# @author AidData
# 
# @description
# \code{expected_aid_ROI} Calculates expected value of aid.
# 
# @usage
# expected_aid_ROI(aidData, roiData, probAidAssume, dollar_set, aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
#  
# @param roiData name of dataframe of ROI data 
# @param aidData name of dataframe of aid project data
# @param aid.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
# @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#  
# @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website].
expected_aid_ROI <- function(aidData, roiData, probAidAssume, dollar_set, aid.precision.code, roi.pc1.name, roi.pc2.name, roi.pc3.name, roi.pc4.name, roi.pc5.name, roi.pc6.name, aid.pc1.centroid.name){
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  aid <- dollar_expected_value(param_set=param_set, dollar_set=dollar_set)
  return(aid)
}

prob_aid_ROI <- function(aidData, roiData, probAidAssume, dollar_set, aid.precision.code, roi.pc1.name, roi.pc2.name, roi.pc3.name, roi.pc4.name, roi.pc5.name, roi.pc6.name, aid.pc1.centroid.name){
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  param_set_probNoAid <- 1 - param_set
  prob.0.prjs <- apply(param_set_probNoAid, 1, function(x) prod(x))
  prob.atLeast1.prjs <- 1 - prob.0.prjs
  
  #prob.atLeast1.prjs[prob.atLeast1.prjs < 0.5] <- 0
  #prob.atLeast1.prjs[prob.atLeast1.prjs >= 0.5] <- 1
  return(prob.atLeast1.prjs)
}

# @title Spatial Uncertainty (Lambda)
# 
# @author AidData
# 
# @description
# \code{calc_lambda} Calculates spatial uncertainty (lambda) of aid project dataset
# 
# @usage
# calc_lambda(aidData, roiData, roi.area="area", aid.precision.code="precision.code", roi.pc1.name="pc1.id", roi.pc2.name="pc2.id", roi.pc3.name="pc3.id", roi.pc4.name="pc4.id", roi.pc5.name="pc5.id", roi.pc6.name="pc6.id", aid.pc1.centroid.name="centroid.pc1.id")
#  
# @param roiData name of dataframe of ROI data 
# @param aidData name of dataframe of aid project data
# @param aid.amount character containing the name of the variable in the aidData dataset which contains aid amounts (e.g., commitment, disbursement). Set value to 1 if interested in number of aid projects rather than dollars.
# @param roi.area character containing the name of the variable in the ROI dataset which contains areas of ROIs. "area" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc1.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc2.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc2.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc3.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc3.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc4.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc4.id" is the default name in datasets produced by AidData's data extraction tool
# @param roi.pc6.name character containing the name of the variable in the ROI dataset which contains names or IDs of the precision code 1 spatial area that each ROI falls within. "pc6.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.pc1.centroid.name character containing the name of the variable in the aidData dataset which contains names or IDs of a precision code 1 spatial area that the aid project falls within. "centroid.pc1.id" is the default name in datasets produced by AidData's data extraction tool
# @param aid.precision.code character containing the name of the variable in the aidData dataset which contains precision codes for each project. "pc1.id" is the default name in datasets produced by AidData's data extraction tool
#  
# @examples 
# lambda.naive <- calc_lambda(aidData=aidData.sub@data, roiData=roiData)
# 
# @note The function is built to work with data from AidData's data extration tool. The extraction tool can be accessed here: [provide website].

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
  
  param_set.a = lapply(1:nrow(aidData), paramCol, aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  
  return(param_set) 
}

dollar_expected_value <- function(param_set, dollar_set){    
  dollar_set[is.na(dollar_set)] <- 0
  paramDollars <- lapply(1:length(dollar_set), function(i) param_set[,i] * dollar_set[i])
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

geoSimulateError <- function(probIncPC, aidData=aidData, roiData=roiData, probAidAssume=probAid, PC_researcherSees=precision.code.original, maxLambda_denom=maxLambda_denom, roi.area=roi.area, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name, aid.amount=aid.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary, sim_pc1=sim_pc1){
    
  # Initializing Results Matrix
  results <- matrix(NA, ncol=2,nrow=1)
  results <- as.data.frame(results)
  names(results) <- c("lambda","aid.expected.coef")
  
  # Simulating Error
  aidData[aid.precision.code] <- PC_researcherSees
  probStayPC <- 1 - probIncPC
  
  # ... Precision Code 4s
  aidData[aid.precision.code][aidData[aid.precision.code] == 4] <- sample(size=length(aidData[aid.precision.code][aidData[aid.precision.code] == 4]), 
                                                                          x=c(4,6), 
                                                                          prob=c(probStayPC,probIncPC), 
                                                                          replace=TRUE)
  # ... Precision Code 3s
  aidData[aid.precision.code][aidData[aid.precision.code] == 3] <- sample(size=length(aidData[aid.precision.code][aidData[aid.precision.code] == 3]), 
                                                                          x=c(3,4,6), 
                                                                          prob=c(probStayPC,probIncPC/2,probIncPC/2), 
                                                                          replace=TRUE)
  
  # ... Precision Code 2s
  aidData[aid.precision.code][aidData[aid.precision.code] == 2] <- sample(size=length(aidData[aid.precision.code][aidData[aid.precision.code] == 2]), 
                                                                          x=c(2,3,4,6), 
                                                                          prob=c(probStayPC,probIncPC/3,probIncPC/3,probIncPC/3), 
                                                                          replace=TRUE)
  
  # ... Precision Code 1s
  if(sim_pc1){
    aidData[aid.precision.code][aidData[aid.precision.code] == 1] <- sample(size=length(aidData[aid.precision.code][aidData[aid.precision.code] == 1]), 
                                                                            x=c(1,2,3,4,6), 
                                                                            prob=c(probStayPC,probIncPC/4,probIncPC/4,probIncPC/4,probIncPC/4), 
                                                                            replace=TRUE)
  }
  
  # Calculating paramSet
  param_set = paramSet(aidData=aidData, roiData=roiData, probAidAssume=probAidAssume, aid.precision.code=aid.precision.code, roi.pc1.name=roi.pc1.name, roi.pc2.name=roi.pc2.name, roi.pc3.name=roi.pc3.name, roi.pc4.name=roi.pc4.name, roi.pc5.name=roi.pc5.name, roi.pc6.name=roi.pc6.name, aid.pc1.centroid.name=aid.pc1.centroid.name)
  
  # Update aid variable
  if(binary){
    roiData[geoSIMEXvariable] <- prob_aid(param_set=param_set)
  } else{
    roiData[geoSIMEXvariable] <- dollar_expected_value(param_set=param_set, dollar_set=as.matrix(aidData[aid.amount]))
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


bootIter <- function(i, geoSimulateError.results.df, geoSimulateError.results.df.se, bins, numFromBin){
  
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
  
  numIter <- (maxLambda - minLambda) / binSize
  for(i in 1:ceiling(numIter)){
    
    results.boot <- rbind(results.boot, geoSimulateError.results.df[(geoSimulateError.results.df$lambda >= binSize_lb) & (geoSimulateError.results.df$lambda <= binSize_ub),][1:numFromBin,])
    results.boot.se <- rbind(results.boot.se, geoSimulateError.results.df.se[(geoSimulateError.results.df.se$lambda >= binSize_lb) & (geoSimulateError.results.df.se$lambda <= binSize_ub),][1:numFromBin,])
    
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
  
  numVars <- ncol(results.boot.se) - 3
  se.geoSIMEX.boot <- matrix(NA, nrow=1, ncol=numVars)
  for(i in 1:numVars){
    se.geoSIMEX.boot[i] <- summary(lm(as.matrix(results.boot.se[i]) ~ lambda + lambda_sq, data = results.boot.se))$coefficients[1]
  }
  se.geoSIMEX.boot <- as.data.frame(se.geoSIMEX.boot)
  names(se.geoSIMEX.boot) <- head(names(results.boot.se), -3)
  
  return(list(coefs.geoSIMEX.boot=coefs.geoSIMEX.boot,
              se.geoSIMEX.boot=se.geoSIMEX.boot))
}

bootIter2 <- function(i, geoSimulateError.results.df, geoSimulateError.results.df.se, average_lambda){
  
  geoSimulateError.results.df$rand <- runif(nrow(geoSimulateError.results.df))
  geoSimulateError.results.df.se$rand <- geoSimulateError.results.df$rand
  
  geoSimulateError.results.df <- geoSimulateError.results.df[order(geoSimulateError.results.df$rand),] 
  geoSimulateError.results.df.se <- geoSimulateError.results.df.se[order(geoSimulateError.results.df.se$rand),] 
  
  results.boot <- matrix(NA,nrow=0,ncol=ncol(geoSimulateError.results.df))
  results.boot <- as.data.frame(results.boot)
  
  results.boot.se <- matrix(NA,nrow=0,ncol=ncol(geoSimulateError.results.df.se))
  results.boot.se <- as.data.frame(results.boot.se)
  
  boot.coefs <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(geoSimulateError.results.df)))
  names(boot.coefs) <- names(geoSimulateError.results.df)
  
  boot.se <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(geoSimulateError.results.df.se)))
  names(boot.se) <- names(geoSimulateError.results.df.se)
  
  # Randomly Pulling from Results Matrix to Get Subsample for Bootstrapping 
  for(b in names(table(geoSimulateError.results.df$bin_12))){
    temp.df <- geoSimulateError.results.df[geoSimulateError.results.df$bin_12 == b,]
    boot.coefs <- rbind(boot.coefs, temp.df[1,])
  }
  
  for(b in names(table(geoSimulateError.results.df.se$bin_12))){
    temp.df <- geoSimulateError.results.df.se[geoSimulateError.results.df.se$bin_12 == b,]
    boot.se <- rbind(boot.se, temp.df[1,])
  }
  
  # Extrapolating 
  if(average_lambda){
    numVars <- ncol(boot.coefs) - 8
    coefs.geoSIMEX.boot <- matrix(NA, nrow=1, ncol=numVars)
    for(i in 1:numVars){
      coefs.geoSIMEX.boot[i] <- summary(lm(as.matrix(boot.coefs[i]) ~ lambda_1 + lambda_1_sq, data = boot.coefs))$coefficients[1]
    }
    coefs.geoSIMEX.boot <- as.data.frame(coefs.geoSIMEX.boot)
    names(coefs.geoSIMEX.boot) <- head(names(boot.coefs), -8)
    
    numVars <- ncol(boot.se) - 8
    se.geoSIMEX.boot <- matrix(NA, nrow=1, ncol=numVars)
    for(i in 1:numVars){
      se.geoSIMEX.boot[i] <- summary(lm(as.matrix(boot.se[i]) ~ lambda_1 + lambda_1_sq, data = boot.se))$coefficients[1]
    }
    se.geoSIMEX.boot <- as.data.frame(se.geoSIMEX.boot)
    names(se.geoSIMEX.boot) <- head(names(boot.se), -8)
    
  } else {
    numVars <- ncol(boot.coefs) - 8
    coefs.geoSIMEX.boot <- matrix(NA, nrow=1, ncol=numVars)
    for(i in 1:numVars){
      coefs.geoSIMEX.boot[i] <- summary(lm(as.matrix(boot.coefs[i]) ~ lambda_1 + lambda_2 + lambda_1_sq + lambda_2_sq, data = boot.coefs))$coefficients[1]
    }
    coefs.geoSIMEX.boot <- as.data.frame(coefs.geoSIMEX.boot)
    names(coefs.geoSIMEX.boot) <- head(names(boot.coefs), -8)
    
    numVars <- ncol(boot.se) - 8
    se.geoSIMEX.boot <- matrix(NA, nrow=1, ncol=numVars)
    for(i in 1:numVars){
      se.geoSIMEX.boot[i] <- summary(lm(as.matrix(boot.se[i]) ~ lambda_1 + lambda_2 + lambda_1_sq + lambda_2_sq, data = boot.se))$coefficients[1]
    }
    se.geoSIMEX.boot <- as.data.frame(se.geoSIMEX.boot)
    names(se.geoSIMEX.boot) <- head(names(boot.se), -8)
  }
  
  return(list(coefs.geoSIMEX.boot=coefs.geoSIMEX.boot,
              se.geoSIMEX.boot=se.geoSIMEX.boot))
}


# Model Averaging Change Probability 


model_rand_prob <- function(j, param_set.bin=param_set.bin, aidData=aidData, roiData=roiData, aid.amount=aid.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary){
  
  # Update aid variable
  if(binary){
    # WILL NEED THE RANDOMLY GENERATING PARAM SET HERE!
  } else{
    roiData[geoSIMEXvariable] <- dollar_expected_value_randProb(param_set.bin=param_set.bin, dollar_set=as.matrix(aidData[aid.amount]))
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



dollar_expected_value_randProb <- function(j, param_set.bin, dollar_set){
  
  dist.type <- sample(1,x=c(1,2),prob=c(1,1))
  if(dist.type==1){
    probAidGuess.current <- runif(nrow(param_set.bin))
  }
  
  if(dist.type==2){
    probAidGuess.current <- rgamma(nrow(param_set.bin), shape=runif(1,1,20))
  }
  
  param_set.bin <- param_set.bin*probAidGuess.current
  param_set.bin.colsums <- colSums(param_set.bin)
  
  paramDollars <- lapply(1:length(dollar_set), function(i) param_set.bin[,i] / param_set.bin.colsums[i])
  paramDollars <- as.data.frame(paramDollars)
  
  paramDollars <- lapply(1:length(dollar_set), function(i) paramDollars[,i] * dollar_set[i])
  paramDollars <- as.data.frame(paramDollars)
  
  return(rowSums(paramDollars))
}


geoSimulate_realization <- function(j, param_set=param_set, roiData=roiData, aid.amount=aid.amount, model=model, geoSIMEXvariable=geoSIMEXvariable, binary=binary, aidData=aidData){
  
  # Update aid variable
  if(binary){
    temp.var <- realization_of_aid(param_set=param_set, dollar_set=as.matrix(aidData[aid.amount]))
    temp.var[temp.var > 0] <- 1
    roiData[geoSIMEXvariable] <- temp.var
  } else{
    roiData[geoSIMEXvariable] <- realization_of_aid(param_set=param_set, dollar_set=as.matrix(aidData[aid.amount]))
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

realization_of_aid <- function(param_set, dollar_set){
  
  param_set_1pjrPerROI <- function(i, param_set=param_set, dollar_set=dollar_set){
    id <- sample(size=1,x=1:length(param_set[,i]), prob=param_set[,i])
    col.temp <- as.data.frame(matrix(0,nrow=length(param_set[,i]), ncol=1))
    col.temp[id,1] <- 1
    return(col.temp)
  }
  
  param_set_realization <- as.data.frame(lapply(1:ncol(param_set), param_set_1pjrPerROI, param_set=param_set, dollar_set=dollar_set))
  
  paramDollars <- lapply(1:length(dollar_set), function(i) param_set_realization[,i] * dollar_set[i])
  paramDollars <- as.data.frame(paramDollars)
  dollars_realization <- rowSums(paramDollars)
  return(dollars_realization)
}

##### SECTION 4: SUBSETTING AID DATA #####

subset.aiddata <- function(json){
  
  # Import dataset. Dataset names must be loaded from Rpackage and be the same that appears in JSON.
  #geo.data <- eval(parse(text=json$dataset))
  
  # For now, have all the datasets loaded on github; will pull from there. 
  geo.data.name <- json$dataset
  if(geo.data.name == "colombiaaims_geocodedresearchrelease_level1_v1_1_1"){
    # geo.data <- COLUMBIA DATA
  }
  
  # By default, just loading UGA data in here now.
  geo.data <- read.csv("https://raw.githubusercontent.com/ramarty/geoSIMEX/master/Example/UgandaAMP_GeocodedResearchRelease_Level1_v1.3/data/level_1a.csv")
  
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




