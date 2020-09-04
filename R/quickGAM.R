#'quickGAM
#'
#'@description
#'\code{quickGAM} is a function to fit a GAM and perform the related analyses described in Simpson (2018). This function essentially applies the code provided in that article's Supplementary Information, to a pair of vectors.
#'
#'@usage quickGAM(x, y, xLabel = NULL, yLabel = NULL, yRange = NULL, skipAuto = FALSE)
#'
#'@param x vector containing the independent data/linear predictor
#'@param y vector containing the dependent data/response variable
#'@param xLabel (defaults to NULL): optional custom string to ensure consistent x-axis labels across the plots
#'@param yLabel (defaults to NULL): optional custom string to ensure consistent y-axis labels across the plots
#'@param yRange (defaults to NULL): optional custom range to ensure consistent y-axis across all plots
#'@param skipAuto (Defaults to FALSE): logical to indicate whether to skip the autocorrelation calculation and plot
#'
#'@details Two input vectors are required (x and y, the independent and dependent variable respectively). The optional arguments allow customization of the x and y axes labels to ensure they remain consistent across all graphs. Similarly the y-range can be manually defined.
#'
#'The function produces six plots:
#'\itemize{
#' \item Basic scatterplot of the data
#' \item Estimated CAR1 process from the GAM fitted to the data ( i.e. autocorrelation check)
#' \item GAM-based trends fitted to the data
#' \item Estimated trends (w 20 random draws of the GAM fit to the data)
#' \item 95\% simultaneous confidence intervals (light grey) and across-the-function confidence intervals (dark grey) on the estimated trends.
#' \item Estimated first derivatives (black lines) and 95\% simultaneous confidence intervals of the GAM trends fitted to the data. Where the simultaneous interval does not include 0, the models detect significant temporal change in the response.
#'   }
#'
#'@examples
#'#Fit a GAM to the example chlorophyll a data in 'vrsChlaInput'
#'data(vrsChlaInput)
#'plotData <- vrsChla(vrsChlaInput)
#'#Values/labels to be preserved across all plots
#'indepVarLabel <- "Core Depth (cm)"
#'depenVarLabel <- expression(
#'  "VRS-Inferred Chlorophyll "~italic("a")~" (mg"%.%"g"^"-1"*textstyle(")"))
#'depenVarRange <- c(0, 0.0825)
#'
#'#NOTE: drop=FALSE when reading in data columns preserves column names, allowing their use in plots.
#'quickGAM(plotData[,1, drop=FALSE],plotData[,2, drop=FALSE], xLabel = indepVarLabel, 
#'yLabel = depenVarLabel, yRange = depenVarRange)
#'
#'@author Adam Jeziorski
#' 
#'@references Simpson GL (2018) Modelling Palaeoecological Time Series Using Generalised Additive Models. Frontiers in Ecology and Evolution [https://dx.doi.org/10.3389/fevo.2018.00149]
#'
#'@export

quickGAM <- function(x, y, xLabel=NULL, yLabel=NULL, yRange=NULL, skipAuto = FALSE){
  #Set default ggplot theme
  theme_set(theme_bw())
  
  plotData <- cbind(x,y)
  plotData <- as.data.frame(plotData)
  
  #If xLabel and yLabel are not explicitly provided, use the column names of the dataframe
  if (is.null(xLabel) == TRUE){
    indepVarLabel <- colnames(plotData)[1]
  } else {indepVarLabel <- xLabel}
  
  if (is.null(yLabel) == TRUE){
    depenVarLabel <- colnames(plotData)[2]
  } else {depenVarLabel <- yLabel}
  
  #If yRange not explicitly provided, determine from min/max values of y
  if (is.null(yRange) == TRUE){
    depenVarRange <- c((min(y) - ((max(y)-min(y))*0.5)), (max(y) + ((max(y)-min(y))*0.5)))
  } else {depenVarRange <- yRange}

  
  #Overwrite columns names of plotData to allow direct references
  colnames(plotData) <- c("indep", "depen")
  
  #Initial Plot
  initialPlot <- ggplot(plotData, aes(x=indep, y=depen)) +
    geom_point() +
    labs(y=depenVarLabel, x= indepVarLabel) +
    ylim(depenVarRange)
  print(initialPlot)
  
  #Fit a GAM to the dataset ####
  #Simple GAM fit the dataset
  #'depen' is the response variable, 'indep' is the linear predictor. The linear predictor part of the model comtains a smooth function of the predictor variable in the data set, controlled by the "s" function, the smoothness selection method is set to 'REML'.
  m <- gam(depen ~ s(indep, k=15), data=plotData, method="REML")
  
  #GAM plus CAR(1) process fit to the dataset using gamm()
  #Model assumes the residuals are independent, so to account for residual temporal autcorrelation, can include a continuous-time first-order autoregressive (CAR(1)) process in the model residuals
  mod <- gamm(depen ~ s(indep, k=15), data=plotData, correlation = corCAR1(form = ~ indep), method="REML")
  

  
  #Summary object
  #The returned object includes both the linear mixed model and GAM elements of the model, need to access these elements individually.
  #Output shows estimated complexity of the fitted smooth, in terms of the effective degreees of freedom of the spline. An associated  F statistic and test of the null hypothesis are also shown.
  summary(mod$gam)
  
  #Plot CAR(1) process (can be bypassed using the skipAuto argument)
  if (skipAuto == F){
    
    #Estimate of phi and confidence interval
    #Estimated value of phi for the CAR(1) can be extracted from the fit model via the $lme component.
    smallPhi <- intervals(mod$lme, which = "var-cov")$corStruct
    smallPhi
    
    #In example 1 (chlorophyll a data), exponential decline in correlation with increasing separation is evident, little dependence once samples are ~5 cm apart
    S <-seq(0,50, length=100)
    car1 <- setNames(as.data.frame(t(outer(smallPhi, S, FUN = `^`)[1, , ])), c("Lower","Correlation","Upper"))
    car1 <- transform(car1, S = S)
  
    car1Plt <- ggplot(car1, aes(x = S, y = Correlation)) +
      geom_ribbon(aes(ymax = Upper, ymin = Lower), fill = "black", alpha = 0.2) +
      geom_line() +
      ylab(expression(italic(c) * (list(h, varphi)))) +
      xlab(bquote(h ~"("*.(indepVarLabel)*")"))
    print(car1Plt)
  }
  
  #Plot the fitted GAM
  #General idea is to predict from the fitted model for a fine grid of points over the range of the linear predictor (i.e. time or depth). The plot is of the trend with an ~95% point-wise confidence interval that assumes asymptotic normality.
  #Define the number of points at which to evaluate the smooth
  N <- 300 
  
  #Create new data to predict at; 200 evenly-spaced values over the independent variable `indep`
  newIndep <- with(plotData, data.frame(indep = seq(min(indep), max(indep), length.out = 200)))
  
  #Predict from the fitted model; note we predict from the $gam part
  newIndep <- cbind(newIndep, data.frame(predict(mod$gam, newIndep, se.fit = TRUE)))
  
  #Create the confidence interval
  crit.t <- qt(0.975, df = df.residual(mod$gam))
  newIndep <- transform(newIndep, upper = fit + (crit.t * se.fit), lower = fit - (crit.t * se.fit))
  
  #Plot the estimated trend
  fittedPlot <- ggplot(newIndep, aes(x = indep, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, x = indep), alpha = 0.2,
                inherit.aes = FALSE, fill = "black") +
    geom_point(data = plotData, mapping = aes(x = indep, y = depen),
               inherit.aes = FALSE) +
    geom_line() +
    labs(y = depenVarLabel, x = indepVarLabel) +
    ylim(depenVarRange)
  print(fittedPlot)
  
  #Posterior simulation ####
  #Samples from the posterior distribution of a GAM can be drawn using the simulate() methods from the 'gratia' package.
  
  set.seed(1) #set the random seed to make this reproducible
  nsim <- 20 #how many simulations to draw
  
  #Do the simulations
  sims <- simulate(mod, nsim = nsim, newdata = newIndep, unconditional = TRUE)
  
  #Rearrange the output into a long/tidy format
  colnames(sims) <- paste0("sim", seq_len(nsim))
  sims <- setNames(stack(as.data.frame(sims)), c("simulated", "run"))
  sims <- transform(sims, indep = rep(newIndep$indep, nsim), simulated = simulated)
  
  #Plot simulated trends
  simulationPlot <- ggplot(newIndep, aes(x = indep, y = fit)) +
    geom_line(data = sims, mapping = aes(y = simulated, x = indep, group = run),
              colour = "grey80") +
    geom_line(lwd = 1) +
    labs(y = depenVarLabel, x = indepVarLabel) +
    ylim(depenVarRange)
  print(simulationPlot)
  
  #Confidence and simultaneous intervals ####
  #Across-the-function and simultaneous confidence intervals computed using confint() method. The confint() methods return data frames suitable for plotting with 'ggplot'. Columns labelled est and se are esitimate values of hte smooth and its standar error. The variables lower and upper contain values of the lower and upper bounds for each interval.
  #What does shift doe? it seems to fix things...
  
  depenCInt <- confint(mod, parm = "indep", newdata = newIndep, type = "confidence", shift=TRUE)
  depenSInt <- confint(mod, parm = "indep", newdata = newIndep, type = "simultaneous", shift=TRUE)
  head(depenCInt)
  
  intervalPlot <- ggplot(depenCInt, aes(x = indep, y = est)) +
    geom_ribbon(data = depenSInt, mapping = aes(ymin = lower, ymax = upper, x = indep), fill = "grey80", inherit.aes = FALSE) +
    geom_ribbon(mapping = aes(ymin = lower, ymax = upper, x = indep), fill = "grey60", inherit.aes = FALSE) +
    geom_line(lwd = 1) +
    labs(y = depenVarLabel, x = indepVarLabel) +
    ylim(depenVarRange)
  print(intervalPlot)
  
  #Derivatives of the estimated trend ####
  #First derivative of the estimated trend calculated using finite differences using fderiv(). So first derivatives and a 95% simultaneous confidence interval for the trend computed and plotted:
  depenDerivative <- fderiv(mod, newdata = newIndep, n = N)
  depenSInt <- with(newIndep, cbind(confint(depenDerivative, nsim = nsim, type = "simultaneous"), indep = indep))
  
  derivativePlot <- ggplot(depenSInt, aes(x = indep, y = est)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "black") +
    geom_line() +
    labs(x = indepVarLabel, y = "First derivative")
  print(derivativePlot)
}
