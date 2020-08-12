#'interpDates
#'
#'@description
#'\code{interpDates} interpolates dates for any undated intervals in a sediment core from the midpoint and age of the dated intervals, as well as the sectioning resolution.
#'
#'Three sets of interpolated dates are returned:
#'
#'\itemize{
#'   \item connectTheDotsDates: are calculated from straight line between sequential date pairs
#'   \item linearDates: are calculated by fitting a straight line through all dated intervals
#'   \item polynomialDates: are calculated by fitting a 2nd order polynomial line through all dated intervals
#'   }
#'   
#'For more involved approaches to the estimation of age-depth relationships consult Blaauw and Heegaard (2012).
#'
#'@usage interpDates(d, intervalWidth=0.5)
#'
#'@param d matrix or data frame containing the input data
#'@param intervalWidth (defaults to a constant value of 0.5): a single numeric values is used to indicate a constant sectioning resolution, but a vector (with interval widths for the full core) can also be used to provide the width of each interval in the core
#'
#'@details the input data from the dated sediment intervals should be arranged in 2 columns:
#'\itemize{
#'   \item Column 1: midpoint depth of the dated interval
#'   \item Column 2: date of the interval
#'   }
#'
#'@return core, the table containing the interpolated dates for each interval.
#'
#'@author Adam Jeziorski
#'@export
#'
#'@references Blaauw M, Heegaard E (2012) 12. Estimation of age-depth relationships. In: Birks HJB, Lotter AF, Juggins S, Smol JP (eds.) Tracking Environmental Change Using Lake Sediments. Volume 5: Data Handling and Numerical Techniques. Springer, Netherlands, pp 379-413
#'
#'@examples
#'#4 dated intervals with a constant sectioning resolution of 0.5cm
#'interpDates.input <- cbind(c(0.25, 4.25, 8.25, 18.25), c(2017, 2000, 1950, 1850))
#'interpDates(interpDates.input)
#'
#'#4 dated intervals in a core with slices 1-10 = 0.5cm, and 11-20 = 1.0cm
#'interpDates.input <- cbind(c(0.25, 4.25, 8.5, 13.5), c(2017, 1980, 1920, 1850))
#'intervalWidth <- c(rep(0.5,10), rep(1.0,10))
#'interpDates(interpDates.input, intervalWidth)

interpDates <- function (d, intervalWidth=0.5){
  #Ensure d is a data.frame
  d <- as.data.frame(d)
  colnames(d) <- c("intervalMidpt", "intervalAge")
  
  #Remove intervals without a CRS_Age before building the model
  if(any(is.na(d[,2]))==TRUE) (d <- d[-which(is.na(d[,2])==TRUE),])
  
  #Create a sequence containing each interval from top to bottom 
  #When intervalWidth is a single value
  if (length(intervalWidth)==1 ){
    intervalTop <- seq(from=0, to=(d[nrow(d),1] - (intervalWidth/2)), by=intervalWidth)
    intervalBot <- intervalTop+intervalWidth
    intervalMidpt <- (intervalTop+intervalBot)/2
  } else {
    #When intervalWidth is a vector
    intervalBot <- cumsum(intervalWidth)
    intervalTop <- intervalBot-intervalWidth
    intervalMidpt <- (intervalTop+intervalBot)/2
  }
  
  core <- cbind(intervalTop, intervalBot, intervalMidpt)
  
  #Merge 'd' onto 'core' according to 'intervalMidpt' while keeping all intervals
  core <- merge(core, d, by="intervalMidpt", all.x=TRUE)
  
  #Identify row positions of dated intervals
  measuredDatesPositions <- which(is.na(core$intervalAge)==FALSE)
  measuredDatesValues <- core$intervalAge[measuredDatesPositions]
  #Remove any rows below the last dated interval
  core <- core[(1:max(measuredDatesPositions)),]
  #Identify the number of emptyCells or dates to interpolate
  emptyCells <- diff(measuredDatesPositions)-1
  
  #Create matrix with columns first date, second date, # of intervening undated intervals
  dateMatrix <- cbind(measuredDatesValues, c(measuredDatesValues[2:length(measuredDatesValues)],NA), c(emptyCells,NA))
  colnames(dateMatrix)=c("firstDate","secondDate","emptyCells")
  
  #Function interpolates between two dates over 'emptyCells' steps, then removes last value to avoid duplicates
  dateInterpolation <- function(e){seq(e[1],e[2], length.out=(e[3]+2))[-length(seq(e[1],e[2], length.out=(e[3]+2)))]}
  
  #Apply 'dateInterplation' to 'dateMatrix' by rows (excluding the last row), unlist to turn into a vector and append the final date
  connectTheDotsDates <- apply(dateMatrix[(1:nrow(dateMatrix)-1),], 1, dateInterpolation)
  connectTheDotsDates <- unlist(connectTheDotsDates)
  connectTheDotsDates <- c(connectTheDotsDates, as.numeric(dateMatrix[nrow(dateMatrix),1]))
  
  #Append 'connectTheDotsDates to 'core'
  core <- cbind(core, connectTheDotsDates)
  
  #Intial plot of age vs. depth with the dated samples
  plot(d$intervalMidpt, d$intervalAge, xlab = "Midpoint Depth (cm)", ylab = "Sediment Age (Year)")
  
  #Perform a first-order or linear regression 'fit1' on intervalAge vs. intervalMidpt
  fit1 <- lm(d$intervalAge ~ d$intervalMidpt)
  summary(fit1)
  abline(fit1, col = "blue")
  
  #Perform a second-order polynomial regression 'fit2' on intervalAge vs. intervalMidpt
  fit2 <- lm(d$intervalAge ~ poly(d$intervalMidpt, 2, raw=TRUE))
  summary(fit2)
  points(d$intervalMidpt, predict(fit2), type="l", col="red", lwd=1)
  
  #Pull the adjusted R^2 values from both models and use the coefficients to add equations of both lines to the plot
  fit1.adj.r.squared <- summary(fit1)$adj.r.squared
  fit2.adj.r.squared <- summary(fit2)$adj.r.squared
  
  #fit1
  text(0.55*(max(d$intervalMidpt)), 0.945*(max(d$intervalAge)), paste("y = ", round(fit1$coefficient[2],2),"x + ", round(fit1$coefficient[1],2)), col="blue", pos = 2, cex=.8)
  text(0.55*(max(d$intervalMidpt)), 0.940*(max(d$intervalAge)), paste("Adjusted R^2 =", round(fit1.adj.r.squared, 3)), col="blue", pos = 2, cex=.8)
  
  #fit2
  text(0.31*(max(d$intervalMidpt)), 0.995*(max(d$intervalAge)), paste("y = ", round(fit2$coefficient[3],2) ,"x^2 + " , round(fit2$coefficient[2],2),"x + ", round(fit2$coefficient[1],2)), col="red", pos = 4, cex=.8)
  text(0.31*(max(d$intervalMidpt)), 0.990*(max(d$intervalAge)), paste("Adjusted R^2 =", round(fit2.adj.r.squared, 3)), col="red", pos = 4, cex=.8)
  
  #Calculate ages for all the 'intervalMidpt' values in 'core' using 'fit1'
  lin1 <- function(x) fit1$coefficient[2]*x + fit1$coefficient[1]
  linearDates <- lin1(core$intervalMidpt)
  core <- cbind(core, linearDates)
  
  #Calculate ages for all the 'intervalMidpt' values in 'core' using 'fit2'
  pol2 <- function(x) fit2$coefficient[3]*x^2 + fit2$coefficient[2]*x + fit2$coefficient[1]
  polynomialDates <- pol2(core$intervalMidpt)
  core <- cbind(core, polynomialDates)
  
  return(core)
}
