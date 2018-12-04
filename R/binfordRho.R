#'binfordRho
#'
#'@description
#'\code{binfordRho} calculates the percent water and dry sediment density for the individual samples of a sediment core freeze-dried in preparation for gamma dating. The values for intervening (i.e. non-prepared samples) are interpolated.
#'
#'@usage binfordRho(d, bagwt=3.25, subsample=TRUE, coreD=7.62)
#'
#'@param d matrix or data frame containing the input data
#'@param bagwt (defaults to 3.25): average weight (g) of your brand/size of "Whirlpak" bags
#'@param subsample (Defaults to TRUE): logical to indicate whether the samples were subsampled before freeze-drying
#'@param coreD (defaults to 7.62):inner tube diameter (cm) of the sediment corer used to collect the core
#'
#'@details \code{binfordRho} interpolates the values of dry sediment density needed for dating a sediment core from the dry masses of the freeze-dried intervals and the bag weights of wet sediment (by calculating and interpolating percent water). The data to be examined should be in the form of a matrix or data frame with the following 7 columns:
#'\itemize{
#'  \item INTTOP: depth (cm) to the top of the analyzed interval
#'  \item INTBOT: depth (cm) to the bottom of the analyzed interval
#'  \item WT_SED+BAG: mass (g) of the wet sediment and the WhirlPak bag
#'  \item WT_VIAL: mass (g) of the empty vial used for freeze drying (disregarded if \code{subsample}=FALSE)
#'  \item WT_VIAL+SED: mass (g) of the vial after adding the subsample of wet sediment (disregarded if \env{subsample}=FALSE)
#'  \item WT_FD = mass (g) of the vial of sediment after freeze drying (if \code{subsample}=FALSE, this value should instead be the mass (g) of the dry sediment and the Whirlpak bag)
#'  \item WT_GAMMA: mass (g) of freeze-dried sediment added to the gamma tube
#'  }
#'
#'NOTE: Columns 1-3 must have values for the entire core length. 
#'Columns 4-7 will only have values for the intervals prepared for dating by freeze-drying, with the values for non-prepared intervals left blank.
#'If the entire "Whirlpak" bag was freeze-dried (i.e. subsample=FALSE), then Columns 4-5 should be left blank and the mass (g) of the bag of sediment after freeze drying entered into Column 6.
#'
#'@return
#'\code{binfordRhoOutput} table of results for use with the \code{binfordDates} function.
#'
#'@examples 
#'#Imports the example rho data for a sediment core subsampled prior to freeze-drying
#'data(binfordRhoInput)
#'binfordRho(binfordRhoInput)
#'@author Adam Jeziorski, Joshua Thienpont
#'@export

binfordRho=function(d, bagwt=3.25, subsample=TRUE, coreD=7.62){

# Input matrix is copied into a data matrix called "Weights"
Weights=d

# Creation of a data matrix called "wt" with 12 columns
wt=matrix(data=NA,nrow=nrow(Weights),ncol=12,dimnames=list(c(),c("INTTOP","INTBOT","WT_SED+BG","WT_SED","PER_WATER","WT_DRY","WT_CUMDRY","MASS_CUMDRY","CUM_TOP","CUM_BOT","RHO","WT_GAMMA")))
wt=as.data.frame(wt)
wt[,c(1:3,12)]=Weights[,c(1:3,7)]
wt[,4]=(Weights[,3]-bagwt)

# Calculation of "PER_WATER" for the values that went through the freeze drying process
# If subsample="TRUE" vial info is required, if subsample="FALSE" the vial columns are disregarded
if (subsample==TRUE) {wt[,5]=(1-((Weights[,6]-Weights[,4])/(Weights[,5]-Weights[,4])))*100}
if (subsample==FALSE) {wt[,5]=(1-((Weights[,6]-bagwt)/(Weights[,3]-bagwt)))*100}

# replace missing values in "PER_WATER" with a zero
wt$PER_WATER[is.na(wt$PER_WATER)]=0

# Check to determine presence of intervals below the last freeze-dried interval, the value of "empty" is the number of extra rows
emptytest=function(d){if(d[1]==0){(which(d>0)[1])-1} else(0)}
empty=emptytest(rev(wt[,5]))

# Removes the bottom 'empty' rows of 'wt
wt=head(wt,nrow(wt)-empty)

# Scan 'wt' to deterimine the max # of consecutive intervals in "5%WATER" with a missing value
# A vector "emp" that is nrow(wt) intervals long is created and populated with zeros
emp=rep(0, times = nrow(wt))

# Counts through the %WATER values "k" is equal to the number of consecutive zeroes at i and places the value in "emp"
# "k" resets when it sees a non-zero value
for (i in 1:nrow(wt)){
 if (wt[i,5]>0){k=0}
  else if (wt[(i-1),5]==0) {k=k+1}
  else {k=k+1}
 {emp[i]=k}
}

# Interpolates "PER_WATER" values for intervals that were not freeze dried 
# The maximum value of "emp" controls the number of subdivisions between measured intervals
for (i in 1:nrow(wt)) {
 n=1
 if (wt[i,5]>0){next}
  else  for (j in 1:max(emp)) {if ((wt[(i+j),5])==0) {n=n+1}
  else break
 }
 for (k in 1:n) {wt [(i+k-1),5]= (wt[(i-1),5])+(((wt[(i+n),5]-(wt[(i-1),5]))/(n+1))*k)
 }
}

# Fills in the dry weight (g) of each interval "6WT_DRY" using the calculated "5%WATER"
wt[,6] = (100-(wt[,5]))/100*(wt[,4])

# Fills in the cumulative dry weight "7WT_CUMDRY" going downcore
wt[1,7]=wt[1,6]
for (i in 2:nrow(wt)) {
 wt[i,7] = wt[(i-1),7]+wt[i,6]
}

# Fills in the cumulative dry mass "8MASS_CUMDRY" going downcore assuming an inner tube diameter of "coreD" in cm (change this variable accordingly)
wt[,8]=wt[,7]/(pi*(coreD/2)^2)

# Fills in the cumulative dry mass to the top of the interval "9CUM_TOP"
wt[1,9]=0
for (i in 2:nrow(wt)) {
 wt[i,9] = wt[(i-1),8]
}

# Fills in the cumulative dry mass to the bottom of the interval "10CUM_BOT" (this is just a duplicate of "8MASS_CUMDRY" renamed for clarity)
wt[,10]=wt[,8]

# Fills in the value "11RHO" (the density) of the sediment interval (note it is dependent on the interval thickness)
wt[,11]=((wt[,10])-(wt[,9]))/(wt[,2]-wt[,1])

# The matrix "ex" is a summary table of the results of the calculations
ex=cbind(wt[,(1:2)], wt[,(9:11)], Weights[(1:nrow(wt)),7],wt[,5])
colnames(ex)[6]=("WT_GAMMA")
colnames(ex)[7]=("%WATER")

# 'binfordRhoOutput'=subset of 'ex' for which there is a gamma tube measurement (i.e. ex[,6] has a value)
binfordRhoOutput=subset(ex, is.na(ex[,6])==F)

# returns the "binfordRhoOutput" dataframe for eventual use by the "binforDates" function
return(binfordRhoOutput)
}