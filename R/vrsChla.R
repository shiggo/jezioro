#'vrsChla
#'
#'@description
#'\code{vrsChla} infers chlorophyll a concentrations of sediments from spectral measurements of absorbance at wavelengths between 650-700 nm, following the approach described in Wolfe et al. (2006) and Michelutti et al. (2010).
#'
#'@usage vrsChla(d, profilePlot = TRUE)
#'
#'@param d matrix or data frame containing the input data
#'@param profilePlot (defaults to TRUE): logical to control whether to produce the plot of sediment depth vs. inferred chlorophyll a values
#'
#'@details The input data should be in the form of a matrix or data frame with 28 columns:
#'\itemize{
#'   \item Column 1: midpoint depth (cm) of the sediment interval
#'   \item Columns 2-27: values from the spectrophotometer for wavelengths 650-700 nm, measured every 2 nm
#'   }
#'NOTE: When using the Model 6500 series Rapid Content Analyzer at PEARL, the necessary values are contained in the 'spectra' tab of the excel file output (although they must be transposed). Ensure cells are formatted to 15 decimal places to avoid small rounding errors.
#'
#'@return
#'\code{chl.a.output} table of the inferred chlorophyll a values (mg/g dry mass) for each interval. 
#'
#'@author Joshua Thienpont, Adam Jeziorski
#'@export
#'
#'@examples 
#'#Infer the chlorophyll a values from the example spectral data
#'data(vrsChlaInput)
#'vrsChla(vrsChlaInput)
#'
#'@references Michelutti N, Blais JM, Cumming BF, Paterson AM, Ruhland K, Wolfe AP, Smol JP (2010) Do spectrally inferred determinations of chlorophyll a reflect trends in lake trophic status? Journal of Paleolimnology 43: 208-217
#'
#'Wolfe AP, Vinebrooke RD, Michelutti N, Rivard B, Das B (2006) Experimental calibration of lake-sediment spectral reflectance to chlorophyll a concentrations: methodology and paleolimnological validation. Journal of Paleolimnology 36: 91-100

vrsChla=function(d, profilePlot=TRUE){
chla.data=d
#Move the first column (interval labels) to row names and remove from table
rownames(chla.data)=chla.data[,1]
chla.data=chla.data[,-1]

numrows <-nrow(chla.data)
rows <- as.numeric(rownames(chla.data))
chla.work <- matrix(data=NA, nrow=numrows, ncol=52, dimnames=list(c(rows),c(650:700,"chla")))
chla.work[,1]=chla.data[,1]
chla.work[,3]=chla.data[,2]
chla.work[,5]=chla.data[,3]
chla.work[,7]=chla.data[,4]
chla.work[,9]=chla.data[,5]
chla.work[,11]=chla.data[,6]
chla.work[,13]=chla.data[,7]
chla.work[,15]=chla.data[,8]
chla.work[,17]=chla.data[,9]
chla.work[,19]=chla.data[,10]
chla.work[,21]=chla.data[,11]
chla.work[,23]=chla.data[,12]
chla.work[,25]=chla.data[,13]
chla.work[,27]=chla.data[,14]
chla.work[,29]=chla.data[,15]
chla.work[,31]=chla.data[,16]
chla.work[,33]=chla.data[,17]
chla.work[,35]=chla.data[,18]
chla.work[,37]=chla.data[,19]
chla.work[,39]=chla.data[,20]
chla.work[,41]=chla.data[,21]
chla.work[,43]=chla.data[,22]
chla.work[,45]=chla.data[,23]
chla.work[,47]=chla.data[,24]
chla.work[,49]=chla.data[,25]
chla.work[,51]=chla.data[,26]
chla.work[,52]=0
for (i in seq(2,51,by=2)){chla.work[,i]=(chla.work[,(i-1)]+chla.work[,(i+1)])/2}
for (i in 1:numrows){chla.work[i,52]=(0.0919*((sum(chla.work[i,]))-((chla.work[i,1]*51)+(((chla.work[i,51]-chla.work[i,1])*51)/2)))+0.0011)}

# For reference, the formulae that are used to fill in the final column (chla) are below:
# area.under.curve <- sum(chla.work[1,])
# area.under.line <- (chla.work[1,1]*51)+(((chla.work[1,51]-chla.work[1,1])*51)/2)
# peak.area <-area.under.curve-area.under.line
# chla.work[1,52] <- (0.0919*peak.area)+0.0011

if (profilePlot==TRUE) {
# Produces a figure of the chla over depth, scaled in a separate window
limit=(max(rows)+0.1*(max(rows)))

#Plot titles
main.title=expression(paste("VRS-inferred chlorophyll ", italic("a"), " profile"))
xaxis.title=expression(paste("Inferred chlorophyll ", italic("a"), " (mg"%.%"g"^{-1}, " dry mass)"))
yaxis.title="Core depth (cm)"

plot(chla.work[,52],rows, type="p", pch=21, bg="black", cex =0.75, lwd=1, main=main.title, xlab=xaxis.title, ylab=yaxis.title, ylim=c(limit,0), xlim=c(0,((max(chla.work[,52])+sd(chla.work[,52]))))) 
# adds a line at 0.01 mg/g, indicating the limit of detection for the method
abline(v=0.01, lwd=1.5, col="red")
text(0.15, limit,"limit of detection \n (0.01 mg/g)", cex=0.8, col="red")
}

chl.a.output <- matrix(data=NA,nrow=numrows, ncol=2, dimnames=list(c(rows),c("Core Depth(cm)","VRS-chla (mg/g dry wt)")))
chl.a.output[,1]=rows
chl.a.output[,2]=chla.work[,52]

# returns the 'chl.a.output' matrix
return(chl.a.output)
}
