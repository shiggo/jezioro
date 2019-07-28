#'dipteranVWHO
#'
#'@description
#'\code{dipteranVWHO} is a function to harmonize appropriately formatted chironomid count data with the calibration set described in Quinlan and Smol (2001, 2010), construct the described VWHO model, use the model to infer VWHO from the chironomid sample data, and finally perform an analog matching analysis between the sample data and the calibration set (that follows the procedures described in the 'analogue_methods' vignette provided by the 'analogue' package). 
#'
#'@usage dipteranVWHO(d, evaluate=FALSE, percentileCut=5, lowCount=50)
#'
#'@param d matrix or data frame containing the input data
#'@param evaluate (defaults to FALSE): logical to indicate whether to display the model information
#'@param percentileCut (defaults to 5): cutoff value for flagging samples in the plots as having poor modern analogs
#'@param lowCount (defaults to 50): cutoff value to flag samples that have low counts (i.e. less than 50 individuals) 
#'
#'@details The required format of the input data is provided in 'dipteranVWHOInput'
#'
#'@return \code{output} list  of 4 sublists 
#' #'\itemize{
#'   \item formattedData contains the sample info, the raw chironomid counts, the raw chaoborid counts, the chironomid counts aggregated by taxon code to match the taxonomy used in vwhoQuinlan 2010 dataset, and vector of the 44 taxon codes used in the model.
#'   \item vwhoModel the model object constructed by the 'WA' function provided by 'rioja'
#'   \item vwhoResults contains model ouput, relative abundance of the dipterans used in the model, relative abundances of all dipterans in the dataset
#'   \item analogResults contains the output from the analog analyis, the close modern analogs, and the minimum dissimilarity values
#'   }
#' Three plots are produced: minimum dissimilarity between each sample interval and the calibration set data, sediment depth vs dipteran-inferred VWHO, and a composite combining data from both graphs.
#'
#'@examples 
#'#Import the example input data for use with dipteranVWHO
#'data(dipteranVWHOInput)
#'dipteranVWHO(dipteranVWHOInput, evaluate=TRUE, percentileCut=5, lowCount=50)
#' 
#'@author Adam Jeziorski
#' 
#'@references Quinlan R, Smol JP (2001) Chironomid-based inference models for estimating end-of-summer hypolimnetic oxygen from south-central Ontario shield lakes. Freshwater Biology 46: 1529-1551
#'
#'Quinlan R, Smol JP (2010) Use of Chaoborus subfossil mandibles in models for inferring past hypolimnetic oxygen. Journal of Paleolimnology 44: 43-50
#'@export

dipteranVWHO <- function (d, evaluate=FALSE, percentileCut=5, lowCount=50){
##Save initial 'par' settings and restore on function exit
old.par <- par(no.readonly = TRUE)
on.exit(par(old.par))

##Format the input data ####
formatInput <- function (d){
rawData <- d
##Separate sample info from actual chironomid/chaoborid count data
sampleInfo <- rawData[which(rawData[,1] == "SampleInfo"),]
colnames(sampleInfo) <- sampleInfo[which(sampleInfo[,3] == "Code"),]
chironomidData <- rawData[-which(rawData[,1] == "SampleInfo"),]
colnames(chironomidData) <- sampleInfo[which(sampleInfo[,3] == "Code"),]
##Replace NA values in count data with zero, ensure counts are numeric data
chironomidData[is.na(chironomidData)] <- 0
chironomidData[,4:ncol(chironomidData)] <- sapply(chironomidData[,4:ncol(chironomidData)], as.numeric)
  
##Isolate the chironomid and chaoborid count data
chaoboridData <- chironomidData [which(chironomidData$SampleInfo == "Chaob"),]
chironomidData <- chironomidData [-which(chironomidData$SampleInfo == "Chaob"),]
##Remove Bytho row if present
if (length(which(chironomidData$SampleInfo == "Bytho")) > 0) {
  chironomidData <- chironomidData [-which(chironomidData$SampleInfo == "Bytho"),]
}
  
#Remove first column from each data frame (contains tribe abbreviation)
chironomidData <- chironomidData[,-1]
chaoboridData <- chaoboridData[,-1]
sampleInfo <- sampleInfo[,-1]
  
##Aggregate chironomid data by species code (for use with VWHO model)
chironomidDataAggregated <- aggregate(as.matrix(chironomidData[,3:ncol(chironomidData)])~chironomidData$Code, FUN=sum)
rownames(chironomidDataAggregated) <- chironomidDataAggregated[,1]
chironomidDataAggregated <- chironomidDataAggregated[,-1]
##Vector containing species codes used in the VWHO model
modelTaxa <- c("TANYSL", "TANYCH", "TANYLU", "CLDTYA", "CLDTYM", "MICROP", "STMPLN", "STMPLL", "PSEDCH", "CHIRON", "CHIRS1", "CLADOP", "CRYPCH", "CYPHOM", "DICROT", "ENDOCH", "GLYPTO", "LAUTER", "MICROT", "PGASTL", "PARACH", "PARALA", "PARATE", "POLYPE", "SERGEN", "STICTO", "PROCLD", "PENTAN", "NILOTA", "PROTAN", "CORYTH", "CRICOR", "EUKIEF", "HTRTRS", "LIMNOP", "NANOCL", "PARAKB", "PARAKA", "PSECTM", "PSECTP", "PSDSMA", "SYNORT", "ZALUTZ", "ZALUTS")
##Identify whether any model taxa are absent from aggregated data set and if so append appropriate rows populated by zeroes 
missingTaxa <- setdiff(modelTaxa, rownames(chironomidDataAggregated))
if (length(missingTaxa) > 0) {
  chironomidDataAggregated <- rbind(chironomidDataAggregated, matrix(0, length(missingTaxa), ncol(chironomidDataAggregated), dimnames=list(missingTaxa, colnames(chironomidDataAggregated))))
}
##Reorder aggregated chironomid taxa to match required input order of VWHO model
chironomidDataAggregated <- chironomidDataAggregated[order(match(rownames(chironomidDataAggregated), modelTaxa)),]
  
##Remove species codes from chironomidData as they are no longer needed
rownames(chironomidData) <- chironomidData[,1]
chironomidData <- chironomidData[,-(1:2)]
  
##Append total number of mandibles to chaoboridData
rownames(chaoboridData) <- chaoboridData$Code
chaoboridData <- chaoboridData[,-(1:2)]
chaoboridData <- rbind(chaoboridData, apply(chaoboridData, MARGIN=2, FUN=sum))
rownames(chaoboridData)[nrow(chaoboridData)] <- "C.TOT"
  
##Tidy up the sampleInfo (and ensure remaining rows are numeric)
rownames(sampleInfo) <- sampleInfo[,2]
sampleInfo <- sampleInfo[,-(1:2)]
sampleInfo <- sampleInfo[-which(rownames(sampleInfo) == "interval"),]
sampleInfo <- sampleInfo[-which(rownames(sampleInfo) == "Code"),]
sampleInfo <- sampleInfo[-which(rownames(sampleInfo) == "sampleType"),]
sampleInfo <- sampleInfo[-which(rownames(sampleInfo) == "complete"),]
  
temp <- rownames(sampleInfo)
sampleInfo <- sapply(sampleInfo[,1:ncol(sampleInfo)], as.numeric)
rownames(sampleInfo) <- temp
  
##Append chironomid head capsule density (#/g drysed) to sampleInfo using sediment weights and total counts
sampleInfo <- rbind(sampleInfo, (apply(chironomidData, MARGIN=2, FUN=sum)/2))
rownames(sampleInfo)[nrow(sampleInfo)] <- "chirIndividuals"
sampleInfo <- rbind(sampleInfo, sampleInfo["chirIndividuals",]/((sampleInfo["wetSed",]*(1-sampleInfo["percWater",])) + sampleInfo["drySed",]))
rownames(sampleInfo)[nrow(sampleInfo)] <- "density_chirInd_drySed"
  
##Prepare list containing formatted data
formattedData <- list("sampleInfo" =sampleInfo, "chironomidCount" = chironomidData, "chaoboridCount" = chaoboridData, "chironomidCountAggregated" = chironomidDataAggregated, "modelTaxa" = modelTaxa)
  
return(formattedData)
}
formattedData <- formatInput(d)

##Build the Quinlan and Smol 2010 VWHO model ####
buildVWHOModel <- function (evaluate){
##'vwhoQuinlan2010' in sysdata.rda contains the VWHO calibration set of 54 sites and %RelAbund of identifiable Diptera
##Inference taxa: 44 chironomids, and total chaoborids

##Use the 'WA' function provided by 'rioja' to build the VWHO model
##The assemblage data (columns 3:47) are use to model the VWHO data (column 2)
vwhoModel <- WA(vwhoQuinlan2010[,3:47], vwhoQuinlan2010[,2], tolDW=TRUE)
  
##Evaluation of the model using 'crossval'
if (evaluate==TRUE){
  ##Compare metrics those given in the Quinlan and Smol papers 
  ##WA.inv.tol_XVal model: RMSE=1.98 mg/L, R2=0.5998
  modelPerformance <- rioja::crossval(vwhoModel, cv.method="loo", verbose=TRUE, ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL)
  print("The 'best' (lowest RMSEP, highest r^2) VWHO inference model described in Quinlan and Smol 2010 (JOPL) was a chironomid and chaoborid model (54-lake training set, weighted-averaging with tolerance downweighting, inverse deshrinking). The reported performance metrics were RMSEP=1.98 and r^2(jack)=0.60.")
  print("The performance metrics of the current model are:")
  print(rioja::performance(modelPerformance)$crossval[3,(1:2)])
}
return(vwhoModel)
}
vwhoModel <- buildVWHOModel(evaluate)

##Infer VWHO from sample data using dipteran relative abundances and 'vwhoModel' ####
inferVWHO <- function (formattedData, vwhoModel){
##Simple functions
percentage=function(d){d/sum(d)*100}
chaobChir <- function(d){d[length(d)]/sum(d)}

dipteranData <- rbind(formattedData$chironomidCountAggregated, formattedData$chaoboridCount["C.TOT",])
##Convert dipteran count datas to number of individuals represented
dipteranData <- dipteranData/2
##Initialize VWHO results table
vwhoResults <- cbind(formattedData$sampleInfo["midpt",], formattedData$sampleInfo["age",], apply(dipteranData, MARGIN=2, FUN=chaobChir))
colnames(vwhoResults) <- c("midpt", "age", "chaobChirRatio")
  
##Convert dipteran data to relative abundances
dipteranData <- apply(X=dipteranData, MARGIN=2, FUN=percentage)
##Export the dipteran data to allow for external analyses
dipteranRelAbAll <- dipteranData

##Remove all the taxa not included in the model calibration set
dipteranData=dipteranData[c(1:44,nrow(dipteranData)),]

## % of sample as inference taxa
vwhoResults <- cbind(vwhoResults, apply(dipteranData, 2, sum))
colnames(vwhoResults)[ncol(vwhoResults)] <- c("propInfTaxa")
  
#Predict VWHO using the Quinlan and Smol 2010 model
#Note there is an option to calculate sample-specific errors using bootstrapping
lakeVWHO <- predict(vwhoModel, t(dipteranData))
lakeVWHO
  
vwhoResults <- cbind(vwhoResults, formattedData$sampleInfo["chirIndividuals",], formattedData$sampleInfo["density_chirInd_drySed",], lakeVWHO$fit[,3])
colnames(vwhoResults)="VWHO.mg/L" <- c("midpt", "age", "chaobChirRatio", "propInfTaxa", "chirTotal", "chirDensity/g.Dry", "VWHO.mg/L")

vwhoResults <- list("vwhoResults" = vwhoResults, "dipteranPercModel" = dipteranData, "dipteranPercAll" = dipteranRelAbAll)

return(vwhoResults)
}
vwhoResults <- inferVWHO(formattedData, vwhoModel)

##Analog Matching ####
matchAnalogs <- function(vwhoResults, percentileCut){
  ##Format model and test data for use with the analogue package
  calibrationSet <- vwhoQuinlan2010
  row.names(calibrationSet) <- calibrationSet[,1]
  calibrationSet <- calibrationSet[,-1]
  calibrationSet <- data.frame(calibrationSet)
  calibrationSet_vwho <- calibrationSet[,1]
  names(calibrationSet_vwho) <- row.names(calibrationSet)
  calibrationSet <- calibrationSet[,-1]
  calibrationNames <- rownames(calibrationSet)
  
  ##The testSet initially contains all taxa in the samples
  testSet <- t(vwhoResults$dipteranPercAll)
  testSet <- data.frame(testSet)
  testNames <- rownames(testSet)
  
  ##Ensure calibrationSet and testSet have all the same taxa in the same order
  ##Note that this uses the plyr version of 'join'
  calibrationSet <- suppressMessages(plyr::join(calibrationSet,testSet, type="left"))
  calibrationSet[is.na(calibrationSet)] <- 0
  calibrationSet <- calibrationSet[,order(names(calibrationSet))]
  rownames(calibrationSet) <- calibrationNames
  
  testSet <- suppressMessages(plyr::join(testSet, calibrationSet, type="left"))
  testSet[is.na(testSet)] <- 0 # fill in empty columns with 0's
  testSet <- testSet[,order(names(testSet))] # order columns alphabetically 
  rownames(testSet) <- testNames
  
  ##Finally convert both data sets from percentages to proportions 
  calibrationSet <- calibrationSet/100
  testSet <- testSet/100
  
  ##Perform analog matching between calibrationSet and testSet using Bray-Curtis dissimilarity coefficients (requires "analogue" package)
  results.analog <- analog(calibrationSet, testSet, method=c("bray"))
  
  ##Identify close modern analogues in calibrationSet for each testSet sample
  ##i.e. samples in calibrationSet as close or closer than some critical threshold of dissimilarity (percentileCut) to each testSet sample
  ##percentileCut defaults to 5th percentile, adjusted via argument to be more strict or loose (1=1%, 2=2%, 3=5%, 4=10%, 5=20%)
  if (any(c(1,2,5,10,20) == percentileCut)) {
    if (percentileCut==1){percentileCut <- 1}
    if (percentileCut==2){percentileCut <- 2}
    if (percentileCut==5){percentileCut <- 3}
    if (percentileCut==10){percentileCut <- 4}
    if (percentileCut==20){percentileCut <- 5}
  } else {
    stop("percentileCut requires a value of 1, 2, 5, 10 or 20")
  }
  results.cma=cma(results.analog, cutoff= summary(results.analog)[[3]][percentileCut])
  
  #List of samples in the training set that are close modern analogues to each test interval given in results.cma$close
  
  #Isolate minimum dissimilarity with calibrationSet for each testSet interval
  results.minDC <- minDC(results.analog)
  minDC <- results.minDC$minDC
  
  analogResults <- list("analogResults" = results.analog, "cmaResults" = results.cma, "minDCResults" = results.minDC)
  return(analogResults)
}
analogResults <- matchAnalogs(vwhoResults, percentileCut)

##Create Analogs Plot ####
plotAnalogs <- function(analogResults, depthAxis=TRUE){
##Check pair-wise distribution of minimum dissimilarity coefficents to determine where the cutoffs will be set 
dissimDistrib <- dissim(analogResults$analogResults, which=c("train"))
#plot(dissimDistrib, xlim=c(0,1), prob=c(0.05, 0.1, 0.9, 0.95))
dissimDistrib <- data.frame(matrix(unlist(dissimDistrib)))
d <- density(as.numeric(dissimDistrib[,1]))

##Plot cumulative probablity of dissimiarity coefficients
if(depthAxis==TRUE){par(mar = c(3.5,5,0.5,5))}
else{par(mar = c(3.5,4,0.5,1.25))}

plot(d, xlim=c(0,1), main="", xlab="", ylab="", las=1)
title(xlab="Dissimilarity (Bray-Curtis)", line=2.4, cex.lab=1.3)
title(ylab="Density", line=2.5, cex.lab=1.3)
##Add lines for 5th and 10th percentiles (pulled from results.analog)
abline(v=summary(analogResults$analogResults)[[3]][3], col="red",lty=2)
text((summary(analogResults$analogResults)[[3]][3]-0.02),0.1,"5th", col="red", cex=0.6)
abline(v=summary(analogResults$analogResults)[[3]][4], col="red",lty=2)
text((summary(analogResults$analogResults)[[3]][4]-0.02),0.1,"10th", col="red", cex=0.6)
  
##Pull in midpt depth from 'formattedData' object
depth <- as.numeric(formattedData$sampleInfo[1,])
depthRange=rev(range(depth))
depthRange=c(ceiling(max(depth)/10)*10,0) #round depth up to nearest 10cm
  
##Add minDC vs depth to the existing plot
par(new=TRUE)
plot(analogResults$minDCResults$minDC, depth, xaxt="n", yaxt="n", xlim=c(0,1), ylim=depthRange, xlab=NA, ylab=NA, type="p", pch=19)
axis(side=4, las=1)
if (depthAxis==TRUE){
mtext(side =4, line=2.5, "Sediment Depth (cm)", cex=1.3)
}
legendtext=c("Min. Dissimilarity w Calibration Set", "Pairwise Dist. of Calibration Set MinDC", "Percentile Cutoffs")
legend(0.6, 0.01, inset=c(.01,0), legendtext, pch=c(19,NA_integer_,NA_integer_), lty=c(0,1,2),col=c("black","black","red"), pt.cex=1, cex=0.7, bty="n", y.intersp=0.8, x.intersp=0.55, seg.len=1)
}
plotAnalogs(analogResults, depthAxis=TRUE)

##Create VWHO Plot####
plotVWHO <- function(formattedData, vwhoResults, analogResults, lowCount){
##Samples with low counts (i.e. < 50 head capsules) given different shape
lowCounts <- vwhoResults$vwhoResults[,5]<lowCount
lowCounts[which(lowCounts == TRUE)]<- 25 #filled triangles
lowCounts[which(lowCounts==FALSE)] <- 21 #filled circle
LegendCategory1=eval(substitute(paste(">", nn, "Head Capsules"), list(nn=lowCount)))
LegendCategory2=eval(substitute(paste("<", nn, "Head Capsules"), list(nn=lowCount)))
#Samples with no close modern analogs in the calibration set plotted grey
noCloseModernAnalogs <- analogResults$cmaResults$n.analogs==0
noCloseModernAnalogs[which(noCloseModernAnalogs == FALSE)] <- "black"
noCloseModernAnalogs[which(noCloseModernAnalogs == TRUE)] <- "grey"

depth <- as.numeric(formattedData$sampleInfo[1,])
depthRange=rev(range(depth))
depthRange=c(ceiling(max(depth)/10)*10,0) #round depth up to nearest 10cm

par(mar = c(3.5,5,0.5,5))
plot(vwhoResults$vwhoResults[,7], vwhoResults$vwhoResults[,1],
     type="p", #point only
     xlab="",
     ylab="",
     xlim=c(floor(min(vwhoResults$vwhoResults[,7])), ceiling(max(vwhoResults$vwhoResults[,7]))),
     ylim=depthRange,
     las=1, #rotate y-axis labels
     pch=lowCounts, #triangle indicates <40 head capsules
     cex=1.5,
     col="black",
     bg=noCloseModernAnalogs #
)
title(xlab="VWHO (mg/L)", line=2.4, cex.lab=1.3)
title(ylab="Sediment Depth (cm)", line=2.5, cex.lab=1.3)
legend("bottomleft", inset=c(.01,0), legend=c(LegendCategory1, LegendCategory2, "Poor Modern Analogs"), pch=c(21,25,21), bty="n", col="black", pt.cex=0.8, pt.bg=c("black", "black", "grey"), cex=0.7, y.intersp=0.8, x.intersp=0.95)

##Add Age Axis
par(new=TRUE)
mtext(side =4, line=3.5, expression(" "^{210}~Pb~Age), cex=1.3)
axis(side = 4, at=vwhoResults$vwhoResults[,1], round(vwhoResults$vwhoResults[,2], digits=0), las=1, cex.axis=0.65)
}
plotVWHO(formattedData, vwhoResults, analogResults, lowCount)

##Composite Plot ####
##Save both plots as functions
plot1 <- function(){plotAnalogs(analogResults, depthAxis=FALSE)}
plot2 <- function(){plotVWHO(formattedData, vwhoResults, analogResults, lowCount)}

layout(matrix(c(1,2), nrow= 1), widths=c(1.25,1))
plot1()
plot2()

##Export inference data for external analyses ####
return(vwhoOutput <- list("formattedData" = formattedData, "vwhoModel" = vwhoModel,"vwhoResults" = vwhoResults, "analogResults" = analogResults))
}
