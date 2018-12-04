#'binfordActivity
#'
#'@description
#'\code{binfordActivity} uses the appropriate Region of Interest (ROI) values from the gamma counter to correct the 210Pb, 214Bi and 137Cs activities for the amount of sediment in the gamma tube as well as the background measurements of the gamma counter.
#'
#'@usage binfordActivity(d, LakeName = "Lake Name", live = 80000, PbBk = 0, BiBk = 0, 
#'CsBk = 0, PbBkErr = 0, BiBkErr = 0, CsBkErr = 0, graph.units="dpm")
#'
#'@param d the matrix or data frame containing the input data.
#'@param LakeName Name of the lake being analyzed. Defaults to "Lake Name".
#'@param live Time in seconds that the samples were counted for. Defaults to 80000.
#'@param PbBk Background 210Pb (cpm). Defaults to 0.
#'@param BiBk Background 214Bi (cpm). Defaults to 0.
#'@param CsBk Background 137Cs (cpm). Defaults to 0.
#'@param PbBkErr Background Error for 210Pb (cpm). Defaults to 0.
#'@param BiBkErr Background Error for 214Bi (cpm). Defaults to 0.
#'@param CsBkErr Background Error for 137Cs (cpm). Defaults to 0.
#'@param graph.units determines whether activities should be plotted as dpm/g ("dpm") or Bq/g ("bq"). Defaults to "dpm".
#'
#'@details \code{binfordActivity} calculates the activity of the 210Pb, 137Cs, and 214Bi isotopes, corrected for the efficiency of the gamma counter device, and the time between sediment core collection and analysis. A plot is generated of the corrected isotope activities vs. sediment core depth. The data to be examined should be in the form of a matrix or data frame with the following 14 columns:
#'\itemize{
#'  \item INTTOP: top of the analyzed interval (cm)
#'  \item INTBOT: bottom of the analyzed interval (cm)
#'  \item DATECORR: number of days between core collection and analysis of the interval
#'  \item ROI1: sum of the activity of the 5 channel region of interest immediately less than the 210Pb photopeak (ROI1 from the Maestro ROI Report)
#'  \item ROI2: sum of the activity of the 10 channel region of interest centered on the 210Pb photopeak (45.5keV; ROI2 from the Maestro ROI Report)
#'  \item ROI3: sum of the activity of the 5 channel region of interest immediately greater than the 210Pb photopeak (ROI3 from the Maestro ROI Report)
#'  \item ROI4: sum of the activity of the 5 channel region of interest immediately less than the 214Bi photopeak (ROI10 from the Maestro ROI Report)
#'  \item ROI5: sum of the activity of the 10 channel region of interest centered on the 214Bi photopeak (609keV; ROI11 from the Maestro ROI Report)
#'  \item ROI6: sum of the activity of the 5 channel region of interest immediately greater than the 214Bi photopeak (ROI12 from the Maestro ROI Report)
#'  \item ROI7: sum of the activity of the 5 channel region of interest immediately less than the 137Cs photopeak (ROI13 from the Maestro ROI Report)
#'  \item ROI8: sum of the activity of the 10 channel region of interest centered on the 137Cs photopeak (662keV; ROI14 from the Maestro ROI Report)
#'  \item ROI9: sum of the activity of the 5 channel region of interest immediately greater than the 137Cs photopeak (ROI15 from the Maestro ROI Report)
#'  \item WTinTube: mass (g) of dried sediment added to the gamma tube (equivalent to WT_GAMMA in \code{binfordRho})
#'  \item HTinTube: height (mm) of the sediment added to the gamma tube, generally measured to 3 decimal places using an accurate set of calipers
#'  }
#'
#'NOTE: Columns 1 and 2 must have values for all samples prepared for gamma analysis (the same values entered into Columns 1 and 2 of \code{binfordRho}), but no values for columns 3-9 in cases where samples have not yet been run or are not required (in cases where background has been reached).
#'
#'\code{binfordActivity} is intended to be run in an iterative function following analysis of each sample, thereby allowing accurate and timely determination of when supported levels of 210Pb have been reached, in order to avoid the superfluous analysis of samples beyond background.
#'
#'It is very important to use the correct background and error measurements (e.g. PbBk, PbBkErr, etc.). That is, the values specific to both the gamma counter used and the time period when the samples were measured.
#'
#'@return
#'\code{binfordActivityOutput} table of results for use with the \code{binfordDates} function.
#'
#'@examples 
#'#Imports the example activity data for a sediment core subsampled prior to freeze-drying
#'data(binfordActivityInput)
#'binfordActivity(binfordActivityInput)
#'@author Adam Jeziorski, Joshua Thienpont
#'@references Schelske CL, Peplow A, Brenner M, Spencer CN (1994) Low-background gamma counting: applications for 210Pb dating of sediments Journal of Paleolimnology 10: 115-128
#'@export

binfordActivity=function(d, LakeName="Lake Name", live=80000, PbBk=0, BiBk= 0, CsBk=0, PbBkErr=0, BiBkErr=0, CsBkErr=0, graph.units="dpm"){
  roi=d
  
  # First remove intervals without any data  
  # Replace any missing values in "DATECORR" with a zero
  roi$DATECORR[is.na(roi$DATECORR)]=0
  
  # The value of "empty" is the number of extra rows in 'roi' to be removed 
  emptytest=function(d){if(d[1]==0){(which(d>0)[1])-1} else(0)}
  empty=emptytest(rev(roi[,3]))
  roi=head(roi,nrow(roi)-empty)
   
  # Creation of a data matrix called "cor" with 18 columns, this matrix is populated by the script
  cor=matrix(data=NA,nrow=nrow(roi),ncol=18,dimnames=list(c(),c("1-210PbCTS","2-214BiCTS","3-137CsCTS","4-210PbError","5-214BiError","6-137CsError","7-210PbHtAdjcpm","8-214BiHtAdjcpm","9-137CsHtAdjcpm","10-210PbErrcpm/g","11-214BiErrcpm/g","12-137CsErrcpm/g","13-Corr210Pbdpm/g","14-Corr214Bidpm/g","15-Corr137Csdpm/g","16-Corr210PbErr","17-Corr214BiErr","18-Corr137CsErr")))
  cor[,1]=(roi[,5]-roi[,4]-roi[,6])/(live/60)
  cor[,2]=(roi[,8]-roi[,7]-roi[,9])/(live/60)
  cor[,3]=(roi[,11]-roi[,10]-roi[,12])/(live/60)
  cor[,4]=sqrt(abs(roi[,5]-roi[,4]-roi[,6]))*cor[,1]/(roi[,5]-roi[,4]-roi[,6])
  cor[,5]=sqrt(abs(roi[,8]-roi[,7]-roi[,9]))*cor[,2]/(roi[,8]-roi[,7]-roi[,9])
  cor[,6]=sqrt(abs(roi[,11]-roi[,10]-roi[,12]))*cor[,3]/(roi[,11]-roi[,10]-roi[,12])
  cor[,7]=((cor[,1]-PbBk)/roi[,13])/((-0.017993*roi[,14])+1.5398)
  cor[,8]=((cor[,2]-BiBk)/roi[,13])/((-0.0116*roi[,14])+1.3477)
  cor[,9]=((cor[,3]-CsBk)/roi[,13])/((-0.012*roi[,14])+1.359)
  cor[,10]=abs(cor[,7]*cor[,4]/cor[,1])
  cor[,11]=abs(cor[,8]*cor[,5]/cor[,2])
  cor[,12]=abs(cor[,9]*cor[,6]/cor[,3])
  
  
  # This section performs the corrections for counting efficiency and the density of the dry sediments
  cor[,13]=(cor[,7]/(0.0285*(1+(0.033873*(1.9228-(2.904/((roi[,14]-2.4786)/9.4781)*roi[,13]))))))
  cor[,14]=cor[,8]/0.0419
  # Cs value (dpm/g) corrected for sampling date
  cor[,15]=(cor[,9]/0.09331)/exp(-0.03114*((roi[,3])/365))
  
    # Cumulative error x Correction Factors
  cor[,16]=cor[,13]*cor[,10]/cor[,7]
  cor[,17]=cor[,14]*cor[,11]/cor[,8]
  cor[,18]=cor[,15]*cor[,12]/cor[,9]
  
  # 'plotdpm' peels from 'cor'the corrected activities and errors (in dpm) to be used in the output plot
  plotdpm=(cor[,(13:18)])
  # 'bq' is the corrected values and errors converted into becquerels/g (i.e. plotdpm/60)
  plotbq=(plotdpm/60)
  colnames(plotbq)[(1:3)]=c('Corr210Pbbq/g', 'Corr214Bibq/g', 'Corr137Csbq/g')
  
  # Produce a graph of the Corrected Activities (with error bars) vs. Interval Depth
  # The graph is not saved by the script (if you wish to keep a copy you can do so using the RGUI
  midpt=(roi[,2]+roi[,1])/2
  plot_colours = c("blue","red","purple")
  
  # Determines whether the plot should be in units of dpm or bq
  graph.units <- match.arg(graph.units, c("dpm", "bq"), several.ok = FALSE)
  if (graph.units=="dpm") {plotdata=plotdpm}
  if (graph.units=="bq") {plotdata=plotbq}
  
  # Determines the ranges of the plot
  maxact = (max(plotdata[,1]))*1.05
  maxdepth = max(midpt)*1.05
  
  plot(midpt,plotdata[,1], type="o", pch=21, bg=plot_colours[1], col=plot_colours[1], ann=FALSE, ylim=c(0,maxact), xlim=c(0,maxdepth))
  lines(midpt, plotdata[,2], type="o", pch=22, bg=plot_colours[2], col=plot_colours[2])
  lines(midpt, plotdata[,3], type="o", pch=23, bg=plot_colours[3], col=plot_colours[3])

  
  
  # Creates the error bar function and adds them to each line of the graph
  # Potential for an error to display if an error bar goes off the bottom of the graph (i.e. below zero)
  # Such an error will appear as a prompt, but not stop the running of the script (safe to ignore)
  
  error.bar <- function(x, y, upper, lower=upper, length=0.06,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }
  
  error.bar(midpt, plotdata[,1], plotdata[,4], col=plot_colours[1])
  error.bar(midpt, plotdata[,2], plotdata[,5], col=plot_colours[2])
  error.bar(midpt, plotdata[,3], plotdata[,6], col=plot_colours[3])
  
  
  # Adds titles and a legend to the graph
  title(main=paste("Corrected Activities vs. Depth for", LakeName), col.lab=rgb(0,0.0,0), font.main=2)
  title(xlab="Interval Midpoint (cm)", col.lab=rgb(0,0.0,0))
  if (graph.units=="dpm"){title(ylab="Corrected Activity (dpm/g)", col.lab=rgb(0,0.0,0))}
  if (graph.units=="bq"){title(ylab="Corrected Activity (Bq/g)", col.lab=rgb(0,0.0,0))}
  
  #Formats the labels used in the legend
  Pbleglabel=expression(paste(""^210, " Pb"))
  Bileglabel=expression(paste(""^214, " Bi"))
  Csleglabel=expression(paste(""^137, " Cs"))
  
  legend("topright", inset= 0.075, legend=c(Pbleglabel,Bileglabel,Csleglabel), y.intersp=1.1, cex=1.1, 
         col=plot_colours, pch=21:23, pt.bg=plot_colours, bty="n");
  
  
  # Creation and population of the export matrix containing the Activity and Error Values
  activity.output=matrix(data=NA,nrow=nrow(roi),ncol=8,dimnames=list(c(),c("1-IntMidPt", "13-Corr210Pbdpm/g","14-Corr214Bidpm/g","15-Corr137Csdpm/g","16-Corr210PbErr","17-Corr214BiErr","18-Corr137CsErr", "8DateCorr")))
  activity.output[,1]=midpt
  activity.output[,(2:7)]=cor[,(13:18)]
  activity.output[,8]=roi[,3]
  
  # In order to match up with the density.output the number of empty rows cut from 'roi' are appended to 'activity.output'
  nogamma =matrix(data=NA, nrow=empty, ncol=8)
  nogamma[,1]=(d[(nrow(d)-(empty-1)):nrow(d),2]+d[(nrow(d)-(empty-1)):nrow(d),1])/2
  activity.output=rbind(activity.output, nogamma)
  binfordActivityOutput <- activity.output
  
# Returns the "binfordActivityOutput" dataframe for eventual use by the "binfordDates" function
  return(binfordActivityOutput)
}
