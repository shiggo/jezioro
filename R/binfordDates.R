#'binfordDates
#'
#'@description
#'\code{binfordDates} uses the output of \code{binfordRho} and \code{binfordActivity} to calculate 210Pb dates using the Constant Rate of Supply (CRS) model.
#'
#'@usage binfordDates(binfordRhoOutput, binfordActivityOutput, LakeName = "Lake Name", 
#'CoreName = "KB1", CoreDate = 2006, CIC = FALSE)
#'
#'@param binfordRhoOutput matrix or data frame containing the output data from \code{binfordRho}
#'@param binfordActivityOutput matrix or data frame containing the output data from \code{binfordActivity}
#'@param LakeName (defaults to "Lake Name"): name of the lake being analyzed
#'@param CoreName (defaults to "KB1"): name of the core being analyzed
#'@param CoreDate (defaults to 2006): year of core collection
#'@param CIC (defaults to FALSE): logical to control whether to also calculate and plot CIC dates (NOTE: this was still being worked on when the functions were rendered obsolete, and remains incomplete)
#'
#'@details \code{binfordDates} determines the unsupported fraction of 210Pb activity using the procedure outlined in Binford (1990) and then calculates dates based on the constant rate of supply model. For input, it requires the output data from \code{binfordRho} and \code{binfordActivity}. 
#'
#'@return
#'\code{binfordDatesOutput} returns a table with 15 columns:
#'\itemize{
#'  \item IntTop: depth (cm) to the top of the analyzed interval
#'  \item IntBot: depth (cm) to the bottom of the analyzed interval
#'  \item IntMid: depth (cm) to the midpoint of the analyzed interval
#'  \item TTop: age (years) at the top of the interval
#'  \item SDTTop: standard error of the age at the top of the interval
#'  \item TBot: age (years) at the bottom of the interval.
#'  \item SDTBot: standard error of the age at the bottom of the interval
#'  \item SedRate: sedimentation rate (g/cm^2/yr).
#'  \item SDSedRate: standard error of the sedimentation rate
#'  \item SumTop: cumulative dry mass from the bottom of the core to the top of the interval
#'  \item DTop: date at the top of the interval
#'  \item DBot: date at the bottom of the interval
#'  \item DMid: date at the midpoint of the interval
#'  \item CICDMid: CIC Date at the midpoint of the interval
#'  \item OrgSedRate: organic sedimentation rate
#'  }
#'
#'@examples 
#'#Import example input data
#'data(binfordRhoInput)
#'data(binfordActivityInput)
#'binfordDates(binfordRho(binfordRhoInput), binfordActivity(binfordActivityInput), 
#'LakeName="Test", CoreName="KB1", CIC=TRUE, CoreDate=2010)
#'
#'@author Adam Jeziorski, Joshua Thienpont
#'@export
#'
#'@references Binford MW (1990) Calculation and uncertainty analysis of 210PB dates for PIRLA project lake sediment cores. Journal of Paleolimnology 3: 253-267

binfordDates=function(binfordRhoOutput, binfordActivityOutput, LakeName="Lake Name", CoreName="KB1", CoreDate=2006, CIC=FALSE){
  density=binfordRhoOutput
  activity=binfordActivityOutput
  
  # A new data matrix "bin" is created that collates the relevant data from the two identified output files
  # It is important that your output files have the same intervals (i.e. all freeze-dried intervals need to be included in your activity sheet, just leave blank the intervals with no counts)
  levels=nrow(density)
  bin=matrix(data=NA,nrow=levels,ncol=9,dimnames=list(c(),c("1-IntTop", "2-IntBot","3-DateCorr","4-210PbCorr","5-210PbCorrDate", "6-CumMassTop", "7-CumMassBot","8-Rho", "9-210PbErr1Stdev")))
  
  bin[,1]=density[,1]
  bin[,2]=density[,2]
  bin[,3]=activity[,8]
  bin[,4]=activity[,2]
  bin[,5]=exp(0.00008532*bin[,3])*bin[,4]
  bin[,6]=density[,3]
  bin[,7]=density[,4]
  bin[,8]=density[,5]
  bin[,9]=activity[,5]
  
  
  # As not every sample from your density matrix was necessarily counted (i.e. some samples that were prepped, weren't counted due to hitting background)
  # the matrix "bin2" reduces the size of "bin" by cutting out any intervals for which there is no value for the corrected 210Pb activity
  
  # First there is a count to determine how many intervals need to be removed "t"
  t=0
  for (i in 1:levels){
    if (is.na(bin[i,4])){t=t+1}
  }
  
  # The matrix bin2 has "t" less rows than bin, and all rows that do not have a NA value for corrected 210Pb activity are copied from bin
  bin2=matrix(data=NA,nrow=(levels-t),ncol=9,dimnames=list(c(),c("1-IntTop", "2-IntBot","3-DateCorr","4-210PbCorr","5-210PbCorrDate","6-CumMassTop", "7-CumMassBot","8-Rho", "9-210PbErr1Stdev")))
  j=1
  for (i in 1:levels) {
    if (is.na(bin[i,4])){next}
    else bin2[j,]=bin[i,]
    j=j+1
  }
  
  # Application of the "Binford Rule"
  # In the first iteration of the loop, the mean and standard deviation of the bottom three corrected 210Pb activities are calculated.
  # If the corrected activity of the next interval up is greater than the calculated mean + 1stdev then loop exits as it is above background
  # However if this condition is not met the loop continues, this time including the bottom four activities
  
  # "k" is a counter that will be added to if the requirements of the binford procedure are not met (it begins at 2 so that the loop begins with 3 intervals)
  k=2
  for (i in (levels-t):1) {
    binmean = mean(bin2[((levels-t-k):(levels-t)),5])
    binsd = sd(bin2[((levels-t-k):(levels-t)),5])
    if (bin2[(levels-t-(k+1)),5]>=(binmean+binsd)){next}
    else k=k+1 
  }
  
  
  # Subtracts the supported 210Pbvalue (binmean for the # of intervals required to satisfy the "Binford Rule") from the corrected activities 
  # to create the unsupported 210Pb activity for each interval "UnsuppPb"
  UnsuppPb=bin2[,5]-binmean
  
  
  # The first step in creating an export matrix for the Binford correction program
  # It is initally populated with zeroes so that "8-PercOrg" has something in it when exported
  binexp=matrix(data=0,nrow=(levels-t),ncol=9, dimnames=list(c(),c("1-IntTop", "2-IntBot","3-210PbCorr","4-UnsuppPb","5-Rho", "6-CumMassTop", "7-CumMassBot","8-PercOrg", "9-SDErrAct")))
  binexp[,(1:2)]=bin2[,(1:2)]
  binexp[,3]=bin2[,5]
  binexp[,4]=UnsuppPb
  binexp[,5]=bin2[,8]
  binexp[,(6:7)]=bin2[,(6:7)]
  binexp[,9]=bin2[,9]
  
  
  # Checks downward through the unsupported 210Pb activity "bin[,4]" once it hits a negative value
  # that value and every value below it are set to zero
  for (i in 1:(levels-t)) {
    if (binexp[i,4]<=0){binexp[(i:(levels-t)),4]=0}
  }
  

  # Establishes the # of intervals with an activity value and assigns that value to "valid"
  for (i in (levels-t):1) {
    if (binexp[i,4]==0.0) {valid=i-1}
  }
  
  
  # Creation of a data matrix called "cal" with the 29 columns listed in the Binford paper, this matrix will be filled in by the various calculations
  # Nine columns are copies from "binexp", note that the matrix has (levels-t) rows
  cal=matrix(data=NA,nrow=(levels-t),ncol=29,dimnames=list(c(),c("01inttop","02intbot","03totact","04unsupact","05rho","06percorg","07rhoxcee","08integseg","09integint","10sumbot","11sumtop","12ttop","13tbot","14sedrate","15orgsedrt","16cummasst","17cummassb","18sdact","19sdrho","20sdunsup","21sdrhocee","22sdintseg","23sdintint","24sdsumbot","25sdsumtop","26sdttop","27sdtbot","28sdsedrt","29sdorgrt")))
  cal[,1]=binexp[,1]
  cal[,2]=binexp[,2]
  cal[,3]=binexp[,3]
  cal[,4]=binexp[,4]
  cal[,5]=binexp[,5]
  cal[,6]=binexp[,8]
  cal[,16]=binexp[,6]
  cal[,17]=binexp[,7]
  cal[,18]=binexp[,9]
  
  
  # Calculation of the standard deviation of unsupported Pb-210
  cal[,20]=sqrt(cal[,18]^2+binsd^2)
  
  
  # Fill the Rho X C(X) column
  cal[,7]=cal[,4]*cal[,5]
  
  
  # The standard deviation of RHO is filled in (This variable is assumed to be zero in the PIRLA study - see Binford, 1990)
  cal[,19]=0
  
  
  # Fill standard deviation of RHOCEE column in matrix
  cal[,21]=sqrt(cal[,19]^2*cal[,4]^2+cal[,20]^2*cal[,5]^2)
  
  
  # Fill integrated segment Pb210 Activity column in matrix
  cal[,8]=cal[,7]*(cal[,2]-cal[,1])
  
  
  # Fill standard deviation of integrated segment column in matrix
  cal[,22]=sqrt(cal[,21]^2*(cal[,2]-cal[,1])^2)
  
  
  # Fill integrated interval (between segments) of the Pb210 Activity column in matrix (entered into column 9 of "cal")
  # Fill standard deviation of integrated Pb210 of interval (entered into column 23 of "cal")
  
  # NOTE: The first equation is derived from the exponential rule of integration
  # cal[i,23]=cal[i,9]*sqrt(((sqrt(cal[i,22]^2+cal[i+1,22]^2)/(cal[i,7]-cal[i+1,7]))^2)+((sqrt((cal[i,22]/cal[i,7])^2+(cal[i+1,22]/cal[i+1,7])^2))/log(cal[i,7]/cal[i+1,7]))^2)
  # The second equation, which is used, is derived from the trapezoidal rule of integration
  
  # Adds zeros in the bottom row for columns 9 and 23 to create a floor for the integration
  cal[(levels-t),9]=0.00
  cal[(levels-t),23]=0.00
  
  for (i in 1:(levels-t-1)) {
    cal[i,9]=((cal[i+1,7]-cal[i,7])/log((cal[i+1,7]+0.0000001)/cal[i,7]))*(cal[i+1,1]-cal[i,2])
    cal[i,23]=sqrt(.25*(cal[i,22]^2+cal[i+1,22]^2)*(cal[i+1,1]-cal[i,2])^2)
  }
  
  
  # Fill cumulative Pb210 columns in matrix ([,10] is sum to bottom of interval, [,11] is sum to top of interval)
  # Assigns a zero to column 10 and 11 to the bottom row to create a starting point for the calculation
  cal[(levels-t),(10:11)]=0.00
  
  for (i in seq((levels-t-1),1,-1)) {
    cal[i,10]=cal[i,9]+cal[i+1,11]
    cal [i,11] =cal[i,10]+cal[i,8]
  }
  
  
  # Fill standard deviation columns in matrix for cumulative Pb210 (24 is sum to the bottom of the interval, 25 to top)
  # Assigns a zero to column 25 for the bottom interval to create a starting point for the calculation
  # Sneaking suspicion that this is where the slight difference in the errors from binford.exe are creeping in
  cal[(levels-t),25]=0.00
  
  for (i in seq((levels-t-1),1,-1)) {
    cal[i,24]=sqrt(cal[i,23]^2+cal[i+1,25]^2)
    cal[i,25]=sqrt(cal[i,22]^2+cal[i,24]^2)
  }
  
  
  # ***********************************************************************************
  # *********************************Correction Factor*********************************
  # ***********************************************************************************
  
  # Add in the small tail of cumulative Pb210 not accounted for in the numerical integration
  # A new matrix "cor" is defined, this is where the geometric calculations necessary to determine the correction factor take place
  E=(levels-t)*2
  cor=matrix(data=0,nrow=E,ncol=2)
  
  for (i in seq(1,E,2)) {
    f=(i/2+0.5)
    cor[i,1]=cal[f,11]
    cor[i,2]=cal[f,16]
    cor[i+1,1]=cal[f,10]
    cor[i+1,2]=cal[f,17]
  }
  
  for (i in 1:3) {
    A80=0.2*cor[1,1]
    A90=0.1*cor[1,1]
    A95=0.05*cor[1,1]
    for (m in 1:E) {
      if ((m+1)>=E){next}
      if (cor[m,1] >=0.000001 && cor[m+1,1] <= 0.000001){T=m}
      if (cor[m,1] ==0) {next}
      if (A80>=cor[m+1,1] && A80<=cor[m,1]){X=m}
      if (A90>=cor[m+1,1] && A90<=cor[m,1]){Y=m}
      if (A95>=cor[m+1,1] && A95<=cor[m,1]){Z=m}
    }
    
    L8=(log(A80/cor[X,1])/log(cor[X+1,1]/cor[X,1]))*(cor[X+1,2]-cor[X,2])+cor[X,2]
    L9=(log(A90/cor[Y,1])/log(cor[Y+1,1]/cor[Y,1]))*(cor[Y+1,2]-cor[Y,2])+cor[Y,2]
    L5=(log(A95/cor[Z,1])/log(cor[Z+1,1]/cor[Z,1]))*(cor[Z+1,2]-cor[Z,2])+cor[Z,2]
    LB=(2.*L5+2.822*(3.*L5-2.*L9-L8))/2.
    ZE=(cor[T,2]-L5)/(LB-L5)
    AN=(cor[1,1])*(sqrt(.001^(1.+ZE)*.05^(1.-ZE)))
    
    for (g in 1:(T+1)) {
      cor[g,1]=(cor[g,1]+AN)
    }
  }
  
  for(p in seq(1,T,2)) { 
    s=p/2.+.5
    cal[s,11]=cor[p,1]+AN
    cal[s,10]=cor[p+1,1]+AN
  }
  
  # ***********************************************************************************
  # *********************************Correction Factor*********************************
  # ***********************************************************************************
  
  
  # Calculate times for each of the tops and bottoms of the intervals
  for (i in 1:(levels-t)) {
    if (cal[i,10]<0.001){next}
    else cal[i,12]=(log(cal[1,11]/cal[i,11]))/.03114
    cal[i,13]=(log(cal[1,11]/cal[i,10]))/.03114       
    
    
    # Calculate standard deviations for the time at both the top and bottom of the intervals
    cal[i,26]=sqrt(((1/(.03114^2*cal[1,11]^2))*cal[1,25]^2)+((1/(.03114^2*cal[i,11]^2))*cal[i,25]^2))
    cal[i,27]=sqrt(((1/(.03114^2*cal[1,11]^2))*cal[1,25]^2)+((1/(.03114^2*cal[i,10]^2))*cal[i,24]^2))
  }
  
  
  # Calculate bulk sedimentation rates. Note that the A(X) value is calculated at the mid-point of the interval,
  # and not top of the interval: ((cal(N,11)+cal(N,10))/2
  for (i in 1:(levels-t)) {
    if (cal[i,4]<0.000001) {next}
    else cal[i,14]=.03114*((cal[i,11]+cal[i,10])/2)/cal[i,4]
    # Calculate standard deviation for bulk sedimentation rates
    cal[i,28]=sqrt(.03114^2*cal[i,11]/cal[i,4]*sqrt((cal[i,25]/cal[i,11])^2+(cal[i,20]/cal[i,4])^2))
  }
  
  
  # Calculate Organic sedimentation rate
  for (i in 1:(levels-t)) {
    cal[i,15]=cal[i,14]*cal[i,6]
  }
  
  
  # Calculate 210Pb dates using the Constant Initial Concentration (CIC) model
  # 'CICmatrix' is kept seperate as it is not part of the original Binford program (plus no errors calculated)
  # These values are only used if the argument "CIC" is set to TRUE
  CICmatrix=matrix(data=NA,nrow=(levels-t),ncol=3,dimnames=list(c(),c("UnsuppPb","CICbp","CICcal")))
  decay=0.03114
  CICmatrix[,1]=binexp[,4]
  CICmatrix[,2]=(decay^-1)*log(binexp[1,4]*(binexp[,4]^-1))
  CICmatrix[,3]=CoreDate-(CICmatrix[,2])
  
  
  # Creation of the Export Matrix (NOTE: This is limited to the number of "valid" intervals (i.e. had a UnSupp Pb210 activity value)
  export=matrix(data=NA,nrow=valid,ncol=15,dimnames=list(c(),c("IntTop","IntBot","IntMid","TTop","SDTTop","TBot","SDTBot","SedRate","SDSedRate","SumTop","DTop","DBot","DMid","CICDMid","OrgSedRate")))
  export[(1:valid),1]=cal[(1:valid),1]
  export[(1:valid),2]=cal[(1:valid),2]
  export[,3]=(export[,2]+export[,1])/2
  export[(1:valid),4]=cal[(1:valid),12]
  export[(1:valid),5]=cal[(1:valid),26]
  export[(1:valid),6]=cal[(1:valid),13]
  export[(1:valid),7]=cal[(1:valid),27]
  export[(1:valid),8]=cal[(1:valid),14]
  export[(1:valid),9]=cal[(1:valid),28]
  export[(1:valid),10]=cal[(1:valid),11]
  export[,11]=CoreDate-export[,4]
  export[,12]=CoreDate-export[,6]
  export[,13]=CoreDate-((export[,4]+export[,6])/2)
  if (CIC == TRUE){
    export[(1:valid),14]=CICmatrix[(1:valid),3]
  } else {
    export[(1:valid),15]=cal[(1:valid),15]
  }

  # Creates the error bar function and adds them to each line of the graph
  # Potential for an error to display if error bar smaller than one pixel, such an error will not stop the the script (safe to ignore)
  error.bar <- function(x, y, upper, lower=upper, length=0.06,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }
  
  
  # Allows the graphs to be displayed in a single window in a 2x2 matrix
  par(mfcol=c(2,2) ) 
  
  # logPb210 Activity vs. CumMass graph
  logunsupfig=matrix(data=NA,nrow=(valid),ncol=3,dimnames=list(c(),c("unsuppPb","logunsuppPb","cumdrymass")))
  logunsupfig[,1]=binexp[(1:valid),4]
  logunsupfig[,2]=log10(logunsupfig[,1])
  logunsupfig[,3]=binexp[(1:valid),7]
  reg1=lm(logunsupfig[,2]~logunsupfig[,3])
  plot(logunsupfig[,3],logunsupfig[,2],type="p", pch=21, las = 1, ann=FALSE)
  title(main=paste("Cumulative Dry Mass vs. \n log (Unsupp 210Pb) for ",LakeName),col.lab=rgb(0,0.0,0), font.main=2)
  title(xlab=expression(paste("log(Unsupporte",d ^210*Pb ," Activity)")), col.lab=rgb(0,0.0,0), ylab=expression(paste("Cumulative Dry Mass (g/",cm^2, ")")), col.lab=rgb(0,0.0,0))
  abline(reg1)
  
  
  # Date vs. Depth graph (CRS)
  plot(export[,3], export[,13],type="p", pch=21, las = 1, ann=FALSE, ylim=c((min(export[,13])-30), (CoreDate+5)))
  graph3title= paste("Date vs. Depth (CRS Model) \n for", LakeName)
  title(main=graph3title, col.lab=rgb(0,0.0,0), font.main=2)
  title(xlab="Interval Midpoint Depth (cm)", col.lab=rgb(0,0.0,0), ylab="Date (yr)", col.lab=rgb(0,0.0,0))
  error.bar(export[,3],export[,13], ((export[,5]+export[,7])/2))
  #trendcrs=lm(export[,13]~poly(export[,3],4, raw=TRUE)
  
  
  # Depth vs. SedRate graph (CRS)
  plot(export[,8],export[,3], type="o", pch=21, ylim = c((max(export[,3])*1.04), 0), las = 1, ann=FALSE)
  graph1title= paste("Depth vs. Sedimentation Rate \n for", LakeName)
  title(main=graph1title, col.lab=rgb(0,0.0,0), font.main=2)
  title(xlab=expression(paste("Sedimentation Rate (g/", cm^2, "/yr)")), col.lab=rgb(0,0.0,0), ylab="Interval Midpoint Depth (cm)", col.lab=rgb(0,0.0,0))
  error.bar(export[,8],export[,3], export[,9])
  
  
  
  # Depth vs. Dates graph (CIC) (only drawn if CIC=TRUE)
  if (CIC == TRUE){
    plot(export[,3], export[,14],type="p", pch=21, las = 1, ann=FALSE)
    graph4title= paste("Date vs. Depth (CIC Model) \n for", LakeName)
    title(main=graph4title, col.lab=rgb(0,0.0,0), font.main=2)
    title(xlab="Interval Midpoint Depth (cm)", col.lab=rgb(0,0.0,0), ylab="Date (yr)", col.lab=rgb(0,0.0,0))
  }
  
  # Returns the 'export' table as  "binfordDatesOutput" 
  binfordDatesOutput <- export
  return(binfordDatesOutput)
}
