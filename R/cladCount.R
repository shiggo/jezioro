#'cladCount
#'
#'@description
#'\code{cladCount} determines the maximum number of individuals represented in raw counts of cladoceran subfossils
#'
#'@usage cladCount(d, percCutoff=2, sampleCutoff=2, outputType="indiv")
#'
#'@param d matrix or data frame containing the input data
#'@param percCutoff (defaults to 2): minimum relative abundance (i.e. \%) required for a taxon to be included in the the reduced subset ('gt' - abbreviation of 'greater than')
#'@param sampleCutoff (defaults to 2): minimum number of samples a taxon must be present in with at least a relative abundance of \code{percCutoff} for inclusion in the reduced subset ('gt')
#'@param outputType (defaults to 'indiv'): the format of the output, either individuals ('indiv'), relative abundance ('perc'), or the reduced subset ('gt')
#'
#'@details The input data should be in the form of a matrix or data frame with at least 3 columns:
#'\itemize{
#'   \item Column 1: taxon name
#'   \item Column 2: subfossil name
#'   \item Column 3+: count of each taxon/subfossil remain type
#'   }
#'NOTE: Special characters in the column headings may interfere with import of the input file. The script uses string matching to identify body parts, anything other than the following list in the body part column (i.e. spelling errors) will be ignored.
#'\itemize{
#'   \item Body part (number of subfossils per individual)
#'   \item Headshield (1)
#'   \item Carapace (2)
#'   \item PA Claw (2)
#'   \item Postabdomen (1)
#'   \item Mandible (2)
#'   \item Caudal Furca (2)
#'   \item Exopodite Segment (2)
#'   \item Basal Exp Segment (2)
#'   \item Exp Segment 2 (2)
#'   \item Exp Segment 3 (2)
#'   \item Antennule (2)
#'   \item Tail Stem (1)     
#'   }
#'
#'@return
#'\code{cladCount} will return one of: 
#'\itemize{
#'  \item \code{indiv}: table listing the minimum number of individuals of each taxa represented in each interval
#'  \item \code{perc}: table listing the relative abundance of each taxa within each interval
#'  \item \code{gt}: relative abundances of those taxa meeting the \code{sampleCutoff} and \code{percCutoff} criteria
#'  }
#'
#'@author Adam Jeziorski, Andrew Labaj
#'@export
#'
#'@examples 
#'#load example cladoceran count data
#'data(cladCountInput)
#'cladCount(cladCountInput)
#'
#'#Return the values for only taxa with greater than 4 percent abundance in at least 2 samples
#'cladCount(cladCountInput, percCutoff=4, sampleCutoff=2)
#'
#'@references Korhola A, Rautio M (2001) 2. Cladocera and other branchiopod crustaceans. In: Smol JP, Birks HJB, Last WM (eds.) Tracking Environmental Change Using Lake Sediments. Volume 4: Zoological Indicators. Kluwer Academic Publishers, Dordrecht, The Netherlands, pp 4-41

cladCount=function(d, percCutoff=2, sampleCutoff=2, outputType="indiv"){
 
data=d
#Ensure the species names were imported as a vector
data[,1]=as.vector(data[,1])
#Identify rows without a species name in the first column and replace empty value with NA
taxa.check=(grepl(" ", data[,1])==FALSE)
data[taxa.check,1]=NA

#Replace all missing values in "data" with a zero
data[is.na(data)]=0

#Creates a vector 'names' that is the subset of the individual species names (i.e. removes the zeros)
#The value 'taxa' is the number of individual species
names=data[,1]
names=subset(names, names!=0)
taxa=length(names)

#Creates a vector "part" of the number of body part types per taxa
#'t' increases by 1 for each new taxa
part=rep(1, taxa)
t=0
for (i in 1:(nrow(data))){
  if (data[i,1]!=0){t=t+1}
  else (part[t]=part[t]+1)
}

#Creates the matrix "indiv", which will be populated by the # of individuals using the raw counts
indiv=matrix(data=NA,nrow=taxa,ncol=(ncol(data)-2))
rownames(indiv)=names
colnames(indiv)=dimnames(data)[[2]][3:ncol(data)]

#The following for loop iterates once for each species in the raw data table
#"start" is used to determine the starting point for each iteration of the loop, it is increased by the appropriate value of the "part" vector each time
start=1

#The "part" vector is used to define the size of a new matrix "sp"" and is populated with the raw counts for each of the species body part types
for(i in 1:taxa){
  sp=matrix(data=0,nrow=(part[i]),ncol=(ncol(data)-2),dimnames=NULL)
  sp=data[(start:(start+(part[i])-1)),(3:ncol(data))]

#A series of vectors is created, one for each potential body part
  HS=vector(mode="integer", length=(ncol(data)-2))
  C=vector(mode="integer", length=(ncol(data)-2))
  PAC=vector(mode="integer", length=(ncol(data)-2))
  PA=vector(mode="integer", length=(ncol(data)-2))
  M=vector(mode="integer", length=(ncol(data)-2))
  CF=vector(mode="integer", length=(ncol(data)-2))
  ES=vector(mode="integer", length=(ncol(data)-2))
  TS=vector(mode="integer", length=(ncol(data)-2))
  BES=vector(mode="integer", length=(ncol(data)-2))
  ES2=vector(mode="integer", length=(ncol(data)-2))
  ES3=vector(mode="integer", length=(ncol(data)-2))
  AN=vector(mode="integer", length=(ncol(data)-2))

#A loop to check the part-type and assign the count values to the appropriate vector
  for (j in 1:(part[i])){
    if(grepl("Headshield", data[(start+j-1),2])==TRUE){HS=sp[j,]}
    if(grepl("Carapace", data[(start+j-1),2])==TRUE){C=sp[j,]}
    if(grepl("PA Claw", data[(start+j-1),2])==TRUE){PAC=sp[j,]}
    if(grepl("Postabdomen", data[(start+j-1),2])==TRUE){PA=sp[j,]}
    if(grepl("Mandible", data[(start+j-1),2])==TRUE){M=sp[j,]}
    if(grepl("Caudal Furca", data[(start+j-1),2])==TRUE){CF=sp[j,]}
    if(grepl("Exopodite Segment", data[(start+j-1),2])==TRUE){ES=sp[j,]}
    if(grepl("Tail Stem", data[(start+j-1),2])==TRUE){TS=sp[j,]}
    if(grepl("Basal Exp Segment", data[(start+j-1),2])==TRUE){BES=sp[j,]}
    if(grepl("Exp Segment 2", data[(start+j-1),2])==TRUE){ES2=sp[j,]}
    if(grepl("Exp Segment 3", data[(start+j-1),2])==TRUE){ES3=sp[j,]}
    if(grepl("Antennule", data[(start+j-1),2])==TRUE){AN=sp[j,]}
  }

#A new data frame "sp2" is made from combining each of the part vectors (after modification to calculate the maximum # of indiviuals represented by each body part)
#Note that it rounds UP
  sp2=rbind(HS, ceiling(C/2), ceiling(PAC/2), PA, ceiling(M/2), ceiling(CF/2), ceiling(ES/2), TS, ceiling(BES/2), ceiling(ES2/2), ceiling(ES3/2), ceiling(AN/2))

#The vector count is used to determine the maximum value of individuals represented by any single body part
  count=vector()
  for (k in 1:(ncol(data)-2)){
    count=c(count, max(sp2[,k]))
  }

#The matrix "indiv" is populated with the count data
  indiv[i,(1:(ncol(data)-2))]=count
  start=start+part[i]
}

#Determines species with no presence in the dataset and removes them from indiv
sumindiv=apply(X=indiv, MARGIN=1, FUN=sum)
indiv=subset(indiv, sumindiv>0)

#A copy of the indiv table 'perc', is filled with % Data
percentage=function(d){d/sum(d)*100}
perc=apply(X=indiv, MARGIN=2, FUN=percentage)

#Removes rare species (creates a subset of 'perc' called 'gt')
#Includes any species where the 2nd highest abundance is >2% (i.e. >2% in 2 intervals)
oneless=function(d){sort(d,partial=(length(d)-(sampleCutoff-1)))[length(d)-(sampleCutoff-1)]}
secondhighest=apply(X=perc, MARGIN=1, FUN=oneless)
gt=perc[-which(secondhighest<percCutoff),]

indiv=t(indiv)
perc=t(perc)
gt=t(gt)

#Determines which data frame to return based upon the value of 'outputType'
if(outputType=="indiv"){
  return(indiv)
} else if (outputType=="perc"){
    return(perc)
} else if (outputType=="gt"){
    return(gt)
} else print("outputType must be one of 'indiv', 'perc', or 'gt'")
  
}