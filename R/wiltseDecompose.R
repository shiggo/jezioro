#'wiltseDecompose
#'
#'@description
#'\code{wiltseDecompose} builds upon the Brien's test performed by \code{wiltseBrien} to identify homogenous and coherent subsets within a correlation matrix.
#'
#'First, a test for homogeneity and coherence within a correlation matrix is performed. If heterogeneity is detected (a homogenous matrix fails to reject the null-hypothesis for the tests of main effects, interactions and equal correlations), then the lake with lowest mean correlation within the matrix is deleted and the test is rerun on the reduced matrix. This process is repeated until the remaining variables are both homogenous and synchronous (i.e. fail to reject the null-hypothesis for the grand mean).
#'
#'To determine multiple sets of homogenous and synchronous variables within the data, the test must be run multiple times. Include all variables in the first run, and in subsequent reuns include only those variables previously rejected.
#'
#'This function is detailed in Brendan Wiltse's PhD thesis and calculates the grand mean, main effects interactions and equal correlation for the correlation matrix and returns the degrees of freedom, chi-squared value and p-value for each test. 
#'
#'@usage wiltseDecompose(my.data, print.detail=F, return.detail=F)
#'
#'@param my.data A data frame containing the variables to be tested. The variables should be arranged in columns with appropriate column names
#'@param print.detail (Defaults to FALSE): logical to indicate whether to display output from \code{wiltseBrien} for each step of the process. It is recommended that this is set to TRUE when first using the function to better understand the process of matrix decomposition.
#'@param return.detail (Defaults to FALSE): logical to indicate whether to return the full list output for further manipulation.
#'
#'@details \code{wiltseDecompose} A complete description of the Brien's test can be found in Brien et al. 1984, and the procedure for matrix decomposition follows the methodology described in Rusak et al. 1999
#'
#'@return
#'\code{output} table of results for use with other functions.
#'
#'@examples 
#'#Imports the example input data for use with Brendan Wiltse's functions
#'data(wiltseInput)
#'wiltseDecompose(wiltseInput[,2:9])
#'@author Brendan Wiltse
#'@references Brien CJ, Venables WN, James AT, Mayo O (1984) An analysis of correlation matrices: equal correlations. Biometrika 71: 545-54
#'
#'Rusak JA, Yan ND, Somers KM, McQueen DJ (1999) The temporal coherence of zooplankton population abundances in neighboring north-temperate lakes. The American Naturalist 153: 46–58
#'
#'Wiltse B (2014) The response of Discostella species to climate change at the Experimental Lakes Area, Canada. PhD Thesis Queen's University
#'@export

wiltseDecompose=function(my.data, print.detail=F, return.detail=F){
  #Run Brien on original data set
  BrienOut=wiltseBrien(my.data, print=print.detail, print.detail=print.detail)
  #Assign data set to new variable
  my.data1=my.data  
  
  #Decompose matrix based on variable with the lowest correlation
  while(BrienOut[[1]]$p.value[2]<=0.05 ||
        BrienOut[[1]]$p.value[3]<=0.05 ||
        BrienOut[[1]]$p.value[4]<=0.05){
    b=which.min(BrienOut$meanr)
    my.data1[,b]=NULL
    BrienOut=wiltseBrien(my.data1, print=print.detail, print.detail=print.detail)
    if(length(my.data1[1,])==3){grandmean=c("Grand Mean:",
                                            BrienOut$grandmean)
    coherentsubset=c("Coherent
Subset:",names(my.data1))
    cat(coherentsubset, '\n', BrienOut$grandmean,
        '\n', '\n', "Significance Tests", '\n')
    print(BrienOut[[1]])}
    if(length(my.data1[1,])==3)stop("Only Three Variables
Remain")
  }
  #Prepare output for printing
  grandmean=c("Grand Mean:", BrienOut$grandmean)
  coherentsubset=c("Coherent Subset:",names(my.data1))
  cat(coherentsubset, '\n', BrienOut$grandmean, '\n', '\n',
      "Significance Tests", '\n')
  print(BrienOut[[1]])
  if(return.detail != F){return(BrienOut)}
}

