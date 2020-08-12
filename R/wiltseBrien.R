#'wiltseBrien
#'
#'@description
#'\code{wiltseBrien} performs a Brien's test on a correlation matrix to test whether changes in the underlying variables are coherent. 
#'
#'This function is detailed in Brendan Wiltse's PhD thesis and calculates the grand mean, main effects interactions and equal correlation for the correlation matrix and returns the degrees of freedom, chi-squared value and p-value for each test.
#'
#'@usage wiltseBrien(my.data, print=T, print.detail=T)
#'
#'@param my.data A data frame containing the variables to be tested. The variables should be arranged in columns with appropriate column names
#'@param print (Defaults to TRUE): logical to indicate whether to display the degrees of freedom, chi-squared, and p-values for the grand mean, main effects, interactions, and equal correlations
#'@param print.detail (Defaults to TRUE): logical to indicate whether to display the variables included in the test and the correlation coefficient for the grand mean
#'
#'@details \code{wiltseBrien} Full details on the Brien's test and the underlying mathematics can be found in Brien et al. 1984
#'
#'@return \code{output} table of results for use with other functions
#'
#'@examples 
#'#Import the example input data for use with Brendan Wiltse's functions
#'data(wiltseInput)
#'wiltseBrien(wiltseInput[,2:9])
#'@author Brendan Wiltse
#'@references Brien CJ, Venables WN, James AT, Mayo O (1984) An analysis of correlation matrices: equal correlations. Biometrika 71: 545-54
#'
#'Wiltse B (2014) The response of Discostella species to climate change at the Experimental Lakes Area, Canada. PhD Thesis Queen's University
#'@export

wiltseBrien=function(my.data, print=T, print.detail=T){
  #Calculate correlation matrix from data set
  my.cor<-cor(my.data)
  #Determine number of variables
  n=length(my.data[1,])
  #Determine number of samples
  s=length(my.data[,1])
  #Calculate the mean correlation for each variable
  meanr=0
    for(i in 1:n){
    meanr[i]=((sum(my.cor[,i])-1)/(n-1))
    }
  
  #Calculate Grand Mean of matrix
  grandmean<-(sum(my.cor)-n)/(n*(n-1))
  
  #Calculate z-transformed correlation matrix with proper formatting
  my.zcor=matrix((0.5*(log(1+my.cor)-log(1-my.cor))),n,n)
  Changed <- as.character(my.zcor)
  for (i in 1:length(Inf)){
    Changed <- replace(Changed, Changed == Inf[i], 0[i])
  }
  my.zcor=as.numeric(Changed)
  my.zcor=data.frame(matrix(my.zcor,n,n))
  
  #Calculate z++
  z1=sum(my.zcor)/2
  
  #Calculate zi+
  z2=0
  for(i in 1:n){
    z2[i]=sum(my.zcor[ ,i])
  }
  
  #Calculate z~i+
  z3=0
  for(i in 1:n){
    z3[i]=z2[i]/(n-2)
  }
  
  #Calculate zi-blah
  z4=0
  for(i in 1:n){
    z4[i]=((z2[i]-(2*z1/n))^2/(n-2))
  }
  
  #Calculate var of z-transform
  z5=0
  for(i in 1:n){
    z5[i]=var(my.zcor[,i])
  }
  
  #Calculate z~++
  z6=(2*z1)/((n-1)*(n-2))
  
  #Calculate z~j+
  z7=0
  for(i in 1:n){
    z7[i]=sum(my.zcor[i, ])/(n-2)
  }
  
  #Calculate Q2 intermediate matrix
  my.zcor1=as.matrix(my.zcor)
  my.zcor1=c(my.zcor1)
  z3m=matrix(z3,n,n)
  z7m=matrix(z7,n,n,byrow=T)
  z8=0
  for(i in 1:n^2){
    z8[i]=(my.zcor1[i]-z3m[i]-z7m[i]+z6)^2
  }
  z8=matrix(z8,n,n)
  for(i in 1:n){
    z8[i,i]=0
  }
  
  #Calculate c^
  c1=(2*z1)/(n*(n-1))
  
  #Calculate z^
  z9=(exp(c1)-exp(-1*c1))/(exp(c1)+exp(-1*c1))
  
  #Calculate n(samples-3)
  s1=s-3
  
  #Calculate Q0, Q1, Q2
  Q0=((2*(z1^2))/(n*(n-1)))
  Q1=sum(z4)
  Q2=sum(z8)/2
  
  #Calculate L0,L1,L2
  L0=((1+(n-1)*0)^2)/(s1*((1+0)^2))
  L1=(n-((n-2)*(1-z9)^2))/((2*s1)*(1+z9)^2)
  L2=1/(s1*(1+z9)^2)
  
  #Calculate d.f. for each test
  grandmeandf=1
  maineffectsdf=n-1
  interactionsdf=0.5*n*(n-3)
  equalcordf=0.5*n*(n-1)-1
  totaldf=0.5*n*(n-1)
  d.f.=c(grandmeandf,maineffectsdf,interactionsdf,equalcordf,
         totaldf)
  
  #Calculate chi-squre for each test
  grandmeanx2=Q0/L0
  maineffectsx2=Q1/L1
  interactionsx2=Q2/L2
  equalcorx2=maineffectsx2+interactionsx2
  chi2=c(grandmeanx2,maineffectsx2,interactionsx2,equalcorx2,
         NA)
  
  #Calculate p-value (one-tail) for each test
  grandmeanp=dchisq(grandmeanx2, grandmeandf)*2
  maineffectsp=dchisq(maineffectsx2, maineffectsdf)*2
  interactionsp=dchisq(interactionsx2, interactionsdf)*2
  equalcorp=dchisq(equalcorx2, equalcordf)*2
  p.value=c(grandmeanp,maineffectsp,interactionsp,equalcorp,
            NA)
  
  #Set up variables for printing
  title=c("Test", "p-value")
  grand=c("Grand Mean:", grandmeanp)
  main=c("Main Effects:", maineffectsp)
  inter=c("Interactions:",interactionsp)
  equal=c("Equal Correlation:",equalcorp)
  titles=c("Grand Mean:","Main Effects:","Interactions:","Equal Correlation:", "Total:")
  
  #Set up data frame for main summary
  output=data.frame(d.f., chi2, p.value, row.names=titles)
  
  #Printing detailed output
  if(isTRUE(print.detail)){
    var=c("Variables:", names(my.data))
    grand=c("Grand Mean:", grandmean)
    cat(var, '\n', grand, '\n', '\n')
  }
  
  #Printing output
  if(isTRUE(print)){
    print(output, justify="left")
  }
  
  #Assigning output to globalenv for use in other functions
  #assign("Brien.output", output, envir=globalenv())
  #AJ: Instead of assigning output, my.cor, meanr, and grandmean to the global environment, everything is merged as a list that is returned invisibly for assignment if necessary
  output <- list("wiltseBrien Output"=output, "variables"= names(my.data), "grandmean"=grandmean, "meanr"=meanr, "Correlation Matrix"=my.cor)
  invisible(output)
}