#### EXACT MLE###########################################
###LOG L FUNCTION###
##First the log-likelihood function should be created based on the equion (4) on the pdf file
### The name of the function: Logll 
### Its arguments: TL: lower end of the intervals/ freq: Frequency on each interval/ Theta: the arguments including mu and sigma
Logll <- function(TL,freq,theta){ #Function Arguments
  m<- length(TL) # the total number of intervals
  
  if( (pnorm(theta[2]*TL[2]-theta[1])) < 1e-16  ){ # as the cdf of normal is zero for minus/plus infinity we define the condition
    ###for the first lower bound and the last upper bound which are minus infinity and plus infinity## 
    ### For the first lower bound here is the condition: if the result based on the first term in the summation of Eq(4) which 
    ### is the frequency of first interval*ln(pnorm(the standardized format of the upper bound)) #as the lower bound is minus inf
    a <- -1e+6} 
  else{
      a <- freq[1]*log(pnorm(theta[2]*TL[2]- theta[1]))} 
  #print(a)
  if( (1-pnorm(theta[2]*TL[m]-theta[1])) < 1e-16  ) {
    ### For the last upper bound here is the condition: if the result based on the secind term in the summation of Eq(4) which 
    ### is the frequency of last interval*ln(pnorm(the standardized format of the lower bound)) #as the lower bound is inf
    b <- -1e+6 }else{
      b <- freq[m]*log(1-pnorm(theta[2]*TL[m]-theta[1]))}
  #print(b)
  c<-0
  for(i in 2:(m-1)){ #Here we are doing the last summation of the equation, but we define a condition to avoid zero and non-solvable
    ### results
    #print(i)
    #print(freq[i])
    if ( (pnorm(theta[2]*TL[i+1]-theta[1]) - pnorm(theta[2]*TL[i]-theta[1])) < 1e-16 ){
      c <- c -1e+6 }else{
        c <- c + freq[i]*(log( pnorm(theta[2]*TL[i+1]-theta[1]) - pnorm(theta[2]*TL[i]-theta[1] ))) }
    #print(c)
    #print(TL[i])
  }
  L <- -(a+b+c) ### We add a, b and c, which are the three summations in the equation and assign the result to the variable L
  return(L)### Return the result of the log likelihood
}



#####FOR SIM DATA 2############################################
###Here the above log-likelihood function should be used on a data set to see the estimates of the parameters
### As we have simulated 30 data sets, we have created two matrices to put the results of mu and sigma on each for each data set
### As 30 data sets was in the array, each time we have to read one of the data sets to see the results, so we are using a for loop
### for the 30 data sets. OF course, if it was only one data set, we do not need to do that.
MUMult2<- rep(0,30)
SDMult2<- rep(0,30)

for(j in 1:30){
  TL<- simdata2[,1,j]
  freq<- simdata2[,3,j]
  ### WE have to use an optimization function (here for this case we are using OPTIM function) to find the MLE estimates 
  ### of the above Log-likelihood. 
  ### Note: For this case, from the set of numerical methods, we have to use L-BFGS-B, as the function was very sensitive 
  ### to the method used.
  ### The arguments of optim are: the initial value of the parameters, fn=Logll (the function that we defined above and
  ### wanted to maximize), TL=TL : the lower bounds of the data set, freq=frequency of the data set
  ###method: explained earlier, lower & upper: are the lower and upper limit of the L-BFGS-B method
  res2 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  #estimsim<- ressim$par
  MUMult2[j]<- res2$par[1]/res2$par[2] ### To see the estimate of mu
  SDMult2[j]<- 1/res2$par[2] ### To see the estimates of sigma
}

#MUMult2
#SDMult2
##############################################################################################
#### RUN the optim function for one data set:
###Galton data###
Galton<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton.csv",header=TRUE) ### The original data###
X<- Galton$Parent ### Select the variable Parent
MP<- mean(X) ### mean of parent
SSP<- ((1/length(X))*sum((X-mean(X))^2)) ### Variance of parent


#Galton<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton.csv",header=TRUE)
Y<- Galton$Children ### Select the variable Children
MeC<- mean(Y)### mean of children
SSC<- ((1/length(Y))*sum((Y-mean(Y))^2)) ### Variance of children

mle1<- c(MP,SSP,MeC,SSC) ### put all the estimates in a vector 
###For Galton Data:Parent###
#### This is  the data set of the grouped data for variable parent
Galton2 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton2.csv",header=TRUE)### Now, this data set is the organized 
### version of the grouped data which should be used for all cases which have three columns: one column for the lower bound 
### of each interval, the second column: is the upper bound of each interval, the third column: the frequencies of each interval
Galton2[1,1]<- -Inf ### we set the first lower bound to -inf
Galton2[11,2]<- Inf ### we set the last upper bound to inf
#Galton2
#### This is  the data set of the grouped data for variable Children
Galton3 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton3.csv",header=TRUE)

Galton3[1,1]<- -Inf ### we set the first lower bound to -inf
Galton3[14,2]<- Inf ### we set the last upper bound to inf
#Galton3

thetaintP<- c(MP,SSP) ### assign the thetaintP for the estimates of the mle ignoring grouping for Parent
blP<- Galton2$TL ### rename and assign the lower bound of the intervals of Parent data to blP
buP<- Galton2$TU ### rename and assign the upper bound of the intervals of Parent data to buP
FrP<- Galton2$Freq ### rename and assign the frequency of the intervals of Parent data to FrP


thetaintC<- c(MeC,SSC) ### assign the thetaintC for the estimates of the mle ignoring grouping for Children
blC<- Galton3$TL ### rename and assign the lower bound of the intervals of Children data to blC
buC<- Galton3$TU ### rename and assign the upper bound of the intervals to buC
FrC<- Galton3$Freq ### rename and assign the frequency of the intervals to FrC

res <- optim(c((68/2),(1/2)),fn=Logll,TL=blP,freq=FrP,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))### 
### WE have to use an optimization function (here for this case we are using OPTIM function) to find the MLE estimates 
### of the above Log-likelihood. 
### Note: For this case, from the set of numerical methods, we have to use L-BFGS-B, as the function was very sensitive 
### to the method used.
### The arguments of optim are: the initial value of the parameters, fn=Logll (the function that we defined above and
### wanted to maximize), TL=blP : the lower bounds of the data set, freq=frequency of the data set (for parent data FrP)
###method: explained earlier, lower & upper: are the lower and upper limit of the L-BFGS-B method


estimate<- res$par ### To see the estimate of the parameters
muP1<- round(res$par[1]/res$par[2],4) ### To see the 4 digit rounded estimate of mu parent
SP1<- round(1/res$par[2],4) ### To see the 4 digit rounded estimate of sigma parent


res1 <- optim(c((68/2),(1/2)),fn=Logll,TL=blC,freq=FrC,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2))
### WE have to use an optimization function (here for this case we are using OPTIM function) to find the MLE estimates 
### of the above Log-likelihood. 
### Note: For this case, from the set of numerical methods, we have to use L-BFGS-B, as the function was very sensitive 
### to the method used.
### The arguments of optim are: the initial value of the parameters, fn=Logll (the function that we defined above and
### wanted to maximize), TL=blC : the lower bounds of the data set, freq=frequency of the data set (for parent data FrC)
###method: explained earlier, lower & upper: are the lower and upper limit of the L-BFGS-B method
estimate1<- res1$par ### To see the estimate of the parameters
muC1<- res1$par[1]/res1$par[2] ### To see the 4 digit rounded estimate of mu children
SC1<- 1/res1$par[2]  ### To see the 4 digit rounded estimate of sigma children

Mleexc<- c(muP1,(SP1)^2,muC1,(SC1)^2)