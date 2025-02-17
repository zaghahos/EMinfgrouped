
## In this file we want to use MCEM algorithm method (section 2.1.4) for estimating the parameters 
### We have to do it in E-step & M-step: For E-step We define three functions in R, one called ZMCEM: to simulate the new samples
###MuMCEM for estimates of mu and
### the other SigmaMCEM for etimate of sigma. For M-step: which is a maximaziation step, I define a function called MCEM.

library(truncnorm) ### We have to use this library for the simulations from truncated normal over each separate interval
###1. Update Both Mu & Sigma based on the Formulas###

###Calculation of updates for Mu###
##E-Step for Mu estimation: Simulating Z's###
ZMCEM<- function(theta,data){ ### The function to generate the simulations from the function mentioned in part 2.1.4 for 
  ### each interval with the arguments: theta: the current state of the parameters, and data 
  k<- 1000 ### set k=1000, the number of samples simulated
  sim<- matrix(rep(0,k*nrow(data)),ncol=nrow(data)) ### define a matrix to put the simulated values on that
  #sim<- matrix(rep(0,k*3),ncol=3)
  
  for(i in 1 :nrow(data)) ### use for loop to simulate the samples over each intervals
    sim[,i]<- rtruncnorm(k,a=data[i,1],b=data[i,2],mean=theta[1],sd=theta[2]) ### Simulate data over each specific intervals 
  ### and assign it to the matrix sim
  
  return(sim) ### return the sim matrix
}  
### E-Step for Mu & Sigma###
### Estimate the parameters of mu & sigma in the E-step using the defined function 
MuMCEM<- function(data,simZ){ ### Function to calculate estimate of mu, takes two arguments: data and simZ: which is the 
  ### the simulated samples over the intervals from the above function ZMCEM
  n<- sum(data[,3]) ### total number of data
  Z<- colMeans(simZ) ### A vector of the colmeans of the simulated samples
  numerator<- rep(0,nrow(data)) ### define a vector to calculate mulptiplications of the colmeans of ZMCEM matrix and the
  ### frequencies for each separate interval
  for (i in 1:nrow(data)) #### using for loop to calculate the mentioned multiplications over the intervals
    numerator[i]<- data[i,3]*Z[i]
  sum(numerator) ### We have to sum the result of the multiplications over all the intervals (the numerator on the formula 7)
  MuN<- (1/n)*sum(numerator) ### divide the summation by the total number (fromula 7) and assign it to the value MuN
  return(MuN) ### return the MuN (updated mu from the E-step)
}

################################################
sigmaMCEM<- function(data,simZZ,mupd){ ### Function to calculate estimate of sigma, takes three arguments: 
  ### data and simZ: which is the the simulated samples over the intervals from the above function ZMCEM,
  ### mupd: the updated mu from the first step of E-step
  n<- sum(data[,3]) ### total number of data
  ZZ<- simZZ ### The simulated samples and assign them to the matrix ZZ
  #print(ZZ)
  NewZ<- (ZZ-mupd)^2 ### claculate the squar of (simulated samples-muupdate)^2 from the formula (8) and assign it to the 
  ### matrix NewZ
  #print(NewZ)
  SZNEW<- colMeans(NewZ) ### Calculate the colmeans of the matrix NewZ (the part after n_i in the formula 8) and assign it to 
  ### SZNEW
  #print(SZNEW)
  numerator<- rep(0,nrow(data)) ### define a vector to calculate mulptiplications of the colmeans of NewZ matrix (SZNEW) 
  ###and the frequencies for each separate interval 
  for (i in 1:nrow(data))#### using for loop to calculate the mentioned multiplications over the intervals
    numerator[i]<- data[i,3]*SZNEW[i] ### We have to sum the result of the multiplications over all the intervals (the numerator on the formula 8)
  
  sigmaNN<- (1/n)*sum(numerator) ### divide the summation by the total number (fromula 8) and assign it to the value sigmaNN
  sig<- sqrt(sigmaNN) ### take square from the sigmaNN and assign it to the value sig
  
  # Z<- colMeans(ZZ)
  #print(Z)
  
  
  #zsqr<- (Z-mupd)^2
  
  #numerator<- rep(0,nrow(data))
  #for (i in 1:nrow(data))
  # numerator[i]<- data[i,3]*zsqr[i]
  
  #sigmaNN<- (1/n)*sum(numerator)
  #sig<- sqrt(sigmaNN)
  return(sig) ### return the sig (updated sigma from the E-step)
}

##########################################
###MONTE CARLO EM###
##########################################################
### This is the maximization step of the MCEM algorithm (M-step) which I have defined it using the function MCEM

MCEM<- function(data,theta_init,maxit=1000,tol1=1e-2,tol2=1e-3){ ### Arguments of the function: 
  ### data,theta_init: the initial value of the parameter, maxit=1000: the maximum number of iteration of the MCEM algorithm,
  ### tol1=1e-2: the stopping criteria for updating mu, tol2=1e-3: the stopping criteria for updating sigma
  flag<- 0 ### The value of flag (which will be used in the condition of the stopping rule), it is set to zero
  Mu_cur<- theta_init[1] ### Assign the current value of the mu equal to the mu from the theta initial
  S_cur<- theta_init[2]  ### Assign the current value of the sigma equal to the sigma from the theta initial
  #theta_cur<- c(Mu_cur,S_cur)
  iter<- rep(0,maxit) ### 
  Svec<- rep(0,maxit)
  Mvec<- rep(0,maxit)
  for (i in 1:maxit){ ### this is the updating process which should be done from 1 up to the maximum iteration number (1000)
    #print(paste("Iteration number=", i))
    cur<- c(Mu_cur,S_cur) ### Set the cur values of the parameter from the above (Mu_cur & S_cur)
    Munew<- MuMCEM(data=mydat,simZ=ZMCEM(theta=cur,data=mydat)) ### etimate the new (updated version) mu using the MuMCEM 
    ###function from the E-step above 
    
    Snew<- sigmaMCEM(data=mydat,simZZ=ZMCEM(theta=cur,data=mydat),
                     mupd=MuMCEM(data=mydat,simZ=ZMCEM(theta=cur,data=mydat))) ### etimate the new (updated version) sigma 
                                   ###using the Mest function from the E-step above 
    
    
    
    Mu_new<- Munew ### set the new (upadeted) values of mu calculated two lines above 
    S_new<- Snew ### set the new (upadeted) values of sigma calculated two lines above 
    
    new_step<- c(Mu_new,S_new) ### put the new values of mu and sigma in a vector called new_step
    
    
    if(abs(cur[1]-new_step[1])<tol1 & abs(cur[2]-new_step[2])<tol2){flag<-1 ;break} ### This is the condition for updating 
    ### the parameters, for each parameter, the abs value of the differences between the new and the current are calculated, 
    ### if they are less than the stopping criteria, if they are not match the condition, the updating on 
    ### the iterations will be continued
    ###until the stopping rule is meet, which means the difference between the current and new value are less than tol1 & tol2 
    ###and then the final parameters are updated to the last new values and we assign them the vector called update
    
    
    Mu_cur<- Mu_new 
    
    S_cur<- S_new
    
    iter[i]<- i
    Svec[i]<- S_new
    Mvec[i]<- Mu_new 
    
  }
  
  if(!flag) warning("Didn't Converge \n") ### The condition when the iteration in not convergence which return the warning
  
  #list(paste("updated Mu of the EM=" , Mu_cur),paste("Updated Sigma of the EM=" ,S_cur),
   #    paste("iteration number=",i)) 
  update=c(Mu_cur,(S_cur)^2) ### assigned the update resultes to this vector
  return(update) ### Return the update result
}

##################################################################
############################################################################################
### 30 SIMULATION RESULT, n=50 & 10 Classes####
### Now we have to run the above MCEM function, on a data set to see the result of the estimates of the parameters using MCEM
###algorithm As in this file we have done it on the 30 simulated data sets, we have to create a matrix with 2 colums 
###(one for mu and one for sigma) and 30 rows (one for the results of each data set) 
outputMCEM2<- matrix(rep(0,2*30),ncol=2) ### producing the matrix of the result of the parameters for 30 simulated data sets
colnames(outputMCEM2)<- c("mean","std")  ### naming the column of the matrix

for(i in 1:30){ ### do the iteration of EM method for each data set
  #i<- 1 
  mydat<- simdata2[,,i]
  outputMCEM2[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) ### Running MCEM function on the 
  ###data sets- here we set the initial values theta_init=c(67,2), they could be any numbers but preferably 
  ###closer values to the MLE estimates
  
}

########################################################################################################
### To run MCEM on Galton data (only one data set) first we have to read the data which is in the form of ignoring the
### grouping (see section 2.1.1), and the results of the estimates are the mle ignoring the grouping
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


##GALTON DATA: Parent##
###For Galton Data###
#### This is  the data set of the grouped data for variable parent
Galton2 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton2.csv",header=TRUE) ### Now, this data set is the organized 
### version of the grouped data which should be used for all cases which have three columns: one column for the lower bound 
### of each interval, the second column: is the upper bound of each interval, the third column: the frequencies of each interval
#Galton2
BL <- Galton2$TL ### rename and assign the lower bound of the intervals of Parent data to BL
BL[1] <- -Inf ### we set the first lower bound to -inf
BU <- Galton2$TU ### rename and assign the upper bound of the intervals of Parent data to BU
BU[11] <- Inf ### we set the last upper bound to inf

Freq <- Galton2$Freq ### rename and assign the frequency of the intervals to Freq
GaltonDat<- cbind(BL,BU,Freq) ### put all the columns of BL,BU, and Freq in the columnbind matrix
GD<- as.data.frame(GaltonDat)### change the matrix to the data frame 

mydat<- GD ### Rename the data frame to a name used in the EM /MCEM function (mydat)
#theta_init<- c(Me,SD)
MCEMParent<- MCEM(data=mydat,theta_init=thetaintP,maxit = 1000,tol1=1e-2,tol2=1e-3) ### Run MCEM algorithm for parent using the
###data=mydat,initial value of the theta for parent found above thetaintP, maxit=1000: max of iteration, and 
### tol1 & tol2: the stopping rule criteria


#############################################################
### GALTON DATA: Children###
#### This is  the data set of the grouped data for variable Children
Galton3 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton3.csv",header=TRUE)
#Galton3
#dim(Galton3)
BL1 <- Galton3$TL ### rename and assign the lower bound of the intervals of children data to BL1
BL1[1] <- -Inf ### we set the first lower bound to -inf
BU1 <- Galton3$TU  ### rename and assign the upper bound of the intervals of children data to BU1
BU1[14] <- Inf ### we set the last upper bound to inf

Freqc <- Galton3$Freq ### rename and assign the frequency of the intervals to Freqc
GaltonDatC<- cbind(BL1,BU1,Freqc)### put all the columns of BL1,BU1, and Freqc in the columnbind matrix
GDC<- as.data.frame(GaltonDatC) ### change the matrix to the data frame 

mydat<- GDC ### Rename the data frame to a name used in the EM /MCEM function (mydat)
MCEMChildren<- MCEM(data=mydat,theta_init=thetaintC,maxit = 1000,tol1=1e-2,tol2=1e-3) ### Run MCEM algorithm for children using the
###data=mydat,initial value of the theta for children found above thetaintC, maxit=1000: max of iteration, and 
### tol1 & tol2: the stopping rule criteria
############################################################################################