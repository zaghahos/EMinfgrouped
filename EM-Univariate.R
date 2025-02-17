

## In this file we want to use EM algorithm method (section 2.1.3) for estimating the parameters 
### We have to do it in E-step & M-step: For E-step We define two functions in R, one called Mest for estimates of mu and
### the other SSest for etimate of sigma. For M-step: which is a maximaziation step, I define a function called EM.

### Calculation of Updates mu###

###E-Step For Mean Estimate###
Mest<- function(theta,bl,bu,Freq){ ### Arguments are: theta (parameters), bl: lower bound of the intervals, 
  ### bu: upper bound of the intervals, Freq: frequencies over each interval
  Aj<- rep(0,length(bl)) ### to produce a vector for the mu's
  astar<- rep(0,length(bl)) ### to produce a vector of the standardized value of the lower bound of the intervals 
  bstar<- rep(0,length(bl)) ### to produce a vector of the standardized value of the upper bound of the intervals 
  for(i in 1:length(bl)){ ### the loop for doing the standardization over all the intervals
    bstar[i]<- (bu[i]-theta[1])/theta[2] ###  Standardized upper bound of the intervals
    astar[i]<- (bl[i]-theta[1])/theta[2] ###  Standardized lower bound of the intervals
    
  }
  
  for(i in 1:length(bl)){ ### the expectations in the equation (5) are calculated for all the intervals, The method of 
    ### how the expectation is calculated could be found in appendix B
    
    Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    
  }
  M<- sum(Aj*Freq)/sum(Freq) ### Now the parameters are calculated using formula (5) and assign it to the variable M 
  #return(Aj)
  return(M) ### return M which are the estimates of M in E-step
  
}

############################################
#############################################
### Calculation of Updates Sigma###
###E-Step For Variance Estimate ###
SSest<- function(theta,bl,bu,muupdate,Freq){ ### Arguments are: theta (parameters), bl: lower bound of the intervals, 
  ### bu: upper bound of the intervals,muupdate: which is the updated estimates of mu which calculated in the above 
  ### function (So this is calling a function in another function)
  ###Freq: frequencies over each interval
  Bj<- rep(0,length(bl)) ### to produce a vector for the sigma's
  bstar<- rep(0,length(bl)) ### to produce a vector of the standardized value of the upper bound of the intervals
  astar<- rep(0,length(bl)) ### to produce a vector of the standardized value of the lower bound of the intervals
  
  for(i in 1:length(bl)){### the loop for doing the standardization over all the intervals
    bstar[i]<- (bu[i]-theta[1])/theta[2] ###  Standardized upper bound of the intervals
    astar[i]<- (bl[i]-theta[1])/theta[2] ###  Standardized lower bound of the intervals
    
  }
  
  astar[1]<- -1000 ### Based on the formula of the resulted expectations in appendix B, for the first value of the vector 
  ### of lower bound (which equals minus inf) we have to put a large value to avoid NaN on the results
  bstar[length(bl)]<- 1000 ### Based on the formula of the resulted expectations in appendix B, for the last value of the vector 
  ### of upper bound (which equals inf) we have to put a large value to avoid NaN on the results
  
  
  for(i in 1:length(bl)){ ### the expectations in the equation (6) are calculated for all the intervals, The method of 
    ### how the expectation is calculated could be found in appendix B
    Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i]))))
    
  }
  SS<- sum(Bj*Freq)/sum(Freq) ### Now the parameters are calculated using formula (6) and assign it to the variable SS
  #S<- sqrt(SS)
  return(SS) ### return SS which are the estimates of sigma in E-step
}

#####################################################
##########################################################
####MAXIMAZATION STEP:######
### This is the maximization step of the EM algorithm (M-step) which I have defined it using the function EM
EM<- function(bl,bu,Freq,theta_init,maxit=1000,tol1=1e-3,tol2=1e-4){ ### Arguments of the function: 
  ### bl: lower bound of the intervals, bu: upper bound of the intervals, Freq: Frequency over the intervals,
  ### theta_init: the initial value of the parameter, maxit=1000: the maximum number of iteration of the EM algorithm,
  ### tol1=1e-3: the stopping criteria for updating mu, tol2=1e-4: the stopping criteria for updating sigma
  flag<- 0 ### The value of flag (which will be used in the condition of the stopping rule), it is set to zero
  Mu_cur<- theta_init[1] ### Assign the current value of the mu equal to the mu from the theta initial
  S_cur<- theta_init[2]  ### Assign the current value of the sigma equal to the sigma from the theta initial
  
  for (i in 1:maxit){ ### this is the updating process which should be done from 1 up to the maximum iteration number (1000)
    #print(paste("Iteration number=", i))
    cur<- c(Mu_cur,S_cur) ### Set the cur values of the parameter from the above (Mu_cur & S_cur)
    
    Munew<- Mest(theta=cur,bl,bu,Freq) ### etimate the new (updated version) mu using the Mest function from the E-step above 
    
    SSnew<- SSest(theta=cur,bl,bu,
                  muupdate=Mest(theta=cur,bl,bu,Freq) ,Freq) ### etimate the new (updated version) sigma using the SSest 
     ###function from the E-step above 
    
    
    Mu_new<- Munew ### set the new (upadeted) values of mu calculated two lines above 
    S_new<- sqrt(SSnew)  ### set the new (upadeted) values of sigma calculated two lines above 
    new_step<- c(Mu_new,S_new) ### put the new values of mu and sigma in a vector called new_step
    
    if(abs(cur[1]-new_step[1])<tol1 & abs(cur[2]-new_step[2])<tol2){flag<-1 ;break} ### This is the condition for updating 
    ### the parameters, for each parameter, the abs value of the differences between the new and the current are calculated, 
    ### if they are less than the stopping criteria, if they are not match the condition, the updating on 
    ### the iterations will be continued
    ###until the stopping rule is meet, which means the difference between the current and new value are less than tol1 & tol2 
    ###and then the final parameters are updated to the last new values and we assign them the vector called updateres 
    
    Mu_cur<- Mu_new
    S_cur<- S_new
  }
  if(!flag) warning("Didn't Converge \n") ### The condition when the iteration in not convergence which return the warning
  #list(paste("updated Mu of the EM=" , Mu_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #    paste("iteration number=",i))
  updateres<- c(Mu_cur,(S_cur)^2) ### assigned the updated resultes to this vector
  return(updateres) ### Return the updated result
}

###########################################################
###########################################################
########################################################
########################################################
### 30 SIMULATION RESULT, n=50 & 10 Classes####
### Now we have to run the above EM function, on a data set to see the result of the estimates of the parameters using EM algorithm
### As in this file we have done it on the 30 simulated data sets, we have to create a matrix with 2 colums (one for mu and 
### one for sigma) and 30 rows (one for the results of each data set) 
output2<- matrix(rep(0,2*30),ncol=2) ### producing the matrix of the result of the parameters for 30 simulated data sets
colnames(output2)<- c("mean","std") ### naming the column of the matrix


for(i in 1:30){ ### do the iteration of EM method for each data set 
  output2[i,]<- EM(bl=simdata2[,1,i],bu=simdata2[,2,i],Freq=simdata2[,3,i],theta_init=c(67,2),
                   maxit = 1000,tol1=1e-3,tol2=1e-4) ### Running EM function on the data sets- here we set the initial values
  ### theta_init=c(67,2), they could be any numbers but preferably closer values to the MLE estimates
  
}
#print (output2)
#simdata1    
####################################################################
##### Note: TO see the results on one data set, we can see the file of the Galton-univariate code  

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
outParent<- EM(bl=blP,bu=buP,Freq=FrP,theta_init=thetaintP,maxit=1000,tol1=1e-3,tol2=1e-4)### Run EM algorithm for parent 
###using the blP,buP,Freq,initial value of the theta for parent found above thetaintP, maxit=1000: max of iteration, and 
### tol1 & tol2: the stopping rule criteria


outChildren<- EM(bl=blC,bu=buC,Freq=FrC,theta_init=thetaintC,maxit=1000,tol1=1e-3,tol2=1e-4)### Run EM algorithm for children 
###using the blP,buP,Freq,initial value of the theta for parent found above thetaintP, maxit=1000: max of iteration, and 
### tol1 & tol2: the stopping rule criteria

Emres<- c(outParent,outChildren) ### put all the EM estimates for Parent and Children in a vector called Emres