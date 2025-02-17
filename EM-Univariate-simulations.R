

## In this file we want to use EM algorithm method (section 2.1.3) for estimating the parameters 
### We have to do it in E-step & M-step: For E-step We define two functions in R, one called Mest for estimates of mu and
### the other SSest for etimate of sigma. For M-step: which is a maximaziation step, I define a function called EM.

### Calculation of Updates mu###

###E-Step For Mean Estimate###
Mest<- function(theta,bl,bu,Freq){ 
  Aj<- rep(0,length(bl)) 
  astar<- rep(0,length(bl)) 
  bstar<- rep(0,length(bl))  
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2] 
    astar[i]<- (bl[i]-theta[1])/theta[2] 
    
  }
  
  for(i in 1:length(bl)){  

    Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    
  }
  M<- sum(Aj*Freq)/sum(Freq) 
  return(M) 
  
}

############################################
#############################################
### Calculation of Updates Sigma###
###E-Step For Variance Estimate ###
SSest<- function(theta,bl,bu,muupdate,Freq){ 
  Bj<- rep(0,length(bl)) 
  bstar<- rep(0,length(bl)) 
  astar<- rep(0,length(bl)) 
  
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2] 
    astar[i]<- (bl[i]-theta[1])/theta[2] 
    
  }
  
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  
  
  for(i in 1:length(bl)){ 
    Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i]))))
    
  }
  SS<- sum(Bj*Freq)/sum(Freq) 
  #S<- sqrt(SS)
  return(SS) 
}

#####################################################
##########################################################
####MAXIMAZATION STEP:######
### This is the maximization step of the EM algorithm (M-step) which I have defined it using the function EM
EM<- function(bl,bu,Freq,theta_init,maxit=1000,tol1=1e-3,tol2=1e-4){ 
  flag<- 0 
  Mu_cur<- theta_init[1] 
  S_cur<- theta_init[2]  
  
  for (i in 1:maxit){ 
    cur<- c(Mu_cur,S_cur) 
    
    Munew<- Mest(theta=cur,bl,bu,Freq) 
    
    SSnew<- SSest(theta=cur,bl,bu,
                  muupdate=Mest(theta=cur,bl,bu,Freq) ,Freq) 
    
    Mu_new<- Munew  
    S_new<- sqrt(SSnew)  
    new_step<- c(Mu_new,S_new) 
    
    if(abs(cur[1]-new_step[1])<tol1 & abs(cur[2]-new_step[2])<tol2){flag<-1 ;break} 
    
    Mu_cur<- Mu_new
    S_cur<- S_new
  }
  if(!flag) warning("Didn't Converge \n") 
  #list(paste("updated Mu of the EM=" , Mu_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #    paste("iteration number=",i))
  updateres<- c(Mu_cur,(S_cur)^2)
  return(updateres) 
}

###########################################################
###########################################################
########################################################
########################################################
### 500 SIMULATION RESULT###############
 
output50m8EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output50m8EM)<- c("mean","var") 


for(i in 1:500){ 
  output50m8EM[i,]<- EM(bl=simdata50m8[,1,i],bu=simdata50m8[,2,i],Freq=simdata50m8[,3,i],theta_init=c(67,2),
                   maxit = 1000,tol1=1e-3,tol2=1e-4) 
}

###############################################################################
output100m8EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output100m8EM)<- c("mean","var") 


for(i in 1:500){ 
  output100m8EM[i,]<- EM(bl=simdata100m8[,1,i],bu=simdata100m8[,2,i],Freq=simdata100m8[,3,i],theta_init=c(67,2),
                      maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output300m8EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output300m8EM)<- c("mean","var") 


for(i in 1:500){ 
  output300m8EM[i,]<- EM(bl=simdata300m8[,1,i],bu=simdata300m8[,2,i],Freq=simdata300m8[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output600m8EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output600m8EM)<- c("mean","var") 


for(i in 1:500){ 
  output600m8EM[i,]<- EM(bl=simdata600m8[,1,i],bu=simdata600m8[,2,i],Freq=simdata600m8[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output1000m8EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output1000m8EM)<- c("mean","var") 


for(i in 1:500){ 
  output1000m8EM[i,]<- EM(bl=simdata1000m8[,1,i],bu=simdata1000m8[,2,i],Freq=simdata1000m8[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output50m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output50m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output50m15EM[i,]<- EM(bl=simdata50m15[,1,i],bu=simdata50m15[,2,i],Freq=simdata50m15[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output100m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output100m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output100m15EM[i,]<- EM(bl=simdata100m15[,1,i],bu=simdata100m15[,2,i],Freq=simdata100m15[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output300m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output300m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output300m15EM[i,]<- EM(bl=simdata300m15[,1,i],bu=simdata300m15[,2,i],Freq=simdata300m15[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output600m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output600m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output600m15EM[i,]<- EM(bl=simdata600m15[,1,i],bu=simdata600m15[,2,i],Freq=simdata600m15[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output1000m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output1000m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output1000m15EM[i,]<- EM(bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output50m30EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output50m30EM)<- c("mean","var") 


for(i in 1:500){ 
  output50m30EM[i,]<- EM(bl=simdata50m30[,1,i],bu=simdata50m30[,2,i],Freq=simdata50m30[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output100m30EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output100m30EM)<- c("mean","var") 


for(i in 1:500){ 
  output100m30EM[i,]<- EM(bl=simdata100m30[,1,i],bu=simdata100m30[,2,i],Freq=simdata100m30[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output300m30EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output300m30EM)<- c("mean","var") 


for(i in 1:500){ 
  output300m30EM[i,]<- EM(bl=simdata300m30[,1,i],bu=simdata300m30[,2,i],Freq=simdata300m30[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output600m30EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output600m30EM)<- c("mean","var") 


for(i in 1:500){ 
  output600m30EM[i,]<- EM(bl=simdata600m30[,1,i],bu=simdata600m30[,2,i],Freq=simdata600m30[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
output1000m30EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output1000m30EM)<- c("mean","var") 


for(i in 1:500){ 
  output1000m30EM[i,]<- EM(bl=simdata1000m30[,1,i],bu=simdata1000m30[,2,i],Freq=simdata1000m30[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################

#############################################################################
out50m8EM<- write.csv( output50m8EM,"C:/Users/sh_za/Desktop/Results/EM Results/out50m8EM.csv")
out50m15EM<- write.csv( output50m15EM,"C:/Users/sh_za/Desktop/Results/EM Results/out50m15EM.csv")
out50m30EM<- write.csv( output50m30EM,"C:/Users/sh_za/Desktop/Results/EM Results/out50m30EM.csv")
out100m8EM<- write.csv( output100m8EM,"C:/Users/sh_za/Desktop/Results/EM Results/out100m8EM.csv")
out100m15EM<- write.csv( output100m15EM,"C:/Users/sh_za/Desktop/Results/EM Results/out100m15EM.csv")
out100m30EM<- write.csv( output100m30EM,"C:/Users/sh_za/Desktop/Results/EM Results/out100m30EM.csv")
out300m8EM<- write.csv( output300m8EM,"C:/Users/sh_za/Desktop/Results/EM Results/out300m8EM.csv")
out300m15EM<- write.csv( output300m15EM,"C:/Users/sh_za/Desktop/Results/EM Results/out300m15EM.csv")
out300m30EM<- write.csv( output300m30EM,"C:/Users/sh_za/Desktop/Results/EM Results/out300m30EM.csv")
out600m8EM<- write.csv( output600m8EM,"C:/Users/sh_za/Desktop/Results/EM Results/out600m8EM.csv")
out600m15EM<- write.csv( output600m15EM,"C:/Users/sh_za/Desktop/Results/EM Results/out600m15EM.csv")
out600m30EM<- write.csv( output600m30EM,"C:/Users/sh_za/Desktop/Results/EM Results/out600m30EM.csv")
out1000m8EM<- write.csv( output1000m8EM,"C:/Users/sh_za/Desktop/Results/EM Results/out1000m8EM.csv")
out1000m15EM<- write.csv( output1000m15EM,"C:/Users/sh_za/Desktop/Results/EM Results/out1000m15EM.csv")
out1000m30EM<- write.csv( output1000m30EM,"C:/Users/sh_za/Desktop/Results/EM Results/out1000m30EM.csv")





#simdata50m8<- write.csv(simulation50m8,"C:/Users/sh_za/Desktop/Results/Univariate Simulated Data/simdata50m8.csv")

#simulateddata<- write.csv(clusters,"C:/Users/sh_za/OneDrive/Desktop/Simulated Data/Mixture Poisson/simulateddata.csv")


#Galton2 <- read.csv("C:/Users/sh_za/Desktop/Galton Data/Galton2.csv",header=TRUE) ### Now, this data set is the organized 



####################################################################
##### Note: TO see the results on one data set, we can see the file of the Galton-univariate code  

##GALTON DATA: Parent##
###For Galton Data###
#### This is  the data set of the grouped data for variable parent
Galton2 <- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton2.csv",header=TRUE) ### Now, this data set is the organized 
output100m8<- matrix(rep(0,2*500),ncol=2) 
colnames(output2)<- c("mean","std") 


for(i in 1:500){ 
  output100m8[i,]<- EM(bl=simdata2[,1,i],bu=simdata2[,2,i],Freq=simdata2[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
#############################################################################
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