### MCEM ALGORITHM For simulated data and Galton data###

library(truncnorm) 
###Calculation of updates for Mu###
##E-Step for Mu estimation: Simulating Z's###
ZMCEM<- function(theta,data){ 
  k<- 1000 
  sim<- matrix(rep(0,k*nrow(data)),ncol=nrow(data)) 
  #sim<- matrix(rep(0,k*3),ncol=3)
  
  for(i in 1 :nrow(data)) {
    sim[,i]<- rtruncnorm(k,a=data[i,1],b=data[i,2],mean=theta[1],sd=theta[2]) }
  return(sim) 
}  
### E-Step for Mu & Sigma###
MuMCEM<- function(data,simZ){ 
  n<- sum(data[,3])
  Z<- colMeans(simZ)
  numerator<- rep(0,nrow(data)) 
  for (i in 1:nrow(data)) {
    numerator[i]<- data[i,3]*Z[i]}
  sum(numerator) 
  MuN<- (1/n)*sum(numerator) 
  return(MuN) 
}

################################################
sigmaMCEM<- function(data,simZZ,mupd){ 
  n<- sum(data[,3])
  ZZ<- simZZ 
  #print(ZZ)
  NewZ<- (ZZ-mupd)^2 
  #print(NewZ)
  SZNEW<- colMeans(NewZ) 
  #print(SZNEW)
  numerator<- rep(0,nrow(data))  
  for (i in 1:nrow(data)){
    numerator[i]<- data[i,3]*SZNEW[i] }
  
  sigmaNN<- (1/n)*sum(numerator) 
  sig<- sqrt(sigmaNN) 
  
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
### This is the Iteration of maximization step of the MCEM algorithm (M-step) which I have defined it using the function MCEM

MCEM<- function(data,theta_init,maxit=1000,tol1=1e-2,tol2=1e-3){ 
  flag<- 0 
  Mu_cur<- theta_init[1] 
  S_cur<- theta_init[2]  
  #theta_cur<- c(Mu_cur,S_cur)
  iter<- rep(0,maxit)  
  Svec<- rep(0,maxit)
  Mvec<- rep(0,maxit)
  for (i in 1:maxit){ 
    #print(paste("Iteration number=", i))
    cur<- c(Mu_cur,S_cur) 
    Munew<- MuMCEM(data=mydat,simZ=ZMCEM(theta=cur,data=mydat))  

    Snew<- sigmaMCEM(data=mydat,simZZ=ZMCEM(theta=cur,data=mydat),
                     mupd=MuMCEM(data=mydat,simZ=ZMCEM(theta=cur,data=mydat))) 
    
    
    Mu_new<- Munew 
    S_new<- Snew  
    
    new_step<- c(Mu_new,S_new) 
    
    
    if(abs(cur[1]-new_step[1])<tol1 & abs(cur[2]-new_step[2])<tol2){flag<-1 ;break} 
    
    Mu_cur<- Mu_new 
    
    S_cur<- S_new
    
    iter[i]<- i
    Svec[i]<- S_new
    Mvec[i]<- Mu_new 
    
  }
  
  if(!flag) warning("Didn't Converge \n") 
  
  #update=list(paste("updated Mu of the EM=" , Mu_cur),paste("Updated Sigma of the EM=" ,S_cur),
   #    paste("iteration number=",i)) 
  update=c(Mu_cur,(S_cur)^2) 
  return(update) 
}

##################################################################
############################################################################################
### 500 SIMULATION RESULT####

##############################################################################################
#rm("mydat")
outputMCEM<- matrix(rep(0,3*50),ncol=3) 
colnames(outputMCEM)<- c("mean","var","iter")  

for(i in 1:50){ 
  #i<- 1 
  mydat<- simdata1000m30[,,i]
  #outputMCEM50m8[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
  outputMCEM[i,]<- unlist(MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3))
}

outputMCEM
##################################################################################################

outputMCEM50m8<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM50m8)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata50m8[,,i]
  outputMCEM50m8[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
  #outputMCEM50m8[i,]<- unlist(MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3))
}

outputMCEM50m8
#unlist(MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3))
#######################################################################
rm("mydat")

outputMCEM50m15<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM50m15)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata50m15[,,i]
  outputMCEM50m15[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM50m30<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM50m30)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata50m30[,,i]
  outputMCEM50m30[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM100m8<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM100m8)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata100m8[,,i]
  outputMCEM100m8[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM100m15<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM100m15)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata100m15[,,i]
  outputMCEM100m15[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM100m30<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM100m30)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata100m30[,,i]
  outputMCEM100m30[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM300m8<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM300m8)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata300m8[,,i]
  outputMCEM300m8[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM300m15<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM300m15)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata300m15[,,i]
  outputMCEM300m15[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM300m30<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM300m30)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata300m30[,,i]
  outputMCEM300m30[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM600m8<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM600m8)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata600m8[,,i]
  outputMCEM600m8[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}

#######################################################################
rm("mydat")

outputMCEM600m15<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM600m15)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata600m15[,,i]
  outputMCEM600m15[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM600m30<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM600m30)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata600m30[,,i]
  outputMCEM600m30[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM1000m8<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM1000m8)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata1000m8[,,i]
  outputMCEM1000m8[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM1000m15<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM1000m15)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata1000m15[,,i]
  outputMCEM1000m15[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}
#######################################################################
rm("mydat")

outputMCEM1000m30<- matrix(rep(0,2*500),ncol=2) 
colnames(outputMCEM1000m30)<- c("mean","var")  

for(i in 1:500){ 
  #i<- 1 
  mydat<- simdata1000m30[,,i]
  outputMCEM1000m30[i,]<- MCEM(data=mydat,theta_init=c(67,2),maxit = 1000,tol1=1e-2,tol2=1e-3) 
}

#############################################################################################################
#############################################################################################################
outMCEM50m8<- write.csv(outputMCEM50m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM50m8.csv")
outMCEM50m15<- write.csv(outputMCEM50m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM50m15.csv")
outMCEM50m30<- write.csv(outputMCEM50m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM50m30.csv")

outMCEM100m8<- write.csv(outputMCEM100m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM100m8.csv")
outMCEM100m15<- write.csv(outputMCEM100m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM100m15.csv")
outMCEM100m30<- write.csv(outputMCEM100m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM100m30.csv")

outMCEM300m8<- write.csv(outputMCEM300m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM300m8.csv")
outMCEM300m15<- write.csv(outputMCEM300m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM300m15.csv")
outMCEM300m30<- write.csv(outputMCEM300m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM300m30.csv")

outMCEM600m8<- write.csv(outputMCEM600m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM600m8.csv")
outMCEM600m15<- write.csv(outputMCEM600m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM600m15.csv")
outMCEM600m30<- write.csv(outputMCEM600m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM600m30.csv")

outMCEM1000m8<- write.csv(outputMCEM1000m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM1000m8.csv")
outMCEM1000m15<- write.csv(outputMCEM1000m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM1000m15.csv")
outMCEM1000m30<- write.csv(outputMCEM1000m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/MCEM/outMCEM1000m30.csv")

##############################################################################################################################

#out50m15MCEM<- write.csv( outputMCEM50m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out50m15MCEM.csv")
#out100m15MCEM<- write.csv( outputMCEM100m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out100m15MCEM.csv")
#out300m15MCEM<- write.csv( outputMCEM300m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out300m15MCEM.csv")
#out600m15MCEM<- write.csv( outputMCEM600m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out600m15MCEM.csv")
#out1000m15MCEM<- write.csv( outputMCEM1000m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out1000m15MCEM.csv")

####################################################################################################################################
