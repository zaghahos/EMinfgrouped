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
  dinom<- NULL
  for(i in 1:length(bl)){ 
    #print(i)
    #print(pnorm(bstar[i]))
    #print(pnorm(astar[i]))
    dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i]))
    #print(dinom[i])
    if(dinom[i]==0) { Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001)}
    else {Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i])}
  }
  #Aj[i]<- theta[1]-theta[2]*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i])     
  #}
  M<- sum(Aj*Freq)/sum(Freq) 
  return(M) 
  
}


############################################)
#############################################
### Calculation of Updates Sigma###
###E-Step For Variance Estimate ###
SSest<- function(theta,bl,bu,muupdate,Freq){ 
  Bj<- rep(0,length(bl)) 
  bstar<- rep(0,length(bl)) 
  astar<- rep(0,length(bl)) 
  
  dinom<- NULL
  
  for(i in 1:length(bl)){
    bstar[i]<- (bu[i]-theta[1])/theta[2] 
    astar[i]<- (bl[i]-theta[1])/theta[2] 
    
  }
  
  astar[1]<- -1000 
  bstar[length(bl)]<- 1000 
  
  
  for(i in 1:length(bl)){ 
    dinom[i]<- (pnorm(bstar[i])-pnorm(astar[i])) 
    if(dinom[i]==0) {Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/0.0001)
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/0.0001))}
    
    else{Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/dinom[i])
    +(muupdate-theta[1])^2+
      (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/dinom[i]))}
    
    
    # Bj[i]<- theta[2]^2*(1-(bstar[i]*dnorm(bstar[i])-astar[i]*dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i])))
    #+(muupdate-theta[1])^2+
    # (2*theta[2]*(muupdate-theta[1])*((dnorm(bstar[i])-dnorm(astar[i]))/(pnorm(bstar[i])-pnorm(astar[i]))))
    
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
    #print(Mu_new)
    S_new<- sqrt(SSnew)  
    #print(S_new)
    new_step<- c(Mu_new,S_new) 
    #print(abs(cur[1]-new_step[1])<tol1)
    #print(abs(cur[2]-new_step[2])<tol2)
    
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
### TEST#####
#output1000m15EM<- EM(bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i],theta_init=c(67,2),
 #                    maxit = 1000,tol1=1e-3,tol2=1e-4) 


#output1000m15EM
########################################################
########################################################


#############################################################################
output50m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output50m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output50m15EM[i,]<- EM(bl=simdata50m15[,1,i],bu=simdata50m15[,2,i],Freq=simdata50m15[,3,i],theta_init=c(67,2),
                         maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
output50m15EM
#############################################################################
output100m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output100m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output100m15EM[i,]<- EM(bl=simdata100m15[,1,i],bu=simdata100m15[,2,i],Freq=simdata100m15[,3,i],theta_init=c(67,2),
                          maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
output100m15EM
#############################################################################
output300m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output300m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output300m15EM[i,]<- EM(bl=simdata300m15[,1,i],bu=simdata300m15[,2,i],Freq=simdata300m15[,3,i],theta_init=c(67,2),
                          maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
output300m15EM
#############################################################################
output600m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output600m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output600m15EM[i,]<- EM(bl=simdata600m15[,1,i],bu=simdata600m15[,2,i],Freq=simdata600m15[,3,i],theta_init=c(67,2),
                          maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
output600m15EM
#############################################################################
output1000m15EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output1000m15EM)<- c("mean","var") 


for(i in 1:500){ 
  output1000m15EM[i,]<- EM(bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i],theta_init=c(67,2),
                           maxit = 1000,tol1=1e-3,tol2=1e-4) 
}
output1000m15EM
#############################################################################
#############################################################################
#############################################################################
out50m15EM<- write.csv( output50m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out50m15EM.csv")
out100m15EM<- write.csv( output100m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out100m15EM.csv")
out300m15EM<- write.csv( output300m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out300m15EM.csv")
out600m15EM<- write.csv( output600m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out600m15EM.csv")
out1000m15EM<- write.csv( output1000m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/out1000m15EM.csv")





#simdata50m8<- write.csv(simulation50m8,"C:/Users/sh_za/Desktop/Results/Univariate Simulated Data/simdata50m8.csv")

#simulateddata<- write.csv(clusters,"C:/Users/sh_za/OneDrive/Desktop/Simulated Data/Mixture Poisson/simulateddata.csv")


#Galton2 <- read.csv("C:/Users/sh_za/Desktop/Galton Data/Galton2.csv",header=TRUE) ### Now, this data set is the organized 



####################################################################
##### Note: TO see the results on one data set, we can see the file of the Galton-univariate code  

##GALTON DATA: Parent##
###For Galton Data###
#### This is  the data set of the grouped data for variable parent
Galton <- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton Data/Galton.csv",header=TRUE)  
Galton

X<- Galton$Parent 
MP<- mean(X)
SSP<- ((1/length(X))*sum((X-mean(X))^2))

Y<- Galton$Children
MeC<- mean(Y)
SSC<- ((1/length(Y))*sum((Y-mean(Y))^2))

mle1<- c(MP,SSP,MeC,SSC)
mle1

Galton2 <- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton Data/Galton2.csv",header=TRUE)
Galton2[1,1]<- -Inf
Galton2[11,2]<- Inf
Galton2

Galton3 <- read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Galton Data/Galton3.csv",header=TRUE)

Galton3[1,1]<- -Inf
Galton3[14,2]<- Inf
Galton3

#############################################################################

BL <- Galton2$TL 
BL
#BL[1] <- -Inf
#BU <- Galton2$TU 
#BU[11] <- Inf 
#BL
#BU

Freq <- Galton2$Freq 
GaltonDat<- cbind(BL,BU,Freq) 
GD<- as.data.frame(GaltonDat)
GD
mydat<- GD ### 
#theta_init<- c(Me,SD)
#MCEMParent<- MCEM(data=mydat,theta_init=thetaintP,maxit = 1000,tol1=1e-2,tol2=1e-3) #


#############################################################
### GALTON DATA: Children###

thetaintP<- c(MP,SSP) ##PARENT
blP<- Galton2$TL 
buP<- Galton2$TU 
FrP<- Galton2$Freq 

#blP
thetaintC<- c(MeC,SSC) ## CHILDREN
blC<- Galton3$TL 
buC<- Galton3$TU 
FrC<- Galton3$Freq 
outParent<- EM(bl=blP,bu=buP,Freq=FrP,theta_init=thetaintP,maxit=1000,tol1=1e-3,tol2=1e-4) 

outChildren<- EM(bl=blC,bu=buC,Freq=FrC,theta_init=thetaintC,maxit=1000,tol1=1e-3,tol2=1e-4)

EM_Result<- c(outParent,outChildren) 
EM_Result
outParent
outChildren
