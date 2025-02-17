### EM ALGORITHM For simulated data and Galton data###
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
####Iterating MAXIMAZATION STEP:######
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
output1000m15EM<- EM(bl=simdata1000m15[,1,i],bu=simdata1000m15[,2,i],Freq=simdata1000m15[,3,i],theta_init=c(67,2),
                        maxit = 1000,tol1=1e-3,tol2=1e-4) 


output1000m15EM
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
                       maxit = 1000,tol1=1e-3,tol2=1e-3) 
}
#############################################################################
output300m30EM<- matrix(rep(0,2*500),ncol=2) 
colnames(output300m30EM)<- c("mean","var") 


for(i in 1:500){ 
  output300m30EM[i,]<- EM(bl=simdata300m30[,1,i],bu=simdata300m30[,2,i],Freq=simdata300m30[,3,i],theta_init=c(67,2),
                       maxit = 1000,tol1=1e-3,tol2=1e-3) 
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
simdata1000m30



#############################################################################

#############################################################################
out50m8EM<- write.csv( output50m8EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out50m8EM.csv")
out50m15EM<- write.csv( output50m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out50m15EM.csv")
out50m30EM<- write.csv( output50m30EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out50m30EM.csv")
out100m8EM<- write.csv( output100m8EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out100m8EM.csv")
out100m15EM<- write.csv( output100m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out100m15EM.csv")
out100m30EM<- write.csv( output100m30EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out100m30EM.csv")
out300m8EM<- write.csv( output300m8EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out300m8EM.csv")
out300m15EM<- write.csv( output300m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out300m15EM.csv")
out300m30EM<- write.csv( output300m30EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out300m30EM.csv")
out600m8EM<- write.csv( output600m8EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out600m8EM.csv")
out600m15EM<- write.csv( output600m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out600m15EM.csv")
out600m30EM<- write.csv( output600m30EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out600m30EM.csv")
out1000m8EM<- write.csv( output1000m8EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out1000m8EM.csv")
out1000m15EM<- write.csv( output1000m15EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out1000m15EM.csv")
out1000m30EM<- write.csv( output1000m30EM,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out1000m30EM.csv")

read.csv("C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/EM/out50m8EM.csv")
