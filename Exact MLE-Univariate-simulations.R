#### EXACT MLE###########################################
###LOG L FUNCTION###

Logll <- function(TL,freq,theta){ 
  m<- length(TL) 
  
  if( (pnorm(theta[2]*TL[2]-theta[1])) < 1e-16  ){ 
    a <- -1e+6} 
  else{
      a <- freq[1]*log(pnorm(theta[2]*TL[2]- theta[1]))} 
  #print(a)
  if( (1-pnorm(theta[2]*TL[m]-theta[1])) < 1e-16  ) {
    b <- -1e+6 }else{
      b <- freq[m]*log(1-pnorm(theta[2]*TL[m]-theta[1]))}
  #print(b)
  c<-0
  for(i in 2:(m-1)){ 
    #print(i)
    #print(freq[i])
    if ( (pnorm(theta[2]*TL[i+1]-theta[1]) - pnorm(theta[2]*TL[i]-theta[1])) < 1e-16 ){
      c <- c -1e+6 }else{
        c <- c + freq[i]*(log( pnorm(theta[2]*TL[i+1]-theta[1]) - pnorm(theta[2]*TL[i]-theta[1] ))) }
    #print(c)
    #print(TL[i])
  }
  L <- -(a+b+c)
  return(L)
}



#####FOR SIM DATA 2############################################

MLEbinned50m8<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned50m8)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata50m8[,1,j]
  freq<- simdata50m8[,3,j]
  
  res50m8 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned50m8[j,1]<- res50m8$par[1]/res50m8$par[2] 
  MLEbinned50m8[j,2]<- (1/res50m8$par[2])^2 
  MLEbinned50m8[j,3]<- (1/res50m8$par[2]) 
  
}


##############################################################################################################
rm("TL","freq")

MLEbinned50m15<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned50m15)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata50m15[,1,j]
  freq<- simdata50m15[,3,j]
  
  res50m15 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned50m15[j,1]<- res50m15$par[1]/res50m15$par[2] 
  MLEbinned50m15[j,2]<- (1/res50m15$par[2])^2 
  MLEbinned50m15[j,3]<- (1/res50m15$par[2]) 
  
}

##############################################################################################################
rm("TL","freq")

MLEbinned50m30<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned50m30)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata50m30[,1,j]
  freq<- simdata50m30[,3,j]
  
  res50m30 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned50m30[j,1]<- res50m30$par[1]/res50m30$par[2] 
  MLEbinned50m30[j,2]<- (1/res50m30$par[2])^2 
  MLEbinned50m30[j,3]<- (1/res50m30$par[2]) 
  
}
#############################################################################################################
##############################################################################################################
rm("TL","freq")

MLEbinned100m8<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned100m8)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata100m8[,1,j]
  freq<- simdata100m8[,3,j]
  
  res100m8 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned100m8[j,1]<- res100m8$par[1]/res100m8$par[2] 
  MLEbinned100m8[j,2]<- (1/res100m8$par[2])^2 
  MLEbinned100m8[j,3]<- (1/res100m8$par[2]) 
  
}


##############################################################################################################
rm("TL","freq")

MLEbinned100m15<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned100m15)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata100m15[,1,j]
  freq<- simdata100m15[,3,j]
  
  res100m15 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned100m15[j,1]<- res100m15$par[1]/res100m15$par[2] 
  MLEbinned100m15[j,2]<- (1/res100m15$par[2])^2 
  MLEbinned100m15[j,3]<- (1/res100m15$par[2]) 
  
}

##############################################################################################################
rm("TL","freq")

MLEbinned100m30<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned100m30)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata100m30[,1,j]
  freq<- simdata100m30[,3,j]
  
  res100m30 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned100m30[j,1]<- res100m30$par[1]/res100m30$par[2] 
  MLEbinned100m30[j,2]<- (1/res100m30$par[2])^2 
  MLEbinned100m30[j,3]<- (1/res100m30$par[2]) 
  
}
#############################################################################################################
##############################################################################################################
rm("TL","freq")

MLEbinned300m8<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned300m8)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata300m8[,1,j]
  freq<- simdata300m8[,3,j]
  
  res300m8 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned300m8[j,1]<- res300m8$par[1]/res300m8$par[2] 
  MLEbinned300m8[j,2]<- (1/res300m8$par[2])^2 
  MLEbinned300m8[j,3]<- (1/res300m8$par[2]) 
  
}


##############################################################################################################
rm("TL","freq")

MLEbinned300m15<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned300m15)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata300m15[,1,j]
  freq<- simdata300m15[,3,j]
  
  res300m15 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned300m15[j,1]<- res300m15$par[1]/res300m15$par[2] 
  MLEbinned300m15[j,2]<- (1/res300m15$par[2])^2 
  MLEbinned300m15[j,3]<- (1/res300m15$par[2]) 
  
}

##############################################################################################################
rm("TL","freq")

MLEbinned300m30<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned300m30)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata300m30[,1,j]
  freq<- simdata300m30[,3,j]
  
  res300m30 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned300m30[j,1]<- res300m30$par[1]/res300m30$par[2] 
  MLEbinned300m30[j,2]<- (1/res300m30$par[2])^2 
  MLEbinned300m30[j,3]<- (1/res300m30$par[2]) 
  
}
#############################################################################################################
##############################################################################################################
rm("TL","freq")

MLEbinned600m8<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned600m8)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata600m8[,1,j]
  freq<- simdata600m8[,3,j]
  
  res600m8 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned600m8[j,1]<- res600m8$par[1]/res600m8$par[2] 
  MLEbinned600m8[j,2]<- (1/res600m8$par[2])^2 
  MLEbinned600m8[j,3]<- (1/res600m8$par[2]) 
  
}


##############################################################################################################
rm("TL","freq")

MLEbinned600m15<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned600m15)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata600m15[,1,j]
  freq<- simdata600m15[,3,j]
  
  res600m15 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned600m15[j,1]<- res600m15$par[1]/res600m15$par[2] 
  MLEbinned600m15[j,2]<- (1/res600m15$par[2])^2 
  MLEbinned600m15[j,3]<- (1/res600m15$par[2]) 
  
}

##############################################################################################################
rm("TL","freq")

MLEbinned600m30<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned600m30)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata600m30[,1,j]
  freq<- simdata600m30[,3,j]
  
  res600m30 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned600m30[j,1]<- res600m30$par[1]/res600m30$par[2] 
  MLEbinned600m30[j,2]<- (1/res600m30$par[2])^2 
  MLEbinned600m30[j,3]<- (1/res600m30$par[2]) 
  
}
#############################################################################################################
##############################################################################################################
rm("TL","freq")

MLEbinned1000m8<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned1000m8)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata1000m8[,1,j]
  freq<- simdata1000m8[,3,j]
  
  res1000m8 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned1000m8[j,1]<- res1000m8$par[1]/res1000m8$par[2] 
  MLEbinned1000m8[j,2]<- (1/res1000m8$par[2])^2 
  MLEbinned1000m8[j,3]<- (1/res1000m8$par[2]) 
  
}


##############################################################################################################
rm("TL","freq")

MLEbinned1000m15<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned1000m15)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata1000m15[,1,j]
  freq<- simdata1000m15[,3,j]
  
  res1000m15 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned1000m15[j,1]<- res1000m15$par[1]/res1000m15$par[2] 
  MLEbinned1000m15[j,2]<- (1/res1000m15$par[2])^2 
  MLEbinned1000m15[j,3]<- (1/res1000m15$par[2]) 
  
}

##############################################################################################################
rm("TL","freq")

MLEbinned1000m30<- matrix(rep(0,3*500),ncol=3)
colnames(MLEbinned1000m30)<- c("Mean","Var","Std")

for(j in 1:500){
  TL<- simdata1000m30[,1,j]
  freq<- simdata1000m30[,3,j]
  
  res1000m30 <- optim(c((67/2),(1/2)),fn=Logll,TL=TL,freq=freq,method="L-BFGS-B", lower=c(0.02,0.2),upper=c(180,2)) 
  MLEbinned1000m30[j,1]<- res1000m30$par[1]/res1000m30$par[2] 
  MLEbinned1000m30[j,2]<- (1/res1000m30$par[2])^2 
  MLEbinned1000m30[j,3]<- (1/res1000m30$par[2]) 
  
}
#############################################################################################################
##############################################################################################################
mlebinned50m8<- write.csv(MLEbinned50m8,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned50m8.csv")
mlebinned50m15<- write.csv(MLEbinned50m15,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned50m15.csv")
mlebinned50m30<- write.csv(MLEbinned50m30,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned50m30.csv")

mlebinned100m8<- write.csv(MLEbinned100m8,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned100m8.csv")
mlebinned100m15<- write.csv(MLEbinned100m15,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned100m15.csv")
mlebinned100m30<- write.csv(MLEbinned100m30,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned100m30.csv")

mlebinned300m8<- write.csv(MLEbinned300m8,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned300m8.csv")
mlebinned300m15<- write.csv(MLEbinned300m15,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned300m15.csv")
mlebinned300m30<- write.csv(MLEbinned300m30,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned300m30.csv")

mlebinned600m8<- write.csv(MLEbinned600m8,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned600m8.csv")
mlebinned600m15<- write.csv(MLEbinned600m15,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned600m15.csv")
mlebinned600m30<- write.csv(MLEbinned600m30,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned600m30.csv")

mlebinned1000m8<- write.csv(MLEbinned1000m8,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned1000m8.csv")
mlebinned1000m15<- write.csv(MLEbinned1000m15,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned1000m15.csv")
mlebinned1000m30<- write.csv(MLEbinned1000m30,"C:/Users/sh_za/Desktop/Results/MLEbinned/mlebinned1000m30.csv")


