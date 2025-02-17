library(tmvtnorm)
library(MASS)


Billik<- function(Data,theta){ 
  
  d<- (ncol(Data)-1)/2
  
  mu<- theta[c(1:d)] 
  
  sig<- matrix(rep(0,d*d),ncol=d) 
  
  subtheta<- theta[-c(1:d)] 
  
  diag(sig)<- subtheta[1:d] 
  subtheta1<- subtheta[-c(1:d)] 
  
  sig[lower.tri(sig)] <- sig[upper.tri(sig)]<- subtheta1
  print(sig)

  ind<- ncol(Data)-1 
  
  lowerb<- NULL
  upperb<- NULL 
  for (i in 1:ind){ 
    if (i%%2!=0) {lowerb<- append(lowerb,i)} 
    else {upperb<- append(upperb,i)} 
    
  }
  
  #print(lowerb)
  #print(upperb)
  Data1<- as.matrix(Data) 
  b<- 0 
  for(i in 1:nrow(Data1)){ 
    a<- pmvnorm(lower=Data1[i,c(lowerb)],upper=Data1[i,c(upperb)],mean=mu,sigma=sig) 
    #print(a)
    if (a==0) {a<- 1e-03} 
    #else {a<- a[1]}
    b<- b+Data[i,ncol(Data)]*log(a) 
    #print(b)
  }
  #print(-b)
  
  return(-b) 
}
#############################################
### Initialize theta for d=2 & d=3###
#mu1 <- c(67,67)
#sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
#2.16/(sqrt(3.1*6.05))

#mu2<- c(67,67,67)
#sigma2<- matrix(c(3.1,2.2,1.6,2.2,4.3,2.45,1.6,2.45,4.8),3,3)


#theta1<- c(mu1,as.vector(sigma1)



#theta2<- c(mu2,as.vector(sigma2))
#i=1
#dataMLE_n100_bin10[i,]<- nlm(Billik,Data=simulateddata100b10[,,i],theta<- c(67,67,3.10,6.05,2.15),hessian=TRUE)

#dataMLE_n100_bin10[i,]<- nlm(Billik,Data=simulateddata100b10[,,i],theta<- c(67,67,3.10,6.05,2.15),hessian=TRUE)$estimate
i=1
#nlm(Billik,Data=simulateddata100b10[,,i],theta<- c(67,67,3.10,6.05,2.15),hessian=TRUE)
#dataMLE_n100_bin10[i,]<- optim(theta<- c(67,67,3.5,6.5,2.5),Billik,Data=simulateddata100b10[,,i],hessian=TRUE)$estimate
#########################################################################################################################
Data<- simulateddata1000b10[,,1]
d=2
d=3
d<- (ncol(Data)-1)/2
print(d)
theta<- c(67,67,3.5,6.5,2.5)
theta<- c(67,67,3.5,6.5,0.5)
theta<- c(67,67,67,3,6,5,0.5,0.65,0.4)
mu<- theta[c(1:d)] 
mu
subtheta<- theta[-c(1:d)] 
subtheta


sig<- matrix(rep(0,d*d),ncol=d) 
sig
subtheta<- theta[-c(1:d)] 

diag(sig)<- subtheta[1:d] 
sig
subtheta1<- subtheta[-c(1:d)] 

sig[lower.tri(sig)] <- sig[upper.tri(sig)]<- subtheta1
print(sig)



###########################################################################################################################
dataMLE_n1000_bin10<- matrix(rep(0,5*20),ncol=5)
colnames(dataMLE_n1000_bin10)<- c("mux1","mux2","Vx1","Vx2","Sx1x2")

for(i in 1:nrow(dataMLE_n1000_bin10)){
  dataMLE_n1000_bin10[i,]<- nlm(Billik,Data=simulateddata1000b10[,,i],theta<- c(67,67,3.5,6.5,2.5),hessian=TRUE)$estimate
  #dataMLE_n1000_bin10[i,]<- optim(theta<- c(67,67,3.5,6.5,2.5),Billik,Data=simulateddata100b10[,,i],hessian=TRUE)$estimate
  

}
dataMLE_n1000_bin10
i=1
nlm(Billik,Data=simulateddata1000b10[,,i],theta<- c(67,67,3.5,6.5,2.5),hessian=TRUE)$estimate
i=6
nlm(Billik,Data=simulateddata1000b10[,,i],theta<- c(67,67,3.5,6.5,2.5),hessian=TRUE)$estimate


simulateddata1000b10[,,6]
##################################################################################################################
#################################################################################################################
##### MAIN CODE fOR BIVARIATE  ##############################################################################################################
library(tmvtnorm)
library(MASS)

Billik<- function(Data,theta){
  #mu<- c(theta[1],theta[2])
  Sigma<- matrix(c(theta[3],theta[5]*sqrt(theta[3])*sqrt(theta[4]),theta[5]*sqrt(theta[3])*sqrt(theta[4])
                   , theta[4]),byrow=T,ncol=2)
  #print(Sigma)
  #Sigma<- matrix(c(theta[3],theta[5],theta[5],theta[4]),byrow=T,ncol=2)
  b<- 0
  for(i in 1:nrow(Data)){
    a<- pmvnorm(lower=c(Data[i,1],Data[i,3]),upper=c(Data[i,2],Data[i,4]),mean=c(theta[1],theta[2]),sigma=Sigma)
    if (a[1]==0) {a<- 1e-03}
    else {a<- a[1]}
    b<- b+Data[i,5]*log(a)
    #print(b)
  }
  ##The joint likelihood
  #print(-b)
  return(-b)
  
}

theta<- c(65,65,3,6,2.12132)

Billik(DataBiv,theta)
################################
#MLEExact<- optim(Billik,Data=DataBiv,theta<- c(65,65,3,6,0.5),hessian=TRUE)

#MLEExact<- nlm(Billik,Data=DataBiv,theta<- c(67,67,3.2,6.2,0.5),hessian=TRUE)
#MLEExact$estimate

MLEExact<- nlm(Billik,Data=DataBiv,theta<- c(65,65,3,6,0.5),hessian=TRUE)
MLEExact$estimate
#system.time(nlm(Billik,Data=DataBiv,theta<- c(67,67,3.2,6.2,0.5),hessian=TRUE))
#print(paste("Parameter Estimates using nlm function",MLnlm$estimate))
#####################################################################################################################
dataMLE_n1000_bin10<- matrix(rep(0,5*30),ncol=5)
colnames(dataMLE_n1000_bin10)<- c("mux1","mux2","Vx1","Vx2","Sx1x2")

for(i in 1:nrow(dataMLE_n1000_bin10)){

  dataMLE_n1000_bin10[i,]<- nlm(Billik,Data=simulateddata1000b10[,,i],theta<- c(67,67,3.2,6.2,0.5),hessian=TRUE)$estimate

  
}

dataMLE_n1000_bin10
#######################################################################################################################
dataMLE_n1000_bin15<- matrix(rep(0,5*30),ncol=5)
colnames(dataMLE_n1000_bin15)<- c("mux1","mux2","Vx1","Vx2","Sx1x2")

for(i in 1:nrow(dataMLE_n1000_bin15)){
  #for(i in 1:2){
  
  dataMLE_n1000_bin15[i,]<- nlm(Billik,Data=simulateddata1000b15[,,i],theta<- c(67,67,3.2,6.2,0.5),hessian=TRUE)$estimate
  #dataMLE_n1000_bin10[i,]<- optim(theta<- c(67,67,3.5,6.5,2.5),Billik,Data=simulateddata100b10[,,i],hessian=TRUE)$estimate
  
  
}

dataMLE_n1000_bin15

####################################################################################################################
dataMLE_n1000_bin30<- matrix(rep(0,5*30),ncol=5)
colnames(dataMLE_n1000_bin30)<- c("mux1","mux2","Vx1","Vx2","Sx1x2")

for(i in 1:nrow(dataMLE_n1000_bin30)){
  #for(i in 1:2){
  
  dataMLE_n1000_bin30[i,]<- nlm(Billik,Data=simulateddata1000b30[,,i],theta<- c(67,67,3.2,6.2,0.5),hessian=TRUE)$estimate
  #dataMLE_n1000_bin10[i,]<- optim(theta<- c(67,67,3.5,6.5,2.5),Billik,Data=simulateddata100b10[,,i],hessian=TRUE)$estimate
  
  
}

dataMLE_n1000_bin30

OUTMLEBINNED1000b10<- write.csv(dataMLE_n100_bin10,"C:/Users/sh_za/Desktop/Results/Bivariate/MLE/dataMLE_n100_bin10.csv")


#######################################################################################

dataMLE_n100_bin15<- matrix(rep(0,5*30),ncol=5)
colnames(dataMLE_n100_bin15)<- c("mux1","mux2","Vx1","Vx2","Sx1x2")

for(i in 1:nrow(dataMLE_n100_bin15)){
  dataMLE_n100_bin15[i,]<- nlm(Billik,Data=simulateddata100b15[,,i],theta<- c(67,67,3.5,6.5,2.5),hessian=TRUE)$estimate
  #dataMLE_n100_bin10[i,]<- optim(theta<- c(67,67,3.5,6.5,2.5),Billik,Data=simulateddata100b10[,,i],hessian=TRUE)$estimate
  
  
}
dataMLE_n100_bin15


OUTMLEBINNED1000b10<- write.csv(dataMLE_n100_bin10,"C:/Users/sh_za/Desktop/Results/Bivariate/MLE/dataMLE_n100_bin10.csv")


#######################################################################################


