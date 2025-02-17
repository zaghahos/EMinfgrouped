library(tmvtnorm)
library(MASS)
set.seed(54321)
###data: n=50##########################
### Bivariate grouped data simulation: to simulate 10 data sets
mm<- c(68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,2,6),2,2) ## The covariance matrix for the simulation

ssdata<- array(rep(0,50*2*500),c(50,2,500)) 
for(i in 1: 500){
  ssdata[,,i]<- mvrnorm(n=50,mm,ss) 
}
#ssdata


x<- matrix(rep(0,50*500),ncol=500)
y<- matrix(rep(0,50*500),ncol=500) 

for(i in 1:500){ 
  x[,i]<- ssdata[,1,i] 
  y[,i]<- ssdata[,2,i] 
}


Freqtable1<- array(rep(0,10*10*500),c(10,10,500)) 

for(i in 1:500){  
  x.cut<- cut(x[,i],breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))  
  y.cut<- cut(y[,i],breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))  
  Freqtable1[,,i]<- table(x.cut,y.cut) 
  
}


simulateddata50b10<- array(rep(0,10*5*10*500),c(10*10,5,500)) 

lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10)

lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10)) 

upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10) 

upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10)) 

for(i in 1:500){ ## use the for loop to build the simulated grouped data over all 30 data set
  simulateddata50b10[,1,i]<- lower.x 
  simulateddata50b10[,3,i]<- lower.y 
  simulateddata50b10[,2,i]<- upper.x 
  simulateddata50b10[,4,i]<- upper.y 
  simulateddata50b10[,5,i]<- c(Freqtable1[,,i]) 
}

simulateddata50b10

#########################################################################################################
simdata50b10<- write.csv(simulateddata50b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/Simulations1/simulateddata50b10.csv")
##############################################################################################################

###data: n=100##########################
### Bivariate grouped data simulation: to simulate 10 data sets
mm<- c(68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,2,6),2,2) ## The covariance matrix for the simulation

ssdata<- array(rep(0,100*2*500),c(100,2,500)) 
for(i in 1: 500){
  ssdata[,,i]<- mvrnorm(n=100,mm,ss) 
}
#ssdata


x<- matrix(rep(0,100*500),ncol=500)
y<- matrix(rep(0,100*500),ncol=500) 

for(i in 1:500){ 
  x[,i]<- ssdata[,1,i] 
  y[,i]<- ssdata[,2,i] 
}


Freqtable1<- array(rep(0,10*10*500),c(10,10,500)) 

for(i in 1:500){  
  x.cut<- cut(x[,i],breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))  
  y.cut<- cut(y[,i],breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))  
  Freqtable1[,,i]<- table(x.cut,y.cut) 
  
}


simulateddata100b10<- array(rep(0,10*5*10*500),c(10*10,5,500)) 

lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10)

lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10)) 

upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10) 

upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10)) 

for(i in 1:500){ ## use the for loop to build the simulated grouped data over all 30 data set
  simulateddata100b10[,1,i]<- lower.x 
  simulateddata100b10[,3,i]<- lower.y 
  simulateddata100b10[,2,i]<- upper.x 
  simulateddata100b10[,4,i]<- upper.y 
  simulateddata100b10[,5,i]<- c(Freqtable1[,,i]) 
}

simulateddata100b10 
#########################################
simdata100b10<- write.csv(simulateddata100b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/Simulations1/simulateddata100b10.csv")
############################################################################################################################

######################################################################################################################################
###data: n=600##########################
### Bivariate grouped data simulation: to simulate 10 data sets
mm<- c(68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,2,6),2,2) ## The covariance matrix for the simulation

ssdata<- array(rep(0,600*2*500),c(600,2,500)) 
for(i in 1: 500){
  ssdata[,,i]<- mvrnorm(n=600,mm,ss) 
}
#ssdata


x<- matrix(rep(0,600*500),ncol=500)
y<- matrix(rep(0,600*500),ncol=500) 

for(i in 1:500){ 
  x[,i]<- ssdata[,1,i] 
  y[,i]<- ssdata[,2,i] 
}


Freqtable1<- array(rep(0,10*10*500),c(10,10,500)) 

for(i in 1:500){  
  x.cut<- cut(x[,i],breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))  
  y.cut<- cut(y[,i],breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))  
  Freqtable1[,,i]<- table(x.cut,y.cut) 
  
}


simulateddata600b10<- array(rep(0,10*5*10*500),c(10*10,5,500)) 

lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10)

lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10)) 

upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10) 

upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10)) 

for(i in 1:500){ ## use the for loop to build the simulated grouped data over all 30 data set
  simulateddata600b10[,1,i]<- lower.x 
  simulateddata600b10[,3,i]<- lower.y 
  simulateddata600b10[,2,i]<- upper.x 
  simulateddata600b10[,4,i]<- upper.y 
  simulateddata600b10[,5,i]<- c(Freqtable1[,,i]) 
}

simulateddata600b10 

simdata600b10<- write.csv(simulateddata600b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/Simulations1/simulateddata600b10.csv")

########################################################################################################################################
###data: n=1000##########################
### Bivariate grouped data simulation: to simulate 10 data sets
mm<- c(68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,2,6),2,2) ## The covariance matrix for the simulation

ssdata<- array(rep(0,1000*2*500),c(1000,2,500)) 
for(i in 1: 500){
  ssdata[,,i]<- mvrnorm(n=1000,mm,ss) 
}
#ssdata


x<- matrix(rep(0,1000*500),ncol=500)
y<- matrix(rep(0,1000*500),ncol=500) 

for(i in 1:500){ 
  x[,i]<- ssdata[,1,i] 
  y[,i]<- ssdata[,2,i] 
}


Freqtable1<- array(rep(0,10*10*500),c(10,10,500)) 

for(i in 1:500){  
  x.cut<- cut(x[,i],breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))  
  y.cut<- cut(y[,i],breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))  
  Freqtable1[,,i]<- table(x.cut,y.cut) 
  
}


simulateddata1000b10<- array(rep(0,10*5*10*500),c(10*10,5,500)) 

lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10)

lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10)) 

upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10) 

upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10)) 

for(i in 1:500){ ## use the for loop to build the simulated grouped data over all 30 data set
  simulateddata1000b10[,1,i]<- lower.x 
  simulateddata1000b10[,3,i]<- lower.y 
  simulateddata1000b10[,2,i]<- upper.x 
  simulateddata1000b10[,4,i]<- upper.y 
  simulateddata1000b10[,5,i]<- c(Freqtable1[,,i]) 
}

simulateddata1000b10 

simdata1000b10<- write.csv(simulateddata1000b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/Simulations1/simulateddata1000b10.csv")

##################################################################################################################################
##################################################################################################################################
#set.seed(6789)

###data: n=1500##########################
### Bivariate grouped data simulation: to simulate 10 data sets
mm<- c(68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,2,6),2,2) ## The covariance matrix for the simulation

ssdata<- array(rep(0,1500*2*500),c(1500,2,500)) 
for(i in 1: 500){
  ssdata[,,i]<- mvrnorm(n=1500,mm,ss) 
}
#ssdata


x<- matrix(rep(0,1500*500),ncol=500)
y<- matrix(rep(0,1500*500),ncol=500) 

for(i in 1:500){ 
  x[,i]<- ssdata[,1,i] 
  y[,i]<- ssdata[,2,i] 
}


Freqtable1<- array(rep(0,10*10*500),c(10,10,500)) 

for(i in 1:500){  
  x.cut<- cut(x[,i],breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))  
  y.cut<- cut(y[,i],breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))  
  Freqtable1[,,i]<- table(x.cut,y.cut) 
  
}


simulateddata1500b10<- array(rep(0,10*5*10*500),c(10*10,5,500)) 

lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10)

lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10)) 

upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10) 

upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10)) 

for(i in 1:500){ ## use the for loop to build the simulated grouped data over all 30 data set
  simulateddata1500b10[,1,i]<- lower.x 
  simulateddata1500b10[,3,i]<- lower.y 
  simulateddata1500b10[,2,i]<- upper.x 
  simulateddata1500b10[,4,i]<- upper.y 
  simulateddata1500b10[,5,i]<- c(Freqtable1[,,i]) 
}

simulateddata1500b10 

simdata1500b10<- write.csv(simulateddata1500b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/Simulations1/simulateddata1500b10.csv")







###data: n=300##########################
### Bivariate grouped data simulation: to simulate 10 data sets
mm<- c(68,68) 
ss<- matrix(c(3,2,2,6),2,2) 

ssdata<- array(rep(0,300*2*500),c(300,2,500)) 
for(i in 1: 500){ 
  ssdata[,,i]<- mvrnorm(n=300,mm,ss) 
}



x<- matrix(rep(0,300*500),ncol=500)
y<- matrix(rep(0,300*500),ncol=500) 

for(i in 1:500){ 
  x[,i]<- ssdata[,1,i] 
  y[,i]<- ssdata[,2,i] 
}


Freqtable1<- array(rep(0,10*10*500),c(10,10,500)) 

for(i in 1:500){ 
  x.cut<- cut(x[,i],breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))  
  
  y.cut<- cut(y[,i],breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))   
  Freqtable1[,,i]<- table(x.cut,y.cut)  
  
}


simulateddata300b10<- array(rep(0,10*5*10*500),c(10*10,5,500)) 


lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10) 

lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10)) 


upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10) 

upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10)) 
for(i in 1:500){ 
  simulateddata300b10[,1,i]<- lower.x 
  simulateddata300b10[,3,i]<- lower.y
  simulateddata300b10[,2,i]<- upper.x 
  simulateddata300b10[,4,i]<- upper.y 
  simulateddata300b10[,5,i]<- c(Freqtable1[,,i]) 
}


simulateddata300b10

simdata300b10<- write.csv(simulateddata300b10,"C:/Users/sh_za/OneDrive/Desktop/Results/Bivariate/Simulations1/simulateddata300b10.csv")
