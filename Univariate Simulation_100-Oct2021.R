############SIMULATIONS#####################################################
rm(list=ls())
set.seed(2222)
#####Univariate Simulations#####

############################################################################
### n=50 & 15 Classes###

sim50m15<- matrix(rep(0,50*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim50m15)){
  sim50m15[,i]<- rnorm(50,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,15*ncol(sim50m15)),ncol=ncol(sim50m15))

for(i in 1:ncol(sim50m15)){
  Fr[,i]<- table(cut(sim50m15[,i],breaks=c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)))
  
  
}
#print(Fr)

simdata50m15<- array(rep(0,15*3*500),c(15,3,500))
med2<- array(rep(0,15*2*500),c(15,2,500))
for(i in 1:ncol(sim50m15)){
  simdata50m15[,1,i]<- c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79)
  simdata50m15[,2,i]<- c(59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)
  simdata50m15[,3,i]<- Fr[,i]
  med2[,1,i]<- (simdata50m15[,1,i]+simdata50m15[,2,i])/2
  med2[1,1,i]<- 58.75
  med2[nrow(simdata50m15),1,i]<- 79.75
  med2[,2,i]<- Fr[,i]
}

#simdata50m15

############################################################################
############################################################################

############################################################################
### n=100 & 15 Classes###

sim100m15<- matrix(rep(0,100*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim100m15)){
  sim100m15[,i]<- rnorm(100,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,15*ncol(sim100m15)),ncol=ncol(sim100m15))

for(i in 1:ncol(sim100m15)){
  Fr[,i]<- table(cut(sim100m15[,i],breaks=c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)))
  
}
#print(Fr)

simdata100m15<- array(rep(0,15*3*500),c(15,3,500))
med5<- array(rep(0,15*2*500),c(15,2,500))
for(i in 1:ncol(sim100m15)){
  simdata100m15[,1,i]<- c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79)
  simdata100m15[,2,i]<- c(59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)
  simdata100m15[,3,i]<- Fr[,i]
  med5[,1,i]<- (simdata100m15[,1,i]+simdata100m15[,2,i])/2
  med5[1,1,i]<- 58.75
  med5[nrow(simdata100m15),1,i]<- 79.75
  med5[,2,i]<- Fr[,i]
}

############################################################################

############################################################################
### n=300 & 15 Classes###

sim300m15<- matrix(rep(0,300*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim300m15)){
  sim300m15[,i]<- rnorm(300,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,15*ncol(sim300m15)),ncol=ncol(sim300m15))

for(i in 1:ncol(sim300m15)){
  Fr[,i]<- table(cut(sim300m15[,i],breaks=c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)))
  
}
#print(Fr)

simdata300m15<- array(rep(0,15*3*500),c(15,3,500))
med8<- array(rep(0,15*2*500),c(15,2,500))
for(i in 1:ncol(sim300m15)){
  simdata300m15[,1,i]<- c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79)
  simdata300m15[,2,i]<- c(59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)
  simdata300m15[,3,i]<- Fr[,i]
  med8[,1,i]<- (simdata300m15[,1,i]+simdata300m15[,2,i])/2
  med8[1,1,i]<- 58.75
  med8[nrow(simdata300m15),1,i]<- 79.75
  med8[,2,i]<- Fr[,i]
}

############################################################################
############################################################################
### n=600 & 15 Classes###

sim600m15<- matrix(rep(0,600*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim600m15)){
  sim600m15[,i]<- rnorm(600,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,15*ncol(sim600m15)),ncol=ncol(sim600m15))

for(i in 1:ncol(sim600m15)){
  Fr[,i]<- table(cut(sim600m15[,i],breaks=c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)))
  
}
#print(Fr)

simdata600m15<- array(rep(0,15*3*500),c(15,3,500))
med11<- array(rep(0,15*2*500),c(15,2,500))
for(i in 1:ncol(sim600m15)){
  simdata600m15[,1,i]<- c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79)
  simdata600m15[,2,i]<- c(59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)
  simdata600m15[,3,i]<- Fr[,i]
  med11[,1,i]<- (simdata600m15[,1,i]+simdata600m15[,2,i])/2
  med11[1,1,i]<- 58.75
  med11[nrow(simdata600m15),1,i]<- 79.75
  med11[,2,i]<- Fr[,i]
}

############################################################################
############################################################################
############################################################################
### n=1000 & 15 Classes###

sim1000m15<- matrix(rep(0,1000*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim1000m15)){
  sim1000m15[,i]<- rnorm(1000,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,15*ncol(sim1000m15)),ncol=ncol(sim1000m15))

for(i in 1:ncol(sim1000m15)){
  Fr[,i]<- table(cut(sim1000m15[,i],breaks=c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)))
  
}
#print(Fr)

simdata1000m15<- array(rep(0,15*3*500),c(15,3,500))
med14<- array(rep(0,15*2*500),c(15,2,500))
for(i in 1:ncol(sim1000m15)){
  simdata1000m15[,1,i]<- c(-Inf,59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79)
  simdata1000m15[,2,i]<- c(59.5,61,62.5,64,65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,Inf)
  simdata1000m15[,3,i]<- Fr[,i]
  med14[,1,i]<- (simdata1000m15[,1,i]+simdata1000m15[,2,i])/2
  med14[1,1,i]<- 58.75
  med14[nrow(simdata1000m15),1,i]<- 79.75
  med14[,2,i]<- Fr[,i]
}
##################################################################################################
##################################################################################################

############################################################################
############################################################################
############################################################################
simulation50m15<- write.csv(simdata50m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/simulation50m15.csv")

simulation100m15<- write.csv(simdata100m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/simulation100m15.csv")

simulation300m15<- write.csv(simdata300m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/simulation300m15.csv")

simulation600m15<- write.csv(simdata600m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/simulation600m15.csv")

simulation1000m15<- write.csv(simdata1000m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulations1/simulation1000m15.csv")

