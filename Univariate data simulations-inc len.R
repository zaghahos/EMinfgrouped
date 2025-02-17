############UNIVARIATE Data SIMULATIONS for n=50,100,300,600,1000 and bins=8,15,30######################################
rm(list=ls())
set.seed(4321)
#####Univariate Simulations#####
### n=50 & 8 Classes###

sim50m8<- matrix(rep(0,50*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim50m8)){
  sim50m8[,i]<- rnorm(50,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,8*ncol(sim50m8)),ncol=ncol(sim50m8))

for(i in 1:ncol(sim50m8)){
  Fr[,i]<- table(cut(sim50m8[,i],breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf)))
  
}
#print(Fr)

simdata50m8<- array(rep(0,8*3*500),c(8,3,500))
med1<- array(rep(0,8*2*500),c(8,2,500))
for(i in 1:ncol(sim50m8)){
  simdata50m8[,1,i]<- c(-Inf,64,65.5,67,68.5,70,71.5,73)
  simdata50m8[,2,i]<- c(64,65.5,67,68.5,70,71.5,73,Inf)
  simdata50m8[,3,i]<- Fr[,i]
  med1[,1,i]<- (simdata50m8[,1,i]+simdata50m8[,2,i])/2
  med1[1,1,i]<- 63.25
  med1[nrow(simdata50m8),1,i]<- 73.75
  med1[,2,i]<- Fr[,i]
}

Mynewarray50m8<- array(rep(0*8*2*500),c(8,2,500))
for (i in 1:500){
  Mynewarray50m8[,1,i]<- (simdata50m8[,1,i]+simdata50m8[,2,i])/2
  Mynewarray50m8[,2,i]<- simdata50m8[,3,i]
  Mynewarray50m8[1,1,i]<- 63.25
  Mynewarray50m8[nrow(simdata50m8),1,i]<- 73.75
}
#Mynewarray

newdata50m8<- matrix(rep(0*500*50),ncol = 500,nrow=50)
for (i in 1:500){
  newdata50m8[,i]<- rep(Mynewarray50m8[,1,i],Mynewarray50m8[,2,i])
}
#newdata2[100:150,]
#head(newdata50m8)
#newdata50m8[1:50,1:5]
Me50m8<- colMeans(newdata50m8)
SS50m8<- rep(0,500)
SDD50m8<- rep(0,500)
for(i in 1:500){
  SS50m8[i]<- (sum((newdata50m8[,i]-Me50m8[i])^2))/nrow(newdata50m8)
  SDD50m8[i]<- sqrt(SS50m8[i])
}
#SS2
#SDD2
MLEunbinned50m8<- cbind(Me50m8,SS50m8,SDD50m8)
colnames(MLEunbinned50m8)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

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

Mynewarray50m15<- array(rep(0*15*2*500),c(15,2,500))
for (i in 1:500){
  Mynewarray50m15[,1,i]<- (simdata50m15[,1,i]+simdata50m15[,2,i])/2
  Mynewarray50m15[,2,i]<- simdata50m15[,3,i]
  Mynewarray50m15[1,1,i]<- 58.75
  Mynewarray50m15[nrow(simdata50m15),1,i]<- 79.75
}
#Mynewarray

newdata50m15<- matrix(rep(0*500*50),ncol = 500,nrow=50)
for (i in 1:500){
  newdata50m15[,i]<- rep(Mynewarray50m15[,1,i],Mynewarray50m15[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me50m15<- colMeans(newdata50m15)
SS50m15<- rep(0,500)
SDD50m15<- rep(0,500)
for(i in 1:500){
  SS50m15[i]<- (sum((newdata50m15[,i]-Me50m15[i])^2))/nrow(newdata50m15)
  SDD50m15[i]<- sqrt(SS50m15[i])
}
#SS2
#SDD2
MLEunbinned50m15<- cbind(Me50m15,SS50m15,SDD50m15)
colnames(MLEunbinned50m15)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
### n=50 & 30 Classes###

sim50m30<- matrix(rep(0,50*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim50m30)){
  sim50m30[,i]<- rnorm(50,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,30*ncol(sim50m30)),ncol=ncol(sim50m30))

for(i in 1:ncol(sim50m30)){
  Fr[,i]<- table(cut(sim50m30[,i],breaks=c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                                           65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                                           83.5,85,86.5,88,89.5,Inf)))
  
}
#print(Fr)

simdata50m30<- array(rep(0,30*3*500),c(30,3,500))
med3<- array(rep(0,30*2*500),c(30,2,500))
for(i in 1:ncol(sim50m30)){
  simdata50m30[,1,i]<- c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                         65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                         83.5,85,86.5,88,89.5)
  
  simdata50m30[,2,i]<- c(47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                         65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                         83.5,85,86.5,88,89.5,Inf)
  
  simdata50m30[,3,i]<- Fr[,i]
  med3[,1,i]<- (simdata50m30[,1,i]+simdata50m30[,2,i])/2
  med3[1,1,i]<- 46.75
  med3[nrow(simdata50m30),1,i]<- 90.25
  med3[,2,i]<- Fr[,i]
}

Mynewarray50m30<- array(rep(0*30*2*500),c(30,2,500))
for (i in 1:500){
  Mynewarray50m30[,1,i]<- (simdata50m30[,1,i]+simdata50m30[,2,i])/2
  Mynewarray50m30[,2,i]<- simdata50m30[,3,i]
  Mynewarray50m30[1,1,i]<- 46.75
  Mynewarray50m30[nrow(simdata50m30),1,i]<- 90.25
}
#Mynewarray

newdata50m30<- matrix(rep(0*500*50),ncol = 500,nrow=50)
for (i in 1:500){
  newdata50m30[,i]<- rep(Mynewarray50m30[,1,i],Mynewarray50m30[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me50m30<- colMeans(newdata50m30)
SS50m30<- rep(0,500)
SDD50m30<- rep(0,500)
for(i in 1:500){
  SS50m30[i]<- (sum((newdata50m30[,i]-Me50m30[i])^2))/nrow(newdata50m30)
  SDD50m30[i]<- sqrt(SS50m30[i])
}
#SS2
#SDD2
MLEunbinned50m30<- cbind(Me50m30,SS50m30,SDD50m30)
colnames(MLEunbinned50m30)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
############################################################################
############################################################################
### n=100 & 8 Classes###

sim100m8<- matrix(rep(0,100*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim100m8)){
  sim100m8[,i]<- rnorm(100,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,8*ncol(sim100m8)),ncol=ncol(sim100m8))

for(i in 1:ncol(sim100m8)){
  Fr[,i]<- table(cut(sim100m8[,i],breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf)))
  
}
#print(Fr)

simdata100m8<- array(rep(0,8*3*500),c(8,3,500))
med4<- array(rep(0,8*2*500),c(8,2,500))
for(i in 1:ncol(sim100m8)){
  simdata100m8[,1,i]<- c(-Inf,64,65.5,67,68.5,70,71.5,73)
  simdata100m8[,2,i]<- c(64,65.5,67,68.5,70,71.5,73,Inf)
  simdata100m8[,3,i]<- Fr[,i]
  med4[,1,i]<- (simdata100m8[,1,i]+simdata100m8[,2,i])/2
  med4[1,1,i]<- 63.25
  med4[nrow(simdata100m8),1,i]<- 73.75
  med4[,2,i]<- Fr[,i]
}

Mynewarray100m8<- array(rep(0*8*2*500),c(8,2,500))
for (i in 1:500){
  Mynewarray100m8[,1,i]<- (simdata100m8[,1,i]+simdata100m8[,2,i])/2
  Mynewarray100m8[,2,i]<- simdata100m8[,3,i]
  Mynewarray100m8[1,1,i]<- 63.25
  Mynewarray100m8[nrow(simdata100m8),1,i]<- 73.75
}
#Mynewarray

newdata100m8<- matrix(rep(0*500*100),ncol = 500,nrow=100)
for (i in 1:500){
  newdata100m8[,i]<- rep(Mynewarray100m8[,1,i],Mynewarray100m8[,2,i])
}
#newdata2[100:150,]
#head(newdata50m8)
#newdata50m8[1:50,1:5]
Me100m8<- colMeans(newdata100m8)
SS100m8<- rep(0,500)
SDD100m8<- rep(0,500)
for(i in 1:500){
  SS100m8[i]<- (sum((newdata100m8[,i]-Me100m8[i])^2))/nrow(newdata100m8)
  SDD100m8[i]<- sqrt(SS100m8[i])
}
#SS2
#SDD2
MLEunbinned100m8<- cbind(Me100m8,SS100m8,SDD100m8)
colnames(MLEunbinned100m8)<- c("Mean","Var","Sd")

MLEunbinned50m8<- cbind(Me50m8,SS50m8,SDD50m8)
colnames(MLEunbinned50m8)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

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

Mynewarray100m15<- array(rep(0*15*2*500),c(15,2,500))
for (i in 1:500){
  Mynewarray100m15[,1,i]<- (simdata100m15[,1,i]+simdata100m15[,2,i])/2
  Mynewarray100m15[,2,i]<- simdata100m15[,3,i]
  Mynewarray100m15[1,1,i]<- 58.75
  Mynewarray100m15[nrow(simdata100m15),1,i]<- 79.75
}
#Mynewarray

newdata100m15<- matrix(rep(0*500*100),ncol = 500,nrow=100)
for (i in 1:500){
  newdata100m15[,i]<- rep(Mynewarray100m15[,1,i],Mynewarray100m15[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me100m15<- colMeans(newdata100m15)
SS100m15<- rep(0,500)
SDD100m15<- rep(0,500)
for(i in 1:500){
  SS100m15[i]<- (sum((newdata100m15[,i]-Me100m15[i])^2))/nrow(newdata100m15)
  SDD100m15[i]<- sqrt(SS100m15[i])
}
#SS2
#SDD2
MLEunbinned100m15<- cbind(Me100m15,SS100m15,SDD100m15)
colnames(MLEunbinned100m15)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
### n=100 & 30 Classes###

sim100m30<- matrix(rep(0,100*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim100m30)){
  sim100m30[,i]<- rnorm(100,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,30*ncol(sim100m30)),ncol=ncol(sim100m30))

for(i in 1:ncol(sim100m30)){
  Fr[,i]<- table(cut(sim100m30[,i],breaks=c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                                            65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                                            83.5,85,86.5,88,89.5,Inf)))
  
}
#print(Fr)

simdata100m30<- array(rep(0,30*3*500),c(30,3,500))
med6<- array(rep(0,30*2*500),c(30,2,500))
for(i in 1:ncol(sim100m30)){
  simdata100m30[,1,i]<- c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                          65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                          83.5,85,86.5,88,89.5)
  
  simdata100m30[,2,i]<- c(47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                          65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                          83.5,85,86.5,88,89.5,Inf)
  
  simdata100m30[,3,i]<- Fr[,i]
  med6[,1,i]<- (simdata100m30[,1,i]+simdata100m30[,2,i])/2
  med6[1,1,i]<- 46.75
  med6[nrow(simdata100m30),1,i]<- 90.25
  med6[,2,i]<- Fr[,i]
}

Mynewarray100m30<- array(rep(0*30*2*500),c(30,2,500))
for (i in 1:500){
  Mynewarray100m30[,1,i]<- (simdata100m30[,1,i]+simdata100m30[,2,i])/2
  Mynewarray100m30[,2,i]<- simdata100m30[,3,i]
  Mynewarray100m30[1,1,i]<- 46.75
  Mynewarray100m30[nrow(simdata100m30),1,i]<- 90.25
}
#Mynewarray

newdata100m30<- matrix(rep(0*500*100),ncol = 500,nrow=100)
for (i in 1:500){
  newdata100m30[,i]<- rep(Mynewarray100m30[,1,i],Mynewarray100m30[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me100m30<- colMeans(newdata100m30)
SS100m30<- rep(0,500)
SDD100m30<- rep(0,500)
for(i in 1:500){
  SS100m30[i]<- (sum((newdata100m30[,i]-Me100m30[i])^2))/nrow(newdata100m30)
  SDD100m30[i]<- sqrt(SS100m30[i])
}
#SS2
#SDD2
MLEunbinned100m30<- cbind(Me100m30,SS100m30,SDD100m30)
colnames(MLEunbinned100m30)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
############################################################################
############################################################################
### n=300 & 8 Classes###

sim300m8<- matrix(rep(0,300*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim300m8)){
  sim300m8[,i]<- rnorm(300,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,8*ncol(sim300m8)),ncol=ncol(sim300m8))

for(i in 1:ncol(sim300m8)){
  Fr[,i]<- table(cut(sim300m8[,i],breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf)))
  
}
#print(Fr)

simdata300m8<- array(rep(0,8*3*500),c(8,3,500))
med7<- array(rep(0,8*2*500),c(8,2,500))
for(i in 1:ncol(sim300m8)){
  simdata300m8[,1,i]<- c(-Inf,64,65.5,67,68.5,70,71.5,73)
  simdata300m8[,2,i]<- c(64,65.5,67,68.5,70,71.5,73,Inf)
  simdata300m8[,3,i]<- Fr[,i]
  med7[,1,i]<- (simdata300m8[,1,i]+simdata300m8[,2,i])/2
  med7[1,1,i]<- 63.25
  med7[nrow(simdata300m8),1,i]<- 73.75
  med7[,2,i]<- Fr[,i]
}

Mynewarray300m8<- array(rep(0*8*2*500),c(8,2,500))
for (i in 1:500){
  Mynewarray300m8[,1,i]<- (simdata300m8[,1,i]+simdata300m8[,2,i])/2
  Mynewarray300m8[,2,i]<- simdata300m8[,3,i]
  Mynewarray300m8[1,1,i]<- 63.25
  Mynewarray300m8[nrow(simdata300m8),1,i]<- 73.75
}
#Mynewarray

newdata300m8<- matrix(rep(0*500*300),ncol = 500,nrow=300)
for (i in 1:500){
  newdata300m8[,i]<- rep(Mynewarray300m8[,1,i],Mynewarray300m8[,2,i])
}
#newdata2[100:150,]
#head(newdata50m8)
#newdata50m8[1:50,1:5]
Me300m8<- colMeans(newdata300m8)
SS300m8<- rep(0,500)
SDD300m8<- rep(0,500)
for(i in 1:500){
  SS300m8[i]<- (sum((newdata300m8[,i]-Me300m8[i])^2))/nrow(newdata300m8)
  SDD300m8[i]<- sqrt(SS300m8[i])
}
#SS2
#SDD2
MLEunbinned300m8<- cbind(Me300m8,SS300m8,SDD300m8)
colnames(MLEunbinned300m8)<- c("Mean","Var","Sd")

#SDD2<- stdmle2

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

Mynewarray300m15<- array(rep(0*15*2*500),c(15,2,500))
for (i in 1:500){
  Mynewarray300m15[,1,i]<- (simdata300m15[,1,i]+simdata300m15[,2,i])/2
  Mynewarray300m15[,2,i]<- simdata300m15[,3,i]
  Mynewarray300m15[1,1,i]<- 58.75
  Mynewarray300m15[nrow(simdata300m15),1,i]<- 79.75
}
#Mynewarray

newdata300m15<- matrix(rep(0*500*300),ncol = 500,nrow=300)
for (i in 1:500){
  newdata300m15[,i]<- rep(Mynewarray300m15[,1,i],Mynewarray300m15[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me300m15<- colMeans(newdata300m15)
SS300m15<- rep(0,500)
SDD300m15<- rep(0,500)
for(i in 1:500){
  SS300m15[i]<- (sum((newdata300m15[,i]-Me300m15[i])^2))/nrow(newdata300m15)
  SDD300m15[i]<- sqrt(SS300m15[i])
}
#SS2
#SDD2
MLEunbinned300m15<- cbind(Me300m15,SS300m15,SDD300m15)
colnames(MLEunbinned300m15)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
### n=300 & 30 Classes###

sim300m30<- matrix(rep(0,300*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim300m30)){
  sim300m30[,i]<- rnorm(300,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,30*ncol(sim300m30)),ncol=ncol(sim300m30))

for(i in 1:ncol(sim300m30)){
  Fr[,i]<- table(cut(sim300m30[,i],breaks=c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                                            65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                                            83.5,85,86.5,88,89.5,Inf)))
  
}
#print(Fr)

simdata300m30<- array(rep(0,30*3*500),c(30,3,500))
med9<- array(rep(0,30*2*500),c(30,2,500))
for(i in 1:ncol(sim300m30)){
  simdata300m30[,1,i]<- c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                          65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                          83.5,85,86.5,88,89.5)
  
  simdata300m30[,2,i]<- c(47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                          65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                          83.5,85,86.5,88,89.5,Inf)
  
  simdata300m30[,3,i]<- Fr[,i]
  med9[,1,i]<- (simdata300m30[,1,i]+simdata300m30[,2,i])/2
  med9[1,1,i]<- 46.75
  med9[nrow(simdata300m30),1,i]<- 90.25
  med9[,2,i]<- Fr[,i]
}

Mynewarray300m30<- array(rep(0*30*2*500),c(30,2,500))
for (i in 1:500){
  Mynewarray300m30[,1,i]<- (simdata300m30[,1,i]+simdata300m30[,2,i])/2
  Mynewarray300m30[,2,i]<- simdata300m30[,3,i]
  Mynewarray300m30[1,1,i]<- 46.75
  Mynewarray300m30[nrow(simdata300m30),1,i]<- 90.25
}
#Mynewarray

newdata300m30<- matrix(rep(0*500*300),ncol = 500,nrow=300)
for (i in 1:500){
  newdata300m30[,i]<- rep(Mynewarray300m30[,1,i],Mynewarray300m30[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me300m30<- colMeans(newdata300m30)
SS300m30<- rep(0,500)
SDD300m30<- rep(0,500)
for(i in 1:500){
  SS300m30[i]<- (sum((newdata300m30[,i]-Me300m30[i])^2))/nrow(newdata300m30)
  SDD300m30[i]<- sqrt(SS300m30[i])
}
#SS2
#SDD2
MLEunbinned300m30<- cbind(Me300m30,SS300m30,SDD300m30)
colnames(MLEunbinned300m30)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
############################################################################
############################################################################
### n=600 & 8 Classes###

sim600m8<- matrix(rep(0,600*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim600m8)){
  sim600m8[,i]<- rnorm(600,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,8*ncol(sim600m8)),ncol=ncol(sim600m8))

for(i in 1:ncol(sim600m8)){
  Fr[,i]<- table(cut(sim600m8[,i],breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf)))
  
}
#print(Fr)

simdata600m8<- array(rep(0,8*3*500),c(8,3,500))
med10<- array(rep(0,8*2*500),c(8,2,500))
for(i in 1:ncol(sim600m8)){
  simdata600m8[,1,i]<- c(-Inf,64,65.5,67,68.5,70,71.5,73)
  simdata600m8[,2,i]<- c(64,65.5,67,68.5,70,71.5,73,Inf)
  simdata600m8[,3,i]<- Fr[,i]
  med10[,1,i]<- (simdata600m8[,1,i]+simdata600m8[,2,i])/2
  med10[1,1,i]<- 63.25
  med10[nrow(simdata600m8),1,i]<- 73.75
  med10[,2,i]<- Fr[,i]
}

Mynewarray600m8<- array(rep(0*8*2*500),c(8,2,500))
for (i in 1:500){
  Mynewarray600m8[,1,i]<- (simdata600m8[,1,i]+simdata600m8[,2,i])/2
  Mynewarray600m8[,2,i]<- simdata600m8[,3,i]
  Mynewarray600m8[1,1,i]<- 63.25
  Mynewarray600m8[nrow(simdata600m8),1,i]<- 73.75
}
#Mynewarray

newdata600m8<- matrix(rep(0*500*600),ncol = 500,nrow=600)
for (i in 1:500){
  newdata600m8[,i]<- rep(Mynewarray600m8[,1,i],Mynewarray600m8[,2,i])
}
#newdata2[100:150,]
#head(newdata50m8)
#newdata50m8[1:50,1:5]
Me600m8<- colMeans(newdata600m8)
SS600m8<- rep(0,500)
SDD600m8<- rep(0,500)
for(i in 1:500){
  SS600m8[i]<- (sum((newdata600m8[,i]-Me600m8[i])^2))/nrow(newdata600m8)
  SDD600m8[i]<- sqrt(SS600m8[i])
}
#SS2
#SDD2

MLEunbinned600m8<- cbind(Me600m8,SS600m8,SDD600m8)
colnames(MLEunbinned600m8)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

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

Mynewarray600m15<- array(rep(0*15*2*500),c(15,2,500))
for (i in 1:500){
  Mynewarray600m15[,1,i]<- (simdata600m15[,1,i]+simdata600m15[,2,i])/2
  Mynewarray600m15[,2,i]<- simdata600m15[,3,i]
  Mynewarray600m15[1,1,i]<- 58.75
  Mynewarray600m15[nrow(simdata600m15),1,i]<- 79.75
}
#Mynewarray

newdata600m15<- matrix(rep(0*500*600),ncol = 500,nrow=600)
for (i in 1:500){
  newdata600m15[,i]<- rep(Mynewarray600m15[,1,i],Mynewarray600m15[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me600m15<- colMeans(newdata600m15)
SS600m15<- rep(0,500)
SDD600m15<- rep(0,500)
for(i in 1:500){
  SS600m15[i]<- (sum((newdata600m15[,i]-Me600m15[i])^2))/nrow(newdata600m15)
  SDD600m15[i]<- sqrt(SS600m15[i])
}
#SS2
#SDD2
MLEunbinned600m15<- cbind(Me600m15,SS600m8,SDD600m15)
colnames(MLEunbinned600m15)<- c("Mean","Var","Sd")


#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
### n=600 & 30 Classes###

sim600m30<- matrix(rep(0,600*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim600m30)){
  sim600m30[,i]<- rnorm(600,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,30*ncol(sim600m30)),ncol=ncol(sim600m30))

for(i in 1:ncol(sim600m30)){
  Fr[,i]<- table(cut(sim600m30[,i],breaks=c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                                            65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                                            83.5,85,86.5,88,89.5,Inf)))
  
}
#print(Fr)

simdata600m30<- array(rep(0,30*3*500),c(30,3,500))
med12<- array(rep(0,30*2*500),c(30,2,500))
for(i in 1:ncol(sim600m30)){
  simdata600m30[,1,i]<- c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                          65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                          83.5,85,86.5,88,89.5)
  
  simdata600m30[,2,i]<- c(47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                          65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                          83.5,85,86.5,88,89.5,Inf)
  
  simdata600m30[,3,i]<- Fr[,i]
  med12[,1,i]<- (simdata600m30[,1,i]+simdata600m30[,2,i])/2
  med12[1,1,i]<- 46.75
  med12[nrow(simdata600m30),1,i]<- 90.25
  med12[,2,i]<- Fr[,i]
}

Mynewarray600m30<- array(rep(0*30*2*500),c(30,2,500))
for (i in 1:500){
  Mynewarray600m30[,1,i]<- (simdata600m30[,1,i]+simdata600m30[,2,i])/2
  Mynewarray600m30[,2,i]<- simdata600m30[,3,i]
  Mynewarray600m30[1,1,i]<- 46.75
  Mynewarray600m30[nrow(simdata600m30),1,i]<- 90.25
}
#Mynewarray

newdata600m30<- matrix(rep(0*500*600),ncol = 500,nrow=600)
for (i in 1:500){
  newdata600m30[,i]<- rep(Mynewarray600m30[,1,i],Mynewarray600m30[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me600m30<- colMeans(newdata600m30)
SS600m30<- rep(0,500)
SDD600m30<- rep(0,500)
for(i in 1:500){
  SS600m30[i]<- (sum((newdata600m30[,i]-Me600m30[i])^2))/nrow(newdata600m30)
  SDD600m30[i]<- sqrt(SS600m30[i])
}
#SS2
#SDD2
MLEunbinned600m30<- cbind(Me600m30,SS600m30,SDD600m30)
colnames(MLEunbinned600m30)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
############################################################################
############################################################################
### n=1000 & 8 Classes###

sim1000m8<- matrix(rep(0,1000*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim1000m8)){
  sim1000m8[,i]<- rnorm(1000,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,8*ncol(sim1000m8)),ncol=ncol(sim1000m8))

for(i in 1:ncol(sim1000m8)){
  Fr[,i]<- table(cut(sim1000m8[,i],breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf)))
  
}
#print(Fr)

simdata1000m8<- array(rep(0,8*3*500),c(8,3,500))
med13<- array(rep(0,8*2*500),c(8,2,500))
for(i in 1:ncol(sim1000m8)){
  simdata1000m8[,1,i]<- c(-Inf,64,65.5,67,68.5,70,71.5,73)
  simdata1000m8[,2,i]<- c(64,65.5,67,68.5,70,71.5,73,Inf)
  simdata1000m8[,3,i]<- Fr[,i]
  med13[,1,i]<- (simdata1000m8[,1,i]+simdata1000m8[,2,i])/2
  med13[1,1,i]<- 63.25
  med13[nrow(simdata1000m8),1,i]<- 73.75
  med13[,2,i]<- Fr[,i]
}

Mynewarray1000m8<- array(rep(0*8*2*500),c(8,2,500))
for (i in 1:500){
  Mynewarray1000m8[,1,i]<- (simdata1000m8[,1,i]+simdata1000m8[,2,i])/2
  Mynewarray1000m8[,2,i]<- simdata1000m8[,3,i]
  Mynewarray1000m8[1,1,i]<- 63.25
  Mynewarray1000m8[nrow(simdata1000m8),1,i]<- 73.75
}
#Mynewarray

newdata1000m8<- matrix(rep(0*500*1000),ncol = 500,nrow=1000)
for (i in 1:500){
  newdata1000m8[,i]<- rep(Mynewarray1000m8[,1,i],Mynewarray1000m8[,2,i])
}
#newdata2[100:150,]
#head(newdata50m8)
#newdata50m8[1:50,1:5]
Me1000m8<- colMeans(newdata1000m8)
SS1000m8<- rep(0,500)
SDD1000m8<- rep(0,500)
for(i in 1:500){
  SS1000m8[i]<- (sum((newdata1000m8[,i]-Me1000m8[i])^2))/nrow(newdata1000m8)
  SDD1000m8[i]<- sqrt(SS1000m8[i])
}
#SS2
#SDD2
MLEunbinned1000m8<- cbind(Me1000m8,SS1000m8,SDD1000m8)
colnames(MLEunbinned1000m8)<- c("Mean","Var","Sd")

#Me2<- Memle2
#SDD2<- stdmle2

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

Mynewarray1000m15<- array(rep(0*15*2*500),c(15,2,500))
for (i in 1:500){
  Mynewarray1000m15[,1,i]<- (simdata1000m15[,1,i]+simdata1000m15[,2,i])/2
  Mynewarray1000m15[,2,i]<- simdata1000m15[,3,i]
  Mynewarray1000m15[1,1,i]<- 58.75
  Mynewarray1000m15[nrow(simdata1000m15),1,i]<- 79.75
}
#Mynewarray

newdata1000m15<- matrix(rep(0*500*1000),ncol = 500,nrow=1000)
for (i in 1:500){
  newdata1000m15[,i]<- rep(Mynewarray1000m15[,1,i],Mynewarray1000m15[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me1000m15<- colMeans(newdata1000m15)
SS1000m15<- rep(0,500)
SDD1000m15<- rep(0,500)
for(i in 1:500){
  SS1000m15[i]<- (sum((newdata1000m15[,i]-Me1000m15[i])^2))/nrow(newdata1000m15)
  SDD1000m15[i]<- sqrt(SS1000m15[i])
}
#SS2
#SDD2
MLEunbinned1000m15<- cbind(Me1000m15,SS1000m15,SDD1000m15)
colnames(MLEunbinned1000m15)<- c("Mean","Var","Sd")


#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
### n=1000 & 30 Classes###

sim1000m30<- matrix(rep(0,1000*500),ncol=500)
#ncol(sim1)

for(i in 1:ncol(sim1000m30)){
  sim1000m30[,i]<- rnorm(1000,68,2.5)
}


#dim(sim2)
Fr<- matrix(rep(0,30*ncol(sim1000m30)),ncol=ncol(sim1000m30))

for(i in 1:ncol(sim1000m30)){
  Fr[,i]<- table(cut(sim1000m30[,i],breaks=c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                                             65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                                             83.5,85,86.5,88,89.5,Inf)))
                                             
                                             
                                        
  
}
#print(Fr)

simdata1000m30<- array(rep(0,30*3*500),c(30,3,500))
med15<- array(rep(0,30*2*500),c(30,2,500))
for(i in 1:ncol(sim1000m30)){
  simdata1000m30[,1,i]<- c(-Inf,47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                           65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                           83.5,85,86.5,88,89.5)
  
  simdata1000m30[,2,i]<- c(47.5,49,50.5,52,53.5,55,56.5,58,59.5,61,62.5,64,
                           65.5,67,68.5,70,71.5,73,74.5,76,77.5,79,80.5,82,
                           83.5,85,86.5,88,89.5,Inf)
  
  simdata1000m30[,3,i]<- Fr[,i]
  med15[,1,i]<- (simdata1000m30[,1,i]+simdata1000m30[,2,i])/2
  med15[1,1,i]<- 46.75
  med15[nrow(simdata1000m30),1,i]<- 90.25
  med15[,2,i]<- Fr[,i]
}

Mynewarray1000m30<- array(rep(0*30*2*500),c(30,2,500))
for (i in 1:500){
  Mynewarray1000m30[,1,i]<- (simdata1000m30[,1,i]+simdata1000m30[,2,i])/2
  Mynewarray1000m30[,2,i]<- simdata1000m30[,3,i]
  Mynewarray1000m30[1,1,i]<- 46.75
  Mynewarray1000m30[nrow(simdata1000m30),1,i]<- 90.25
}
#Mynewarray

newdata1000m30<- matrix(rep(0*500*1000),ncol = 500,nrow=1000)
for (i in 1:500){
  newdata1000m30[,i]<- rep(Mynewarray1000m30[,1,i],Mynewarray1000m30[,2,i])
}
#newdata2[100:150,]
#head(newdata50m15)
#newdata50m15[1:50,1:5]
Me1000m30<- colMeans(newdata1000m30)
SS1000m30<- rep(0,500)
SDD1000m30<- rep(0,500)
for(i in 1:500){
  SS1000m30[i]<- (sum((newdata1000m30[,i]-Me1000m30[i])^2))/nrow(newdata1000m30)
  SDD1000m30[i]<- sqrt(SS1000m30[i])
}
#SS2
#SDD2
MLEunbinned1000m30<- cbind(Me1000m30,SS1000m30,SDD1000m30)
colnames(MLEunbinned1000m30)<- c("Mean","Var","Sd")
##################################################################################################
##################################################################################################
#### N=50 ##########################################################
BR50m8<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR50m8)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR50m8[,1]<- Me50m8-68
BR50m8[,2]<- SS50m8-6.25
BR50m8[,3]<- (Me50m8-68)^2
BR50m8[,4]<- (SS50m8-6.25)^2

BR50m15<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR50m15)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR50m15[,1]<- Me50m15-68
BR50m15[,2]<- SS50m15-6.25
BR50m15[,3]<- (Me50m15-68)^2
BR50m15[,4]<- (SS50m15-6.25)^2

BR50m30<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR50m30)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR50m30[,1]<- Me50m30-68
BR50m30[,2]<- SS50m30-6.25
BR50m30[,3]<- (Me50m30-68)^2
BR50m30[,4]<- (SS50m30-6.25)^2

#### N=100 ##########################################################
BR100m8<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR100m8)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR100m8[,1]<- Me100m8-68
BR100m8[,2]<- SS100m8-6.25
BR100m8[,3]<- (Me100m8-68)^2
BR100m8[,4]<- (SS100m8-6.25)^2

BR100m15<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR100m15)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR100m15[,1]<- Me100m15-68
BR100m15[,2]<- SS100m15-6.25
BR100m15[,3]<- (Me100m15-68)^2
BR100m15[,4]<- (SS100m15-6.25)^2

BR100m30<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR100m30)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR100m30[,1]<- Me100m30-68
BR100m30[,2]<- SS100m30-6.25
BR100m30[,3]<- (Me100m30-68)^2
BR100m30[,4]<- (SS100m30-6.25)^2
####N= 300#########################################################
BR300m8<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR300m8)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR300m8[,1]<- Me300m8-68
BR300m8[,2]<- SS300m8-6.25
BR300m8[,3]<- (Me300m8-68)^2
BR300m8[,4]<- (SS300m8-6.25)^2

BR300m15<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR300m15)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR300m15[,1]<- Me300m15-68
BR300m15[,2]<- SS300m15-6.25
BR300m15[,3]<- (Me300m15-68)^2
BR300m15[,4]<- (SS300m15-6.25)^2

BR300m30<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR300m30)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR300m30[,1]<- Me300m30-68
BR300m30[,2]<- SS300m30-6.25
BR300m30[,3]<- (Me300m30-68)^2
BR300m30[,4]<- (SS300m30-6.25)^2
####N= 600#########################################################
BR600m8<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR600m8)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR600m8[,1]<- Me600m8-68
BR600m8[,2]<- SS600m8-6.25
BR600m8[,3]<- (Me600m8-68)^2
BR600m8[,4]<- (SS600m8-6.25)^2

BR600m15<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR600m15)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR600m15[,1]<- Me600m15-68
BR600m15[,2]<- SS600m15-6.25
BR600m15[,3]<- (Me600m15-68)^2
BR600m15[,4]<- (SS600m15-6.25)^2

BR600m30<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR600m30)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR600m30[,1]<- Me600m30-68
BR600m30[,2]<- SS600m30-6.25
BR600m30[,3]<- (Me600m30-68)^2
BR600m30[,4]<- (SS600m30-6.25)^2

####N= 1000#########################################################
BR1000m8<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR1000m8)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR1000m8[,1]<- Me1000m8-68
BR1000m8[,2]<- SS1000m8-6.25
BR1000m8[,3]<- (Me1000m8-68)^2
BR1000m8[,4]<- (SS1000m8-6.25)^2

BR1000m15<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR1000m15)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR1000m15[,1]<- Me1000m15-68
BR1000m15[,2]<- SS1000m15-6.25
BR1000m15[,3]<- (Me1000m15-68)^2
BR1000m15[,4]<- (SS1000m15-6.25)^2

BR1000m30<- matrix(rep(0*500*4),ncol=4,nrow = 500) 
colnames(BR1000m30)<- c("Mean Bias","var Bias","Mean MSE","Var MSE")
BR1000m30[,1]<- Me1000m30-68
BR1000m30[,2]<- SS1000m30-6.25
BR1000m30[,3]<- (Me1000m30-68)^2
BR1000m30[,4]<- (SS1000m30-6.25)^2
#Me2<- Memle2
#SDD2<- stdmle2

############################################################################
############################################################################
############################################################################
simulation50m8<- write.csv(simdata50m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation50m8.csv")
simulation50m15<- write.csv(simdata50m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation50m15.csv")
simulation50m15<- write.csv(simdata50m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation50m30.csv")

simulation100m8<- write.csv(simdata100m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation100m8.csv")
simulation100m15<- write.csv(simdata100m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation100m15.csv")
simulation100m30<- write.csv(simdata100m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation100m30.csv")

simulation300m8<- write.csv(simdata300m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation300m8.csv")
simulation300m15<- write.csv(simdata300m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation300m15.csv")
simulation300m30<- write.csv(simdata300m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation300m30.csv")

simulation600m8<- write.csv(simdata600m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation600m8.csv")
simulation600m15<- write.csv(simdata600m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation600m15.csv")
simulation600m30<- write.csv(simdata600m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation600m30.csv")

simulation1000m8<- write.csv(simdata1000m8,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation1000m8.csv")
simulation1000m15<- write.csv(simdata1000m15,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation1000m15.csv")
simulation1000m30<- write.csv(simdata1000m30,"C:/Users/sh_za/OneDrive/Desktop/Results/Univariate/Simulation/simulation1000m30.csv")

