
###############################################################################################
###############################################################################################
###############################################################################################
###data4: n=1000##########################
### Trivariate grouped data simulation: to simulate 30 data sets
mm<- c(68,68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,1.5,2,4,2.5,1.5,2.5,5),3,3) ## The covariance matrix for the simulation

ssdata<- array(rep(0,3000*3*30),c(3000,3,30)) ## create an initial array of zero's
for(i in 1:30){ ## use the for loop to create 30 data sets
  ssdata[,,i]<- mvrnorm(n=3000,mm,ss) ## simulate data from Multivariate normal distribtion for 3 variables 
}
#ssdata


x<- matrix(rep(0,3000*30),ncol=30) ## initialize matrix of X 
y<- matrix(rep(0,3000*30),ncol=30) ## Initial matrix of Y
z<- matrix(rep(0,3000*30),ncol=30) ## initialize matrix of Z
for(i in 1:30){## use the for loop to do the steps over 30 data sets
  x[,i]<- ssdata[,1,i] ## assign the first column of each simulated array to the variable X
  y[,i]<- ssdata[,2,i] ## assign the second column of each simulated array to the variable Y
  z[,i]<- ssdata[,3,i] ## assign the third column of each simulated array to the variable Z
}


Freqtable1<- array(rep(0,8*8*8*30),c(8,8,8,30)) ## create an array of the initial frequencies which is set to be zero
for(i in 1:30){ ## use the for loop to do the steps for 30 data sets
  x.cut<- cut(x[,i],breaks=c(-Inf,65,66,67,68,69,70,71,Inf)) ## define the break points of the variable X for building 8 intervals 
  y.cut<- cut(y[,i],breaks=c(-Inf,63,65,67,69,71,73,75,Inf)) ## define the break points of the variable Y for building 8 intervals
  z.cut<- cut(z[,i],breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf)) ## define the break points of the variable Z for building 8 intervals
  
  Freqtable1[,,,i]<- table(x.cut,y.cut,z.cut) ## Set the cut points to the frequency array, to count how many data points are 
  ### falling into the surfaces
  
}

#Freqtable1
#Freqtable1[,,,5]
simulateddata<- array(rep(0,30*7*512),c(512,7,30)) ## To Create a trivariate grouped data, first create an initialized array
### with the values of zero's
### This array is included 30 matrix (for each data sets) of 512 (rows) and 7 (columns)- 512 is the number of total intervals
### 8 (intervals of X)* 8 (intervals of Y)* 8 (intervals of Z)- Note that the number of intervals could be different for 
### the variables, but for simplicity we consider the equal number of intervals 
### 7 is the number of columns of each dimension of this array: the column are as follows:
### Columns 1,3,5 (odd numbers): are the lower bounds of the variables X, Y , Z
### Columns 2,4,6 (even numbers): are the upper bounds of the variables X, Y , Z
### Columns 7: are the Frequency of the variables X, Y , Z
lower.x<- rep(c(-Inf,65,66,67,68,69,70,71),64) ## the lower bounds of X
ly<- c(rep(-Inf,8),rep(63,8),rep(65,8),rep(67,8),rep(69,8),rep(71,8),rep(73,8),rep(75,8)) ## to create lower bounds of Y
lower.y<- rep(ly,8) ## lower bounds of Y
lower.z<- c(rep(-Inf,64),rep(64,64),rep(65.5,64),rep(67,64),rep(68.5,64),rep(70,64),rep(71.5,64),rep(73,64)) ## Lower bounds of Z

upper.x<- rep(c(65,66,67,68,69,70,71,Inf),64) ## Upper bounds of X
uy<- c(rep(63,8),rep(65,8),rep(67,8),rep(69,8),rep(71,8),rep(73,8),rep(75,8),rep(Inf,8)) ## To create upper bounds of Y
upper.y<- rep(uy,8) ## Upper bounds of Y
upper.z<- c(rep(64,64),rep(65.5,64),rep(67,64),rep(68.5,64),rep(70,64),
            rep(71.5,64),rep(73,64),rep(Inf,64)) ## upper bounds of Z


for(i in 1:30){ ## use the for loop to build the simulated grouped data over all 30 data sets 
  simulateddata[,1,i]<- lower.x ## Set the first column to the lower bound of X
  simulateddata[,3,i]<- lower.y ## Set the third column to the lower bound of Y
  simulateddata[,5,i]<- lower.z ## Set the fifth column to the lower bound of Z 
  simulateddata[,2,i]<- upper.x ## Set the second column to the upper bound of X
  simulateddata[,4,i]<- upper.y ## Set the forth column to the upper bound of Y
  simulateddata[,6,i]<- upper.z ## Set the sixth column to the upper bound of Z
  simulateddata[,7,i]<- c(Freqtable1[,,,i]) ## Set the seventh column to the frequencies
}

#simulateddata

###################################################################################
### To calculate the MLE's of the parameteres ignore the grouping :
### Here we use the midpoints of the intervals as observations
medx<- (simulateddata[,1,1]+simulateddata[,2,1])/2  ## to calcualte the midpoints of the variable X

for (i in 1:length(medx)){ ## to replace the -Inf & Inf with the proper number for X
  if (medx[i]==-Inf) medx[i]<- 64.5 ## replace -Inf with 64.5
  else if (medx[i]==Inf) medx[i]<- 71.5 ## replace Inf with 71.5
}


medy<- (simulateddata[,3,1]+simulateddata[,4,1])/2 ## to calcualte the midpoints of the variable Y

for (i in 1:length(medy)){ ## to replace the -Inf & Inf with the proper number for Y
  if (medy[i]==-Inf) medy[i]<- 62 ## replace -Inf with 62
  else if (medy[i]==Inf) medy[i]<- 76 ## replace Inf with 76
}



medz<- (simulateddata[,5,1]+simulateddata[,6,1])/2 ## to calcualte the midpoints of the variable Z

for (i in 1:length(medz)){ ## to replace the -Inf & Inf with the proper number for X
  if (medz[i]==-Inf) medz[i]<- 63.25 ## replace -Inf with 63.25
  else if (medz[i]==Inf) medz[i]<- 73.75 ## replace Inf with 73.75
}


midpt<- array(rep(0, 512*4*30),c(512,4,30)) ## create an initial array of zero's for the midpoints of all intervals of 30 datasets
for (i in 1:30){ ## use for loop for doing the process over 30 data sets
 midpt[,,i]<- cbind(medx,medy,medz,simulateddata[,7,i]) ## create each third dimension of array of the columnbind matrix of 
 ### the midpoints of the variables and the frequencies
 colnames(midpt)<- c("medx","medy","medz","Freq")## Rename the columns of the matrix
 }


MLEest<- matrix(rep(0,9*30),ncol=9,nrow=30) ## To calculate the MLE of the parameters, initialize a matrix of zero's with 
### 9 columns, and 30 rows: one columns for each parameter and one row for each data set
colnames(MLEest)<- c("meanx","meany","meanz","sxx","syy","szz","sxy","sxz","syz") ## name the columns of the matrix
for(i in 1:30){## use for loop to do the calculations over 30 data sets

  MLEest[i,1]<- sum(midpt[,1,i]*midpt[,4,i])/(sum(midpt[,4,i])) ## calculate the mean of X using the  formula on page (11)
  MLEest[i,4]<- sum(((midpt[,1,i]-MLEest[i,1])^2)*midpt[,4,i])/(sum(midpt[,4,i])) ## calculate the variance of X using the  formula on page (11)
  
  MLEest[i,2]<- sum(medy*simulateddata[,7,i])/(sum(simulateddata[,7,i])) ## calculate the mean of Y using the  formula on page (11)
  MLEest[i,5]<- sum(((medy-MLEest[i,2])^2)*simulateddata[,7,i])/(sum(simulateddata[,7,i])) ## calculate the variance of Y using the  formula on page (11)
  
  MLEest[i,3]<- sum(medz*simulateddata[,7,i])/(sum(simulateddata[,7,i])) ## calculate the mean of Z using the  formula on page (11)
  MLEest[i,6]<- sum(((medz-MLEest[i,3])^2)*simulateddata[,7,i])/(sum(simulateddata[,7,i])) ## calculate the variance of Z using the  formula on page (11)
  
  
  
  MLEest[i,7]<- sum(((medx-MLEest[i,1])*(medy-MLEest[i,2]))*simulateddata[,7,i])/(sum(simulateddata[,7,i]))
  ### calculate the covariances of X & Y using the  formula on page (11)
  MLEest[i,8]<- sum(((medx-MLEest[i,1])*(medz-MLEest[i,3]))*simulateddata[,7,i])/(sum(simulateddata[,7,i]))
  ### calculate the covariances of X & Z using the  formula on page (11)
  
  MLEest[i,9]<- sum(((medz-MLEest[i,3])*(medy-MLEest[i,2]))*simulateddata[,7,i])/(sum(simulateddata[,7,i]))
  ### calculate the covariances of Y & Z using the  formula on page (11)
  
  
}

MLEest ## the MLE estimation matrix 
colMeans(MLEest) ## calculate the column means to see if the average results iver the 30 data sets are close to the default values
################################################################################
### Here for each parameter, we have calculated its mean, sd, and mean squarred error (MSE)  over the 30 simulated data sets
### and put the results of each patrameter on a vector, that means we have 9 parameters, so we will build 9 vectors of the 
### values of mean, sd and MSE of each parameters over the 30 data sets
### Then we produce a matrix of rbind (row bind) of these 9 vectors, which is a 9*3 matrix- 9 rows for each parameter, 
### 3 columns: for the calculated results of mean, sd and MSE 



MMX<- mean(MLEest[,1])
MMX
MSDX<- sd(MLEest[,1])
MMSEX<- (sum((MLEest[,1]-68)^2))/nrow(MLEest)
MX<- c(MMX,MSDX,MMSEX)
###################################################
MMY<- mean(MLEest[,2])
MSDY<- sd(MLEest[,2])
MMSEY<- (sum((MLEest[,2]-68)^2))/nrow(MLEest)
MY<- c(MMY,MSDY,MMSEY)
###################################################
MMz<- mean(MLEest[,3])
MSDz<- sd(MLEest[,3])
MMSEz<- (sum((MLEest[,3]-68)^2))/nrow(MLEest)
MZ<- c(MMz,MSDz,MMSEz)

###################################################
MVX<- mean(MLEest[,4])
VSDX<- sd(MLEest[,4])
VMSEX<- (sum((MLEest[,4]-3)^2))/nrow(MLEest)
VX<- c(MVX,VSDX,VMSEX)

###################################################
MVY<- mean(MLEest[,5])
VSDY<- sd(MLEest[,5])
VMSEY<- (sum((MLEest[,5]-6)^2))/nrow(MLEest)
VY<- c(MVY,VSDY,VMSEY)

###################################################
MVz<- mean(MLEest[,6])
VSDz<- sd(MLEest[,6])
VMSEz<- (sum((MLEest[,6]-5)^2))/nrow(MLEest)
VZ<- c(MVz,VSDz,VMSEz)

###################################################
MVxY<- mean(MLEest[,7])
VSDxY<- sd(MLEest[,7])
VMSExY<- (sum((MLEest[,7]-2)^2))/nrow(MLEest)
VXY<- c(MVxY,VSDxY,VMSExY)
###################################################
MVxz<- mean(MLEest[,8])
VSDxz<- sd(MLEest[,8])
VMSExz<- (sum((MLEest[,8]-1.5)^2))/nrow(MLEest)
VXZ<- c(MVxz,VSDxz,VMSExz)
###################################################
MVYz<- mean(MLEest[,9])
VSDYz<- sd(MLEest[,9])
VMSEYz<- (sum((MLEest[,9]-2.5)^2))/nrow(MLEest)
VYZ<- c(MVYz,VSDYz,VMSEYz)
###################################################

#Parameters<- c("Mean X","Mean Y","Var X","Var Y","Rho")
#criteria<- c("mean","std","MSE")
MLETRI<- rbind(MX,MY,MZ,VX,VY,VZ,VXY,VXZ,VYZ)
colnames(MLETRI)<- c("mean","std","MSE")
rownames(MLETRI)<- c("Mean X","Mean Y","Mean Z","Var X","Var Y","Var Z","S XY","S XZ","S YZ")
print(MLETRI)


#simulateddata
#library(xtable)
#xtable(MLETRI,digits=7,caption = "Estimate of MLE Parameters for Bivariate Simulation")


### THE END###_ Do not see the rests
########################################################################################################################
############################################################################################
############################################################################################
#simulateddata
#medx
#medy

newsim<- array(rep(0,2*30*1000),c(1000,2,30))
for(i in 1:30){
  newsim[,1,i]<- rep(medx,simulateddata[,5,i])
  newsim[,2,i]<- rep(medy,simulateddata[,5,i])
}
#newsim

#newsim[,,1]

MLEB1<- matrix(rep(0,5*30),ncol=5,nrow=30)
colnames(MLEB1)<- c("MX","MY","VX","VY","rXY")
for(i in 1:30){
  MLEB1[i,1]<- mean(newsim[,1,i])
  MLEB1[i,2]<- mean(newsim[,2,i])
  MLEB1[i,3]<- sum((newsim[,1,i]-mean(newsim[,1,i]))^2)/nrow(newsim[,,i])
  MLEB1[i,4]<- sum((newsim[,2,i]-mean(newsim[,2,i]))^2)/nrow(newsim[,,i])
  MLEB1[i,5]<- ((sum((newsim[,1,i]-mean(newsim[,1,i]))*(newsim[,2,i]-mean(newsim[,2,i]))))/nrow(newsim[,,i]))*(1/sqrt(MLEB1[i,3]*MLEB1[i,4]))
  
}
#MLEB1