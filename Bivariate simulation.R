library(tmvtnorm)
library(MASS)

###data: n=1000##########################
### Bivariate grouped data simulation: to simulate 30 data sets
mm<- c(68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,2,6),2,2) ## The covariance matrix for the simulation

ssdata<- array(rep(0,1000*2*30),c(1000,2,30)) ## create an initial array of zero's
for(i in 1: 30){ ## use the for loop to create 30 data sets
  ssdata[,,i]<- mvrnorm(n=1000,mm,ss) ## simulate data from Multivariate normal distribtion for 2 variables
}
#ssdata


x<- matrix(rep(0,1000*30),ncol=30) ## initialize matrix of X
y<- matrix(rep(0,1000*30),ncol=30) ## Initial matrix of Y

for(i in 1:30){ ## use the for loop to do the steps over 30 data sets
  x[,i]<- ssdata[,1,i] ## assign the first column of each simulated array to the variable X
  y[,i]<- ssdata[,2,i] ## assign the second column of each simulated array to the variable Y
}


Freqtable1<- array(rep(0,10*10*30),c(10,10,30)) ## create an array of the initial frequencies which is set to be zero
for(i in 1:30){ ## use the for loop to do the steps for 30 data sets
  x.cut<- cut(x[,i],breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))  ## define the break points of the variable X for building 10 intervals 
  y.cut<- cut(y[,i],breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))  ## define the break points of the variable Y for building 10 intervals 
  Freqtable1[,,i]<- table(x.cut,y.cut) ## Set the cut points to the frequency array, to count how many data points are 
  ### falling into the rectangles
  
}


simulateddata<- array(rep(0,30*5*100),c(100,5,30)) ## To Create a Bivariate grouped data, first create an initialized array
### with the values of zero's
### This array is included 30 matrix (for each data sets) of 100 (rows) and 5 (columns)- 100 is the number of total intervals
### 10 (intervals of X)* 10 (intervals of Y)- Note that the number of intervals could be different for 
### the variables, but for simplicity we consider the equal number of intervals 
### 5 is the number of columns of each dimension of this array: the column are as follows:
### Columns 1,3 (odd numbers): are the lower bounds of the variables X, Y 
### Columns 2,4 (even numbers): are the upper bounds of the variables X, Y 
### Columns 5: are the Frequency of the variables X, Y 
lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10) ## the lower bounds of X
lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10)) ## the lower bounds of Y
upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10) ## Upper bounds of X
upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10)) ## Upper bounds of Y
for(i in 1:30){ ## use the for loop to build the simulated grouped data over all 30 data sets
  simulateddata[,1,i]<- lower.x ## Set the first column to the lower bound of X
  simulateddata[,3,i]<- lower.y ## Set the third column to the lower bound of Y
  simulateddata[,2,i]<- upper.x ## Set the second column to the upper bound of X
  simulateddata[,4,i]<- upper.y ## Set the forth column to the upper bound of Y
  simulateddata[,5,i]<- c(Freqtable1[,,i]) ## Set the fifth column to the frequencies
}

#simulateddata
##############################################
### To calculate the MLE's of the parameteres ignore the grouping :
### Here we use the midpoints of the intervals as observations
medx<- (simulateddata[,1,1]+simulateddata[,2,1])/2 ## to calcualte the midpoints of the variable X



for (i in 1:length(medx)){ ## to replace the -Inf & Inf with the proper number for X
  if (medx[i]==-Inf) medx[i]<- 63.5 ## replace -Inf with 63.5
  else if (medx[i]==Inf) medx[i]<- 72.5 ## replace Inf with 72.5
}


medy<- (simulateddata[,3,1]+simulateddata[,4,1])/2 ## to calcualte the midpoints of the variable Y


for (i in 1:length(medy)){ ## to replace the -Inf & Inf with the proper number for Y
  if (medy[i]==-Inf) medy[i]<- 63.7 ## replace -Inf with 63.4
  else if (medy[i]==Inf) medy[i]<- 72.7 ## replace Inf with 72.7
}


midpt<- array(rep(0, 100*3*30),c(100,3,30)) ## create an initial array of zero's for the midpoints of all intervals of 30 datasets
for (i in 1:30){ ## use for loop for doing the process over 30 data sets
  midpt[,,i]<- cbind(medx,medy,simulateddata[,5,i]) ## create each third dimension of array of the columnbind matrix of 
  ### the midpoints of the variables and the frequencies
  colnames(midpt)<- c("medx","medy","Freq")## Rename the columns of the matrix
}



MLEest<- matrix(rep(0,6*30),ncol=6,nrow=30) ## To calculate the MLE of the parameters, initialize a matrix of zero's with 
### 6 columns, and 30 rows: one columns for each parameter and one row for each data set
colnames(MLEest)<- c("meanx","meany","sxx","syy","sxy","rho") ## name the columns of the matrix
for(i in 1:30){## use for loop to do the calculations over 30 data sets
  
  MLEest[i,1]<- sum(midpt[,1,i]*midpt[,3,i])/(sum(midpt[,3,i])) ## calculate the mean of X using the  formula on page (10)
  MLEest[i,3]<- sum(((midpt[,1,i]-MLEest[i,1])^2)*midpt[,3,i])/(sum(midpt[,3,i])) ## calculate the variance of X using the  formula on page (10)
  
  MLEest[i,2]<- sum(midpt[,2,i]*midpt[,3,i])/(sum(midpt[,3,i])) ## calculate the mean of Y using the  formula on page (10)
  MLEest[i,4]<- sum(((midpt[,2,i]-MLEest[i,2])^2)*midpt[,3,i])/(sum(midpt[,3,i])) ## calculate the variance of Y using the  formula on page (10)
  
   
  MLEest[i,5]<- sum(((midpt[,1,i]-MLEest[i,1])*(midpt[,2,i]-MLEest[i,2]))*midpt[,3,i])/(sum(midpt[,3,i]))
  ### calculate the covariances of X & Y using the  formula on page (10)
  MLEest[i,6]<- MLEest[i,5]/sqrt(MLEest[i,3]*MLEest[i,4])
  ### calculate the correlation of X & Y using the  formula on page (10)
  
  
}
  
MLEest ## the MLE estimation matrix 
colMeans(MLEest) ## calculate the column means to see if the average results iver the 30 data sets are close to the default values
################################################################################
### Here for each parameter, we have calculated its mean, sd, and mean squarred error (MSE)  over the 30 simulated data sets
### and put the results of each patrameter on a vector, that means we have 5 parameters, so we will build 5 vectors of the 
### values of mean, sd and MSE of each parameters over the 30 data sets
### Then we produce a matrix of rbind (row bind) of these 5 vectors, which is a 5*3 matrix- 5 rows for each parameter, 
### 3 columns: for the calculated results of mean, sd and MSE 


#MLEest
MMX<- mean(MLEest[,1])
MSDX<- sd(MLEest[,1])
MMSEX<- (sum((MLEest[,1]-68)^2))/nrow(MLEest)
MX<- c(MMX,MSDX,MMSEX)
###################################################
MMY<- mean(MLEest[,2])
MSDY<- sd(MLEest[,2])
MMSEY<- (sum((MLEest[,2]-68)^2))/nrow(MLEest)
MY<- c(MMY,MSDY,MMSEY)

###################################################
MVX<- mean(MLEest[,3])
VSDX<- sd(MLEest[,3])
VMSEX<- (sum((MLEest[,3]-3)^2))/nrow(MLEest)
VX<- c(MVX,VSDX,VMSEX)

###################################################
MVY<- mean(MLEest[,4])
VSDY<- sd(MLEest[,4])
VMSEY<- (sum((MLEest[,4]-6)^2))/nrow(MLEest)
VY<- c(MVY,VSDY,VMSEY)

###################################################
MR<- mean(MLEest[,6])
RSD<- sd(MLEest[,6])
RMSE<- (sum((MLEest[,6]-0.4714045)^2))/nrow(MLEest)
RR<- c(MR,RSD,RMSE)

#Parameters<- c("Mean X","Mean Y","Var X","Var Y","Rho")
#criteria<- c("mean","std","MSE")
MLEBIV<- rbind(MX,MY,VX,VY,RR)
colnames(MLEBIV)<- c("mean","std","MSE")
rownames(MLEBIV)<- c("Mean X","Mean Y","Var X","Var Y","Rho")
print(MLEBIV)
#simulateddata
#library(xtable)
#xtable(MLEBIV,digits=7,caption = "Estimate of MLE Parameters for Bivariate Simulation")
#######################################################################################################################
### SIMULATE ONLY ONE BIVARIATE GROUPED DATA SET###
mm<- c(68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,2,6),2,2) ## The covariance matrix for the simulation

ssdata<- matrix(rep(0,1000*2),c(1000,2))## create an initial array of zero's

ssdata<- mvrnorm(n=1000,mm,ss) ## simulate data from Multivariate normal distribtion for 2 variables


x<- rep(0,1000) ## Initialize vector of X
y<- rep(0,1000) ## Initialize vector of Y


x<- ssdata[,1] ## assign the first column of each simulated array to the variable X
y<- ssdata[,2] ## assign the second column of each simulated array to the variable Y


Freqtable1<- array(rep(0,10*10),c(10,10)) ## create an array of the initial frequencies which is set to be zero

x.cut<- cut(x,breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf)) ## define the break points of the variable X for building 10 intervals 
y.cut<- cut(y,breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf)) ## define the break points of the variable Y for building 10 intervals 
Freqtable1<- table(x.cut,y.cut) ## Set the cut points to the frequency array, to count how many data points are 
### falling into the rectangles

simulateddata2<- array(rep(0,5*100),c(100,5)) ## To Create a Bivariate grouped data, first create an initialized array
### with the values of zero's
### This array is included a matrix (for the data set) of 100 (rows) and 5 (columns)- 100 is the number of total intervals
### 10 (intervals of X)* 10 (intervals of Y)- Note that the number of intervals could be different for 
### the variables, but for simplicity we consider the equal number of intervals 
### 5 is the number of columns of each dimension of this array: the column are as follows:
### Columns 1,3 (odd numbers): are the lower bounds of the variables X, Y 
### Columns 2,4 (even numbers): are the upper bounds of the variables X, Y 
### Columns 5: are the Frequency of the variables X, Y 
lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10) ## the lower bounds of X
lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10)) ## the lower bounds of Y
upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10) ## Upper bounds of X
upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10)) ## Upper bounds of Y
simulateddata2[,1]<- lower.x # Set the first column to the lower bound of X
simulateddata2[,3]<- lower.y # Set the third column to the lower bound of Y
simulateddata2[,2]<- upper.x # Set the second column to the upper bound of X
simulateddata2[,4]<- upper.y # Set the forth column to the upper bound of Y
simulateddata2[,5]<- c(Freqtable1)

#simulateddata2
###########################################################################################################
### TO FIND MLE OF THE PARAMETERS: Method 1
### To calculate the MLE's of the parameteres ignore the grouping :
### Here we use the midpoints of the intervals as observations
medx<- (simulateddata2[,1]+simulateddata2[,2])/2 ## to calcualte the midpoints of the variable X



for (i in 1:length(medx)){ ## to replace the -Inf & Inf with the proper number for X
  if (medx[i]==-Inf) medx[i]<- 63.5 ## replace -Inf with 63.5
  else if (medx[i]==Inf) medx[i]<- 72.5 ## replace Inf with 72.5
}


medy<- (simulateddata2[,3]+simulateddata2[,4])/2 ## to calcualte the midpoints of the variable Y


for (i in 1:length(medy)){ ## to replace the -Inf & Inf with the proper number for Y
  if (medy[i]==-Inf) medy[i]<- 63.7 ## replace -Inf with 63.4
  else if (medy[i]==Inf) medy[i]<- 72.7 ## replace Inf with 72.7
}


midpt<- array(rep(0, 100*3),c(100,3)) ## create an initial array of zero's for the midpoints of all intervals of the dataset
  midpt<- cbind(medx,medy,simulateddata2[,5]) ## create a columnbind matrix of 
  ### the midpoints of the variables and the frequencies
  colnames(midpt)<- c("medx","medy","Freq")## Rename the columns of the matrix

#midpt

MLEest<- rep(0,6) ## To calculate the MLE of the parameters, initialize a vecor of zero's with 
  ### 6 elements one for each parameter 
  names(MLEest)<- c("meanx","meany","sxx","syy","sxy","rho") ## name the elements of the vector

    MLEest[1]<- sum(midpt[,1]*midpt[,3])/(sum(midpt[,3])) ## calculate the mean of X using the  formula on page (10)
    MLEest[3]<- sum(((midpt[,1]-MLEest[1])^2)*midpt[,3])/(sum(midpt[,3])) ## calculate the variance of X using the  formula on page (10)
    
    MLEest[2]<- sum(midpt[,2]*midpt[,3])/(sum(midpt[,3])) ## calculate the mean of Y using the  formula on page (10)
    MLEest[4]<- sum(((midpt[,2]-MLEest[2])^2)*midpt[,3])/(sum(midpt[,3])) ## calculate the variance of Y using the  formula on page (10)
    
    
    MLEest[5]<- sum(((midpt[,1]-MLEest[1])*(midpt[,2]-MLEest[2]))*midpt[,3])/(sum(midpt[,3]))
    ### calculate the covariances of X & Y using the  formula on page (10)
    MLEest[6]<- MLEest[5]/sqrt(MLEest[3]*MLEest[4])
    ### calculate the correlation of X & Y using the  formula on page (10)
    
MLEest    
#########################################################################################
### TO FIND MLE OF THE PARAMETERS: Method 2

newsim<- array(rep(0,2*1000),c(1000,2)) ## Create an other matrix of zero's with 2 columns (number of variables) and 
### 1000 rows (number of observations)
### First column for the first variable (let's say X) and the second column for the other variable (Y)
  newsim[,1]<- rep(medx,simulateddata2[,5]) ### The first column for variable X, is compromising the replications of 
  ### midpoints of the intervals- the frequency of each midpoint equals the column of frequency of the simulated 
  ### bivariate grouped data created above  
  newsim[,2]<- rep(medy,simulateddata2[,5]) ### The second column for variable Y, is compromising the replications of 
  ### midpoints of the intervals- the frequency of each midpoint equals the column of frequency of the simulated 
  ### bivariate grouped data created above  
newsim ### If you run this code, this new version of our simulated data set, is very similar to the bivariate Galton data 
### when you read its original form of the data before organizing

#newsim[,,1]

MLEB1<- rep(0,5) ## create a new vector for the MLE's that we want to find from these new version of data
names(MLEB1)<- c("MX","MY","VX","VY","rXY") ## Name the elements of the vector
  
MLEB1[1]<- mean(newsim[,1]) ## set the first element of the vector equals mean of the variable X
MLEB1[2]<- mean(newsim[,2]) ## set the second element of the vector equals mean of the variable Y
MLEB1[3]<- sum((newsim[,1]-mean(newsim[,1]))^2)/nrow(newsim) ## set the third element of the vector equals variance of the
### variable X
### See how the calculation of the variance is different from the previous part, here we only find the difference between
### the values (which are midpoints, here as individual observations) and the mean of the variable, then squared it, and
### devide it by the total number of frequencies
###In method one, as these calculations is done for each interval, we have
### to multiply the result of sum of squares by the frequency of each interval and then devide it by the total number of 
### frequencies
### In method 2 as we replicate the midpoints of each interval with the numbers of frequency for the relevant interval, 
### we have use the origil mean method and do not need to multiply it by the number of frequecy of each interval.
### ****If this explanation is not clear, I can send you an screen shot of a data from frequency tables calculations of 
### the mean & variances

MLEB1[4]<- sum((newsim[,2]-mean(newsim[,2]))^2)/nrow(newsim) ## set the forth element of the vector equals variance of the
### variable X
MLEB1[5]<- ((sum((newsim[,1]-mean(newsim[,1]))*(newsim[,2]-mean(newsim[,2]))))/nrow(newsim))*(1/sqrt(MLEB1[3]*MLEB1[4])) ## set 
###the fifth element of the vector equals correlation between the variables  X and Y
  

MLEB1 ### See the result of this part **** (you can compare it with result of method 1, they should be the same)

MLEest

#######################################################################################################################
#### The metod 2 process for the 30 simulated data sets#### 
### *** I have described the process above, so just I have keep it to have all the parts related to bivariate simulations 
### in one specific R file

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