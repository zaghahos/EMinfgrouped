rm(list=ls())
set.seed(123456)
library(tmvtnorm) ## The libraries needed for the functions to be implemented: tmvtnorm, MASS
library(MASS) 

#####E-STEP####
###SIMULATE Z###
#############################################
ZMCEM1<- function(Data,mu,sigma){ ## This function to created to generate samples of data for the MCEM algorithm, 
  ### It has 3 arguments: Data=the data set we are working with, mu=the mean vector, sigme=the covariance matrix
  
  #k=5
  k<- 1000 ## Number of simulations that should be done by MCEM 
  n<- nrow(Data) ## to find the number of all rectangles
  d<- (ncol(Data)-1)/2 ## To find the number of dimension of data 
  
  ind<- ncol(Data)-1 ## select the number of columns of data needed to create the lower & upper bounds for the variables, to 
  ### select the index of which columns shuold be used for the lower bounds and which ones are for the upper bound
  lowerb<- NULL ## The initialized lowerbound incides of the variables (the column numbers that contain the lowerbound of the variables)
  upperb<- NULL ##The initialized upperbound incides of the variables (the column numbers that contain the upperbound of the variables)
  for (i in 1:ind){ ## set the loop to select the columns for the variables lower & upper bounds
    if (i%%2!=0) {lowerb<- append(lowerb,i)} ## to select the lower bounds columns
    else {upperb<- append(upperb,i)} ## to select the upper bound columns
    
  }
  Data1<- as.matrix(Data) ## As mtmvnorm function, needed the numerics for the lower & upper bounds and the results of the previous
  ### lines of code were in the format of lists, so we have to transform data to the matrix to have the results as numeric/vectors
  ### otherwise, we will get error when running the pmvnorm function.
  
  m<- array(rep(0,k*d*n),c(k,d,n)) ## Create an initial array of zero's for the simulated data 
  #print(m)
  
  for(i in 1:n){ ## Use the for loop to do the simulation over all rows
    m[,,i]<- rtmvnorm(k, mean=mu,
                      sigma=sigma, lower=Data1[i,c(lowerb)], upper=Data1[i,c(upperb)],
                      algorithm="gibbs") ## Simulate the data using K: number of data to be generated, mean: the mean vector, 
    ### sigma: Sigma vector, lower: lower bounds, upper: upper bounds, algorithm: which here "gibbs" should be selected 
    #  print(m)
  }
  #print(m)
  
  #mm<- matrix(rep(0,2*n),ncol=2)
  #print(mm)
  #for(i in 1:n){
  #mm[i,]<- apply(m[,,i],2,mean)
  #}
  #print(mm)
  return(m)## Return the array of simulations
  
}

## d=2##
mu1 <- c(67,67)
sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
## d=3##
mu2<- c(67,67,67)
sigma2<- matrix(c(3.1,2.2,1.6,2.2,4.3,2.45,1.6,2.45,4.8),3,3)
#sigma2

out1<- ZMCEM1(Data=DataBiv,mu1,sigma1)
out2<- ZMCEM1(Data=simulateddata2,mu1,sigma1)
out3<- ZMCEM1(Data=simulateddata3,mu2,sigma2)

###E-step: Mean Calculation###
######################################################
Mu<- function(Data,Y){ ## This is the function for calculation mean vector from the generated samples in the previous  part
  ###The arguments of this function are Data (original data that we have) and Y= which is the 
  ### simulated (generated) data or samples we have done in the previous function ZMCEM1 
  n<- nrow(Data) ## set n equals number of row of data (equals number of rectangles/surfaces) 
  F<- ncol(Data) ## number of column of dataset
  T<- sum(Data[,F]) ## Total sum of observations (Data)
  #print(T)
  d<- (ncol(Data)-1)/2 ## To find the number of dimension of data 
  mm<- matrix(rep(0,d*n),ncol=d) ## create & initialize a matrix of zero's with d(=number of variables) columns 
  #print(mm)
  for(i in 1:n){ ## use for loop to implement over all rectangles/surfaces
    mm[i,]<- apply(Y[,,i],2,mean) ## using apply function to sum over the columns of Y's (which is the array generated 
    ### using the ZMCEM1 function above), and assign the results to each row of the mm matrix created before  
  }
  #print(mm)
  Ym<- mm*Data[,F] ## multiply the sumation of each rows by the frequencies of each relevant row and assign it to Ym
  #print(Ym)
  sm<- rep(0,d) ## create and initialize a vector called sm of zero's with d(=number of dimensions/variables)
  for(i in 1:d){ ## use for loop to implement over each element of the vector
    sm[i]<- sum(Ym[,i]) ## sum over the columns of Ym and assign it to each elements of sm
  }
  #print(sm)
  mu<- (sm/T) ## devide each element of sm by the total number of observations to calculate the mean of each variable
  return(mu) ## return the vector of means
}

#out1<- ZMCEM1(Data=DataBiv,mu1,sigma1)
#out2<- ZMCEM1(Data=simulateddata2,mu1,sigma1)
#out3<- ZMCEM1(Data=simulateddata3,mu2,sigma2)

out4=Mu(Data=DataBiv,Y=ZMCEM1(Data=DataBiv,mu1,sigma1))
out5=Mu(Data=simulateddata2,Y=ZMCEM1(Data=simulateddata2,mu1,sigma1))
out6=Mu(Data=simulateddata3,Y=ZMCEM1(Data=simulateddata3,mu2,sigma2))

### E-Step: Sigma Calculation###
#######################################################
#######################################################################
### SIGMA###
Sigma1<- function(Data,Y,M){## Using the above two functions, in this function to find the covariance matrix, 
  ### The arguments of this function are: Data= the data set we are working with, Y= the generated (simulated) samples 
  ### from the function ZMCEM1, M= the mean vector from the Mu function above
  n<- nrow(Data) ## set n equals number of row of data (equals number of rectangles/surfaces) 
  F<- ncol(Data) ## number of column of dataset
  total<- sum(Data[,F]) ## Total sum of observations (Data)
  #print(Y)
  #k=5
  k=1000 ## number of generated samples
  d<- (ncol(Data)-1)/2 ## To find the number of dimension of data 
  myarray<- array(rep(0,d*d*k*n),c(d,d,k,n)) ## create and initialize an array of zero's with the dims: d,d,k,n to create 
  ### d*d covariance matrices over the generated samples and the rectangles/surfaces
  for(i in 1:n){ ## use for loop to implement over all rows (rectangles/surfaces)
    for(j in 1:k){ ## use another for loop to implement over all generated samples
      myarray[,,j,i]<- (Y[j,,i]-M)%*%t(Y[j,,i]-M) ## over the array of generated samples, for each row, and each generated 
      ### samples, we have to find the difference between the generated samples and the mean vectors and square it by 
      ### multiply it to its transpose
    }
    
  }
  #print(myarray)
  #m<- apply(myarray[,,,1],c(1,2),mean)
  #print(m)
  Newarray<- array(rep(0,d*d*n),c(d,d,n)) ## create and initialize an array of zero's for the covariances over the 
  ### rectangles/surfaces
  #print(Newarray)
  for(i in 1:n){ ## use for loop to do the function over all rectangles/surfaces
    Newarray[,,i]<- apply(myarray[,,,i],c(1,2),mean) ## take mean over the first and second dimensions of array from last
    ### step (myarray) and assign the results for each rectangles of Newarray  
  }
  #print(Newarray)
  
  NNarray<- array(rep(0,d*d*n),c(d,d,n))  ## create and initialize an array of zero's for the covariances over the 
  ### rectangles/surfaces
  for(i in 1:n){ ## use for loop to do the function over all rectangles/surfaces
    NNarray[,,i]<- Newarray[,,i]*Data[i,F] ## Multiply the covariance matrices over the rectangles/surfaces
    ### from above (Newaaray) by the frequency of each surface
  }
  #print(NNarray)
  S1<- apply(NNarray,c(1,2),sum) ## sum over all the rectangles/surfaces
  #print(S1)
  mcov<- S1/total ## devide the result by total number of observations to find the covariance matrix
  #print(mcov)
  #print(myarray)
  #m1<- array(rep(0,2*2*n),c(2,2,n))
  #m2<- array(rep(0,2*2*n),c(2,2,n))
  #for(i in 1:n){
  # m1[,,i]<- (Y[i,]-M)%*%t(Y[i,]-M)
  #m2[,,i]<- m1[,,i]*Data[i,5]
  
  #  }
  #print(m1)
  #print(m2)
  # m3<- apply(m2,1:2,sum)
  #print(m3)
  #mcov<- m3/total
  #print(mcov)
  return(mcov) ## return the covariance matrix
}


#out4=Mu(Data=DataBiv,Y=ZMCEM1(Data=DataBiv,mu1,sigma1))
#out5=Mu(Data=simulateddata2,Y=ZMCEM1(Data=simulateddata2,mu1,sigma1))
#out6=Mu(Data=simulateddata3,Y=ZMCEM1(Data=simulateddata3,mu2,sigma2))


out7=Sigma1(Data=DataBiv,Y=ZMCEM1(Data=DataBiv,mu1,sigma1),
            M=Mu(Data=DataBiv,Y=ZMCEM1(Data=DataBiv,mu1,sigma1)))

out8=Sigma1(Data=simulateddata2,Y=ZMCEM1(Data=simulateddata2,mu1,sigma1),
            M=Mu(Data=simulateddata2,Y=ZMCEM1(Data=simulateddata2,mu1,sigma1)))

out9=Sigma1(Data=simulateddata3,Y=ZMCEM1(Data=simulateddata3,mu2,sigma2),
            M=Mu(Data=simulateddata3,Y=ZMCEM1(Data=simulateddata3,mu2,sigma2)))

### M-Step###
MCEMBiv<- function(data,mu_init,sigma_init,maxit=1000,tol1=1e-3,tol2=1e-3){## This is the function of M-step, to iterate 
  ### over the results of the previuos parts and find the maximum updating estimates. It has 5 arguments, 
  ### data=the data set we are working with, mu_init=the initialized mean vector, sigma_init=the initialize covariance matrix,
  ### tol1 & tol 2: the stopping criterias for means vectors and covariance matrix 
  flag<- 0 ## initlize flag at zero for the stopping rule 
  M_cur<- mu_init ## rename the initialized mean vector
  S_cur<- sigma_init ## rename the initialized covariance matrix
  #theta_cur<- c(Mu_cur,S_cur)
  iter<- rep(0,maxit) ## create a vector of iter to record the number of iteration 
  #Svec<- rep(0,maxit)
  #Mvec<- rep(0,maxit)
  for (i in 1:maxit){ ## use the for loop to do all the steps by maximum number of the iterations
    #print(paste("Iteration number=", i))
    #cur<- c(Mu_cur,S_cur)
    MU_cur<- M_cur ## rename the mean vector inside the for loop
    SS_cur<- S_cur ## rename the covariance matrix inside the for loop
    
    MuNew<- Mu(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur)) ## using the Mu fuction to create the updated 
    ### mean vector and assign it to MuNew
    
    SigmaNew<- Sigma1(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur),
                      M=Mu(Data=mydat,Y=ZMCEM1(Data=mydat,mu=MU_cur,sigma=SS_cur) ) ) ## using the Sigma1 fuction to create
    ###the updated covariance matrix and assign it to SigmaNew
    
    Mu_new<- MuNew ## rename the MuNew
    Sigma_new<- SigmaNew ## rename the SigmaNew
    
    diff1<- MU_cur-Mu_new ## find the differences between the current mean vector and the new estimate of mean vector 
    diff2<- SS_cur-Sigma_new ## find the differences between the current covariance matrix and the new estimate of covariance matrix
   
    D1<- abs(mean(diff1)) ## calculate the absolute value of the mean of differences in means
    D2<- abs(mean(diff2)) ## calculate the absolute value of the mean of differences 
    ##in the elements of the covariance matrix
    
    #print(D1)
    #print(D2)
    if(D1<tol1 & D2<tol2){flag<- 1 ;break} ## define the condition for the absolute values of the differences according to the 
    ## stopping criteria, continue updating until the values of the parameter estimates meet both conditions, then stop and 
    ## take the final new estimates as the updated estimates of the parameters 
    
    M_cur<- Mu_new ## update the current value of mean vector with the new estimate
   
    S_cur<- Sigma_new ## update the current value of covariance matrix with the new estimate
  }
  
  if(!flag) warning("Didn't Converge \n") ## the condition if there is no convergence in the stopping rule
  #list(paste("updated Mu of the EM=" , M_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #    paste("iteration number=",i))
  
  
  updatedest=list(M_cur,(S_cur)) ## the list of the updated estimates of the parameters
  return(updatedest) ## return the list of updated estimates of the parameters
  
}
#####################################################################################
### Check if the function works on data ###

###INITIAL VALUES for dim=2 & dim =3####
## d=2##
mu1 <- c(67,67)
sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
## d=3##
mu2<- c(67,67,67)
sigma2<- matrix(c(3.1,2.2,1.6,2.2,4.3,2.45,1.6,2.45,4.8),3,3)
#sigma2
################################################################
### Implement the function on Galton data###
###DATA###
DataBiv<- read.csv("C:/Users/Zahra/Desktop/Galton Data/Galton-Bivariate.csv") ## read the data in the form of grouped data (bivariate)

DataBiv[1:14,1]<- -Inf ## set the first lower 14 limits of parent data equals -inf (you can see the table of 
##data to understand about the intervals and see the gropued form of data using that table )   
DataBiv[141:154,2]<- Inf ## set the last 14 limilts of children equals -inf 
DataBiv[c(1,15,29,43,57,71,85,99,113,127,141),3]<- -Inf ##set the last upper limit of the last intervals (parent data) equals inf  
DataBiv[c(14,28,42,56,70,84,98,112,126,140,154),4]<- Inf ##set the last upper limit of the last intervals (children data) equals inf

#########################################################################
mydat<- DataBiv ## Assign mydat to the Galton data

out10=MCEMBiv(data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-3,tol2=1e-3)

rm("mydat", "DataBiv")

#############################################################################################################
### Implement the function on a smulation of Bivariate grouped data 
### Simulate the data d=2###
mm<- c(68,68)
ss<- matrix(c(3,2,2,6),2,2)

ssdata<- matrix(rep(0,1000*2),c(1000,2))

ssdata<- mvrnorm(n=1000,mm,ss)

#ssdata

x<- rep(0,1000)
y<- rep(0,1000)


x<- ssdata[,1]
y<- ssdata[,2]

Freqtable1<- array(rep(0,10*10),c(10,10))

x.cut<- cut(x,breaks=c(-Inf,64,65,66,67,68,69,70,71,72,Inf))
y.cut<- cut(y,breaks=c(-Inf,64.2,65.2,66.2,67.2,68.2,69.2,70.2,71.2,72.2,Inf))
Freqtable1<- table(x.cut,y.cut)


simulateddata2<- array(rep(0,5*100),c(100,5))
lower.x<- rep(c(-Inf,64,65,66,67,68,69,70,71,72),10)
lower.y<- c(rep(-Inf,10),rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10)
            ,rep(71.2,10),rep(72.2,10))
upper.x<- rep(c(64,65,66,67,68,69,70,71,72,Inf),10)
upper.y<- c(rep(64.2,10),rep(65.2,10),rep(66.2,10),rep(67.2,10),rep(68.2,10),rep(69.2,10),rep(70.2,10),
            rep(71.2,10),rep(72,2,10),rep(Inf,10))
simulateddata2[,1]<- lower.x
simulateddata2[,3]<- lower.y
simulateddata2[,2]<- upper.x
simulateddata2[,4]<- upper.y
simulateddata2[,5]<- c(Freqtable1)

#simulateddata2
##########################################################################################
### run the EM function###
mydat<- simulateddata2 ## Assign mydat to the simulated data

out11=MCEMBiv(data=simulateddata2,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-3,tol2=1e-3)

out11 ## see the output
unlist(out11) ## unlist the output

rm("mydat","simulateddata2")
##############################################################################################
### Implement the function on a smulation of Trivariate grouped data 
### Simulate the data d=3###
#####################################################################################
mm<- c(68,68,68) ## The vector of mu that we want to simulate from
ss<- matrix(c(3,2,1.5,2,4,2.5,1.5,2.5,5),3,3) ## The covariance matrix for the simulation

ssdata<- array(rep(0,3000*3),c(3000,3)) 
ssdata<- mvrnorm(n=3000,mm,ss) 

x<- rep(0,3000) 
y<- rep(0,3000) 
z<- rep(0,3000) 

x<- ssdata[,1] 
y<- ssdata[,2] 
z<- ssdata[,3] 

Freqtable1<- array(rep(0,8*8*8),c(8,8,8)) 
x.cut<- cut(x,breaks=c(-Inf,65,66,67,68,69,70,71,Inf))  
y.cut<- cut(y,breaks=c(-Inf,63,65,67,69,71,73,75,Inf)) 
z.cut<- cut(z,breaks=c(-Inf,64,65.5,67,68.5,70,71.5,73,Inf)) 

Freqtable1<- table(x.cut,y.cut,z.cut) 


simulateddata3<- array(rep(0,7*512),c(512,7)) 
lower.x<- rep(c(-Inf,65,66,67,68,69,70,71),64) 
ly<- c(rep(-Inf,8),rep(63,8),rep(65,8),rep(67,8),rep(69,8),rep(71,8),rep(73,8),rep(75,8)) 
lower.y<- rep(ly,8) 
lower.z<- c(rep(-Inf,64),rep(64,64),rep(65.5,64),rep(67,64),rep(68.5,64),rep(70,64),rep(71.5,64),rep(73,64)) 

upper.x<- rep(c(65,66,67,68,69,70,71,Inf),64) 
uy<- c(rep(63,8),rep(65,8),rep(67,8),rep(69,8),rep(71,8),rep(73,8),rep(75,8),rep(Inf,8)) 
upper.y<- rep(uy,8) 
upper.z<- c(rep(64,64),rep(65.5,64),rep(67,64),rep(68.5,64),rep(70,64),
            rep(71.5,64),rep(73,64),rep(Inf,64)) 

simulateddata3[,1]<- lower.x 
simulateddata3[,3]<- lower.y 
simulateddata3[,5]<- lower.z  
simulateddata3[,2]<- upper.x 
simulateddata3[,4]<- upper.y 
simulateddata3[,6]<- upper.z 
simulateddata3[,7]<- c(Freqtable1) 
#simulateddata3

#############################################################################
### Run the function###
mydat<- simulateddata3 ## Assign mydat to the trivariate simulated data

out12=MCEMBiv(data=mydat,mu_init=mu2,sigma_init=sigma2,maxit=1000,tol1=1e-3,tol2=1e-3)

out12 ## see the output
unlist(out12) ## unlist the output


######################################################################################
#### TO SEE THE RESULT of MCEM FUNCTION OVER 30 SIMULATED DATA SET

mu1 <- c(67,67)
sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
############################################################
####n=50####
MCEMoutsim50<- matrix(rep(0*5*30),ncol=5,nrow=30)
colnames(MCEMoutsim50)<- c("mux1","mux2","Sx1","Sx2","rhox1x2")
for(i in 1:nrow(MCEMoutsim50)){
  mydat<- simulateddata50[,,i]
  MCEMoutsim50[i,]<- MCEMBiv(data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-3,tol2=1e-3)
  
  
}

###########################################################
####n=150####
MCEMoutsim150<- matrix(rep(0*5*30),ncol=5,nrow=30)
colnames(MCEMoutsim150)<- c("mux1","mux2","Sx1","Sx2","rhox1x2")
for(i in 1:nrow(MCEMoutsim150)){
  mydat<- simulateddata150[,,i]
  MCEMoutsim150[i,]<- MCEMBiv(data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-3,tol2=1e-3)
  
  
}

###########################################################
####n=500####
MCEMoutsim500<- matrix(rep(0*5*30),ncol=5,nrow=30)
colnames(MCEMoutsim500)<- c("mux1","mux2","Sx1","Sx2","rhox1x2")
for(i in 1:nrow(MCEMoutsim500)){
  mydat<- simulateddata500[,,i]
  MCEMoutsim500[i,]<- MCEMBiv(data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-3,tol2=1e-3)
  
  
}

##########################################################
####n=1000####
MCEMoutsim<- matrix(rep(0*5*30),ncol=5,nrow=30)
colnames(MCEMoutsim)<- c("mux1","mux2","Sx1","Sx2","rhox1x2")
for(i in 1:nrow(MCEMoutsim)){
  mydat<- simulateddata[,,i]
  MCEMoutsim[i,]<- MCEMBiv(data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-3,tol2=1e-3)
  
  
}
#MCEMoutsim
#############################################################################
###n=50###
#MLEest
MMX50<- mean(MLEest50[,1])
MSDX50<- sd(MLEest50[,1])
MMSEX50<- (sum((MLEest50[,1]-68)^2))/nrow(MLEest50)
MX50<- c(MMX50,MSDX50,MMSEX50)
###################################################
MMY50<- mean(MLEest50[,2])
MSDY50<- sd(MLEest50[,2])
MMSEY50<- (sum((MLEest50[,2]-68)^2))/nrow(MLEest50)
MY50<- c(MMY50,MSDY50,MMSEY50)

###################################################
MVX50<- mean(MLEest50[,3])
VSDX50<- sd(MLEest50[,3])
VMSEX50<- (sum((MLEest50[,3]-3)^2))/nrow(MLEest50)
VX50<- c(MVX50,VSDX50,VMSEX50)

###################################################
MVY50<- mean(MLEest50[,5])
VSDY50<- sd(MLEest50[,5])
VMSEY50<- (sum((MLEest50[,5]-6)^2))/nrow(MLEest50)
VY50<- c(MVY50,VSDY50,VMSEY50)

#####################################################################
MR50<- mean(MLEest50[,6])
RSD50<- sd(MLEest50[,6])
RMSE50<- (sum((MLEest50[,6]-0.4714045)^2))/nrow(MLEest50)

RR50<- c(MR50,RSD50,RMSE50)

#Parameters<- c("Mean X","Mean Y","Var X","Var Y","Rho")
#criteria<- c("mean","std","MSE")
MLEBIV50<- rbind(MX50,MY50,VX50,VY50,RR50)
colnames(MLEBIV50)<- c("mean","std","MSE")
rownames(MLEBIV50)<- c("Mean X","Mean Y","Var X","Var Y","Rho")
######################################################################
###n=150####
#MLEest
MMX150<- mean(MLEest150[,1])
MSDX150<- sd(MLEest150[,1])
MMSEX150<- (sum((MLEest150[,1]-68)^2))/nrow(MLEest150)
#MX150<- c(MMX150,MSDX150,MMSEX150)
###################################################
MMY150<- mean(MLEest150[,2])
MSDY150<- sd(MLEest150[,2])
MMSEY150<- (sum((MLEest150[,2]-68)^2))/nrow(MLEest150)
MY150<- c(MMY150,MSDY150,MMSEY150)

###################################################
MVX150<- mean(MLEest150[,3])
VSDX150<- sd(MLEest150[,3])
VMSEX150<- (sum((MLEest150[,3]-3)^2))/nrow(MLEest150)
VX150<- c(MVX150,VSDX150,VMSEX150)

###################################################
MVY150<- mean(MLEest150[,5])
VSDY150<- sd(MLEest150[,5])
VMSEY150<- (sum((MLEest150[,5]-6)^2))/nrow(MLEest150)
VY150<- c(MVY150,VSDY150,VMSEY150)

###################################################
MR150<- mean(MLEest150[,6])
RSD150<- sd(MLEest150[,6])
RMSE150<- (sum((MLEest150[,6]-0.4714045)^2))/nrow(MLEest150)

RR150<- c(MR150,RSD150,RMSE150)

#Parameters<- c("Mean X","Mean Y","Var X","Var Y","Rho")
#criteria<- c("mean","std","MSE")
MLEBIV150<- rbind(MX150,MY150,VX150,VY150,RR150)
colnames(MLEBIV150)<- c("mean","std","MSE")
rownames(MLEBIV150)<- c("Mean X","Mean Y","Var X","Var Y","Rho")
#print(MLEBIV50)

#######################################################################
###n=500###
#MLEest
MMX500<- mean(MLEest500[,1])
MSDX500<- sd(MLEest500[,1])
MMSEX500<- (sum((MLEest500[,1]-68)^2))/nrow(MLEest500)
MX500<- c(MMX500,MSDX500,MMSEX500)
###################################################
MMY500<- mean(MLEest500[,2])
MSDY500<- sd(MLEest500[,2])
MMSEY500<- (sum((MLEest500[,2]-68)^2))/nrow(MLEest500)
MY500<- c(MMY500,MSDY500,MMSEY500)

###################################################
MVX500<- mean(MLEest500[,3])
VSDX500<- sd(MLEest500[,3])
VMSEX500<- (sum((MLEest500[,3]-3)^2))/nrow(MLEest500)
VX500<- c(MVX500,VSDX500,VMSEX500)

###################################################
MVY500<- mean(MLEest500[,5])
VSDY500<- sd(MLEest500[,5])
VMSEY500<- (sum((MLEest500[,5]-6)^2))/nrow(MLEest500)
VY500<- c(MVY500,VSDY500,VMSEY500)

###################################################
MR500<- mean(MLEest500[,6])
RSD500<- sd(MLEest500[,6])
RMSE500<- (sum((MLEest500[,6]-0.4714045)^2))/nrow(MLEest500)

RR500<- c(MR500,RSD500,RMSE500)

#Parameters<- c("Mean X","Mean Y","Var X","Var Y","Rho")
#criteria<- c("mean","std","MSE")
MLEBIV500<- rbind(MX500,MY500,VX500,VY500,RR500)
colnames(MLEBIV500)<- c("mean","std","MSE")
rownames(MLEBIV500)<- c("Mean X","Mean Y","Var X","Var Y","Rho")
#print(MLEBIV50)
#simulateddata

