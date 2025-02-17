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
######################################################################################
#### TO SEE THE RESULT of MCEM FUNCTION OVER 30 SIMULATED DATA SET

mu1 <- c(67,67)
sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
############################################################
####n=1000####
#### MCEM1000b10: MEEM RESULTS FOR N=1000 & Binns=10
MCEM1000b10<- matrix(rep(0*6*10),ncol=6,nrow=10)
colnames(MCEM1000b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")
for(i in 1:nrow(MCEM1000b10)){
  mydat<-  simulateddata1000b10[,,i] ###JOHN: HERE NOTE THAT THE PROPER SIMULATED DATA SHOULD BE REPLACED
  MCEM1000b10[i,]<- unlist(MCEMBiv(data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-3,tol2=1e-3))
  
  
}

OUTMCEM1000b10<- write.csv(MCEM1000b10,"C:/Users/sh_za/Desktop/Results/Bivariate/MCEM/Sample size1000/MCEM1000b10.csv")

###########################################################
