## This part is based on section 2.2.3 & 2.2.5 EM algorithm
## In the first part, we have find and define the functions for the estimate of the parameters based on formulas 19-23,
##and the equations on page 15,
## then in the second part, we iterate over the current estimates of the parameters to find the most updated estimates based on 
##M-step for the parameters

###Parameter Estimation Using EM Algorithm####
#######################################################################
###E-STEP:####
###Moments###
MEXI<- function(Data,mu,sigma){ ## the function to find the first moment of the truncated bivariate normal distributions for 
  ## the rectangles which will be used to calculate the estimate of the covariance matrix 
  ##it has three arguments:  Data: in the form of bivariate grouped data, mu: mean vector, sigma: covariance matrix
  n<- nrow(Data) ## to find the number of all rectangles
  d<- (ncol(Data)-1)/2 ## To find the number of dimension of data 
  
  ind<- ncol(Data)-1 ## select the number of columns of data needed to create the lower & upper bounds for the variables, to 
  ### select the index of which columns should be used for the lower bounds and which ones are for the upper bound
  lowerb<- NULL ## The initialized lower bound incides of the variables (the column numbers that contain the lowerbound of the variables)
  upperb<- NULL ##The initialized upper bound incides of the variables (the column numbers that contain the upperbound of the variables)
  for (i in 1:ind){ ## set the loop to select the columns for the variables lower & upper bounds
    if (i%%2!=0) {lowerb<- append(lowerb,i)} ## to select the lower bounds columns
    else {upperb<- append(upperb,i)} ## to select the upper bound columns
    
  }
  Data1<- as.matrix(Data) ## As mtmvnorm function, needed the numerics for the lower & upper bounds and the results of the previous
  ### lines of code were in the format of lists, so we have to transform data to the matrix to have the results as numeric/vectors
  ### otherwise, we will get error when running the pmvnorm function.
  
  
  
  EX<- matrix(rep(0,n*d),ncol=d,nrow=n) ## define an initialized matrix of zeros to put the mean vectors in it
  for(i in 1:nrow(Data)){ ## use for loop to do the iterations over all rectangles/surfaces
    mnts<- mtmvnorm(mean=mu, sigma=sigma, 
                    lower=Data1[i,c(lowerb)],
                    upper=Data1[i,c(upperb)]) ## to find the first moments vectors for each rectangle, we have to use the mtmvnorm
    ## function from the library of tmvtnorm, the lower: is the lower limit of rectangle, the upper: is the upper limit of rectangle
    ## for each separated interval
    EX[i,]<- mnts$tmean ## we select the tmeans (to get only the first moment vectors results from the mnts' output) and put the 
    ## resulted mean vectors for each rectangle on each raw of the matrix
  }
  return(EX) ## return the first moment vector matrix of the rectagles (the expectation part of mu on page 13, 16 
  ### for the rectangles/surfaces)
  
}


mu1 <- c(67,67)
sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)



#############################################################
McovXI<- function(Data,mu,sigma){ ## the function to find the covariance matrices using the current theta of the truncated 
  ##bivariate normal distributions for the rectangles which will be used to calculate the estimate of the covariance matrix 
  ##it has three arguments:  Data: in the form of bivariate grouped data, mu: mean vector, sigma: covariance matrix
  #### Attention: in this part, the covariance is only found from the mtmvnorm function for the current state of theta, 
  #### however, in the EM algorithm, to calculate the estimate of the covariance matrix, we have to use the 
  #### updated means/mean vectors to find the new estimate for covariances, so the result of this part, will be used
  #### indirectly, to calculate the estimate of the covariances (If you have any question regarding this part, I can explain it
  #### to Camila, and she will explain it simply in your language to be easily understandable for you )/
  #### also see the detail on the following part related to function with name: EXXest
  n<- nrow(Data) ## to caclulate the number of all rectangles/surfaces
  d<- (ncol(Data)-1)/2 ## for the following lines, see the previous part
  
  ind<- ncol(Data)-1
  lowerb<- NULL
  upperb<- NULL
  for (i in 1:ind){ 
    if (i%%2!=0) {lowerb<- append(lowerb,i)} 
    else {upperb<- append(upperb,i)} 
    
  }
  Data1<- as.matrix(Data)
  
  s1<- array(rep(0,n*d*d),c(d,d,n)) ## define and assigned an initialized array for the covariance matrices using the current 
  ##theta for all the rectangles/ surfaces
  for(i in 1:nrow(Data)){ ## use for loop to iterate the calcultion over the rectanlges
    mnts<- mtmvnorm(mean=mu, sigma=sigma,
                    lower=Data1[i,c(lowerb)],
                    upper=Data1[i,c(upperb)]) ## again to find the covariance matrices using the current theta for 
    ##each rectangle, we have to use the mtmvnorm
    ## function from the library of tmvtnorm, the lower: is the lower limit of rectangle, the upper: is the upper limit of rectangle
    ## for each separated interval
    s1[,,i]<- mnts$tvar ## assing the tvar output from the above results to the i element of the defined array for the rectangles/surfaces
  }
  return(s1) ## return the s1, the array of the covariance matrices (using the current theta) results (This is used in the expectation in the formula of 
  ## the estimate of the covariance matrix on page 13,16)
  
}



#############################################################
###MU:E-step###
MEM<- function(Data,EXi){ ## the function to calculate the  mean vectors , formula is on page 15,
  ##or the for the bivariate case each element of the vector are eaxctly according to the formulas 19,20 
  ## this function has two arguments: Data in the form of bivariate grouped data, EXi: the first moments matrix of the
  ## rectangles from above
  F<- ncol(Data) ## To find the column number of the frequencies
  Freq<- Data[,F] ## The vector of frequencies over the rectangles/surfaces
  A<- EXi*Data[,F] ## the numerator of the estimate of mean vector (see page 13 (also 16), also for the bivariate case, the elementes
  ## of the vecor are exactly the numerator of formulas 19, 20
     
  mupd<- apply(A,2,sum)/sum(Data[,F]) ## the estimate of mean vector (see page 13, 16, also for the bivariate case, the elementes
  ## of the vecor are exactly the estimate of the means of formulas 19, 20)
  return(mupd) ## return the mean vector
  
}

################################################################################
###E(XX) estimate###
EXXest<- function(Data,EXi,ss1){ ## the function to calculate the second moments, using the formulas is on page 15,
  ##or the for the bivariate case each element of the vector are eaxctly according to the formulas 20,21,22 
  ## this function has three arguments: Data in the form of bivariate grouped data, EXi: the first moments matrix of the
  ## rectangles from above, ss1: the array of covariances of the rectangles from above 
  n<- nrow(Data)## the number of all rectangles/surfaces
  d<- (ncol(Data)-1)/2 ## To find the number of variables inorder to illustrate the covariance matrices for the same dimensions
  ExxNew<- array(rep(0,n*d*d),c(d,d,n)) ## define and assigne an initialized array, where the elements of the array, are the 
  ## second moments for the rectangles
  for(i in 1:n){ ## use for loop to do all the computations over all rectangles
    ExxNew[,,i]<- ss1[,,i]+EXi[i,]%*%t(EXi[i,]) ## using the estimated first moments and the covariances (for the current theta)
    ## to find the second moments matrices
    
    #### Explanation to clarify this part: If you see the formula of the expectation for calculating the estimate of the 
    #### covariance matrix (page 13 (16), also for bivariate: formulas 21,22,23), we have to extend the expectation exactly, 
    #### and E[(x-muupdate)^t * (x-muupdate)]=E[x^t*x]-muupdate^t*E[x]-E[x^t]*muupdate+muupdate^t*muupdate, where the
    #### expectation is done over the current state of theta
    #### in the extended form above, E[x^t*x]: is the second moment that is used in this formula (using the current state of 
    #### theta) and replace in this summation to find the required expectation in the formulas (page 13 (16),formulas 21,22,23)
    #### in the functions defined above, from mtmvnorm output over the rectangles, we can find the first moments and the 
    #### covariances where the covariances are not updated version, however, we use those covariances to find the 
    #### second moments (on the extended form of the expectations) under the current state of theta's (parameters)
    
  }
  return(ExxNew) ## return the array where its elements are the second moment matrics for the current theta 
  
}




######################################################################
Mysig<- function(Data,EXX,EX,UPDMU){ ## the function to calculate the estimate of the covariances (formula on page 13 & 16, also 
  ## for bivariate case, each resulted elements are in the form of formulas 21,22, 23)
  ## this function has four arguments: Data in the form of bivariate grouped data, EX: the first moments matrix of the
  ## rectangles from above, EXX: the array of second moments matrices over the rectangles from above, UPDMU: the updated mu
  ## from the above functions
  n<- nrow(Data) ## the number of rectangles/surfaces 
  d<- (ncol(Data)-1)/2 ## the number of dimensions
  s2<- array(rep(0,n*d*d),c(d,d,n)) ## define and assign an array to put the estimate of the covariance matrices
  
  for(i in 1:n){ ## use for loop to do the computations over all rectangles/surfaces

  s2[,,i]<- (EXX[,,i]-UPDMU%*%t(EX[i,])-EX[i,]%*%t(UPDMU)+(UPDMU%*%t(UPDMU)))*Data[i,ncol(Data)] ## Here, we use the extended form 
  ## of the expectations for each rectangle and multiply it by the frequency of that relevant rectangles/surfaces
  
  }
  
  SigmaNew<- apply(s2,c(1,2),sum)/sum(Data[,ncol(Data)]) ## use the fuction apply to find the summation in the numerator of the 
  ## formulas and devide it by total frequency in the bivariate grouped data
  
  if (ncol(SigmaNew)>2){ ## As the covariances produced from the the output of the function tmvtnorm$var, in the 
    ### function McovxI (EXX: argumant here) is not symmetric exactly for the dimension of more than one (maybe not exact
    ### in the very last digits), so we have to make sure if the matrix is symmetric, we use the condition, to check it for
    ### dimensions above d>2
    for (i in 1:ncol(SigmaNew)){ ## use for loop to check over all elements
      if (SigmaNew[upper.tri(SigmaNew)][i]!=SigmaNew[lower.tri(SigmaNew)][i]){ ## check if the elements of the upper diagonal
        ### are not equal with the lower diagonal, then set (replace) all the elements of the lower diagonal with the 
        ### elements of the upper diagonal. As they are not equal in the vary last digits mostly after the second digits, 
        ### so their exact difference are very little, however, as the output of this step will be used in M-step for iterate
        ### and update the covariance matrix for the EM algorithm, we need to make sure this is a symmetric matrix to avoid
        ### any further errors
        SigmaNew[lower.tri(SigmaNew)][i]<- SigmaNew[upper.tri(SigmaNew)][i] ## replace the values of the lower diagonal 
        ### elements with the values of the upper diagonal elements  
      }
    }
    
  }
  
  
  return(SigmaNew) ## Return the estimate of the covariance matrix
}


 
  
########################################################################
########################################################################
###M-STEP:###
### M-Step###
EMBivG<- function(Data,mu_init,sigma_init,maxit=1000,tol1=1e-4,tol2=1e-3){ ## The EM function, to find the estimate of 
  ## the mean vector/means and covariance matrix
  ## This function has 6 elements: 
  ## Data: the grouped bivariate data, mu_iniit: the initial values of mu vector, sigma_init: the initial covariance matrix,
  ## tol1: stopping criteria for means, tol2: stopping criteria for variances
  flag<- 0 ## the initial value of flag which will be used in the stopping criteria conditions
  M_cur<- mu_init ## assign the current value of mean equal to the initial point 
  S_cur<- sigma_init ## assing the current value of covariance equal to the initial point 
  #iter<- rep(0,maxit) 
  
  for (i in 1:maxit){ ## define a for loop to do update the estimates up to the maximum of 1000 times
    
    
    MU_cur<- M_cur ## define a for loop to do update the estimates up to the maximum of 1000 times
    SS_cur<- S_cur ## the current state of covariance calling inside the for loop
    
    
    MuNew<- MEM(Data=mydat,EXi=MEXI(Data=mydat,mu= MU_cur,sigma=SS_cur)) ## updating step for mean
    
    
    SigmaNew<- Mysig(Data=mydat,
                   EXX=EXXest(Data=mydat,
                              EXi=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur),
                              ss1=McovXI(Data=mydat,mu=MU_cur,sigma=SS_cur)),
                   EX=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur),
                   UPDMU=MEM(Data=mydat,EXi=MEXI(Data=mydat,mu=MU_cur,sigma=SS_cur)))  ## updating step for covariance matrix
    
    
    
    Mu_new<- MuNew ## reassign the new estimate of mean (renaming)
    Sigma_new<- SigmaNew ## reassign the new estimate of covariance (renaming)
    
    
    diff1<- MU_cur-Mu_new ## find the differences between the current mean vector and the new estimate of mean vector 
    diff2<- SS_cur-Sigma_new ## find the differences between the current covariance matrix and the new estimate of covariance matrix
    #print(diff1)
    #print(diff2)
    
    D1<- abs(mean(diff1)) ## calculate the absolute value of the mean of differences in means
    D2<- abs(mean(diff2)) ## calculate the absolute value of the mean of differences 
    ##in the elements of the covariance matrix
  
    if(D1<tol1 & D2<tol2){flag<- 1 ;break} ## define the condition for the absolute values of the differences according to the 
    ## stopping criteria, continue updating until the values of the parameter estimates meet both conditions, then stop and 
    ## take the final new estimates as the updated estimates of the parameters 
    
    M_cur<- Mu_new ## update the current value of mean vector with the new estimate
    
    S_cur<- Sigma_new ## update the current value of covariance matrix with the new estimate
    
    #iter[i]<- i
    
  }
  
  if(!flag) warning("Didn't Converge \n") ## the condition if there is no convergence in the stopping rule
  
  #list(paste("updated Mu of the EM=" , M_cur),paste("Updated Sigma of the EM=" ,S_cur),
  #    paste("iteration number=",i))
  #print(M_cur,S_cur)
  #print(S_cur)
  update=list(M_cur,(S_cur)) ## the list of the updated estimates of the parameters
  return(update) ## return the list of updated estimates of the parameters
  
  }
#############################################################################
##########################################################################
### Check if the function works on data ###

###INITIAL VALUES for dim=2 & dim =3####
## d=2##
mu1 <- c(67,67)
sigma1 <- matrix(c(3.1, 2.16, 2.16, 6.05), 2, 2)
## d=3##
#mu2<- c(67,67,67)
#sigma2<- matrix(c(3.1,2.2,1.6,2.2,4.3,2.45,1.6,2.45,4.8),3,3)
#sigma2
################################################################
####n=1000####
### TO produce the estimates of EM over 10 samples from simulated data
####JOHN: EM1000b10: means EM for n=1000 & binns=10 (it should be changed for other cases similarly)
EM1000b10<- matrix(rep(0*5*10),ncol=5,nrow=10) ## define and assign a matrix to put the resulted estimates of 
## of the parameters using EM algorithm for 10 simulated datasets
colnames(EM1000b10)<- c("mux1","mux2","Varx1","COVx1x2","COVx2x1","Varx2")## name the columns of the matrix
for(i in 1:nrow(EM1000b10)){ ## use for loop to find the estimates for all 10 simulated datasets 
  mydat<- simulateddata1000b10[,,i]## assign the new mydat (each of the 10 datasets) to mydat
  ######### JIHN: PLEASE NOTE THAT HERE, EACH TIME, simulateddata1000b10 SHOULD BE REPLACED BASED ON THE CASE
  
  EM1000b10[i,]<- unlist(EMBivG(Data=mydat,mu_init=mu1,sigma_init=sigma1,maxit=1000,tol1=1e-4,tol2=1e-5))## implement the EMBivG function
  ## and put the results on the relevant raw of the matrix
  
}

OUTEM1000b10<- write.csv(EM1000b10,"C:/Users/sh_za/Desktop/Results/Bivariate/EM/Sample size1000/EM1000b10.csv")

#colMeans(myoutsim50) 









