# Simulate data with truncation X<Y
# Parameters: 
# n - sample size 
# dependence_type - distribution (copula, normal, ..)
# prms - parameters of distribution
#
# Output: 
# data - an n*2 array with (x,y) values
# 
simulate_biased_sample <- function(n, dependence_type, bias_type, prms)
{
  # rejection sampling   
  data <- matrix(0, n, 2)
  prms$W_max <- 1.0 # temp. Should be input    
  k = 1
  while(k<=n)
  {
    switch(dependence_type, # First sample from F_XY
           'Gaussian' ={ library(mvtnorm)
             xy<-rmvnorm(1, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2))         
           },
           'LogNormal'={ library(mvtnorm)
             xy<-exp(rmvnorm(1, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2)))         
           },
           'Gumbel'= { # here rho must be > 1 
             gumbel.cop <- gumbelCopula(prms$rho)#(1.2)
             ranks<- rCopula(1, gumbel.cop)
             xy = c(qnorm(ranks[1]), qnorm(ranks[2]))
           }, 
           'Clayton'={ library('copula')
             clayton.cop <- claytonCopula(prms$rho)
             u<- rCopula(n, clayton.cop)
             ranks<- rCopula(1, clayton.cop)
             xy = c(qnorm(ranks[1]), qnorm(ranks[2]))
           }, 
           'LogNormal'={ library(mvtnorm)
             xy<-exp(rmvnorm(1, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2)))         
           },
           'UniformStrip'={
             while((xy<-abs(runif(2)))>pems$rho)
             {
             }
           }
    ) # end switch 
    
    switch(bias_type, # Next decide if to keep point based on W
           'truncation' = { # w(x,y)=1_{x<y}
             keep <- xy[1] <= xy[2]
           },
           'strictly_positive' = { # w(x,y)>0 , use rejection sampling 
             keep <- rbinom(1, 1, biased.sampling.w(xy[1], xy[2], 'bias_type')/prms$W_max)
           })
    
    if(keep) 
    {
      data[k,] <- xy
      k <- k+1
    }
  }    
  return (data)  
}


#######################################################################
# Compute bias function used 
#
########################################################################
biased.sampling.w <- function(x, y, bias_type)
{
  r <- switch(bias_type, 
              'truncation'={(x<y)},
              'Hyperplane_Truncation'={(x<y)},
              'exp'= { exp((-abs(x)-abs(y))/4)},
              'huji'={pmax(pmin(66-x-y,18),0)},  # changed length bias to 66 (from back to 65)
              'sum'={x+y}
  )
  return(r)
}



##########################################################################
#  Compute the N*N matrix of sampling weights:
# Parameters: 
# data - n*2 matrix with (x,y) sample
# N - sample size 
# biased_method - method for matrix computation
# bias_params - parameters
#
#########################################################################
get.biased.sampling.weights <- function(data, N, biased_method, bias_params)
{
  W = matrix(0,N,N)
  
  # New: just call function w for each row 
  for(i in 1:N)
  {
    W[i,] <- biased.sampling.w(data[i,1], data[,2], biased_method)
  }
  return (W); 
}
