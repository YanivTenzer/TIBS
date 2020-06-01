#########################################################################
# Simulate data with biased-sampling weighting function W
# Parameters: 
# n - sample size 
# dependence.type - distribution (copula, normal, ..)
# w.fun - function W(x,y) to use
# prms - parameters of distribution
#
# Output: 
# data - an n*2 array with (x,y) values
#########################################################################
SimulateBiasedSample <- function(n, dependence.type, w.fun, prms)
{
  library('copula')
  
  # rejection sampling   
  data <- matrix(0, n, 2)
  if(!('keep.all' %in% names(prms)))
    prms$keep.all <- 0
  if(prms$keep.all) # keep also data-points that we through away in biased samples
  {
    all.data <- matrix(0, 2*n, 2)
    all.k <- 1
  }
  
  if(!('W.max' %in% names(prms)))
    prms$W.max <- 1.0 # temp. W.max should be input    
  k = 1
  while(k<=n)
  {
    switch(dependence.type, # First sample from Fxy
           'Gaussian' ={ library(mvtnorm)
             xy<-rmvnorm(1, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2))         
           },
           'LogNormal'={ library(mvtnorm)
             xy<-exp(rmvnorm(1, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2)))         
           },
           'LD'={ library(copula)  # y ~ Weibull, x ~ Exponential 
             GaussianCop<- normalCopula(param=prms$rho, dim = 2, dispstr = "ex") # if needed
             ranks<- rCopula(1, GaussianCop)
             xy <- rep(0, 2)
             xy[2]<-qweibull(ranks[,1], shape = 3, scale = 8.5, lower.tail = TRUE, log.p = FALSE)
             xy[1]<-qexp(ranks[,2], rate = 0.2)
           }, 
           'nonmonotone_nonexchangeable'={ library(copula)  # y ~ weibull, x ~ Gaussian copula 
             GaussianCop<- normalCopula(param=prms$rho, dim = 2, dispstr = "ex") # if needed
             ranks<- rCopula(1, GaussianCop)
             xy <- rep(0, 2)
             xy[2]<-qweibull(ranks[,1], shape = 0.5, scale = 2, lower.tail = TRUE, log.p = FALSE)
             xy[1]<-0.5 * (ranks[,2] * sample(c(-1,1), 1)+ 1)
             # need to set a copula for the dependency between x and y
           },
           'Gumbel'= { # here rho must be > 1 
             xy <- qnorm(rCopula(1, gumbelCopula(prms$rho)))
           }, 
           'Clayton'={ library('copula')
             xy <- qnorm(rCopula(1, claytonCopula(prms$rho)))
           }, 
           'CLmix'={ library('copula')
             xy <- rCopula(1, claytonCopula(0.5-rbinom(1, 1, 0.5)))
           }, 
           'strictly_positive'={ library('mvtnorm')  # w(x,y) = exp( -(|x|+|y|)/4 ) < 1 
             xy <- rmvnorm(1, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2))
           },
           'UniformStrip'={
             xy.abs.diff <- 2
             while(xy.abs.diff>prms$rho)
             {
               xy <- runif(2)
               xy.abs.diff <- abs(xy[2]-xy[1])
             }
           }
    ) # end switch 
    
    # Next decide if to keep point based on W
    if(w.fun %in% c('truncation'))  # w(x,y)=1_{x<y}
    {
      keep <- xy[1] <= xy[2]
    } else
    {
      # w(x,y)>0 , use rejection sampling 
      keep <- rbinom(1, 1, w_fun_eval(xy[1], xy[2], w.fun)/prms$W.max)
    }
    if(keep) 
    {
      data[k,] <- xy
      k <- k+1
    }
    if(prms$keep.all)
    {
      if(all.k <= dim(all.data)[1])
        all.data[all.k,] <- xy
      else
        all.data <- rbind(all.data, xy)
      all.k <- all.k+1
    }
  }    
#  return(data)  
  if(prms$keep.all)
  {
    return(list(data=data, all.data=all.data[1:(all.k-1),]))
  }
  else
    return(list(data=data))
}

#######################################################################
# A set of biased sampling functions to be used 
# Input: 
# x, y - data 
# w.fun - string indicating W type 
# 
# Output: 
# The values of w evaluated at the (x,y) array 
########################################################################
w_fun_eval <- function(x, y, w.fun) {
  if(typeof(w.fun)=="character") {
    r <- switch(w.fun, 
                'truncation'={x<y},
                'Hyperplane_Truncation'={(x<y)},
                'exp'= { exp((-abs(x)-abs(y))/4)},
                'exponent_minus_sum_abs'= { exp((-abs(x)-abs(y))/4)},
                'huji'={pmax(pmin(65-x-y,18),0)},  # changed length bias to 65 (from back to 66)
                'stritcly_positive'={exp((-abs(x)-abs(y))/4)}, # like exp? 
                'sum'={x+y},
                'naive'={1})
  } else {  # here w.fun is a function. Apply it to the array 
    r <- w.fun(x,y)
  }
  return (r)
}


#########################################################################
#  Compute the n*n matrix of sampling weights:
# Parameters: 
# data - n*2 matrix with (x,y) sample
# w.fun - biased sampling function W
#########################################################################
w_fun_to_mat <- function(data, w.fun)
{
  n <- dim(data)[1]  # get sample size 
  w.mat = matrix(0,n,n)
  for(i in 1:n)
    w.mat[i,] <- w_fun_eval(data[i,1], data[,2], w.fun)
  return (w.mat)
}
