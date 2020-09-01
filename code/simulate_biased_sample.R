#########################################################################
# Simulate data with biased-sampling weighting function W
# Parameters: 
# input.sample (optional) - fixed large sample from F_XY
# n - sample size 
# dependence.type - distribution (copula, normal, ..). Currently only string - will allow a function as input 
# w.fun - function W(x,y) to use
# prms - parameters of distribution
#
# Output: 
# data - an n*2 array with (x,y) values
#########################################################################
SimulateBiasedSample <- function(n, dependence.type, w.fun, prms, input.sample)
{
  library('copula')
  
  # rejection sampling   
  data <- matrix(0, n, 2)
  if(!('keep.all' %in% names(prms)))
    prms$keep.all <- 0
  if(prms$keep.all) # keep also data-points that we through away in biased samples
  {
    all.data <- matrix(0, 2*n, 2)
    all.k <- 0
  }
  
  
  if(!('w.max' %in% names(prms)))  # set w.max 
  {
    prms$w.max <- set_w_max(2*n, dependence.type, w.fun, prms)
    print(paste0("Setting w.max=", round(prms$w.max, 4)))    
  }
  
  
  if(!('sample.by.bootstrap' %in% names(prms)))  # new sampling method
    prms$sample.by.bootstrap = 0
  if(prms$sample.by.bootstrap) # here use the method for drawing many samples from F_XY, treat the empirical distribution as F_XY and sample with wegihts w from it
  {
    if(!exists('input.sample') | (is.na(input.sample)))
      input.sample <- SimulateSample(10^6, dependence.type, prms)  # draw a lot 
    w.vec <- w_fun_eval(input.sample[,1], input.sample[,2], w.fun)
    idx <- sample(dim(input.sample)[1], n, replace = TRUE, w.vec / sum(w.vec))   # sample with probabilities as normalized weights 
    data <- input.sample[idx,]
    all.data <- input.sample # if we want to return all
  } else # do rejection sampling
  {
    k = 0
    while(k<n)  # sample one by one (should change to sampling a vector)
    {
      if(is.function(dependence.type))  # new: sample from a general density function
        xy <- dependence.type(n-k)  # dependence.type returns a sample vector of length (n-k)*2
      else  
        xy <- SimulateSample(n-k, dependence.type, prms) # new! sample vector !
      
      # Next decide if to keep point based on W
      if(w.fun %in% c('truncation'))  # w(x,y)=1_{x<y}
        keep <- which(xy[,1] <= xy[,2])
      else       # w(x,y)>0 , use rejection sampling 
        keep <- which(rand(n-k,1) < w_fun_eval(xy[,1], xy[,2], w.fun)/prms$w.max)  # rbinom(1, 1, w_fun_eval(xy[1], xy[2], w.fun)/prms$w.max)
      n.keep <- length(keep)
      if(prms$keep.all)
      {
        if(all.k + n-k <= dim(all.data)[1])
          all.data[(all.k+1):(all.k+n-k),] <- xy
        else
          all.data <- rbind(all.data[1:all.k,], xy)
        all.k <- all.k+n-k
      }
      if(n.keep > 0) 
      {
        data[(k+1):(k+n.keep),] <- xy[keep,]
        k <- k+n.keep
      }
    }  # if bootstrap/rejection sampling
#    if(prms$keep.all)
#      print(paste0("Sampling Ratio = ", all.k / n))
    #    print(k)
  } # end while   k <= n
  #  print("return")
  if(prms$keep.all)
    return(list(data=data, all.data=all.data[1:all.k,], w.max=prms$w.max))
  else
    return(list(data=data, w.max=prms$w.max))
}



# Simulate sample without bias W 
SimulateSample <- function(n, dependence.type, prms)
{
  switch(dependence.type, # First sample from Fxy
         'Gaussian' ={ library(mvtnorm)
           xy.mat <- rmvnorm(n, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2))    
#           if(prms$rho == 0)
#           {
#             print("Randomize 2!!!")
#             xy.mat[,1] <- rnorm(n)
#             xy.mat[,2] <- rnorm(n)
#           }
             
         },
         'LogNormal'={ library(mvtnorm)
           xy.mat <- exp(rmvnorm(n, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2)))         
         },
         'strictly_positive'={ library(mvtnorm)  # another name for LogNormal  
           xy.mat <- exp(rmvnorm(n, c(0,0), matrix(c(1, prms$rho, prms$rho,1),2,2)))         
         },
         # Copula distributions  from here          
         'LD'={ library(copula)  # y ~ Weibull, x ~ Exponential 
           GaussianCop <- normalCopula(param=prms$rho, dim = 2, dispstr = "ex") # if needed
           ranks<- rCopula(n, GaussianCop)
           xy.mat <- cbind(qweibull(ranks[,1], shape = 3, scale = 8.5, lower.tail = TRUE, log.p = FALSE), 
                           qexp(ranks[,2], rate = 0.2))
         }, 
         'nonmonotone_nonexchangeable'={ library(copula)  # y ~ weibull, x ~ Gaussian copula 
           GaussianCop<- normalCopula(param=prms$rho, dim = 2, dispstr = "ex") # if needed
           ranks<- rCopula(n, GaussianCop)
           xy.mat <- cbind(qweibull(ranks[,1], shape = 0.5, scale = 2, lower.tail = TRUE, log.p = FALSE), 
                           0.5 * (ranks[,2] * sample(c(-1,1), 1)+ 1))
         },
         'Clayton'={ library('copula')
           xy <- qnorm(rCopula(n, claytonCopula(prms$rho)))
         }, 
         'Gumbel'= { # here rho must be > 1 
           xy.mat <- qnorm(rCopula(n, gumbelCopula(prms$rho)))
         }, 
  )
  if(!exists('xy.mat')) # for these distributions we simulate one by one
  {
    xy.mat <- matrix(0, n, 2)
    for(i in c(1:n))    
    {
      switch(dependence.type, # First sample from Fxy
             'CLmix'={ library('copula') # choose from mixture for each point independently 
               xy <- rCopula(n, claytonCopula(0.5-rbinom(1, 1, 0.5)))
             }, 
             
             'UniformStrip'={ # use rejection sampling 
               xy.abs.diff <- 2
               while(xy.abs.diff>prms$rho)
               {
                 xy <- runif(2)
                 xy.abs.diff <- abs(xy[2]-xy[1])
               }
             }
      ) # end switch dependence.type
      xy.mat[i,] <- xy 
    } # end loop on i  
  } # end if one-by-one 
  return(xy.mat)
}

  #######################################################################
# A set of biased sampling functions to be used 
# Input: 
# x, y - data 
# w.fun - string indicating W type . We assume this is a vectorized function !!! 
# 
# Output: 
# The values of w evaluated at the (x,y) array 
########################################################################
w_fun_eval <- function(x, y, w.fun) {
  if(typeof(w.fun)=="character") {
    r <- switch(w.fun, 
                'truncation'={as.numeric(x<y)},
                'Hyperplane_Truncation'={as.numeric(x<y)},
                'gaussian'= {exp((-x**2-y**2)/2)},
                'exp'= {exp((-abs(x)-abs(y))/4)},
                'exponent_minus_sum_abs'= { exp((-abs(x)-abs(y))/4)},
                'huji'={pmax(pmin(65-x-y,18),0)},  # changed length bias to 65 (from back to 66)
                'stritcly_positive'={exp((-abs(x)-abs(y))/4)}, # like exp? 
                'step'={as.numeric(x<y)+0.1},
                'sum'={x+y},
                'naive'={rep(1, length(x))}, 
                'const'={rep(1, length(x))})
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


#######################################################################
# New: return a function A set of biased sampling functions to be used 
# Input: 
# w.fun - string indicating W type 
# 
# Output: 
# w.fun - a real nonnegative function of two variables  
########################################################################
w_str_to_fun <- function(w.str)
{
  if(is.function(w.str))
    return(w.str)
  w.fun=function(x,y){
    return(w_fun_eval(x, y, w.str))
  }
  return(w.fun)
}


# Determine empirically  w_max for functions where we don't know it (this is approximate and may fail)
set_w_max <- function(n=1000, dependence.type, w.fun, prms)
{
  xy <- SimulateSample(n, dependence.type, prms)
  return(5 * (max(w_fun_eval(xy[,1], xy[,2], w.fun)) + 5 * std(w_fun_eval(xy[,1], xy[,2], w.fun)) + 1.0))     
} 

# Set maximum value for a specific sample 
set_w_max_sample <- function(data, w.fun)
{
    return(max(w_fun_to_mat(data, w.fun)))
}


# Determine if weighing function is positive 
is_pos_w <- function(w.fun, data, mat.flag)
{
  if(!is.function(w.fun))
  {
    if(w.fun %in% c('sum', 'sum_coordinates', 'exponent_minus_sum_abs', 'const', 'naive'))
      return(TRUE)
    if(w.fun %in% c('truncation', 'Hyperplane_Truncation'))
      return(FALSE)
  } 
  if(!mat.flag)  # run and compute values 
      return(min(w_fun_eval(data[,1], data[,2], w.fun)) > 0) # test only on sample points x_i, y_i
    else
      return(min(w_fun_to_mat(data, w.fun)) > 0) # test all pairs x_i, y_j
}




