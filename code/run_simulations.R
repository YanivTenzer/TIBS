path = 'C:/Users/Or Zuk/Documents/GitHub/TIBS/code'  # change to your path

setwd(path)
args=commandArgs(trailingOnly = TRUE)
# Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 

source('simulate_and_test.R')
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
source('import_samp.R')
library(foreach)
library(doSNOW)
library(parallel)
library(doParallel)
library(gdata)
library(mvtnorm)
library(ggplot2)  
library(pracma)

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

hyperplane.prms<-c(-1,1,0)
alpha <- 0.05 # significane threshold 
plot.flag <- 0 # plot and save figures 
run.flag <- 1 # 1: run simulations inside R. -1: run simulations from outside command line.  0: load simulations results from file if they're available
sequential.stopping <- 1 # New! use early stopping to save time ! 


# Vectors with different dependency settings 
dependence.type <- c('UniformStrip', 'Gaussian', 'Clayton', 'Gumbel', 
                     'LD', 'nonmonotone_nonexchangeable', 'CLmix', 'strictly_positive')
w.fun <- c('truncation', 'truncation', 'truncation', 'truncation', 
               'truncation', 'truncation', 'truncation', 'exponent_minus_sum_abs') # not good that we have only one simulation with positive W. Should add X+Y?
monotone.type <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) # is monotone
exchange.type <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) # is exchangeable
sample.size <- c(500, 100, 100, 100, 100, 100, 100, 100) # set all sample.sizes to 100 
prms.rho <- list(0.3, seq(-0.9, 0.9, 0.1), 0.5, 1.6, c(0, 0.4),
                 seq(-0.9, 0.9, 0.1), 0.5, seq(-0.9, 0.9, 0.1)) # Parameters for each sampling type 
test.type <- c('tsai', 'minP2', # other's tests 
               'permutations', 'bootstrap', 'fast-bootstrap', 'naive-bootstrap', 'naive-permutations') # different tests to run - new: add 
test.type <- c('permutations', 'bootstrap') #  'importance.sampling')  # for fast simulations . Add new importance sampling test 

if(run.flag == 1)
{
  run.dep <- c(5) # 2:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
  iterations = 4   # for minP2 which is very slow  # 00 # 500  # Number of simulated dataset. Shared by all simulations
  B = 150  # number of permtuations or bootstrap samples. Shared by all simulations 
  
} else  # run from command line 
{
  run.dep <- as.integer(args[1]) #  c(7) # 2:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
  iterations = as.integer(args[2])  # 4  # 10 for minP2 which is very slow  # 00 # 500  # Number of simulated dataset. Shared by all simulations
  B =  as.integer(args[3])  # 10^2  # number of permtuations or bootstrap samples. Shared by all simulations 
  
}
num.sim <- length(dependence.type)
if(isempty(intersect(run.dep,c(3,4,5,6,7)))) # %in% )
  library(copula) # needed for most simulations 

for(s in run.dep) # Run all on the farm  
{
  prms = list(B=100, sample.size=100, iterations=50, plot.flag=0, alpha=0.05, sequential.stopping=0, use.cpp=0) # NEW! Use RCPP 
  
  if(run.flag != 1)
    prms.rho[[s]] = as.numeric(args[4]) # temp for loading from user 
  print(paste0("s=", s))
  print(paste0("rho=", prms.rho[[s]]))
  set.seed(1) # set for randomization 
  # Call function. # run simulations function 
  print(paste("n=", sample.size[s]))
  T.OUT <- simulate_and_test(dependence.type[s], prms.rho[[s]], w.fun[s], test.type, prms) # run all tests 

} # end loop on dependency types


