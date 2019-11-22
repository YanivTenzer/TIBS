# path = 'C:\\Users\\Or Zuk\\Dropbox\\BiasedSampling\\Code'  # change to your path
# setwd(path)
args=commandArgs(trailingOnly = TRUE)

source('simulate_and_test.R')
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
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
run.flag <- 1 # 1: run simulations again. 0: load simulations results from file if they're available

# Vectors with different dependency settings 
dependence.type <- c('UniformStrip', 'Gaussian', 'Clayton', 'Gumbel', 
                     'LD', 'nonmonotone_nonexchangeable', 'CLmix', 'strictly_positive')
bias.type <- c('truncation', 'truncation', 'truncation', 'truncation', 
               'truncation', 'truncation', 'truncation', 'exponent_minus_sum_abs') # not good that we have only one simulation with positive W. Should add X+Y?
monotone.type <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) # is monotone
exchange.type <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) # is exchangeable
sample.size <- c(500, 100, 100, 100, 100, 100, 100, 100) # set all sample.sizes to 100 
prms.rho <- list(0.3, seq(-0.9, 0.9, 0.1), 0.5, 1.6, c(0, 0.4),
                 seq(-0.9, 0.9, 0.1), 0.5, seq(-0.9, 0.9, 0.1)) # Parameters for each sampling type 
test.type <- c('tsai', 'minP2', # other's tests 
               'permutations', 'bootstrap', 'fast-bootstrap', 'naive-bootstrap', 'naive-permutations') # different tests to run - new: add 
run.dep <- as.integer(args[1]) #  c(7) # 2:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
iterations = as.integer(args[2])  # 4  # 10 for minP2 which is very slow  # 00 # 500  # Number of simulated dataset. Shared by all simulations
B =  as.integer(args[3])  # 10^2  # number of permtuations or bootstrap samples. Shared by all simulations 
num.sim <- length(dependence.type)
if(isempty(intersect(run.dep,c(3,4,5,6,7)))) # %in% )
  library(copula) # needed for most simulations 

#print(paste0("sample-size=", sample.size))
for(s in run.dep) # Run all on the farm  
{
  prms.rho[[s]] = as.integer(args[4]) # temp for loading from user 
  print(paste0("s=", s))
  print(paste0("rho=", prms.rho[[s]]))
  set.seed(1) # set for randomization 
  # Call function. # run simulations function 
  print(paste("n=", sample.size[s]))
  T.OUT <- simulate_and_test(dependence.type[s], prms.rho[[s]], bias.type[s], test.type, # run all tests 
                             B, sample.size[s], iterations, plot.flag)  
  
} # end loop on dependency types

