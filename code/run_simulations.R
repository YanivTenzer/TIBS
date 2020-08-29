rm(list=ls())
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # get path 
path = getwd()
args=commandArgs(trailingOnly = TRUE)
Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 

source('simulate_and_test.R')
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
source('import_samp.R')
library(foreach)
library(doSNOW)
#library(parallel)
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
const.seed <- 1 # set constant seed 

# Vectors with different dependency settings 
dependence.type <- c('UniformStrip', 'Gaussian', 'Clayton', 'Gumbel', 
                     'LD', 'nonmonotone_nonexchangeable', 'CLmix', 'Gaussian',
                     'strictly_positive', 'Gaussian') # The last one is log-normal 
w.fun <- c('truncation', 'truncation', 'truncation', 'truncation', 
               'truncation', 'truncation', 'truncation', 
           'exponent_minus_sum_abs', 'sum', 'const') # not good that we have only one simulation with positive W. Should add X+Y?
monotone.type <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE) # is monotone
exchange.type <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE) # is exchangeable
# sample.size <- c(500, 100, 100, 100, 100, 100, 100, 100) # set all sample.sizes to 100 
prms.rho <- list(0.3, seq(-0.9, 0.9, 0.1), 0.5, 1.6, c(0, 0.4),
                 seq(-0.9, 0.9, 0.1), 0.5, c(0), c(0), c(0))#seq(-0.9, 0.9, 0.1), 
                 #c(0)) # Parameters for each sampling type 

test.type<- c('permutations', 'uniform_importance_sampling', 'permutations_inverse_weighting', 'uniform_importance_sampling_inverse_weighting') # ,'bootstrap')#c( 'permutations','permutations_inverse_weighting',
            #  #'uniform_importance_sampling',
            #  'uniform_importance_sampling_inverse_weighting',
            #  'bootstrap', 
            #  'bootstrap_inverse_weighting', 
            #  'min_P2', 'Tsai')
num.sim <- length(dependence.type)
if(run.flag == 1)
{
  run.dep <- c(10)#(8:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
} else  # run from command line 
{
  run.dep <- as.integer(args[1]) #  c(7) # 2:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
  iterations = as.integer(args[2])
} # 4  # 10 for minP2 which is very slow  # 00 # 500  # Number of simulated dataset. Shared by all simulations

#test.type<-c( 'permutations','permutations_inverse_weighting',
#              #'uniform_importance_sampling',
#              'uniform_importance_sampling_inverse_weighting',
#              'bootstrap', 
#              'bootstrap_inverse_weighting')

if(isempty(intersect(run.dep, c(3,4,5,6,7)))) # %in% )
  library(copula) # needed for most simulations 

for(s in run.dep) # Run all on the farm  
{
  for(num_of_observations in c(100))#seq(250, 400, 50))
  {
    prms = list(B=100, sample.size=num_of_observations, iterations=100, plot.flag=0, alpha=0.1, sequential.stopping=0, 
                use.cpp=0, keep.all=0) # , sample.by.bootstrap=1) # set running parameters here ! 
    prms$w.max = 1
    if(run.flag != 1)
      prms.rho[[s]] = as.numeric(args[4]) # temp for loading from user 
    print(paste0("s=", s))
    print(paste0("rho=", prms.rho[[s]]))
    # Call function. # run simulations function 
    print(paste("n=", prms$sample.size))
    if(const.seed)
      prms$seed <- 38322
    T.OUT <- simulate_and_test(dependence.type[s], prms.rho[[s]], w.fun[s], test.type, prms) # run all tests 
  }
} # end loop on dependency types

col.vec <- c("blue", "black", "green", "orange")
plot(c(0, prms$iterations), c(0,1), col="red", type="l")
for(i in 1:length(test.type))
  points(sort(T.OUT$test.pvalue[1,i,]), col=col.vec[i], pch=20)
legend(0, 1, test.type, lwd=c(2,2), col=col.vec[1:length(test.type)], y.intersp=0.8, cex=0.6)


