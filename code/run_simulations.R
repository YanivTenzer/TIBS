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

run.flag <- 1 # 1: run simulations inside R. -1: run simulations from outside command line.  0: load simulations results from file if they're available
const.seed <- 1 # set constant seed 

# Vectors with different dependency settings 
dependence.type <- c('UniformStrip', 'Gaussian', 'Clayton', 'Gumbel', 
                     'LD', 'nonmonotone_nonexchangeable', 'CLmix', 'Gaussian',
                     'LogNormal', 'Gaussian') # The last one is log-normal 
w.fun <- c('truncation', 'truncation', 'truncation', 'truncation', 
               'truncation', 'truncation', 'truncation', 
           'exponent_minus_sum_abs', 'sum', 'truncation') # not good that we have only one simulation with positive W. Should add X+Y?
monotone.type <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE) # is monotone
exchange.type <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE) # is exchangeable
# sample.size <- c(500, 100, 100, 100, 100, 100, 100, 100) # set all sample.sizes to 100 
prms.rho <- list(0.3, seq(-0.9, 0.9, 0.1), 0.5, 1.6, c(0, 0.4),
                 seq(-0.9, 0.9, 0.1), 0.5, c(0), c(0), c(0))#seq(-0.9, 0.9, 0.1), 
                 #c(0)) # Parameters for each sampling type 

# test.type<- c('uniform_importance_sampling') # ,'bootstrap')#c( 'permutations','permutations_inverse_weighting',
test.type <- c('permutations', 'uniform_importance_sampling', 'permutations_inverse_weighting', 'uniform_importance_sampling_inverse_weighting', 
            'bootstrap', 'bootstrap_inverse_weighting') #c( 'permutations','permutations_inverse_weighting',
            #  #'uniform_importance_sampling',
            #  'uniform_importance_sampling_inverse_weighting',
            #  'bootstrap', 
            #  'bootstrap_inverse_weighting', 
            #  'min_P2', 'Tsai')
num.sim <- length(dependence.type)
num.tests <- length(test.type)

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
    prms = list(B=100, sample.size=num_of_observations, iterations=500, plot.flag=0, alpha=0.05, sequential.stopping=0, 
                use.cpp=1, keep.all=1, perturb.grid=1) # , sample.by.bootstrap=1) # set running parameters here ! 
    prms$w.max = 1
    if(run.flag != 1)
      prms.rho[[s]] = as.numeric(args[4]) # temp for loading from user 
    print(paste0("s=", s))
    print(paste0("rho=", prms.rho[[s]]))
    # Call function. # run simulations function 
    print(paste("n=", prms$sample.size))
    if(const.seed)
      prms$seed <- 92669484 # 4524553
    T.OUT <- simulate_and_test(dependence.type[s], prms.rho[[s]], w.fun[s], test.type, prms) # run all tests 
  }
} # end loop on dependency types


test.legend <- paste(test.type, as.character(T.OUT$test.power))
col.vec <- c("blue", "black", "green", "orange", "gray", "pink", "yellow")
i=1
plot(T.OUT$test.stat[1,i,]- rowMeans(T.OUT$test.null.stat[1,i,,]), col=col.vec[i], pch=20, main='differences')
for(i in 2:length(test.type))
  points(T.OUT$test.stat[1,i,]- rowMeans(T.OUT$test.null.stat[1,i,,]), col=col.vec[i], pch=20)
legend(0, 200, test.legend, lwd=c(2,2), col=col.vec[1:length(test.type)], y.intersp=0.8, cex=0.6)
# points(T.OUT$test.stat[1,,]- rowMedians(T.OUT$test.null.stat[1,1,,]), col="blue")


jpeg(paste0("../figs/check_valid_", dependence.type[s], "_",  w.fun[s],  "_perturb_grid_", prms$perturb.grid, ".jpg"), width = 400, height = 400)
plot(c(0, prms$iterations), c(0,1), col="red", type="l", 
     main=paste0("Tests pvals and power, n=", num_of_observations, ", alpha=", prms$alpha, " pert=", prms$perturb.grid))
valid.tests <- rep(0, num.tests)
for(i in 1:num.tests)
{
  if(max(T.OUT$test.pvalue[1,i,])> -1)
  {
    points(sort(T.OUT$test.pvalue[1,i,]), col=col.vec[i], pch=20)
    valid.tests[i] <- 1
  }
}
legend(0, 1, test.legend[which(valid.tests>0)], lwd=c(2,2), col=col.vec[which(valid.tests>0)], y.intersp=0.8, cex=0.6)
dev.off()


#library(matrixStats)
#plot(T.OUT$test.stat[1,,], rowMeans(T.OUT$test.null.stat[1,1,,]))
#plot(T.OUT$test.stat[1,,], rowMedians(T.OUT$test.null.stat[1,1,,]), col="blue")
#lines(c(15000, 50000), c(15000, 50000), col="red")

