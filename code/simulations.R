cat("\014")
# path = 'D:/cond_ind_2019' # change path
path = 'C:\\Users\\Or Zuk\\Dropbox\\BiasedSampling\\Code'
setwd(path)
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
library(foreach)
library(doSNOW)
library(parallel)
library(doParallel)
library(gdata)
library(copula)
library(mvtnorm)
library(xtable)
library(Matrix)
library(ggplot2) # need to install it 

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
Estimate_Marginals = 1
alpha <- 0.05 # significane threshold 


prms = c()

# Vectors with different dependency settings 
dependence_type <- c('UniformStrip', 'Gaussian', 'Clayton', 'Gumbel', 'LD', 'nonmonotone_nonexchangeable', 'CLmix', 'strictly_positive')
bias_type <- c('truncation', 'truncation', 'truncation', 'truncation', 'truncation', 'truncation', 'truncation', 'exponent_minus_sum_abs')
monotone_type <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) # is monotone
exchange_type <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) # is exchangeable
sample_size <- c(500, 100, 100, 100, 500, 500, 100, 500) # sample_size=100
prms_rho <- list(0.3, seq(-0.9, 0.9, 0.1), 0.5, 1.6, c(0, 0.4),
                 seq(-0.9, 0.9, 0.1), NULL, seq(-0.9, 0.9, 0.1)) # Parameters for each sampling type 

test_type <- c('bootstrap', 'permutations', 'tsai', 'minP2') # different tests to run 
iterations = 20 # 500  # shared by all simulations
B = 10^2 # 10^4 # number of permtuations or bootstrap samples. Shared by all simulations 
num_sim <- length(dependence_type)

for(s in 1:num_sim)  # Loop on different dependency types 
{
  print(paste0("s=", s))
  set.seed(1)
  num_prms <- length(prms_rho[[s]])
  test_Pvalue <- test_time <- array(-1, c(num_prms, 4, iterations)) 
  for(iter in 1:num_prms) # loop on different simulation parameters 
  {
    prms$rho = prms_rho[[s]][iter]
    prms$W_max <- 1 # temp: 1 for all weights
    
    ## Parallel on multiple cores 
    ##    results_table<- foreach(i=seq(iterations), .combine=rbind) %dopar%{ 
    ##      iteration_result <- matrix(0, 4, B+1)
    for(i in 1:iterations){
      biased_data<-simulate_biased_sample(sample_size[s], dependence_type[s], bias_type[s], prms) 
      for(j in 1:3) # loop on statistical tests . Last test (4) minP2 is very slow 
      {
        # Skip irrelevant tests 
        if((!(bias_type[s] %in% c('truncation', 'Hyperplane_Truncation'))) & (test_type[j] %in% c("tsai", 'minP2')))
          next  # these tests run only for truncation 
        if((test_type[j] == 'bootstrap') & (bias_type[s] %in% c('truncation', 'Hyperplane_Truncation', 'huji')))
          next # can't run bootstrap because w can be zero 
        
        print(paste0(dependence_type[s], ', rho=', prms$rho, ' test: ', test_type[j], ' iter=', i, ' of ', iterations))
        start_time <- Sys.time()
        test_results<-TIBS(biased_data, bias_method, bias_params,   # change to bias_type[s]
                           B, test_type[j], prms)
        test_time[iter, j, i] <- Sys.time() - start_time
        test_Pvalue[iter, j, i] <- test_results$Pvalue
        print(Sys.time() - start_time)
      } # end loop on tests 
    }  # end parallel collection. Simulation and testing for one parameter (iter)
    
  } # end simulation and testing for one dependency type (loop over iter)
  
  # Compute power (this is also type-1-error alpha under the null)
  test_Power <- apply(test_Pvalue<alpha, 1, rowMeans) 
  # Modify to save in ONE file: (need to change output_dir / file name to indicate dependence type)
  save(test_Pvalue, test_time, test_Power, prms_rho, sample_size, s, file=paste0('results/', dependence_type[s], '_four_tests_results', '.Rdata'))
#  print(xtable(results_table, type = "latex"), file = paste0(output_dir ,'/Gaussian_rho_',toString(prms$rho), '_four_tests_results', '.tex')) # new! save in latex format 
} # end loop on dependency types 


