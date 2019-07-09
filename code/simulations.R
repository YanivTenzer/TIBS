cat("\014")
# path = 'D:/cond_ind_2019' # change path
path = 'C:\\Users\\Or Zuk\\Dropbox\\BiasedSampling\\Code'
setwd(path)
source('simulate_biased_sample.R')
source('runner_single_iteration.R')
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

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

output_dir = 'monotone_dependence_simulation'
dir.create(output_dir)
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
Estimate_Marginals = 1

prms = c()
output_dir = 'Gaussian_dependence_simulation'
dir.create(output_dir)
dependence_type = 'Gaussian'
bias_type = 'truncation'

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
  print(c("s=", s))
  set.seed(1)
  num_prms <- length(prms_rho[[s]])
  
  for(iter in 1:num_prms)
  {
    prms$rho = prms_rho[[s]][iter]
    prms$W_max <- 1 # temp: 1 for all weights
    
    results_table<- foreach(i=seq(iterations), .combine=rbind) %dopar%{ 
      iteration_result <- matrix(0, 4, B+1)
      truncated_data<-simulate_biased_sample(sample_size[s], dependence_type[s], bias_type[s], prms) 
      for(j in 1:4) # test_type in c('bootstrap', 'permutations', 'tsai', 'minP2'))
      {
        print(c("j=", j, test_type[j]))
        start_time <- Sys.time()
        test_results<-runner_single_iteration(truncated_data, bias_method, bias_params,   # change to bias_type[s]
                                                      B, path, test_type[j], prms)
        print(Sys.time() - start_time)
        if("Pvalue" %in% names(test_results))
        {
          iteration_result[j,1:(B+1)] <- test_results$Pvalue
        } else
        {
          iteration_result[j,1:B] <- test_results$statistics_under_null
          iteration_result[j,B+1] <- test_results$True_T
        }
      } 
      iteration_result
    }  # end parallel collection

    # Modify to save in ONE file: (need to change output_dir / file name to indicate dependence type)
    save(results_table, file=paste0(output_dir ,'/Gaussian_rho_',toString(prms$rho), '_four_tests_results', '.Rdata'))
    print(xtable(results_table, type = "latex"), file = paste0(output_dir ,'/Gaussian_rho_',toString(prms$rho), '_four_tests_results', '.tex')) # new! save in latex format 
  } # end simulation and testing for one parameter (iter)
} # end loop on dependency types 


