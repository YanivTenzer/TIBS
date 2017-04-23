cat("\014")
path = 'D:/cond_ind_local'
setwd(path)
source('Simulate _Data.R')
source('runner_single_iteration.R')
library(foreach)
library(doSNOW)
library(parallel)
library(doParallel)
library(gdata)
require(gdata)
#######################################################################################################################################
####Permutations with known Gaussian marginals:
cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

test_type = 'permutations'
statistic_type = 'HHG'

#The size of the original sample (i.e., before truncation)
sample_size=700
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
num_of_iterations = 200
#How many permutations to use/how many bootstrap samples to use:
num_of_repetitions_per_iterations = 300

rho =0.9
prms<-list(rho)
names(prms)<-c("rho")

results_table<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
  #Simulate bivariate Gaussian data:
  data<-simulate_data(sample_size, 'Gaussian', prms) 
  truncated_data<-Create_Bias(data, bias_method, bias_params)
  #Run a permutations test (assume marginals are known):
  results<-runner_single_iteration(truncated_data, 'Gaussian', bias_method, bias_params,num_of_repetitions_per_iterations, path, prms, test_type,statistic_type)
  switch(test_type,
         'bootstrap'={temp<-c(results$True_T, results$statistics_bootstrap)},
         'permutations'={temp<-c(results$True_T, results$statistics_permutations)})
    
    temp
  }
  
Pvalues<-apply(results_table, 1, function(x) length(which(x[-1]>x[1]))/length(x))
Power<-length(which(Pvalues<alpha))/num_of_iterations 
stopCluster(cl)
#################################################################################################################################
####Permutations with *unknown* (Gaussian) marginals:
test_type = 'permutations'
statistic_type = 'HHG'

#The size of the original sample (i.e., before truncation)
sample_size=700
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
num_of_iterations = 200
#How many permutations to use:
num_of_repetitions_per_iterations = 300

rho = 0
prms<-list(rho)
names(prms)<-c("rho")

cores=detectCores()
cl<-makeCluster(cores[1]-1) 
registerDoParallel(cl) 

results_table<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
  #Simulate bivariate Gaussian data:
  data<-simulate_data(sample_size, 'Gaussian', prms) 
  truncated_data<-Create_Bias(data, bias_method, bias_params)
  #Run the Permutations test (assume marginals are not known):
  results<-runner_single_iteration(truncated_data, 'Unknown', bias_method, bias_params,num_of_repetitions_per_iterations, path, prms, test_type,statistic_type)
  switch(test_type,
         'bootstrap'={temp<-c(results$True_T, results$statistics_bootstrap)},
         'permutations'={temp<-c(results$True_T, results$statistics_permutations)})
  
  temp
}

Pvalues<-apply(results_table, 1, function(x) length(which(x[-1]>x[1]))/length(x))
Power<-length(which(Pvalues<alpha))/num_of_iterations 
stopCluster(cl)
##################################################################################################################################
test_type = 'bootstrap'
statistic_type = 'HHG'

#The size of the original sample (i.e., before truncation)
sample_size=700
dependence_type='Gaussian'
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
num_of_iterations = 200
#How many permutations to use:
num_of_repetitions_per_iterations = 300

rho = 0
prms<-list(rho)
names(prms)<-c("rho")

cores=detectCores()
cl<-makeCluster(cores[1]-1) 
registerDoParallel(cl) 

results_table<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
  #Simulate bivariate Gaussian data:
  data<-simulate_data(sample_size, 'Gaussian', prms) 
  truncated_data<-Create_Bias(data, bias_method, bias_params)
  #Run a bootstrap test (assume marignals are known):
  results<-runner_single_iteration(truncated_data, 'Gaussian', bias_method, bias_params,num_of_repetitions_per_iterations, path, prms, test_type,statistic_type)
  switch(test_type,
         'bootstrap'={temp<-c(results$True_T, results$statistics_bootstrap)},
         'permutations'={temp<-c(results$True_T, results$statistics_permutations)})
  
  temp
}

Pvalues<-apply(results_table, 1, function(x) length(which(x[-1]>x[1]))/length(x))
Power<-length(which(Pvalues<alpha))/num_of_iterations 
stopCluster(cl)
##################################################################################################################################
test_type = 'bootstrap'
statistic_type = 'HHG'

#The size of the original sample (i.e., before truncation)
sample_size=700
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
num_of_iterations = 200
#How many permutations to use:
num_of_repetitions_per_iterations = 300

rho =0.01
prms<-list(rho)
names(prms)<-c("rho")

cores=detectCores()
cl<-makeCluster(cores[1]-1) 
registerDoParallel(cl) 

results_table<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
  #Simulate bivariate Gaussian data:
  data<-simulate_data(sample_size, 'Gaussian', prms) 
  truncated_data<-Create_Bias(data, bias_method, bias_params)
  #Run a bootstrap test (assume marignals are not known):
  results<-runner_single_iteration(truncated_data, 'Unknown', bias_method, bias_params,num_of_repetitions_per_iterations, path, prms, test_type,statistic_type)
  switch(test_type,
         'bootstrap'={temp<-c(results$True_T, results$statistics_bootstrap)},
         'permutations'={temp<-c(results$True_T, results$statistics_permutations)})
  
  temp
}

Pvalues<-apply(results_table, 1, function(x) length(which(x[-1]>x[1]))/length(x))
Power<-length(which(Pvalues<alpha))/num_of_iterations 
stopCluster(cl)
##################################################################################################################################
cores=detectCores()
cl<-makeCluster(cores[1]-1) 
registerDoParallel(cl) 

num_of_iterations = 200
sample_size=300
statistic_type = 'HHG'
test_type = 'permutations'
dependence_type='Sin'

results_table<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
  data<-simulate_data(sample_size, dependence_type, NULL) 
  truncated_data<-Create_Bias(data, bias_method, bias_params)
  results<-runner_single_iteration(truncated_data, 'Unknown', bias_method, bias_params,num_of_repetitions_per_iterations, 
                                   path, NULL, test_type,statistic_type)
  switch(test_type,
         'bootstrap'={temp<-c(results$True_T, results$statistics_bootstrap)},
         'permutations'={temp<-c(results$True_T, results$statistics_permutations)})
  
  temp }

Pvalues<-apply(results_table, 1, function(x) length(which(x[-1]>x[1]))/length(x))
outfile<-'Permutations_vs_Tsai_sin_example.Rdata'
save(Pvalues,results_table, file=outfile)
stopCluster(cl)