cat("\014")
path = 'D:/cond_ind_local'
setwd(path)
source('Simulate _Data.R')
source('runner_single_iteration.R')
library(gdata)
require(gdata)

#The size of the original sample (i.e., before truncation)
sample_size=500
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
#How many permutations/bootstrap samples to use:
num_of_repetitions_per_iterations = 200

#######Example 1: weighted-permutations test for Gaussian data with correlation parameter \rho=0.1
dependence_type='Gaussian'
switch(dependence_type,
       'Gaussian'={rho = 0
       prms<-list(rho)
       names(prms)<-c("rho")},
       'Exponential'={lambda=20
       prms<-list(lambda)
       names(marginal_prms)<-c('lambda')})       

test_type = 'permutations'
data<-simulate_data(sample_size, dependence_type, prms) 
truncated_data<-Create_Bias(data, bias_method, bias_params)
results<-runner_single_iteration(truncated_data, bias_method, bias_params,num_of_repetitions_per_iterations, 
                                  path, test_type,statistic_type)
Pvalue<-length(which(results$statistics_under_null>=results$True_T))/num_of_repetitions_per_iterations
cat("Pvalue = ", Pvalue, '\n')
######################################
#######Example 2: bootstrap test for Gumbel data:
dependence_type='Gumbel'
test_type = 'bootstrap'
data<-simulate_data(sample_size, dependence_type, NULL) 
truncated_data<-Create_Bias(data, bias_method, bias_params)
results<-runner_single_iteration(truncated_data, bias_method, bias_params,num_of_repetitions_per_iterations, 
                                 path, test_type,statistic_type)
Pvalue<-length(which(results$statistics_under_null>=results$True_T))/num_of_repetitions_per_iterations
cat("Pvalue = ", Pvalue, '\n')
#######################################
