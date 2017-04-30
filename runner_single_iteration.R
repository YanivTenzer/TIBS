# Run single iteration of ... 
# 
runner_single_iteration<-function(data, dependence_type, bias_method, bias_params, TargetSampleSize, input_dir, prms, test_type)
{  
  library(pracma)
  source('Compute_Statistic.R')
  source('Hyperplane_Truncated_Bivariate_Gaussian.R')
  source('estimate_marginals.R')
  source('get_null_distribution.R')
  source('Bootstrap_according_to_null.R')
  source('Get_Dummy_Sample.R')
  statistic_type = 'HHG'
  
  #1.################Creating a fixed grid:#######################
  switch(dependence_type,
         'Gaussian'={t<-meshgrid(seq(-2,2,0.1), seq(-2,2,0.1))},
         'Exponential'={t<-meshgrid(seq(0,1.5,0.05), seq(0,1.5,0.05))})
  
  x<-as.vector(t$X)
  y<-as.vector(t$Y)  
  grid_points<-cbind(x,y)
  
  switch (bias_method,
          'Hyperplane_Truncation' = {signs<-bias_params$Hyperplane_prms[1]*x+bias_params$Hyperplane_prms[2]*y+bias_params$Hyperplane_prms[3]
          idx<-which(signs>0)
          grid_points<-grid_points[idx,]})
  #################################################################
  #2. Compute weights matrix W:    
  W=get.biased.sampling.weights(data, dim(data)[1], bias_method, bias_params)
  #################################################################
  switch(test_type,
         'bootstrap'={
           #3. Estimate the marginals (note that at this stage we assume that we know the marginals, that are standard normal):
           bEstimate_Marginals = 1
           marginals<-estimate_marginals(data, dependence_type, bias_method, bias_params, bEstimate_Marginals, prms)
           pdfs<-marginals$PDF
           cdfs<-marginals$CDF
           
           #4. Estimate W(x,y)*F_X*FY/normalizing_factor
           null_distribution<-get_null_distribution(data,pdfs,W)
           normalizing_factor<-null_distribution$normalizing_factor
           null_distribution<-null_distribution$null_distribution
           
           #5. Compute the test statistics:
           if(statistic_type == 'HHG')
           {
             #1. First compute the statistics based on the original data set:
             True_T=Compute_Statistic(data, statistic_type, grid_points, data, null_distribution)
             #2. Compute statistic for dummy sample:
             statistics_bootstrap=matrix(0,TargetSampleSize,1)
             for(ctr in 1:TargetSampleSize) 
             {
               dummy_sample<-Bootstrap_according_to_null(data, null_distribution)
               tryCatch(
                 {
                   marginals_bootstrap<-estimate_marginals(dummy_sample, dependence_type, bias_method, bias_params, bEstimate_Marginals)
                   pdfs_bootstrap<-marginals_bootstrap$PDF
                   cdfs_bootstrap<-marginals_bootstrap$CDF
                   
                   #3. Compute weights matrix W:    
                   W_bootstrap=get.biased.sampling.weights(dummy_sample, dim(dummy_sample)[1], bias_method, bias_params)
                   
                   #4. Estimate W(x,y)*F_X*FY/normalizing_factor
                   null_distribution_bootstrap<-get_null_distribution(dummy_sample,pdfs_bootstrap,W_bootstrap)
                   normalizing_factor_bootstrap<-null_distribution_bootstrap$normalizing_factor
                   null_distribution_bootstrap<-null_distribution_bootstrap$null_distribution
                   
                   statistics_bootstrap[ctr]<-Compute_Statistic(dummy_sample, statistic_type, grid_points,dummy_sample, null_distribution_bootstrap)},
                 error= function(e){statistics_bootstrap[ctr]=Inf})
               browser()
             }
           }
           
           output<-list(True_T,statistics_bootstrap)
           names(output)<-c('True_T', 'statistics_bootstrap')
         },
         ########################
         'permutations'={
           source('MCMC_Permutations.R')  
           source('prepare_null_mass_table_given_grid.R')
           source('compute_statistic_based_permutaions.R')  
           
           bEstimate_Marginals <- 1; # 1 - estimate expected, 0 - take analytic expectation   
           
           number_of_permutations <- TargetSampleSize #   100 # why do we need a new variable here? 
           Permutations=MCMC_Permutations(W, number_of_permutations)
           if(file.exists('resources/permutations_test_mass_table.Rdata'))
           {
             load('resources/permutations_test_mass_table.Rdata')
           }else
           {
             mass_table<-prepare_null_mass_table_given_grid(bias_params, grid_points) # This assumes known distribution 
             dir.create('resources',showWarnings = FALSE)
             save(mass_table, file='resources/permutations_test_mass_table.Rdata')
           }
           
           
           if(bEstimate_Marginals) # New: compute expected from data !!! (assume unknown distribution)
           {
             perm_weights <- matrix(1, number_of_permutations, 1); # take uniform weights for sampled permutations. 
             Q <- compute_permutations_marginal_table(Permutations, perm_weights); # Compute marginal table
             mass_table <- compute_statistic_based_permutaions_quartiles(data, Q, grid_points) # Compute expected 
           }
           
           
           True_T=compute_statistic_based_permutaions(data, statistic_type, dim(data)[1], bias_method, dependence_type, grid_points, mass_table)
           statistics_permutations=matrix(0,TargetSampleSize,1)
           for(ctr in 1:TargetSampleSize) 
           {
             print(ctr)
             statistics_permutations[ctr]=compute_statistic_based_permutaions(cbind(data[,1],data[Permutations[,ctr],2]), statistic_type, dim(data)[1],bias_method,dependence_type, grid_points, mass_table)
           }
           output<-list(True_T,statistics_permutations)
           names(output)<-c('True_T', 'statistics_permutations')
         })
  return(output)
}









