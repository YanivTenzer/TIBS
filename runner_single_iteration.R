runner_single_iteration<-function(data, bias_method, bias_params,TargetSampleSize, input_dir, test_type,statistic_type)
{  
  library(pracma)
  source('Compute_Statistic.R')
  source('Hyperplane_Truncated_Bivariate_Gaussian.R')
  source('get_null_distribution.R')
  source('Bootstrap.R')
  source('Get_Dummy_Sample.R')
  source('Tsai_test.R')
  source('Yanivs_estimate_marginals.R')
  
  #################################################################
  #1. Compute weights matrix W:    
  W=get.biased.sampling.weights(data, dim(data)[1], bias_method, bias_params)
  #1.################Creating a fixed grid:#######################
  permutation<-MCMC_Permutations(data,W,1,dim(data)[1])
  temp_data<-cbind(data[,1], data[permutation,2])
  
  idx_x_min<-which(temp_data[,1]==min(temp_data[,1]))
  idx_x_max<-which(temp_data[,1]==max(temp_data[,1]))
  
  idx_y_min<-which(temp_data[,2]==min(temp_data[,2]))
  idx_y_max<-which(temp_data[,2]==max(temp_data[,2]))
  
  t<-list(temp_data[-c(idx_x_min,idx_x_max,idx_y_min,idx_y_max),1],temp_data[-c(idx_x_min,idx_x_max,idx_y_min,idx_y_max),2])
  names(t)<-c('X', 'Y')
  x<-as.vector(t$X)
  y<-as.vector(t$Y)  
  grid_points<-cbind(x,y)
  
  switch (bias_method,
          'Hyperplane_Truncation' = {signs<-bias_params$Hyperplane_prms[1]*x+bias_params$Hyperplane_prms[2]*y+bias_params$Hyperplane_prms[3]
          idx<-which(signs>0)
          grid_points<-grid_points[idx,]})
  #####################################################################
 switch(test_type,
  'bootstrap'={
  #3. Estimate the marginals
  marginals<-Yanivs_estimate_marginals(data, 1)
  cdfs<-marginals$CDFs
  pdfs<-marginals$PDFs

  #4. Estimate W(x,y)*F_X*FY/normalizing_factor
  null_distribution<-get_null_distribution(data,pdfs,W)
  normalizing_factor<-null_distribution$normalizing_factor
  null_distribution<-null_distribution$null_distribution
  
  #5. Compute the test statistics:
  if(statistic_type == 'HHG' ||statistic_type == 'kendall')
  {
    #1. First compute the statistics based on the original data set:
    True_T=Compute_Statistic(data, statistic_type, grid_points, data, null_distribution)
    
    #2. Compute statistic for dummy sample:
    statistics_under_null=matrix(0,TargetSampleSize,1)
    for(ctr in 1:TargetSampleSize) 
    {
      dummy_sample<-Bootstrap(data, null_distribution)
      tryCatch(
      {
        marginals_bootstrap<-Yanivs_estimate_marginals(dummy_sample, 1)
        cdfs_bootstrap<-marginals_bootstrap$CDFs
        pdfs_bootstrap<-marginals_bootstrap$PDFs
        
        #3. Compute weights matrix W:    
        W_bootstrap=get.biased.sampling.weights(dummy_sample, dim(dummy_sample)[1], bias_method, bias_params)

        #4. Estimate W(x,y)*F_X*FY/normalizing_factor
        null_distribution_bootstrap<-get_null_distribution(dummy_sample,pdfs_bootstrap,W_bootstrap)
        normalizing_factor_bootstrap<-null_distribution_bootstrap$normalizing_factor
        null_distribution_bootstrap<-null_distribution_bootstrap$null_distribution
        
        statistics_under_null[ctr]<-Compute_Statistic(dummy_sample, statistic_type, grid_points,dummy_sample, null_distribution_bootstrap)},
        error= function(e){statistics_under_null[ctr]=Inf})
    }
  }
  output<-list(True_T,statistics_under_null)
  names(output)<-c('True_T', 'statistics_under_null')
  },
  ########################
  'permutations'={
    source('MCMC_Permutations.R')  
    source('compute_statistic_based_permutaions.R')  
    source('Evaluate_mass_table_by_permutations.R')
   
    number_of_permutations =  TargetSampleSize 
    Permutations=MCMC_Permutations(data,W,number_of_permutations,dim(data)[1])
  
    expectations_table<-Evaluate_mass_table_by_permutations(data, Permutations, number_of_permutations, grid_points)
    
    True_T=compute_statistic_based_permutaions(data, statistic_type, dim(data)[1], bias_method, grid_points, expectations_table)
    #Compute the statistics value for each permutation:
    statistics_under_null=matrix(0,TargetSampleSize,1)
    statistics_under_null_cor=matrix(0,TargetSampleSize,1)
    for(ctr in 1:TargetSampleSize) 
    {
      statistics_under_null[ctr]<-compute_statistic_based_permutaions(cbind(data[,1],data[Permutations[,ctr],2]),statistic_type,dim(data)[1], bias_method, grid_points, expectations_table)
    }
    output<-list(True_T,statistics_under_null)
    names(output)<-c('True_T', 'statistics_under_null')})
    
  return(output)
}









