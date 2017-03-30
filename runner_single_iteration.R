runner_single_iteration<-function(data, marginals_type, bias_method, bias_params,TargetSampleSize, input_dir, prms, test_type,statistic_type)
{  
  library(pracma)
  source('Compute_Statistic.R')
  source('Compute_Hoeffding_Statistic.R')
  source('Hyperplane_Truncated_Bivariate_Gaussian.R')
  source('estimate_marginals.R')
  source('get_null_distribution.R')
  source('Bootstrap_according_to_null.R')
  source('Get_Dummy_Sample.R')
  source('Tsai_test.R')
 
  #1. Compute weights matrix W:    
  W=get.biased.sampling.weights(data, dim(data)[1], bias_method, bias_params)
  bEstimate_Marginals = ifelse(identical(marginals_type, 'Unknown'),1,0)
  #2.################Creating a fixed grid:#######################
  switch(marginals_type,
         'Gaussian'={t<-meshgrid(seq(-2,2,0.1), seq(-2,2,0.1))
          x<-as.vector(t$X)
          y<-as.vector(t$Y)
          grid_points<-cbind(x,y)},
         'Exponential'={t<-meshgrid(seq(0,1.5,0.05), seq(0,1.5,0.05))},
         'Unknown'={
           grid_permutation<-MCMC_Permutations(data,W,1,dim(data)[1])
           grid_points<-cbind(data[,1], data[grid_permutation,2])
           
           idx_x_min<-which(grid_points[,1]==min(grid_points[,1]))
           idx_x_max<-which(grid_points[,1]==max(grid_points[,1]))      
           idx_y_min<-which(grid_points[,2]==min(grid_points[,2]))
           idx_y_max<-which(grid_points[,2]==max(grid_points[,2]))
           grid_points<-grid_points[-c(idx_x_min,idx_x_max, idx_y_min, idx_y_max),]
           x<-grid_points[,1]
           y<-grid_points[,2]})
  
    switch (bias_method,
          'Hyperplane_Truncation' = {signs<-bias_params$Hyperplane_prms[1]*x+bias_params$Hyperplane_prms[2]*y+bias_params$Hyperplane_prms[3]
          idx<-which(signs>0)
          grid_points<-grid_points[idx,]})
  #################################################################
 switch(test_type,
  'bootstrap'={
  #3. Estimate the marginals:
  marginals<-estimate_marginals(data, marginals_type, bias_method, bias_params, bEstimate_Marginals, prms)
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
        marginals_bootstrap<-estimate_marginals(dummy_sample, marginals_type, bias_method, bias_params, bEstimate_Marginals)
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
    source('Evaluate_mass_table_by_permutations.R')
  
    number_of_permutations =  TargetSampleSize 
    Permutations=MCMC_Permutations(data,W,number_of_permutations,dim(data)[1])
    switch(toString(bEstimate_Marginals), 
    '1'={expectations_table<-Evaluate_mass_table_by_permutations(data, Permutations, number_of_permutations, grid_points)},
    '0'={if(file.exists('resources/permutations_test_mass_table.Rdata')){
        load('resources/permutations_test_mass_table.Rdata')
       }else{
         expectations_table<-prepare_null_mass_table_given_grid(bias_params, grid_points)
          dir.create('resources',showWarnings = FALSE)
          save(expectations_table, file='resources/permutations_test_mass_table.Rdata')
        }})
    #Compute the "true" statistics value:
    True_T=compute_statistic_based_permutaions(data, statistic_type, dim(data)[1], bias_method, marginals_type, grid_points, expectations_table)
    #Compute the statistics value for each permutation:
    #Permutations=MCMC_Permutations(data,W,number_of_permutations,dim(data)[1])
    statistics_permutations=matrix(0,TargetSampleSize,1)
    for(ctr in 1:TargetSampleSize) 
    {
      statistics_permutations[ctr]=compute_statistic_based_permutaions(cbind(data[,1],data[Permutations[,ctr],2]), statistic_type, 
                                                                       dim(data)[1],bias_method,marginals_type, grid_points, expectations_table)
    }
   browser()  
   output<-list(True_T,statistics_permutations)
   names(output)<-c('True_T', 'statistics_permutations')})
    
  return(output)
}









