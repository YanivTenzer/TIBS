# Parameters: 
# data - n*2 array of (x,y) samples
# bias_method - string indicating W 
# bias_params - parameters of W 
# TargetSampleSize - number of bootstrap/permutation samples to perform 
# input_dir - where to save results (?)
# test_type - test to perform 
# prms - additional parameters (needed?)
runner_single_iteration<-function(data, bias_method, bias_params, TargetSampleSize, input_dir, test_type, prms)
{  
  library(pracma)
  source('utilities.R')
  source('Tsai_test.R')
  #source('marginal_estimation.R')
  #################################################################
  #1. Compute weights matrix W:  
  W=get.biased.sampling.weights(data, dim(data)[1], bias_method, bias_params)
  #2.Create a grid of points, based on the data:
  permutation<-MCMC_Permutations(W, 1, dim(data)[1]) # why always sample a permutation using MCMC? 
  temp_data<-cbind(data[,1], data[permutation,2])
  
  #discard the extremum points to avoid numerical issues
  idx_minmax <-which( (temp_data[,1] %in% c(min(temp_data[,1],max(temp_data[,1])))) | 
                        (temp_data[,2] %in% c(min(temp_data[,2],max(temp_data[,2])))) )
  
  #  x<-as.vector(temp_data[-idx_minmax,1])
  #  y<-as.vector(temp_data[-idx_minmax,2])  
  grid_points<-temp_data[-idx_minmax,] # cbind(x,y)
  #  switch (bias_method, # keep only grid points with w>0? but MCMC sampling guarantees this anyway
  #          'Hyperplane_Truncation' = {signs<-bias_params$Hyperplane_prms[1]*x+bias_params$Hyperplane_prms[2]*y+bias_params$Hyperplane_prms[3]
  #          idx<-which(signs>0)
  #          grid_points<-grid_points[idx,]},
  #          
  #          'michas_example' = {abs_diff = apply(grid_points,1,function(x) abs(x[1]-x[2]))
  #          idx<-which(abs_diff<0.1)
  #          grid_points<-grid_points[idx,]},
  #          
  #          'huji'={signs<- ifelse(min(65-rowSums(grid_points), 18)>=0, 1, 0)})
  #####################################################################
  switch(test_type,
         'bootstrap'={
           fast_bootstrap <- 0 # 0 # New: allow not to reestimate nulls for each bootstrap sample 
           #3. Estimate the marginals
           marginals<-main_estimate_marginals(data, bias_method)
           pdfs<-marginals$PDFs  
           print("get null dist.")
           #4. Estimate W(x,y)*F_X*FY/normalizing_factor
           null_distribution<-get_null_distribution(data,pdfs,W)
           normalizing_factor<-null_distribution$normalizing_factor
           null_distribution<-null_distribution$null_distribution
           
           print("compute statistic")
           #1. First compute the statistics based on the original data set:
           True_T=Compute_Statistic(data, grid_points, null_distribution)
           #2. Compute statistic for dummy sample:
           statistics_under_null=matrix(0,TargetSampleSize,1)
           for(ctr in 1:TargetSampleSize) 
           {
             if(mod(ctr,10)==0)
               print(c("Run Boots=", ctr))
             dummy_sample<-Bootstrap(data, pdfs, bias_method, prms)
             
             if(!fast_bootstrap)
             {
               marginals_bootstrap<-main_estimate_marginals(dummy_sample, bias_method)  # Why are the marginals estimated each time? 
               pdfs_bootstrap<-marginals_bootstrap$PDFs
               
               print("Compute weights")
               #3. Compute weights matrix W:    
               W_bootstrap=get.biased.sampling.weights(dummy_sample, dim(dummy_sample)[1], bias_method, bias_params)
               
               print("Get Null distribution")
               #4. Estimate W(x,y)*F_X*FY/normalizing_factor
               null_distribution_bootstrap<-get_null_distribution(dummy_sample,pdfs_bootstrap,W_bootstrap)
               normalizing_factor_bootstrap<-null_distribution_bootstrap$normalizing_factor
               null_distribution_bootstrap<-null_distribution_bootstrap$null_distribution
             } else # just copy the distribution 
             {
               null_distribution_bootstrap<-null_distribution
             }
             print("compute statistic boots")
             statistics_under_null[ctr]<-Compute_Statistic(dummy_sample, grid_points, null_distribution_bootstrap)
           }
           output<-list(True_T=True_T,statistics_under_null=statistics_under_null)
         },
         ########################
         'permutations'={
           number_of_permutations = TargetSampleSize 
           print("Sample Perms")
           Permutations=MCMC_Permutations(W, number_of_permutations,dim(data)[1])
           
           expectations_table<-evaluate_mass_table_by_permutations(data, Permutations, grid_points)
           True_T=Compute_Statistic(data, grid_points, expectations_table)
           
           #Compute the statistics value for each permutation:
           statistics_under_null=matrix(0,TargetSampleSize,1)
           for(ctr in 1:TargetSampleSize) 
           {
             if(mod(ctr,100)==0)
               print(c("Comp. Stat. Perm=", ctr))
             statistics_under_null[ctr] <- Compute_Statistic(
               cbind(data[,1],data[Permutations[,ctr],2]), grid_points, expectations_table)
           }
           output<-list(True_T=True_T,statistics_under_null=statistics_under_null)
         },
         
         #Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
         'tsai' = { result<-tsai.test.ties(data[,1],data[,2])
         output<-list(Pvalue=rep(result[2], 1 + TargetSampleSize))
         },
         
         #MinP2 test, relevant only for truncation W(x,y)=1_{x<=y}
         'minP2' = {library(permDep)
           dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           B = TargetSampleSize
           results<-with(dat, permDep(trun, obs, B, delta, nc = 1, minp2Only = TRUE, sampling = 'conditional', kendallOnly = FALSE))
           
           output<-list(Pvalue=rep(results$p.valueMinp2, 1+TargetSampleSize))
         }
  )
  output$permuted_data <- temp_data # add example of permuted data 
  return(output)
}









