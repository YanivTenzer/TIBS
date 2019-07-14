# Perform Test for Independence under general Biased Sampling (TIBS)
# 
# Parameters: 
# data - n*2 array of (x,y) samples
# bias_method - string indicating biased sampling function W 
# bias_params - parameters of W 
# B - number of bootstrap/permutation samples to perform 
# test_type - test to perform 
# prms - additional parameters (needed for bootstrap)
TIBS <- function(data, bias_method, bias_params, B, test_type, prms)
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
  grid_points<-temp_data[-idx_minmax,] # cbind(x,y)
  
  switch(test_type,
         'bootstrap'={
           fast_bootstrap <- 1 # 0 # New: allow not to reestimate nulls for each bootstrap sample 
           #3. Estimate the marginals
           marginals<-main_estimate_marginals(data, bias_method)
           pdfs<-marginals$PDFs  
           #4. Estimate W(x,y)*F_X*FY/normalizing_factor
           null_distribution<-get_null_distribution(pdfs, W)
           normalizing_factor<-null_distribution$normalizing_factor
           null_distribution<-null_distribution$null_distribution
           
           #1. First compute the statistics based on the original data set:
           True_T=Compute_Statistic(data, grid_points, null_distribution)
           #2. Compute statistic for dummy sample:
           statistics_under_null=matrix(0,B,1)
           for(ctr in 1:B) 
           {
             if(mod(ctr,100)==0)
               print(c("Run Boots=", ctr))
             dummy_sample<-Bootstrap(data, pdfs, bias_method, prms)
             if(!fast_bootstrap)
             {
               marginals_bootstrap<-main_estimate_marginals(dummy_sample, bias_method)  # Why are the marginals estimated each time? 
               pdfs_bootstrap<-marginals_bootstrap$PDFs
               
               #3. Compute weights matrix W:    
               W_bootstrap=get.biased.sampling.weights(dummy_sample, dim(dummy_sample)[1], bias_method, bias_params)
               
               #4. Estimate W(x,y)*F_X*FY/normalizing_factor
               null_distribution_bootstrap<-get_null_distribution(pdfs_bootstrap, W_bootstrap)
               normalizing_factor_bootstrap<-null_distribution_bootstrap$normalizing_factor
               null_distribution_bootstrap<-null_distribution_bootstrap$null_distribution
             } else # just copy the distribution 
             {
               null_distribution_bootstrap<-null_distribution
             }
             statistics_under_null[ctr]<-Compute_Statistic(dummy_sample, grid_points, null_distribution_bootstrap)
           }
           output<-list(True_T=True_T,statistics_under_null=statistics_under_null)
         },
         ########################
         'permutations'={
           Permutations=MCMC_Permutations(W, B, dim(data)[1])
           print("Compute expectation table")
           expectations_table<-evaluate_mass_table_by_permutations(data, Permutations, grid_points)
           True_T=Compute_Statistic(data, grid_points, expectations_table)
           
           #Compute the statistics value for each permutation:
           statistics_under_null=matrix(0,B,1)
           for(ctr in 1:B) 
           {
             if(mod(ctr,100)==0)
               print(c("Comp. Stat. Perm=", ctr))
             statistics_under_null[ctr] <- Compute_Statistic(
               cbind(data[,1], data[Permutations[,ctr],2]), grid_points, expectations_table)
           }
           output<-list(True_T=True_T,statistics_under_null=statistics_under_null)
         },
         
         #Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
         'tsai' = { result<-tsai.test.ties(data[,1],data[,2])
         output<-list(Pvalue=result[2])
         },
         
         #MinP2 test, relevant only for truncation W(x,y)=1_{x<=y}
         'minP2' = {library(permDep)
           dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           results<-with(dat, permDep(trun, obs, B, delta, nc = 1, minp2Only = TRUE, sampling = 'conditional', kendallOnly = FALSE))
           
           output<-list(Pvalue=results$p.valueMinp2)
         }
  )
  output$permuted_data <- temp_data # add example of permuted data 
  if(!("Pvalue" %in% names(output))) # Compute empirical P-value
    output$Pvalue <- length(which(output$statistics_under_null>=output$True_T))/B
  return(output)
}

