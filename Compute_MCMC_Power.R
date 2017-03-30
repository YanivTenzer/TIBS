# Compute power for the case where we know the effect size (i.e. joint distribution of X,Y)
compute.mcmc.permutation.power <- function(n, num_datasets, alpha,dependce_type, dependence_params,num_permutations,statistic_type,bias_method, bias_params)
{
  power <- 0; 
  for(i in 1:num_datasets) # Simulate many times (X,Y data set)
  {
    Data <- simulate.xy.data(n, dependce_type, dependence_params)
    BiasedData = Biased.xy.data(Data, n, bias_method, bias_params)
    
    d=dim(BiasedData)
    new_n=d[1]
    True_T<-compute.dependence.statistic(BiasedData, statistic_type,new_n)
    
    W=get.biased.sampling.weights(BiasedData, new_n, bias_method, bias_params)
    Permutations=MCMC_Permutations(W,num_permutations,new_n)
    
    #Compute statistic for each permutation:
    T = matrix(0,num_permutations,1);
    for(ctr in 1:num_permutations) 
    {
      T[ctr] <- compute.dependence.statistic(cbind(BiasedData[,1],BiasedData[Permutations[,ctr],2]), statistic_type,new_n)
    }
    
    pval <- sum(T > True_T)/num_permutations # compute empirical p-value
    
    #pval <- compute.mcmc.permutation.pval(Data, W, num_permutations) # compute empirical p-value
    if(pval < alpha)
    {
      power <- power+1; 
    }
  }
  
  power <- power / num_datasets;
  return (power); 
}