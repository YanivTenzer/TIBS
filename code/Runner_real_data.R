RUNNER_real_data<-function(data, Hyperplane_prms)
{  
  path = 'C:/Users/User/Dropbox/Independence_Test'
  setwd(paste(path, '/src', sep=""))
  
  source('GET_QUARTER_PROBABILITY.r')
  source('MCMC_Permutations.R')
  source('Compute_Power.R')
  source('Estimate_Marginals.R')
  source('Simulate_Data.R')
  source('Get_Quarter_Approximated_Mass.R')
  source('compute_statistic_utilizing_estimated_marginals.R')
  
  bUSE_FIXED_GRID = 0
  dependence_type='Gaussian'
  bUSE_CORRECTION = 1
  #dependence_params<-list(1,0.99)
  #names(dependence_params)<-c("noise", "rho")
  
  TargetSampleSize = 100
  maxIter = 50
  #Which tests to use:
  bUSE_PERMUTATIONS_TEST = 1 
  biased_method = 'Hyperplane_Truncation'  
  bCOMPUTE_POWER = 1;
  statistic_type = 'HHG'

  dir.create(file.path('D:/','Independence_Test'), showWarnings=FALSE )
  dir.create(file.path('D:/Independence_Test', 'Real_Data' ),showWarnings=FALSE)
  dir.create(file.path('D:/Independence_Test/Real_Data/', 'Data' ),showWarnings=FALSE)
  
  Root_Dir ='D:/Independence_Test/Real_Data' 
  DataDir=paste(Root_Dir, '/Data', sep="")
  ResultsDir = paste(Root_Dir, '/Results', sep="")
  
  dir.create(file.path(Root_Dir, 'Results' ),showWarnings=FALSE)
  
  True_T=matrix(0,maxIter,1)
  True_T_b=matrix(0,maxIter,1)
  Perm_pval=matrix(0,maxIter,1)
  Perm_pval_b=matrix(0,maxIter,1)
  bREJECT = matrix(0,maxIter,1)
  bREJECT_b = matrix(0,maxIter,1)
  ###The Main loop:
  for (iter in 1:maxIter)
  { 
      mu=c(0,0)
      rho = 0.99
      Sigma=matrix(c(1, rho, rho,1),2,2)
      data<-mvrnorm(n = 700, mu, Sigma)
      browser() 
      # data<-read.table('C:/Users/User/Dropbox/Independence_Test/data/simple data.txt')
      # data<-data[-1,1:2]
      # x<-suppressWarnings(as.numeric(data[,1]))
      # y<-suppressWarnings(as.numeric(data[,2]))
      # data<-cbind(x,y)
      # browser()
      
      Hyperplane_prms<-c(-1, 2,0.5)
      bias_params<-list(Hyperplane_prms)
      names(bias_params)<-c("Hyperplane_prms")
      
      BiasedData = Biased.xy.data(data,dim(data)[1],biased_method,bias_params)
      d=dim(BiasedData)
      n_biased=d[1]
      
      browser()
      marginals<-Estimate_Marginals(BiasedData, biased_method, bias_params)
      pdfs<-marginals$PDF
      cdfs<-marginals$CDF
      
      if(bUSE_PERMUTATIONS_TEST)
      {
        #Towards MCMC - get W matrix:
        W=get.biased.sampling.weights(BiasedData, n_biased, biased_method, bias_params) 
        #Sample permutations:
        Permutations=MCMC_Permutations(W, TargetSampleSize)
        #Now test for independency using the permutations/LRT tests:
        if(statistic_type == 'HHG')
        {
          True_T[iter]=compute_statistic_utilizing_estimated_marginals(BiasedData, statistic_type, n_biased, pdfs,
                                                                       biased_method, bias_params)
        }
        #Compute statistic for each permutation:
        T = matrix(0,TargetSampleSize,1);
        T_b = matrix(0,TargetSampleSize,1);
        for(ctr in 1:TargetSampleSize) 
        {
          
          T[ctr]=compute_statistic_utilizing_estimated_marginals(cbind(BiasedData[,1],BiasedData[Permutations[,ctr],2]), 
                                                                 statistic_type, n_biased,cbind(pdfs[,1],pdfs[Permutations[,ctr],2]),biased_method, bias_params)
        }
            
        Perm_pval[iter] <- length(which(T > True_T[iter]))/TargetSampleSize # compute  p-value (under the null hypothesis)
        Perm_pval_b[iter] <- length(which(T_b > True_T_b[iter]))/TargetSampleSize

        browser()
        if(bCOMPUTE_POWER)
        {
          #Should null hypothesis be rejected ?
          alpha=0.05

          SORTED_STATISTICS = sort(c(True_T[iter],T), index.return=TRUE)
          RANK=which(SORTED_STATISTICS$ix==1)
          if(TargetSampleSize-RANK<=floor(alpha*TargetSampleSize)-1)
          {
            bREJECT[iter]=1
          }
          
          # SORTED_STATISTICS = sort(c(True_T_b[iter],T_b), index.return=TRUE)
          # RANK=which(SORTED_STATISTICS$ix==1)
          # if(TargetSampleSize-RANK<=floor(alpha*TargetSampleSize)-1)
          # {
          #   bREJECT_b[iter]=1
          # }
        }
            
      }
                    
  }
              
  if(bUSE_PERMUTATIONS_TEST)
  {
    Power = mean(bREJECT)
    Pvalue = mean(Perm_pval)
    
    # Power_b = mean(bREJECT_b)
    # Pvalue_b = mean(Perm_pval_b)
  }
  
  browser()
  return(c(Power,Pvalue, Power_b, Pvalue_b))
}









