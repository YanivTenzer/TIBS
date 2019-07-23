# Perform Test for Independence under general Biased Sampling (TIBS)
# 
# Parameters: 
# data - n*2 array of (x,y) samples
# bias.method - string indicating biased sampling function W 
# bias.params - parameters of W 
# B - number of bootstrap/permutation samples to perform 
# test.type - test to perform 
# prms - additional parameters (needed for bootstrap)
#
TIBS <- function(data, bias.method, bias.params, B, test.type, prms)
{  
  library(pracma)
  source('utilities.R')
  source('Tsai_test.R')
  #################################################################
  #1. Compute weights matrix W:  
  W=GetBiasedSamplingWeights(data, dim(data)[1], bias.method, bias.params)
  #2.Create a grid of points, based on the data:
  permutation<-PermutationsMCMC(W, 1, dim(data)[1]) # why always sample a permutation using MCMC? 
  temp.data<-cbind(data[,1], data[permutation,2])
  #discard the extremum points to avoid numerical issues
  idx.minmax <-which( (temp.data[,1] %in% c(min(temp.data[,1],max(temp.data[,1])))) | 
                        (temp.data[,2] %in% c(min(temp.data[,2],max(temp.data[,2])))) )
  grid.points<-temp.data[-idx.minmax,]
  
  switch(test.type,
         'bootstrap'={
           fast.bootstrap <- 1 # 0 # New: allow not to reestimate nulls for each bootstrap sample 
           #3. Estimate the marginals
#           print("estimate marginals")
           marginals<-EstimateMarginals(data, bias.method)
           pdfs<-marginals$PDFs  
           #4. Estimate W(x,y)*Fx*Fy/normalizing.factor
#           print("Get null")
           null.distribution<-GetNullDistribution(pdfs, W)
           normalizing.factor<-null.distribution$normalizing.factor
           null.distribution<-null.distribution$null.distribution
           
           #1. First compute the statistics based on the original data set:
#           print("Comp. stat")
           TrueT=ComputeStatistic(data, grid.points, null.distribution)
           #2. Compute statistic for dummy sample:
           statistics.under.null=matrix(0,B,1)
           null.distribution.bootstrap<-null.distribution
           for(ctr in 1:B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Run Boots=", ctr))
             dummy.sample<-Bootstrap(data, pdfs, bias.method, prms)
             if(!fast.bootstrap)
             {
               marginals.bootstrap<-EstimateMarginals(dummy.sample, bias.method)  # Why are the marginals estimated each time? 
               pdfs_bootstrap<-marginals.bootstrap$PDFs
               
               #3. Compute weights matrix W:    
               W_bootstrap=GetBiasedSamplingWeights(dummy.sample, dim(dummy.sample)[1], bias.method, bias.params)
               
               #4. Estimate W(x,y)*Fx*FY/normalizing.factor
               null.distribution.bootstrap<-GetNullDistribution(pdfs_bootstrap, W_bootstrap)
               normalizing.factor_bootstrap<-null.distribution.bootstrap$normalizing.factor
               null.distribution.bootstrap<-null.distribution.bootstrap$null.distribution
             } # else # just copy the distribution 
             #             {
             #               null.distribution.bootstrap<-null.distribution
             #             }
             # MOST TIME SPENT HERE (STAT COMPUTE!) TEMP for profiling
             statistics.under.null[ctr]<-ComputeStatistic(dummy.sample, grid.points, null.distribution.bootstrap)
           }
           output<-list(TrueT=TrueT,statistics.under.null=statistics.under.null)
         },
         
         'permutations'={
           Permutations=PermutationsMCMC(W, B, dim(data)[1])
           expectations.table<-QuarterProbFromPermutations(data, Permutations, grid.points)
           TrueT=ComputeStatistic(data, grid.points, expectations.table)
           
           #Compute the statistics value for each permutation:
           statistics.under.null=matrix(0,B,1)
           for(ctr in 1:B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Comp. Stat. Perm=", ctr))
             statistics.under.null[ctr] <- ComputeStatistic(
               cbind(data[,1], data[Permutations[,ctr],2]), grid.points, expectations.table)
           }
           output<-list(TrueT=TrueT,statistics.under.null=statistics.under.null)
         },
         
         'tsai' = { result<-TsaiTestTies(data[,1],data[,2]) #Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
         output<-list(Pvalue=result[2])
         },
         
         'minP2' = {library(permDep)  #MinP2 test, relevant only for truncation W(x,y)=1_{x<=y}
           dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           results<-with(dat, permDep(trun, obs, B, delta, nc = 1, minp2Only = TRUE, 
                                      sampling = 'conditional', kendallOnly = FALSE))
           output<-list(Pvalue=results$p.valueMinp2)
         }
  )
  output$permuted.data <- temp.data # add example of permuted data 
  if(!("Pvalue" %in% names(output))) # Compute empirical P-value
    output$Pvalue <- length(which(output$statistics.under.null>=output$TrueT))/B
  return(output)
}

