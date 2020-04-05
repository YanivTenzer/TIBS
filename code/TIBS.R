########################################################################
# Perform Test for Independence under general Biased Sampling (TIBS)
# 
# Parameters: 
# data - n*2 array of (x,y) samples
# bias.type - string indicating biased sampling function W 
# test.type - test to perform 
# prms - additional parameters (needed for bootstrap) # B - number of bootstrap/permutation samples to perform 
########################################################################
TIBS <- function(data, bias.type, test.type, prms)
{  
  library(PerMallows) # distance between permutations 
  library(pracma)
  source('utilities.R')
  source('Tsai_test.R')
  source('simulate_biased_sample.R')
  source('marginal_estimation.R')
  
  # Set defaults
  if(!('fast.bootstrap' %in% names(prms)))
    prms$fast.bootstrap <- 0
  if(!('minp.eps' %in% names(prms)))
    prms$minp.eps <- NULL # default: let permDep algorithm select minp.eps
  if(!('PL.expectation' %in% names(prms)))  # what's that?
    prms$PL.expectation <- FALSE
  if(!('naive.expectation' %in% names(prms)))
    prms$naive.expectation <- 0
  
  #################################################################
  # 1.Compute weights matrix W:  
  W=GetBiasedSamplingWeights(data, dim(data)[1], bias.type)
  # 2.Create a grid of points, based on the data:
#  permutation<-PermutationsMCMC(W, dim(data)[1], list(B=1)) # why always sample a permutation using MCMC? 
###  permutation <- rev(dim(data)[1]:1)  # always use the same permutation and grid - TEMP
###  temp.data<-cbind(data[,1], data[permutation,2])
###  permutation<-PermutationsMCMC(W, 1, dim(data)[1], cyc) # why always sample a permutation using MCMC? 
  #temp.data<-cbind(data[,1], data[permutation,2])
  grid.points <- cbind(data[,1], data[,2])  # keep original points 
  grid.points <- unique.matrix(grid.points)  # why set unique? we want to increase the density/weight for ties.

  
    # Discard the extremum points to avoid numerical issues
  # idx.minmax <- which( (temp.data[,1] %in% c(min(temp.data[,1],max(temp.data[,1])))) | 
  #                      (temp.data[,2] %in% c(min(temp.data[,2],max(temp.data[,2])))) )
  # grid.points <- temp.data[-idx.minmax,]
  # grid.points<- unique.matrix(temp.data[-idx.minmax,])
  
  switch(test.type,
         'bootstrap'={
           #3. Estimate the marginals
           marginals <- EstimateMarginals(data, bias.type)
           W = GetBiasedSamplingWeights(marginals$xy, dim(marginals$xy)[1], bias.type) # compute W again for augmented data
           
           #4. Estimate W(x,y)*Fx*Fy/normalizing.factor
           if(prms$naive.expectation) # here we ignore W (using statistic for unbiased sampling)
           {
             marginals.naive <- EstimateMarginals(data, 'naive')
             null.distribution <- GetNullDistribution(marginals.naive$PDF, 1)
             expectations.table <- QuarterProbFromBootstrap(
               marginals.naive$xy, null.distribution$null.distribution, grid.points)
           } else
           {
             null.distribution <- GetNullDistribution(marginals$PDF, W)
             expectations.table <- QuarterProbFromBootstrap(
               marginals$xy, null.distribution$null.distribution, grid.points)             
           }

           #1. First compute the statistics based on the original data set:
           #           expectations.table<-QuarterProbFromBootstrap(marginals$xy, null.distribution, grid.points) # data
           TrueT=ComputeStatistic(data, grid.points, expectations.table)
           obs.table <- TrueT$obs.table
           TrueT <- TrueT$Statistic

           #2. Compute statistic for bootstrap sample:
           statistics.under.null=matrix(0, prms$B, 1)
           null.distribution.bootstrap<-null.distribution
           for(ctr in 1:prms$B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Run Boots=", ctr))
             bootstrap.sample <- Bootstrap(marginals$xy, marginals$PDF, bias.type, prms, dim(data)[1]) # draw new sample. Problem: which pdf and data? 
             if(!prms$fast.bootstrap) # re-estimate marginals for null expectation for each bootstrap sample
             {
               marginals.bootstrap <- EstimateMarginals(bootstrap.sample, bias.type)  # Why are the marginals estimated each time? 
               #3. Compute weights matrix W:    
               W.bootstrap=GetBiasedSamplingWeights(marginals.bootstrap$xy, dim(marginals.bootstrap$xy)[1], bias.type)
               #4. Estimate W(x,y)*Fx*FY/normalizing.factor
               ifelse(prms$naive.expectation,  # here we ignore W (using statistic for unbiased sampling)
                 null.distribution.bootstrap <- GetNullDistribution(marginals.bootstrap$PDFs, 1),
                 null.distribution.bootstrap <- GetNullDistribution(marginals.bootstrap$PDFs, W.bootstrap))
               
               expectations.table <- QuarterProbFromBootstrap(
                 marginals.bootstrap$xy, null.distribution.bootstrap$null.distribution, grid.points)
             } 
             NullT <- ComputeStatistic(bootstrap.sample, grid.points, expectations.table)
             null.obs.table <- NullT$obs.table
             statistics.under.null[ctr] <- NullT$Statistic
           }
           output<-list(TrueT=TrueT,statistics.under.null=statistics.under.null)
         },
         'permutations'={
           Permutations=PermutationsMCMC(W, dim(data)[1], prms) # burn.in=prms$burn.in, Cycle=prms$Cycle)
           P=Permutations$P
           Permutations=Permutations$Permutations
           if(prms$naive.expectation) # here we ignore W (using statistic for unbiased sampling)
           {
             marginals <- EstimateMarginals(data, 'naive')           
             null.distribution <- GetNullDistribution(marginals$PDF, 1)
             expectations.table <- QuarterProbFromBootstrap(marginals$xy, null.distribution$null.distribution, grid.points) # data
           } else
           {
              if(prms$PL.expectation)  # what's this? 
              {
                  marginals <- EstimateMarginals(data, bias.type)
                  W = GetBiasedSamplingWeights(marginals$xy, dim(marginals$xy)[1], bias.type)
                  null.distribution <- GetNullDistribution(marginals$PDF, W)
                  expectations.table <- QuarterProbFromBootstrap(marginals$xy, null.distribution$null.distribution, grid.points)
              } else
                  expectations.table <- QuarterProbFromPermutations(data, P, grid.points)  # Permutations
           }
    #             expectations.table <- QuarterProbFromPermutations(data, Permutations, grid.points)
           TrueT = ComputeStatistic(data, grid.points, expectations.table)$Statistic
           
           #Compute the statistics value for each permutation:
           statistics.under.null = matrix(0, prms$B, 1)
           for(ctr in 1:prms$B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Comp. Stat. Perm=", ctr))
             statistics.under.null[ctr] <- ComputeStatistic(
               cbind(data[,1], data[Permutations[,ctr],2]), grid.points, expectations.table)$Statistic
           }
           output<-list(TrueT=TrueT,statistics.under.null=statistics.under.null, Permutations=Permutations)
         },  # end permutations test 
         'tsai' = {result <- TsaiTestTies(data[,1],data[,2]) # Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
         output <- list(Pvalue=result[2])
         },
         'minP2' = { library(permDep)  #MinP2 test, relevant only for truncation W(x,y)=1_{x<=y}
           require(survival)  
           # library(permDep) # should use new library installed from github: https://github.com/stc04003/permDep (not CRAN)
           # QU: what should be min.eps?
###           dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           if(!is.na(prms$delta))
           {
              dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = prms$delta))
           }
           else{
              dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           }
           
           
           results<- permDep(dat$trun, dat$obs, prms$B, dat$delta, nc = 4, minp2Only = TRUE, kendallOnly = FALSE) # set number of cores 
##                                      sampling = 'conditional', kendallOnly = FALSE) #  minp.eps= prms$minp.eps) # ,  new! set also min epsilon
           output<-list(Pvalue=results$p.valueMinp2)
         }
  )
  output$permuted.data <- NA #temp.data # add example of permuted data 
  if(!("Pvalue" %in% names(output))) # Compute empirical P-value
    output$Pvalue <- length(which(output$statistics.under.null>=output$TrueT))/prms$B
  return(output)
}
