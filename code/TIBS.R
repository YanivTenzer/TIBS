Rcpp::sourceCpp("C/ToR.cpp")  # new: replace functions by their c++ version

########################################################################
# Perform Test for Independence under general Biased Sampling (TIBS)
# 
# Parameters: 
# data - n*2 array of (x,y) samples
# w.fun - biased sampling function W 
# test.type - test to perform 
# prms - additional parameters (needed for bootstrap) including B - number of bootstrap/permutation samples to perform 
# 
# Output:
# TrueT - test statistic for the data
# statistics.under.null - vector of statistics under the null 
########################################################################
TIBS <- function(data, w.fun, test.type, prms)
{  
  library(PerMallows) # distance between permutations 
  library(pracma)
  source('utilities.R')
  source('Tsai_test.R')
  source('simulate_biased_sample.R')
  source('marginal_estimation.R')
  
  use_cpp <- 0 # new: a flag for using c++ code 
  
  # Set defaults
  if(!('fast.bootstrap' %in% names(prms)))
    prms$fast.bootstrap <- 0
  if(!('minp.eps' %in% names(prms)))
    prms$minp.eps <- NULL # default: let permDep algorithm select minp.eps
  if(!('PL.expectation' %in% names(prms)))  # get expectation form the bootstrap: product-limit estimator
    prms$PL.expectation <- FALSE
  if(!('naive.expectation' %in% names(prms)))
    prms$naive.expectation <- 0
  
  #################################################################
  # 1.Compute weights matrix W: (not needed here, just for permutations test)
#  w.mat=w_fun_to_mat(data, w.fun)
  # 2.Create a grid of points, based on the data:
  grid.points <- cbind(data[,1], data[,2])  # keep original points 
  grid.points <- unique.matrix(grid.points)  # set unique for ties? for discrete data
 
  switch(test.type,
         'bootstrap_inverse_weighting'={
           marginals <- EstimateMarginals(data, w.fun)
           w.mat = w_fun_to_mat(marginals$xy, w.fun) 
           null.distribution <- GetNullDistribution(marginals$PDF, w.mat)
           TrueT <- ComputeStatistic_inverse_weighting(data, grid.points, w.mat)$Statistic
           
           statistics.under.null=matrix(0, prms$B, 1)
           for(ctr in 1:prms$B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Run Boots=", ctr))
             
             bootstrap.sample <- Bootstrap(marginals$xy, marginals$PDF, w.fun, prms, dim(data)[1]) # draw new sample. Problem: which pdf and data? 
             w.mat.bootstrap <- w_fun_to_mat(bootstrap.sample, w.fun)
             NullT <- ComputeStatistic_inverse_weighting(bootstrap.sample, grid.points, w.mat.bootstrap)
             statistics.under.null[ctr] <- NullT$Statistic
           }
           output<-list(TrueT=TrueT,statistics.under.null=statistics.under.null)
         },
         'bootstrap'={
           #3. Estimate the marginals
           marginals <- EstimateMarginals(data, w.fun)
           w.mat = w_fun_to_mat(marginals$xy, w.fun) # compute W again for augmented data
           
           #4. Estimate W(x,y)*Fx*Fy/normalizing.factor
           if(prms$naive.expectation) # here we ignore W (using statistic for unbiased sampling)
           {
             marginals.naive <- EstimateMarginals(data, 'naive')
             null.distribution <- GetNullDistribution(marginals.naive$PDF, 1)
             expectations.table <- QuarterProbFromBootstrap(
               marginals.naive$xy, null.distribution$null.distribution, grid.points)
           } else
           {
             null.distribution <- GetNullDistribution(marginals$PDF, w.mat)
             expectations.table <- QuarterProbFromBootstrap(
               marginals$xy, null.distribution$null.distribution, grid.points)             
           }

           #1. First compute the statistic based on the original data set:
           if(use_cpp)  # new: use cpp
           {
             print("USE C")
             TrueT=ComputeStatistic_cpp(n, data, grid.points, expectations.table)
           }
           
           else
           {
             print("USE R")
             TrueT=ComputeStatistic(data, grid.points, expectations.table)
           }
           # obs.table <- TrueT$obs.table
           TrueT <- TrueT$Statistic

           #2. Compute statistic for bootstrap sample:
           statistics.under.null=matrix(0, prms$B, 1)
           null.distribution.bootstrap<-null.distribution
           for(ctr in 1:prms$B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Run Boots=", ctr))
             bootstrap.sample <- Bootstrap(marginals$xy, marginals$PDF, w.fun, prms, dim(data)[1]) # draw new sample. Problem: which pdf and data? 
             if(!prms$fast.bootstrap) # re-estimate marginals for null expectation for each bootstrap sample
             {
               marginals.bootstrap <- EstimateMarginals(bootstrap.sample, w.fun)  # Why are the marginals estimated each time? 
               #3. Compute weights matrix W:    
               w.mat.bootstrap=w_fun_to_mat(marginals.bootstrap$xy, w.fun)
               #4. Estimate W(x,y)*Fx*FY/normalizing.factor
               ifelse(prms$naive.expectation,  # here we ignore W (using statistic for unbiased sampling)
                 null.distribution.bootstrap <- GetNullDistribution(marginals.bootstrap$PDFs, 1),
                 null.distribution.bootstrap <- GetNullDistribution(marginals.bootstrap$PDFs, w.mat.bootstrap))
               
               expectations.table <- QuarterProbFromBootstrap(
                 marginals.bootstrap$xy, null.distribution.bootstrap$null.distribution, grid.points)
             } 
             NullT <- ComputeStatistic(bootstrap.sample, grid.points, expectations.table)
             # null.obs.table <- NullT$obs.table
             statistics.under.null[ctr] <- NullT$Statistic
           }
           output<-list(TrueT=TrueT, statistics.under.null=statistics.under.null)
         },
         'permutations'={
           w.mat=w_fun_to_mat(data, w.fun)
           Permutations=PermutationsMCMC(w.mat, dim(data)[1], prms) # burn.in=prms$burn.in, Cycle=prms$Cycle)
           P=Permutations$P
           Permutations=Permutations$Permutations
           if(prms$naive.expectation) # here we ignore W (using statistic for unbiased sampling)
           {
             marginals <- EstimateMarginals(data, 'naive')           
             null.distribution <- GetNullDistribution(marginals$PDF, 1)
             expectations.table <- QuarterProbFromBootstrap(marginals$xy, null.distribution$null.distribution, grid.points) # data
           } else
           {
             if(prms$PL.expectation)  # get expectations from the bootstrap estimator
             {
               marginals <- EstimateMarginals(data, w.fun)
               w.mat = w_fun_to_mat(marginals$xy, w.fun)
               null.distribution <- GetNullDistribution(marginals$PDF, W)
               expectations.table <- QuarterProbFromBootstrap(marginals$xy, null.distribution$null.distribution, grid.points)
             } else
               expectations.table <- QuarterProbFromPermutations(data, P, grid.points)  # Permutations
           }
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
           output<-list(TrueT=TrueT, statistics.under.null=statistics.under.null, Permutations=Permutations)
         },  # end permutations test 
         'permutations_inverse_weighting'={
           w.mat = w_fun_to_mat(data, w.fun)
           TrueT = ComputeStatistic_inverse_weighting(data, grid.points, w.mat)$Statistic
          
           Permutations=PermutationsMCMC(W, dim(data)[1], prms) # burn.in=prms$burn.in, Cycle=prms$Cycle)
           Permutations=Permutations$Permutations
          #Compute the statistics value for each permutation:
           statistics.under.null = matrix(0, prms$B, 1)
           for(ctr in 1:prms$B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Comp. Stat. Perm=", ctr))
            
             permuted.sample =  cbind(data[,1], data[Permutations[,ctr],2])
             w.mat.permutation=w_fun_to_mat(permuted.sample, w.fun)
             NullT <- ComputeStatistic_inverse_weighting(permuted.sample, grid.points, w.mat.permutation)
             statistics.under.null[ctr] <- NullT
           }
           output<-list(TrueT=TrueT, statistics.under.null=statistics.under.null, Permutations=Permutations)
         },  # end permutations with inverse weighting test 
         'tsai' = {result <- TsaiTestTies(data[,1],data[,2]) # Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
         output <- list(Pvalue=result[2])
         },
         'minP2' = { library(permDep)  #MinP2 test, relevant only for truncation W(x,y)=1_{x<=y}
           require(survival)  
           # library(permDep) # use new library installed from github: https://github.com/stc04003/permDep (not CRAN)
           if(!is.na(prms$delta))
           {
              dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = prms$delta))
           }
           else{
              dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           }
           
           results<- permDep(dat$trun, dat$obs, prms$B, dat$delta, nc = 4, minp2Only = TRUE, kendallOnly = FALSE) # set number of cores 
##                                      sampling = 'conditional') #  minp.eps= prms$minp.eps) # ,  new! set also min epsilon
           output <-list(Pvalue=results$p.valueMinp2)
         }, 
         'importance.sampling' = {  # new importance sampling permutations test (not working yet)
           results <- IS.permute(dat, prms$B, w.fun) # W)  # ComputeStatistic.W(dat, grid.points, w.fun)
         }
  )
  output$permuted.data <- NA #temp.data # add example of permuted data 
  if(!("Pvalue" %in% names(output))) # Compute empirical P-value
    output$Pvalue <- length(which(output$statistics.under.null>=output$TrueT))/prms$B
  return(output)
}
