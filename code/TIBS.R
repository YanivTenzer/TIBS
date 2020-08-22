
#' @export
add.one.sqrt <- function(x) {
  add_one_sqrt(x)
}



# Special function for arranging bootstrap marginals in the indices of the original sample 
MarginalsBootstrapOrganized <- function(data.marginals, bootstrap, w.fun, prms)
{
  if(w.fun %in% c('truncation', 'Hyperplane_Truncation'))
    marginals.n <- 2*prms$sample.size
  else
    marginals.n <- prms$sample.size
  
  marginals.bootstrap <- EstimateMarginals(bootstrap$sample, w.fun)
  marginals.bootstrap.organized <- c()
  marginals.bootstrap.organized$xy <- data.marginals$xy
  marginals.bootstrap.organized$PDFs <- matrix(0, marginals.n, 2) #   TrueT$marginals; marginals.bootstrap.organized$PDFs[] <-0 # reorder marginals
  if(w.fun %in% c('truncation', 'Hyperplane_Truncation'))
  {
    for(i in c(1:prms$sample.size))
      for(j in c(1:2))
      {
        marginals.bootstrap.organized$PDFs[bootstrap$indices[i,j], 1] <- 
          marginals.bootstrap.organized$PDFs[bootstrap$indices[i,j], 1] + 
          marginals.bootstrap$PDFs[i+(j-1)*prms$sample.size, 1]
        marginals.bootstrap.organized$PDFs[bootstrap$indices[i,j], 2] <- 
          marginals.bootstrap.organized$PDFs[bootstrap$indices[i,j], 2] + 
          marginals.bootstrap$PDFs[i+(j-1)*prms$sample.size, 1]
        
      } 
  } else 
    for(i in c(1:prms$sample.size))
      for(j in c(1:2))
        marginals.bootstrap.organized$PDFs[bootstrap$indices[i,j], j] <- 
          marginals.bootstrap.organized$PDFs[bootstrap$indices[i,j], j] + 
            marginals.bootstrap$PDFs[i, j] # NullT$marginals$PDFs[i,j]  # works for "Standard Marginals"
  
  
  marginals.bootstrap.organized$CDFs <- PDFToCDFMarginals(data.marginals$xy, marginals.bootstrap.organized$PDFs)  # here what happens if we duplicate? 
  return(marginals.bootstrap.organized)  
}

########################################################################
# 
# All steps of TIBS test: 
# 1. estimate marginals
# 2. compute null distribution
# 3. compute expectation table
# 4. compute statistic
# 
########################################################################
TIBS.steps <- function(data, w.fun, w.mat, grid.points, expectations.table, prms)
{
  marginals <- c()
  null.distribution <- c()
  if(prms$naive.expectation) # here we ignore W (using statistic for unbiased sampling)
  {
    use.w <- 'naive'
    w.mat <- 1
  }
  else 
    use.w <- w.fun
  if(prms$use.cpp)  # run in cpp 
  {
    if(missing(expectations.table) | isempty(expectations.table))
    {
      marginals <- EstimateMarginals_rcpp(data, use.w)
      if((!prms$naive.expectation) & (missing(w.mat) | isempty(w.mat)))
        w.mat = w_fun_to_mat_rcpp(marginals$xy, w.fun) # compute W again for augmented data
      null.distribution <- GetNullDistribution_rcpp(marginals$PDF, as.matrix(w.mat))
      expectations.table <- prms$sample.size * QuarterProbFromBootstrap_rcpp(marginals$xy, null.distribution$distribution, grid.points)          
    }
    T <- ComputeStatistic_rcpp(data, grid.points, expectations.table)  # keep same grid points for bootstrap sample?
  } else # use R
  {
    if(missing(expectations.table) | isempty(expectations.table))
    {
      marginals <- EstimateMarginals(data, use.w)
      if((!prms$naive.expectation) & (missing(w.mat) | isempty(w.mat)))
        w.mat = w_fun_to_mat(marginals$xy, w.fun) # compute W again for augmented data
      null.distribution <- GetNullDistribution(marginals$PDF, w.mat)
      expectations.table <- prms$sample.size * QuarterProbFromBootstrap(marginals$xy, null.distribution$distribution, grid.points)             
    }
    T <- ComputeStatistic(data, grid.points, expectations.table)$Statistic  # keep same grid points for bootstrap sample?
  }  
  return(list(Statistic=T, expectations.table=expectations.table, marginals=marginals, w.mat=w.mat, null.distribution=null.distribution))
}

########################################################################
# Perform Test for Independence under general Biased Sampling (TIBS)
# 
# Parameters: 
# data - n*2 array of (x,y) samples (NOT LIST!)
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
  
  if(!is.numeric(data))   # unlist and keep dimensions for data 
    data <- array(as.numeric(unlist(data)), dim(data))  
  
  
  # Set defaults
  if(!('use.cpp' %in% names(prms)))  # new: a flag for using c++ code 
    prms$use.cpp <- 0
  if(!('fast.bootstrap' %in% names(prms)))
    prms$fast.bootstrap <- FALSE
  if(!('minp.eps' %in% names(prms)))
    prms$minp.eps <- NULL # default: let permDep algorithm select minp.eps
  if(!('PL.expectation' %in% names(prms)))  # get expectation form the bootstrap: product-limit estimator
    prms$PL.expectation <- FALSE
  if(!('naive.expectation' %in% names(prms)))
    prms$naive.expectation <- FALSE
  if(!('delta' %in% names(prms)))
    prms$delta <- NA
  n <- dim(data)[1]
  if(!('w.max' %in% names(prms)))
    prms$w.max <- max(w_fun_to_mat(data, w.fun)) # update max
  
  if(prms$use.cpp) # use cpp code 
  {
    library(Rcpp)
    library(RcppArmadillo)
    Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 
  }    
  
  
  ############### ##################################################
  # 1.Compute weights matrix W: (not needed here, just for permutations test)
  #  w.mat=w_fun_to_mat(data, w.fun)
  # 2.Create a grid of points, based on the data:
  #  grid.points <- cbind(data[,1], data[,2])  # keep original points 
  grid.points <- unique.matrix(data)  # set unique for ties? for discrete data
  
  print("START TIBS")
  
  switch(test.type,
         'bootstrap'={
           TrueT <- TIBS.steps(data, w.fun, c(), grid.points, c(),  prms)  # compute statistic 
           
           #2. Compute statistic for bootstrap sample:
           statistics.under.null = matrix(0, prms$B, 1)
           #           null.distribution.bootstrap <- null.distribution
           for(ctr in 1:prms$B) # heavy loop: run on bootstrap 
           {
             if(prms$use.cpp)
               bootstrap <- Bootstrap_rcpp(TrueT$marginals$xy, TrueT$marginals$CDFs, w.fun, prms, dim(data)[1]) # draw new sample. Problem: which pdf and data? 
             else
               bootstrap <- Bootstrap(TrueT$marginals$xy, TrueT$marginals$PDFs, w.fun, prms, dim(data)[1]) 
             if(!prms$fast.bootstrap) # re-estimate marginals for null expectation for each bootstrap sample
             {
               new.bootstrap=1
               if(!new.bootstrap)
               {
                 ttt <- TIBS.steps(bootstrap$sample, w.fun, c(), grid.points, c(), prms)
                 statistics.under.null[ctr] <- ttt$Statistic
#                 print("OLD")  # TEMP 
               }
               else # do all steps except recomputing w_mat 
               {
                 marginals.bootstrap.new <- MarginalsBootstrapOrganized(TrueT$marginals, bootstrap, w.fun, prms)
                 
                 null.distribution.bootstrap.new <- GetNullDistribution(marginals.bootstrap.new$PDFs, TrueT$w.mat) # keep w_mat of ORIGINAL DATA! 
                 expectations.table.new <- prms$sample.size * QuarterProbFromBootstrap(
                   marginals.bootstrap.new$xy, null.distribution.bootstrap.new$distribution, grid.points)
                 statistics.under.null[ctr] <- ComputeStatistic(bootstrap$sample, grid.points, expectations.table.new)$Statistic # NEW! Compute null statistic without recomputing the entire matrix !!                
                 
#                 z0 <- max(abs(expectations.table.new - ttt$expectations.table))
#                 if(z0>0.00001)
#                   print(paste0("Should be zero expected table:", z0))
#                 z1 <- statistics.under.null[ctr] - ttt$Statistic
#                 if(z1>0.00001)
#                   print(paste0("Should be zero statistic:", z1))
               }               
             }
             else # use same expectation as original sample 
             {
               if(prms$use.cpp)
                 statistics.under.null[ctr] <- ComputeStatistic_rcpp(bootstrap$sample, grid.points, TrueT$expectations.table)
               else
               {
                 statistics.under.null[ctr] <- ComputeStatistic(bootstrap$sample, grid.points, TrueT$expectations.table)$Statistic  # keep same grid points for bootstrap sample?
                 #                 NullT.cpp <- ComputeStatistic_rcpp(bootstrap$sample, grid.points, TrueT$expectations.table)
                 #                 print(paste0("should be zero ComputeStatistic: ", max(abs(statistics.under.null[ctr] - NullT.cpp ))))
               }
             }
           }  # end loop on B bootstrap samples 
           output<-list(TrueT=TrueT$Statistic, statistics.under.null=statistics.under.null)
         },
         'bootstrap_inverse_weighting'={
           if(prms$use.cpp)
           {
             marginals <- EstimateMarginals_rcpp(data, w.fun)
             null.distribution <- GetNullDistribution_rcpp(marginals$PDFs, w.mat)
           } else {
             marginals <- EstimateMarginals(data, w.fun)
             null.distribution <- GetNullDistribution(marginals$PDFs, w.mat)
           }
           TrueT <- ComputeStatistic.W(data, grid.points, w.fun)$Statistic
           
           statistics.under.null=matrix(0, prms$B, 1)
           for(ctr in 1:prms$B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Run Boots=", ctr))
             bootstrap <- Bootstrap(marginals$xy, marginals$PDFs, w.fun, prms, dim(data)[1]) # draw new sample. Problem: which pdf and data? 
             NullT <- ComputeStatistic.W(bootstrap$sample, grid.points, w.fun)
             statistics.under.null[ctr] <- NullT$Statistic
           }
           output<-list(TrueT=TrueT,statistics.under.null=statistics.under.null)
         },
         'permutations'={
           w.mat=w_fun_to_mat(data, w.fun)
           if(prms$use.cpp)
             Permutations <- PermutationsMCMC_rcpp(w.mat, prms)
           else
             Permutations <- PermutationsMCMC(w.mat, prms)
           P=Permutations$P
           Permutations=Permutations$Permutations
           if((!prms$naive.expectation) & (!prms$PL.expectation))
           {
             if(prms$use.cpp)
               expectations.table <- QuarterProbFromPermutations_rcpp(data, P, grid.points)  # Permutations
             else
               expectations.table <- QuarterProbFromPermutations(data, P, grid.points)
           } else
             expectations.table <- c()
           #           print(expectations.table[1:5,])
           #           print(P[1:5,1:5]) # show P_ij 
           TrueT <- TIBS.steps(data, w.fun, w.mat, grid.points, expectations.table, prms)  # compute statistic. Use permutations for expected table 
           permuted.data <- cbind(data[,1], data[Permutations[,1],2]) # save one example 
           
           #Compute the statistics value for each permutation:
           statistics.under.null = matrix(0, prms$B, 1) # no need for new expectations.table 
           for(ctr in 1:prms$B)  
           {
             if(prms$use.cpp)  # new: use cpp
               statistics.under.null[ctr] <- ComputeStatistic_rcpp( # need to modify to calculate grid.points inside function!
                 cbind(data[,1], data[Permutations[,ctr],2]), grid.points, TrueT$expectations.table)
             else 
               statistics.under.null[ctr] <- ComputeStatistic( # here calculate grid.points inside function ! 
                 cbind(data[,1], data[Permutations[,ctr],2]), grid.points, TrueT$expectations.table)$Statistic
           }
           output<-list(TrueT=TrueT$Statistic, statistics.under.null=statistics.under.null, Permutations=Permutations)
         },  # end permutations test 
         'permutations_inverse_weighting'={ # statistic Our MCMC permutations 
           #             w.mat = w_fun_to_mat(data, w.fun)
           TrueT = ComputeStatistic.W(data, grid.points, w.fun)$Statistic
           w.mat=w_fun_to_mat(data, w.fun)
           Permutations=PermutationsMCMC(w.mat, prms) # burn.in=prms$burn.in, Cycle=prms$Cycle)
           Permutations=Permutations$Permutations
           #Compute the statistics value for each permutation:
           statistics.under.null = matrix(0, prms$B, 1)
           for(ctr in 1:prms$B) 
           {
             statistics.under.null[ctr] <- ComputeStatistic.W(cbind(data[,1], data[Permutations[,ctr],2]), grid.points, w.fun)$Statistic # grid.points calculated inside function
             #             statistics.under.null.cpp <- ComputeStatistic_w_rcpp(cbind(data[,1], data[Permutations[,ctr],2]), grid.points, w.fun) # grid.points calculated inside function
             #             print(paste0("should be zero ComputeStatisticW: ", max(abs(statistics.under.null.cpp - statistics.under.null[ctr] ))))
           }
           output<-list(TrueT=TrueT, statistics.under.null=statistics.under.null, Permutations=Permutations)
         },  # end permutations with inverse weighting test 
         'uniform_importance_sampling' = {  # new uniform importance sampling permutations test with the weighted Hoeffding statistic - only for positive W 
           w.mat=w_fun_to_mat(data, w.fun)
           if(prms$use.cpp)
             Permutations <- PermutationsMCMC_rcpp(w.mat, prms)
           else
             Permutations <- PermutationsMCMC(w.mat, prms)
           P=Permutations$P
           Permutations=Permutations$Permutations
           if((!prms$naive.expectation) & (!prms$PL.expectation))
           {
             if(prms$use.cpp)
               expectations.table <- QuarterProbFromPermutations_rcpp(data, P, grid.points)  # Permutations
             else
               expectations.table <- QuarterProbFromPermutations(data, P, grid.points)
           } else
             expectations.table <- c()
           TrueT <- TIBS.steps(data, w.fun, w.mat, grid.points, expectations.table, prms)  # compute statistic. Use permutations for expected table 
           permuted.data <- cbind(data[,1], data[Permutations[,1],2]) # save one example 
           
           output <- IS.permute(data, prms$B, w.fun, TrueT$expectations.table) # W)  # ComputeStatistic.W(dat, grid.points, w.fun)
         },
         'uniform_importance_sampling_inverse_weighting' = {  # new uniform importance sampling permutations test with our Hoeffding statistic - only for positive W
           output <- IS.permute(data, prms$B, w.fun) # W)  # ComputeStatistic.W(dat, grid.points, w.fun)
         }, 
         'tsai' = {result <- Tsai.test(data[,1],data[,2]) # TsaiTestTies  Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
         output <- list(Pvalue=result[4])
         },
         'minP2' = { library(permDep)  #MinP2 test, relevant only for truncation W(x,y)=1_{x<=y}
           require(survival)  
           # library(permDep) # use new library installed from github: https://github.com/stc04003/permDep (not CRAN)
           if(!is.na(prms$delta)) # when running minP we must have a delta parameter 
             dat <- data.frame(list(trun = data[,1], obs = data[,2], delta = prms$delta))
           else
             dat <- data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           results <- permDep(dat$trun, dat$obs, prms$B, dat$delta, nc = 4, minp2Only = TRUE, kendallOnly = FALSE) # set number of cores 
           ##                                      sampling = 'conditional') #  minp.eps= prms$minp.eps) # ,  new! set also min epsilon
           output <- list(Pvalue=results$p.valueMinp2)
         }
  )
  if(exists("permuted.data"))
    output$permuted.data <- permuted.data
  else
    output$permuted.data <- NA
  if(!("Pvalue" %in% names(output))) # Compute empirical P-value
    output$Pvalue <- length(which(output$statistics.under.null>=output$TrueT))/prms$B
  return(output)
}
