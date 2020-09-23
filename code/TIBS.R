
#' @export
add.one.sqrt <- function(x) {
  add_one_sqrt(x)
}



# Special function for arranging bootstrap marginals in the indices of the original sample 
MarginalsBootstrapOrganized <- function(data.marginals, bootstrap, w.fun, prms)
{
  if(!is.function(w.fun) && (w.fun %in% c('truncation', 'Hyperplane_Truncation')))
    marginals.n <- 2*prms$sample.size
  else
    marginals.n <- prms$sample.size
  
  marginals.bootstrap <- EstimateMarginals(bootstrap$sample, w.fun)
  marginals.bootstrap.organized <- c()
  marginals.bootstrap.organized$xy <- data.marginals$xy
  marginals.bootstrap.organized$PDFs <- matrix(0, marginals.n, 2) #   TrueT$marginals; marginals.bootstrap.organized$PDFs[] <-0 # reorder marginals
  if(!is.function(w.fun) && (w.fun %in% c('truncation', 'Hyperplane_Truncation')))
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
  
  
  ## NOT NEEDED!!!  marginals.bootstrap.organized$CDFs <- PDFToCDFMarginals(data.marginals$xy, marginals.bootstrap.organized$PDFs)  # here what happens if we duplicate? 
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
  if(!('naive.expectation' %in% names(prms)))
    prms$naive.expectation <- FALSE
  
  marginals <- c()
  null.distribution <- c()
  if(prms$naive.expectation) # here we ignore W (using statistic for unbiased sampling)
  {
    use.w <- 'naive'
    w.mat <- 1
  }
  else 
    use.w <- w.fun
  if(prms$use.cpp && !(is.function(use.w)))  # run in cpp 
  {
    if(missing(expectations.table) || isempty(expectations.table))
    {
      marginals <- EstimateMarginals_rcpp(data, use.w)
      if((!prms$naive.expectation) & (missing(w.mat) | isempty(w.mat)))
        w.mat = w_fun_to_mat_rcpp(marginals$xy, w.fun) # compute W again for augmented data
      null.distribution <- GetNullDistribution_rcpp(marginals$PDFs, as.matrix(w.mat))
      expectations.table <- prms$sample.size * QuarterProbFromBootstrap_rcpp(marginals$xy, null.distribution$distribution, grid.points)          
    }
    T <- ComputeStatistic_rcpp(data, grid.points, expectations.table)  # keep same grid points for bootstrap sample?
  } else # use R
  {
    if(missing(expectations.table) || isempty(expectations.table))
    {
      marginals <- EstimateMarginals(data, use.w)
      if((!prms$naive.expectation) & (missing(w.mat) | isempty(w.mat)))
        w.mat = w_fun_to_mat(marginals$xy, w.fun) # compute W again for augmented data
      null.distribution <- GetNullDistribution(marginals$PDFs, w.mat)
      expectations.table <- prms$sample.size * QuarterProbFromBootstrap(marginals$xy, null.distribution$distribution, grid.points)             
    }
    T <- ComputeStatistic(data, grid.points, expectations.table)$Statistic  # keep same grid points for bootstrap sample?
  }  
  
  #  print("expectations sum inside tibs.steps")
  #  print(rowSums(expectations.table))
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
TIBS <- function(data, w.fun, prms, test.method, test.stat)
{  
  #  print("Start TIBS")
  
  ##  library(PerMallows) # distance between permutations 
  ##  library(pracma)
  ##  source('utilities.R')
  ##  source('Tsai_test.R')
  ##  source('simulate_biased_sample.R')
  ##  source('marginal_estimation.R')
  
  epsilon <- 0.0000000001 # tolerance 
  
  #  print("called libraries TIBS")
  
  
  
  # Set defaults
  # new! we separate test statistic and test method 
  if(test.method %in% c('tsai', 'minP2'))
    test.stat <- test.method 
  if(test.stat %in% c('tsai', 'minP2'))
    test.method <- test.stat 
  
  if(missing(test.stat) || is.na(test.stat) || (test.stat==""))  # new: a flag for the test statistic
    test.stat <- "adjusted_w_hoeffding"
  if(missing(test.method) || is.na(test.method) || (test.method==""))  # new: a flag for method for assigning p-value (and statistic) 
    test.method <- "permutations"
  
  
  
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
  else
  {
    if(!(test.stat %in% c('tsai', 'minP2')))  # remove cenosred points
      data <- data[prms$delta==1,]
  }
  if(!is.numeric(data))   # unlist and keep dimensions for data 
    data <- array(as.numeric(unlist(data)), dim(data))  
  
  #  print("SET GRID")  
  
  if(!('perturb.grid' %in% names(prms)))  # default: perturb grid to avoid ties 
    prms$perturb.grid <- TRUE
  if(!('counts.flag' %in% names(prms)))  # new: a flag for using counts (not frequencies) in test statistic 
    prms$counts.flag <- 1 # default is counts statistic 
  if(!('include.ID' %in% names(prms)))  # new: a flag for including ID permutation in p-value calculation  
    prms$include.ID <- 1 # default is include 
  
  
  prms$sample.size <- n <- dim(data)[1]
  
  if(('grid.points' %in% names(prms)))  # new! enable input grid !  
  {
    grid.points <- prms$grid.points
    prms$perturb.grid <- FALSE
  } else
  {
    grid.points <- data # no unique !!! unique.matrix(data)  # set unique for ties? for discrete data
    if(prms$perturb.grid)  # new! add a small perturbation to avoid ties with data 
      grid.points = grid.points + epsilon * matrix( rnorm(2*n), n, 2)
  }
  
  
  
  #  print(paste0("RUN n=", n))
  
  #  if(!('w.max' %in% names(prms)))
  #    prms$w.max <- max(w_fun_to_mat(data, w.fun)) # update max
  
  if(prms$use.cpp && is.character(w.fun) && (!(test.stat %in% c('tsai', 'minP2')))) # use cpp code for our hard-coded tests and w 
  {
    #    print("Run TIBS CPP!!")
    #    library(Rcpp) # for some reason calling library causes sometimes crush of code 
    #    library(RcppArmadillo)
    #    Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 
    
    #    print("Call TIBS_RCPP")
    TR <- TIBS_rcpp(data, w.fun, prms, test.method, test.stat)
    #    print("Finished TIBS_RCPP")
    #    print("PVALUE:")
    #    print(TR$Pvalue)
    #    print("output names")
    #    print(names(TR))
    return(TR)
    #    return(TIBS_rcpp(data, w.fun, test.type, prms)) # run all in cpp 
  } # else
  # if(test.type == 'bootstrap')
  #  prms$use.cpp = 0 # make sure everything runs in R (only for esitmate marginals)
  
  
  
  switch(test.method,
         'bootstrap'={
           
           switch(test.stat, 
                  "adjusted_w_hoeffding"={
                    TrueT <- TIBS.steps(data, w.fun, c(), grid.points, c(),  prms)  # compute statistic 
                  }, 
                  "inverse_w_hoeffding"={
                    TrueT <- ComputeStatistic.W(data, grid.points, w.fun, prms$counts.flag)
                    if(prms$use.cpp)
                      marginals <- EstimateMarginals_rcpp(data, w.fun)
                    else 
                      marginals <- EstimateMarginals(data, w.fun)
                    w.mat <- w_fun_to_mat(marginals$xy, w.fun)
                  })  # end switch on test stat          
           
           statistics.under.null = matrix(0, prms$B, 1)
           #           null.distribution.bootstrap <- null.distribution
           for(ctr in 1:prms$B) # heavy loop: run on bootstrap 
           {
             switch(test.stat, 
                    "adjusted_w_hoeffding"={
                      
                      if(prms$use.cpp && !(is.function(w.fun)))  # here we must use sorted marginals 
                      {
                        bootstrap <- Bootstrap_rcpp(TrueT$marginals$xy, TrueT$marginals$CDFs, 
                                                    TrueT$w.mat, prms, dim(data)[1]) # TrueT$w.mat draw new sample. Problem: which pdf and data? 
                        bootstrap$indices <- bootstrap$indices + 1 # adjust from cpp to R
                      }
                      else
                        bootstrap <- Bootstrap(TrueT$marginals$xy, TrueT$marginals$PDFs, w.fun, prms, dim(data)[1]) 
                      
                      if(!prms$fast.bootstrap) # re-estimate marginals for null expectation for each bootstrap sample
                      {
                        marginals.bootstrap.new <- MarginalsBootstrapOrganized(TrueT$marginals, bootstrap, w.fun, prms)
                        
                        if(prms$use.cpp)
                        {
                          null.distribution.bootstrap.new <- GetNullDistribution_rcpp(marginals.bootstrap.new$PDFs, TrueT$w.mat) # TrueT$w.mat) # keep w_mat of ORIGINAL DATA! 
                          expectations.table.new <- prms$sample.size * QuarterProbFromBootstrap_rcpp(
                            marginals.bootstrap.new$xy, null.distribution.bootstrap.new$distribution, grid.points)
                          statistics.under.null[ctr] <- ComputeStatistic_rcpp(bootstrap$sample, grid.points, expectations.table.new) # NEW! Compute null statistic without recomputing the entire matrix !!              
                        }  else
                        {
                          null.distribution.bootstrap.new <- GetNullDistribution(marginals.bootstrap.new$PDFs, TrueT$w.mat) # keep w_mat of ORIGINAL DATA! 
                          expectations.table.new <- prms$sample.size * QuarterProbFromBootstrap(
                            marginals.bootstrap.new$xy, null.distribution.bootstrap.new$distribution, grid.points)
                          statistics.under.null[ctr] <- ComputeStatistic(bootstrap$sample, grid.points, expectations.table.new)$Statistic # NEW! Compute null statistic without recomputing the entire matrix !!          
                        }
                      }
                      else # use same expectation as original sample 
                      {
                        if(prms$use.cpp)
                          statistics.under.null[ctr] <- ComputeStatistic_rcpp(bootstrap$sample, grid.points, TrueT$expectations.table)
                        else
                          statistics.under.null[ctr] <- ComputeStatistic(bootstrap$sample, grid.points, TrueT$expectations.table)$Statistic  # keep same grid points for bootstrap sample?
                      }
                    }, 
                    "inverse_w_hoeffding"={
                      bootstrap <- Bootstrap(marginals$xy, marginals$PDFs, w.fun, prms, dim(data)[1]) # draw new sample. Problem: which pdf and data? 
                      statistics.under.null[ctr] <- ComputeStatistic.W(bootstrap$sample, grid.points, w.fun, prms$counts.flag)$Statistic
                    }) # end switch test stat
             
           }  # end loop on B bootstrap samples 
           
           
           output<-list(TrueT=TrueT$Statistic, statistics.under.null=statistics.under.null)
         },
         'permutationsMCMC'={
           w.mat = w_fun_to_mat(data, w.fun)
           if(prms$use.cpp)
             Permutations <- PermutationsMCMC_rcpp(w.mat, prms)
           else
             Permutations <- PermutationsMCMC(w.mat, prms)
           P=Permutations$P
           Permutations=Permutations$Permutations
           switch(test.stat, 
                  "adjusted_w_hoeffding"={
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
                  }, 
                  "inverse_w_hoeffding"={
                    TrueT <- ComputeStatistic.W(data, grid.points, w.fun, prms$counts.flag)
                  }) # end switch test.stat
           
           #Compute the statistics value for each permutation:
           statistics.under.null = matrix(0, prms$B, 1) # no need for new expectations.table 
           for(ctr in 1:prms$B)  
           {
             switch(test.stat, 
                    "adjusted_w_hoeffding"={
                      if(prms$use.cpp)  # new: use cpp
                        statistics.under.null[ctr] <- ComputeStatistic_rcpp( # need to modify to calculate grid.points inside function!
                          cbind(data[,1], data[Permutations[,ctr],2]), grid.points, TrueT$expectations.table)
                      else 
                        statistics.under.null[ctr] <- ComputeStatistic( # here calculate grid.points inside function ! 
                          cbind(data[,1], data[Permutations[,ctr],2]), grid.points, TrueT$expectations.table)$Statistic
                    }, 
                    "inverse_w_hoeffding"={
                      statistics.under.null[ctr] <- ComputeStatistic.W(cbind(data[,1], data[Permutations[,ctr],2]), grid.points, w.fun, prms$counts.flag)$Statistic # grid.points calculated inside function
                    }) # end switch test stat
           }
           output<-list(TrueT=TrueT$Statistic, statistics.under.null=statistics.under.null, Permutations=Permutations)
         },  # end permutations test 
         "permutationsIS" = {  # importance sampling. The IS distribution, and test statistic are encoded inside prms : prms$importance.sampling.dist, and test.stat 
           output <- IS.permute(data, grid.points, w.fun, prms, test.stat)
         },
         'tsai' = {
           if(!('delta' %in% names(prms)) || any(is.na(prms$delta)))  # new: Tsai with censoring
             result <- Tsai.test(data[,1], data[,2]) # TsaiTestTies  Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
           else
             result <- Tsai.test(data[,1], data[,2], prms$delta) # TsaiTestTies  Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
           output <- list(Pvalue=result[4], TrueT=result[3])
         },
         'minP2' = { library(permDep) #MinP2 test, relevant only for truncation W(x,y)=1_{x<=y}
           require(survival)  
           # use new library installed from github: https://github.com/stc04003/permDep (not CRAN)
           if(!any(is.na(prms$delta))) # when running minP we must have a delta parameter 
             dat <- data.frame(list(trun = data[,1], obs = data[,2], delta = prms$delta))
           else
             dat <- data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           results <- permDep(dat$trun, dat$obs, prms$B, dat$delta, nc = 4, minp2Only = TRUE, kendallOnly = FALSE) # set number of cores 
           ##                                      sampling = 'conditional') #  minp.eps= prms$minp.eps) # ,  new! set also min epsilon
           output <- list(Pvalue=results$p.valueMinp2)
         }
  )
  
  if(!("permuted.data" %in% names(output)) && exists("permuted.data"))
    output$permuted.data <- permuted.data
  if(!("permuted.data" %in% names(output)))
    output$permuted.data <- NA
  prms$include.ID = as.numeric(prms$include.ID>0)  # set to zero or one
  if(!("Pvalue" %in% names(output))) # Compute empirical P-value
    output$Pvalue <- (length(which(output$statistics.under.null>=output$TrueT))+prms$include.ID) / (prms$B+prms$include.ID) # new! include also id permutation. Pval between 1/(N+1) to N/(N+1)
  return(output)
}
