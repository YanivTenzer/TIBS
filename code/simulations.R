path = 'C:\\Users\\Or Zuk\\Dropbox\\BiasedSampling\\Code'  # change to your path
setwd(path)
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
library(foreach)
library(doSNOW)
library(parallel)
library(doParallel)
library(gdata)
library(copula)
library(mvtnorm)
library(xtable)
library(Matrix)
library(ggplot2)  
library(latex2exp)

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

hyperplane.prms<-c(-1,1,0)
alpha <- 0.05 # significane threshold 
plot.flag <- 1 # plot and save figures 
run.flag <- 1 # 1: run simulations again. 0: load simulations results from file if they're available

prms = c()

# Vectors with different dependency settings 
dependence.type <- c('UniformStrip', 'Gaussian', 'Clayton', 'Gumbel', 
                     'LD', 'nonmonotone_nonexchangeable', 'CLmix', 'strictly_positive')
bias.type <- c('truncation', 'truncation', 'truncation', 'truncation', 
               'truncation', 'truncation', 'truncation', 'exponent_minus_sum_abs') # not good that we have only one simulation with positive W. Should add X+Y?
monotone.type <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) # is monotone
exchange.type <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE) # is exchangeable
# sample.size <- c(500, 500, 100, 100, 500, 500, 100, 500) 
sample.size <- c(500, 100, 100, 100, 100, 100, 100, 100) # set all sample.sizes to 100 
prms.rho <- list(0.3, seq(-0.9, 0.9, 0.1), 0.5, 1.6, c(0, 0.4),
                 seq(-0.9, 0.9, 0.1), 0.5, seq(-0.9, 0.9, 0.1)) # Parameters for each sampling type 

# test.type <- c('bootstrap', 'permutations', 'tsai', 'minP2') # different tests to run - new: add 
test.type <- c('bootstrap', 'permutations', 'fast-bootstrap', 'naive-bootstrap', 'naive-permutations', 
               'tsai', 'minP2') # different tests to run - new: add 
num.tests <- length(test.type)
iterations = 500 # 500  # Number of simulated dataset. Shared by all simulations
B = 10^2 # 10^4 # number of permtuations or bootstrap samples. Shared by all simulations 
num.sim <- length(dependence.type)

for(s in 2:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
{
  output.file <- paste0('results/', dependence.type[s], '_4tests_results_B_', 
                        B, '_', 'iters_', iterations) # set output file name
  print(paste0("s=", s))
  set.seed(1)
  num.prms <- length(prms.rho[[s]])
  if((!run.flag) & file.exists(paste0(output.file, '.Rdata')))
    load(file=paste0(output.file, '.Rdata'))
  else
    test.pvalue <- test.time <- array(-1, c(num.prms, num.tests, iterations)) 
  
  for(i.prm in 1:num.prms) # only 0.9 correlation num.prms) # loop on different simulation parameters - should plot for each of them? 
  {
    prms$rho = prms.rho[[s]][i.prm]
    prms$W.max <- 1 # temp: 1 for all weights
    prms$minp.eps <- "in"
    prms$keep.all <- 0 # set to 1 for plotting
    
    ## Parallel on multiple cores 
    ##    results.table<- foreach(i=seq(iterations), .combine=rbind) %dopar%{ 
    ##      iteration_result <- matrix(0, 4, B+1)
    for(i in 1:iterations) # run simulations again
    { 
      biased.data <- SimulateBiasedSample(sample.size[s], dependence.type[s], bias.type[s], prms) 
      all.data <- biased.data$all.data
      biased.data <- biased.data$data 
      if((!run.flag) & file.exists(paste0(output.file, '.Rdata'))) # simulate only once 
        break
      for(t in 1:6) # loop on statistical tests . Last test (7) minP2 is very slow 
      {
        prms$fast.bootstrap <- 0 # method for computing null expectations 
        prms$naive.expectation <- 0 # Check test with standard expectations
        # Skip irrelevant tests 
        if((!(bias.type[s] %in% c('truncation', 'Hyperplane_Truncation'))) & 
           (test.type[t] %in% c("tsai", 'minP2')))
          next  # these tests run only for truncation 
        if(test.type[t] == 'bootstrap')
        {
          if(bias.type[s] == 'huji')
            next
          if((bias.type[s] %in% c('truncation', 'Hyperplane_Truncation')) & (!exchange.type[s]))
            next # can't run bootstrap because w can be zero, unless we assume exchangability !!! 
        }      
        if((t %in% c(3,4,5)) & (s != 2))
          next # try comparison with naive settings only for Gaussian case (?) 
        cur.test.type <- test.type[t]
        switch(test.type[t], # Set test type
               'fast-bootstrap'={prms$fast.bootstrap <- 1
               cur.test.type <- 'bootstrap'},
               'naive-bootstrap'={prms$fast.bootstrap <- 1
               prms$naive.expectation <- 1
               cur.test.type <- 'bootstrap'},
               'naive-permutations'={prms$naive.expectation <- 1
               cur.test.type <- 'permutations'}
        ) 
        start.time <- Sys.time()
        test.results<-TIBS(biased.data, bias.type[s], B, cur.test.type, prms)
        test.time[i.prm, t, i] <- difftime(Sys.time(), start.time, units='secs')
        test.pvalue[i.prm, t, i] <- test.results$Pvalue
        print(paste0(dependence.type[s], ', rho=', prms$rho, '. Test: ', test.type[t], 
                     ' i=', i, ' of ', iterations, '. Pval=', round(test.pvalue[i.prm, t, i], 3), '. Time (sec.)=', round(test.time[i.prm, t, i], 2)))
      } # end loop on tests 
    }  # end loop on iterations (parallel collection). Simulation and testing for one parameter (i.prm)
    prms$title <- as.integer(s>1) 
    if(plot.flag) # plot example 
      PlotBiasedData(dependence.type[s], biased.data, prms)
  } # end simulation and testing for one dependency type (loop over i.prm)
  # Compute power (this is also type-1-error alpha under the null)
  test.power <- apply(test.pvalue<alpha, 1, rowMeans) 
  # Save in ONE file
  save(test.pvalue, test.time, test.power, prms.rho, sample.size, s, file=paste0(output.file, '.Rdata'))
  test.output <- t(cbind(test.power, colSums(rowSums(test.time, dims=2))))
  colnames(test.output) <- test.type
  rownames(test.output) <- c(prms.rho[[s]], 'time')
  test.output <- test.output[, test.time[1,,1]>=0] # take only relevant tests 
  print(xtable(test.output, type = "latex"), file = paste0(output.file, '.tex'), size="\\tiny") # new! save in latex format 
} # end loop on dependency types 

