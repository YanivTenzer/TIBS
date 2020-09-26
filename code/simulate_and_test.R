rm(list=ls())
gc()

library(latex2exp)
library(xtable)
library(binom)  # for binomial confidence intervals 

library(stringr)
library(foreach)
library(doSNOW)
library(doParallel)
library(gdata)
library(mvtnorm)
library(ggplot2)  
library(pracma)
library(matrixStats)
library(PerMallows) # distance between permutations 
library(Matrix)

args=commandArgs(trailingOnly = TRUE)
Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 


source('simulate_biased_sample.R')
source('TIBS.R')
source("import_samp.R")
source('utilities.R')
## source('simulate_biased_sample.R')
## source('TIBS.R')
source('marginal_estimation.R')
## source('utilities.R')
## source('import_samp.R')
source('Tsai_test.R')



# Function running for one dataset: simulate data and perform weighted independent test 
# 
# Parameters:
# sequential.stopping (new!) - stop drawing permutations/bootstrap samples if we see that the p-value is very far from alpha 
# 
# Output:
# test.power - estimated power of the test
# test.output - structure containing power and running time 
# 
simulate_and_test <- function(dependence.type='Gaussian', prms.rho=c(0.0), w.fun='truncation',
                              # new: separate test.stat from test.method 
                              test.comb = rbind(c("tsai", ""), c("minP2", ""), c("permutationsMCMC", ""), c("bootstrap", ""), c("fast-bootstrap", ""), c("bootstrap", "hoeffding"), c("permutationsMCMC", "hoeffding")),
                              #                              test.type=c('tsai', 'minP2', 'permutations', 'bootstrap', 'fast-bootstrap', 'naive-bootstrap', 'naive-permutations'), 
                              prms) 
{
  if(is.character(prms)) # new: enable reading parameters from a parameters file: 
    load(prms)  # load parameters structure from file. Can also contain test.comb 

  if(prms$use.cpp)
  {
    library(Rcpp)
    library(RcppArmadillo)
  }
  
  
  if('seed' %in% names(prms))
    set.seed(prms$seed)
  
  if(!('simulate.once' %in% names(prms)))  # default: simulate a new dataset for each iteration. 
    prms$simulate.once = 0
  if(!('run.flag' %in% names(prms)))  # default: run test for each iteration . 
    prms$run.flag = 1
  if(!('run.sim' %in% names(prms)))  # default: simulate a new dataset for each iteration. 
    prms$run.sim = 1
  
  if(!('B' %in% names(prms)))  # default: simulate using 
    prms$B = 1000
  
  if(is.character(w.fun))
    w.fun.str <- w.fun
  else # w.fun given as a function 
    w.fun <- 'function'
  test.method <- test.comb[,1]
  test.stat <- test.comb[,2]
  test.params <- test.comb[,3]
  
  num.tests <- dim(test.comb)[1]
  #browser()
  output.file <- paste0('results/', dependence.type, '_w_', w.fun.str, '_all_tests_results_B_', 
                        prms$B, '_iters_', prms$iterations, '_n_', prms$sample.size, '_rho_', paste(prms.rho, collapse = '_')) # set output file name
  
  num.prms <- length(prms.rho)
  
  for(i.prm in 1:num.prms) # loop on different simulation parameters - should plot for each of them? 
  {
    prms$rho = prms.rho[i.prm] # [[s]]
    prms$w.rho = -prms$rho  # new: set negative rho for w. (Should speciffy from outside if you want a different value) 
    sim.file <- paste0('sim/', dependence.type, '_w_', w.fun.str, '_sim_iters_', 
                       prms$iterations, '_n_', prms$sample.size, '_rho_', prms.rho[i.prm], '.Rdata') # set simulated data file name for each rho 
    
    
    prms$minp.eps <- "in"
    if(!('w.max' %in% names(prms))) # set for each parameter separately 
      prms$w.max <- set_w_max(2*prms$sample.size, dependence.type, w.fun, prms)
    
    ## Parallel on multiple cores 
    ##    results.table<- foreach(i=seq(prms$iterations), .combine=rbind) %dopar%{ 
    ##      iteration_result <- matrix(0, 4, B+1)
    if(prms$run.sim == 0)   # here we load pre-existing data 
    {
      save.run.flag <- prms$run.flag
      save.use.cpp <- prms$use.cpp
      if(file.exists(sim.file))
        load(sim.file)  # load simulations data file
      else
      {
        all.biased.data <- array(0, c(prms$sample.size, 2, prms$iterations)) # save a matrix for all data 
        for(i in 1:prms$iterations) # run simulations multiple times 
        {
          if(prms$simulate.once && (i>1)) # all iterations are with the same simulated dataset 
            sim.data <- all.data
          else
            sim.data <- NA
          
          if( i%%10==0)
            print(paste0("sim i=", i, " out of ", prms$iterations))
          biased.data <- SimulateBiasedSample(prms$sample.size, 
                                              dependence.type, w.fun, 
                                              prms, sim.data) 
          if(i==1)
            all.data <- biased.data$all.data # not used 
          all.biased.data[,,i] <- biased.data$data
        }
        save(all.biased.data, prms, file=sim.file)  
      }
      prms$run.flag <- save.run.flag
      prms$use.cpp <- save.use.cpp
    } 
    
    
    # run flag: 1 - always run, 0 - never run, -1 - run only if file doesn't exist
    test.true.stat <- array(-1, c(num.prms, num.tests, prms$iterations)) 
    test.null.stat <- array(-1, c(num.prms, num.tests, prms$iterations, prms$B)) 
    if((prms$run.flag != 1) && file.exists(paste0(output.file, '.Rdata')))
    {
      save.prms <- prms
      load(file=paste0(output.file, '.Rdata'))
      if(!("run.sim" %in% names(prms)))
        prms$run.sim = save.prms$run.sim
      if(!("run.flag" %in% names(prms)))
        prms$run.flag = save.prms$run.flag
      if(!("w.rho" %in% names(prms)))
        prms$w.rho = save.prms$w.rho
      print(setdiff(names(save.prms), names(prms)))
    }
    else
    {
      test.pvalue <- test.time <- test.true.stat <- array(-1, c(num.prms, num.tests, prms$iterations)) 
    }
    
    
    for(i in 1:prms$iterations) # run simulations multiple times 
    { 
      if(prms$run.sim == 0) # load simulated data from file 
      {
#        if(i==1)  # no need to load again
#          load(sim.file)
        biased.data = all.biased.data[,,i]
      } else  # new simulation
      {
        if(prms$simulate.once && (i>1)) # all iterations are with the same simulated dataset 
          sim.data <- all.data
        else
          sim.data <- NA
        biased.data <- SimulateBiasedSample(prms$sample.size, 
                                            dependence.type, w.fun, 
                                            prms, sim.data) 
        if(i==1)
          all.data <- biased.data$all.data # not used 
        biased.data <- biased.data$data 
      }      
      #      print(paste0("sampling ratio: ", dim(all.data)[1] / dim(biased.data)[1]))      
      
      if((prms$run.flag==0) || ((prms$run.flag==-1) && file.exists(paste0(output.file, '.Rdata'))))#      if((!prms$run.flag) && file.exists(paste0(output.file, '.Rdata'))) # simulate only once 
        break
      
      for(t in c(1:num.tests))  # ALL TESTs !    #  c(1,2,6,7)) # loop on statistical tests . Last test (7) minP2 is very slow 
      {
        if(is.na(test.method[t]))
          test.method[t] <- ""
        if(is.na(test.stat[t]))
          test.stat[t] <- ""
        
        prms$fast.bootstrap <- 0 # method for computing null expectations 
        prms$naive.expectation <- 0 # Check test with standard expectations
        # Skip irrelevant tests 
        if((!(w.fun %in% c('truncation', 'Hyperplane_Truncation'))) & 
           (test.method[t] %in% c('tsai', 'minP2')))
          next  # these tests run only for truncation 
        
        if(((w.fun %in% c('truncation', 'Hyperplane_Truncation'))) && (test.method[t] == "permutationsIS") && !(test.params[t] %in% c("Tsai", "KouMcculough.w")))  # no importance sampling for truncation (should add Chen&Liu method)
          next
        
        if( (!is_pos_w(w.fun)) && test.method[t] == "inverse_w_hoeffding" )  # run only for positive distributions 
          next
        if(test.method[t] == 'bootstrap')
        {
          if(w.fun == 'huji')
            next
          # Run bootstrap even for truncation! (although the results are wrong here)
          #          if((w.fun[s] %in% c('truncation', 'Hyperplane_Truncation')) & (!exchange.type[s]))
          #            next # can't run bootstrap because w can be zero, unless we assume exchangability !!! 
        }      
        if(((test.method[t] %in% c('fast-bootstrap') || (test.stat[t] == "hoeffding")) ) && !(dependence.type %in% c('Gaussian', 'Clayton', 'Gumbel', 'LD')))
          next # try comparison with naive settings only for Gaussian case and other cases in Table 1 
        cur.test.method <- test.method[t]
        prms$importance.sampling.dist <- test.params[t]
        if(test.method[t] == "fast-bootstrap")
        {
          prms$fast.bootstrap <- 1
          cur.test.method <- 'bootstrap'
        }
        if(test.stat[t] == "hoeffding")
          prms$naive.expectation <- 1
        start.time <- Sys.time()
        #if(i==24 &&test.type[t]=='uniform_importance_sampling')
        #browser()
        # new! run sequencial tests and stop early for bootstrap/permutations
        if(prms$sequential.stopping)
        {
          #browser()
          block.size <- 100  # must divide prms$B !!! 
          prms$B <- ceil(prms$B / block.size) * block.size  # set a multiplication of block size 
          # prms$B <- block.size # wtf?
          prms$gamma <- 0.0001  # chance that we miss being below/above alpha at an early stop - 1/10000 for each block
          #          prms$z_gamma <- abs(qnorm(prms$gamma/2))
          cur.pvalue <- 0
          stop.flag <- 0
          for(b in 1:(prms$B/block.size))  # run block.size permutations each time 
          {
            cur.test.results <- TIBS(biased.data, w.fun, prms, cur.test.method, test.stat[t]) # we get to rcpp version inside TIBS if needed
            cur.pvalue <- cur.test.results$Pvalue + cur.pvalue
            
            cur.conf.int <- binom.confint(cur.pvalue*block.size, b*block.size, conf.level = 1-prms$gamma, method = "wilson")
            
            if((prms$alpha < cur.conf.int$lower) ||  (prms$alpha > cur.conf.int$upper))   # here stop early
            {
              print(paste0('early stopping after ', b*block.size, ' ', cur.test.method, ' out of ', prms$B, ' saving: ', round(100*(1.0-(b*block.size)/prms$B), 1), '% of work!'))
              stop.flag <- 1
            }
            if(stop.flag || (b == prms$B/block.size)) # reached last value
            {
              test.results <- c()
              test.results$Pvalue <- cur.pvalue / b
              break
            }
          }
        } else   # run full test: just simulate all B permutations 
        {
          #browser()
          test.results <- TIBS(biased.data, w.fun, prms, cur.test.method, test.stat[t])  # get into rcpp inside TIBS if needed
        }
        test.time[i.prm, t, i] <- difftime(Sys.time(), start.time, units='secs')
        test.pvalue[i.prm, t, i] <- test.results$Pvalue
        test.true.stat[i.prm, t, i] <- test.results$TrueT
        if('statistics.under.null' %in% names(test.results)) # not all tests have null statistics 
          test.null.stat[i.prm, t, i, ] <- test.results$statistics.under.null
        # optional printing:
        print(paste0("Dep: ", dependence.type, ', rho=', prms$rho, ', w: ', w.fun.str, '. Test: ', cur.test.method, " ", test.stat[t], 
                     ' i=', i, ' of ', prms$iterations, '. Pval=', round(test.pvalue[i.prm, t, i], 3), '. Time(sec.)=', round(test.time[i.prm, t, i], 2)))
      } # end loop on tests 
    }  # end loop on iterations (parallel collection). Simulation and testing for one parameter (i.prm)
    
    prms$title <- as.integer(dependence.type != 'UniformStrip') # s>1) 
    if(prms$plot.flag) # plot example 
      PlotBiasedData(dependence.type, biased.data, prms)
    
    # New: save intermediate results also for latex format 
    # Compute power (this is also type-1-error alpha under the null)
    test.power <- apply(test.pvalue < prms$alpha, 1, rowMeans) 
    # Save results in one file per dataset 
    test.output <- t(cbind(test.power, colSums(rowSums(test.time, dims=2))))
    colnames(test.output) <- test.method
    rownames(test.output) <- c(prms.rho, 'time') # [[s]]
    test.output <- test.output[, test.time[1,,1]>=0, drop = FALSE] # take only relevant tests 
    
    if(prms$run.flag!=0)
    {
      print('save results:') # save partial results for cases of script crashing
      if(i.prm < num.prms) # intermediate loops 
      {
        save(test.pvalue, test.time, prms.rho, prms, file=paste0(output.file, '.partial.Rdata'))
        print(xtable(test.output[c(1:i.prm, i.prm+2),], type = "latex", digits=3), 
              file = paste0(output.file, '.partial.tex'), size="\\tiny") # save in latex format
      } else #    if(i.prm == num.prms) # final loop 
      {
        save(test.pvalue, test.time, test.power, prms.rho, prms, file=paste0(output.file, '.Rdata'))  
        if(dim(test.output)[2]>0)
          print(xtable(test.output, type = "latex", digits=3), 
                file = paste0(output.file, '.tex'), size="\\tiny") # save in latex format 
      }
    }    
  } # end simulation and testing for one dependency type (loop over i.prm)

  return(list(test.power=test.power, test.output=test.output, test.pvalue=test.pvalue, 
              test.time=test.time, test.true.stat=test.true.stat, test.null.stat=test.null.stat))
  
} # end function

