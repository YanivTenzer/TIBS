library(latex2exp)
library(xtable)
library(binom)  # for binomial confidence intervals 


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
                              test.type=c('tsai', 'minP2', 'permutations', 'bootstrap', 'fast-bootstrap', 'naive-bootstrap', 'naive-permutations'), 
                              B=100, sample.size=100, iterations=50, plot.flag=0, alpha=0.05, sequential.stopping=0)
{
  print(paste0("rho.inside=", prms.rho))
  prms = list(B = B)
  run.flag = 1
  num.tests <- length(test.type)
  
  output.file <- paste0('results/', dependence.type, '_all_tests_results_B_', 
                        B, '_iters_', iterations, '_n_', sample.size, '_rho_', paste(prms.rho, collapse = '_')) # set output file name
  num.prms <- length(prms.rho)
  if((!run.flag) & file.exists(paste0(output.file, '.Rdata')))
    load(file=paste0(output.file, '.Rdata'))
  else
    test.pvalue <- test.time <- array(-1, c(num.prms, num.tests, iterations)) 
  
  for(i.prm in 1:num.prms) # loop on different simulation parameters - should plot for each of them? 
  {
    prms$rho = prms.rho[i.prm] # [[s]]
    prms$W.max <- 1 # temp: 1 for all weights
    prms$minp.eps <- "in"
    prms$keep.all <- 0 # set to 1 for plotting
    
    ## Parallel on multiple cores 
    ##    results.table<- foreach(i=seq(iterations), .combine=rbind) %dopar%{ 
    ##      iteration_result <- matrix(0, 4, B+1)
    for(i in 1:iterations) # run simulations multiple times 
    { 
      biased.data <- SimulateBiasedSample(sample.size, dependence.type, w.fun, prms) 
      all.data <- biased.data$all.data
      biased.data <- biased.data$data 
      if((!run.flag) & file.exists(paste0(output.file, '.Rdata'))) # simulate only once 
        break
      for(t in c(1:length(test.type)))  # ALL TESTs !    #  c(1,2,6,7)) # loop on statistical tests . Last test (7) minP2 is very slow 
      {
        prms$fast.bootstrap <- 0 # method for computing null expectations 
        prms$naive.expectation <- 0 # Check test with standard expectations
        # Skip irrelevant tests 
        if((!(w.fun %in% c('truncation', 'Hyperplane_Truncation'))) & 
           (test.type[t] %in% c("tsai", 'minP2')))
          next  # these tests run only for truncation 
        if(test.type[t] == 'bootstrap')
        {
          if(w.fun == 'huji')
            next
          # Run bootstrap even for truncation! (although the results are wrong here)
          #          if((w.fun[s] %in% c('truncation', 'Hyperplane_Truncation')) & (!exchange.type[s]))
          #            next # can't run bootstrap because w can be zero, unless we assume exchangability !!! 
        }      
        if((t %in% c(5:7)) & !(dependence.type %in% c('Gaussian', 'Clayton', 'Gumbel', 'LD')))
          next # try comparison with naive settings only for Gaussian case and other cases in Table 1 
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

        # new! run sequencial tests and stop early for bootstrap/permutations
        if(sequential.stopping)
        {
          block.size <- 50
          prms$B <- block.size
          prms$gamma <- 0.0001  # chance that we miss being below/above alpha at an early stop - 1/10000 for each block
#          prms$z_gamma <- abs(qnorm(prms$gamma/2))
          cur.pvalue <- 0
          stop.flag <- 0
          for(b in 1:(B/block.size))  # run block.size permutations each time 
          {
            cur.test.results<-TIBS(biased.data, w.fun, cur.test.type, prms)
            cur.pvalue <- cur.test.results$Pvalue + cur.pvalue
            
            cur.conf.int <- binom.confint(cur.pvalue*block.size, b*block.size, conf.level = 1-prms$gamma, method = "wilson")
            
            # test if we should stop! 
#            print(paste('cur.pval=', round(cur.pvalue/b, 3), 'conf.int=(', round(cur.conf.int$lower, 3), ', ', round(cur.conf.int$upper, 3), ') alpha=', alpha))
            if((alpha < cur.conf.int$lower) ||  (alpha > cur.conf.int$upper))   # here stop early !!! 
            {
              print(paste0('early stopping after ', b*block.size, ' ', cur.test.type, ' out of ', B, ' saving: ', round(100*(1-b*block.size/B), 1), '% of work!'))
              stop.flag <- 1
            }
            if(stop.flag || (b == B/block.size)) # reached last value!
            {
              test.results <- c()
              test.results$Pvalue <- cur.pvalue / b
              break
            }
          }
          
        } else   # run full test: just simulate all B permutations 
            test.results<-TIBS(biased.data, w.fun, cur.test.type, prms)
        test.time[i.prm, t, i] <- difftime(Sys.time(), start.time, units='secs')
        test.pvalue[i.prm, t, i] <- test.results$Pvalue
        print(paste0(dependence.type, ', rho=', prms$rho, '. Test: ', test.type[t], 
                     ' i=', i, ' of ', iterations, '. Pval=', round(test.pvalue[i.prm, t, i], 3), '. Time (sec.)=', round(test.time[i.prm, t, i], 2)))
      } # end loop on tests 
    }  # end loop on iterations (parallel collection). Simulation and testing for one parameter (i.prm)
    prms$title <- as.integer(dependence.type != 'UniformStrip') # s>1) 
    if(plot.flag) # plot example 
      PlotBiasedData(dependence.type, biased.data, prms)
    
    # New: save intermediate results also for latex format 
    
    # Compute power (this is also type-1-error alpha under the null)
    test.power <- apply(test.pvalue<alpha, 1, rowMeans) 
    # Save results in one file per dataset 
    test.output <- t(cbind(test.power, colSums(rowSums(test.time, dims=2))))
    colnames(test.output) <- test.type
    rownames(test.output) <- c(prms.rho, 'time') # [[s]]
    test.output <- test.output[, test.time[1,,1]>=0, drop = FALSE] # take only relevant tests 
    
    print('save results:') # save partial results for cases of script crashing
    if(i.prm < num.prms) # intermediate loops 
    {
      save(test.pvalue, test.time, prms.rho, sample.size, B, iterations, file=paste0(output.file, '.partial.Rdata'))
      print(xtable(test.output[c(1:i.prm, i.prm+2),], type = "latex", digits=3), 
            file = paste0(output.file, '.partial.tex'), size="\\tiny") # save in latex format
    } else #    if(i.prm == num.prms) # final loop 
    {
      save(test.pvalue, test.time, test.power, prms.rho, sample.size, B, iterations, file=paste0(output.file, '.Rdata'))  
      print(xtable(test.output, type = "latex", digits=3), 
            file = paste0(output.file, '.tex'), size="\\tiny") # save in latex format 
    }
  } # end simulation and testing for one dependency type (loop over i.prm)
  
  return(list(test.power=test.power, test.output=test.output))
  
} # end function

