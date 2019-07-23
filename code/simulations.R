cat("\014")
# path = 'D:/cond_ind_2019' # change path
path = 'C:\\Users\\Or Zuk\\Dropbox\\BiasedSampling\\Code'
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
library(ggplot2) # need to install it 

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

hyperplane.prms<-c(-1,1,0)
bias.params<-list(hyperplane.prms)
names(bias.params)<-c("hyperplane.prms")
bias_method = 'Hyperplane_Truncation'  
Estimate_Marginals = 1
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
sample.size <- c(500, 100, 100, 100, 500, 500, 100, 500) # should all sample.sizes be 500 ?
prms.rho <- list(0.3, seq(-0.9, 0.9, 0.1), 0.5, 1.6, c(0, 0.4),
                 seq(-0.9, 0.9, 0.1), 0.5, seq(-0.9, 0.9, 0.1)) # Parameters for each sampling type 

test.type <- c('bootstrap', 'permutations', 'tsai', 'minP2') # different tests to run 
iterations = 10 # 500  # Number of simulated dataset. Shared by all simulations
B = 10^2 # 10^4 # number of permtuations or bootstrap samples. Shared by all simulations 
num_sim <- length(dependence.type)

for(s in 1:num_sim)  # 1 # Loop on different dependency types 
{
  output.file <- paste0('results/', dependence.type[s], '_4tests_results_B_', 
                        B, '_', 'iters_', iterations, '.Rdata') # set output file name
  print(paste0("s=", s))
  if((!run.flag) & file.exists(output.file))
  {  # try loading
    load(file=output.file)
    next
  }
  set.seed(1)
  num.prms <- length(prms.rho[[s]])
  test.pvalue <- test.time <- array(-1, c(num.prms, 4, iterations)) 
  for(iter in 1:num.prms) # loop on different simulation parameters 
  {
    prms$rho = prms.rho[[s]][iter]
    prms$W_max <- 1 # temp: 1 for all weights
    
    ## Parallel on multiple cores 
    ##    results.table<- foreach(i=seq(iterations), .combine=rbind) %dopar%{ 
    ##      iteration_result <- matrix(0, 4, B+1)
    for(i in 1:iterations){
      biased.data<-SimulateBiasedSample(sample.size[s], dependence.type[s], bias.type[s], prms) 
      for(j in 1:3) # loop on statistical tests . Last test (4) minP2 is very slow 
      {
        # Skip irrelevant tests 
        if((!(bias.type[s] %in% c('truncation', 'Hyperplane_Truncation'))) & (test.type[j] %in% c("tsai", 'minP2')))
          next  # these tests run only for truncation 
        if(test.type[j] == 'bootstrap')
        {
          if(bias.type[s] == 'huji')
            next
          if((bias.type[s] %in% c('truncation', 'Hyperplane_Truncation')) & (!exchange.type[s]))
            next # can't run bootstrap because w can be zero, unless we assume exchangability !!! 
        }      
        start.time <- Sys.time()
        test.results<-TIBS(biased.data, bias.type[s], bias.params,   # change to bias.type[s]
                           B, test.type[j], prms)
        test.time[iter, j, i] <- difftime(Sys.time(), start.time, units='secs')
        test.pvalue[iter, j, i] <- test.results$Pvalue
        print(paste0(dependence.type[s], ', rho=', prms$rho, ' test: ', test.type[j], 
                     ' iter=', i, ' of ', iterations, ' time (sec.)=', round(test.time[iter, j, i], 2)))
      } # end loop on tests 
    }  # end parallel collection. Simulation and testing for one parameter (iter)
  } # end simulation and testing for one dependency type (loop over iter)
  
  # Compute power (this is also type-1-error alpha under the null)
  test.power <- apply(test.pvalue<alpha, 1, rowMeans) 
  # Modify to save in ONE file: (need to change output_dir / file name to indicate dependence type)
  save(test.pvalue, test.time, test.power, prms.rho, sample.size, s, file=output.file)

      #  print(xtable(results.table, type = "latex"), file = paste0(output_dir ,'/Gaussian_rho_',toString(prms$rho), '_four_tests_results.tex')) # new! save in latex format 
  if(plot.flag) # plot example 
  {
    output.scatter.file <- paste0(path, '/../Figures/simulations/', 
                                  dependence.type[s], '.bmp')  # set output figure file name
    bmp(output.scatter.file, units="in", width=5, height=5, res=300)
    
    xy <- data.frame(biased.data)   # plot and save to image file 
    jpeg() # set output file name 
        ggplot(xy, aes(x=xy[,1], y=xy[,2])) + 
      geom_point() + 
      ggtitle(paste(as.character(dependence.type[s]), "(", expression(theta), "=", prms.rho[s], ")")) +
      xlab("X") + ylab("Y") +
      scale_y_continuous(breaks=seq(-2, 2, 2)) +
      scale_x_continuous(breaks=seq(-2, 2, 2)) +
      theme(
        plot.title = element_text(size=14, face="bold.italic", hjust=0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size = 14),
        axis.text.x = element_text(face="bold", size=12), 
        axis.text.y = element_text(face="bold", size=12)
      )
#    save(figure, file=paste0('../Figures/', dependence.type[s], '_scatter.jpeg'))
    dev.off()    
    
    if(s==1) # Show example with Kaplan Meyer marginal estimation  
    {
      marginals <-EstimateMarginals(biased.data, 'survival')
      plot(-marginals$CDFs[,1]$time, marginals$CDFs[,1]$surv, type='s',col='blue',xlab="t",ylab="KM(t)")
      lines(marginals$CDFs[,2]$time, 1-marginals$CDFs[,2]$surv, type='s',col='red')
      vec.marg <- Vectorize(FUN = UniformStripMarginal,vectorize.args = "t")
      t <- seq(0,1,len=200)
      F <- vec.marg(t,0.3)
      lines(t,F)
      
      # Show data sampled from marginals
    }
    
  } # end if plot.flag
  
} # end loop on dependency types 


