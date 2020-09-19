rm(list=ls())
gc()

library(stringr)
library(foreach)
library(doSNOW)
#library(parallel)
library(doParallel)
library(gdata)
library(mvtnorm)
library(ggplot2)  
library(pracma)
library(matrixStats)
library(PerMallows) # distance between permutations 
library(Matrix)


isRStudio <- Sys.getenv("RSTUDIO") == "1" # check if we run interactively or inside a script
run.flag <- isRStudio # 1: run simulations inside R. -1: run simulations from outside command line.  0: load simulations results from file if they're available
if(isRStudio)
{
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # get path 
  path = getwd()
}
args=commandArgs(trailingOnly = TRUE)
Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 

source('simulate_and_test.R')
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
source('import_samp.R')
source('Tsai_test.R')

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

const.seed <- 1 # set constant seed 

# Vectors with different dependency settings : dependence-type, w, monotone, exchangable, rho
run.params.mat <- t(matrix(c('UniformStrip', 'truncation', TRUE, TRUE, list(0.3), 
                             'Gaussian', 'truncation', TRUE, TRUE, list(seq(-0.9,0.9,0.1)),
                             'Clayton','truncation', TRUE, TRUE, list(0.5),
                             'Gumbel', 'truncation', TRUE, TRUE, list(1.6),
                             'LD', 'truncation', TRUE, FALSE, list(c(0, 0.4)),
                             'nonmonotone_nonexchangeable', 'truncation', FALSE, FALSE, list(seq(-0.9, 0.9, 0.1)),
                             'CLmix','truncation', FALSE, TRUE, list(0.5), 
                             'LogNormal', 'sum', TRUE, TRUE,  list(c(0)),
                             'Gaussian', 'gaussian', TRUE, TRUE, list(seq(-0.9, 0.9, 0.1)) ), 5, 9)) # replace by CLmix / non-monotone and centered at zero 
#  'Gaussian','exponent_minus_sum_abs', TRUE, TRUE, # not needed (w(x,y)=w(x)*w(y), not interesting)

run.params.mat

dependence.type <- run.params.mat[,1]
w.fun <- run.params.mat[,2]
monotone.type <- run.params.mat[,3]
exchange.type <- run.params.mat[,4]
prms.rho <- run.params.mat[,5]

#dependence.type <- c('UniformStrip', 'Gaussian', 'Clayton', 'Gumbel', 
#                     'LD', 'nonmonotone_nonexchangeable', 'CLmix', 'Gaussian',
#                     'LogNormal', 'Gaussian') # The last one is log-normal 
#w.fun <- c('truncation', 'truncation', 'truncation', 'truncation', 
#           'truncation', 'truncation', 'truncation', 
#           'exponent_minus_sum_abs', 'sum', 'gaussian') # not good that we have only one simulation with positive W. Should add X+Y?
#monotone.type <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE) # is monotone
#exchange.type <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE) # is exchangeable
# sample.size <- c(500, 100, 100, 100, 100, 100, 100, 100) # set all sample.sizes to 100 
#prms.rho <- list(0.3, seq(-0.9, 0.9, 0.1), 0.5, 1.6, c(0, 0.4),
#                 seq(-0.9, 0.9, 0.1), 0.5, c(0), c(0), c(0))#seq(-0.9, 0.9, 0.1), 
#c(0)) # Parameters for each sampling type 

### test.type<- c('permutations','bootstrap') #c( 'permutations','permutations_inverse_weighting',
## test.type <- c('uniform_importance_sampling', 'uniform_importance_sampling_inverse_weighting') #c( 'permutations','permutations_inverse_weighting',

test.type <- c("uniform_importance_sampling_inverse_weighting", "uniform_importance_sampling", 'match_importance_sampling', 'monotone_importance_sampling')
# Official tests:
#    test.type <- c('permutations',  'permutations_inverse_weighting', 'bootstrap', 'bootstrap_inverse_weighting',  'tsai', 
#                   'monotone_importance_sampling', 'uniform_importance_sampling', 'match_importance_sampling', 'uniform_importance_sampling_inverse_weighting')  # official testing
#c( 'permutations','permutations_inverse_weighting', # everything except minP2
#  #'uniform_importance_sampling',
#  'uniform_importance_sampling_inverse_weighting',
#  'bootstrap', 
#  'bootstrap_inverse_weighting', 
#  'min_P2', 'Tsai')
num.sim <- length(dependence.type)
num.tests <- length(test.type)

if(run.flag == 1)
{
  iterations = 200 # official: 500
  B = 100 # official:  1000
  sample.size = 301 #  official:  100
  run.dep <- c(8) #  official: 1:9 # c(8:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
  
} else  # run from command line 
{
  run.dep <- as.integer(args[1]) #  c(7) # 2:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
  iterations = as.integer(args[2])
  B = as.integer(args[3])
  sample.size = as.integer(args[4])
  print(paste0("run.dep:", run.dep, " iters:", iterations))
  
} # 4  # 10 for minP2 which is very slow  # 00 # 500  # Number of simulated dataset. Shared by all simulations

#test.type<-c( 'permutations','permutations_inverse_weighting',
#              #'uniform_importance_sampling',
#              'uniform_importance_sampling_inverse_weighting',
#              'bootstrap', 
#              'bootstrap_inverse_weighting')

if(!isempty(intersect(run.dep, c(3,4,5,6,7)))) # %in% )
  library(copula) # needed for most simulations 

overall.start.time <- Sys.time()

for(s in run.dep) # Run all on the farm  
{
  for(n in c(sample.size)) #seq(250, 400, 50))
  {
    prms = list(B=B, sample.size=n, iterations=iterations, plot.flag=0, alpha=0.05, sequential.stopping=0, # pilot study 
                use.cpp=0, keep.all=0, perturb.grid=1, simulate.once=0, new.bootstrap=1) # , sample.by.bootstrap=1) # set running parameters here ! 
    #    if(run.flag != 1)
    #      prms.rho[[s]] = as.numeric(args[4]) # temp for loading from user 
    print(paste0("s=", s))
    print(paste0("rho=", prms.rho[[s]]))
    # Call function. # run simulations function 
    print(paste("n=", prms$sample.size))
    if(const.seed)
      prms$seed <- 1141248 # 4524553
    T.OUT <- simulate_and_test(dependence.type[[s]], prms.rho[[s]], w.fun[[s]], test.type, prms) # run all tests 
    
    # New: just create jobs strings 
    for(k in c(1:length(prms.rho[[s]])))  # run each parameter separately:
      run_str <- paste0("T.OUT <- simulate_and_test(", dependence.type[[s]], ", ", prms.rho[[s]][k], ", ", w.fun[[s]], ", ", prms)
  }
} # end loop on dependency types



if(isRStudio)  # plot results in interactive mode
{
  test.legend <- paste(test.type, as.character(T.OUT$test.power))
  col.vec <- c("blue", "black", "green", "orange", "gray", "pink", "yellow", "purple", "cyan")
  i=1
  plot(T.OUT$test.stat[1,i,]- rowMeans(T.OUT$test.null.stat[1,i,,]), col=col.vec[i], pch=20, main='differences')
  for(i in 2:length(test.type))
    points(T.OUT$test.stat[1,i,]- rowMeans(T.OUT$test.null.stat[1,i,,]), col=col.vec[i], pch=20)
  legend(0, 200, test.legend, lwd=c(2,2), col=col.vec[1:length(test.type)], y.intersp=0.8, cex=0.6)
  # points(T.OUT$test.stat[1,,]- rowMedians(T.OUT$test.null.stat[1,1,,]), col="blue")
  
  
  jpeg(paste0("../figs/check_valid_n_", n, "_B_", prms$B, "_dep_", dependence.type[s], "_w_",  w.fun[s], 
              "_perturb_grid_", prms$perturb.grid, ".jpg"), width = 400, height = 400)
  plot(c(0, prms$iterations), c(0,1), col="red", type="l", 
       main=paste0("Tests pvals and power, n=", n, ", alpha=", prms$alpha, " pert=", prms$perturb.grid))
  valid.tests <- rep(0, num.tests)
  for(i in 1:num.tests)
  {
    if(max(T.OUT$test.pvalue[1,i,])> -1)
    {
      points(sort(T.OUT$test.pvalue[1,i,]), col=col.vec[i], pch=20)
      valid.tests[i] <- 1
    }
  }
  legend(prms$iterations*0.45, 0.25, test.legend[which(valid.tests>0)], lwd=c(2,2), col=col.vec[which(valid.tests>0)], y.intersp=0.8, cex=0.6)
  dev.off()
}

#library(matrixStats)
#plot(T.OUT$test.stat[1,,], rowMeans(T.OUT$test.null.stat[1,1,,]))
#plot(T.OUT$test.stat[1,,], rowMedians(T.OUT$test.null.stat[1,1,,]), col="blue")
#lines(c(15000, 50000), c(15000, 50000), col="red")


# small mcmc test: 
run.mcmc <- 0
if(run.mcmc)
{
  w = c(2,3,5,10, 1, 1, 0.1, 55)
  k = length(w)
  z = sum(w)
  p_w = w / sum(w)
  B = 1000000
  x = sample(k, B, prob=p_w, replace=TRUE)
  u = sample(k, B, replace=TRUE)
  
  w[x]
  z_hat = B*k / (sum(1 / w[x]))
  z_hat_u = (k*sum(w[u])) / (B)
  
  print(z)
  print(z_hat)
  print(z_hat_u)
  
}







overall.time <-  difftime(Sys.time() , overall.start.time, units="secs") 
print(paste0("Overall simulations time (no min.p), B= ", prms$B, ", iters=", prms$iterations, ":"))
print(overall.time)
