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

print("set R studio")
isRStudio <- Sys.getenv("RSTUDIO") == "1" # check if we run interactively or inside a script
print("set run flag")
run.flag <- isRStudio # 1: run simulations inside R. -1: run simulations from outside command line.  0: load simulations results from file if they're available
if(isRStudio)
{
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # get path 
  path = getwd()
}
run.flag <- 1 # temp, just generate running scripts 
args=commandArgs(trailingOnly = TRUE)
print("Compile C docode ")
Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 

print("Compiled C docode ")

source('simulate_and_test.R')
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
source('import_samp.R')
source('Tsai_test.R')

print("Included sources ")
isRStudio <- Sys.getenv("RSTUDIO") == "1" # check if we run interactively or inside a script
run.flag <- 1 # set again 
#print(paste0("Now Rstudio=", isRStudio))

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
                             'LogNormal', 'sum', TRUE, TRUE,  list(c(0, 0.2)),  # added also a signal 
                             'Gaussian', 'gaussian', TRUE, TRUE, list(seq(-0.9, 0.9, 0.1)) ), 5, 9)) # -0.9 - 0.9 replace by CLmix / non-monotone and centered at zero 
#  'Gaussian','exponent_minus_sum_abs', TRUE, TRUE, # not needed (w(x,y)=w(x)*w(y), not interesting)


#print(paste0("Again Rstudio=", isRStudio))

dependence.type <- run.params.mat[,1]
w.fun <- run.params.mat[,2]
monotone.type <- run.params.mat[,3]
exchange.type <- run.params.mat[,4]
prms.rho <- run.params.mat[,5]


#print(paste0("Again2 Rstudio=", isRStudio))
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


#################################################################################
# Official parameters for long run:
test.stat <- c("adjusted_w_hoeffding", "inverse_w_hoeffding", "tsai") # possible test statistics # "hoeffding", , "tsai", "minP2" "adjusted_w_hoeffding", 
test.method <- c("permutationsIS", "permutationsMCMC", "bootstrap", "tsai") # possible methods for computing the test statistic "fast-bootstrap", "bootstrap",  
IS.methods <- c("tsai", "KouMcculough.w", "uniform", "monotone.w", "monotone.grid.w", "match.w") #  different methods for importance sampling of permutations
prms.file <- "sim/prms.sim"
run.script.file <- "run.all.sim.sh"
#################################################################################
# Temp parameters for experimentation
## test.stat <- c("adjusted_w_hoeffding")
## test.method <- c("permutationsMCMC")

#################################################################################


#print(paste0("Again3 Rstudio=", isRStudio))
##test.type <- c("uniform_importance_sampling_inverse_weighting", "uniform_importance_sampling", 'match_importance_sampling', 'monotone_importance_sampling')
# Official tests:
#    test.type <- c('permutations',  'permutations_inverse_weighting', 'bootstrap', 'bootstrap_inverse_weighting',  'tsai', 
#                   'monotone_importance_sampling', 'uniform_importance_sampling', 'match_importance_sampling', 'uniform_importance_sampling_inverse_weighting')  # official testing
#c( 'permutations','permutations_inverse_weighting', # everything except minP2
#  #'uniform_importance_sampling',
#  'uniform_importance_sampling_inverse_weighting',
#  'bootstrap', 
#  'bootstrap_inverse_weighting', 
#  'min_P2', 'tsai')
num.sim <- length(dependence.type)
# num.tests <- length(test.type)

print("Setting parameters")
run.rho <- prms.rho
  
if(isRStudio == 1)
{
  iterations = 500 # 00 # official: 500
  B = 1000 # 0 # official:  1000
  sample.size = 100 #  official:  100
  run.dep <- c(9) #  official: 1:9 # c(8:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
  
} else  # run from command line 
{
  print("inside else!!")
  args <- commandArgs(TRUE)
  
  print(as.integer(args[1]))
  run.dep <- as.integer(args[1]) #  c(7) # 2:num.sim) # 2 is only Gaussians (to compare to minP2 power) # 1 # Loop on different dependency types 
  print("set arg1")
  
  iterations = as.integer(args[2])
  print("set arg2")
  B = as.integer(args[3])
  print("set arg3")
  sample.size = as.integer(args[4])
  
  run.rho[[run.dep]] = as.numeric(args[5]) # new: run only on a single rho 
  print(paste0("run.dep:", run.dep, " iters:", iterations))
} # 4  # 10 for minP2 which is very slow  # 00 # 500  # Number of simulated dataset. Shared by all simulations

print("Finished Setting parameters")

#test.type<-c( 'permutations','permutations_inverse_weighting',
#              #'uniform_importance_sampling',
#              'uniform_importance_sampling_inverse_weighting',
#              'bootstrap', 
#              'bootstrap_inverse_weighting')

if(!isempty(intersect(run.dep, c(3,4,5,6,7)))) # %in% )
  library(copula) # needed for most simulations 

overall.start.time <- Sys.time()

run_str <- c()
ctr=1
for(s in run.dep) # Run all on the farm  
{
  for(n in c(sample.size)) #seq(250, 400, 50))
  {
    prms = list(B=B, sample.size=n, iterations=iterations, plot.flag=0, alpha=0.05, sequential.stopping=0, # pilot study 
                use.cpp=1, keep.all=0, perturb.grid=1, simulate.once=0, new.bootstrap=1, diagnostic.plot=0, 
                IS.methods=IS.methods, include.ID=1, run.sim=0, run.flag=1) # , sample.by.bootstrap=1) # set running parameters here ! 
    #    if(run.flag != 1)
    #      prms.rho[[s]] = as.numeric(args[4]) # temp for loading from user 
    print(paste0("s=", s))
    print(paste0("rho=", prms.rho[[s]]))
    # Call function. # run simulations function 
    print(paste("n=", prms$sample.size))
    if(const.seed)
      prms$seed <- 1234567890 # 4524553
    
    # New: set applicible tests: 
    test.comb <- GetTestCombinations(prms, w.fun[[s]], dependence.type[[s]], test.stat, test.method)
    
    if(s < 8) # get rid of all IS methods 
    {
      test.comb <- test.comb[test.comb[,1] != "permutationsIS",]
      test.comb <- test.comb[test.comb[,2] != "inverse_w_hoeffding",]
    }
    num.tests <- dim(test.comb)[1] # can change with s !! 

    
    save(prms, test.comb, file=paste0(prms.file, '.', s, '.Rdata'))
    
    
    # new: separate into different rho values: 
    if(run.flag==1)
      T.OUT <- simulate_and_test(dependence.type[[s]], run.rho[[s]], w.fun[[s]], test.comb, prms) # run all tests on one type of simulatee data 
    else  # prepare job strungs 
    {
      # New: just create jobs strings 
      for(k in c(1:length(run.rho[[s]])))  # run each parameter separately:
      {
        run_str[ctr] <- paste0("Rscript run_simulations.R ", run.dep, " ", iterations, " ",  B, " ", sample.size, " ",  run.rho[[s]][k],  " > out/run.sim.s.", s, ".rho.",  run.rho[[s]][k], ".out ", " &")
                          
##        run_str[ctr] <- paste0("Rscript simulate_and_test.R ", dependence.type[[s]], " ", run.rho[[s]][k], " ", w.fun[[s]], " ", "\"c()\"",  " ", 
##                               paste0(prms.file, '.', s, '.Rdata'), " > out/run.sim.s.", s, ".rho.",  prms.rho[[s]][k], ".out ", " &")  # test.comb,
        print(run_str[ctr])  # save all run strungs into a script file 
        ctr = ctr+1
      }
      
    }
  }
} # end loop on dependency types

if(run.flag==0)# Create running script
{
  run.script.ptr <- file(run.script.file)
  writeLines(c("#!/bin/sh", run_str), run.script.file)
  close(run.script.ptr)
}

run.flag=0
if(run.flag && isRStudio)  # plot results in interactive mode
{
  test.legend <- apply(cbind(test.comb[,c(1,3)], T.OUT$test.power), 1, paste0, collapse=" ")#    paste(test.type, as.character(T.OUT$test.power))
  col.vec <- c("orange", "cyan", "navy", "dodgerblue3", "pink", "green", "black", "blue", "black", "green", "orange", "gray", "pink",  "purple", "brown", "yellow", "magenta", "darkgreen", "gold")
  i=1
  plot(T.OUT$test.true.stat[1,i,]- rowMeans(T.OUT$test.null.stat[1,i,,]), col=col.vec[i], pch=20, main='differences')
  for(i in 2:dim(test.comb)[1])
    points(T.OUT$test.true.stat[1,i,]- rowMeans(T.OUT$test.null.stat[1,i,,]), col=col.vec[i], pch=20)
  legend(0, 200, test.legend, lwd=c(2,2), col=col.vec[1:num.tests], y.intersp=0.8, cex=0.6)
  # points(T.OUT$test.true.stat[1,,]- rowMedians(T.OUT$test.null.stat[1,1,,]), col="blue")
  
  
  ##jpeg(paste0("../figs/check_valid_n_", n, "_B_", prms$B, "_dep_", dependence.type[s], "_w_",  w.fun[s], 
  ##            "_perturb_grid_", prms$perturb.grid, ".jpg"), width = 400, height = 400)
  plot(c(0, prms$iterations), c(0,1), col="red", type="l", 
       main=TeX(paste0("Tests Cumulative Pvalues, $n=", n, ", \\alpha =", prms$alpha)), # , "$ pert=", prms$perturb.grid)), 
       xlab="Rank", ylab="Pvalue")
  valid.tests <- rep(0, num.tests)
  for(i in c(1:num.tests)) # c(1:4, 6:num.tests))
  {
    if(max(T.OUT$test.pvalue[1,i,])> -1)
    {
      points(sort(T.OUT$test.pvalue[1,i,]), col=col.vec[i], pch=20)
      valid.tests[i] <- 1
    }
  }
  for(i in 1:num.tests)
  {
    test.legend[i] <- str_remove(test.legend[i], "permutations")
    test.legend[i] <- str_replace(test.legend[i], ".w", "")
  }
  grid(NULL,NULL, lwd=1)
  legend(prms$iterations*0.75, 0.18, test.legend[c(1:num.tests)], lwd=c(2,2), col=col.vec[c(1:num.tests)], 
         y.intersp=0.8, cex=0.6, box.lwd = 0,box.col = "white",bg = "white")
#  legend(prms$iterations*0.75, 0.18, test.legend[c(1:4,6:num.tests)], lwd=c(2,2), col=col.vec[c(1:4,6:num.tests)], 
#         y.intersp=0.8, cex=0.6, box.lwd = 0,box.col = "white",bg = "white")
  ##  dev.off()
}

#library(matrixStats)
#plot(T.OUT$test.true.stat[1,,], rowMeans(T.OUT$test.null.stat[1,1,,]))
#plot(T.OUT$test.true.stat[1,,], rowMedians(T.OUT$test.null.stat[1,1,,]), col="blue")
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





# Temp: samll test of product of Gaussian distribution 
test.gaussian <- 0
if(test.gaussian)
{
  rho <- 0.95
  rho2 <- 0.7
  Sigma.p <- matrix(c(1, rho, rho, 1), 2, 2)
  Sigma.m <- matrix(c(1, -rho2, -rho2, 1), 2, 2)
  mu <- c(0,0)
  
  
  
  x <- seq(-1,1,0.05)
  y <- seq(-1,1,0.05)
  n <- length(x)
  xy.mat <- matrix(0, n, n)
  for(i in 1:n)
    for(j in 1:n)
      xy.mat[i,j] <- dmvnorm(c(x[i], y[j]), mu, Sigma.m) * dmvnorm(c(x[i], y[j]), mu, Sigma.p)
  
  persp(x, y, xy.mat)
  
  
  
  rho <- 0.95
  n <- 500
  g.prms <- c()
  g.prms$rho <- rho
  g.prms$w.rho <- -rho
  g.prms$w.max <- 1
  xy <- SimulateBiasedSample(n, "Gaussian", "gaussian", g.prms)
  g.prms$rho <- -rho
  g.prms$w.rho <- rho
  xy2 <- SimulateBiasedSample(n, "Gaussian", "gaussian", g.prms)
  
  
  
  plot(xy$data[,1], xy$data[,2])
  cor(xy$data[,1], xy$data[,2])
  points(xy2$data[,1], xy2$data[,2], col="red")
  cor(xy2$data[,1], xy2$data[,2])
  
}


all.test.comb <- c()
for(s in run.dep)
{
  load(paste0("sim/prms.sim.", s, ".Rdata"))
  all.test.comb[[s]] <- rbind(c(dependence.type[[s]], w.fun[[s]], ""),  test.comb)
}
