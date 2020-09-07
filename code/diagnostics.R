# Run diagnostic for weighted permutations tests 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # get path 
path = getwd()
setwd(path)
args=commandArgs(trailingOnly = TRUE)
Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 
source('simulate_and_test.R')
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
source('import_samp.R')
library(mvtnorm)
library(pracma)
library(matrixStats)
library(rcpp)
library(RcppArmadillo)


# 1. Check overdispersion of variance of p-vals calculatoin
B_vec <- seq(200, 2000, 200) # try different numbers of permutations
iters <- 25 # number of times to estimate for the same dataset: 
test.type <- 'bootstrap'
prms <- c()
rho.max <- -0.05 # small correlations 

n_vec <- c(100, 200, 500)
prms$sample.size <- n_vec[1]  # sample size 
prms$use.cpp <- 0

prms$sample.by.bootstrap <- 0 # enable new sampling method for general w 
w.fun <- 'sum'
dependence.type <- 'LogNormal'
pvals <- matrix(0, iters, length(B_vec))

prms$rho <- 0
biased.data <- SimulateBiasedSample(prms$sample.size, dependence.type, w.fun, prms) 
# FIX SO W>MAX ISN"T SET TO ! IN BIASED DATA. 
# THEN CHECK BOTH OPTIONS SUM AND TRUNCATION ON R AND CPP 

temp.data <- biased.data
temp.data$data <- biased.data$data+0.0000001
temp.data.sorted <- biased.data
temp.data.sorted$data[,1]= sort(temp.data.sorted$data[,1])
temp.data.sorted$data[,2]= sort(temp.data.sorted$data[,2])
if(!('w.max' %in% names(prms)))
  prms$w.max <- biased.data$w.max
prms$B <- 200
plot(biased.data$data[,1], biased.data$data[,2])
points(temp.data$data[,1], temp.data$data[,2], col='red')
points(temp.data.sorted$data[,1], temp.data.sorted$data[,2], col='blue')
prms$w.max <- max(prms$w.max, set_w_max_sample(biased.data$data, w.fun))

test.type = "permutations"
data <- biased.data$data
grid.points <- data
p0 <-  TIBS.steps(data, w.fun, c(), grid.points, c(),  prms)
prms$new.bootstrap = FALSE
p0.cpp <- TIBS_steps_rcpp(data, w.fun, matrix(1,1,1), grid.points, matrix(1,1,1), prms)
prms$new.bootstrap = TRUE
p0.new.cpp <- TIBS_steps_rcpp(data, w.fun, matrix(1,1,1), grid.points, matrix(1,1,1), prms)


p0$Statistic
p0.cpp$Statistic
p0.new.cpp$Statistic


#test.type = 'permutations'
p.zero <- TIBS(biased.data$data, w.fun, test.type, prms)
prms$new.bootstrap = FALSE
p.zero.cpp <- TIBS_rcpp(biased.data$data, w.fun, test.type, prms)
prms$new.bootstrap = TRUE
p.zero.cpp.new <- TIBS_rcpp(biased.data$data, w.fun, test.type, prms)

#points(biased.data$data[,1], biased.data$data[,2], col='green') # Problem: C code CHANGES the data !!!! 
#points(temp.data$data[,1], temp.data$data[,2], col='orange') # Problem: C code CHANGES the data !!!! 

p.zero$TrueT
p.zero.cpp$TrueT
p.zero.cpp.new$TrueT

p.zero$Pvalue
p.zero.cpp$Pvalue
p.zero.cpp.new$Pvalue

mean(p.zero$statistics.under.null)
mean(p.zero.cpp$statistics.under.null)
mean(p.zero.cpp.new$statistics.under.null)


# Next, run permutations test, and verify that it is valud
iters <- 100
p.val <- matrix(0, iters)
p.val.cpp <- matrix(0, iters)
prms$B=200

#for(t in c(1:2))
#{
#  for(iter in c(1:iters))
#  {
#    print(paste0("Run Iter ", iter))
#    biased.data <- SimulateBiasedSample(prms$sample.size, dependence.type, w.fun, prms) 
#    p.val[iter] <- TIBS(biased.data$data, w.fun, test.type, prms)$Pvalue
#    p.val.cpp[iter] <- TIBS_rcpp(biased.data$data, w.fun, test.type, prms)$Pvalue
#  }
#  plot(sort(p.val))
#  points(sort(p.val.cpp), col="red")#
#
#  sum(p.val < 0.05)  # power at alpha=0.05
#  sum(p.val.cpp < 0.05)
#}
  



test.type = c("permutations" , "permutations_inverse_weighting") # check scaling for both 
num.tests <- length(test.type)
B_vec <- round(logspace(2, 4, 10)) # round(logspace(2, 4, 10))
n_vec <- c(200) # , 200, 500, 1000) # small n to find memory leak fast 

iters <- 25 # to get reasonable variance 
prms$use.cpp=1
prms$perturb.grid <- 0 # CHECK IF THIS ADDS VARIANCE!!
prms$counts.flag <- 1 # New statistic for inverse weighting (use counts and not probabilities)
pvals <- matrix(0, iters, length(B_vec))
empirical.std <- matrix(0, num.tests, length(B_vec))
epsilon <- 0.0000000001 # tolerance 
for(n in n_vec[1])
{
  prms$rho <- rho.max * min(n_vec) / n  # change rho to keep power the same 
  prms$sample.size <- n
  biased.data <- SimulateBiasedSample(prms$sample.size, dependence.type, w.fun, prms) 
  prms$w.max <- max(biased.data$w.max, set_w_max_sample(biased.data$data, w.fun))
  prms$grid.points <- biased.data$data  # fix the grid for all sample sizes and iterations 
  prms$grid.points <- prms$grid.points + epsilon * matrix( rnorm(2*n), n, 2) # perturb   
  for(t in c(1:num.tests))
  {
    for(i in c(1:iters))  # iters can be much lower - we just need st.d. !!! 
    {
      print(paste0(test.type[t], " run i=", i, " out of ", iters))
      for(b in 1:length(B_vec))
      {
#        print("Start again")
        prms$B <- B_vec[b]
#        print("B=")
#        print(prms$B)
        start.time <- Sys.time()
#        print("took time again!")
        TR <- TIBS(biased.data$data, w.fun, test.type[t], prms) #  TIBS(biased.data, w.fun, cur.test.type, prms)
#        print("Finished TIBS")
#        print("PVAL=")
#        print(TR$Pvalue)
        pvals[i,b] <-TR$Pvalue #  TIBS(biased.data$data, w.fun, test.type[t], prms)$Pvalue #  TIBS(biased.data, w.fun, cur.test.type, prms)
        test.time <- difftime(Sys.time(), start.time, units='secs')
#        print("Total time:")
        print(test.time)
      }
#      print(paste0("Finished run iter i=", i, " out of ", iters))
    }
#    print("Pvals:")
#    print(pvals)
    true.pval <- mean(pvals)
#    print("True pval mean:")
#    print(true.pval)
        
#    print("Theo St.d.:")
    theoretical.std <- sqrt(true.pval * (1-true.pval) / B_vec)
#    print("Empiric:")
    empirical.std[t,] <- sqrt(colVars(pvals))
  }
#  jpeg(paste0("../figs/var_p_hat_dep_", dependence.type, "_w_", w.fun, "_n_", n, ".jpg"), width = 400, height = 400)
  plot(B_vec, log10(theoretical.std), type="l", col="green", lwd=5, xlab="B", ylab="log_{10}(std)", 
       ylim=log10(range(empirical.std, theoretical.std)),
       main=paste0(" std(p.hat) for n=", prms$sample.size, " p.val=", round(true.pval, 3) ))
  color.vec <- c("red", "blue", "orange")
  for(t in c(1:num.tests))
    lines(B_vec, log10(empirical.std[t,]), col=color.vec[t])
  
#  legend(max(B_vec)*0.7, max(log10(empirical.std)), c("theoretical", test.type), lwd=c(5,2,2), col=c("green", "red", "blue"))
#  legend(0, min(log10(theoretical.std)+0.3), c("theoretical", test.type), 
#         lwd=c(5,2,2), cex=0.7, col=c("green", "red", "blue", "orange" ))
  legend(5000, max(log10(theoretical.std)+0.1), c("theoretical", test.type), 
         lwd=c(5,2,2,2), cex=0.6, col=c("green", "red", "blue", "orange" ))
#  dev.off()
}



# Compute auto-correlation function of MCMC


# Add test for reproduciability of the pvalue: (Yaniv)
# Take B=100,500,1000,2000
# For each B run 20 iterations and compute the pvalues for the SAME dataset, P_1,...P_20.
# Check that P_1,..,P_20 look like they were P_i ~ Binom(B, p) / B   i.i.d. and don't have a higher variance 
# Here check also minP (for one time this is feasible)
test.pvalues.variance = 1; 
if(test.pvalues.variance)  
{
  # To complete ...
  
}









