# Run diagnostic for weighted permutations tests 
path = 'C:/Users/Or Zuk/Documents/GitHub/TIBS/code'  # change to your path
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


# 1. Check overdispersion of variance of p-vals calculatoin
B_vec <- seq(200, 2000, 200) # try different numbers of permutations
iters <- 25 # number of times to estimate for the same dataset: 
test.type <- 'bootstrap'
prms <- c()
rho.max <- -0.05 # small correlations 

n_vec <- c(100, 200, 500)
prms$n <- n_vec[1]  # sample size 
prms$use.cpp <- 1 
w.fun <- 'truncation'
dependence.type <- 'LogNormal'
pvals <- matrix(0, iters, length(B_vec))

biased.data <- SimulateBiasedSample(prms$n, dependence.type, w.fun, prms) 
if(!('w.max' %in% names(prms)))
  prms$w.max <- biased.data$w.max


# p.zero <- TIBS(biased.data$data, w.fun, test.type, prms)
# p.zero.cpp <- TIBS_rcpp(biased.data$data, w.fun, test.type, prms)

for(n in n_vec)
{
  prms$rho <- rho.max * min(n_vec) / n  # change rho to keep power the same 
  prms$n <- n
  biased.data <- SimulateBiasedSample(prms$n, dependence.type, w.fun, prms) 
  for(i in c(1:iters))
  {
    print(paste0("run i=", i, " out of ", iters))
    for(b in 1:length(B_vec))
    {
      prms$B <- B_vec[b]
      pvals[i,b] <- TIBS(biased.data$data, w.fun, test.type, prms)$Pvalue #  TIBS(biased.data, w.fun, cur.test.type, prms)
    }
  }
  true.pval <- mean(pvals)
  
  theoretical.std <- sqrt(true.pval * (1-true.pval) / B_vec)
  empirical.std <- sqrt(colVars(pvals))
  
  jpeg(paste0("../figs/var_p_hat_", test.type, "_", w.fun, "_n_", n, ".jpg"), width = 400, height = 400)
  plot(B_vec, log10(theoretical.std), type="l", col="green", lwd=5, xlab="B", ylab="log_{10}(std)", 
       ylim=log10(range(empirical.std, theoretical.std)),
       main=paste0(test.type, " std(p.hat) for n=", prms$n, " p.val=", round(true.pval, 3) ))
  lines(B_vec, log10(empirical.std), col="red")
  legend(max(B_vec)*0.7, max(log10(empirical.std)), c("theoretical", "empirical"), lwd=c(5,2), col=c("green", "red"))
  dev.off()
}



# Compute auto-correlation function of MCMC
