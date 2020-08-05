#############################################################
# Calculate p-value using importance sampling
# (If product of weights is very small/large, we can
#  multiply the weight function by a constant so the weights
#  are in some sense centered at 1)
#############################################################

IS.permute <- function(data,B,w=function(x){1}){
  n <- dim(data)[1]
  TrueT <- ComputeStatistic.W(data,data,w=w)$Statistic # no unique? 
  reject <- 0
  sum.p <- 0
  for (b in 1:B){
    perm <- sample(n)
    T.b <- ComputeStatistic.W(cbind(data[,1], data[perm,2]), cbind(data[,1], data[perm,2]), w=w)$Statistic # grid depends on permuted data
    W <- w_fun_eval(data[,1], data[,2], w) #     W <- apply(data,1,w)
    p.w <- prod(W)   # could cause overflow 
    reject <- reject+(T.b>=TrueT)/p.w
    sum.p <- sum.p+1/p.w
  }
  return(list(Pvalue=reject/sum.p, TrueT=TrueT))
}

#############################################################
# simulation for w(x,y)=x+y
# rep - number of replicaions
# B - number of permutations
# r - correlation (must be non-negative)
#############################################################
simul <- function(rep=100,B=2000,r=0){
  prep <- c()
  for (i in 1:rep){
    Z0 <- runif(200000)
    Z1 <- sqrt(r)*Z0 + sqrt(1-r)*runif(200000)
    Z2 <- sqrt(r)*Z0 + sqrt(1-r)*runif(200000)
    samp <- sample(x = 200000,size = 100,replace = TRUE,prob = Z1+Z2)
    x <- Z1[samp]
    y <- Z2[samp]
    dat <- data.frame(x,y)
    res <- IS.permute(data = dat,B = 2000,w = function(x){sum(x)})
    prep <- c(prep,res$p.val)
    print(i)
  }
  return(prep)
}

# prep <- simul(rep=100,B=2000,r=0.25)
# summary(prep)
# mean(prep<0.05)
# plot(ecdf(prep))
# abline(0,1)


