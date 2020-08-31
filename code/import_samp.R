#############################################################
# Calculate p-value using importance sampling
# (If product of weights is very small/large, we can
#  multiply the weight function by a constant so the weights
#  are in some sense centered at 1)
#############################################################
IS.permute <- function(data, grid.points, B, w.fun=function(x){1}, expectations.table=c())
{
  n <- dim(data)[1]
  if(missing(grid.points) || isempty(grid.points))
    grid.points = data #  grid.points <- matrix(rnorm(2*n, 0, 1), n, 2) # TEMP DEBUG!!
  inverse.weight = missing(expectations.table) || isempty(expectations.table) # default is using inverse weighting 
  if(inverse.weight)
    TrueT <- ComputeStatistic.W(data, grid.points, w.fun)$Statistic # no unique in grid-points 
  else
    TrueT <- ComputeStatistic(data, grid.points, expectations.table)$Statistic # no unique in grid-points !!
  
  p.w <- matrix(0, B, 1)
  T.b <- matrix(0, B, 1) # statistics under null 
  for(b in 1:B)
  {
    perm <- sample(n)  # uniformly sampled permutation 
    if(inverse.weight)
      T.b[b] <- ComputeStatistic.W(cbind(data[,1], data[perm,2]), grid.points, w.fun)$Statistic # grid depends on permuted data
    else
      T.b[b] <- ComputeStatistic(cbind(data[,1], data[perm,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
    p.w[b] <- sum(log(w_fun_eval(data[,1], data[perm,2], w.fun))) # prod of w[i, pi(i)] 
  }
  p.w <- exp(p.w - max(p.w)) # shift max to prevent overflow 
  return(list(Pvalue=sum((T.b>=TrueT) * p.w) / sum(p.w), TrueT=TrueT, statistics.under.null=T.b))
}

#############################################################
# simulation for w(x,y)=x+y
# rep - number of replicaions
# B - number of permutations
# r - correlation (must be non-negative)
#############################################################
#simul <- function(rep=100,B=2000,r=0){
#  prep <- c()
#  for (i in 1:rep){
#    Z0 <- runif(200000)
#    Z1 <- sqrt(r)*Z0 + sqrt(1-r)*runif(200000)
#    Z2 <- sqrt(r)*Z0 + sqrt(1-r)*runif(200000)
#   samp <- sample(x = 200000,size = 100,replace = TRUE,prob = Z1+Z2)
#    x <- Z1[samp]
#    y <- Z2[samp]
#    dat <- data.frame(x,y)
#    res <- IS.permute(data = dat,B = 2000,w = function(x){sum(x)})
#    prep <- c(prep,res$p.val)
#    print(i)
#  }
#  return(prep)
#}
#
# prep <- simul(rep=100,B=2000,r=0.25)
# summary(prep)
# mean(prep<0.05)
# plot(ecdf(prep))
# abline(0,1)


