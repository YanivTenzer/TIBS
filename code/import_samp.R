##################################################################################################
# Compute the modified Hoeffding's test statistic corrected for w:
# Observed and Expected calculated using inverse weighting
# GOOD ONLY FOR A STRICTLY POSITIVE w
# Parameters: 
# data - n*2 matrix with (x,y) sample
# grid.points - all possible (x_i,y_j) points  
# w - weight function (default w(x,y)=1) 
# 
#  Quardants convension:
#   4 | 1
#   ------
#   3 | 2
##################################################################################################

ComputeStatistic.W<- function(data, grid.points,w=function(x){1}){
  W <- apply(data,1,w)
  n.w <- sum(1/W)
  obs.table<-exp.table <- matrix(0, dim(grid.points)[1], 4)
  Obs<-Exp<-matrix(0,4,1) # observed & expected
  Statistic <- 0 
  for (i in 1:dim(grid.points)[1])  # Slow loop on grid points 
  {
    Rx <- data[,1]>grid.points[i,1]
    Ry <- data[,2]>grid.points[i,2]
    Exp[1] <- sum(Rx/W)*sum(Ry/W)/n.w^2
    Exp[2] <- sum(Rx/W)*sum((!Ry)/W)/n.w^2
    Exp[4] <- sum((!Rx)/W)*sum(Ry/W)/n.w^2
    Exp[3] <- sum((!Rx)/W)*sum((!Ry)/W)/n.w^2
    Obs[1] <- sum(Rx*Ry/W)/n.w
    Obs[2] <- sum(Rx*(!Ry)/W)/n.w
    Obs[4] <- sum((!Rx)*Ry/W)/n.w
    Obs[3] <- sum((!Rx)*(!Ry)/W)/n.w
    obs.table[i,] <- Obs
    exp.table[i,] <- Exp
    if (min(Exp)>(1/dim(data)[1])) {
      Statistic <-  Statistic + sum((Obs-Exp)^2 / Exp) # set valid statistic when expected is 0 or very small 
    } 
  } # end loop on grid points 
  
  return(list(Statistic=Statistic, obs.table=obs.table,exp.table=exp.table)) 
  # returns also expected and observed tables for diagnostics
}


#############################################################
# Calculate p-value using importance sampling
# (If product of weights is very small/large, we can
#  multiply the weight function by a constant so the weights
#  are in some sense centered at 1)
#############################################################

IS.permute <- function(data,B,w=function(x){1}){
  n <- dim(data)[1]
  T.obs <- ComputeStatistic.W(data,data,w=w)$Statistic 
  reject <- 0
  sum.p <- 0
  for (b in 1:B){
    perm <- sample(n)
    dat.b <- data.frame(x=data[1:n,1],y=data[perm,2]) # permuted data
    T.b <- ComputeStatistic.W(dat.b,dat.b,w=w)$Statistic # grid depends on permuted data
    W <- apply(data,1,w)
    p.w <- prod(W)   
    reject <- reject+(T.b>=T.obs)/p.w
    sum.p <- sum.p+1/p.w
  }
  return(list(p.val=reject/sum.p, T.obs=T.obs))
}

# simulation for w(x,y)=x+y
# rep - number of replicaions
# B - number of permutations
# r - correlation (must be non-negative)

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

prep <- simul(rep=100,B=2000,r=0.25)
summary(prep)
mean(prep<0.05)
plot(ecdf(prep))
abline(0,1)


