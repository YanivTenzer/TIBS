# This code compares two bootstrap methods for truncated data
# usnig the Tsai's test
# Both methods estimate the CDFs using the Aalen-Nelson estimate
# undrer the independence assumption, and then apply the 
# bootstrap by sampling from this CDF under truncation.
# 
# method 1 - sample from both CDFs keeping datum satysfing
#            the truncation constraint
# method 2 - sample only one coordinate (y) from the conditional
#            distribution given the 2nd (Y|Y>x)


# Tsai's test for uncensored data
Tsai.test <- function(L,Y){
  n <- length(Y)
  Risk <- vector(mode = "list", length = n)
  risk.size <- NULL
  for (i in 1:n) {
    set.i <- which(L<=Y[i] & Y>=Y[i])
    Risk[[i]] <- set.i
    risk.size <- c(risk.size,length(set.i))
  }
  S.i <- NULL
  for (i in 1:n){
    S.i[i] <- sum(sign(L[Risk[[i]]]-L[i]))
  }
  K.star <- sum(S.i)
  return(K.star)
}

# marginal estimationion 
# we use the Nelson-Aalen estimate to eliminate 0 probabilities
PL.estim <- function(L,Y,both=TRUE){
  require(survival)
  n <- length(Y)
  srv.S <- survfit(Surv(L,Y,rep(1,n))~1)
  sup.Y <- srv.S$time
  S.Y <- exp(-srv.S$cumhaz)
  p.Y <- c(1,S.Y[-n])-S.Y
  if (both==FALSE){
    return(list(sup.Y=sup.Y,p.Y=p.Y))
  } else {
    srv.S1 <- survfit(Surv(-Y,-L,rep(1,n))~1)
    sup.L <- -srv.S1$time
    S.L <- exp(-srv.S1$cumhaz)
    p.L <- c(1,S.L[-n])-S.L
    return(list(sup.Y=sup.Y,p.Y=p.Y,sup.L=sup.L,p.L=p.L))
  }
}

## bootstrap - sampling both coordinates

boot.both <- function(L,Y,B=200){
  estmarg <- PL.estim(L,Y)  # estimate marginals
  res <- c()
  n <- length(L)
  for (b in 1:B){
    # sample from truncated dist
    k <- 0
    samp.L <- c()
    samp.Y <- c()
    while (k<n){
      samp.LL <- sample(x = estmarg$sup.L, size = n-k, replace = TRUE, prob = estmarg$p.L)
      samp.YY <- sample(x = estmarg$sup.Y, size = n-k, replace = TRUE, prob = estmarg$p.Y)
      ok <- (samp.LL<samp.YY)
      samp.L <- c(samp.L,samp.LL[ok])
      samp.Y <- c(samp.Y,samp.YY[ok])
      k <- length(samp.L)
    }
    res.b <- Tsai.test(L = samp.L,Y = samp.Y) # test statistic
    res <- c(res,res.b)
  }
  res
}

## bootstrap - sampling one coordinate

boot.one <- function(L,Y,B=200){
  estmarg <- PL.estim(L,Y,both=TRUE)  # estimate marginal of Y
  res <- c()
  n <- length(L)
  for (b in 1:B){
    # sample from truncated dist
    samp.L <- L
    samp.Y <- c()
    for (i in 1:n){
      ind.i <- which(Y>L[i])
      if (min(estmarg$p.Y[ind.i])<=0) {return(list(samp.L=L,samp.Y=Y))}
      if (length(ind.i)==1) {
        samp.i <- ind.i
      } else{
        samp.i <- sample(x = ind.i, size = 1, replace = FALSE, prob = estmarg$p.Y[ind.i])
      }
      samp.Y[i] <- Y[samp.i]
    }
    res.b <- Tsai.test(L = samp.L,Y = samp.Y) # test statistic
    res <- c(res,res.b)
  }
  res
}

# Simulation to compare the two bootstrap methods

## sampling truncated data - bivariate normal
samp.trunc <- function(n,rho=0){
  library(mvtnorm)
  k <- 0
  L <- c()
  Y <- c()
  while(k<n){
    X <- rmvnorm(n-k,c(0,0),matrix(c(1,rho,rho,1),2,2))
    ok <- (X[,1]<X[,2])
    L <- c(L,X[ok,1])
    Y <- c(Y,X[ok,2])
    k <- length(Y)
  }
  list(L=L,Y=Y)
}

# running the simulation. 
# once in a while (rare), there are sampling problems,
# probably due to outliers in the sample that cuase
# the risk set to vanish (known problem in truncation).
# These are ignored using the tryCatch function
compare.boot <- function(n,rho,B=500,rep=100){
  p.val <- c()
  for (i in 1:rep){
    skipnext <- FALSE
    tryCatch({dat <- samp.trunc(n = n, rho=rho);
    boot1 <- boot.one(L = dat$L,Y=dat$Y,B = B);
    boot2 <- boot.both(L = dat$L,Y=dat$Y,B = B);
    tsai <- Tsai.test(L = dat$L, Y=dat$Y);
    p.v1 <- mean(abs(tsai)<=abs(boot1));
    p.v2 <- mean(abs(tsai)<=abs(boot2));
    p.val <- rbind(p.val,c(p.v1,p.v2));
    print(i)}, error = function(e) {skipnext <- TRUE})
    if (skipnext) {next}
  }
  p.val
}

p0 <- compare.boot(n = 200,rho = 0,B = 200,rep = 100)
plot(p0)
abline(0,1)
colMeans(p0<0.05)

p2 <- compare.boot(n = 200,rho = -0.2,B = 200,rep = 100)
plot(p2)
abline(0,1)
colMeans(p2<0.05)

p4 <- compare.boot(n = 200,rho = -0.4,B = 200,rep = 100)
plot(p4)
abline(0,1)
colMeans(p4<0.05)

p6 <- compare.boot(n = 200,rho = -0.6,B = 200,rep = 100)
plot(p6)
abline(0,1)
colMeans(p6<0.05)

p8 <- compare.boot(n = 200,rho = -0.8,B = 200,rep = 100)
plot(p8)
abline(0,1)
colMeans(p8<0.05)

#dat <- samp.trunc(1000,rho=-0.8)
#est <- PL.estim(dat$L,dat$Y)
