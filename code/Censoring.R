# TEST FOR LEFT-TRUNCATED RIGHT-CENSORED DATA
# THE NELSON-AALEN ESTIMATOR IS USED FOR THE CENSORING SURVIVAL FUNCTION
# THIS PREVENTS 0 AT THEE RIGHT TAIL

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(mvtnorm)
library(ggplot2)
library(survival)
source("TIBS.R")

# sampling left truncated right censored data
# log normal with correlation rho

samp.LTRC <- function(n=200,rho=0){
  k <- 0 
  X <- NULL; Y<- NULL
  while(k<n){
    XX <- rmvnorm(n = n-k,mean = c(0,0),sigma = matrix(c(1,rho,rho,1),2,2))
    keep <- I(XX[,1]<XX[,2])
    X <- c(X,exp(XX[keep,1]))
    Y <- c(Y,exp(XX[keep,2]))
    k <- length(X)
  }
  C <- rgamma(n,shape=4,scale=0.5)
  X.obs <- X
  Y.obs <- pmin(Y,X+C)
  delta <- (Y<= X+C)
  dat.LTRC <- data.frame(X.obs=X.obs,Y.obs=Y.obs,delta=delta)
  return(dat.LTRC)
}


### Simulation - comparing true weight fun to two estimates
n.sim <- 100
n.samp <- 200
prms <- list(B = 200)
rho <- -0.5
res <- c()
for (i in 1:n.sim){
  dat <- samp.LTRC(n=n.samp,rho=rho)
  KM <- survfit(Surv(Y.obs-X.obs,!delta) ~ 1,data=dat)
  Srv.C1 <- stepfun(KM$time,c(1,exp(-KM$cumhaz)))
  Srv.C2 <- stepfun(KM$time,c(1,KM$surv))
  Srv.C3 <- function(x){pgamma(x,shape=4,scale=0.5,lower.tail = FALSE)}
  only.uncens <- dat[dat$delta,1:2]
  w.fun1 <- function(x,y){(x<y)*Srv.C1(y-x)}
  w.fun2 <- function(x,y){(x<y)*Srv.C2(y-x)}
  w.fun3 <- function(x,y){(x<y)*Srv.C3(y-x)}
  tibs1 <- TIBS(data=only.uncens, bias.type=w.fun1, test.type='permutations',prms=prms)
  tibs2 <- TIBS(data=only.uncens, bias.type=w.fun2, test.type='permutations',prms=prms)
  tibs3 <- TIBS(data=only.uncens, bias.type=w.fun3, test.type='permutations',prms==prms)
  res <- rbind(res,c(tibs1$Pvalue,tibs2$Pvalue,tibs3$Pvalue,dim(only.uncens)[1]))
}

summary(res)
c(mean(res[,1]<0.05),mean(res[,2]<0.05),mean(res[,3]<0.05))


#############################################################
#############################################################
# a simple example 

dat <- samp.LTRC(n=200,rho=-0.3)
# descriptive
table(dat$delta)
ggplot(dat, aes(x = X.obs, y = Y.obs, colour = as.factor(delta))) +
  geom_point() + geom_abline(intercept = 0, slope = 1)

KM <- survfit(Surv(Y.obs-X.obs,!delta) ~ 1,data=dat)
# to eliminate 0 values, we use Nelson AAlen instead of Kaplan-Meier
Srv.C <- stepfun(KM$time,c(1,exp(-KM$cumhaz)))
curve(Srv.C(x),col="green",0,max(dat$Y.obs))
Srv.C <- function(x){pgamma(x,shape=4,scale=0.5,lower.tail = FALSE)}

only.uncens <- dat[dat$delta,1:2]
w.fun <- function(x,y){
  (x<y)*Srv.C(y-x)
}

w.dat <- w.fun(only.uncens$X.obs,only.uncens$Y.obs)
summary(w.dat)

TIBS(data=only.uncens, bias.type=w.fun,list(B = 1000), test.type='permutations',prms=c())

