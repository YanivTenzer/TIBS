# TEST FOR LEFT-TRUNCATED RIGHT-CENSORED DATA
# THE NELSON-AALEN ESTIMATOR IS USED FOR THE CENSORING SURVIVAL FUNCTION
# THIS PREVENTS 0 AT THEE RIGHT TAIL

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(mvtnorm)
library(ggplot2)
library(survival)
source("TIBS.R")
source("simulate_biased_sample.R")
Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 

source('simulate_and_test.R')
source('simulate_biased_sample.R')
source('TIBS.R')
source('marginal_estimation.R')
source('utilities.R')
source('import_samp.R')
library(foreach)
library(doSNOW)
#library(parallel)
library(doParallel)
library(gdata)
library(mvtnorm)
library(ggplot2)  
library(pracma)

cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 


# sampling left truncated right censored data
samp.LTRC <- function(n=200,rho=0, dependence.type){
  k <- 0 
  X <- NULL; Y<- NULL
  while(k<n){
    
    if(dependence.type=='micha_example')
    {
      XX <- rmvnorm(n = n-k,mean = c(0,0),sigma = matrix(c(1,rho,rho,1),2,2))
      keep <- I(XX[,1]<XX[,2])
      X <- c(X,exp(XX[keep,1]))
      Y <- c(Y,exp(XX[keep,2]))
      k <- length(X)
    }
    
    if(dependence.type=='Gaussian')
    {
      library(mvtnorm)
      xy<-rmvnorm(n=n-k, c(0,0), matrix(c(1, rho, rho,1),2,2))  
      keep <- I(xy[,1]<xy[,2])
      X_new = xy[keep==1,1]
      Y_new = xy[keep==1,2]
      X <- c(X,X_new)
      Y <- c(Y,Y_new)
      k <- length(X)
    }
    if(dependence.type=='LD')
    {
      library(copula)  
      GaussianCop<- normalCopula(param=rho, dim = 2, dispstr = "ex") # if needed
      ranks<- rCopula(n-k, GaussianCop)
      X_1<-qexp(ranks[,2], rate = 0.2)
      Y_1<-qweibull(ranks[,1], shape = 3, scale = 2, lower.tail = TRUE, log.p = FALSE)
      keep <- I(X_1<Y_1)
      X_new = X_1[keep==1]
      Y_new = Y_1[keep==1]
      X <- c(X,X_new)
      Y <- c(Y,Y_new)
      k <- length(X)
    }
    
    if(dependence.type=='nonmonotone_nonexchangeable')
    {
      library(copula)  # y ~ weibull, x ~ Gaussian copula 
      GaussianCop<- normalCopula(param=rho, dim = 2, dispstr = "ex") # if needed
      ranks<- rCopula(n-k, GaussianCop)
      Y_1<-qweibull(ranks[,1], shape = 0.5, scale = 2, lower.tail = TRUE, log.p = FALSE)
      X_1<-0.5 * (ranks[,2] * sample(c(-1,1), 1)+ 1)
      keep <- I(X_1<Y_1)
      X_new = X_1[keep==1]
      Y_new = Y_1[keep==1]
      X <- c(X,X_new)
      Y <- c(Y,Y_new)
      k <- length(X)
    }
    if(dependence.type=='CLmix')
    {
      library(copula)
      xy <- rCopula(n-k, claytonCopula(0.5-rbinom(1, 1, 0.5)))
      X_1<-xy[,1]
      Y_1<-xy[,2]
      keep <- I(X_1<Y_1)
      X_new = X_1[keep==1]
      Y_new = Y_1[keep==1]
      X <- c(X,X_new)
      Y <- c(Y,Y_new)
      k <- length(X)
    }
    if(dependence.type=='Gumbel')
    {
      library(copula)
      xy <- qnorm(rCopula(n-k, gumbelCopula(rho)))
      X_1<-xy[,1]
      Y_1<-xy[,2]
      keep <- I(X_1<Y_1)
      X_new = X_1[keep==1]
      Y_new = Y_1[keep==1]
      X <- c(X,X_new)
      Y <- c(Y,Y_new)
      k <- length(X)
    }
  }

  if(dependence.type %in% c("Gumbel", "CLmix"))
  {
    C <- rgamma(n,shape=1,scale=1.5)
  }else{
    if(dependence.type %in% c('LD'))
    {
      C <- rgamma(n,shape=1.5,scale=1.5)
    }else
    {
      C <- rgamma(n,shape=5,scale=1.5)
    }
  }
  X.obs <- X
  Y.obs <- pmin(Y,X+C)
  delta <- (Y<= X+C)
  dat.LTRC <- data.frame(X.obs=X.obs,Y.obs=Y.obs,delta=delta)
  return(dat.LTRC)
}


### Simulation - comparing true weight fun to two estimates
#n.sim <- 1000
#n.samp <-200
#B <- 100
#Prms<-list()
#Prms$naive.expectation = 0
#Prms$PL.expectation = 0
#Prms$B = B

Prms<-list(B=100, sample.size=200, iterations=1000, 
           naive.expectation = 0, PL.expectation=0,
           alpha=0.05, sequential.stopping=0, 
           use.cpp=0, keep.all=0, perturb.grid=1, 
           simulate.once=0, new.bootstrap=1) 

dependence.type = c('nonmonotone_nonexchangeable', 'LD', 'Gumbel', 'CLmix')
prms <- list(seq(0,0.2,0.1), seq(0,0.2, 0.1), c(1.6), c(0.5))
names(prms) <- c('nonmonotone_nonexchangeable', 'LD', 'Gumbel', 'CLmix')

output_dir <- '~/Documents/TIBS/code/LTRC' #(set your output dir)
########################################################################
# Main loop
#######################################################################
for(dependence in dependence.type)
{
  idx = 0
  res <- c()
  for(prm in prms[[dependence]])
  {
    idx = idx+1
    for(rep_idx in seq(Prms$iterations))
    {
      dat <- samp.LTRC(n=Prms$sample.size,rho=prm, dependence)
      only.uncens <- dat[dat$delta==1,1:2]
      KM <- survfit(Surv(Y.obs-X.obs,!delta) ~ 1,data=dat)
      Srv.C1 <- stepfun(KM$time,c(1,exp(-KM$cumhaz)))
      w.fun1 <- function(x,y){(x<y)*Srv.C1(y-x)}
      
      
      #TIBS, permutations:
      tibs1 <- TIBS(only.uncens, w.fun=w.fun1, test.type='permutations',Prms)
      #TIBS, importance_sampling:
      tibs2<-TIBS(data=only.uncens, w.fun=w.fun1, 
                  test.type='uniform_importance_sampling',Prms)
      
      tibs3<-TIBS(data=only.uncens, w.fun=w.fun1, 
                  test.type='permutations_inverse_weighting',prms=Prms)
      
      
      #minP2:
      Prms_minP2<-Prms
      Prms_minP2$delta <- dat$delta
      minP2 <- TIBS(data=dat[,1:2], w.fun1, test.type='minP2',Prms_minP2)
      
      #Tsai:
      tsai <- TIBS(data=dat[,1:2], NA, test.type='tsai',NA)
      
      res <- rbind(res,c(tibs1$Pvalue, tibs2$Pvalue, tibs3$Pvalue, minP2$Pvalue, tsai$Pvalue, dim(only.uncens)[1]))
    }
    res = as.data.frame(res)
    colnames(res)<-c('tibs1', 'tibs2', 'tibs3', 'minP2', 'tsai', 'n')
    save(dependence.type, prm, res,file=paste0(output_dir, '/LTRC_', dependence, toString(idx),'.Rdata'))
  }
}
rm(list=ls())
#collect results:
Power <- c()
for(dependence in dependence.type)
{
  idx = 0
  for(prm in prms[[dependence]])
  {
    idx = idx+1
    load(paste0(output_dir, '/LTRC_', dependence, toString(idx),'.Rdata'))
    Power<-rbind(Power,colMeans(apply(res[,seq(5)],2,function(x) ifelse(x<0.05, 1, 0))))
  }
  Power<-data.frame(Power)
}
################################################################
#  end of main loop
#############################################################
power<-c()
for(idx in seq(1,10))
{
  #File = paste0('~/Dropbox/cond_ind/code/censoring/Gaussian_with_censoring/','Gaussian_', toString(idx), '.Rdata')
  File = paste0('~/Dropbox/cond_ind/code/LD', toString(idx), '.Rdata')
  load(File)
  power<-rbind(power, c(mean(res[,1]<0.05),mean(res[,2]<0.05),mean(res[,3]<0.05), mean(res[,4]<0.05)))
}

power = data.frame(power)
colnames(power)<-c('tibs1', 'tibs2', 'tibs3', 'minp2')
library(cast)
power_melt<-melt(power)
power_melt<-data.frame(power_melt)
colnames(power_melt)<-c("test", 'power')
rho = rep(seq(0, 0.9, 0.1), 4)
power_melt<-cbind(power_melt, rho)

jpeg('LD_35_with_censoring.jpeg')
p<-ggplot(power_melt, aes(x=rho, y=power, group=test, linetype=test)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.ticks.y = element_blank())+#,
  #legend.position = "none")+
  scale_y_continuous(limits = c(0, 1), breaks=c(0, 0.3, 0.6, 0.9))+
  geom_line(aes(color=test),size=1.2)+
  geom_point(aes(color=test, shape=test), stroke = 2)+
  scale_shape_manual(values = c(5, 3,4,19,25,7,8,9))+
  scale_size_manual(values=rep(20,7)) 
print(p)
dev.off()
#############################################################