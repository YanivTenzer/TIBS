# Set path to your working directory: 
# The script contains real data analysis 
# Due to privacy issues, the data files for the ICU datasets are not part of the released package. 
# Please contact the authors if you are interested in them. 

rm(list=ls())
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path = getwd()
Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 

source('TIBS.R')
source('simulate_biased_sample.R')
source('marginal_estimation.R')
source('Tsai_test.R')
source('utilities.R') 
source("import_samp.R")
library(lubridate)
library(stringr)
library(mvtnorm)
library(matrixStats)
library(permDep)
library(tictoc)
library(xtable)
library(ggplot2)
library(ggmap)
library(pracma)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
load("Channing_run_IS.Rdata")


datasets <- c('huji', 'AIDS', 'ICU', 'Infection',  'ChanningHouse')  # Run Channing and skipping Dementia. 
w.fun <- c('huji', 'Hyperplane_Truncation', 'Hyperplane_Truncation', 'sum', 'Hyperplane_Truncation_censoring') # last one is dementia (sum?)

# Read real dataset 
ReadDataset <- function(data_str)
{
  naive.w.fun = "const"
  switch(data_str,  # read dataset 
         'huji'={
           load('../data/APage_APtime.RData')
           input.data <- HUJI.dat
           w.fun <- 'huji'
         },
         'AIDS'={
           library(DTDA)
           data("AIDS")
           input.data <- AIDS[,2:3]
           w.fun <- 'Hyperplane_Truncation'
         },
         'ICU'={  # this dataset is not part of the released package
           input.data <- read.table("../data/ICU_data.txt", header = TRUE)
           input.data <- input.data[,c(1:2)]  # discard third column 
           input.data[,2] <- input.data[,2]-0.02 # correction: reduce by 1: the truncation criterion is L<TU (L==TU is removed)
           w.fun <- 'Hyperplane_Truncation'
         }, 
         'Infection'={ # this dataset is not part of the released package
           load('../data/ICU_INF.Rdata')
           input.data <- cbind(X,Y)
#           w.max <- max(X)+max(Y) # for truncation we need 
           w.fun <- 'sum'
         }, 
         'Dementia'={ # this dataset is not part of the released package
           # DEMENTIA DATA
           first_time <- 0
           if(file.exists('../data/Dementia.Rdata')) ### Get data.
           {
             load('../data/Dementia.Rdata')  # should be added 
           } else
           {
             prevdata<-read.csv("C:/Users/mm/Dropbox/marco/Nonparametric bivariate estimation/Old folder/Data analysis/cshaforRCSV.csv"); # change to relative path 
             ### Clean data.
             badindex<-which(prevdata$duration-prevdata$truncation<=0);
             prevdata <- prevdata[-badindex,];
             
             ### Order data according to disease duration.
             x<-prevdata$AAO;
             v<-prevdata$duration;
             w<-prevdata$AAO+prevdata$truncation;
             delta<-prevdata$death;
             order.v<-order(v);
             x<-x[order.v];
             v<-v[order.v];
             w<-w[order.v];
             delta<-delta[order.v];
             
             ### Create data frame.
             cshadata <- data.frame(list("x"=x,"v"=v,"w"=w,"delta"=delta))  # cshadata
             
             # the risk set vanishes before the last observation. 
             # remove the largest 3 v's to take care of that.
             id.rs <- which(cshadata$v>25)
             cshadata1 <- cshadata[-id.rs,] 
             
             # before ommiting the largest 3
             KM <- survfit(Surv(v-(w-x),!delta) ~ 1,data=cshadata)
             Srv.C2 <- stepfun(KM$time,c(1,KM$surv))
             
             w.fun2 <- function(x,y){(x<y)*Srv.C2(y-x)}
             x.csha <- cshadata$w-cshadata$x
             y.csha <- cshadata$v
             delta.csha <- cshadata$delta
             input.data <- data.frame(x.csha[delta.csha],y.csha[delta.csha])  # csha.delta1
             
             # Save data frame (next time we can load only this)
             save(input.data, KM, file='../data/Dementia.Rdata')          
           }                
           Srv.C1 <- stepfun(KM$time,c(1,exp(-KM$cumhaz)))
           w.fun <- function(x,y){(x<y)*Srv.C1(y-x)}  # modify w.fun 
         }, # end Dementia
         'ChanningHouse'={ # Read Channing House dataset : FILL CODE
           library(survival)
           data("channing", package = "boot")
           ok <- which(channing$exit>channing$entry) # keep only these 
           ### Create data frame.
           input.data <- data.frame(list("x"=channing$entry[ok], "y"=channing$exit[ok], "delta"=channing$cens[ok]))
           # filter the uncensored observations and apply TIBS
#           input.data <- data.frame(x=ch.data$x[ch.data$delta==1], y=ch.data$y[ch.data$delta==1])  # keep all values (including censored!)
           
           # estimating censoring distribution and calculating W function
           KM <- survfit(Surv(y-x,!delta) ~ 1, data=input.data)
           Srv.C <- stepfun(KM$time, c(1,KM$surv))
           w.fun <- function(x,y){(x<y)*Srv.C(y-x)}
           naive.w.fun <- function(x,y){Srv.C(y-x)}
#           return(list(input.data=input.data, w.fun=w.fun, delta=ch.data$delta)) # return also delta for censoring
           
         }
  ) # end switch 
  
  if(!is.numeric(input.data))   # unlist and keep dimensions for data 
    input.data <- array(as.numeric(unlist(input.data)), dim(input.data))  
  
  return(list(input.data=input.data, w.fun=w.fun, naive.w.fun=naive.w.fun)) # new! return list with also w (also update for other datasets).  
}



# example of a cpp function
cppFunction('NumericVector timesTwo(NumericVector x) {
  return x * 2;
}')

########################################################
# "Official parameters"
test.method <- c('bootstrap', 'permutationsMCMC', "permutationsIS", 'tsai', 'minP2') # different tests to run 
test.stat <- c("adjusted_w_hoeffding", "adjusted_w_hoeffding", "adjusted_w_hoeffding", "tsai", "minP2")
IS.method <- c("tsai", "KouMcculough.w")
exchange.type <- c(FALSE, FALSE, FALSE, FALSE, FALSE) # no reason to assume real data is exchangeable
B = 100000 # number of permutations/bootstrap samples 
minP2.B = 10000 # for minP take a smaller value due to running time 
plot.flag <- 0
########################################################

# parameters for experimenting
test.method <- c('permutationsMCMC') # different tests to run 
test.stat <- c("adjusted_w_hoeffding")
B = 99
minP2.B = 5
plot.flag <- 1
########################################################

# w.fun <- c('huji', 'Hyperplane_Truncation', 'Hyperplane_Truncation', 'sum', 'sum') # last one is dementia (sum?)
# w.max <- c(65, 1, 1, -1) # -1 denotes calculate max of w from data  
n.datasets <- length(datasets)
n.tests <- length(test.method)
Hyperplane.prms<-c(-1,1,0)
test.pvalue <- matrix(-1, n.datasets, n.tests) # -1 denotes test wasn't performed
test.time <- matrix(-0.01, n.datasets, n.tests) # -1 denotes test wasn't performed

overll.time <-  Sys.time()
for(d in 5:5) # temp just check channing 1:n.datasets) # loop on datasets (last is dementia)
{
#    if(datasets[d] %in% c('ICU', 'Infection'))  # these datasets are available by request. Uncomment this line if you don't have them 
#      next
  dat <- ReadDataset(datasets[d]) # read dataset 

  prms <- list(B = B)  # NOTE: for minP we may need a lower number of permutations
  prms$w.max <- set_w_max_sample(dat$input.data, dat$w.fun) # set maximum value of w for the dataset
    
  if(dim(dat$input.data)[2]==3)  # enable censoring (delta)
  {
      prms$delta <- dat$input.data[,3]
      dat$input.data <- dat$input.data[,1:2]
  }
  
#    max(w.max[d], max(w_fun_to_mat(dat$input.data, dat$w.fun, prms))) # update max 
  prms$use.cpp <- 1 # New! enable one to run with c++ code (faster)
  for(t in 1:n.tests) # run all tests 
  {
    if(test.method[t] == 'minP2')  # smaller sample size for minP2 (slower)
      prms$B = minP2.B
    else
      prms$B = B
    if(is.character(dat$w.fun))
      if((!(dat$w.fun %in% c('truncation', 'Hyperplane_Truncation'))) & (test.method[t] %in% c("tsai", 'minP2')))
        next  # these tests run only for truncation 
    if(is.character(dat$w.fun))
      if((test.method[t] == 'bootstrap') & (dat$w.fun %in% c('truncation', 'Hyperplane_Truncation', 'huji')))
        next # can't run bootstrap because w can be zero 
    if((test.method[t] == 'bootstrap') & (min(w_fun_eval(dat$input.data[,1], dat$input.data[,2], dat$w.fun, prms))==0))  # check for icu that we can run it
      next # can't run bootstrap because w can be zero 


    set.seed(100)
    
    print(paste0("Running ", datasets[d], ", ", test.method[t], ":"))
    start.time <- Sys.time()
    
    if(test.method[t] == "permutationsIS")
    {
      if(!is_pos_w(dat$w.fun, dat$input.data, 1))
        prms$importance.sampling.dist = "tsai" # for truncation x2>=x1 or truncation + censoring
      else
        prms$importance.sampling.dist = "KouMcculough.w"
    }
    
    save(dat, prms, test.method, test.stat, t, file="Channing_run_IS.Rdata")
    results.test <- TIBS(dat$input.data, dat$w.fun, prms, test.method[t], test.stat[t])  # can also be cpp 
    naive.results.test <- TIBS(dat$input.data, dat$naive.w.fun, prms, test.method[t], test.stat[t])  # see the effect of w
    difftime(Sys.time() , start.time, units="secs")
    test.time[d,t] <- as.numeric(difftime(Sys.time() , start.time, units="secs")) # Sys.time() - start.time
    print(paste0("test time: ", test.time[d,t]))
    test.pvalue[d,t] <- results.test$Pvalue 
    cat(datasets[d], ', ', test.method[t], ', Pvalue:', test.pvalue[d,t], '\n')
    if(plot.flag & (test.method[t] == "permutationsMCMC")) # permutations
    {
      if(("delta" %in% names(prms)) && (length(prms$delta) == dim(dat$input.data)[1]))
        xy <- cbind(data.frame(dat$input.data[which(prms$delta==1),]), as.data.frame(results.test$permuted.data))
      else
        xy <- cbind(data.frame(dat$input.data), as.data.frame(results.test$permuted.data))
      
      x.vec <- seq(min(xy[,1]), max(xy[,1]), (max(xy[,1])-min(xy[,1]))/100)
      y.vec <- seq(min(xy[,2]), max(xy[,2]), (max(xy[,2])-min(xy[,2]))/100)
      n <- length(x.vec)
      z.mat <- matrix(0, n,n)
      for(i in c(1:n))
        for(j in c(1:n))
          z.mat[i,j] <- w_fun_eval(x.vec[i], y.vec[j], dat$w.fun)
      
      image(x.vec, y.vec, z.mat, col=gray.colors(33), xlab="X", ylab="Y")
      levelplot(z.mat)
#      contour(x.vec, y.vec, z.mat)
      image.plot(x.vec, y.vec, z.mat, col = gray(100:0/100), xlab="X", ylab="Y")
      points(xy[,1], xy[,2], col="red", pch=16)
      points(xy[,3], xy[,4], col="green", pch=3, cex=0.5)
      
      ggplot(xy, aes(x=xy[,1], y=xy[,2])) + 
        geom_point(aes(x=xy[,1], y=xy[,2]), col="red" , size=2) + 
        geom_point(shape=3, aes(x=xy[,3], y=xy[,4]), col="black" , size=2) + 
#        ggtitle(datasets[d]) +
        xlab("X") + ylab("Y") +
        theme( # plot.title = element_text(size=14, face="bold.italic", hjust=0.5),
              axis.title.y = element_text(face="bold", size=14),
              axis.title.x = element_text(face="bold", size = 14),
              axis.text.x = element_text(face="bold", size=12), 
              axis.text.y = element_text(face="bold", size=12), 
              legend.position = "none")
##      ggsave(paste(path, '/../figs/real_data/', datasets[d], '.png', sep=''),
##          width=5, height=5, dpi=300)
    }
  } # end loop on tests 
} # end loop on datasets   

# save table with p-values and running times  
results.table <- cbind(test.pvalue, test.time)
rownames(test.pvalue) <- datasets
rownames(test.time) <- datasets
colnames(test.pvalue) <- test.method
colnames(test.time) <- test.method

rownames(results.table) <- datasets
# colnames(results.table) <- test.method
save(test.pvalue, test.time, test.method, 
     file=paste0(path ,'/../docs/Tables/real_datasets_B_', B, '.Rdata'))
# print(xtable(results.table, type = "latex"), 
#      file = paste0(path ,'/../docs/Tables/real_datasets.tex')) # save also in latex format 


print(paste0("Overall real-life-datasets time (sec.)  : ", sum(test.time)))


#####################################################

