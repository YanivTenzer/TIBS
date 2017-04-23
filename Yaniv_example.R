cat("\014")
rm(list=ls())
path = 'C:/Users/oorzu/Documents/GitHub/test_for_independence_under_biased_sampling' # (Or's path)  # 'D:/cond_ind_local'; # (Yaniv's path)
setwd(path)
source('MCMC_Permutations.R')  
library(mvtnorm)

start.time <- Sys.time()
########################
#This is a toy example, simulating bivariate vector (x,y) \in {0,1,2}x{0,1,2}):
#1. Define P(X,Y):
probability<-matrix(0,3,3)
probability[1,1]=1/8
probability[2,2]=1/8
probability[3,3]=1/4
probability[1,2]=1/20
probability[1,3]=1/20
probability[2,1]=1/20
probability[2,3]=2/20
probability[3,1]=4/20
probability[3,2]=1/20
#2. W*P(x,y):
biased_probability<-matrix(0,3,3) #w(0,0)=w(0,1)=w(0,2)=w(1,2)=w(2,1)=0:
S<-sum(colSums(probability))-(sum(probability[1,])+probability[2,3] +0*probability[3,2])
#Normalize the resulting distribution:
biased_probability[2:3,1:3]<-probability[2:3,1:3]/S
biased_probability[2,3]=0
biased_probability[3,2]=probability[3,2] # 0 # Add another point to destroy quasi-independence
#3. Sample from  W*P(x,y):
n=3000 # num. permutations

idx <- sample(1:3^2, n, prob=as.vector(biased_probability), replace = TRUE) # weighted sampling from un-truncated values 

Sample<-matrix(0,n,2) # Sample truncated data 
for(i in 1:3^2)
{
  Sample[idx==i,1] <- (i-1)%%3
  Sample[idx==i,2] <- floor((i-1)/3)
}  

#############################
#4. Define the weight matrix:
Weights_Matrix_1<-matrix(0, dim(Sample)[1],dim(Sample)[1])
Weights_Matrix_independent<-matrix(0, dim(Sample)[1],dim(Sample)[1])
for(i in seq(1,dim(Sample)[1],1))
{
  temp<-cbind(rep(Sample[i,1], dim(Sample)[1]), Sample[,2])  
  Weights_Matrix_1[i,]<-biased_probability[Sample[i,1]+1, temp[,2]+1] 
  Weights_Matrix_independent[i,]<-ifelse(biased_probability[Sample[i,1]+1, temp[,2]+1]>0, 1,0) # take W 
}

#5. Sample valid permutations:
TargetSampleSize = 3000
Permutations<-MCMC_Permutations(NULL,Weights_Matrix_1,TargetSampleSize,dim(Weights_Matrix_1)[1])
Permutations_independent<-MCMC_Permutations(NULL,Weights_Matrix_independent,TargetSampleSize,dim(Weights_Matrix_independent)[1])
#6. We are set:
temp<-matrix(0,TargetSampleSize,1)
temp_independent<-matrix(0,TargetSampleSize,1)
####################################
p_test<-c(1,0)#c(1,1)#c(2,2)#c(2,0)

for(i in seq(1,TargetSampleSize,1))
{
  permuted_data<-cbind(Sample[,1], Sample[Permutations[,i],2])
  temp[i]<-length(which(permuted_data[,1]==p_test[1] & permuted_data[,2]==p_test[2] ))/dim(permuted_data)[1]
}
cat('Sampling with the induced distribution: ', toString(mean(temp)))

for(i in seq(1,TargetSampleSize,1))
{
  permuted_data<-cbind(Sample[,1], Sample[Permutations_independent[,i],2])
  temp_independent[i]<-length(which(permuted_data[,1]==p_test[1] & permuted_data[,2]==p_test[2] ))/dim(permuted_data)[1]
}
cat('Sampling with the P_W: ', toString(mean(temp_independent)))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#########################################################
temp<-matrix(0, dim(Permutations)[2],1)
for(i in seq(1, dim(Permutations)[2],1))
{
  permuted_data<-cbind(data[,1], data[Permutations[,i],2])
  temp[i]<-length(which(permuted_data[,1]<0.5 & permuted_data[,2]<0.5))/dim(permuted_data)[1]
}
mean(temp)
######Now compare with this:
rho=0.6
sigma<-matrix(c(1, rho, rho , 1),2,2)
mu=c(0,0)
pmvnorm(upper=c(0.5, 0.5),mean=c(0,0), corr=sigma)