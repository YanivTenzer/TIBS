########################################################################
# Perform Test for Independence under general Biased Sampling (TIBS)
# 
# Parameters: 
# data - n*2 array of (x,y) samples
# bias.type - string indicating biased sampling function W 
# B - number of bootstrap/permutation samples to perform 
# test.type - test to perform 
# prms - additional parameters (needed for bootstrap)
########################################################################
TIBS <- function(data, bias.type, B, test.type, prms)
{  
  library(pracma)
  source('utilities.R')
  source('Tsai_test.R')
  #################################################################
  #1. Compute weights matrix W:  
  W=GetBiasedSamplingWeights(data, dim(data)[1], bias.type)
  #2.Create a grid of points, based on the data:
  permutation<-PermutationsMCMC(W, 1, dim(data)[1]) # why always sample a permutation using MCMC? 
  temp.data<-cbind(data[,1], data[permutation,2])
  #discard the extremum points to avoid numerical issues
  idx.minmax <-which( (temp.data[,1] %in% c(min(temp.data[,1],max(temp.data[,1])))) | 
                        (temp.data[,2] %in% c(min(temp.data[,2],max(temp.data[,2])))) )
  grid.points<-temp.data[-idx.minmax,]
  
  switch(test.type,
         'bootstrap'={
           fast.bootstrap <- 0 # 0 # New: allow not to reestimate nulls for each bootstrap sample 
           #3. Estimate the marginals
           marginals<-EstimateMarginals(data, bias.type)
           W=GetBiasedSamplingWeights(marginals$xy, dim(marginals$xy)[1], bias.type) # compute W again for augmented data
           
           # Temp: set real marginals (Gaussians) for debug           
           analytic_marginals <- 0
           if(analytic_marginals)
           {
             marginals$PDFs[,1] <- dnorm(marginals$xy[,1]) / sum(dnorm(marginals$xy[,1]))           
             marginals$PDFs[,2] <- dnorm(marginals$xy[,2]) / sum(dnorm(marginals$xy[,2]))   
             Px <- sort(marginals$xy[,1], index.return=TRUE)  # Permute to order x_i, y_i
             marginals$CDFs[Px$ix,1] <- cumsum(marginals$PDFs[Px$ix,1])
             Py <- sort(marginals$xy[,2], index.return=TRUE)  # Permute to order x_i, y_i
             marginals$CDFs[Py$ix,2] <- cumsum(marginals$PDFs[Py$ix,2])
           }
           #4. Estimate W(x,y)*Fx*Fy/normalizing.factor
           null.distribution<-GetNullDistribution(marginals$PDF, W)
           normalizing.factor<-null.distribution$normalizing.factor
           null.distribution<-null.distribution$null.distribution
           
           #1. First compute the statistics based on the original data set:
           expectations.table<-QuarterProbFromBootstrap(marginals$xy, null.distribution, grid.points) # data
           TrueT=ComputeStatistic(data, grid.points, expectations.table)
           obs.table <- TrueT$obs.table
           TrueT <- TrueT$Statistic
           #2. Compute statistic for dummy sample:
           statistics.under.null=matrix(0, B, 1)
           null.distribution.bootstrap<-null.distribution
           for(ctr in 1:B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Run Boots=", ctr))
             dummy.sample<-Bootstrap(marginals$xy, marginals$PDF, bias.type, prms, dim(data)[1]) # draw new sample. Problem: which pdf and data? 
             if(!prms$fast.bootstrap)
             {
               marginals.bootstrap<-EstimateMarginals(dummy.sample, bias.type)  # Why are the marginals estimated each time? 
               pdfs.bootstrap<-marginals.bootstrap$PDFs
               #3. Compute weights matrix W:    
               W.bootstrap=GetBiasedSamplingWeights(marginals.bootstrap$xy, dim(marginals.bootstrap$xy)[1], bias.type)
               #4. Estimate W(x,y)*Fx*FY/normalizing.factor
               null.distribution.bootstrap<-GetNullDistribution(pdfs.bootstrap, W.bootstrap)
               normalizing.factor.bootstrap<-null.distribution.bootstrap$normalizing.factor
               null.distribution.bootstrap<-null.distribution.bootstrap$null.distribution
               expectations.table<-QuarterProbFromBootstrap(data, null.distribution.bootstrap, grid.points)
             } # MOST TIME SPENT HERE (STAT COMPUTE!) TEMP for profiling
             NullT <- ComputeStatistic(dummy.sample, grid.points, expectations.table)
             null.obs.table <- NullT$obs.table
             statistics.under.null[ctr] <- NullT$Statistic
           }
           output<-list(TrueT=TrueT,statistics.under.null=statistics.under.null)
         },
         'permutations'={
           Permutations=PermutationsMCMC(W, B, dim(data)[1])
           expectations.table<-QuarterProbFromPermutations(data, Permutations, grid.points)
           TrueT=ComputeStatistic(data, grid.points, expectations.table)$Statistic
           
           #Compute the statistics value for each permutation:
           statistics.under.null=matrix(0, B, 1)
           for(ctr in 1:B) 
           {
             if(mod(ctr,100)==0)
               print(paste0("Comp. Stat. Perm=", ctr))
             statistics.under.null[ctr] <- ComputeStatistic(
               cbind(data[,1], data[Permutations[,ctr],2]), grid.points, expectations.table)$Statistic
           }
           output<-list(TrueT=TrueT,statistics.under.null=statistics.under.null)
         },
         'tsai' = { result<-TsaiTestTies(data[,1],data[,2]) #Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
         output<-list(Pvalue=result[2])
         },
         'minP2' = {library(permDep)  #MinP2 test, relevant only for truncation W(x,y)=1_{x<=y}
           dat<-data.frame(list(trun = data[,1], obs = data[,2], delta = rep(1, dim(data)[1])))
           results<-with(dat, permDep(trun, obs, B, delta, nc = 1, minp2Only = TRUE, 
                                      sampling = 'conditional', kendallOnly = FALSE))
           output<-list(Pvalue=results$p.valueMinp2)
         }
  )
  output$permuted.data <- temp.data # add example of permuted data 
  if(!("Pvalue" %in% names(output))) # Compute empirical P-value
    output$Pvalue <- length(which(output$statistics.under.null>=output$TrueT))/B
  return(output)
}
