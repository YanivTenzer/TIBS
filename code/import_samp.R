#############################################################
# Calculate p-value using importance sampling
# (If product of weights is very small/large, we can
#  multiply the weight function by a constant so the weights
#  are in some sense centered at 1)
#############################################################
IS.permute <- function(data, grid.points, w.fun=function(x){1}, prms, test.type)
#                       expectations.table=c(), counts.flag, test.type)
{
  n <- dim(data)[1]
  if(missing(grid.points) || isempty(grid.points))
    grid.points = data #  grid.points <- matrix(rnorm(2*n, 0, 1), n, 2) # TEMP DEBUG!!
  if(!('counts.flag' %in% names(prms))) # set default
    counts.flag <- 1  # default: use counts in the statistic (not probabilities)
  if(!('importance.sampling.dist' %in% names(prms))) # set default uniform distribution 
    prms$importance.sampling.dist <- "uniform"
  
  Permutations <- PermutationsIS(w_fun_to_mat(data, w.fun), prms)  # sample permutations 
  P.MCMC <- PermutationsMCMC(w_fun_to_mat(data, w.fun), prms)   # debug for comparison 
  orig.IS.dist <- prms$importance.sampling.dist
  
  prms$importance.sampling.decreasing <- TRUE
  Permutations.D <- PermutationsIS(w_fun_to_mat(data, w.fun), prms)  # sample permutations 
  
  prms$importance.sampling.dist <- "uniform"
  Permutations.U <- PermutationsIS(w_fun_to_mat(data, w.fun), prms)  # sample permutations   
#  inverse.weight = missing(expectations.table) || isempty(expectations.table) # default is using inverse weighting 
  inverse.weight = str_detect(test.type, "inverse_weight")
  if(inverse.weight)
    TrueT <- ComputeStatistic.W(data, grid.points, w.fun, prms$counts.flag)$Statistic # no unique in grid-points 
  else    # new! here we also need to compute the expectatuibs table !!!  
  {
    # Compute expectations.table from permutations   
    expectations.table <- QuarterProbFromPermutations(data, Permutations$P, grid.points) #Permutations
    expectations.MCMC <- QuarterProbFromPermutations(data, P.MCMC$P, grid.points) #Permutations
    plot(expectations.MCMC, expectations.table, main=orig.IS.dist)
    expectations.table <- expectations.MCMC # TEMP FOR DEBUG! DOES MCMC GIVE GOOD EXPECTATIONS AND CHANGES PVALUE?
    print(sum((expectations.table-expectations.MCMC)^2))
    TrueT <- ComputeStatistic(data, grid.points, expectations.table)$Statistic # no unique in grid-points !!
  }
  
  T.b <- matrix(0, prms$B, 1) # statistics under null 
  T.b.U <- matrix(0, prms$B, 1)
  T.b.D <- matrix(0, prms$B, 1)
  T.b.MC <- matrix(0, prms$B, 1)
  for(b in 1:prms$B)
  {
    perm <- Permutations$Permutations[,b] # take sampled permutation sample(n)  # uniformly sampled permutation 
    if(inverse.weight)
      T.b[b] <- ComputeStatistic.W(cbind(data[,1], data[perm,2]), grid.points, w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
    else
      T.b[b] <- ComputeStatistic(cbind(data[,1], data[perm,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
    perm.U <- Permutations.U$Permutations[,b]
    T.b.U[b] <- ComputeStatistic(cbind(data[,1], data[perm.U,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
    perm.D <- Permutations.D$Permutations[,b]
    T.b.D[b] <- ComputeStatistic(cbind(data[,1], data[perm.D,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
    perm.MC <- P.MCMC$Permutations[,b]
    T.b.MC[b] <- ComputeStatistic(cbind(data[,1], data[perm.MC,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
    
  }
  Pvalue = (0+0*Permutations$P.W.IS0 + sum((T.b>=TrueT) * Permutations$P.W.IS)) / 
    (0+0*Permutations$P.W.IS0 + sum(Permutations$P.W.IS))  # New: give weight 1 to the identity permutation 
  plot(Permutations$log.P.W-Permutations.U$log.P.W0, T.b, 
       xlim = c(min(Permutations$log.P.W-Permutations$log.P.W0, Permutations.D$log.P.W-Permutations$log.P.W0), 
                max(Permutations$log.P.W-Permutations$log.P.W0, Permutations.U$log.P.W-Permutations$log.P.W0, P.MCMC$log.P.W-Permutations$log.P.W0)), 
       ylim = c(0, max(T.b, T.b.U, T.b.D, TrueT)), main=orig.IS.dist, ylab="Statistic")
  points(Permutations.D$log.P.W-Permutations.D$log.P.W0, T.b.D, col="green") 
  points(Permutations.U$log.P.W-Permutations.U$log.P.W0, T.b.U, col="blue") 
  points(P.MCMC$log.P.W-Permutations$log.P.W0, T.b.MC, col="orange") 
  points(0, TrueT, col="red", pch=10)
  legend(1, 1200, legend=c(orig.IS.dist, 'monotone.d', 'uniform', 'MCMC', 'True'), col=c('black', 'green', 'blue', 'orange', 'red'), lwd=c(2,2,2,2,2), cex=0.6)
  return(list(Pvalue= Pvalue, TrueT=TrueT, statistics.under.null=T.b))
  # WTF?
}


# New: Importance sampling for permutations
# orms$importance.sampling.dist - determines the distribution 
PermutationsIS <- function(w.mat, prms) # burn.in=NA, Cycle=NA)  # New: allow non-default burn-in 
{ 
    n <- dim(w.mat)[1];
    P <- matrix(0, n, n) # New! matrix with P[i]=j estimate

    # Set mcmc default sampling parameters 
    if(!('B' %in% names(prms))) # set default
      prms$B <- 1000
    if(!('importance.sampling.dist' %in% names(prms))) # set default uniform distribution 
      prms$importance.sampling.dist <- "uniform"
    Permutations = matrix(0, n, prms$B)
    
    switch(prms$importance.sampling.dist,
           'uniform'={
             for(b in 1:prms$B)  # sample permutations 
               Permutations[,b]=sample(n);
             P.IS <- zeros(prms$B,1)  # get importance probabilities (up to a constant)
             P.IS0 <- 0 # these are on log-scale 
           },    
           'monotone.w'={
             if(!('importance.sampling.decreasing' %in% names(prms)))
               prms$importance.sampling.decreasing = FALSE  # default is increasing! 
#             w.order <- order(  rowVars(log(w.mat)), decreasing=TRUE  )  # order by variance: first set as high the ones with high variance!! 
             w.order <- order(rowSums(w.mat), decreasing=prms$importance.sampling.decreasing)  # order x values based on sum of w. We wnat Small with Large! (at least for w(x,y)=x+y)
             P.IS <- zeros(prms$B, 1)
             P.IS0 <- 0
             for(j in 1:n) # Caculate prob. of identity under importance sampling 
             {
               weights <- w.mat[w.order[j],]
               if(j>1)
                 weights[w.order[1:(j-1)]] <- 0    
#               weights <- weights / sum(weights)
               P.IS0 <- P.IS0 + log(weights[w.order[j]] / sum(weights))
             }
                 
             for(b in 1:prms$B)  # sample permutations 
             {
               for(j in 1:n) # sample permutation according to w
               {
                 weights <- w.mat[w.order[j],]
                 if(j>1) # remove the ones we already occupied
                   weights[Permutations[w.order[1:(j-1)]  , b]] <- 0                   
                 weights <- weights / sum(weights)
                 Permutations[w.order[j],b] = sample(n, 1, prob = weights)
                 P.IS[b] <- P.IS[b] + log(weights[Permutations[w.order[j],b]]) # Need to update also P_IS here
               }
             }
           }
    ) # end switch on importance sampling distribution 

    # Next, compute importance weights, from P_I and P_W, and then use them to estimate P_ij 
    log.w.mat <- log(w.mat)
    # Comptue also the id permutation weight 
    P.W.IS0 <- -P.IS0 + sum(diag(log.w.mat))

    P.W.IS <- -P.IS # Computing P_W / P_IS 
    for(b in 1:prms$B)
      P.W.IS[b] <- P.W.IS[b] + sum(log.w.mat[cbind(1:n, Permutations[,b])])
#    print("P.W.IS.LOG:")
#    print(t(P.W.IS))
    
        max.log <- max(max(P.W.IS), P.W.IS0)
    P.W.IS <- exp(P.W.IS - max.log)  # compute P_W/P_I (up to a constant)
    P.W.IS0 <- exp(P.W.IS0 - max.log) # forgot exponent here 

#    print("P.IS:")
#    print(P.IS)
    
#    print("P.W.IS:")
#    print(t(P.W.IS))
#    print("P.W.IS0:")
#    print(P.W.IS0)
    log.P.W0 <- sum(diag(log.w.mat))
    log.P.W <- zeros(prms$B, 1)
    for(b in 1:prms$B)
      log.P.W[b] <- sum(log.w.mat[cbind(1:n, Permutations[,b])])
    max.w <- max(max(log.P.W), log.P.W0)
#    P.W <- exp(P.W - max.w)  # compute P_W/P_I (up to a constant)
#    P.W0 <- exp(P.W0 - max.w) # forgot exponent here         
    print(paste0(prms$importance.sampling.dist, " rand-perms: P_W=", sum(exp(log.P.W)), ", P_W/P_IS=", sum(P.W.IS), " id: P_W=", log.P.W0, ", P_W/P_IS=", P.W.IS0))      

    for(b in 1:prms$B)  # next, compute expectations P[i,j]
      P[cbind(1:n, Permutations[,b])] =  P[cbind(1:n, Permutations[,b])] + P.W.IS[b] 
    P[cbind(1:n, 1:n)] <- P[cbind(1:n, 1:n)] + 0*P.W.IS0 # temp: ignore contribution from identity
    P <- P / (0*P.W.IS0 + sum(P.W.IS)) # normalize P    
      
#    print(P[1:5,1:5])
    
#    print(t(rowSums(P)))
#    print(t(colSums(P)))
    return(list(Permutations=Permutations, P=P, P.W.IS=P.W.IS, P.W.IS0=P.W.IS0, log.P.W=log.P.W, log.P.W0=log.P.W0, max.w=max.w)) # New: return also P, a matrix with Pr(pi(i)=j)
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


