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

  Permutations <- PermutationsIS(w_fun_to_mat(data, w.fun), prms)  # sample permutations 
#  P.MCMC <- PermutationsMCMC(w_fun_to_mat(data, w.fun), prms)   # debug for comparison 
  
  
#  inverse.weight = missing(expectations.table) || isempty(expectations.table) # default is using inverse weighting 
  inverse.weight = str_detect(test.type, "inverse_weight")
  if(inverse.weight)
    TrueT <- ComputeStatistic.W(data, grid.points, w.fun, prms$counts.flag)$Statistic # no unique in grid-points 
  else    # new! here we also need to compute the expectatuibs table !!!  
  {
    # Compute expectations.table from permutations   
    expectations.table <- QuarterProbFromPermutations(data, Permutations$P, grid.points) #Permutations
  #  expectations.MCMC <- QuarterProbFromPermutations(data, P.MCMC$P, grid.points) #Permutations
    TrueT <- ComputeStatistic(data, grid.points, expectations.table)$Statistic # no unique in grid-points !!
  }
  
  T.b <- matrix(0, prms$B, 1) # statistics under null 
  for(b in 1:prms$B)
  {
    perm <- Permutations$Permutations[,b] # take sampled permutation sample(n)  # uniformly sampled permutation 
    if(inverse.weight)
      T.b[b] <- ComputeStatistic.W(cbind(data[,1], data[perm,2]), grid.points, w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
    else
      T.b[b] <- ComputeStatistic(cbind(data[,1], data[perm,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
  }
  Pvalue = (1+0*Permutations$P.W.IS0 + sum((T.b>=TrueT) * Permutations$P.W.IS)) / 
    (1+0*Permutations$P.W.IS0 + sum(Permutations$P.W.IS))  # New: give weight 1 to the identity permutation 
  return(list(Pvalue= Pvalue, TrueT=TrueT, statistics.under.null=T.b))
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
             P_IS <- zeros(prms$B,1)  # get importance probabilities (up to a constant)
             P_IS0 <- 0 # these are on log-scale 
           },    
           'monotone.w'={
##             w.order <- order(  rowVars(log(w.mat)), decreasing=TRUE  )  # order by variance: first set as high the ones with high variance!! 
             w.order <- order(rowSums(w.mat), decreasing=FALSE)  # order x values based on sum of w. We wnat Small with Large! (at least for w(x,y)=x+y)
             P_IS <- zeros(prms$B, 1)
             P_IS0 <- 0
             for(j in 1:n) # Caculate prob. of identity under importance sampling 
             {
               weights <- w.mat[w.order[j],]
               if(j>1)
                 weights[w.order[1:(j-1)]] <- 0    
               weights <- weights / sum(weights)
               P_IS0 <- P_IS0 + log(weights[w.order[j]])
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
                 P_IS[b] <- P_IS[b] + log(weights[Permutations[w.order[j],b]]) # Need to update also P_IS here
               }
             }
           }
    ) # end switch on importance sampling distribution 

    # Next, compute importance weights, from P_I and P_W, and then use them to estimate P_ij 
    log.w.mat <- log(w.mat)
    # Comptue also the id permutation weight 
    P.W.IS0 <- -P_IS0 + sum(diag(log.w.mat))

    P.W.IS <- -P_IS # Computing P_W / P_IS 
    for(b in 1:prms$B)
      P.W.IS[b] <- P.W.IS[b] + sum(log.w.mat[cbind(1:n, Permutations[,b])])
    max.log <- max(max(P.W.IS), P.W.IS0)
    P.W.IS <- exp(P.W.IS - max.log)  # compute P_W/P_I (up to a constant)
    P.W.IS0 <- exp(P.W.IS0 - max.log) # forgot exponent here 

    P.W0 <- sum(diag(log.w.mat))
    P.W <- zeros(prms$B, 1)
    for(b in 1:prms$B)
      P.W[b] <- sum(log.w.mat[cbind(1:n, Permutations[,b])])
    max.w <- max(max(P.W), P.W0)
    P.W <- exp(P.W - max.w)  # compute P_W/P_I (up to a constant)
    P.W0 <- exp(P.W0 - max.w) # forgot exponent here         
    print(paste0(prms$importance.sampling.dist, " rand-perms: P_W=", sum(P.W)/P.W0, ", P_W/P_IS=", sum(P.W.IS)/P.W.IS0))      

    for(b in 1:prms$B)  # next, compute expectations P[i,j]
      P[cbind(1:n, Permutations[,b])] =  P[cbind(1:n, Permutations[,b])] + P.W.IS[b] 
    P[cbind(1:n, 1:n)] <- P[cbind(1:n, 1:n)] + 0*P.W.IS0 # temp: ignore contribution from identity
    P <- P / (0*P.W.IS0 + sum(P.W.IS)) # normalize        
      
    return(list(Permutations=Permutations, P=P, P.W.IS=P.W.IS, P.W.IS0=P.W.IS0)) # New: return also P, a matrix with Pr(pi(i)=j)
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


