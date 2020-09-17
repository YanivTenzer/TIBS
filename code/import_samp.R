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
  Permutations.MC <- PermutationsMCMC(w_fun_to_mat(data, w.fun), prms)   # debug for comparison 
  orig.IS.dist <- prms$importance.sampling.dist
  
  prms$importance.sampling.decreasing <- TRUE
  Permutations.D <- PermutationsIS(w_fun_to_mat(data, w.fun), prms)  # sample permutations 
  
  prms$importance.sampling.dist <- "uniform"
  Permutations.U <- PermutationsIS(w_fun_to_mat(data, w.fun), prms)  # sample permutations   
  
  prms$importance.sampling.dist <- "match.w"
  Permutations.MATCH <- PermutationsIS(w_fun_to_mat(data, w.fun), prms)  # sample permutations   
  plot.flag=1   
  #  inverse.weight = missing(expectations.table) || isempty(expectations.table) # default is using inverse weighting 
  inverse.weight = str_detect(test.type, "inverse_weight")
  if(inverse.weight)
    TrueT <- ComputeStatistic.W(data, grid.points, w.fun, prms$counts.flag)$Statistic # no unique in grid-points 
  else    # new! here we also need to compute the expectatuibs table !!!  
  {
    # Compute expectations.table from permutations   
    expectations.table <- QuarterProbFromPermutations(data, Permutations$P, grid.points) #Permutations
    expectations.U <- QuarterProbFromPermutations(data, Permutations.U$P, grid.points) #Permutations
    expectations.D <- QuarterProbFromPermutations(data, Permutations.D$P, grid.points) #Permutations
    expectations.MCMC <- QuarterProbFromPermutations(data, Permutations.MC$P, grid.points) #Permutations
    expectations.MATCH <- QuarterProbFromPermutations(data, Permutations.MATCH$P, grid.points) #Permutations
    if(plot.flag)
    {
    plot(expectations.MCMC, expectations.table, main=orig.IS.dist)
    lines(sort(expectations.MCMC), sort(expectations.MCMC), col="red", lwd=2)
    }    
    expect.str <- "Expectations by Monotone.MATCH" # This choice affects pvalue greatly 
    expectations.table <- expectations.MATCH # MCMC # TEMP FOR DEBUG! DOES MCMC GIVE GOOD EXPECTATIONS AND CHANGES PVALUE?
    ###    print(sum((expectations.table-expectations.MCMC)^2))
    TrueT <- ComputeStatistic(data, grid.points, expectations.table)$Statistic # no unique in grid-points !!
  }
  
  T.b <- matrix(0, prms$B, 1) # statistics under null 
  T.b.U <- matrix(0, prms$B, 1)
  T.b.D <- matrix(0, prms$B, 1)
  T.b.MC <- matrix(0, prms$B, 1)
  T.b.MATCH <- matrix(0, prms$B, 1)
  
  
  for(b in 1:prms$B)
  {
    perm <- Permutations$Permutations[,b] # take sampled permutation sample(n)  # uniformly sampled permutation 
    if(inverse.weight)
    {
      expect.str <- "inverse_weight" # This choice affects pvalue greatly 
      T.b[b] <- ComputeStatistic.W(cbind(data[,1], data[perm,2]), grid.points, w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
      perm.U <- Permutations.U$Permutations[,b]
      T.b.U[b] <- ComputeStatistic.W(cbind(data[,1], data[perm.U,2]), grid.points,  w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
      perm.D <- Permutations.D$Permutations[,b]
      T.b.D[b] <- ComputeStatistic.W(cbind(data[,1], data[perm.D,2]), grid.points,  w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
      perm.MC <- Permutations.MC$Permutations[,b]
      T.b.MC[b] <- ComputeStatistic.W(cbind(data[,1], data[perm.MC,2]), grid.points,  w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
      perm.MATCH <- Permutations.MATCH$Permutations[,b]
      T.b.MATCH[b] <- ComputeStatistic.W(cbind(data[,1], data[perm.MATCH,2]), grid.points,  w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
      
    }
    else
    {
      T.b[b] <- ComputeStatistic(cbind(data[,1], data[perm,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
      perm.U <- Permutations.U$Permutations[,b]
      T.b.U[b] <- ComputeStatistic(cbind(data[,1], data[perm.U,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
      perm.D <- Permutations.D$Permutations[,b]
      T.b.D[b] <- ComputeStatistic(cbind(data[,1], data[perm.D,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
      perm.MC <- Permutations.MC$Permutations[,b]
      T.b.MC[b] <- ComputeStatistic(cbind(data[,1], data[perm.MC,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
      perm.MATCH <- Permutations.MATCH$Permutations[,b]
      T.b.MATCH[b] <- ComputeStatistic(cbind(data[,1], data[perm.MATCH,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
    }
  }
  
  CV2 = var(Permutations$P.W.IS) / mean(Permutations$P.W.IS)^2 # new: add Coefficient of variation for the importance weights for diagnostics
  CV.D2 = var(Permutations.D$P.W.IS) / mean(Permutations.D$P.W.IS)^2 # new: add Coefficient of variation for the importance weights for diagnostics
  CV.U2 = var(Permutations.U$P.W.IS) / mean(Permutations.U$P.W.IS)^2 # new: add Coefficient of variation for the importance weights for diagnostics
  CV.W2 = var(Permutations$P.W.IS) / mean(Permutations$P.W.IS)^2 # new: add Coefficient of variation for the importance weights for diagnostics
  CV.MC2 = 1 # var(Permutations.MC$P.W.IS) / mean(Permutations.MC$P.W.IS)^2 # new: add Coefficient of variation for the importance weights for diagnostics
  CV.MATCH2 = var(Permutations.MATCH$P.W.IS) / mean(Permutations.MATCH$P.W.IS)^2 # new: add Coefficient of variation for the importance weights for diagnostics
  

  
      
  Pvalue = (Permutations$P.W.IS0 + sum((T.b>=TrueT) * Permutations$P.W.IS)) / 
    (Permutations$P.W.IS0 + sum(Permutations$P.W.IS))  # New: give weight 1 to the identity permutation 
  Pval.W = (Permutations$P.W.IS0 + sum((T.b>=TrueT) * Permutations$P.W.IS)) / 
    (Permutations$P.W.IS0 + sum(Permutations$P.W.IS))  # New: give weight 1 to the identity permutation 
  Pval.D = (Permutations.D$P.W.IS0 + sum((T.b.D>=TrueT) * Permutations.D$P.W.IS)) / 
    (Permutations.D$P.W.IS0 + sum(Permutations.D$P.W.IS))  # New: give weight 1 to the identity permutation 
  Pval.U = (Permutations.U$P.W.IS0 + sum((T.b.U>=TrueT) * Permutations.U$P.W.IS)) / 
    (Permutations.U$P.W.IS0 + sum(Permutations.U$P.W.IS))  # New: give weight 1 to the identity permutation 
  Pval.MC = (1 + sum((T.b.MC>=TrueT) )) / (B+1)  # New: give weight 1 to the identity permutation 
  Pval.MATCH = (Permutations.MATCH$P.W.IS0 + sum((T.b.MATCH>=TrueT) * Permutations.MATCH$P.W.IS)) / 
    (Permutations.MATCH$P.W.IS0 + sum(Permutations.MATCH$P.W.IS))  # New: give weight 1 to the identity permutation 
  Pval.MATCH0 = (0*Permutations.MATCH$P.W.IS0 + sum((T.b.MATCH>=TrueT) * Permutations.MATCH$P.W.IS)) / 
    (0*Permutations.MATCH$P.W.IS0 + sum(Permutations.MATCH$P.W.IS))  # New: give weight 1 to the identity permutation 
  Pval.MATCH.no.weights = mean(T.b.MATCH>=TrueT)
  
  legend.vec <- apply(rbind(c(orig.IS.dist, 'monotone.d', 'uniform', 'MCMC', 'Match', 'True'), rep(" CV=", 6),
                      round(c(CV2, CV.D2, CV.U2, CV.MC2, CV.MATCH2, 0), 3), rep(" Pval=", 6),
                            round(c(Pval.W, Pval.D, Pval.U, Pval.MC, Pval.MATCH, 0), 3)), 2, paste0, collapse="")
#  for(j in c(1:length(legend.vec)))
#    legend.vec[j] <- c(legend.vec[j], ", pval=", pval.method[j], ", CV=", cv.method[j])
  
  plot(Permutations.MATCH$log.P.W, Permutations.MATCH$log.P.IS)
  points(Permutations.MATCH$log.P.W[T.b.MATCH>=TrueT], Permutations.MATCH$log.P.IS[T.b.MATCH>=TrueT], col="green")
  points(Permutations.MATCH$log.P.W0, Permutations.MATCH$log.P.IS0, col="red", pch=19, cex=2 )
  points(Permutations.MATCH$log.P.W[89], Permutations.MATCH$log.P.IS[89], col="blue")
  
  plot(Permutations.MATCH$log.P.W - Permutations.MATCH$log.P.IS)
  points(which(T.b.MATCH>=TrueT), Permutations.MATCH$log.P.W[T.b.MATCH>=TrueT] - Permutations.MATCH$log.P.IS[T.b.MATCH>=TrueT], col="green")
  points(Permutations.MATCH$log.P.W0 - Permutations.MATCH$log.P.IS0, col="red", pch=19, cex=2 )
  
  if(plot.flag)
  {
    y.lim <- c(0, max(T.b, T.b.U, T.b.D, T.b.MC, T.b.MATCH, TrueT))
    x.lim <- c(min(Permutations$log.P.W-Permutations$log.P.W0, Permutations.D$log.P.W-Permutations$log.P.W0, 
                   Permutations.MATCH$log.P.W-Permutations$log.P.W0), 
               max(Permutations$log.P.W-Permutations$log.P.W0, Permutations.U$log.P.W-Permutations$log.P.W0, 
                   Permutations.MC$log.P.W-Permutations$log.P.W0, Permutations.MATCH$log.P.W-Permutations$log.P.W0))
    
    plot(Permutations$log.P.W-Permutations.U$log.P.W0, T.b, 
         xlim = x.lim, ylim = y.lim, main=expect.str, ylab="Statistic") # with/without inverse weighting stat 
    points(Permutations.D$log.P.W-Permutations.D$log.P.W0, T.b.D, col="green") 
    points(Permutations.U$log.P.W-Permutations.U$log.P.W0, T.b.U, col="blue") 
    points(Permutations.MC$log.P.W-Permutations$log.P.W0, T.b.MC, col="orange") 
    points(Permutations.MATCH$log.P.W-Permutations$log.P.W0, T.b.MATCH, col="pink") 
    points(0, TrueT, col="red", pch=19, cex=2)
      
    
    legend(-10, y.lim[2], legend=legend.vec, col=c('black', 'green', 'blue', 'orange', 'pink', 'red'), lwd=c(2,2,2,2,2), cex=0.6)
  }
  xxx = 132
  
  return(list(Pvalue= Pvalue, TrueT=TrueT, CV2=CV2, statistics.under.null=T.b, permuted.data = cbind(data[,1], data[Permutations$Permutations[,1],2])))
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
  log.w.mat <- log(w.mat)
  
  switch(prms$importance.sampling.dist,
         'uniform'={
           for(b in 1:prms$B)  # sample permutations 
             Permutations[,b]=sample(n);
           log.P.IS <- zeros(prms$B,1)  # get importance probabilities (up to a constant)
           log.P.IS0 <- 0 # these are on log-scale 
         },    
         'monotone.w'={
           if(!('importance.sampling.decreasing' %in% names(prms)))
             prms$importance.sampling.decreasing = FALSE  # default is increasing! 
           #             w.order <- order(  rowVars(log(w.mat)), decreasing=TRUE  )  # order by variance: first set as high the ones with high variance!! 
           w.order <- order(rowSums(w.mat), decreasing=prms$importance.sampling.decreasing)  # order x values based on sum of w. We wnat Small with Large! (at least for w(x,y)=x+y)
           log.P.IS <- zeros(prms$B, 1)
           log.P.IS0 <- 0
           for(j in 1:n) # Caculate prob. of identity under importance sampling 
           {
             weights <- w.mat[w.order[j],]
             if(j>1)
               weights[w.order[1:(j-1)]] <- 0    
             #               weights <- weights / sum(weights)
             log.P.IS0 <- log.P.IS0 + log(weights[w.order[j]] / sum(weights))
           }
           
           for(b in 1:prms$B)  # sample permutations 
           {
             for(j in 1:n) # sample permutation according to w
             {
               weights <- w.mat[w.order[j],]
               if(j>1) # remove the ones we already occupied
                 weights[Permutations[w.order[1:(j-1)]  , b]] <- 0                   
               #                print("weights not-normalized")
               #                                  print(weights)
               
               weights <- weights / sum(weights)
               #                 print("weights normalized")
               #                 print(weights)
               #                 print(c(j, n))
               Permutations[w.order[j],b] = sample(n, 1, prob = weights)
               log.P.IS[b] <- log.P.IS[b] + log(weights[Permutations[w.order[j],b]]) # Need to update also P_IS here
             }
           }
         }, 
         'match.w'={  # here the goal is to sample a permutation with P_W similar to the ID permutation
           w.order <- order(  rowVars(log(w.mat)), decreasing=TRUE  )  # order by variance: first set as high the ones with high variance!! 
           log.w.sum.id.vec <- cumsum(diag(log.w.mat)[w.order]) # we want to get this sum 
           log.P.IS <- zeros(prms$B, 1)
#           P.IS0 <- -sum(log(1:n)) # temp, use uniform distribution 
           sigma <- 2 # choose width of Gaussian for importance sampling. Can change to be preplexity
           preplex <- 8 #  effective number of possible choices at each step (not implemented yet)
           log.P.IS0 <- 0
           log.w.sum.rand <- 0
           for(j in 1:n) # Caculate prob. of identity under importance sampling 
           {
             weights <- exp(-(log.w.mat[w.order[j],] + log.w.sum.rand - log.w.sum.id.vec[j] )^2/(2*sigma^2))
             if(j>1)
               weights[w.order[1:(j-1)]] <- 0    
             weights <- weights / sum(weights)
             log.P.IS0 <- log.P.IS0 + log(weights[w.order[j]])
             log.w.sum.rand <- log.w.sum.rand + log.w.mat[w.order[j], w.order[j]]
           }
           
           for(b in 1:prms$B)  # sample permutations 
           {
             log.w.sum.rand <- 0
             for(j in 1:n) # sample permutation according to w
             {
               weights <- exp(-(log.w.mat[w.order[j],] + log.w.sum.rand - log.w.sum.id.vec[j] )^2/(2*sigma^2)) # penalize deviations from probability
               if(j>1) # remove the ones we already occupied
                 weights[Permutations[w.order[1:(j-1)] , b]] <- 0                   
               weights <- weights / sum(weights)
               #                print(weights)
               #                print(paste0("j=", j))
               Permutations[w.order[j],b] = sample(n, 1, prob = weights)
               log.P.IS[b] <- log.P.IS[b] + log(weights[Permutations[w.order[j],b]]) # Need to update also P_IS here
               log.w.sum.rand <- log.w.sum.rand + log.w.mat[w.order[j],Permutations[w.order[j],b]]
               #                print(paste0("sum.rand=", log.w.sum.rand, " sum.id=",  log.w.sum.id.vec[j], " diff=", log.w.sum.rand-log.w.sum.id.vec[j]))
             }
             ###              print(paste0("log.w.sum.rand=", log.w.sum.rand))
             ###              print(paste0("P_W=", sum(log.w.mat[cbind(1:n, Permutations[,b])])))
             
           }     
           
           
         }
  ) # end switch on importance sampling distribution 
  
  # Next, compute importance weights, from P_I and P_W, and then use them to estimate P_ij 
  # Comptue also the id permutation weight 
  P.W.IS0 <- -log.P.IS0 + sum(diag(log.w.mat))
  
  P.W.IS <- -log.P.IS # Computing P_W / P_IS 
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
  ###    print(paste0(prms$importance.sampling.dist, " rand-perms: P_W=", sum(exp(log.P.W)), ", P_W/P_IS=", sum(P.W.IS), " id: P_W=", log.P.W0, ", P_W/P_IS=", P.W.IS0))      
  
  for(b in 1:prms$B)  # next, compute expectations P[i,j]
    P[cbind(1:n, Permutations[,b])] =  P[cbind(1:n, Permutations[,b])] + P.W.IS[b] 
  P[cbind(1:n, 1:n)] <- P[cbind(1:n, 1:n)] + 0*P.W.IS0 # temp: ignore contribution from identity
  P <- P / (0*P.W.IS0 + sum(P.W.IS)) # normalize P    
  
  #    print(P[1:5,1:5])
  
  #    print(t(rowSums(P)))
  #    print(t(colSums(P)))
  return(list(Permutations=Permutations, P=P, P.W.IS=P.W.IS, P.W.IS0=P.W.IS0, log.P.W=log.P.W, log.P.W0=log.P.W0, max.w=max.w, 
              log.P.IS=log.P.IS, log.P.IS0=log.P.IS0)) # New: return also P, a matrix with Pr(pi(i)=j)
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


