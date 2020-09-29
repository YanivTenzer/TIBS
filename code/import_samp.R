#############################################################
# Calculate p-value using importance sampling
# (If product of weights is very small/large, we can
#  multiply the weight function by a constant so the weights
#  are in some sense centered at 1)
#############################################################
IS.permute <- function(data, grid.points, w.fun=function(x){1}, prms, test.stat)  #                       expectations.table=c(), counts.flag, test.type)
{
  epsilon <- 0.00000000000001
  col.vec <- c("orange", "cyan", "navy", "dodgerblue3", "green", "red", "yellow",  "gray", "pink",  "purple", "cyan")
  pch.vec <- c(20, 1, 4, 13, 3, 17 )
  
  
  n <- dim(data)[1]
  if(missing(grid.points) || isempty(grid.points))
    grid.points = data #  grid.points <- matrix(rnorm(2*n, 0, 1), n, 2) # TEMP DEBUG!!
  if(!('counts.flag' %in% names(prms))) # set default
    counts.flag <- 1  # default: use counts in the statistic (not probabilities)
  if(!('importance.sampling.dist' %in% names(prms))) # set default KouMcculough.w distribution (works also for truncation!)
    prms$importance.sampling.dist <- "KouMcculough.w"
  
  if(!('diagnostic.plot' %in% names(prms))) # set default
    prms$diagnostic.plot <- 0
  
  w.mat <- w_fun_to_mat(data, w.fun, prms)
#  prob.typical <- estimate_log_prob_typical(w.mat)
  
  orig.IS.dist <- prms$importance.sampling.dist
  Permutations <- PermutationsIS(w.mat, prms, data)  # sample permutations 
  Permutations.cpp <- PermutationsIS_rcpp(w.mat, prms, data) # TEMP - compare with cpp 
  if(prms$diagnostic.plot)
  {
    IS.methods <- c("KouMcculough.w", "uniform", "monotone.w", "monotone.grid.w", "MCMC") #   "sqrt.w", ) "match.w", 
    IS.methods <- setdiff(IS.methods, prms$importance.sampling.dist)
    num.IS <- length(IS.methods)
    perm.IS <- c()
    T.b.IS <- zeros(num.IS, prms$B)
    
    for(i in c(1:num.IS))
    {
      if(IS.methods[i] == "MCMC")
      {
        perm.IS[[i]] <- PermutationsMCMC(w_fun_to_mat(data, w.fun, prms), prms)   # debug for comparison 
        
      } else
      {
        prms$importance.sampling.dist <- IS.methods[i]
        if(IS.methods[i] == "monotone.w")
          prms$importance.sampling.decreasing <- FALSE
        else
          prms$importance.sampling.decreasing <- FALSE
        perm.IS[[i]] <- PermutationsIS(w_fun_to_mat(data, w.fun, prms), prms) 
      }
    }
    
  }
  
  inverse.weight = str_detect(test.stat, "inverse")
  if(inverse.weight)
    TrueT <- ComputeStatistic.W(data, grid.points, w.fun, prms$counts.flag)$Statistic # no unique in grid-points 
  else    # new! here we also need to compute the expectatuibs table !!!  
  {
    # Compute expectations.table from permutations   
    expectations.table <- QuarterProbFromPermutations(data, Permutations$P, grid.points) #Permutations
    
    if(prms$diagnostic.plot)
    {
      expectations.IS <- c()
      for(i in c(1:num.IS))
        expectations.IS[[i]] <- QuarterProbFromPermutations(data, perm.IS[[i]]$P, grid.points) #Permutations
      i.MCMC <- which(IS.methods=="MCMC")
      plot(expectations.IS[[i.MCMC]], expectations.table, main=orig.IS.dist)
      lines(sort(expectations.IS[[i.MCMC]]), sort(expectations.IS[[i.MCMC]]), col="red", lwd=2)
    }
    TrueT <- ComputeStatistic(data, grid.points, expectations.table)$Statistic # no unique in grid-points !!
    
  }    
  #    expect.str <- "Expectations by Monotone.MATCH" # This choice affects pvalue greatly 
  #    expectations.table <- expectations.MATCH # MCMC # TEMP FOR DEBUG! DOES MCMC GIVE GOOD EXPECTATIONS AND CHANGES PVALUE?
  ###    print(sum((expectations.table-expectations.MCMC)^2))
  
  T.b <- matrix(0, prms$B, 1) # statistics under null 
  # T.b.U <- matrix(0, prms$B, 1)
  # T.b.D <- matrix(0, prms$B, 1)
  # T.b.MC <- matrix(0, prms$B, 1)
  # T.b.MATCH <- matrix(0, prms$B, 1)
  
  
  for(b in 1:prms$B)
  {
    
    perm <- Permutations$Permutations[,b] # take sampled permutation sample(n)  # uniformly sampled permutation 
    if(inverse.weight)
    {
      expect.str <- "Inverse Weighted Hoeffding Null Statistics" # This choice affects pvalue greatly 
      T.b[b] <- ComputeStatistic.W(cbind(data[,1], data[perm,2]), grid.points, w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
      if(prms$diagnostic.plot)
      {
        
        for(i in c(1:num.IS))
        {
          perm <- perm.IS[[i]]$Permutations[,b]
          T.b.IS[i,b] <- ComputeStatistic.W(cbind(data[,1], data[perm,2]), grid.points, w.fun, prms$counts.flag)$Statistic # grid depends on permuted data
        }
        
      }
    } else
    {
      expect.str <- paste0("Expectations by ", orig.IS.dist)
      T.b[b] <- ComputeStatistic(cbind(data[,1], data[perm,2]), grid.points, expectations.table)$Statistic # grid depends on permuted data
      if(prms$diagnostic.plot)
      {
        for(i in c(1:num.IS))
        {
          perm <- perm.IS[[i]]$Permutations[,b]
          T.b.IS[i,b] <- ComputeStatistic(cbind(data[,1], data[perm,2]), grid.points, expectations.table)$Statistic # here use the expectation table of original ! 
        }
      }
    }
  }  # end for on permutations
  CV = std(Permutations$P.W.IS) / mean(Permutations$P.W.IS) # new: add Coefficient of variation for the importance weights for diagnostics
  
  if(prms$include.ID==2) # new: count the ID permutation as 1, all others as P_W(pi)/P_IS(pi). We don't know the normalizing constant so take max (P_W(pi)/P_IS(pi))
    Pvalue = (max(Permutations$P.W.IS) + sum((T.b>=TrueT) * Permutations$P.W.IS)) / 
    (max(Permutations$P.W.IS) + sum(Permutations$P.W.IS))  # New: give weight P_W(pi)/P_IS(pi) to the identity permutation 
  else # here inclue.ID is 0 or 1    
    Pvalue = (Permutations$P.W.IS0*prms$include.ID + sum((T.b>=TrueT) * Permutations$P.W.IS)) / 
      (Permutations$P.W.IS0*prms$include.ID + sum(Permutations$P.W.IS))  # New: give weight 1 to the identity permutation 
  Pvalue0 = (sum((T.b>=TrueT) * Permutations$P.W.IS)) /   (sum(Permutations$P.W.IS))  # New: give weight 1 to the identity permutation 
  
  
#  if(Pvalue < 0.01)
#  {
#    xxxx = 31245435
#    print("low pval - why?")
#  }
  
    
    
  if(prms$diagnostic.plot)
  {

    x.lim <- range(c(Permutations$log.P.W, Permutations$log.P.W0))  
    y.lim <- range(c(Permutations$log.P.IS, Permutations$log.P.IS0)) 
    plot(Permutations$log.P.W, Permutations$log.P.IS, col="black", pch=20, xlim=x.lim, ylim=y.lim, xlab="log(P_w)", ylab="log(P_IS)", main=orig.IS.dist)
    points(Permutations$log.P.W[T.b>=TrueT], Permutations$log.P.IS[T.b>=TrueT], col=col.vec[5], pch=19)
    
    #  points(Permutations$log.P.W[J], Permutations$log.P.IS[J], col="green", pch=5, cex=2)
    points(Permutations$log.P.W0, Permutations$log.P.IS0, col="blue", pch=9, cex=2)
    
    plot(Permutations$log.P.W - Permutations$log.P.IS, T.b, col="black", main=orig.IS.dist)
    lines(range(Permutations$log.P.W - Permutations$log.P.IS), c(TrueT, TrueT), col="red")      
    points(Permutations$log.P.W0 - Permutations$log.P.IS0, TrueT, col="red", pch=9, cex=2)      
    
    
    CV.IS <- zeros(num.IS,1)
    PVAL.IS <- zeros(num.IS,1)
    PVAL.IS0 <- zeros(num.IS,1)
    for(i in c(1:num.IS))
    {
      if(IS.methods[i] == "MCMC")
      {
        CV.IS[i] = 1
        PVAL.IS[i] = (prms$include.ID + sum(T.b.IS[i,]>=TrueT) ) /  (prms$include.ID+prms$B)
      } else
      {
        CV.IS[i] =  std(perm.IS[[i]]$P.W.IS) / mean(perm.IS[[i]]$P.W.IS)
        if(prms$include.ID==2) # new: count the ID permutation as 1, all others as P_W(pi)/P_IS(pi). We don't know the normalizing constant so take max (P_W(pi)/P_IS(pi))
          PVAL.IS[i] = (max(perm.IS[[i]]$P.W.IS) + sum((T.b.IS[i,]>=TrueT) * perm.IS[[i]]$P.W.IS)) / 
          (max(perm.IS[[i]]$P.W.IS) + sum(perm.IS[[i]]$P.W.IS))  # New: give weight 1 to the identity permutation 
        else
          PVAL.IS[i] = (perm.IS[[i]]$P.W.IS0*prms$include.ID + sum((T.b.IS[i,]>=TrueT) * perm.IS[[i]]$P.W.IS)) / 
          (perm.IS[[i]]$P.W.IS0*prms$include.ID + sum(perm.IS[[i]]$P.W.IS))  # New: give weight 1 to the identity permutation 
        PVAL.IS0[i] = ( sum((T.b.IS[i,]>=TrueT) * perm.IS[[i]]$P.W.IS)) / 
          ( sum(perm.IS[[i]]$P.W.IS))  # New: give weight 1 to the identity permutation
      }
    }    
    legend.vec <- apply(rbind(c(orig.IS.dist, IS.methods, "$T_0$"), c(rep(" $CV=", num.IS+1), ""),
                              c(round(c(CV, CV.IS), 3), ""), c(rep("$ $Pval=", num.IS+1), ""),
                              c(round(c(Pvalue, PVAL.IS), 3), ""), c(rep("$", num.IS+1), "")), 2, paste0, collapse="")
    
    for(i in c(1:(num.IS+2)))
    {
      legend.vec[i] <- str_remove(legend.vec[i], ".w")
      legend.vec[i] <- str_replace(legend.vec[i], "e\\.g", "e g")
      legend.vec[i] <- TeX(paste0(legend.vec[i]))
    }      
    plot(Permutations$log.P.W, Permutations$log.P.IS, col=col.vec[1])
    for(i in c(1:num.IS))
    {
      if(IS.methods[i] == "MCMC")
      {
        perm.IS[[i]]$P.W.IS = rep(1, prms$B)
        perm.IS[[i]]$P.W.IS0 = 1
        next
      }
      
      plot(perm.IS[[i]]$log.P.W, perm.IS[[i]]$log.P.IS, col="black", pch=20, main=IS.methods[i])
      points(perm.IS[[i]]$log.P.W[T.b.IS[i,]>=TrueT], perm.IS[[i]]$log.P.IS[T.b.IS[i,]>=TrueT], col=col.vec[i+1], pch=20)
      
#      points(perm.IS[[i]]$log.P.W[J], perm.IS[[i]]$log.P.IS[J], col="green", pch=5, cex=2)
      points(perm.IS[[i]]$log.P.W0, perm.IS[[i]]$log.P.IS0, col="blue", pch=9, cex=2)
      
      
      plot(perm.IS[[i]]$log.P.W - perm.IS[[i]]$log.P.IS, T.b.IS[i,], col="black", main=IS.methods[i])
      lines(range(perm.IS[[i]]$log.P.W - perm.IS[[i]]$log.P.IS), c(TrueT, TrueT), col="red")      
    }

    # Supp. method: Plot diff of log(P_W) - log(P_IS)
    x.lim <- range(log(Permutations$P.W.IS) - log(Permutations$P.W.IS0))
    y.lim <- range(c(TrueT, T.b))
    for(i in c(1:num.IS))
    {
      x.lim[1] <- min(x.lim[1], min(log(perm.IS[[i]]$P.W.IS) - log(perm.IS[[i]]$P.W.IS0)))
      x.lim[2] <- max(x.lim[2], max(log(perm.IS[[i]]$P.W.IS) - log(perm.IS[[i]]$P.W.IS0)))
      y.lim[1] <- min(y.lim[1], min(T.b.IS[i,]))
      y.lim[2] <- max(y.lim[2], max(T.b.IS[i,]))
    }      
    plot(log(Permutations$P.W.IS) - log(Permutations$P.W.IS0), T.b, col=col.vec[1], main="Importance weights vs. Statistic", 
         xlim=x.lim, ylim=y.lim, xlab=TeX("$log(\frac{P_W(\\pi_i)}{P_{IS}(\\pi_i)}) - log(\frac{P_W(\\pi_i)}{P_{IS}(\\pi_i)})$"), ylab=TeX("$T_i$"), pch=pch.vec[1], cex=0.8 )
    for(i in c(1:num.IS))
    {
      if(!("log.P.W.IS.minus.P.W.IS0" %in% names(perm.IS[[4]])) || is.na(perm.IS[[i]]$log.P.W.IS.minus.P.W.IS0))
        log.ratio.diff <- log(perm.IS[[i]]$P.W.IS) - log(perm.IS[[i]]$P.W.IS0)
      else
        log.ratio.diff <- perm.IS[[i]]$log.P.W.IS.minus.P.W.IS0
      points(log.ratio.diff, T.b.IS[i,], col=col.vec[i+1], pch=pch.vec[i+1], cex=0.8)
    }
    lines(range(x.lim), c(TrueT, TrueT), col="red")      
    grid(NULL,NULL, lwd=1)
    legend(quantile(x.lim, 0.25), y.lim[2], legend=legend.vec, col=col.vec[1:(num.IS+2)] ,  lwd=c(rep(1, num.IS+2), 10), 
           lty=c(rep(NA, num.IS+1), 1), pch=c(pch.vec[1:(num.IS+1)], NA) ,  cex=rep(0.8), # c(rep(20, num.IS+1), 17)
           box.lwd = 0,box.col = "white",bg = "white") # bty="n")
    
        #    legend.vec <- apply(rbind(c(orig.IS.dist, 'monotone.d', 'uniform', 'MCMC', 'Match', 'True'), rep(" CV=", 6),
    #                        round(c(CV2, CV.D2, CV.U2, CV.MC2, CV.MATCH2, 0), 3), rep(" Pval=", 6),
    #                            round(c(Pval.W, Pval.D, Pval.U, Pval.MC, Pval.MATCH, 0), 3)), 2, paste0, collapse="")
    #  for(j in c(1:length(legend.vec)))
    #    legend.vec[j] <- c(legend.vec[j], ", pval=", pval.method[j], ", CV=", cv.method[j])
    
    #    plot(Permutations.MATCH$log.P.W, Permutations.MATCH$log.P.IS)
    #    points(Permutations.MATCH$log.P.W[T.b.MATCH>=TrueT], Permutations.MATCH$log.P.IS[T.b.MATCH>=TrueT], col="green")
    #    points(Permutations.MATCH$log.P.W0, Permutations.MATCH$log.P.IS0, col="red", pch=19, cex=2 )
    #    points(Permutations.MATCH$log.P.W[89], Permutations.MATCH$log.P.IS[89], col="blue")
    #  
    #    plot(Permutations.MATCH$log.P.W - Permutations.MATCH$log.P.IS)
    #    points(which(T.b.MATCH>=TrueT), Permutations.MATCH$log.P.W[T.b.MATCH>=TrueT] - Permutations.MATCH$log.P.IS[T.b.MATCH>=TrueT], col="green")
    #    points(Permutations.MATCH$log.P.W0 - Permutations.MATCH$log.P.IS0, col="red", pch=19, cex=2 )

    
     
    # Supp. Methods. Figure: Plot log probs vs. statistic        
    y.lim <- c(0, max(max(T.b), max(T.b.IS))) # , T.b.U, T.b.D, T.b.MC, T.b.MATCH, TrueT))
    x.lim <- range(Permutations$log.P.W-Permutations$log.P.W0)
    for(i in c(1:num.IS))
    {
      x.lim[1] <- min(x.lim[1],  min(perm.IS[[i]]$log.P.W-perm.IS[[i]]$log.P.W0))      
      x.lim[2] <- max(x.lim[2],  max(perm.IS[[i]]$log.P.W-perm.IS[[i]]$log.P.W0))      
    }
    plot(Permutations$log.P.W-Permutations$log.P.W0, T.b, col=col.vec[1], pch=pch.vec[1], cex=0.8,
         xlim = x.lim, ylim = y.lim, main=expect.str, xlab=TeX("$log(P_W(\\pi_i))-log(P_{W}(\\pi_0))$"), ylab=TeX("$T_i$")) # with/without inverse weighting stat 
    grid(NULL,NULL, lwd=1)
    for(i in c(1:num.IS))
      points(perm.IS[[i]]$log.P.W-perm.IS[[i]]$log.P.W0, T.b.IS[i,], col=col.vec[i+1], pch=pch.vec[i+1], cex=0.8) 
    points(0, TrueT, col=col.vec[num.IS+2], pch=17, cex=2 ) # plot true point 
    
    legend(quantile(x.lim, 0.25), y.lim[2], legend=legend.vec, col=col.vec[1:(num.IS+2)] ,  lwd=c(rep(1, num.IS+2), 10), 
           lty=rep(NA, num.IS+2), pch=pch.vec[1:(num.IS+2)] ,  cex=rep(0.8), # c(rep(20, num.IS+1), 17)
           box.lwd = 0,box.col = "white",bg = "white") # bty="n")
    
#    prob.typical <- estimate_log_prob_typical(w.mat)
    
    
  } # if diagnostic plots 
  
  return(list(Pvalue= Pvalue, TrueT=TrueT, CV=CV, statistics.under.null=T.b, permuted.data = cbind(data[,1], data[Permutations$Permutations[,1],2])))
} 


# New: Importance sampling for permutations
# orms$importance.sampling.dist - determines the distribution 
# data - (optional). The x,y values 
PermutationsIS <- function(w.mat, prms, data) # burn.in=NA, Cycle=NA)  # New: allow non-default burn-in 
{ 
  epsilon <- 0.000000000000000001
  n <- dim(w.mat)[1];
  P <- matrix(0, n, n) # New! matrix with P[i]=j estimate
  
  # Set mcmc default sampling parameters 
  if(!('B' %in% names(prms))) # set default
    prms$B <- 1000
  if(!('importance.sampling.dist' %in% names(prms))) # set default uniform distribution 
    prms$importance.sampling.dist <- "uniform"
  Permutations = matrix(0, n, prms$B)
  log.w.mat <- log(w.mat)
  
  
  #  print(paste0("Sample IS Dist: ", prms$importance.sampling.dist))
  
  switch(prms$importance.sampling.dist,
         'uniform'={
           for(b in 1:prms$B)  # sample permutations 
             Permutations[,b]=sample(n)
           log.P.IS <- zeros(prms$B,1)  # get importance probabilities (up to a constant)
           log.P.IS0 <- 0 # these are on log-scale 
         },    
         'monotone.w'={
           if(!('importance.sampling.decreasing' %in% names(prms)))
             prms$importance.sampling.decreasing = FALSE  # default is increasing! 
           w.order <- order(  rowVars(log(w.mat)), decreasing=TRUE  )  # order by variance: first set as high the ones with high variance!! 
  #         w.order <- order(rowSums(w.mat), decreasing=prms$importance.sampling.decreasing)  # order x values based on sum of w. We wnat Small with Large! (at least for w(x,y)=x+y)
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
         'monotone.grid.w'={  # order vars in a decreasing manner. 
           if(!('importance.sampling.decreasing' %in% names(prms)))  
             prms$importance.sampling.decreasing = FALSE  # default is increasing! match small w with large w (can also take monotone by variance)
           w.order <- order(  rowVars(log(w.mat)), decreasing=TRUE  )  # order by variance: first set as high the ones with high variance!! 
##           w.order <- order(rowSums(w.mat), decreasing=prms$importance.sampling.decreasing)  # order x values based on sum of w. We wnat Small with Large! (at least for w(x,y)=x+y)
           log.P.IS <- zeros(prms$B, 1)
           
           num.grid.points <- prms$B / 10 # take 10 permutations at each grid point 
           grid.points <- (0:num.grid.points) /  num.grid.points
           B.per.point <- 10
           log.P.IS0.vec <- zeros(num.grid.points, 1)
           log.P.W.IS.minus.P.W.IS0 <- zeros(prms$B, 1)
           
           for(g in 1:num.grid.points)
           {
             cur.w.mat <- w.mat^grid.points[g]  # take power 
             for(j in 1:n) # Caculate prob. of identity under importance sampling 
             {
               weights <- cur.w.mat[w.order[j],]
               if(j>1)
                 weights[w.order[1:(j-1)]] <- 0    
               #               weights <- weights / sum(weights)
               log.P.IS0.vec[g] <- log.P.IS0.vec[g] + log(weights[w.order[j]] / sum(weights))
             }
             
             for(b in ((g-1)*B.per.point+1):(g*B.per.point))  # sample permutations in this grid point 
             {
               for(j in 1:n) # sample permutation according to w
               {
                 weights <- cur.w.mat[w.order[j],]
                 if(j>1) # remove the ones we already occupied
                   weights[Permutations[w.order[1:(j-1)]  , b]] <- 0                   
                 
                 weights <- weights / sum(weights)
                 Permutations[w.order[j],b] = sample(n, 1, prob = weights)
                 log.P.IS[b] <- log.P.IS[b] + log(weights[Permutations[w.order[j],b]]) # Need to update also P_IS here
               }
               
               log.P.W.IS.minus.P.W.IS0[b] <- log.P.IS[b] + sum(log.w.mat[cbind(1:n, Permutations[,b])]) - 
                 log.P.IS0.vec[g] - sum(diag(log.w.mat))
             }
           }  # end loop on grid points 
           log.P.IS0 <- mean(log.P.IS0.vec) #  / num.grid.points # take average 
         }, 
         # new: monotone grid : interpolate between monotone and uniform 
         
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
         }, 
         "sqrt.w"={
           log.P.IS <- zeros(prms$B, 1)
           log.P.IS0 <- 0
           w.col.sums <- colSums(w.mat)
           for(j in 1:n)
           {
             weights <- w.mat[j,]
             if(j>1)
               weights[1:(j-1)] <- 0  
             
             weights[j:n] <- weights[j:n] / sqrt(sum(weights[j:n]))
             weights[j:n] <- weights[j:n] / sqrt( w.col.sums[j:n]  )
             weights <- weights / sum(weights)
             
             
             log.P.IS0 <- log.P.IS0 + log(weights[j]) # Need to update also P_IS here                       
             w.col.sums <- pmax(0, w.col.sums - w.mat[j,])
             w.col.sums[j] <- 0
           }
           
           for(b in 1:prms$B)  # sample permutations 
           {
             w.col.sums <- colSums(w.mat)
             for(j in 1:n)
             {
               weights <- w.mat[j,]
               if(j>1)
                 weights[Permutations[1:(j-1),b]] <- 0  
               
               weights[Permutations[j:n,b]] <- weights[Permutations[j:n,b]] / sqrt(sum(weights[Permutations[j:n,b]]))
               weights[Permutations[j:n,b]] <- weights[Permutations[j:n,b]] / sqrt( w.col.sums[Permutations[j:n,b]]  )
               weights <- weights / sum(weights)
               
               Permutations[j,b] = sample(n, 1, prob = weights)
               log.P.IS[b] <- log.P.IS[b] + log(weights[Permutations[j,b]]) # Need to update also P_IS here                       
               
               w.col.sums <- pmax(0, w.col.sums - w.mat[j,])
               w.col.sums[Permutations[j,b]] <- 0
               
             }
           }
         },
         "KouMcculough.w"={  # no ordering of variables 
           log.P.IS <- zeros(prms$B, 1)
           log.P.IS0 <- 0
           w.col.sums <- colSums(w.mat)
           for(j in 1:n)
           {
             weights <- w.mat[j,]
             if(j>1)
               weights[1:(j-1)] <- 0  
             
             weights[j:n] <- weights[j:n] / max(epsilon, w.col.sums[j:n] -weights[j:n]) 
             weights <- weights / sum(weights)
             
             
             log.P.IS0 <- log.P.IS0 + log(weights[j]) # Need to update also P_IS here                       
             w.col.sums <- pmax(0, w.col.sums - w.mat[j,])
             w.col.sums[j] <- 0
           }

           b = 1           
           while(b <= prms$B)  # sample permutations . Sometimes can fail for truncation 
           {
             w.col.sums <- colSums(w.mat)
             for(j in 1:n)
             {

               weights <- w.mat[j,]
               if(j>1)
                 weights[Permutations[1:(j-1),b]] <- 0  

               if(max(weights)==0)  # here we failed to sample (can happen for truncation) 
                 break

                              # How do we know the permutations? they are not set yet!
               weights <- weights / max(epsilon, w.col.sums - weights)

#               weights[Permutations[j:n,b]] <- weights[Permutations[j:n,b]] / max(epsilon, w.col.sums[Permutations[j:n,b]] - weights[Permutations[j:n,b]] )
               weights <- weights / sum(weights)

               Permutations[j,b] = sample(n, 1, prob = weights) # need to fix for ties 

               if(weights[Permutations[j,b]]==0)
                 print("Error! sampled zero weight!!")
               log.P.IS[b] <- log.P.IS[b] + log(weights[Permutations[j,b]]) # Need to update also P_IS here                       
               
               w.col.sums <- pmax(0, w.col.sums - w.mat[j,])
               w.col.sums[Permutations[j,b]] <- 0
               
             }
             if((j == n) && (Permutations[j,n]>0))
               b = b+1
           }
         },
         "tsai"={  # new: sample EXACTLY from the truncation distribution using Tsai's algorithm . We need the ordering of x,y !! x[i] <= y[i] (i.e. data[,2] >= data[,1])
           w.order <- order(data[,2])
           R.mat <- matrix(0, n,n)
            for(i in 1:n)  # compute the set R_i for each i
            {
              for(j in 1:n)
                if((data[j,1] <= data[i,2]) && (data[i,2] <= data[j,2]))
                  R.mat[i,j] <- 1
            }
            log.P.IS0 <- sum(log(rowSums(R.mat))) # these are on log-scale. Uniform probability over all legal permutations
            log.P.IS <- rep(log.P.IS0, prms$B)  # get importance probabilities (up to a constant)
            
            R.mat.ordered <- R.mat[w.order, w.order]
            print(R.mat.ordered[1:10,1:10])
            
            for(b in c(1:prms$B))  # go over permutations
            {
              Permutations[,b] <- 1:n # need to initialize to deal with ties 
              Q <- w.order
              for(i in 1:n) # next, sample sequentially using R
              {
                j <- sample(n, 1, prob = R.mat.ordered[i,]/sum(R.mat.ordered[i,])) 
                Permutations[Q[j],b] <- w.order[i]                      # set pre-image 
                # swap
                temp <- Q[i]
                Q[i] <- Q[j]
                Q[j] <- temp
              }
              if( min(data[Permutations[,b],2]-data[,1]) <0 )  # check Tsai's sampling
                print("Erorr! permutation doesn't satisfy truncation!!!!")
            }
         }  # end Tsai 
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
  P[cbind(1:n, 1:n)] <- P[cbind(1:n, 1:n)] + min(1, prms$include.ID)*P.W.IS0 # new: add contribution from identity
  P <- P / (min(1, prms$include.ID)*P.W.IS0 + sum(P.W.IS)) # normalize P    
  
  #    print(P[1:5,1:5])
  
  #    print(t(rowSums(P)))
  #    print(t(colSums(P)))
  if(!exists("log.P.W.IS.minus.P.W.IS0"))
    log.P.W.IS.minus.P.W.IS0 = NA
  
  return(list(Permutations=Permutations, P=P, P.W.IS=P.W.IS, P.W.IS0=P.W.IS0, log.P.W=log.P.W, log.P.W0=log.P.W0, max.w=max.w, max.log=max.log,
              log.P.IS=log.P.IS, log.P.IS0=log.P.IS0, log.P.W.IS.minus.P.W.IS0=log.P.W.IS.minus.P.W.IS0)) # New: return also P, a matrix with Pr(pi(i)=j)
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


