##################################################################################################
# Compute the modified Hoeffding's test statistic, for the permutation test
# Parameters: 
# data - n*2 matrix with (x,y) sample
# grid.points - all possible (x_i,y_j) points  
# null.distribution - n*n matrix with Fx*Fy * W estimate OR 4*n mass table with pre-computed mass estimates 
# 
#  Quardants convension:
#   4 | 1
#   ------
#   3 | 2
##################################################################################################
ComputeStatistic<- function(data, grid.points, null.distribution)
{
  epsilon = 0.00000001
  num.samples=dim(data)[1]
  Obs<-matrix(0,4,1) # observed
  Exp<-matrix(0,4,1) # expected
  perm.flag <- !(min(dim(null.distribution))==num.samples) # check type of null. Doesn't work for n=4
  if(!perm.flag)
  { # New! compute CDF for botstrap. Current implementation ignores ties !
    null.distribution.CDF <- PDFToCDF2d(null.distribution, data) # Could be slow !! why depends on data? 
  }
  
  Statistic <- 0 # = rep(0, dim(grid.points)[1])
  for (i in 1:dim(grid.points)[1] )  # Slow loop on grid points 
  {
    if(perm.flag) # here we already computed the null in a 4*n table 
    {
      Exp = pmax(null.distribution[i,], epsilon)
    }
    else {        #UP right (Expected)
      for(j in 1:3)
      {
        #        Exp[j] = GetQuarterExpectedProb(grid.points[i,], j, data, null.distribution) # old version uses PDF 
        Exp[j] <- GetQuarterExpectedProb(grid.points[i,], j, data, null.distribution.CDF)
      }
      Exp[4] = max(1-sum(Exp[1:3]), epsilon)
      Exp <- num.samples*Exp
    }
    #Up right quarter (Observed) - this is the computationally costly part for permutations test 
    Rx <- data[,1]>=grid.points[i,1]
    Ry <- data[,2]>=grid.points[i,2]
    Obs[1] <- sum(Rx*Ry)
    Obs[2] <- sum(Rx)-Obs[1]
    Obs[4] <- sum(Ry)-Obs[1]
    Obs[3] <- num.samples-sum(Obs[c(1,2,4)])
    
#    if(any(is.na(Exp)))
#      print(paste("i=", i, " Exp=", min(Exp), ' Obs=', min(Obs)))
    Statistic <-  Statistic + sum((Obs-Exp)^2 / pmax(Exp, 0.0001)) # set valid statistic when expected is 0 or very small 
  } # end loop on grid points 
  
  return(Statistic) # Excluded Nans #   sum(Statistic[!(is.nan(Statistic) | is.infinite(Statistic))]))
}

#########################################################################################
# sample permutations, using MCMC, over the set of valid permutations, 
# with respect to the distribution appears in Eq 8
# Parameters: 
# W - matrix with weights 
# B - number of permutations to draw
# N - sample size (can be read from data or W?)
#########################################################################################
PermutationsMCMC<-function(W, B, N)
{ 
  # Set mcmc sampling parameters 
  burn.in = 2*N
  Cycle = N
  ctr = 1
  Idx = 1
  PermutationsTable = matrix(0,N,B)
  Perm = 1:N
  while(Idx<=B)
  {
    #A Metropolis Hastings algorithm with target stationary distribution \pi
    #Choose the two indices to be switched
    switchIdx = sample(1:N, 2, replace = FALSE)  
    i = switchIdx[1];
    j = switchIdx[2];
    ratio = W[i,Perm[j]]*W[j,Perm[i]]/(W[i,Perm[i]]*W[j,Perm[j]]) # New! Big-fix (?) denomerator didn't have parenthesis
    
    #    print(paste0("ratio=", ratio, ", W=", W[i,Perm[j]], " ", W[j,Perm[i]], " ", W[i,Perm[i]], " ", W[j,Perm[j]]))
    if(rbinom(1, 1, min(1,ratio))) #we accept the transition with probability min(1, ratio)
    {
      temp <- Perm[i] # SWAP 
      Perm[i] <- Perm[j]
      Perm[j] <- temp
      if(ctr==burn.in || (ctr%%Cycle==0 && ctr>burn.in))
      {
        PermutationsTable[,Idx]=Perm;
        Idx = Idx+1;
        if(mod(Idx,100)==0)
          print(c("Sample Perm=", Idx))
      }
      ctr <- ctr+1;
    }
  }
  return(PermutationsTable)
}

###################################################################################
# Estimate the null distribution fx*fy*W (given the estimated PDFs f_x, f_y)
# 
# Parameters: 
# pdfs - marginal distributions fx, fy probabilities 
# w - n*n matrix with weights W[i,j] = w(x_i, y_j) 
###################################################################################
GetNullDistribution <- function(pdfs, W)
{
  #compute the normalizing factor under the null:
  null.distribution <- W * (pdfs[,1] %*% t(pdfs[,2]))
  Z <- sum(null.distribution)
  null.distribution<-null.distribution/Z
  return( list(null.distribution=null.distribution, normalizing.factors=Z) )
}

############################################################################################
# New: draw a bootstrap sample, given the estimated null distribution Fx, Fy, W
# Use rejection sampling. No need for computing n*n table 
# (Problem: what if prob.(rejection) close to 1?)
# Parameters: 
# data - n*2 array with (X,Y) samples
# pdfs - fx and fy  
# bias.type - W
# prms - for w max 
############################################################################################
Bootstrap <- function(data, pdfs, bias.type, prms)
{
  n = dim(data)[1]
  boot_sample<-matrix(0,n,2)
  k <- 0
  xy <- c(0,0)  
  while(k<n) # could be slow - try to not sample 1-by-1
  {
    # New: faster sampling n-k together
    x <- data[sample(n, n-k, prob=pdfs[,1], replace=TRUE),1] # Sample X from Fx
    y <- data[sample(n, n-k, prob=pdfs[,2], replace=TRUE),2] # Sample Y sample from Fy
    keep <- which(as.logical(rbinom(n-k, 1, BiasedSamplingW(x, y, bias.type)/prms$W_max)))
    boot_sample[keep+k,] <- cbind(x[keep],y[keep]) 
    k <- k+length(keep)
    
    # Old: one-by-one    
    #     xy[1] <- data[sample(n, 1, prob=pdfs[,1]),1] # Sample X from Fx
    # xy[2] <- data[sample(n, 1, prob=pdfs[,2]),2] # Sample Y sample from Fy
    # keep <- rbinom(1, 1, BiasedSamplingW(xy[1], xy[2], bias.type)/prms$W_max)
    # if(keep) 
    # {
    #   boot_sample[k,] <- xy
    #   k <- k+1
    # }
  }    
  return(boot_sample)
}


###################################################################################################
# given a quadrant, evaluate the mass function within it
# New version - uses the two-dimensional CDF (not PDF)
# 
# Point - (x,y) value defining quadrants 
# QId - which quadrant 1,2,3,4
# data - 2*n array of (x,y)
# null.distribution.CDF - n*n table of 2d cumulative of Fx*Fy*w (problem with ties)
###################################################################################################
GetQuarterExpectedProb <- function(Point, QId, data, null.distribution.CDF)
{
  if(QId %in% c(1,2))
  {
    idx.x <- which(data[,1]>=Point[1])
    idx.x <- idx.x[which.min(data[idx.x,1])]
  } else
  {
    idx.x <- which(data[,1]<Point[1])
    idx.x <- idx.x[which.max(data[idx.x,1])]
  }
  if(QId %in% c(1,4))
  {
    idx.y <- which(data[,2]>=Point[2])
    idx.y <- idx.y[which.min(data[idx.y,2])]
  } else
  {
    idx.y <- which(data[,2]<Point[2])
    idx.y <- idx.y[which.max(data[idx.y,2])]
  }
  
  if(isempty(idx.x) | isempty(idx.y))
    return(0)    
  
  m <- which.max(data[,1])
  n <- which.max(data[,2])
  
  if(QId == 1)
    S <- null.distribution.CDF[m, n] + null.distribution.CDF[idx.x, idx.y] - 
    null.distribution.CDF[idx.x, n]  - null.distribution.CDF[m, idx.y]
  if(QId == 2)
    S <- null.distribution.CDF[m, idx.y]  - null.distribution.CDF[idx.x, idx.y]  
  if(QId == 3)    
    S <- null.distribution.CDF[idx.x, idx.y]  
  if(QId == 4)
    S <- null.distribution.CDF[idx.x, n]  - null.distribution.CDF[idx.x, idx.y]  
  
  return(S)       
}



#computing Expect[Q_i(p_j)] for 1<=i<=4, and all j, given a grid of points
# (New, faster implementation)
# Parameters: 
# data - 2*n array (X,Y)
# Permutations - set of permutations
# grid.points - centers of partitions
#
# Output: 
# mass.table - a 4*#grid-points table with quadrants expectation
QuarterProbFromPermutations <- function(data, Permutations, grid.points)
{
  num.permutations <- dim(Permutations)[2]
  n <- dim(data)[1]
  P <- matrix(0, n, n) # matrix with P[i]=j estimate
  for(i in 1:num.permutations)
    P[cbind(1:n, Permutations[,i])] <-  P[cbind(1:n, Permutations[,i])]+1 # need to vector indices here  
  P <- P / num.permutations # normalize 
  
  #next each point has its probability of being selected we evaluate Expect[Q_i(p_j)] by summation
  mass.table<-matrix(0, dim(grid.points)[1], 4)
  for(i in seq(1, dim(grid.points)[1],1))
  {
    x<-grid.points[i,]
    mass.table[i,1]<-sum(P[data[,1]>=x[1], data[,2]>=x[2]])
    mass.table[i,2]<-sum(P[data[,1]>=x[1], data[,2]<x[2]])
    mass.table[i,3]<-sum(P[data[,1]<x[1], data[,2]<x[2]])
    mass.table[i,4]<-sum(P[data[,1]<x[1], data[,2]>=x[2]])
  } 
  mass.table<-mass.table+0.000001
  return(mass.table)
}




# New: cumulative 2d (for faster calculation of test statistic)
# When we have ties we need to correct 
PDFToCDF2d <- function(pdf.2d, data)
{
  Px <- sort(data[,1], index.return=TRUE)  # Permute to order x_i, y_i 
  Py <- sort(data[,2], index.return=TRUE)
  
  cdf.2d <- apply(pdf.2d[Px$ix, Py$ix], 2, cumsum)
  cdf.2d <- apply(cdf.2d, 1, cumsum)
  
  # Deal with ties (not working yet)
  #  ties.x <- which(data[-1,1] == head(data[,1], -1))
  #  ties.y <- which(data[-1,2] == head(data[,2], -1))
  #  print(paste0('ties.x=', ties.x))
  #  print(paste0('ties.y=', ties.y))
  #  for(i in rev(ties.y))
  #    cdf.2d[,i] <- cdf.2d[,i+1]
  #  for(i in rev(ties.x))
  #    cdf.2d[i,] <- cdf.2d[i+1,]
  
  return( t(cdf.2d[invPerm(Py$ix), invPerm(Px$ix)]))  
  #    return( t(apply(apply(pdf.2d, 2, cumsum), 1, cumsum)) )
}


# The true marginals (same for x and y) for uniform strip of width a
UniformStripMarginal <- function(t,a){
  if (a>0.5 | a<0) stop("a must be in (0,0.5]")
  c <- a*(2-a)
  if (t<a){
    val <- t*(t+2*a)/(2*c)
  } else if (t<1-a) {
    val <- 3*a^2/(2*c) + 2*a*(t-a)/c
  } else {
    val <- 3*a^2/(2*c) + 2*a*(1-2*a)/c + (1+3*a-t)*(t+a-1)/(2*c)
  }
  return(val)
}




