library(Matrix)

##################################################################################################
# Compute the modified Hoeffding's test statistic, for the permutation test
# Parameters: 
# data - n*2 matrix with (x,y) sample
# grid.points - n*2 all possible (x_i,y_j) points  (can be the same as data, or a subset)
# null.expectations.table - n*4 mass table with pre-computed mass estimates for grid.points 
# 
#  Quardants convension:
#   4 | 1
#   ------
#   3 | 2
##################################################################################################
ComputeStatistic <- function(data, grid.points, null.expectations.table)
{
  # treat grid.points 
  if(missing(grid.points) || isempty(grid.points))
    grid.points <- data  # default: set grid as data 
  
  obs.table <- matrix(0, dim(grid.points)[1], 4)
  Obs <- Exp <- matrix(0,4,1) # observed & expected
  
  #  print("DIM NULL TABLE")
  #  print(dim(null.expectations.table))
  Statistic <- 0 
  for (i in 1:dim(grid.points)[1])  # Slow loop on grid points 
  {
    Exp <- null.expectations.table[i,]
    Rx <- data[,1] > grid.points[i,1]
    Ry <- data[,2] > grid.points[i,2]
    # new: deal also with ties 
    Eqx <- data[,1] == grid.points[i,1]
    Eqy <- data[,2] == grid.points[i,2]
#    if((Eqx > 0) || (Eqy > 0))
#      print(paste0("Equal X,Y, Xy:", sum(Eqx), " ", sum(Eqy), " ", sum(Eqx*Eqy)))
    Obs[1] <- sum(Rx*Ry)
    Obs[2] <- sum(Rx)-Obs[1]-sum(Eqx)
    Obs[4] <- sum(Ry)-Obs[1]-sum(Eqy)
    Obs[3] <- dim(data)[1]-sum(Obs[c(1,2,4)]) - sum(Eqx) - sum(Eqy) + sum(Eqx*Eqy) # new! don't count points equal 
    obs.table[i,] <- Obs
    
#    if(i < 5)
#      print(paste0("i:", i, " Sum-Exp:", sum(Exp), " Sum-Obs:", sum(Obs)))
    if (min(Exp)>1) {
      #      print("Add To Expected")
      Statistic <-  Statistic + sum((Obs-Exp)^2 / Exp) # set valid statistic when expected is 0 or very small 
      #      print(Statistic)
    }
  } # end loop on grid points 
  
  return(list(Statistic=Statistic, obs.table=obs.table)) # return also observed table for diagnostics
}


##################################################################################################
# Compute the modified Hoeffding's test statistic corrected for w:
# Observed and Expected calculated using inverse weighting
# GOOD ONLY FOR A STRICTLY POSITIVE w
# Parameters: 
# data - n*2 matrix with (x,y) sample
# grid.points - all possible (x_i,y_j) points  
# w - weight function (default w(x,y)=1) 
# 
#  Quardants convension:
#   4 | 1
#   ------
#   3 | 2
##################################################################################################
ComputeStatistic.W <- function(data, grid.points, w.fun=function(x){1})
{
  # treat grid.points 
  if(missing(grid.points) | isempty(grid.points))
    grid.points <- unique.matrix(data)  # default: set unique for ties? for discrete data
  
  w.vec <- w_fun_eval(data[,1], data[,2], w.fun) # w_fun_to_mat(data, w.fun) # calculate w n*n matrix 
  n.w <- sum(1/w.vec)
  obs.table <- exp.table <- matrix(0, dim(grid.points)[1], 4)
  Obs <- Exp <- matrix(0,4,1) # observed & expected
  Statistic <- 0 
  for (i in 1:dim(grid.points)[1])  # Slow loop on grid points 
  {
    Rx <- data[,1] > grid.points[i,1] # ties are ignored in both observed and expected 
    Ry <- data[,2] > grid.points[i,2]
    Lx <- data[,1] < grid.points[i,1] 
    Ly <- data[,2] < grid.points[i,2]
    Exp[1] <- sum(Rx/w.vec)*sum(Ry/w.vec)/n.w^2
    Exp[2] <- sum(Rx/w.vec)*sum(Ly/w.vec)/n.w^2
    Exp[4] <- sum(Lx/w.vec)*sum(Ry/w.vec)/n.w^2
    Exp[3] <- sum(Lx/w.vec)*sum(Ly/w.vec)/n.w^2
    Obs[1] <- sum(Rx*Ry/w.vec)/n.w
    Obs[2] <- sum(Rx*Ly/w.vec)/n.w
    Obs[4] <- sum(Lx*Ry/w.vec)/n.w
    Obs[3] <- sum(Lx*Ly/w.vec)/n.w
    obs.table[i,] <- Obs
    exp.table[i,] <- Exp
    if (min(Exp)>(1/dim(data)[1])) {
      Statistic <-  Statistic + sum((Obs-Exp)^2 / Exp) # set valid statistic when expected is 0 or very small 
    } 
  } # end loop on grid points 
  
  return(list(Statistic=Statistic, obs.table=obs.table, exp.table=exp.table))  
  # returns also expected and observed tables for diagnostics
}

##################################################################################################
#  Same as ComputeStatistic.W, but with test statistic multiplied by n.w
#  (this is the standard form of (exp-obs)^2/exp, ComputeStatistic.W uses (exp/n.w-obs/n.w)^2/(exp/n.w)
#   
##################################################################################################
ComputeStatistic.W1 <- function(data, grid.points, w.fun=function(x){1})
{
  # treat grid.points 
  if(missing(grid.points) | isempty(grid.points))
    grid.points <- unique.matrix(data)  # default: set unique for ties? for discrete data
  
  w.vec <- w_fun_eval(data[,1], data[,2], w.fun) # w_fun_to_mat(data, w.fun) # calculate w n*n matrix 
  n.w <- sum(1/w.vec)
  obs.table <- exp.table <- matrix(0, dim(grid.points)[1], 4)
  Obs <- Exp <- matrix(0,4,1) # observed & expected
  Statistic <- 0 
  for (i in 1:dim(grid.points)[1])  # Slow loop on grid points 
  {
    Rx <- data[,1] > grid.points[i,1] # ties are ignored in both observed and expected 
    Ry <- data[,2] > grid.points[i,2]
    Lx <- data[,1] < grid.points[i,1] 
    Ly <- data[,2] < grid.points[i,2]
    Exp[1] <- sum(Rx/w.vec)*sum(Ry/w.vec)/n.w
    Exp[2] <- sum(Rx/w.vec)*sum(Ly/w.vec)/n.w
    Exp[4] <- sum(Lx/w.vec)*sum(Ry/w.vec)/n.w
    Exp[3] <- sum(Lx/w.vec)*sum(Ly/w.vec)/n.w
    Obs[1] <- sum(Rx*Ry/w.vec)
    Obs[2] <- sum(Rx*Ly/w.vec)
    Obs[4] <- sum(Lx*Ry/w.vec)
    Obs[3] <- sum(Lx*Ly/w.vec)
    obs.table[i,] <- Obs
    exp.table[i,] <- Exp
    if (min(Exp)>(1/dim(data)[1])) {
      Statistic <-  Statistic + sum((Obs-Exp)^2 / Exp) # set valid statistic when expected is 0 or very small 
    } 
  } # end loop on grid points 
  
  return(list(Statistic=Statistic, obs.table=obs.table, exp.table=exp.table))  
  # returns also expected and observed tables for diagnostics
}



#################################################################
# sample permutations, using MCMC, over the set of valid permutations, 
# with respect to the distribution appears in Eq 8
# Parameters: 
# w.mat - matrix with weights 
# prms - parameters, including B, number of permutations to draw
# 
# Output: 
# Permutations - A matrix representing the sampled permutations 
# P - An n*n matrix with P(i,j) = Pr(pi(i)=j) over the sampled permutations 
#########################################################################################
PermutationsMCMC <- function(w.mat, prms) # burn.in=NA, Cycle=NA)  # New: allow non-default burn-in 
{ 
  n <- dim(w.mat)[1];
  P <- matrix(0, n, n) # New! matrix with P[i]=j estimate
  #  for(i in 1:num.permutations) 
  #    P[cbind(1:n, Permutations[,i])] <- P[cbind(1:n, Permutations[,i])]+1 # need to vector indices here  
  
  # Set mcmc default sampling parameters 
  if(!('B' %in% names(prms)))
    prms$B <- 1000
  if(!('burn.in' %in% names(prms)))
    prms$burn.in <- 0  # 2*n  # set to n without burn.in, or 0 without burn.in to include the identity permutation 
  if(!('Cycle' %in% names(prms)))
    prms$Cycle <- n

#  print(paste0("burn.in=", prms$burn.in, " cycle=", prms$Cycle))    
  Idx <- ctr <- 1
  Permutations = matrix(0, n, prms$B)
  Perm = 1:n # start with the identity 
  while(Idx <= prms$B)
  {
    # A Metropolis Hastings algorithm with target stationary distribution \pi
    # Choose the two indices to be switched
    switchIdx = sample(1:n, 2, replace = FALSE)  
    i = switchIdx[1]
    j = switchIdx[2]
    ratio = w.mat[i,Perm[j]]*w.mat[j,Perm[i]] / (w.mat[i,Perm[i]]*w.mat[j,Perm[j]]) 
    
    if(rand() < ratio) #     rbinom(1, 1, min(1,ratio))) #we accept the transition with probability min(1, ratio)
    {
      temp <- Perm[i] # SWAP 
      Perm[i] <- Perm[j]
      Perm[j] <- temp
      if(ctr==prms$burn.in || (ctr%%prms$Cycle==0 && ctr>prms$burn.in)) # save permutation
      {
        Permutations[,Idx]=Perm;
        Idx = Idx+1;
        #          if(mod(Idx,100)==0)
        #            print(c("Sample Perm=", Idx))
      }
      P[cbind(1:n, Perm)] <- P[cbind(1:n, Perm)]+1 # Update table P
      ctr <- ctr+1;  # update counter only after swap 
    }
  }  # end while
  P <- P / (ctr-1) # normalize 
  
  return(list(Permutations=Permutations, P=P)) # New: return also P, a matrix with Pr(pi(i)=j)
}

###################################################################################
# Estimate the null distribution fx*fy*W (given the estimated PDFs f_x, f_y)
# 
# Parameters: 
# pdfs - 2*n table with marginal distributions fx, fy probabilities 
# w - n*n matrix with weights W[i,j] = w(x_i, y_j) 
###################################################################################
GetNullDistribution <- function(pdfs, w_mat)
{
  # Compute the normalizing factor under the null:
  null.distribution <- w_mat * (pdfs[,1] %*% t(pdfs[,2]))
  Z <- sum(null.distribution)
  null.distribution <- null.distribution/Z
  
  #  print("NULL DIM")
  #  print(dim(null.distribution))
  return( list(distribution=null.distribution, Z=Z) )
}

############################################################################################
# Draw a bootstrap sample, given the estimated null distribution Fx, Fy, W
# Use rejection sampling. No need for computing n*n table 
# (Problem: what if prob.(rejection) close to 1?)
# Parameters: 
# data - n*2 array with (X,Y) samples
# pdfs - fx and fy  
# w.fun - Biased sampling function w 
# prms - for w max 
# n - allow a different sample size 
############################################################################################
Bootstrap <- function(data, pdfs, w.fun, prms, n=NULL)
{
  #  print("Inside Bootstrap")
  if(is.null(n))
    n = dim(data)[1]
  #  print(dim(data))
  boot.sample <- matrix(-1, n, 2) # wtf
  indices <- matrix(0, n, 2)
  k <- 0
  ctr <- 0
  while(k<n) 
  {   # sampling n-k together
    I <- sample(dim(pdfs)[1], n-k, prob=pdfs[,1], replace=TRUE)
    J <- sample(dim(pdfs)[1], n-k, prob=pdfs[,2], replace=TRUE)
    x <- data[I,1] # Sample X ~ Fx
    y <- data[J,2] # Sample Y ~ Fy
    if(max(w_fun_eval(x, y, w.fun)/prms$w.max)>1)
    {
      print(paste0("WTF MAX: ", max(w_fun_eval(x, y, w.fun)/prms$w.max)))
      print(cbind(I, J))
    }
      
    keep <- which(as.logical(rbinom(n-k, 1, w_fun_eval(x, y, w.fun)/prms$w.max)))
    if(isempty(keep))
      next
    boot.sample[(1:length(keep))+k,] <- cbind(x[keep],y[keep]) 
    indices[(1:length(keep))+k,] <- cbind(I[keep], J[keep]) # maybe order is a problem? 
    ctr <- ctr + n-k
    k <- k+length(keep)
#    print(x[keep] - data[indices[1:k,1],1])
#    print(y[keep] - data[indices[1:k,2],2])
#    print(boot.sample[1:k,1] - data[indices[1:k,1],1])
#    print(boot.sample[1:k,2] - data[indices[1:k,2],2])
  }    
  #  print(paste0("Sampled in total k=", ctr, " to get n=", n, " samples"))
  return(list(sample=boot.sample, indices=indices))  # new: return also indices! (they're all the information we need!)
}

###################################################################################################
# given a quadrant, evaluate the mass function within it
# New version - uses the two-dimensional CDF (not PDF)
# 
# Point - (x,y) value defining quadrants 
# QId - which quadrant 1,2,3,4
# data - 2*n array of (x,y)
# null.distribution.CDF - n*n table of 2d cumulative of Fx*Fy*w (problem with ties)
# Output: 
# S - expected probability at the quardant QId centered at Point for the distribution defined by null.distribution.CDF and data 
###################################################################################################
GetQuarterExpectedProb_old <- function(Point, QId, data, null.distribution.CDF)
{
  if(QId %in% c(1,2))
    idx.x <- which(data[,1] <= Point[1])  # >= 
  else
    idx.x <- which(data[,1] < Point[1]) # <=
  idx.x <- idx.x[which.max(data[idx.x,1])] # min
#  if(isempty(idx.x))  # new: take care of edges - why take the min? 
#    idx.x <- which.min(data[,1])

  if(QId %in% c(1,4))
    idx.y <- which(data[,2] <= Point[2]) # >= 
  else
    idx.y <- which(data[,2] < Point[2])
  idx.y <- idx.y[which.max(data[idx.y,2])] # min 
#  if(isempty(idx.y))  # new: take care of edges 
#    idx.y <- which.min(data[,2])

#  print("IDX X,Y:")
#  print(idx.x)
#  print(idx.y)
  if(isempty(idx.x) | isempty(idx.y)) # didn't find any - why zero? 
    cdf.point <- 0  #    return(0)     
  else
    cdf.point <- null.distribution.CDF[idx.x, idx.y]
#   print("CDFPONIT")
#   print(cdf.point)

  m <- which.max(data[,1])
  n <- which.max(data[,2])
  if(isempty(idx.x))
    cdf.x <- 0
  else
    cdf.x <- null.distribution.CDF[idx.x, n]
  if(isempty(idx.y))
    cdf.y <- 0
  else
    cdf.y <- null.distribution.CDF[m, idx.y]
        
#  print(paste0("CDF2D: x,y: ", null.distribution.CDF[idx.x, idx.y], 
#               " x, n: ", null.distribution.CDF[idx.x, n], 
#               " m, y:", null.distribution.CDF[m, idx.y], 
#               " m, n:", null.distribution.CDF[m, n]))
  
  switch(QId, # First sample from Fxy
         {S <- 1 + cdf.point - cdf.x - cdf.y}, # 1
         {S <- cdf.y - cdf.point}, # 2
         {S <- cdf.point}, # 3
         {S <- cdf.x - cdf.point}) # 4
  return(S)       
}

# New function version: using ecdf
GetQuarterExpectedProb <- function(Point, QId, data, null.distribution.CDF)
{
  Point.minus <-  Point - .Machine$double.eps*100 # need to have lower tolerance! 
  switch(QId, # First sample from Fxy
         {s <- 1 + ecdf2(Point, null.distribution.CDF, data)  - 
           ecdf2(c(Point[1], max(data[,2])), null.distribution.CDF, data) - 
           ecdf2(c(max(data[,1]), Point[2]), null.distribution.CDF, data)}, # 1
         {s <- ecdf2(c(max(data[,1]), Point.minus[2]), null.distribution.CDF, data) - 
           ecdf2(c(Point[1], Point.minus[2]), null.distribution.CDF, data)}, # 2
         {s <- ecdf2(Point.minus, null.distribution.CDF, data)}, # 3
         {s <-  ecdf2(c(Point.minus[1], max(data[,2])), null.distribution.CDF, data) - 
           ecdf2(c(Point.minus[1], Point[2]), null.distribution.CDF, data)}) # 4
  return(s)       
}  

###################################################################################################
# Compute Expect[Qi(p_j)] for 1<=i<=4, and all j, given a grid of points and bootstrap null distribution
# Parameters: 
# data - 2*n array (X,Y)
# null.distribution - a two-dim array n*n of probabilities under the null 
# Permutations - set of permutations
# grid.points - centers of partitions
#
# Output: 
# mass.table - a 4*#grid-points table with quadrants expectation
###################################################################################################
QuarterProbFromBootstrap <- function(data, null.distribution, grid.points)
{
  mass.table <- matrix(0, dim(grid.points)[1], 4)
  
  #  print("NULL-DIST-DIM:")
  #  print(dim(null.distribution))
  null.distribution.CDF <- PDFToCDF2d(null.distribution, data) 
  
  for(i in seq(1, dim(grid.points)[1],1)) # find empty indices 
  {
    for(j in 1:4) # print index ? 
    {
      mass.table[i,j] <- GetQuarterExpectedProb(grid.points[i,], j, data, null.distribution.CDF)  # new! try using ecdf2 
#      temp.cpp <- GetQuarterExpectedProb_rcpp(grid.points[i,], j-1, data, null.distribution.CDF)
#      temp.cpp2 <- GetQuarterExpectedProb_rcpp2(grid.points[i,], j-1, data, null.distribution.CDF)
#      print(paste0("j=", j, " R=", mass.table[i,j], " cpp=", temp.cpp, " cpp2=", temp.cpp2))
#      x.diff = mass.table[i,j] - temp.cpp + 0.0000000000000001
      
#      if(j == 3)  # just take ecdf2
#      {
#        print(paste0("ecdf2: ", ecdf2(grid.points[i,], null.distribution.CDF, data)))
#        print(paste0("ecdf2.cpp: ", ecdf2_rcpp(grid.points[i,], null.distribution.CDF, data)))
#      }
        
      
#      tmp_debug <-   GetQuarterExpectedProb(grid.points[i,], j, data, null.distribution.CDF)
#      if(abs(tmp_debug - mass.table[i,j]) > 0.000001)
#      {
#        print("Quarter Mismatch!!!")
#        save(grid.points, data, null.distribution.CDF, file='BadQuarter.Rdata')
#        print(paste0("i,j=", i, ", ", j))        
#      }
    }      
      
#    mass.table[i,4] = 1-sum(mass.table[i,1:3]) # , epsilon)
  }
#  mass.table <- dim(data)[1]*mass.table # NEW! Do NOT normalize to counts - this should be sample size n, NOT the data size! 
  
  return(mass.table)
}

###################################################################################################
# compute Expect[Qi(p_j)] for 1<=i<=4, and all j, given a grid of points using permutations
# Parameters: 
# data - 2*n array (X,Y)
# P - n*n table representing permutations probabilities Pr(pi(i)=j) (# Permutations - set of permutations)
# grid.points - centers of partitions 2*k array 
#
# Output: 
# mass.table - a 4*#grid-points table with quadrants expectation
###################################################################################################
QuarterProbFromPermutations <- function(data, P, grid.points) #Permutations
{
  #  num.permutations <- dim(Permutations)[2]
  #  n <- dim(data)[1]
  #  P <- matrix(0, n, n) # matrix with P[i]=j estimate
  #  for(i in 1:num.permutations) 
  #    P[cbind(1:n, Permutations[,i])] <- P[cbind(1:n, Permutations[,i])]+1 # need to vector indices here  
  #  P <- P / num.permutations # normalize 
  
  #next each point has its probability of being selected we evaluate Expect[Qi(p_j)] by summation
  mass.table<-matrix(0, dim(grid.points)[1], 4)
  for(i in seq(1, dim(grid.points)[1],1))
  {
    x<-grid.points[i,]
    mass.table[i,1] <- sum(P[data[,1]>x[1], data[,2]>x[2]])
    mass.table[i,2] <- sum(P[data[,1]>x[1], data[,2]<=x[2]])
    mass.table[i,3] <- sum(P[data[,1]<=x[1], data[,2]<=x[2]])
    mass.table[i,4] <- sum(P[data[,1]<=x[1], data[,2]>x[2]])
  } 
  #  mass.table<-mass.table+0.000001
  return(mass.table)
}

###################################################################################################
# Compute quarter probabilities for the standard bivariate Gaussians 
# with truncation y>x. We assume grid.points satisfy y>x
#
# grid.points - where to evaluate probabilitites 
###################################################################################################
QuarterProbGaussianAnalytic <- function(grid.points, w.fun)
{
  mass.table<-matrix(0, dim(grid.points)[1], 4)
  x <- grid.points[,1]
  y <- grid.points[,2]

  if(w.fun %in% c("truncation", "Hyperplane_truncation")) 
  {
    mass.table[,1] <- (1-pnorm(y))*(1+pnorm(y)-2*pnorm(x))
    mass.table[,2] <- (pnorm(x)-pnorm(y))^2
    mass.table[,3] <- pnorm(x)*(2*pnorm(y)-pnorm(x))
    mass.table[,4] <- 2*pnorm(x)*(1-pnorm(y))
  }
  if(w.fun %in% c('naive', 'const'))
  {
    mass.table[,1] <- (1-pnorm(x))*(1-pnorm(y))
    mass.table[,2] <- (1-pnorm(x))*pnorm(y)
    mass.table[,3] <- pnorm(x)*pnorm(y)
    mass.table[,4] <- pnorm(x)*(1-pnorm(y))
  }
  return(mass.table)
}


##############################################################################
# Convert marginal PDF to PDF from data 
# Parameters: 
# CDF.table - matrix with every column a vector of CDF of each variable 
# Output: 
# PDF.table - matrix with every column a vector of PDF of each variable  
##############################################################################
CDFToPDFMarginals <- function(CDF.table)
{
  n<-dim(CDF.table)[1]  # number of samples 
  PDF.table <- array(0L, dim(CDF.table))  # matrix(0, num.samples, num.variables)
  for(i in 1:dim(CDF.table)[2])  # loop on variables 
  {
    Px <- sort(CDF.table[,i], index.return=TRUE)
    PDF.table[Px$ix,i] <- c(Px$x[1], Px$x[-1]-Px$x[-n])
  }
  return(PDF.table)
}



##############################################################################
# Convert marginal PDF to CDF from data 
# Parameters: 
# data - k*n array of variables 
# PDF.table - matrix with every column a vector of PDF of each variable  
# Output: 
# CDF.table - matrix with every column a vector of CDF of each variable 
##############################################################################
PDFToCDFMarginals <- function(data, PDF.table)
{
  n<-dim(PDF.table)[1]  # number of samples 
  CDF.table <- array(0L, dim(PDF.table))  # matrix(0, num.samples, num.variables)
  for(i in 1:dim(PDF.table)[2])  # loop on variables 
  {
    Px <- sort(data[,i], index.return=TRUE)  # Permute to order x_i, y_i 
    CDF.table[Px$ix,i] <- cumsum(PDF.table[Px$ix,i])
    for(j in seq(n, 2))
      if(data[Px$ix[j-1],i] == data[Px$ix[j],i])
        CDF.table[Px$ix[j-1],i] = CDF.table[Px$ix[j],i]     # new: set CDF for ties
  }
  return(CDF.table)
}


###################################################################################################
# Compute 2d cumulative distribution. When we have ties we need to correct this
# Input: 
# pdf.2d - a two-dim array n*n of probabilities 
# data - xy points with probabilities, array size n*2 (used for sorting)
# Output: 
# cdf.2d - a two-dim array n*n of cumulative probabilities 
###################################################################################################
PDFToCDF2d <- function(pdf.2d, data)
{
  Px <- sort(data[,1], index.return=TRUE)  # Permute to order x_i, y_i 
  Py <- sort(data[,2], index.return=TRUE)
  cdf.2d <- apply(apply(pdf.2d[Px$ix, Py$ix], 1, cumsum), 1, cumsum)  # cumsum on rows and columns 
  
  n <- dim(pdf.2d)[1]
  for(i in seq(n, 2))
    if(data[Px$ix[i-1],1] == data[Px$ix[i],1])
      cdf.2d[i-1,] = cdf.2d[i,]     # new: set CDF for ties
  for(j in seq(n, 2))
    if(data[Py$ix[j-1],2] == data[Py$ix[j],2])
      cdf.2d[, j-1] = cdf.2d[, j]     # new: set CDF for ties
        
  return( cdf.2d[invPerm(Px$ix), invPerm(Py$ix)] )  # why Py first and then Px? 
}


# Compute one-dimensional empirical cdf
ecdf1 <- function(x, cdf, data)
{
  i <- which(data <= x)
  i <- i[which.max(data[i])]
  if(isempty(i))
    return(0.0)
  return(cdf[i])  
}


# Compute two-dimensional empirical cdf
ecdf2 <- function(xy, cdf.2d, data)
{
  i <- which(data[,1] <= xy[1])
  i <- i[which.max(data[i,1])]
  j <- which(data[,2] <= xy[2])
  j <- j[which.max(data[j,2])]
  if(isempty(i) | isempty(j))
    return(0.0)
  
#  print(paste0("MAX.INDEX i,j: ", i, ", ", j))
  return(cdf.2d[i,j])  
}

###################################################################################################
# The true marginals (same for x and y) for uniform distribution on strip of width a
# t - where to evaluate CDF
# a - width of strip 
###################################################################################################
UniformStripMarginalCDF <- function(t,a)
{
  if (a>0.5 | a<0) stop("a must be in (0,0.5]")
  c <- a*(2-a)
  val <- 3*a^2/(2*c) + 2*a*(1-2*a)/c + (1+3*a-t)*(t+a-1)/(2*c) # all indices 
  val[t<a] <- t[t<a]*(t[t<a]+2*a)/(2*c) # for t<a
  val[(t<1-a) & (t>=a)] <- 3*a^2/(2*c) + 2*a*(t[(t<1-a) & (t>=a)]-a)/c # others with t<1-a
  return(val)
}


###################################################################################################
# Plot scatter and estimated mariginals 
# dependence.type - dependency label
# biased.data - 2*n array with sampled data
# prms - structure with parameters 
###################################################################################################
PlotBiasedData <- function(dependence.type, biased.data, prms)
{
  n <- dim(biased.data)[1]
  output.scatter.file <- paste0(path, '/../Figures/simulations/', dependence.type, '_rho_', prms$rho)  # set output figure file name
  for(plot.type in c('scatter', 'KM_marginal', 'marginal_scatter'))
  {
    switch(plot.type, # First sample from Fxy
           'scatter' ={
             xy <- data.frame(biased.data)   # plot and save to image file 
           },
           'KM_marginal' ={
             marginals <- EstimateMarginals(biased.data, 'survival')
             marginals.CDFs <- data.frame(cbind(marginals$xy, marginals$CDFs, 
                                                seq(0, 1, len=n), UniformStripMarginalCDF(seq(0, 1, len=n), prms$rho)))   # estimated marginals 
           },
           'marginal_scatter' ={
             Px <- c(marginals$CDFs[-1,1], 1)-marginals$CDFs[,1]
             Py <- marginals$CDFs[,2]-c(0, marginals$CDFs[-n,2])
             x.ind <- sample(x = marginals$xy[,1], size = 10000, replace = TRUE, prob = Px)
             y.ind <- sample(x = marginals$xy[,2], size = 10000, replace = TRUE, prob = Py)
             w.ind <- which(x.ind < y.ind) 
             xy <- data.frame(cbind(x.ind[w.ind], y.ind[w.ind]))   # sample from estiamted marginals
           }, 
           'all_scatter' = { # here plot both biased and original samples 
             xy <- data.frame(biased.data)             
           }
    ) # end switch 
    if(plot.type=='KM_marginal')
      print( ggplot(marginals.CDFs, aes(x=marginals.CDFs[,1], y=marginals.CDFs[,3], color= )) + 
               geom_line(aes(x=marginals.CDFs[,1], y=marginals.CDFs[,3], col="\\hat{F}_X"), size=1.5) + 
               geom_line(aes(x=marginals.CDFs[,2], y=marginals.CDFs[,4], col="\\hat{F}_y"), size=1.5) + 
               geom_line(aes(x=marginals.CDFs[,5], y=marginals.CDFs[,6], col='F_X===F_Y'), size=1.5) + # change line width
               ggtitle(TeX(rep(paste0("$", gsub("_", '-', dependence.type), " (\\theta = ", prms$rho, " )$"), prms$title))) + 
               xlab("X") + ylab("Y") + # keep defaults 
               theme(plot.title = element_text(size=14, face="bold.italic", hjust=0.5),
                     axis.title.y = element_text(face="bold", size=14),
                     axis.title.x = element_text(face="bold", size = 14),
                     axis.text.x = element_text(face="bold", size=12), 
                     axis.text.y = element_text(face="bold", size=12), 
                     legend.background = element_rect(fill = alpha("lightgray", 0)),
                     legend.title = element_text(), 
                     legend.position =  c(0.85, 0.2)) + 
               labs(color = "") + 
               scale_color_discrete(labels=lapply(sprintf(c("\\hat{F}_X", "\\hat{F}_Y", "F_X=F_Y")), TeX)) )
    else # new: plot two scatters on same datasets 
      print( ggplot(xy, aes(x=xy[,1], y=xy[,2])) + 
               #               geom_point(aes(x=xy[,1], y=xy[,2], col="original")) + 
               #               geom_point(shape=3, aes(x=xy[,3], y=xy[,4], col="biased")) + 
               geom_point() + 
               ggtitle(TeX(rep(paste0( gsub("_", '-', dependence.type), " $(\\theta = ", prms$rho, " )$"), prms$title))) + 
               xlab("X") + ylab("Y") +
               theme(plot.title = element_text(size=14, face="bold.italic", hjust=0.5),
                     axis.title.y = element_text(face="bold", size=14),
                     axis.title.x = element_text(face="bold", size = 14),
                     axis.text.x = element_text(face="bold", size=12), 
                     axis.text.y = element_text(face="bold", size=12), 
                     legend.position = "none") )
    
    ggsave(paste0(output.scatter.file, '_', plot.type, '.png'), units="in", width=5, height=5, dpi=300)
    if(dependence.type != 'UniformStrip') # last two plots are only for first dataset  
      break
  } # end loop on plot type
}


###################################################################################################
# Compute density of product of two Gaussian densities  
# Parameters:
# mu1, mu2 - mean vectors 
# sigma1, sigma2 - covariance matrices
# 
# Output:
# mu - mean vector of product distribution
# sigma - covariance matrix of product distribution
###################################################################################################
GaussianDensityProduct <- function(mu1, mu2, sigma1, sigma2)
{
  sigma <- solve(solve(sigma1) + solve(sigma2))
  mu <- sigma * (solve(sigma1) * mu1 + solve(sigma2) * mu2)
  return(list(mu=mu, sigma=sigma))
}  


# Use to keep marginals always sorted by x and y
sort_marginals <- function(marginals)
{
  marginals_sorted <- marginals
  marginals_sorted$xy = apply(marginals$xy, 2, sort)
  marginals_sorted$CDFs = apply(marginals$CDFs, 2, sort)
  marginals_sorted$PDFs[,1] = marginals_sorted$PDFs[order(marginals$CDFs[,1]),1] 
  marginals_sorted$PDFs[,2] = marginals_sorted$PDFs[order(marginals$CDFs[,2]),2] 

  return(marginals_sorted)  
}  
  
  
