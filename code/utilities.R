##################################################################################################
# Compute the modified Hoeffding's test statistic, for the permutation test
# Parameters: 
# data - n*2 matrix with (x,y) sample
# grid_points - all possible (x_i,y_j) points  
# null_distribution - n*n matrix with F_x*F_Y * W estimate OR 4*n mass table with pre-computed mass estimates 
# 
#  Quardants convension:
#   4 | 1
#   ------
#   3 | 2
##################################################################################################
Compute_Statistic<- function(data, grid_points, null_distribution)
{
  epsilon = 0.00000001
  num_samples=dim(data)[1]
  Obs<-matrix(0,4,1) # observed
  Exp<-matrix(0,4,1) # expected
  perm_flag <- !(min(dim(null_distribution))==num_samples) # check type of null. Doesn't work for n=4
  if(!perm_flag)
  { # New! compute CDF for botstrap. Current implementation ignores ties !
    null_distribution_CDF <- PDF_to_CDF_2d(null_distribution, data)
  }
  
  Statistic <- 0 # = rep(0, dim(grid_points)[1])
  for (i in 1:dim(grid_points)[1] ) 
  {
    if(perm_flag) # here we already computed the null in a 4*n table 
    {
      Exp = pmax(null_distribution[i,], epsilon)
    }
    else {        #UP right (Expected)
      for(j in 1:3)
      {
#        Exp[j] = Get_Quarter_Approximated_Mass(grid_points[i,], j, data, null_distribution) # old version uses PDF 
        Exp[j] <- Get_Quarter_Approximated_Mass(grid_points[i,], j, data, null_distribution_CDF)
      }
      Exp[4] = max(1-sum(Exp[1:3]), epsilon)
      Exp <- num_samples*Exp
    }
    #Up right quarter (Observed) - this is the computationally costly part for permutations test 
    R_X <- data[,1]>=grid_points[i,1]
    R_Y <- data[,2]>=grid_points[i,2]
    Obs[1] <- sum(R_X*R_Y)
    Obs[2] <- sum(R_X)-Obs[1]
    Obs[4] <- sum(R_Y)-Obs[1]
    Obs[3] <- num_samples-sum(Obs[c(1,2,4)])
    
    if(any(is.na(Exp)))
      print(paste("i=", i, " Exp=", min(Exp), ' Obs=', min(Obs)))
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
MCMC_Permutations<-function(W, B, N)
{ 
  # Set mcmc sampling parameters 
  burn_in = 2*N
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
    
    if(rbinom(1, 1, min(1,ratio))) #we accept the transition with probability min(1, ratio)
    {
      temp <- Perm[i] # SWAP 
      Perm[i] <- Perm[j]
      Perm[j] <- temp
      if(ctr==burn_in || (ctr%%Cycle==0 && ctr>burn_in))
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
# Estimate the null distribution f_x*f_y*W (given the estimated PDFs f_x, f_y)
# 
# Parameters: 
# pdfs - marginal distributions f_x, f_y probabilities 
# w - n*n matrix with weights W[i,j] = w(x_i, y_j) 
###################################################################################
get_null_distribution <- function(pdfs, W)
{
  #compute the normalizing factor under the null:
  null_distribution <- W * (pdfs[,1] %*% t(pdfs[,2]))
  Z <- sum(null_distribution)
  null_distribution<-null_distribution/Z
  return( list(null_distribution=null_distribution, normalizing_factors=Z) )
}

############################################################################################
# New: draw a bootstrap sample, given the estimated null distribution F_x, F_Y, W
# Use rejection sampling. No need for computing n*n table 
# (Problem: what if prob.(rejection) close to 1?)
# Parameters: 
# data - n*2 array with (X,Y) samples
# pdfs - f_x and f_y  
# bias_type - W
# prms - for w max 
############################################################################################
Bootstrap <- function(data, pdfs, bias_type, prms)
{
  n = dim(data)[1]
  boot_sample<-matrix(0,n,2)
  k <- 0
  xy <- c(0,0)  
  while(k<=n)
  {
    xy[1] <- data[sample(n, 1, prob=pdfs[,1]),1] # Sample X from F_X
    xy[2] <- data[sample(n, 1, prob=pdfs[,2]),2] # Sample Y sample from F_Y
    keep <- rbinom(1, 1, biased.sampling.w(xy[1], xy[2], bias_type)/prms$W_max)
    if(keep) 
    {
      boot_sample[k,] <- xy
      k <- k+1
    }
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
# null_distribution_CDF - n*n table of 2d cumulative of F_x*F_y*w (problem with ties)
###################################################################################################
Get_Quarter_Approximated_Mass<-function(Point, QId, data, null_distribution_CDF)
{
  if(QId %in% c(1,2))
  {
    idx_x <- which(data[,1]>=Point[1])
    idx_x <- idx_x[which.min(data[idx_x,1])]
  } else
  {
    idx_x <- which(data[,1]<Point[1])
    idx_x <- idx_x[which.max(data[idx_x,1])]
  }
  if(QId %in% c(1,4))
  {
    idx_y <- which(data[,2]>=Point[2])
    idx_y <- idx_y[which.min(data[idx_y,2])]
  } else
  {
    idx_y <- which(data[,2]<Point[2])
    idx_y <- idx_y[which.max(data[idx_y,2])]
  }

  if(isempty(idx_x) | isempty(idx_y))
    return(0)    
    
  m <- which.max(data[,1]) # dim(null_distribution_CDF)[1]
  n <- which.max(data[,2]) # dim(null_distribution_CDF)[2]
  
  if(QId == 1)
    S <- null_distribution_CDF[m, n] + null_distribution_CDF[idx_x, idx_y] - 
      null_distribution_CDF[idx_x, n]  - null_distribution_CDF[m, idx_y]
  if(QId == 2)
    S <- null_distribution_CDF[m, idx_y]  - null_distribution_CDF[idx_x, idx_y]  
  if(QId == 3)    
    S <- null_distribution_CDF[idx_x, idx_y]  
  if(QId == 4)
    S <- null_distribution_CDF[idx_x, n]  - null_distribution_CDF[idx_x, idx_y]  

  return(S)       
}



#computing Expect[Q_i(p_j)] for 1<=i<=4, and all j, given a grid of points
# (New, faster implementation)
# Parameters: 
# data - 2*n array (X,Y)
# Permutations - set of permutations
# grid_points - centers of partitions
#
# Output: 
# mass_table - a 4*#grid-points table with quadrants expectation
evaluate_mass_table_by_permutations <- function(data, Permutations, grid_points)
{
  num_permutations <- dim(Permutations)[2]
  n <- dim(data)[1]
  Pij <- matrix(0, n, n) # matrix with P[i]=j estimate
  for(i in 1:num_permutations)
    Pij[cbind(1:n, Permutations[,i])] <-  Pij[cbind(1:n, Permutations[,i])]+1 # need to vector indices here  
  Pij <- Pij / num_permutations # normalize 

  #next each point has its probability of being selected we evaluate Expect[Q_i(p_j)] by summation
  mass_table<-matrix(0, dim(grid_points)[1], 4)
  for(i in seq(1, dim(grid_points)[1],1))
  {
    x<-grid_points[i,]
    mass_table[i,1]<-sum(Pij[data[,1]>=x[1], data[,2]>=x[2]])
    mass_table[i,2]<-sum(Pij[data[,1]>=x[1], data[,2]<x[2]])
    mass_table[i,3]<-sum(Pij[data[,1]<x[1], data[,2]<x[2]])
    mass_table[i,4]<-sum(Pij[data[,1]<x[1], data[,2]>=x[2]])
  } 
  
  mass_table<-mass_table+0.000001
  return(mass_table)
}




# New: cumulative 2d (for faster calculation of test statistic)
# When we have ties we need to correct 
PDF_to_CDF_2d <- function(pdf_2d, data)
{
  P_x <- sort(data[,1], index.return=TRUE)  # Permute to order x_i, y_i 
  P_y <- sort(data[,2], index.return=TRUE)
  
  cdf_2d <- apply(pdf_2d[P_x$ix, P_y$ix], 2, cumsum)
  cdf_2d <- apply(cdf_2d, 1, cumsum)
  
# Deal with ties (not working yet)
#  ties_x <- which(data[-1,1] == head(data[,1], -1))
#  ties_y <- which(data[-1,2] == head(data[,2], -1))
#  print(paste0('ties_x=', ties_x))
#  print(paste0('ties_y=', ties_y))
#  for(i in rev(ties_y))
#    cdf_2d[,i] <- cdf_2d[,i+1]
#  for(i in rev(ties_x))
#    cdf_2d[i,] <- cdf_2d[i+1,]
  
  
  return( t(cdf_2d[invPerm(P_y$ix), invPerm(P_x$ix)]))  
#    return( t(apply(apply(pdf_2d, 2, cumsum), 1, cumsum)) )
    
}




