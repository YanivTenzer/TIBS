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
  #  print("Start Compute")
  epsilon = 0.00000001
  num_samples=dim(data)[1]
  Obs<-matrix(0,4,1) # observed
  Exp<-matrix(0,4,1) # expected
  perm_flag <- !(min(dim(null_distribution))==num_samples) # check type of null. Doesn't work for n=4
  if(!perm_flag)
    null_distribution_CDF <- PDF_to_CDF_2d(null_distribution)
  
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
        Exp[j] = Get_Quarter_Approximated_Mass(grid_points[i,], j, data, null_distribution)
#        Exp2 <- Get_Quarter_Approximated_Mass(grid_points[i,], j, data, null_distribution_CDF)
#        cat(j, 'Exp=', Exp[j], ' Exp2=', Exp2)
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

    # #Up right quarter (Observed)
    # Obs[1] = sum(data[,1]>=grid_points[i,1] & data[,2]>=grid_points[i,2])
    # #Down right quarter
    # Obs[2] = sum(data[,1]>=grid_points[i,1] & data[,2]<grid_points[i,2])
    # #Down left quarter
    # Obs[3] = sum(data[,1]<grid_points[i,1] & data[,2]<grid_points[i,2])
    # #Up left quarter
    # Obs[4] = num_samples-sum(Obs[1:3]) #    sum(data[,1]<grid_points[i,1] & data[,2]>=grid_points[i,2])
    
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
# TargetSampleSize - number of permutations to draw
# N - sample size (can be read from data or W?)
#########################################################################################
MCMC_Permutations<-function(W, TargetSampleSize, N)
{ 
  # Set mcmc sampling parameters 
  burn_in = 2*N;
  Cycle = N;
  
  ctr = 1;
  Idx = 1;
  PermutationsTable = matrix(0,N,TargetSampleSize);
  Perm = 1:N;
  while(Idx<=TargetSampleSize)
  {
    #A MHS algorithm with target stationary distribution \pi
    #Choose the two indices to be switched
    #We are choosing a neighbor at random with probability of 1/(n*(n-1))
    switchIdx = sample(1:N, 2, replace = FALSE) # why replace = TRUE? 
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
      ctr = ctr+1;
    }
  }
  return(PermutationsTable)
}

###################################################################################
# Estimate the null distribution (given the estimated PDFs)
# 
# Parameters: 
# data - n*2 array (x,y)
# pdfs - marginal distributions f_x, f_y
# w - n*n matrix with W[i,j] = w(x_i, y_j) 
###################################################################################
get_null_distribution<-function(data,pdfs,W)
{
  num_samples<-dim(pdfs)[1]
  null_distribution<-matrix(0,num_samples,num_samples)
  #################################
  dx<-matrix(0,num_samples,1)
  dy<-matrix(0,num_samples,1)
  epsilon=0.0001
  x<-c(min(data[,1])-epsilon,unique(data[,1]),max(data[,1])+epsilon)
  y<-c(min(data[,2])-epsilon, data[,2],max(data[,2])+epsilon) # why no unique for y? 
  
  sorted_x<-sort(x,index.return = TRUE)
  idx_x<-sorted_x$ix
  sorted_x<-sorted_x$x
  sorted_y<-sort(y,index.return = TRUE)
  idx_y<-sorted_y$ix
  sorted_y<-sorted_y$x
  
  for(i in seq(2, length(x)-1,1))
  {
    dx[i]<-0.5*(sorted_x[i+1]-sorted_x[i])+0.5*(sorted_x[i]-sorted_x[i-1])
    dy[i]<-0.5*(sorted_y[i+1]-sorted_y[i])+0.5*(sorted_y[i]-sorted_y[i-1])
  }
  dx<-dx[idx_x[-c(1,num_samples+2)]]
  dy<-dy[idx_y[-c(1,num_samples+2)]]
  #################################
  #compute the normalizing factor under the null:
  Z=0
  temp<-dim(pdfs)[1]
  for(i in 1:temp)
  {
    null_distribution[i,]<-W[i,]*pdfs[i,1]*pdfs[,2]
    Z<-Z+sum(W[i,]*pdfs[i,1]*pdfs[,2])
  }
  null_distribution<-null_distribution/Z
  output<-list(null_distribution,Z)
  names(output)<-c("null_distribution", "normalizing_factor")
  return(output)
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
Bootstrap<-function(data, pdfs, bias_type, prms)
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
# 
# Point - 
# QId - 1,2,3,4
# data - 2*n array of (x,y)
# null_distirbution - n*n table of F_x*F_y*w
#
###################################################################################################
Get_Quarter_Approximated_Mass<-function(Point, QId, data, null_distribution)
{
  S<-0
  if(QId %in% c(1,2))
    idx_x <- which(data[,1]>=Point[1])
  else
    idx_x <- which(data[,1]<Point[1])
  if(QId %in% c(1,4))
    idx_y <- which(data[,2]>=Point[2])
  else
    idx_y <- which(data[,2]<Point[2])
  for(i in 1:length(idx_x))
  {
    if(is.na(S))
      print(paste("S=", S))
    S <- S + sum(null_distribution[idx_x[i], idx_y])  
  } 
  
  return(S)       
}

Get_Quarter_Approximated_Mass2<-function(Point, QId, data, null_distribution_CDF)
{
  if(QId %in% c(1,2))
    idx_x <- which.min(data[,1]>=Point[1])
  else
    idx_x <- which.max(data[,1]<Point[1])
  if(QId %in% c(1,4))
    idx_y <- which.min(data[,2]>=Point[2])
  else
    idx_y <- which.max(data[,2]<Point[2])
  
  if(QId == 1)
    S <- null_distribution_CDF[dim(data)[1], dim(data)[2]] + null_distribution_CDF[idx_x, idx_y] - 
      null_distribution_CDF[idx_x, dim(data)[2]]  - null_distribution_CDF[dim(data)[1], idx_y]
  if(QId == 2)
    S <- null_distribution_CDF[dim(data)[1], idx_y]  - null_distribution_CDF[idx_x, idx_y]  
  if(QId == 3)    
    S <- null_distribution_CDF[idx_x, idx_y]  
  if(QId == 4)
    S <- null_distribution_CDF[idx_x, dim(data)[2]]  - null_distribution_CDF[idx_x, idx_y]  
  
  return(S)       
}



############################################################################################################
#computing Expect[Q_i(p_j)] for 1<=i<=4, and all j, given a grid of points
# Parameters: 
# data - 2*n array (X,Y)
# Permutations - set of permutations
# grid_points - centers of partitions
############################################################################################################
evaluate_mass_table_by_permutations<-function(data, Permutations, grid_points)
{
  num_of_permutations <- dim(Permutations)[2]
  t<-meshgrid(unique(data[,1]), unique(data[,2]))
  x<-as.vector(t$X)
  y<-as.vector(t$Y)  
  temp_grid<-cbind(x,y)
  Prob<-matrix(0,dim(temp_grid)[1],1)
  for(k in 1:dim(temp_grid)[1]) # loop on unique x values 
  {
    p<-temp_grid[k,]
    i<-which(data[,1]==p[1])
    j<-which(data[,2]==p[2])
    #The probability of this point to be selected when choosing randomly N points
    if(length(i)>1 || length(j)>1)
    {
      Prob[k]<-sum(apply(Permutations, 2, function(x) ifelse(x[i] %in% j,1,0)))/(num_of_permutations)
    }else
    {
      Prob[k]<-sum(apply(Permutations, 2, function(x) ifelse(x[i]==j,1,0)))/(num_of_permutations)
    }
  }
  #next each point has its probability of being selected we evaluate Expect[Q_i(p_j)] by summation
  mass_table<-matrix(0,dim(grid_points)[1],4)
  for(i in seq(1, dim(grid_points)[1],1))
  {
    x<-grid_points[i,]
    P1<-sum(Prob[which(temp_grid[,1]>=x[1] & temp_grid[,2]>=x[2])])
    P2<-sum(Prob[which(temp_grid[,1]>=x[1] & temp_grid[,2]<=x[2])])
    P3<-sum(Prob[which(temp_grid[,1]<=x[1] & temp_grid[,2]<=x[2])])
    P4<-sum(Prob[which(temp_grid[,1]<=x[1] & temp_grid[,2]>=x[2])])
    mass_table[i,]<-c(P1,P2,P3,P4)
  } 
  
  mass_table<-mass_table+0.000001
  return(mass_table)
}

# New: cumulative null (for faster calculation of test statistic)
PDF_to_CDF_2d <- function(pdf_2d)
{
  return( t(apply(apply(pdf_2d, 2, cumsum), 1, cumsum)) )
}




