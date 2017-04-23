compute_statistic_based_permutaions <- function(data, statistic_type, num_of_observation, bias_method, dependence_type, grid, mass_table)
{
  if(statistic_type == 'HHG') #compute the HHG test statistic
  {
    d=dim(mass_table)
    Statistic <- matrix(0,d[1],1);
    P_expected <- matrix(0, 4, 1); 
    
    for (i in 1:d[1] ) # loop on data points 
    {
      switch(dependence_type,
             "Gaussian"={
               for (j in 1:4)
                 P_expected[j] <- mass_table[i, j]; 
             })
      
      P_empiric <- matrix(0, 4, 1);
      for (j in 1:4)
      {
        j1 <- (-1)^ceil(j/2); j2 <- (-1)^round(abs((j-2.5)/2)); # use to flip sign
        P_empiric[j] <- length(which(j1*data[,1]<=j1*grid[i,1] & j2*data[,2]<=j2*grid[i,2]));
        Statistic[i] <- Statistic[i] + (P_empiric[j] - num_of_observation*P_expected[j])^2 / (num_of_observation*P_expected[j]); 
      }
    }
    
  }
  T = sum(Statistic)
  return(T)
}


# New: Compute expected statistics using all permutations 
# Output: 
# Q - an nXn matrix such that Q[i,j] is the probability for a random permutation P to have P[i] = j. 
# 
compute_permutations_marginal_table <- function(PermutationsTable, perm_weights)
{
  n <- dim(PermutationsTable)[1]; # get number of data points 
  num_perms <- dim(PermutationsTable)[2]; # get number of permutations
  Q <- matrix(0, n, n); 

    # loop over all permutations and compute marginals 
  for (p in 1:num_perms)
    for (i in 1:n ) 
    {
      Q[i,PermutationsTable[i,p]] = Q[i,PermutationsTable[i,p]] + perm_weights[p];  # update marginal probabilities 
    }
  Q <- Q/rowSums(Q); # divide each row by its sum. 
  
  return(Q); 
}


# Compute the expected marginal in each quartile 
# 
# Input: 
# data - data set (matrix of size nX2)
# Q - matrix of marginal probabilities 
# X - pivots points which seperates plane into 4 quartiles 
#
compute_statistic_based_permutaions_quartiles <- function(data, Q, X)
{
  n_grid <- dim(X)[1]; 
  n <- dim(data[1]);
  
  P_expected <- matrix(0, n_grid, 4); 
  for (i in 1:n_grid)
    for(j in 1:4)
    {
      j1 <- (-1)^ceil(j/2); j2 <- mod(j,2); # use to flip sign
      Quartile_indices1 <- which(j1*data[,1]<=j1*X[i,1]); # find indices in quartile 
      Quartile_indices2 <- which(j2*data[,2]<=j2*X[i,2]); 
      
      for( i1 in Quartile_indices1)
        for (i2 in Quartile_indices2)
          P_expected[i,j] <- P_expected[i,j] + Q[Quartile_indices1[i1], Quartile_indices2[i2]]; 
    }
  return(P_expected / n); # Normalize to get total expectation 1 (not n)
}




