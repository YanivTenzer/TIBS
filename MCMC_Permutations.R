# Draw permutations at random accoring to distirubtion determined by matrix W 
# Input: 
# W - an N*N weight matrix 
# num_perms <- number of permtuations to sample
# 
# Output: 
# PermutationsTable - table of sampled permutations
#
MCMC_Permutations<-function(W,num_perms)
{  

  N <- dim(W)[1]
  BurningTime = 1500; # Set default parameters 
  Cycle = 100;

  Counter = 1;
  Idx = 1;
  PermutationsTable = matrix(0,N,num_perms);
  Perm = 1:N; # init with ID permutation 

  while(Idx<=num_perms)
  {
    #Here we implement a MHS algorithm with target stationary distribution \pi
    #Choose the two indices to be switched
    #note that this way we are choosing a neighboring permutation at random with probability of 2/(n*(n-1))
    switchIdx = sample(1:N, 2, replace = FALSE)
  
    #Should we accept the new permutation ?
    i = switchIdx[1];
    j = switchIdx[2];
    ratio = W[i,Perm[j]]*W[j,Perm[i]]/(W[i,Perm[i]]*W[j,Perm[j]]);

    if(runif(1,0,1)<min(1,ratio)) # accept the transition: Toss a coin with prob. min(1,ratio)
    {
      temp <- Perm[i]; Perm[i] <- Perm[j]; Perm[j] <- temp;# swap i and j
      if(Counter==BurningTime || (Counter%%Cycle==0 && Counter>BurningTime))
      {
        PermutationsTable[,Idx]=Perm;
        Idx = Idx+1;
      }
      Counter = Counter+1;
    }
  }
  return(PermutationsTable)
}