MCMC_Permutations<-function(Dat,W,TargetSampleSize,N)
{  
  BurningTime = 2*dim(Dat)[1];
  Cycle = dim(Dat)[1];

  Counter = 1;
  Idx = 1;
  PermutationsTable = matrix(0,N,TargetSampleSize);
  oldPerm = 1:N;

  while(Idx<=TargetSampleSize)
  {
    #Here we implement a MHS algorithm with target stationary distribution \pi
    #Choose the two indices to be switched
    #note that this way we are choosing a neighbor at random with probability of 1/(n*(n-1))
    newPerm = oldPerm;
    switchIdx = sample(1:N, 2, replace = TRUE)
  
    #Should we accept the new permutation ?
    i = switchIdx[1];
    j = switchIdx[2];
    ratio = W[i,oldPerm[j]]*W[j,oldPerm[i]]/W[i,oldPerm[i]]*W[j,oldPerm[j]];
    p_old_to_new = min(1,ratio)
    #Toss a coin with probability p_old_to_new:
    u = runif(1,0,1)
  
    if(u<p_old_to_new)#we accept the transition:
    {
      newPerm[switchIdx[1]] = oldPerm[switchIdx[2]];
      newPerm[switchIdx[2]] = oldPerm[switchIdx[1]];
    
      if(Counter==BurningTime || (Counter%%Cycle==0 && Counter>BurningTime))
      {
        PermutationsTable[,Idx]=newPerm;
        Idx = Idx+1;
      }
    
      Counter = Counter+1;
      oldPerm = newPerm;
      
    }
  }
  return(PermutationsTable)
}