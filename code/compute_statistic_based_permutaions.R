compute_statistic_based_permutaions <- function(data, statistic_type, num_of_observation,bias_method, grid, mass_table)
{
  epsilon = 0.000001
  if(statistic_type == 'HHG') #compute the HHG test statistic
  {
    d=dim(mass_table)
    Statistic = matrix(0,d[1],1);
    
    for (i in 1:d[1] ) 
    {
      #UP right
      P1 = mass_table[i,1]
      #Down right  
      P2 = mass_table[i,2]
      #Down left  
      P3 = mass_table[i,3]
      #Up left  
      P4 = ifelse(mass_table[i,4] >0,mass_table[i,4], epsilon) 
               
      #Up right quarter
      EmpP1 = length(which(data[,1]>=grid[i,1] & data[,2]>=grid[i,2]));
      #Down right quarter
      EmpP2 = length(which(data[,1]>=grid[i,1] & data[,2]<=grid[i,2]));
      #Down left quarter
      EmpP3 = length(which(data[,1]<=grid[i,1] & data[,2]<=grid[i,2]));
      #Up left quarter
      EmpP4 = length(which(data[,1]<=grid[i,1] & data[,2]>=grid[i,2]));
      Statistic[i] = (EmpP1-P1)^2/(P1) + (EmpP2-P2)^2/(P2)+(EmpP3-P3)^2/(P3) + (EmpP4-P4)^2/(P4)
    }
  }
  T = sum(Statistic)
  return(T)
}
