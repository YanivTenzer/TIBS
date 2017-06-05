Compute_Statistic<- function(data, statistic_type, grid_points, null_data, null_distribution)
{
  source('Get_Quarter_Approximated_Mass.R')
  num_of_samples=dim(data)[1]
  P<-matrix(0,4,1)
 
  if(statistic_type == 'HHG') #compute the HHG test statistic
  {
    test_idx<-seq(1, dim(grid_points)[1],1)
    Statistic = matrix(0,length(test_idx),1)
    counter = 0
    for (i in 1:length(test_idx) ) 
    {
      idx<-test_idx[i]
      #UP right
      P[1] = Get_Quarter_Approximated_Mass(grid_points[idx,], 1, null_data, null_distribution)
      #Down right 
      P[2] = Get_Quarter_Approximated_Mass(grid_points[idx,], 2, null_data, null_distribution)
      #Down left  
      P[3] = Get_Quarter_Approximated_Mass(grid_points[idx,], 3, null_data, null_distribution)
      #Up left  
      P[4] = 1-(P[1]+P[2]+P[3])
      if(P[4]<=0)
      {
        P4=0.000000001
        counter = counter +1
      }
      
      #Up right quarter
      EmpP1 = length(which(data[,1]>=grid_points[idx,1] & data[,2]>=grid_points[idx,2]));
      #Down right quarter
      EmpP2 = length(which(data[,1]>grid_points[idx,1] & data[,2]<=grid_points[idx,2]));
      #Down left quarter
      EmpP3 = length(which(data[,1]<grid_points[idx,1] & data[,2]<grid_points[idx,2]));
      #Up left quarter
      EmpP4 = length(which(data[,1]<grid_points[idx,1] & data[,2]>=grid_points[idx,2]));
      
        Statistic[i] = (EmpP1-num_of_samples*P[1])^2/(num_of_samples*P[1]) + (EmpP2-num_of_samples*P[2])^2/(num_of_samples*P[2])+ 
                       (EmpP3-num_of_samples*P[3])^2/(num_of_samples*P[3]) + (EmpP4-num_of_samples*P[4])^2/(num_of_samples*P[4]);
      
      
  }
  
    idx1<-which(apply(Statistic, 1, function(x) is.nan(x)|| (x==Inf))==FALSE)
    T = sum(Statistic[idx1])
  
  }
  
  if(identical(statistic_type,'kendall')) #compute the HHG test statistic
  {
    T=cor(data, method='kendall')[1,2]
  }

  return(T)
}
