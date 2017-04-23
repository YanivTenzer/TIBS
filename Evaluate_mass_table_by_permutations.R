#This function get a list of points(grid_points) and compute Expect[Q_i(p_j)] for 1<=i<=4, and all j
#For more detalis see Or's comment in this paper
Evaluate_mass_table_by_permutations<-function(data, Permutations, num_of_permutations, grid_points)
{
  t<-meshgrid(unique(data[,1]), unique(data[,2]))
  x<-as.vector(t$X)
  y<-as.vector(t$Y)  
  temp_grid<-cbind(x,y)
  Prob<-matrix(0,dim(temp_grid)[1],1)
  for(k in seq(1,dim(temp_grid)[1],1))
  {
    p<-temp_grid[k,]
    i<-which(data[,1]==p[1])
    j<-which(data[,2]==p[2])
    #The probability of this point to be selected when choosing randomly N points
    Prob[k]<-sum(apply(Permutations, 2, function(x) ifelse(x[i]==j,1,0)))/(num_of_permutations)
  }
  #Now, when each point has its probability of being selected we evaluate Expect[Q_i(p_j)] by summation
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

