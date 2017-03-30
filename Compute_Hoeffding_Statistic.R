Compute_Hoeffding_Statistic<- function(data, null_distribution)
{
  d<-dim(data)[1]
  T<-0
  for(i in seq(1,d,1))
  {
    for(j in seq(1,d,1))
    {
      point<-c(data[i,1], data[j,2])
      idx_x<-which(data[,1]<=point[1])
      idx_y<-which(data[,2]<=point[2])
    
      s = 0
      for(i in seq(1,length(idx_x),1))
      {
        s<-s+sum(null_distribution[idx_x[i],idx_y])
      }
      
      empirical_cdf<-length(which(data[,1]<=point[1] & data[,2]<=point[2]))/d
      T<-T+(empirical_cdf-s)^2
    }
  }
  browser()
  return(T)
}