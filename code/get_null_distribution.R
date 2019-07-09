get_null_distribution<-function(data,pdfs,W)
{
  num_of_samples<-dim(pdfs)[1]
  null_distribution<-matrix(0,num_of_samples,num_of_samples)
  #################################
  dx<-matrix(0,num_of_samples,1)
  dy<-matrix(0,num_of_samples,1)
  epsilon=0.0001
  x<-c(min(data[,1])-epsilon,unique(data[,1]),max(data[,1])+epsilon)
  y<-c(min(data[,2])-epsilon, data[,2],max(data[,2])+epsilon)
  
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
  dx<-dx[idx_x[-c(1,num_of_samples+2)]]
  dy<-dy[idx_y[-c(1,num_of_samples+2)]]
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