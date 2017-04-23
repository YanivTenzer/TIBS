Bootstrap_according_to_null<-function(data, null_distribution)
{
  library(gdata)
  n = dim(data)[1]
  dummy_sample<-matrix(0,n,2)
  #next we sample n observations from the pdf (multinomial) distribution:
  idx<- rmultinom(n=dim(data)[1], size=1, prob=unmatrix(null_distribution,byrow=T))
  d <- matrix(idx, nrow = dim(data)[1], byrow = TRUE)
  y_idx<-apply(d, 1, function(x) which(x==1)%%dim(data)[1])
  y_idx[y_idx==0]=dim(data)[1]
  x_idx<-ceiling((apply(d, 1, function(x) which(x==1))-y_idx)/dim(data)[1])
  x_idx=x_idx+1
  dummy_sample<-cbind(data[x_idx,1], data[y_idx,2])
  return(dummy_sample)
}