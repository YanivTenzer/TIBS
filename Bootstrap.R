Bootstrap<-function(data, null_distribution)
{
  n = dim(data)[1]
  dummy_sample<-matrix(0,n,2)
  
  prob=as.vector(t(null_distribution))
  non_zero_idx<-which(prob!=0)
  prob_2<-prob[non_zero_idx]
  
  a<-sort(prob_2, index.return=TRUE)
  temp_2<-a$x
  key_table<-cbind(prob_2[a$ix],non_zero_idx[a$ix])
  prob2_cum<-cumsum(temp_2)
  prob2_cum<-c(-Inf,prob2_cum, Inf )
  
  rand_idx<-runif(n)
  idx<-matrix(0,length(rand_idx),1)

  for(i in seq(2, length(prob2_cum)-1,1))
  {
    idx[which(prob2_cum[i-1]<rand_idx & rand_idx<prob2_cum[i+1])]=key_table[i-1,2]
  }
  
  #next we sample n observations from the pdf (multinomial) distribution:
  y_idx<-idx%%dim(data)[1]
  y_idx[y_idx==0]=dim(data)[1]
  x_idx<-ceiling((idx-y_idx)/dim(data)[1])
  x_idx=x_idx+1
  dummy_sample<-cbind(data[x_idx,1], data[y_idx,2])
  return(dummy_sample)
}