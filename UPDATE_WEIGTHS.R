UPDATE_WEIGTHS<-function(data, marginals_PDF,bias_params)
{
  Hyperplane_prms = bias_params$Hyperplane_prms
  
  num_of_variables<-dim(marginals_PDF)[2]
  num_of_samples<-dim(marginals_PDF)[1]

  W<-matrix(0,num_of_samples,num_of_variables)
  
  for(i in 1:2)
  {
    for( j in 1:num_of_samples)
    {
      temp<-cbind(rep.int(data[j,i],num_of_samples),data[,-i])
      S<-ifelse(Hyperplane_prms[i]*temp[,1]+Hyperplane_prms[setdiff(c(1,2),i)]*temp[,2]+Hyperplane_prms[3]>0,1,0)
      W[j,i]= sum(S*marginals_PDF[,-i])
    }
  }
  
  return(W)
}