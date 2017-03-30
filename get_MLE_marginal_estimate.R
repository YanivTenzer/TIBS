get_MLE_marginal_estimate<-function(data, w_type,bias_params)
{
  num_of_variables<-dim(data)[2]
  num_of_samples<-dim(data)[1]
  #browser()
  marginals_cdf_0<-matrix(0,num_of_samples, num_of_variables)
  
  if(identical(w_type, 'Hyperplane_Truncation'))
  {
    Hyperplane_prms = bias_params$Hyperplane_prms
    S = apply(data,1,function(x) (Hyperplane_prms[1]*x[1]+Hyperplane_prms[2]*x[2]+Hyperplane_prms[3]))
    W<-vector(mode = 'double', length=num_of_samples)
    W[which(sign(S)>0)] = 1
    W[which(sign(S)<=0)] = 0
  }
  
  Z<-sum(W^(-1))
  for(i in seq(1, num_of_variables,1))
  {
    for(j in seq(1,num_of_samples,1))
    {
      marginals_cdf_0[j,i]<-sum(ifelse(data[,i]<=data[j,i], 1, 0)*W^(-1))/Z  
    }
  }
  
  return(marginals_cdf_0)       
}
##############################################################################
get_marginals_PDF<-function(data, CDF_table)
{
  num_of_variables<-dim(CDF_table)[2]
  num_of_samples<-dim(CDF_table)[1]
  
  PDF_table<-matrix(0, dim(CDF_table)[1], dim(CDF_table)[2])
  density_table<-matrix(0, dim(CDF_table)[1], dim(CDF_table)[2])
  
  for( i in seq(1,num_of_variables,1))
  {
    sorted_CDF<-sort(CDF_table[,i], index.return=TRUE)
    d<-sorted_CDF$x[-1]-sorted_CDF$x[-num_of_samples]
    temp<-c(d[1],sorted_CDF$x[-1]-sorted_CDF$x[-num_of_samples])
    for(j in seq(1,dim(data)[1],1))
    {
      PDF_table[j,i]<-temp[which(sorted_CDF$ix==j)] 
    }
  }
  return(PDF_table)
}