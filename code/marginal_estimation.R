# Estimate marginals F_X, F_Y from biased sample with density proportional to F_XY * W
# Parameters: 
# data - 2*n sample (X,Y)
# bias_method - w function 
# 
main_estimate_marginals<-function(data, bias_method)
{
  if(bias_method %in% c('sum', 'sum_coordinates', 'exponent_minus_sum_abs')) # for w(x,y)>0 cases 
  { #case 1: strictly positive W, use ML estimator
#    print("Estimate positive marginals")
    w_inv <- 1/biased.sampling.w(data[,1], data[,2], bias_method)
    f_x <- f_y <- w_inv / sum(w_inv) # normalize 
    PDF_table = as.data.frame(cbind(f_x, f_y))
    CDF_table = NULL
  }
  else{ #case 2: left truncation, use the estimator of proposition 1 in the paper for exchangable distributions
#    print("Estimate marginals with left truncation and exchangable distribution")
    F_1<-ecdf(data[,1])
    F_2<-ecdf(data[,2])
    F_x<-(F_1(data[,1])+F_2(data[,1]))/2  # F_X, F_Y are the same CDFs evaluated at different data x,y
    F_Y<-(F_1(data[,2])+F_2(data[,2]))/2   
    CDF_table<-cbind(F_x,F_Y)
    PDF_table<-get_marginals_PDF(data, CDF_table)
  }
  
  return( list(CDFs=CDF_table, PDFs=PDF_table) )
}

##############################################################################
# 
# Parameters: 
# data - 2*n array of (X,Y)
# CDF_table - vector of CDF of X and Y
get_marginals_PDF<-function(data, CDF_table)
{
  num_variables<-dim(CDF_table)[2]
  num_samples<-dim(CDF_table)[1]
  PDF_table<-matrix(0, num_samples, num_variables)

#  print(paste0("num_samples=", num_samples, " and ", dim(data)[1]))
  for(i in 1:num_variables)
  {
    sorted_CDF<-sort(CDF_table[,i], index.return=TRUE)
    temp<-c(sorted_CDF$x[1], sorted_CDF$x[-1]-sorted_CDF$x[-num_samples])
    for(j in 1:num_samples)
    {
      PDF_table[j,i]<-temp[which(sorted_CDF$ix==j)]
#      tryCatch({PDF_table[j,i]<-temp[which(sorted_CDF$ix==j)]}, error = function(e){browser()} )
    }
  }
  return(PDF_table)
}
