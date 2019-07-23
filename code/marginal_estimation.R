# Estimate marginals Fx, Fy from biased sample with density proportional to FxY * W
# Parameters: 
# data - 2*n sample (X,Y)
# bias.method - string indicating w function 
# 
EstimateMarginals<-function(data, bias.method)
{
  if(bias.method %in% c('sum', 'sum_coordinates', 'exponent_minus_sum_abs')) # for w(x,y)>0 cases 
  { #case 1: strictly positive W, use ML estimator
    #    print("Estimate positive marginals")
    w_inv <- 1/BiasedSamplingW(data[,1], data[,2], bias.method)
    Fx <- Fy <- w_inv / sum(w_inv) # normalize 
    PDF.table = as.data.frame(cbind(Fx, Fy))
    CDF.table = NULL
  }
  else{ 
    if(bias.method %in% c('survival')) 
    {
      # New: estimate marginals using Kaplan-Meier estimator (equivalent to case 2?)
      require(survival)
      y.srv <- Surv(time=data[,1], time2=data[,2], event = rep(1,n))
      Fy <- survfit(y.srv~1)
      x.srv <- Surv(time=-data[,1], time2=-data[,2], event = rep(1,n))
      Fx <- survfit(x.srv~1)
      CDF.table<-cbind(Fx,Fy)
      PDF.table<-GetMarginalsPDF(data, CDF.table)
    } else {
      
      #case 2: left truncation, use the estimator of proposition 1 in the paper for exchangable distributions
      #    print("Estimate marginals with left truncation and exchangable distribution")
      F1<-ecdf(data[,1])
      F2<-ecdf(data[,2])
      Fx<-(F1(data[,1])+F2(data[,1]))/2  # Fx, Fy are the same CDFs evaluated at different data x,y
      Fy<-(F1(data[,2])+F2(data[,2]))/2   
      CDF.table<-cbind(Fx,Fy)
      PDF.table<-GetMarginalsPDF(data, CDF.table)
    }
    
    
    return( list(CDFs=CDF.table, PDFs=PDF.table) )
  }
  
  ##############################################################################
  # 
  # Parameters: 
  # data - 2*n array of (X,Y)
  # CDF.table - vector of CDF of X and Y
  GetMarginalsPDF<-function(data, CDF.table)
  {
    num.variables<-dim(CDF.table)[2]
    num.samples<-dim(CDF.table)[1]
    PDF.table<-matrix(0, num.samples, num.variables)
    
    #  print(paste0("num.samples=", num.samples, " and ", dim(data)[1]))
    for(i in 1:num.variables)
    {
      sorted.CDF<-sort(CDF.table[,i], index.return=TRUE)
      temp<-c(sorted.CDF$x[1], sorted.CDF$x[-1]-sorted.CDF$x[-num.samples])
      for(j in 1:num.samples)
        PDF.table[j,i]<-temp[which(sorted.CDF$ix==j)]
    }
    return(PDF.table)
  }
  