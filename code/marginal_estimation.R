##############################################################################
# Estimate marginals Fx, Fy from biased sample with density proportional to Fxy * W
# Parameters: 
# data - 2*n sample (X,Y)
# w.fun - biased sampling w function 
#
# Output: 
# xy - data
# CDFs - cdf of the estimated marginals
# PDFs - pdfs of the estimated marginals 
##############################################################################
EstimateMarginals<-function(data, w.fun)
{
  if(w.fun %in% c('sum', 'sum_coordinates', 'exponent_minus_sum_abs')) # for w(x,y)>0 cases 
  { #case 1: strictly positive W, use ML estimator
    w.inv <- 1/w_fun_eval(data[,1], data[,2], w.fun)
    Fx <- Fy <- w.inv / sum(w.inv) # normalize 
    PDF.table = as.data.frame(cbind(Fx, Fy))
    CDF.table = NULL
  } else{ 
    if(w.fun %in% c('survival'))  # what is the definition of w here? 
    {
      # Estimate marginals using Kaplan-Meier estimator 
      n <- dim(data)[1]
      require(survival)
      y.srv <- Surv(time=data[,1], time2=data[,2], event = rep(1,n))
      Fy <- survfit(y.srv~1)
      x.srv <- Surv(time=-data[,2], time2=-data[,1], event = rep(1,n))
      Fx <- survfit(x.srv~1)
      CDF.table<-cbind(rev(Fx$surv), 1-Fy$surv)
      data <- cbind(rev(-Fx$time), Fy$time)
      PDF.table<-GetMarginalsPDF(data, CDF.table)
    } else {  
      if(w.fun %in% c('naive', 'no_bias')) # no bias (assume W(x,y)=1)
      {
        Fx<-ecdf(data[,1])
        Fy<-ecdf(data[,2])
        CDF.table<-cbind(Fx(data[,1]),Fy(data[,2]))
        PDF.table<-GetMarginalsPDF(data, CDF.table)
      } else {  # general biased sampling function (w.fun can be a function not a string)
        #case 2: left truncation, use the estimator of proposition 1 in the paper for exchangable distributions
        # Augment data to include both x and y values for each axis (due to symmetry)
        augment.data <- 1 # new: add xy values 
        if(augment.data) # duplicate x and y values 
          data <- cbind(union(data[,1], data[,2]), union(data[,1], data[,2]))
        F1<-ecdf(data[,1])
        F2<-ecdf(data[,2])
        Fx<-(F1(data[,1])+F2(data[,1]))/2  # Fx, Fy are the same CDFs evaluated at different data x,y
        Fy<-(F1(data[,2])+F2(data[,2]))/2   
        CDF.table<-cbind(Fx,Fy)
        PDF.table<-GetMarginalsPDF(data, CDF.table)
      }
    } # end if 
  }
  return( list(xy=data, CDFs=CDF.table, PDFs=PDF.table) ) # new: return also x,y (might be different than original)
}

##############################################################################
# Estimate marginal PDF from data 
# Parameters: 
# data - 2*n array of (X,Y)
# CDF.table - vector of CDF of X and Y
##############################################################################
GetMarginalsPDF<-function(data, CDF.table)
{
  num.variables<-dim(CDF.table)[2]
  num.samples<-dim(CDF.table)[1]
  PDF.table<-matrix(0, num.samples, num.variables)
  
  for(i in 1:num.variables)
  {
    sorted.CDF<-sort(CDF.table[,i], index.return=TRUE)
    temp<-c(sorted.CDF$x[1], sorted.CDF$x[-1]-sorted.CDF$x[-num.samples])
    for(j in 1:num.samples)
      PDF.table[j,i]<-temp[which(sorted.CDF$ix==j)]
  }
  return(PDF.table)
}
