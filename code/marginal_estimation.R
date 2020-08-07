##############################################################################
# Estimate marginals Fx, Fy from biased sample with density proportional to Fxy * W
# Parameters: 
# data - 2*n sample (X,Y)
# w.fun - biased sampling w function 
# prms - new (optional): what estimation algorithm to use 
#
# Output: 
# xy - data
# CDFs - cdf of the estimated marginals
# PDFs - pdfs of the estimated marginals 
##############################################################################
EstimateMarginals <- function(data, w.fun, prms=c())
{
  if(!missing(prms))     # get marginal CDFs from PDFs
    PDF.table <- iterative_marginal_estimation(data, w.fun)

  if(w.fun %in% c('sum', 'sum_coordinates', 'exponent_minus_sum_abs')) # for w(x,y)>0 cases 
  { #case 1: strictly positive W, use ML estimator
    w.inv <- 1/w_fun_eval(data[,1], data[,2], w.fun)
    Fx <- Fy <- w.inv / sum(w.inv) # normalize
    PDF.table = cbind(Fx, Fy) # PDF for x and y is the same 
    CDF.table = PDFToCDFMarginals(data, PDF.table)
  } else { 
    if(w.fun %in% c('survival'))  # what is the definition of w here? 
    { # Estimate marginals using Kaplan-Meier estimator 
      n <- dim(data)[1]
      require(survival)
      y.srv <- Surv(time=data[,1], time2=data[,2], event = rep(1,n))
      x.srv <- Surv(time=-data[,2], time2=-data[,1], event = rep(1,n))
      Fx <- survfit(x.srv~1)
      Fy <- survfit(y.srv~1)
#      print("Survival CDF!!")
      CDF.table <- cbind(rev(Fx$surv), 1-Fy$surv)
      data <- cbind(rev(-Fx$time), Fy$time)
    } else {  
      if(w.fun %in% c('naive', 'no_bias')) # no bias (assume W(x,y)=1)
      {
        Fx <- ecdf(data[,1])
        Fy <- ecdf(data[,2])
        CDF.table <- cbind(Fx(data[,1]),Fy(data[,2]))
      } else {  # general biased sampling function (w.fun can be a function not a string) with exchangable distributions
        #case 2: left truncation, use the estimator of Proposition 1 in the paper for exchangable distributions
        # Augment data to include both x and y values for each axis (due to symmetry)
#        print("DIM DATA BEFORE")
#        print(dim(data))
        augment.data <- 1 # new: add xy values 
        if(augment.data) # duplicate x and y values 
          data <- cbind(union(data[,1], data[,2]), union(data[,1], data[,2])) # change data size: n-> 2*n 
        
#        print("DIM DATA AFTER")
#        print(dim(data))
        F1 <- ecdf(data[,1])
        F2 <- ecdf(data[,2])
        Fx <- (F1(data[,1])+F2(data[,1]))/2  # Fx, Fy are the same CDFs evaluated at different data x,y
        Fy <- (F1(data[,2])+F2(data[,2]))/2   
        CDF.table<-cbind(Fx,Fy)
      }
    } # end if 
#    print('Estimating PDF marginal! returning!')
    PDF.table <- CDFToPDFMarginals(CDF.table)
  }  # else on w.fun type 
  save(data, PDF.table, CDF.table, file='cdfpdf.Rdata')
  return( list(xy=data, CDFs=CDF.table, PDFs=PDF.table) ) # new: return also x,y (might be different than original)
}


##############################################################################
# Iterative Algorithm for marginal estimation from samples of F_{XY}^{(w)}
# Parameters: 
# data - 2*n array of (X,Y)
# CDF.table - vector of CDF of X and Y
# Output: 
# 
##############################################################################
iterative_marginal_estimation <- function(data, w.fun)  
{
  epsilon <- 0.00000001
  iters <- 1000 
  f_x <- rep(0, n)
  f_y <- rep(0, n)
  w.mat <- w_fun_to_mat(data[,1], data[,2], w.fun) # was w_fun_eval?
  
  for(t in 1:iters)
  {
    f_x_prev <- f_x
    f_y_prev <- f_y
    f_x <- 1 / (w.mat * f_y)
    f_y <- 1 / (t(w.mat) * f_x)
    if(sum(abs(f_x-f_x_prev)) + sum(abs(f_x-f_x_prev)) < epsilon)
      break
  }    
  return(c(f_x, f_y)) # PDF.table 
}  
  
  