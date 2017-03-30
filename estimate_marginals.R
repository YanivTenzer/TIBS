estimate_marginals<-function(Truncated_Data, dependence_type, Bias_Method, bias_params, bEstimate_Marginals, prms)
{
  source('get_MLE_marginal_estimate.R')
  source('UPDATE_WEIGTHS.R')
  #browser()
  if(bEstimate_Marginals)
  {
    marginals_MLE_CDF<-get_MLE_marginal_estimate(Truncated_Data,'Hyperplane_Truncation',bias_params)
    oldF<-marginals_MLE_CDF
    newF<-oldF
    marginals_PDF<-get_marginals_PDF(Truncated_Data,marginals_MLE_CDF)
    oldPDF<-marginals_PDF
    newPDF<-oldPDF
    
    num_of_iterations = 1000
    epsilon = 10^(-3)
    ######Main loop#######
    for(i in 1:num_of_iterations)
    {
      new_W<-UPDATE_WEIGTHS(Truncated_Data,oldPDF,bias_params)
      #Update CDFs:
      for(j in seq(1,dim(Truncated_Data)[2],1))
      {
        z=sum(new_W[,j]^(-1))
        for(k in seq(1,dim(Truncated_Data)[1],1))
        {
          idx<-which(Truncated_Data[,j]<=Truncated_Data[k,j])
          newF[k,j]= sum(new_W[idx,j]^(-1))/z
        }
      }
      
      newPDF<-get_marginals_PDF(Truncated_Data,newF)
      delta=abs(newF-oldF)
      if(mean(delta)<epsilon)
      {
        break
      }else{
        oldF=newF
        oldPDF=newPDF
      }
    }
    #####End of main loop#####
    newPDF<-get_marginals_PDF(Truncated_Data,newF)
    marginals<-list(newPDF, newF)
    names(marginals)<-c("PDF", "CDF")
    
  }else
  {
    #Here are assume that we know that the CDF-s are standard normal
    switch (dependence_type,
            'Gaussian' = {newF=pnorm(Truncated_Data)},
            'Exponential'={newF=cbind(1-pexp(0.6-Truncated_Data[,1], rate = prms$lambda),pexp(Truncated_Data[,2], rate = prms$lambda))})
    
    newPDF<-get_marginals_PDF(Truncated_Data,newF)
    marginals<-list(newPDF, newF)
    names(marginals)<-c("PDF", "CDF")
  }
  #browser()
  return(marginals)
}