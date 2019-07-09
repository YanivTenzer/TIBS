Yanivs_estimate_marginals<-function(data, bEstimate_Marginals)
{
  source('get_MLE_marginal_estimate.R')
  if(bEstimate_Marginals)
  {
    F_1<-ecdf(data[,1])
    F_2<-ecdf(data[,2])
    F_x<-(F_1(data[,1])+F_2(data[,1]))/2
    F_Y<-(F_1(data[,2])+F_2(data[,2]))/2
    CDF_table<-cbind(F_x,F_Y)
    
    PDF_table<-get_marginals_PDF(data, CDF_table)
  }
  
  output<-list(CDF_table,PDF_table)
  names(output)<-c('CDFs', 'PDFs')
  return(output)
}
