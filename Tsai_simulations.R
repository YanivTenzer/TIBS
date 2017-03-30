cores=detectCores()
cl<-makeCluster(cores[1]-1) #not to overload your computer registerDoSNOW(cl)
registerDoParallel(cl) 

sample_size=500
Hyperplane_prms<-c(-1,1,0)
bias_params<-list(Hyperplane_prms)
names(bias_params)<-c("Hyperplane_prms")
bias_method = 'Hyperplane_Truncation'  
num_of_iterations = 200
dependence_type = 'Gaussian'

for(iter in seq(1,9,1))
{
  switch(dependence_type,
         'Gaussian'={rho =0.1*iter
         prms<-list(rho)
         names(prms)<-c("rho")})       
  #The main loop:
  Pvalues<- foreach(i=1:num_of_iterations, .combine=rbind) %dopar%{ 
   data<-simulate_data(sample_size,'Gaussian', prms) 
   truncated_data<-Create_Bias(data, bias_method, bias_params)
   result<-tsai.test.ties(truncated_data[,1], truncated_data[,2])
   temp<-result[2]
    
   temp
  }
  save(Pvalues, file=paste('Tsai_simulations_results/Pvalues_rho',toString(iter), '.Rdata', sep='' ))
  rm(Pvalues)
}

