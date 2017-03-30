Get_Hyperplane_Truncated_BG_Likelihood <- function(Data, rho, Hyperplane_PRMS)
{
  d=dim(Data)
  Likelihood<-matrix(0,d[1],1)
  b<-Hyperplane_PRMS[3]
  w_1<-Hyperplane_PRMS[1]
  w_2<-Hyperplane_PRMS[2]
  c_rho <- (-b)/sqrt((1+rho)/2*(w_1+w_2)^2+(1-rho)/2*(w_1-w_2)^2)
  Z=1-pnorm(c_rho)
  
  sigma<-matrix(1,2,2)
  sigma[1,2]=rho
  sigma[2,1]=rho
  
  Likelihood<-dmvnorm(Data, c(0,0), sigma, log=FALSE)/Z
  
  return(Likelihood)
}
######Estimate ML PRMS:
Get_Hyperplane_Truncated_BG_ML_PRMS <- function(Data,Hyperplane_PRMS,MAX_ITERATIONS,Epsilon,step_size)
{
  b<-Hyperplane_PRMS[3]
  w_1<-Hyperplane_PRMS[1]
  w_2<-Hyperplane_PRMS[2]
  
  Hyperplane_PRMS=c(w_1,w_2,b)
  RhoGrid = seq(-0.99, 0.99, 0.01)
  LOG_LIKELIHOOD = matrix(0,length(RhoGrid),1)
  
  for(i in seq(1,length(RhoGrid),1))
  {
    LOG_LIKELIHOOD[i] = sum(log(Get_Hyperplane_Truncated_BG_Likelihood(Data,RhoGrid[i],Hyperplane_PRMS)))
  }
 
  Idx = which(LOG_LIKELIHOOD==max(LOG_LIKELIHOOD))

  return(c(RhoGrid[Idx],LOG_LIKELIHOOD[Idx]))
  # OldRho=runif(1)
  # DataLikelihood = Get_Hyperplane_Truncated_BG_Likelihood(Data,OldRho,Hyperplane_PRMS)
  # DataLikelihood = max(DataLikelihood, exp(-500))
  # oldLOG_LIKELIHOOD = sum(log(DataLikelihood))
  # 
  # Delta = Inf
  # iter = 1
  # 
  # while(iter<=MAX_ITERATIONS && Delta>=Epsilon)
  # {
  # 
  #   rho = OldRho
  #   PRMS<-list(rho,Hyperplane_PRMS)
  #   names(PRMS)<-c("rho","Hyperplane_PRMS")
  # 
  #   Gradient<-Get_Gradient(Data,'Hyperplane_Truncated_BG', PRMS)
  #   NewRho<-OldRho+step_size*Gradient
  #   
  #   DataLikelihood = Get_Hyperplane_Truncated_BG_Likelihood(Data,NewRho,Hyperplane_PRMS)
  #   DataLikelihood = max(DataLikelihood, exp(-500))
  #   newLOG_LIKELIHOOD=sum(log(DataLikelihood))
  #   
  #   #browser()
  #   GAIN = newLOG_LIKELIHOOD - oldLOG_LIKELIHOOD
  #   Delta=abs(NewRho-OldRho)
  #   if(is.na(GAIN))
  #   {
  #     browser()
  #   }
  #   if(GAIN<=0)
  #   {
  #     break;
  #   }
  #   else
  #   {
  #     oldLOG_LIKELIHOOD = newLOG_LIKELIHOOD
  #     OldRho = NewRho
  #     iter=iter+1
  #   }
  # 
  # }
  #print(paste('iter=',toString(iter)))
  #return(c(OldRho,oldLOG_LIKELIHOOD))
}
######Calculate gradient:
Get_Gradient<-function(Data,Type,PRMS)
{
  switch(Type,
   
   'Hyperplane_Truncated_BG'=
    {
     d=dim(Data)    
     
     Statistics<-matrix(0,d[1],3)
     Statistics[,1]=Data[,1]^2
     Statistics[,2]=Data[,2]^2
     Statistics[,3]=Data[,1]*Data[,2]
     
     rho = PRMS$rho
     w_1=PRMS$Hyperplane_PRMS[1]
     w_2=PRMS$Hyperplane_PRMS[2]
     b=PRMS$Hyperplane_PRMS[3]
     
     c_rho <- (-b)/sqrt((1+rho)/2*(w_1+w_2)^2+(1-rho)/2*(w_1-w_2)^2)
     
     Gradient<-sum(1/(sqrt(2*pi))*(rho/(1-rho)^2 + Statistics[,3]/sqrt((1-rho^2)) - rho*(Statistics[,1]+Statistics[,2])/(2*(1-rho^2)^(3/2))
     +rho^2*Statistics[,3]/(1-rho^2)^(3/2)
     +b/4*(dnorm(c_rho)/pnorm(c_rho))*( ((w_1+w_2)^2+(w_1-w_2)^2)/((1+rho)/2*(w_1+w_2)^2+(1-rho)/2*(w_1-w_2)^2)^(3/2)) ))
     },
   {
     print('Error: unrecognaized density function specified')
    })#END_OF_SWITCH
}