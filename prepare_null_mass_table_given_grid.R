prepare_null_mass_table_given_grid<-function(bias_params, grid)
{
  library(pracma)
  
  d<-dim(grid)
  mass_table<-matrix(0, d[1],4)

  for(i in seq(1,d[1],1))
  {
    for(j in seq(1,4,1))
    {
        mass_table[i,j]=get_quarter_null_probability(grid[i,], j,'Hyperplane_Truncation',bias_params)
    }
  }
  return(mass_table)    
}
##########################################
get_quarter_null_probability<-function(Point, QId,bias_method,bias_params)
{
  epsilon = 0.01
  
  switch (bias_method,
          'Hyperplane_Truncation' = {
            HyperplanePRMS<-bias_params$Hyperplane_prms
            alpha1<-HyperplanePRMS[1]
            alpha2<-HyperplanePRMS[2]
            b<-HyperplanePRMS[3]
            
            rho = 0
            c_rho <- (-b)/sqrt((1+rho)/2*(alpha1+alpha2)^2+(1-rho)/2*(alpha1-alpha2)^2)
            Z=1-pnorm(c_rho)
            
            switch(QId,
                   #Q1
                   "1"={
                      t<-meshgrid(seq(-5,5,epsilon), seq(-5,5,epsilon))
                      x<-as.vector(t$X)
                      y<-as.vector(t$Y)  
                      grid_points<-cbind(x,y)
                      idx<-which(alpha1*grid_points[,1]+alpha2*grid_points[,2]+b>0 & grid_points[,1]>=Point[1] & grid_points[,2]>=Point[2])
                      Probability<-sum(dnorm(grid_points[idx,1])*dnorm(grid_points[idx,2])*epsilon^2)
                     },
                   #Q2
                   "2"={
                      t<-meshgrid(seq(-5,5,epsilon), seq(-5,5,epsilon))
                      x<-as.vector(t$X)
                      y<-as.vector(t$Y)  
                      grid_points<-cbind(x,y)
                      idx<-which(alpha1*grid_points[,1]+alpha2*grid_points[,2]+b>0 & grid_points[,1]>=Point[1] & grid_points[,2]<=Point[2])
                      Probability<-sum(dnorm(grid_points[idx,1])*dnorm(grid_points[idx,2])*epsilon^2)
                     },
                   
                   "3"={
                     t<-meshgrid(seq(-5,5,epsilon), seq(-5,5,epsilon))
                     x<-as.vector(t$X)
                     y<-as.vector(t$Y)  
                     grid_points<-cbind(x,y)
                     idx<-which(alpha1*grid_points[,1]+alpha2*grid_points[,2]+b>0 & grid_points[,1]<=Point[1] & grid_points[,2]<=Point[2])
                     Probability<-sum(dnorm(grid_points[idx,1])*dnorm(grid_points[idx,2])*epsilon^2)
                     },
                   "4"={
                      t<-meshgrid(seq(-5,5,epsilon), seq(-5,5,epsilon))
                      x<-as.vector(t$X)
                      y<-as.vector(t$Y)  
                      grid_points<-cbind(x,y)
                      idx<-which(alpha1*grid_points[,1]+alpha2*grid_points[,2]+b>0 & grid_points[,1]<=Point[1] & grid_points[,2]>=Point[2])
                      Probability<-sum(dnorm(grid_points[idx,1])*dnorm(grid_points[idx,2])*epsilon^2)
                      })})
  if(Probability==0)
  {
    Probability=0.000001
  }
  return(Probability/Z)       
}