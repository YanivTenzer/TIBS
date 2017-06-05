Get_Quarter_Approximated_Mass<-function(Point, QId, data, null_distribution)
{
  #epsilon = 0.00001
  switch(QId,
  #Q1
  "1"={
        S=0
        idx_x<-which(data[,1]>=Point[1])
        idx_y<-which(data[,2]>=Point[2])
        if(length(idx_y)>0 && length(idx_x)>0)
        {  
          for(i in 1:length(idx_x))
          {
            S=S+sum(null_distribution[idx_x[i], idx_y])  
          } 
        }
        
        },
                   
    #Q2
    "2"={
          S=0
          idx_x<-which(data[,1]>=Point[1])
          idx_y<-which(data[,2]<Point[2])
          if(length(idx_y)>0)
          {
            for(i in 1:length(idx_x)&& length(idx_x)>0)
            {
              S=S+sum(null_distribution[idx_x[i], idx_y])  
            }
          }
        },
      #Q3
      "3"={
          S=0
          idx_x<-which(data[,1]<Point[1])
          idx_y<-which(data[,2]<Point[2])
          if(length(idx_y)>0 && length(idx_x)>0)
          {
            for(i in 1:length(idx_x))
            {
              S=S+sum(null_distribution[idx_x[i], idx_y])  
            
            }
          }
        },
      
      #Q4
      "4"={
           S=0
           idx_x<-which(data[,1]<Point[1])
           idx_y<-which(data[,2]>=Point[2])
           if(length(idx_y)>0 && length(idx_x)>0)
           {
            for(i in 1:length(idx_x))
            {
                if(length(idx_y)>0)
                {
                  S=S+sum(null_distribution[idx_x[i], idx_y])   
                } 
            }
            }
          }         #End_of_switch_QId
        )
    
  
  #S=max(epsilon,S)
  return(S)       
}