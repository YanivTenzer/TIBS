library(copula)
library(ggplot2)
###################################################################
# scatter plot of monotone-exchangeable example, Gumbel(\theta=1.6) 
###################################################################
n=400
prms<-list(rho=1.6)
gumbel.cop <- gumbelCopula(prms$rho)#(1.2)
ranks<- rCopula(n, gumbel.cop)
xy = data.frame(cbind(qnorm(ranks[,1]), qnorm(ranks[,2])))

ggplot(xy, aes(x=xy[,1], y=xy[,2])) + 
geom_point() + 
ggtitle(expression(paste("GC(", theta, ")=1.6"))) +
xlab("X") + ylab("Y") +
scale_y_continuous(breaks=seq(-2, 2, 2))+
scale_x_continuous(breaks=seq(-2, 2, 2))+
theme(
  plot.title = element_text(size=14, face="bold.italic", hjust=0.5),
  axis.title.y = element_text(face="bold", size=14),
  axis.title.x = element_text(face="bold", size = 14),
  axis.text.x = element_text(face="bold", size=12), 
  axis.text.y = element_text(face="bold", size=12),
)



