source("simulate_biased_sample.R")
library(copula)
library(ggplot2)
library(tikzDevice)
library(latex2exp)
library(dplyr)



# Same parameters of simulation, but here only subset that we want to plot
run.params.mat <- t(matrix(c('UniformStrip', 'truncation', TRUE, TRUE, list(0.3), 
                             'Gaussian', 'truncation', TRUE, TRUE, list(c(-0.9,-0.5,0.5,0.9)),
                             'Clayton','truncation', TRUE, TRUE, list(0.5),
                             'Gumbel', 'truncation', TRUE, TRUE, list(1.6),
                             'LD', 'truncation', TRUE, FALSE, list(c(0, 0.4)),
                             'nonmonotone_nonexchangeable', 'truncation', FALSE, FALSE, list(c(-0.9, -0.5, 0.0, 0.5, 0.9)),
                             'CLmix','truncation', FALSE, TRUE, list(0.5), 
                             'LogNormal', 'sum', TRUE, TRUE,  list(c(0, 0.2, 0.5, 0.9)),  # added also a signal 
                             'Gaussian', 'gaussian', TRUE, TRUE, list(c(0.0, 0.5, 0.9)) ), 5, 9)) # no need for negatives # -0.9 - 0.9 replace by CLmix / non-monotone and centered at zero 

dependence.type <- run.params.mat[,1]
w.fun <- run.params.mat[,2]
monotone.type <- run.params.mat[,3]
exchange.type <- run.params.mat[,4]
prms.rho <- run.params.mat[,5]

run.dep <- c(2:9) # (1:8)

to.sim = "LogNormal"
###################################################################
# scatter plot of monotone-exchangeable example, Gumbel(\theta=1.6) 
###################################################################
n=400

prms$keep.all=1

for(d in run.dep)
{
  
  for(r in c(1:length(prms.rho[[d]])))
  {
    prms$rho <- prms.rho[[d]][r]  
    prms$w.rho <- -prms$rho
    
    dat <- SimulateBiasedSample(n, dependence.type[[d]], w.fun[[d]], prms)
    xy <- as.data.frame(dat$data)
    xy.all <- as.data.frame(dat$all.data[1:n,])
    r <- range(xy.all)
    # New add another column to enable a line
    xy.all$S1 <- seq(r[1], r[2], (r[2]-r[1])/(n-1))
    xy.all$S2 <- seq(r[1], r[2], (r[2]-r[1])/(n-1))  
    title.str = c()  # expression(italic("LogNormal(0.5)"))
    
    if(w.fun[[d]] == "truncation")  # for truncation 
    {
      bad.points <- which(xy.all[,1] > xy.all[,2])
      good.points <- setdiff(c(1:n), bad.points)
      xy.good <- as.data.frame(xy.all[good.points,]) #  %>% filter(xy.all$V1 <= xy.all$V2))
      xy.bad <- as.data.frame(xy.all[bad.points,])
      trans <- 1
      x.lim <- range(xy.all[,1])
      y.lim <- range(xy.all[,2])
    } else
    {
      xy.good <- xy  # here show all points 
      xy.bad <- xy.all
      trans <- 0 
      x.lim <- c(min(xy.all[,1], xy[,1]), max(xy.all[,1], xy[,1]))
      y.lim <- c(min(xy.all[,2], xy[,2]), max(xy.all[,2], xy[,2]))
    }
    
    ggplot(xy.all, aes(x=xy.all[,1], y=xy.all[,2])) + 
      geom_point(data=xy.good, aes(x=xy.good[,1], y=xy.good[,2]), colour="red", size=2) + 
      geom_point(data=xy.bad, aes(x=xy.bad[,1], y=xy.bad[,2]), colour="gray40", size=1) + 
      geom_line(aes(x=S1, y=S2), colour="black", alpha=trans) +
      ggtitle(title.str) +
      xlab("X") + ylab("Y") +
      xlim(x.lim) + ylim(y.lim) + 
      theme(
        plot.title = element_text(size=14, face="italic", hjust=0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size = 14),
        axis.text.x = element_text(face="bold", size=12), 
        axis.text.y = element_text(face="bold", size=12),
      )
    ggsave(paste('../figs/', dependence.type[[d]], '_', w.fun[[d]], '_', prms$rho, '.jpg', sep=''))
  }  
} # loop on datasets 


#plot(xy[,1], xy[,2])
#points(xy.all[,1], xy.all[,2], col="red")




#  if(to.sim == "Gumbel")
#  {
#    prms<-list(rho=1.6)
#    gumbel.cop <- gumbelCopula(prms$rho)#(1.2)
#    ranks<- rCopula(n, gumbel.cop)
#    xy = data.frame(cbind(qnorm(ranks[,1]), qnorm(ranks[,2])))
#    title.str = expression(paste("GC(", theta, ")=1.6"))
#  } else 
#  {
#    # New: show also scatter plot of positive biased sampling with log normal and sum 
#    
#    prms = c()
#    prms$rho=0.0
#    prms$keep.all=1
#    dat <- SimulateBiasedSample(n, "LogNormal", "sum", prms)
#    xy <- as.data.frame(dat$data)
#    xy.all <- as.data.frame(dat$all.data[1:n,])
