source("simulate_biased_sample.R")
source("TIBS.R")
library(tidyverse)
library(copula)
library(ggplot2)
library(tikzDevice)
library(latex2exp)
library(dplyr)
library(mvtnorm)
Rcpp::sourceCpp("C/utilities_ToR.cpp")  # all functions are here 



# Which one changed: LD (5), nonmonotone_nonexxhagable (6)  (0.9, 0.7, ... -0.9), , CLmix (7) - should re-run them 


# Same parameters of simulation, but here only subset that we want to plot
run.params.mat <- t(matrix(c('UniformStrip', 'truncation', TRUE, TRUE, list(0.3), 
                             'Gaussian', 'truncation', TRUE, TRUE, list(c(-0.9,-0.5,0.5,0.9)),
                             'Clayton','truncation', TRUE, TRUE, list(0.5),
                             'Gumbel', 'truncation', TRUE, TRUE, list(1.6),
                             'LD', 'truncation', TRUE, FALSE, list(c(0, 0.4)),
                             'nonmonotone_nonexchangeable', 'truncation', FALSE, FALSE, list(c(-0.9,-0.5,0,0.5,0.9)), # list(c(-0.9, -0.5, 0.0, 0.5, 0.9)),
                             'CLmix','truncation', FALSE, TRUE, list(0.5), 
                             'LogNormal', 'sum', TRUE, TRUE,  list(c(0, 0.2, 0.5, 0.9)),  # added also a signal 
                             'Gaussian', 'gaussian', TRUE, TRUE, list(seq(0.0, 0.9, 0.1)) ), 5, 9)) # no need for negatives # -0.9 - 0.9 replace by CLmix / non-monotone and centered at zero 

dependence.type <- run.params.mat[,1]
w.fun <- run.params.mat[,2]
monotone.type <- run.params.mat[,3]
exchange.type <- run.params.mat[,4]
prms.rho <- run.params.mat[,5]

run.dep <- c(9) # (1:8)

to.sim = "LogNormal"
###################################################################
# scatter plot of monotone-exchangeable example, Gumbel(\theta=1.6) 
###################################################################
n=400

prms <- c()
prms$keep.all=1
prms$use.cpp=1

all.cor.vec <- rep(0, length(prms.rho[[d]]))
for(d in run.dep)
{
  
  for(r in c(1:length(prms.rho[[d]])))
  {
    prms$rho <- prms.rho[[d]][r]  
    prms$w.rho <- -prms$rho
    
    
    x.perm <- c()
    y.perm <- c()
    for(i in c(1:50))
    {
      dat <- SimulateBiasedSample(n, dependence.type[[d]], w.fun[[d]], prms)

      prms$B = 5
      cor.vec <- rep(0, prms$B)
      T <- TIBS(dat$data, w.fun[[d]], prms, "permutationsMCMC", "adjusted_w_hoeffding")  # New: generate a permuted dataset 
      for(b in 1:prms$B)
      {
        cor.vec[b] = cor(dat$data[,1], dat$data[T$Permutations[,b],2])
        x.perm <- c(x.perm, dat$data[,1])
        y.perm <- c(y.perm, dat$data[T$Permutations[,b],2])
      }
    }
    
    df.perm <- data.frame(x=x.perm, y=y.perm)
          print(r)
    print(mean(cor.vec))
    all.cor.vec[r] = mean(cor.vec)
    
    ggplot(df.perm, aes(x=x, y=y) ) +
      stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
#      geom_density_2d() + # geom_hex() +
      theme_bw()
    
        
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
    if(d==6)  # Cut tail of Weibull distribution 
    {
      x.lim[2] <- 40
    }


        
    ggplot(xy.all, aes(x=xy.all[,1], y=xy.all[,2])) + 
      stat_density_2d(data=df.perm, aes(x=df.perm[,1], y=df.perm[,2], fill = ..level..), geom = "polygon", colour="white", alpha=0.7) +
      geom_point(data=xy.good, aes(x=xy.good[,1], y=xy.good[,2]), colour="red", size=1) + 
      geom_point(data=xy.bad, aes(x=xy.bad[,1], y=xy.bad[,2]), colour="gray40", size=0.5) + 
      geom_line(aes(x=S1, y=S2), colour="black", alpha=trans) +
#      ggtitle(title.str) +
      xlab("X") + ylab("Y") +
      xlim(c(-3,3)) + ylim(c(-3,3)) +  #      xlim(x.lim) + ylim(y.lim) + 
      theme(
        plot.title = element_text(size=14, face="italic", hjust=0.5),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size = 14),
        axis.text.x = element_text(face="bold", size=12), 
        axis.text.y = element_text(face="bold", size=12),
      )
    ggsave(paste('../figs/perm_dist_', dependence.type[[d]], '_', w.fun[[d]], '_', prms$rho, '.jpg', sep=''))
  }  
} # loop on datasets 


plot(xy.all[,1], xy.all[,2], col="red")
points(xy[,1], xy[,2])
points(T$permuted.data[,1], T$permuted.data[,2], co="green")
abline(0,1)

plot(xy[,1], xy[,2])
abline(0,1)


plot(xy[,1], xy[,2])
points(T$permuted.data[,1], T$permuted.data[,2], co="green")
cor(xy[,1], xy[,2])
cor(T$permuted.data[,1], T$permuted.data[,2])
abline(0,1)



#r.vec <- prms$rho / (1+sqrt(1+prms$rho^2))

rho.vec <- seq(0,100,1)/100
r.vec <-  -rho.vec / (1+sqrt(1+rho.vec^2)) #  2*(-1+sqrt(5-4*rho.vec^2)) / ( 5 -4*rho.vec^2 + sqrt(5-4*rho.vec^2) )
sigma.x.vec <-  1 / (1+sqrt(1+rho.vec^2))
plot(rho.vec, r.vec)
points(prms.rho[[d]], all.cor.vec, col="red")

sigma2.vec2 <- (-rho.vec^2 + sqrt(rho.vec^2+4*(1-rho.vec^2)^2)) / 2
r.vec2 <- rho.vec * sigma2.vec2 / (sigma2.vec2 + 1 - rho.vec^2)
points(rho.vec, r.vec2, col="green")

#plot(rho.vec, sigma.x.vec)


prms$rho / (1+sqrt(1+prms$rho^2))


# New calculation: 
rho=0.5
(1-2*rho^2)^2 - 4*(rho^2-1)*(4*rho^2+2)
new.sigma2.vec <- (2*rho.vec^2-1 + sqrt( (2*rho.vec^2-1)^2 - 4 * (rho.vec^2-1)*(4*rho.vec^2+2)  )) / (4*(2*rho.vec^2+1))
new.r.vec <- -rho.vec*new.sigma2.vec / (1-rho.vec^2+new.sigma2.vec)
plot(rho.vec, new.sigma2.vec, ylim=c(-1,1))
points(rho.vec, new.r.vec, col="green")

new.var.vec.from.sigma2.vec <- new.sigma2.vec*(1-rho.vec^2)*(1-rho.vec^2+new.sigma2.vec) / 
  ((1-rho.vec^2+new.sigma2.vec)-4*rho.vec^2*new.sigma2.vec^2)

plot((1-rho.vec^2)/2 -  new.var.vec.from.sigma2.vec)

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



if(plot_gaussian_example)
{
  rho <- 0.9
  df <- data.frame(x = seq(-3, 3,0.02), y = seq(-3, 3,0.02))
#  y = x
  n <- length(df$x)
  

  Sigma.perm <- (1-rho^2)/2 * matrix(c(1, -rho/(1+sqrt(1+rho^2)), -rho/(1+sqrt(1+rho^2)), 1), nrow=2, ncol=2) 
#    Sigma.perm <- matrix(c(1, -rho/(1+sqrt(1+rho^2)), -rho/(1+sqrt(1+rho^2)), 1), nrow=2, ncol=2) / (1+sqrt(1+rho^2))
  Sigma.f <- matrix(c(1, rho, rho, 1), nrow=2, ncol=2)
  Sigma.f.w <- matrix(c((1-rho^2)/2, 0, 0, (1-rho^2)/2), nrow=2, ncol=2)
  z.f <- matrix(0, n, n)
  z.f.w <-  matrix(0, n, n)
  z.perm <- matrix(0, n, n) 
  for(i in c(1:n))
    for(j in c(1:n))
    {
      z.f[i,j] = dmvnorm(c(df$x[i],df$y[j]), rep(0,2), Sigma.f)  #  log(y+1) + log(y+1)/log(x+1))
      z.f.w[i,j] = dmvnorm(c(df$x[i],df$y[j]), rep(0,2), Sigma.f.w)  #  log(y+1) + log(y+1)/log(x+1))
      z.perm[i,j] = dmvnorm(c(df$x[i],df$y[j]), rep(0,2), Sigma.perm)  #  log(y+1) + log(y+1)/log(x+1))
    }
  
  L <- 1.2
  contour(df$x, df$y, z.f, xlab="x", ylab="y", xlim=c(-L,L), ylim=c(-L,L))
  contour(df$x, df$y, z.f.w, col="red", add=TRUE)
  contour(df$x, df$y, z.perm, col="green", add=TRUE)
  
  colnames(z) = y
  rownames(z) = x
  
  # Convert to long format, convert row and column values to numeric, and create groups for colors
  df = as.data.frame(z) %>% 
    rownames_to_column(var="x") %>% 
    gather(y, value, -x) %>% 
    mutate(y=as.numeric(y), 
           x=as.numeric(x),
           value_range = cut(value, 8))
  
  
  ggplot(df, aes(x, y, fill=value_range)) + 
    geom_raster() +
    scale_fill_manual(values=colorRampPalette(c("red","orange","yellow"))(8)) +
    theme_classic() +
    guides(fill=guide_legend(reverse=TRUE))
  library(plotly)
  
  plot_ly(df, x = ~x, y = ~y, z = ~z) %>% 
    add_surface()
  
  
}  
  

  
  
  
  
  
  