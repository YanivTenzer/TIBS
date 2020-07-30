library(Matrix)

##################################################################################################
# Compute the modified Hoeffding's test statistic, for the permutation test
# Parameters: 
# data - n*2 matrix with (x,y) sample
# grid.points - n*2 all possible (x_i,y_j) points  (can be the same as data, or a subset)
# null.expectations.table - n*4 mass table with pre-computed mass estimates for grid.points 
# 
#  Quardants convension:
#   4 | 1
#   ------
#   3 | 2
##################################################################################################
ComputeStatistic <- function(data, grid.points, null.expectations.table)
{
  obs.table <- matrix(0, dim(grid.points)[1], 4)
  Obs <- Exp <- matrix(0,4,1) # observed & expected
  
  #  print("DIM NULL TABLE")
  #  print(dim(null.expectations.table))
  Statistic <- 0 
  for (i in 1:dim(grid.points)[1])  # Slow loop on grid points 
  {
    Exp <- null.expectations.table[i,]
    Rx <- data[,1] > grid.points[i,1]
    Ry <- data[,2] > grid.points[i,2]
    Obs[1] <- sum(Rx*Ry)
    Obs[2] <- sum(Rx)-Obs[1]
    Obs[4] <- sum(Ry)-Obs[1]
    Obs[3] <- dim(data)[1]-sum(Obs[c(1,2,4)]) 
    obs.table[i,] <- Obs
#    print(Exp)
#    print(Obs)
    if (min(Exp)>1) {
#      print("Add To Expected")
      Statistic <-  Statistic + sum((Obs-Exp)^2 / Exp) # set valid statistic when expected is 0 or very small 
#      print(Statistic)
    }
  } # end loop on grid points 
  
  return(list(Statistic=Statistic, obs.table=obs.table)) # return also observed table for diagnostics
}

#########################################################################################
ComputeStatistic_inverse_weighting<- function(data, grid.points, w_mat)
{
  obs.table<-matrix(0, dim(grid.points)[1], 4)
  Obs<-Exp<-matrix(0,4,1) # observed & expected
  n_tilde <- sum(1/diag(w_mat))
  Statistic <- 0  
  min_Exp<-1/dim(data)[1]
  
  for (i in 1:dim(grid.points)[1])  # Slow loop on grid points 
  {
    Rx <- data[,1]>grid.points[i,1]
    Ry <- data[,2]>grid.points[i,2]
    
    idx1 <- which(Rx*Ry==1)
    Obs[1] <- sum(1/diag(w_mat)[idx1])
    idx2 <- which(Rx*(!Ry)==1)
    Obs[2] <- sum(1/diag(w_mat)[idx2])
    idx3 <- which((!Rx)*(!Ry)==1)
    Obs[3] <- sum(1/diag(w_mat)[idx3])
    idx4 <- which(Ry*(!Rx)==1)
    Obs[4] <- sum(1/diag(w_mat)[idx4])
    Obs <- Obs*(n_tilde^(-1))
    
    Exp[1] <- sum(1/diag(w_mat)[Rx])*sum(1/diag(w_mat)[Ry])
    Exp[2] <- sum(1/diag(w_mat)[Rx])*sum(1/diag(w_mat)[which(!Ry==1)])
    Exp[3] <- sum(1/diag(w_mat)[which(!Rx==1)])*sum(1/diag(w_mat)[which(!Ry==1)])
    Exp[4] <- sum(1/diag(w_mat)[which(!Rx==1)])*sum(1/diag(w_mat)[Ry])
    Exp <- Exp*(n_tilde^(-2))
    
    obs.table[i,] <- Obs
    
    if (min(Exp)>min_Exp){
      Statistic <-  Statistic + sum((Obs-Exp)^2 / Exp) # set valid statistic when expected is 0 or very small 
    }
  } # end loop on grid points 
  
  return(list(Statistic=Statistic, obs.table=obs.table)) # return also observed table for diagnostics
}

#########################################################################################
ComputeStatistic.W <- function(data, grid.points,w=function(x){1}){
  W <- apply(data,1,w)
  n.w <- sum(1/W)
  obs.table <- exp.table <- matrix(0, dim(grid.points)[1], 4)
  Obs <- Exp <- matrix(0,4,1) # observed & expected
  Statistic <- 0 
  for (i in 1:dim(grid.points)[1])  # Slow loop on grid points 
  {
    Rx <- data[,1]>grid.points[i,1]
    Ry <- data[,2]>grid.points[i,2]
    Exp[1] <- sum(Rx/W)*sum(Ry/W)/n.w^2
    Exp[2] <- sum(Rx/W)*sum((!Ry)/W)/n.w^2
    Exp[4] <- sum((!Rx)/W)*sum(Ry/W)/n.w^2
    Exp[3] <- sum((!Rx)/W)*sum((!Ry)/W)/n.w^2
    Obs[1] <- sum(Rx*Ry/W)/n.w
    Obs[2] <- sum(Rx*(!Ry)/W)/n.w
    Obs[4] <- sum((!Rx)*Ry/W)/n.w
    Obs[3] <- sum((!Rx)*(!Ry)/W)/n.w
#    print(Exp)
#    print(Obs)
    obs.table[i,] <- Obs
    exp.table[i,] <- Exp
    if (min(Exp)>(1/dim(data)[1])) {
      Statistic <-  Statistic + sum((Obs-Exp)^2 / Exp) # set valid statistic when expected is 0 or very small 
    } 
  } # end loop on grid points 
  
  return(list(Statistic=Statistic, obs.table=obs.table,exp.table=exp.table)) 
  # returns also expected and observed tables for diagnostics
}

  #################################################################
  # sample permutations, using MCMC, over the set of valid permutations, 
  # with respect to the distribution appears in Eq 8
  # Parameters: 
  # W - matrix with weights 
  # B - number of permutations to draw
  # N - sample size (can be read from data or W?)
  # 
  # Output: 
  # PermutationsTable - A matrix representing the sampled permutations 
  # P - An n*n matrix with P(i,j) = Pr(pi(i)=j) over the sampled permutations 
  #########################################################################################
  PermutationsMCMC <- function(w.mat, prms) # burn.in=NA, Cycle=NA)  # New: allow non-default burn-in 
  { 
    n <- dim(w.mat)[1];
    P <- matrix(0, n, n) # New! matrix with P[i]=j estimate
    #  for(i in 1:num.permutations) 
    #    P[cbind(1:n, Permutations[,i])] <- P[cbind(1:n, Permutations[,i])]+1 # need to vector indices here  
    
    
    
    # Set mcmc default sampling parameters 
    if(!('B' %in% names(prms)))
      prms$B <- 1000
    if(!('burn.in' %in% names(prms)))
      prms$burn.in <- 0  # 2*n  # set to n without burn.in, or 0 without burn.in to include the identity permutation 
    if(!('Cycle' %in% names(prms)))
      prms$Cycle <- n
    
    Idx <- ctr <- 1
    PermutationsTable = matrix(0, n, prms$B)
    Perm = 1:n # start with the identity 
    while(Idx<=prms$B)
    {
      
      # A Metropolis Hastings algorithm with target stationary distribution \pi
      # Choose the two indices to be switched
      switchIdx = sample(1:n, 2, replace = FALSE)  
      i = switchIdx[1]
      j = switchIdx[2]
      ratio = w.mat[i,Perm[j]]*w.mat[j,Perm[i]]/(w.mat[i,Perm[i]]*w.mat[j,Perm[j]]) 
      
      if(rbinom(1, 1, min(1,ratio))) #we accept the transition with probability min(1, ratio)
      {
        temp <- Perm[i] # SWAP 
        Perm[i] <- Perm[j]
        Perm[j] <- temp
        if(ctr==prms$burn.in || (ctr%%prms$Cycle==0 && ctr>prms$burn.in))
        {
          PermutationsTable[,Idx]=Perm;
          Idx = Idx+1;
#          if(mod(Idx,100)==0)
#            print(c("Sample Perm=", Idx))
        }
        P[cbind(1:n, Perm)] <- P[cbind(1:n, Perm)]+1 # Update table P
        ctr <- ctr+1;  # update counter only after swap 
      }
    }  # end while
    P <- P / (ctr-1) # normalize 
    
    return(list(PermutationsTable=PermutationsTable, P=P)) # New: return also P, a matrix with Pr(pi(i)=j)
  }

###################################################################################
# Estimate the null distribution fx*fy*W (given the estimated PDFs f_x, f_y)
# 
# Parameters: 
# pdfs - 2*n table with marginal distributions fx, fy probabilities 
# w - n*n matrix with weights W[i,j] = w(x_i, y_j) 
###################################################################################
GetNullDistribution <- function(pdfs, w_mat)
{
  # Compute the normalizing factor under the null:
  null.distribution <- w_mat * (pdfs[,1] %*% t(pdfs[,2]))
  Z <- sum(null.distribution)
  null.distribution <- null.distribution/Z
  
#  print("NULL DIM")
#  print(dim(null.distribution))
  return( list(distribution=null.distribution, Z=Z) )
}

############################################################################################
# Draw a bootstrap sample, given the estimated null distribution Fx, Fy, W
# Use rejection sampling. No need for computing n*n table 
# (Problem: what if prob.(rejection) close to 1?)
# Parameters: 
# data - n*2 array with (X,Y) samples
# pdfs - fx and fy  
# w.fun - Biased sampling function w 
# prms - for w max 
# n - allow a different sample size 
############################################################################################
Bootstrap <- function(data, pdfs, w.fun, prms, n=NULL)
{
  #  print("Inside Bootstrap")
  if(is.null(n))
    n = dim(data)[1]
#  print(dim(data))
  boot.sample <- matrix(-1, n, 2)
  k <- 0
  while(k<n) 
  {   # sampling n-k together
#       print("Inside Bootstrap sample x")
    x <- data[sample(dim(pdfs)[1], n-k, prob=pdfs[,1], replace=TRUE),1] # Sample X ~ Fx
#       print("Inside Bootstrap sample y")
    y <- data[sample(dim(pdfs)[1], n-k, prob=pdfs[,2], replace=TRUE),2] # Sample Y ~ Fy
#        print("Inside Bootstrap keep")
    keep <- which(as.logical(rbinom(n-k, 1, w_fun_eval(x, y, w.fun)/prms$W.max)))
#        print(keep)
    if(isempty(keep))
      next
    boot.sample[1:length(keep)+k,] <- cbind(x[keep],y[keep]) 
    k <- k+length(keep)
#     print(k)
  }    
  return(boot.sample)
}

###################################################################################################
# given a quadrant, evaluate the mass function within it
# New version - uses the two-dimensional CDF (not PDF)
# 
# Point - (x,y) value defining quadrants 
# QId - which quadrant 1,2,3,4
# data - 2*n array of (x,y)
# null.distribution.CDF - n*n table of 2d cumulative of Fx*Fy*w (problem with ties)
# Output: 
# S - expected probability at the quardant QId centered at Point for the distribution defined by null.distribution.CDF and data 
###################################################################################################
GetQuarterExpectedProb <- function(Point, QId, data, null.distribution.CDF)
{
  if(QId %in% c(1,2))
  {
    idx.x <- which(data[,1]>Point[1])
    idx.x <- idx.x[which.min(data[idx.x,1])]
  } else
  {
    idx.x <- which(data[,1]<=Point[1])
    idx.x <- idx.x[which.max(data[idx.x,1])]
  }
  if(QId %in% c(1,4))
  {
    idx.y <- which(data[,2]>Point[2])
    idx.y <- idx.y[which.min(data[idx.y,2])]
  } else
  {
    idx.y <- which(data[,2]<=Point[2])
    idx.y <- idx.y[which.max(data[idx.y,2])]
  }
  
  if(isempty(idx.x) | isempty(idx.y))
    return(0)    
  m <- which.max(data[,1])
  n <- which.max(data[,2])
  
  switch(QId, # First sample from Fxy
         {S <- null.distribution.CDF[m, n] + null.distribution.CDF[idx.x, idx.y] - 
           null.distribution.CDF[idx.x, n] - null.distribution.CDF[m, idx.y]}, # 1
         {S <- null.distribution.CDF[m, idx.y] - null.distribution.CDF[idx.x, idx.y]}, # 2
         {S <- null.distribution.CDF[idx.x, idx.y]}, # 3
         {S <- null.distribution.CDF[idx.x, n] - null.distribution.CDF[idx.x, idx.y]}) # 4
  return(S)       
}

###################################################################################################
# Compute Expect[Qi(p_j)] for 1<=i<=4, and all j, given a grid of points and bootstrap null distribution
# Parameters: 
# data - 2*n array (X,Y)
# null.distribution - a two-dim array n*n of probabilities under the null 
# Permutations - set of permutations
# grid.points - centers of partitions
#
# Output: 
# mass.table - a 4*#grid-points table with quadrants expectation
###################################################################################################
QuarterProbFromBootstrap <- function(data, null.distribution, grid.points)
{
  mass.table <- matrix(0, dim(grid.points)[1], 4)
  
#  print("NULL-DIST-DIM:")
#  print(dim(null.distribution))
  null.distribution.CDF <- PDFToCDF2d(null.distribution, data) 
  
  for(i in seq(1, dim(grid.points)[1],1))
  {
    for(j in 1:3)
      mass.table[i,j] <- GetQuarterExpectedProb(grid.points[i,], j, data, null.distribution.CDF)
    mass.table[i,4] = 1-sum( mass.table[i,1:3]) # , epsilon)
  }
  mass.table <- dim(data)[1]*mass.table # normalize to counts 
  
  return(mass.table)
}

###################################################################################################
# compute Expect[Qi(p_j)] for 1<=i<=4, and all j, given a grid of points using permutations
# Parameters: 
# data - 2*n array (X,Y)
# P - n*n table representing permutations probabilities Pr(pi(i)=j) (# Permutations - set of permutations)
# grid.points - centers of partitions 2*k array 
#
# Output: 
# mass.table - a 4*#grid-points table with quadrants expectation
###################################################################################################
QuarterProbFromPermutations <- function(data, P, grid.points) #Permutations
{
  #  num.permutations <- dim(Permutations)[2]
  #  n <- dim(data)[1]
  #  P <- matrix(0, n, n) # matrix with P[i]=j estimate
  #  for(i in 1:num.permutations) 
  #    P[cbind(1:n, Permutations[,i])] <- P[cbind(1:n, Permutations[,i])]+1 # need to vector indices here  
  #  P <- P / num.permutations # normalize 
  
  #next each point has its probability of being selected we evaluate Expect[Qi(p_j)] by summation
  mass.table<-matrix(0, dim(grid.points)[1], 4)
  for(i in seq(1, dim(grid.points)[1],1))
  {
    x<-grid.points[i,]
    mass.table[i,1] <- sum(P[data[,1]>x[1], data[,2]>x[2]])
    mass.table[i,2] <- sum(P[data[,1]>x[1], data[,2]<=x[2]])
    mass.table[i,3] <- sum(P[data[,1]<=x[1], data[,2]<=x[2]])
    mass.table[i,4] <- sum(P[data[,1]<=x[1], data[,2]>x[2]])
  } 
  #  mass.table<-mass.table+0.000001
  return(mass.table)
}

###################################################################################################
# Compute quarter probabilities for standard bivariate Gaussians 
# with truncation y>x. We assume grid.points satisfy y>x
#
# grid.points - where to evaluate probabilitites 
###################################################################################################
QuarterProbGaussianAnalytic <- function(grid.points)
{
  mass.table<-matrix(0, dim(grid.points)[1], 4)
  x <- grid.points[,1]
  y <- grid.points[,2]
  
  mass.table[,1]<-(1-pnorm(y))*(1+pnorm(y)-2*pnorm(x))
  mass.table[,2]<-(pnorm(x)-pnorm(y))^2
  mass.table[,3]<-pnorm(x)*(2*pnorm(y)-pnorm(x))
  mass.table[,4]<-2*pnorm(x)*(1-pnorm(y))
  return(mass.table)
}


##############################################################################
# Convert marginal PDF to PDF from data 
# Parameters: 
# CDF.table - matrix with every column a vector of CDF of each variable 
# Output: 
# PDF.table - matrix with every column a vector of PDF of each variable  
##############################################################################
CDFToPDFMarginals <- function(CDF.table)
{
  n<-dim(CDF.table)[1]  # number of samples 
  PDF.table <- array(0L, dim(CDF.table))  # matrix(0, num.samples, num.variables)
#  print("DIM CDF -> PDF:")
#  print(dim(PDF.table))
  for(i in 1:dim(CDF.table)[2])  # loop on variables 
  {
    sorted.CDF<-sort(CDF.table[,i], index.return=TRUE)
#    print("sorted.CDF:")
#    print(sorted.CDF)
    PDF.table[sorted.CDF$ix,i] <- c(sorted.CDF$x[1], sorted.CDF$x[-1]-sorted.CDF$x[-n])
  }
  return(PDF.table)
}



##############################################################################
# Convert marginal PDF to CDF from data 
# Parameters: 
# data - k*n array of variables 
# PDF.table - matrix with every column a vector of PDF of each variable  
# Output: 
# CDF.table - matrix with every column a vector of CDF of each variable 
##############################################################################
PDFToCDFMarginals <- function(data, PDF.table)
{
  n<-dim(PDF.table)[1]  # number of samples 
  CDF.table <- array(0L, dim(PDF.table))  # matrix(0, num.samples, num.variables)
#  print("DIM PDF -> CDF:")
#  print(dim(PDF.table))
  for(i in 1:dim(PDF.table)[2])  # loop on variables 
  {
    Px <- sort(data[,i], index.return=TRUE)  # Permute to order x_i, y_i 
#    print("Px=")
#    print(Px)
    CDF.table[Px$ix,i] <- cumsum(PDF.table[Px$ix,i])
#    print("idx:")
#    print(Px$ix)
#    print("cumsum:")
#    print(CDF.table[Px$ix,i])
        
  }
  return(CDF.table)
}


###################################################################################################
# Compute 2d cumulative distribution. When we have ties we need to correct this
# Input: 
# pdf.2d - a two-dim array n*n of probabilities 
# data - xy points with probabilities, array size n*2 (used for sorting)
# Output: 
# cdf.2d - a two-dim array n*n of cumulative probabilities 
###################################################################################################
PDFToCDF2d <- function(pdf.2d, data)
{
#  print("DIM PDF.2D:")
#  print(dim(pdf.2d))
#  print("DIM DATA:")
#  print(dim(data))
  Px <- sort(data[,1], index.return=TRUE)  # Permute to order x_i, y_i 
  Py <- sort(data[,2], index.return=TRUE)
  cdf.2d <- apply(apply(pdf.2d[Px$ix, Py$ix], 1, cumsum), 1, cumsum)  # cumsum on rows and columns 
  
  # Use data to deal with ties (not working yet)
  #  ties.x <- which(data[-1,1] == head(data[,1], -1))
  #  ties.y <- which(data[-1,2] == head(data[,2], -1))
  #  for(i in rev(ties.y))
  #    cdf.2d[,i] <- cdf.2d[,i+1]
  #  for(i in rev(ties.x))
  #    cdf.2d[i,] <- cdf.2d[i+1,]

  # This part is equivalent to the returned matrix 
#  n = length(Px$ix)  
#  cdf.3d <- matrix(0, nrow=n, ncol=n)
#  cdf.3d[Px$ix, Py$ix] = cdf.2d
#  return (cdf.3d)
#  return(cdf.2d)
  return( cdf.2d[invPerm(Px$ix), invPerm(Py$ix)] )  # why Py first and then Px? 
}

###################################################################################################
# The true marginals (same for x and y) for uniform distribution on strip of width a
# t - where to evaluate CDF
# a - width of strip 
###################################################################################################
UniformStripMarginalCDF <- function(t,a)
{
  if (a>0.5 | a<0) stop("a must be in (0,0.5]")
  c <- a*(2-a)
  val <- 3*a^2/(2*c) + 2*a*(1-2*a)/c + (1+3*a-t)*(t+a-1)/(2*c) # all indices 
  val[t<a] <- t[t<a]*(t[t<a]+2*a)/(2*c) # for t<a
  val[(t<1-a) & (t>=a)] <- 3*a^2/(2*c) + 2*a*(t[(t<1-a) & (t>=a)]-a)/c # others with t<1-a
  return(val)
}


###################################################################################################
# Plot scatter and estimated mariginals 
# dependence.type - dependency label
# biased.data - 2*n array with sampled data
# prms - structure with parameters 
###################################################################################################
PlotBiasedData <- function(dependence.type, biased.data, prms)
{
  n <- dim(biased.data)[1]
  output.scatter.file <- paste0(path, '/../Figures/simulations/', dependence.type, '_rho_', prms$rho)  # set output figure file name
  for(plot.type in c('scatter', 'KM_marginal', 'marginal_scatter'))
  {
    switch(plot.type, # First sample from Fxy
           'scatter' ={
             xy <- data.frame(biased.data)   # plot and save to image file 
           },
           'KM_marginal' ={
             marginals <- EstimateMarginals(biased.data, 'survival')
             marginals.CDFs <- data.frame(cbind(marginals$xy, marginals$CDFs, 
                                                seq(0, 1, len=n), UniformStripMarginalCDF(seq(0, 1, len=n), prms$rho)))   # estimated marginals 
           },
           'marginal_scatter' ={
             Px <- c(marginals$CDFs[-1,1], 1)-marginals$CDFs[,1]
             Py <- marginals$CDFs[,2]-c(0, marginals$CDFs[-n,2])
             x.ind <- sample(x = marginals$xy[,1], size = 10000, replace = TRUE, prob = Px)
             y.ind <- sample(x = marginals$xy[,2], size = 10000, replace = TRUE, prob = Py)
             w.ind <- which(x.ind < y.ind) 
             xy <- data.frame(cbind(x.ind[w.ind], y.ind[w.ind]))   # sample from estiamted marginals
           }, 
           'all_scatter' = { # here plot both biased and original samples 
             xy <- data.frame(biased.data)             
           }
    ) # end switch 
    if(plot.type=='KM_marginal')
      print( ggplot(marginals.CDFs, aes(x=marginals.CDFs[,1], y=marginals.CDFs[,3], color= )) + 
               geom_line(aes(x=marginals.CDFs[,1], y=marginals.CDFs[,3], col="\\hat{F}_X"), size=1.5) + 
               geom_line(aes(x=marginals.CDFs[,2], y=marginals.CDFs[,4], col="\\hat{F}_y"), size=1.5) + 
               geom_line(aes(x=marginals.CDFs[,5], y=marginals.CDFs[,6], col='F_X===F_Y'), size=1.5) + # change line width
               ggtitle(TeX(rep(paste0("$", gsub("_", '-', dependence.type), " (\\theta = ", prms$rho, " )$"), prms$title))) + 
               xlab("X") + ylab("Y") + # keep defaults 
               theme(plot.title = element_text(size=14, face="bold.italic", hjust=0.5),
                     axis.title.y = element_text(face="bold", size=14),
                     axis.title.x = element_text(face="bold", size = 14),
                     axis.text.x = element_text(face="bold", size=12), 
                     axis.text.y = element_text(face="bold", size=12), 
                     legend.background = element_rect(fill = alpha("lightgray", 0)),
                     legend.title = element_text(), 
                     legend.position =  c(0.85, 0.2)) + 
               labs(color = "") + 
               scale_color_discrete(labels=lapply(sprintf(c("\\hat{F}_X", "\\hat{F}_Y", "F_X=F_Y")), TeX)) )
    else # new: plot two scatters on same datasets 
      print( ggplot(xy, aes(x=xy[,1], y=xy[,2])) + 
               #               geom_point(aes(x=xy[,1], y=xy[,2], col="original")) + 
               #               geom_point(shape=3, aes(x=xy[,3], y=xy[,4], col="biased")) + 
               geom_point() + 
               ggtitle(TeX(rep(paste0( gsub("_", '-', dependence.type), " $(\\theta = ", prms$rho, " )$"), prms$title))) + 
               xlab("X") + ylab("Y") +
               theme(plot.title = element_text(size=14, face="bold.italic", hjust=0.5),
                     axis.title.y = element_text(face="bold", size=14),
                     axis.title.x = element_text(face="bold", size = 14),
                     axis.text.x = element_text(face="bold", size=12), 
                     axis.text.y = element_text(face="bold", size=12), 
                     legend.position = "none") )
    
    ggsave(paste0(output.scatter.file, '_', plot.type, '.png'), units="in", width=5, height=5, dpi=300)
    if(dependence.type != 'UniformStrip') # last two plots are only for first dataset  
      break
  } # end loop on plot type
}


###################################################################################################
# Compute density of product of two Gaussian densities  
# Parameters:
# mu1, mu2 - mean vectors 
# sigma1, sigma2 - covariance matrices
# 
# Output:
# mu - mean vector of product distribution
# sigma - covariance matrix of product distribution
###################################################################################################
GaussianDensityProduct <- function(mu1, mu2, sigma1, sigma2)
{
  sigma <- solve(solve(sigma1) + solve(sigma2))
  mu <- sigma * (solve(sigma1) * mu1 + solve(sigma2) * mu2)
  return(list(mu=mu, sigma=sigma))
}  



# Read real dataset 
ReadDataset <- function(data_str)
{
  switch(data_str,  # read dataset 
         'huji'={
           load('../data/APage_APtime.RData')
           input.data <- HUJI.dat
         },
         'AIDS'={
           library(DTDA)
           data("AIDS")
           input.data <- AIDS[,2:3]
         },
         'ICU'={  # this dataset is not part of the released package
           input.data <- read.table("../data/ICU_data.txt", header = TRUE)
           input.data <- input.data[,c(1:2)]  # discard third column 
           input.data[,2] <- input.data[,2]-0.02 # correction: reduce by 1: the truncation criterion is L<TU (L==TU is removed)
         }, 
         'Infection'={ # this dataset is not part of the released package
           load('../data/ICU_INF.Rdata')
           input.data <- cbind(X,Y)
           W.max <- max(X)+max(Y)
         }, 
         'Dementia'={ # this dataset is not part of the released package
           # DEMENTIA DATA
           
           if(first_time) ### Get data.
           {
             prevdata<-read.csv("C:/Users/mm/Dropbox/marco/Nonparametric bivariate estimation/Old folder/Data analysis/cshaforRCSV.csv"); # change to relative path 
             ### Clean data.
             badindex<-which(prevdata$duration-prevdata$truncation<=0);
             prevdata<-prevdata[-badindex,];
             
             ### Order data according to disease duration.
             x<-prevdata$AAO;
             v<-prevdata$duration;
             w<-prevdata$AAO+prevdata$truncation;
             delta<-prevdata$death;
             order.v<-order(v);
             x<-x[order.v];
             v<-v[order.v];
             w<-w[order.v];
             delta<-delta[order.v];
             
             ### Create data frame.
             cshadata <- data.frame(list("x"=x,"v"=v,"w"=w,"delta"=delta))  # cshadata
             
             # the risk set vanishes before the last observation. 
             # remove the largest 3 v's to take care of that.
             id.rs <- which(cshadata$v>25)
             cshadata1 <- cshadata[-id.rs,] 
             
             # before ommiting the largest 3
             KM <- survfit(Surv(v-(w-x),!delta) ~ 1,data=cshadata)
             Srv.C1 <- stepfun(KM$time,c(1,exp(-KM$cumhaz)))
             Srv.C2 <- stepfun(KM$time,c(1,KM$surv))
             w.fun <- function(x,y){(x<y)*Srv.C1(y-x)}
             w.fun2 <- function(x,y){(x<y)*Srv.C2(y-x)}
             x.csha <- cshadata$w-cshadata$x
             y.csha <- cshadata$v
             delta.csha <- cshadata$delta
             input.data <- data.frame(x.csha[delta.csha],y.csha[delta.csha])  # csha.delta1
             
             # Save data frame (next time we can load only this)
             save(input.data, KM, file='../data/Dementia.Rdata')          
           }                
           else
           {
             load('../data/Dementia.Rdata')
             Srv.C1 <- stepfun(KM$time,c(1,exp(-KM$cumhaz)))
             w.fun[d] <- function(x,y){(x<y)*Srv.C1(y-x)}  # modify w.fun 
           }
           #           TIBS(data=csha.delta1, w.fun=w.fun1, B=1000, test.type='permutations',prms=c())
         }
  ) # end switch 

  if(!is.numeric(input.data))   # unlist and keep dimensions for data 
    input.data <- array(as.numeric(unlist(input.data)), dim(input.data))  
  
  return(input.data) 
}
