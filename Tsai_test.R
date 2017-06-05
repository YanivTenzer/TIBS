# Tsai's test for quasi-independence
# test for truncation independent of total time

#ICU.dat <- read.table("C:\\Users\\mmandel\\Dropbox\\students\\Bella\\StatMed\\simple data.txt",header=TRUE)
#head(ICU.dat)
#ICU.dat[ICU.dat$delta==0,] # censoring is always at 30

# X - truncation time, Y - lifetime (no censoring)
# the criterion of truncation is X<=Y
tsai.test.ties <- function(X,Y){
  ii <- order(Y)
  X.sort <- X[ii]
  Y.sort <- Y[ii]
  n <- length(X)
  S <- NULL
  R <- NULL
  V <- NULL
  for (i in 1:n){
    risk.i <- which(X.sort[i:n]<=Y.sort[i])+i-1
    r.i <- length(risk.i)
    S.i <- sum(X.sort[risk.i]>X.sort[i])-sum(X.sort[risk.i]<X.sort[i])
    t.risk <- table(X.sort[risk.i])
    var.i <- (r.i^2-1)/3-sum(t.risk^3-t.risk)/(3*r.i)
    S[i] <- S.i
    R[i] <- r.i
    V[i] <- var.i
  }
  #browser()
  sumS <- sum(S)
  var.H0 <- sum(V[which(is.nan(V)==FALSE)])
  tsai.stat <- sumS^2/var.H0
  tsai.p <- 1-pchisq(tsai.stat,1)
  return(c(tsai.stat,tsai.p))
}

# in our data the truncation criterion is L<TU
# L=TU is not observed!
# thus, the risk sets should be adjusted for that
# therfore, we subtract 1 from TU, making the criterion L<=TU

#tsai.test.ties(ICU.dat$L,ICU.dat$TU-1)
