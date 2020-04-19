# Tsai's test for quasi-independence [Tsai, 1990].
# P-val by normal approximation - OK censoring and ties. 
# L - truncation time
# Y - lifetime (no censoring)
# del - optional censoring indicator (1=failure, 0=censoring)
# NOTE: The criterion of truncation is L<=Y so:
# Y_j is in risk set of Y_i if L_j=Y_i
# For strict truncation L<Y add 'epsilon' to all Y's

Tsai.test <- function(L,Y,del=rep(1,length(Y))){
  n <- length(Y)
  Risk <- vector(mode = "list", length = n)
  for (i in 1:n) {
    set.i <- which(L<=Y[i] & Y>=Y[i])
    Risk[[i]] <- set.i
  }
  S.i <- NULL
  V.i <- NULL
  for (i in 1:n){
    if (del[i]==1) {
      S.i[i] <- sum(sign(L[Risk[[i]]]-L[i]))
      ties <- table(L[Risk[[i]]])
      risk.size <- length(Risk[[i]])
      V.i[i] <- (risk.size^2-1)/3-sum(ties^3-ties)/(3*risk.size)
    } else {
      S.i[i] <- 0
      V.i[i] <- 0
    }
  }
  K.star <- sum(S.i)
  V <- sum(V.i)
  test.stat <- K.star/sqrt(V)
  p <- 2*pnorm(-abs(test.stat))
  return(rbind(K.star,sqrt(V),test.stat,p))
}

