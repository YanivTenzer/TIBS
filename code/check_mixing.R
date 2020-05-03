# Test permutations mixing
setwd("C:/Users/Or Zuk/Documents/GitHub/TIBS/code")  # Set to your directory 
source('TIBS.R')
library(PerMallows)

# Problem: power distribution doesn't look binomial 
set.seed(123456)
X <- rnorm(300)
Y <- rnorm(300)
trnc <- (X<Y)
x <- X[trnc]
y <- Y[trnc]
dat <- data.frame(x,y)
prep <- c()

prms = c()  # set parameters 
#prms$burn.in = 1000 # get default values 
#prms$Cycle = 400
prms$B = 300
prms$naive.expectation <- 0 # for 1 we get more stable statistic 

for (i in 1:50){  # Test 100 times the same data with a different permutation test randomization 
  print(c('Run', i, 'out of', 100))
  T <- TIBS(data=dat, w.fun="truncation", test.type='permutations',prms=prms)
  prep <- c(prep,T$Pvalue)
}


# Check permutations distances for last run 
#first.distance=first.distance, mixing.distance=mixing.distance, random.distance=random.distance)  # templ for debug
S <- TIBS(data=dat, w.fun="truncation", test.type='permutations',prms=prms) # sample again

mixing.distance <- rep(0, prms$B)  # diagnostics for MCMC mixing
first.distance <- rep(0, prms$B)
random.distance <- rep(0, prms$B)
tworuns.distance <- rep(0, prms$B)
for(i in 1:(prms$B-1)) # For Debugging! Calculate the distances between permutations
{
  #             print(c('dist ', i))
  if(mod(i,100)==0)
    print(c("Compute Distance Perm=", i))
  mixing.distance[i] <- distance(T$Permutations[,i], T$Permutations[,i+1])
  first.distance[i] <- distance(T$Permutations[,1], T$Permutations[,i+1])
  random.distance[i] <- distance(sample(1:dim(dat)[1]), sample(1:dim(dat)[1]))
  tworuns.distance[i] <- distance(T$Permutations[,i], S$Permutations[,i])
}

jpeg("../figs/perm_mixing.jpg", width = 350, height = 350)
plot(mixing.distance)
points(first.distance, col='red')
points(random.distance, col='green')
points(tworuns.distance, col='blue')
legend(0, 2000, legend=c('consecutive', 'from-first', 'random', 'two-runs'), col=c('black', 'red', 'green', 'blue'), 
       lty=c(1, 1, 1, 1), cex=0.8)
dev.off()


# New: plot the statistic: 
plot(T$statistics.under.null)
points(S$statistics.under.null, col='red')
lines(c(0,prms$B), c(T$TrueT, T$TrueT), col='green')


plot(sign(T$statistics.under.null-T$TrueT))
points(sign(S$statistics.under.null-T$TrueT), col='red')



# Consecutive differences vs. random differences
plot(T$statistics.under.null[-1] - head(T$statistics.under.null, -1))

diff <- matrix(c(T$statistics.under.null[-1] - head(T$statistics.under.null, -1), 
  sample(T$statistics.under.null[-1]) - head(T$statistics.under.null, -1), 
  S$statistics.under.null[-1] - head(S$statistics.under.null, -1), 
  sample(S$statistics.under.null[-1]) - head(S$statistics.under.null, -1)), 
  ncol=4)
boxplot(diff)



# Check if rejection counts for different simulations follow binomial distribution
p.sort <- sort(prep)
jpeg("../figs/fixed_grid_estimate_Pij.jpg", width = 350, height = 350)
ks.stat <- ks.test(p.sort*prms$B, "pbinom", prms$B, mean(prep))
plot(ecdf(prep), main=c('Binomial KS Pvalue:', round(ks.stat$p.value, 3)))
lines(p.sort,pbinom(p.sort*prms$B,prms$B,mean(prep)), col='blue')
summary(prep)
dev.off()


