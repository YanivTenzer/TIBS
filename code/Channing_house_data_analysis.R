# TESTING INDEPENDENCE - CHANNING HOUSE DATA 
# Should unite this with real_life_datasets 

library(survival)
library(permDep)
if(Sys.getenv('RSTUDIO') == '1')
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("TIBS.R")
source("Tsai_test.R")

#######################################################
#########################################################
# CHANNING HOUSE DATA
###########################################

### Get data.
data("channing",package = "boot")
x <- channing$entry
y <- channing$exit
cens <- channing$cens

# descriptive
plot(x,y,col=2-cens)
abline(0,1)
# correcting errors in data
ok <- which(y>x)
channing[which(y<=x),]
xx <- x[ok]
yy <- y[ok]
delta <- cens[ok]

### Create data frame.
ch.data <- data.frame(list("x"=xx,"y"=yy,"delta"=delta))
# chacking the risk set doesn't vanish before last obs
KM <- survfit(Surv(time=x,time2=y,event = delta) ~ 1,data=ch.data)
which(KM$n.risk==(KM$n.event+KM$n.censor))
length(KM$time)

# estimating censoring distribution and calculating W function
KM <- survfit(Surv(y-x,!delta) ~ 1,data=ch.data)
Srv.C <- stepfun(KM$time,c(1,KM$surv))
w.fun <- function(x,y){(x<y)*Srv.C(y-x)}

###################################
##   TESTS
###################################

set.seed(4820)
# filter the uncensored observations and apply TIBS
ch.delta1 <- data.frame(x=xx[delta==1],y=yy[delta==1])
for(i in 1:10){
tibs.res <- TIBS(data=ch.delta1, w.fun=w.fun, test.type='permutations',
                 prms=list(B=10000))
print(tibs.res$Pvalue)
}
# 0.8229

Tsai.test(xx,yy,delta)
# p         9.910721e-02

set.seed(4820)
print(permDep(trun = xx, obs = yy, permSize = 1000, cens = delta,
        sampling = "conditional",
        nc = 4, minp2Only = TRUE, kendallOnly = FALSE))
# Hypothesis test of quasi-independence based on minp2 statistics
# Conditional permutation size =  1000 
# p-value = 0.1409 

