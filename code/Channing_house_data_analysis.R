# TESTING INDEPENDENCE - CHANNING HOUSE DATA 
# Should unite this with real_life_datasets 

library(survival)
if(Sys.getenv('RSTUDIO') == '1')
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("TIBS.R")

#######################################################
#########################################################
# CHANNING HOUSE DATA
###########################################

### Get data.
data("channing",package = "boot")
x <- channing$entry
y <- channing$exit
cens <- channing$cens

# dewcriptive
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
# chacking risk set doesn't vanish before last obs
KM <- survfit(Surv(time=x,time2=y,event = delta) ~ 1,data=ch.data)
which(KM$n.risk==(KM$n.event+KM$n.censor))
length(KM$time)

# estimating censoring distribution and calculate W function
KM <- survfit(Surv(y-x,!delta) ~ 1,data=ch.data)
Srv.C <- stepfun(KM$time,c(1,KM$surv))
w.fun <- function(x,y){(x<y)*Srv.C(y-x)}

# filter the uncensored observations and apply TIBS
ch.delta1 <- data.frame(x=xx[delta==1],y=yy[delta==1])
test.res <- TIBS(data=ch.delta1, w.fun=w.fun, test.type='permutations',prms=list(B=100))



