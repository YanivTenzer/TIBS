# TESTING INDEPENDENCE - DEMENTIA DATA 

library(survival)
source("TIBS.R")

#######################################################
#########################################################
# DEMENTIA DATA
###########################################

### Get data.
prevdata<-read.csv("C:/Users/mm/Dropbox/marco/Nonparametric bivariate estimation/Old folder/Data analysis/cshaforRCSV.csv");
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
cshadata <- data.frame(list("x"=x,"v"=v,"w"=w,"delta"=delta))

# the risk set vanishes before the last observation. 
# remove the largest 3 v's to take care of that.
id.rs <- which(cshadata$v>25)
cshadata1 <- cshadata[-id.rs,] 

# before ommiting the largest 3
KM <- survfit(Surv(v-(w-x),!delta) ~ 1,data=cshadata)
Srv.C1 <- stepfun(KM$time,c(1,exp(-KM$cumhaz)))
Srv.C2 <- stepfun(KM$time,c(1,KM$surv))
w.fun1 <- function(x,y){(x<y)*Srv.C1(y-x)}
w.fun2 <- function(x,y){(x<y)*Srv.C2(y-x)}
x.csha <- cshadata$w-cshadata$x
y.csha <- cshadata$v
delta.csha <- cshadata$delta
csha.delta1 <- data.frame(x.csha[delta.csha],y.csha[delta.csha])
TIBS(data=csha.delta1, bias.type=w.fun1, B=1000, test.type='permutations',prms=c())



