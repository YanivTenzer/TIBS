# Set path to your working directory: 
# The script contains real data analysis 
# Due to privacy issues, the data files for the ICU datasets are not part of the released package. 
# Please contact the authors if you are interested in them. 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('TIBS.R')
source('simulate_biased_sample.R')
source('marginal_estimation.R')
source('Tsai_test.R')
source('utilities.R') 
library(lubridate)
library(permDep)
library(tictoc)
library(xtable)
library(ggplot2)
library(Matrix)

B = 100 # 10000 # number of permutations/bootstrap samples 
plot.flag <- 1

test.type <- c('bootstrap', 'permutations', 'tsai', 'minP2') # different tests to run 
datasets <- c('huji', 'AIDS', 'ICU', 'Infection')
exchange.type <- c(FALSE, FALSE, FALSE, FALSE) # no reason to assume real data is exchangeable
w.fun <- c('huji', 'Hyperplane_Truncation', 'Hyperplane_Truncation', 'sum')
prms = c()
W.max <- c(65, 1, 1, -1) # -1 denotes calculate max of w from data  
n.datasets <- length(datasets)
n.tests <- length(test.type)
Hyperplane.prms<-c(-1,1,0)
test.pvalue <- matrix(-1, n.datasets, n.tests) # -1 denotes test wasn't performed
test.time <- matrix(-1, n.datasets, n.tests) # -1 denotes test wasn't performed

for(d in 1:n.datasets) # loop on datasets.
{
#    if(datasets[d] %in% c('ICU', 'Infection'))  # these datasets are available by request. Uncomment this line if you don't have them 
#      next
    switch(datasets[d],  # read dataset 
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
           input.data <- input.data[,c(1:2)]
           input.data[,2] <- input.data[,2]-0.02 # correction: reduce by 1: the truncation criterion is L<TU (L==TU is removed)
         }, 
         'Infection'={ # this dataset is not part of the released package
           load('../data/ICU_INF.Rdata')
           input.data <- cbind(X,Y)
           W.max[d] <- max(X)+max(Y)
         }
  ) # end switch 
  
  prms <- list(B = 100)
  prms$W.max <- W.max[d] 
  for(t in 1:n.tests) # run all tests 
  {
    if((!(w.fun[d] %in% c('truncation', 'Hyperplane_Truncation'))) & (test.type[t] %in% c("tsai", 'minP2')))
      next  # these tests run only for truncation 
    if((test.type[t] == 'bootstrap') & (w.fun[d] %in% c('truncation', 'Hyperplane_Truncation', 'huji')))
      next # can't run bootstrap because w can be zero 
    set.seed(1)
    
    print(paste0(datasets[d], ", ", test.type[t], ":"))
    start.time <- Sys.time()
    results.test<-TIBS(input.data, w.fun[d], test.type[t], prms)
    test.time[d,t] <- Sys.time() - start.time
    test.pvalue[d,t] <- results.test$Pvalue 
    cat(datasets[d], ', ', test.type[t], ', Pvalue:', test.pvalue[d,t], '\n')
    if(plot.flag & (t==2)) # permutations
    {
      xy <- cbind(data.frame(input.data), as.data.frame(results.test$permuted.data))
      ggplot(xy, aes(x=xy[,1], y=xy[,2])) + 
        geom_point(aes(x=xy[,1], y=xy[,2], col="original")) + 
        geom_point(shape=3, aes(x=xy[,3], y=xy[,4], col="permuted")) + 
        ggtitle(datasets[d]) +
        xlab("X") + ylab("Y") +
        theme(plot.title = element_text(size=14, face="bold.italic", hjust=0.5),
              axis.title.y = element_text(face="bold", size=14),
              axis.title.x = element_text(face="bold", size = 14),
              axis.text.x = element_text(face="bold", size=12), 
              axis.text.y = element_text(face="bold", size=12), 
              legend.position = "none")
      ggsave(paste(path, '/../Figures/real_data/', datasets[d], '.png', sep=''),
          width=5, height=5, dpi=300)
    }
  } # end loop on tests 
} # end loop on datasets   

# save table with p-values and running times  
results.table <- cbind(test.pvalue, test.time)
rownames(test.pvalue) <- datasets
rownames(test.time) <- datasets
colnames(test.pvalue) <- test.type
colnames(test.time) <- test.type

rownames(results.table) <- datasets
# colnames(results.table) <- test.type
save(test.pvalue, test.time, test.type, 
     file=paste0(path ,'/../docs/Tables/real_datasets_B_', B, '.Rdata'))
# print(xtable(results.table, type = "latex"), 
#      file = paste0(path ,'/../docs/Tables/real_datasets.tex')) # save also in latex format 

