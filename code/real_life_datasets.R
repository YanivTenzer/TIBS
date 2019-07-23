# Set path to your working directory: 
# The script contains real data analysis 
# Due to privacy issues, the data files for the ICU datasets are not part of the released package. 
# Please contact the authors if you are interested in them. 
path = 'C:\\Users\\Or Zuk\\Dropbox\\BiasedSampling\\Code' # modify to your path  # path = 'D:/cond_ind_2019'
setwd(path)

source('TIBS.R')
source('simulate_biased_sample.R')
source('marginal_estimation.R')
source('Tsai_test.R')
source('utilities.R') 
library(lubridate)
library(permDep)
library(tictoc)
library(xtable)

num.statistics = 100 # 10000
plot_flag <- 1

test.type <- c('bootstrap', 'permutations', 'tsai', 'minP2') # different tests to run 
datasets <- c('huji', 'AIDS', 'ICU', 'Infection')
exchange_type <- c(FALSE, FALSE, FALSE, FALSE) # no reason to assume real data is exchangeable
bias.method <- c('huji', 'Hyperplane_Truncation', 'Hyperplane_Truncation', 'sum')
prms = c()
W_max <- c(65, 1, 1, -1) # -1 denotes calculate max of w from data  
n_datasets <- length(datasets)
n_tests <- length(test.type)

Hyperplane_prms<-c(-1,1,0)
bias.params<-list(NA, list(Hyperplane_prms=Hyperplane_prms), list(Hyperplane_prms=Hyperplane_prms))

test.pvalue <- matrix(-1, n_datasets, n_tests) # -1 denotes test wasn't performed
test.time <- matrix(-1, n_datasets, n_tests) # -1 denotes test wasn't performed

for(d in 1:n_datasets) # loop on datasets.
{
    if(datasets[d] %in% c('ICU', 'Infection'))  # datasets available by request. Remove this line if you have them 
      next
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
         }, 
         'Infection'={ # this dataset is not part of the released package
           load('../data/ICU_INF.Rdata')
           input.data <- cbind(X,Y)
           W_max[d] <- max(X)+max(Y)
         }
  ) # end switch 
  
  prms <- c()
  prms$W_max <- W_max[d] 
  for(t in 1:n_tests) # run all tests 
  {
    if((!(bias.method[d] %in% c('truncation', 'Hyperplane_Truncation'))) & (test.type[t] %in% c("tsai", 'minP2')))
      next  # these tests run only for truncation 
    if((test.type[t] == 'bootstrap') & (bias.method[d] %in% c('truncation', 'Hyperplane_Truncation', 'huji')))
      next # can't run bootstrap because w can be zero 
    
    set.seed(1)
    
    print(paste0(datasets[d], ", ", test.type[t], ":"))
    start.time <- Sys.time()
    results.test<-TIBS(input.data, bias.method[d], bias.params[[d]], num.statistics,                                       
                       test.type[t], prms)
    test.time[d,t] <- Sys.time() - start.time
    test.pvalue[d,t] <- length(which(results.test$statistics.under.null>=results.test$True_T))/num.statistics
    cat(datasets[d], ', ', test.type[t], ', Pvalue:', test.pvalue[d,t], '\n')
    if(plot_flag & (t==2)) # permutations
    {
      bmp(paste(path, '/../Figures/real_data/', datasets[d], '.bmp', sep=''),
          units="in", width=5, height=5, res=300)
      plot(input.data[,1], input.data[,2], main=datasets[d], xlab='x', ylab='y')
      points(results.test$permuted_data[,1], results.test$permuted_data[,2], col='red', pch=3, cex=0.75)
      #      xy_lim <- par("usr")
      #      alpha<-0.4
      #      legend(xy_lim[1]*alpha+xy_lim[2]*(1-alpha), xy_lim[3]*alpha+xy_lim[4]*(1-alpha), legend=c('orig.', 'permuted'),
      #             col=c("black", "red"),  pch=c(1,3)) # No need for legend
      dev.off()      
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
     file=paste0(path ,'/../docs/Tables/real_datasets_B_', num.statistics, '.Rdata'))
#print(xtable(results.table, type = "latex"), 
#      file = paste0(path ,'/../docs/Tables/real_datasets.tex')) # new! save in latex format 

