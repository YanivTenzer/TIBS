path = 'C:\\Users\\Or Zuk\\Dropbox\\BiasedSampling\\Code' # modify to your path  # path = 'D:/cond_ind_2019'
setwd(path)

source('runner_single_iteration.R')
source('simulate_biased_sample.R')
source('marginal_estimation.R')
source('Tsai_test.R')
source('utilities.R')
library(lubridate)
library(permDep)
library(tictoc)
library(xtable)

num_of_statistics = 100 # 10000
plot_flag <- 1

test_type <- c('bootstrap', 'permutations', 'tsai', 'minP2') # different tests to run 
datasets <- c('huji', 'AIDS', 'ICU', 'Infection')
exchange_type <- c(FALSE, FALSE, FALSE, FALSE) # no reason to assume real data is exchangeable
bias_method <- c('huji', 'Hyperplane_Truncation', 'Hyperplane_Truncation', 'sum')
prms = c()
W_max <- c(65, 1, 1, -1) # -1 denotes calculate max of w from data  
n_datasets <- length(datasets)
n_tests <- length(test_type)

Hyperplane_prms<-c(-1,1,0)
bias_params<-list(NA, list(Hyperplane_prms=Hyperplane_prms), list(Hyperplane_prms=Hyperplane_prms))

test_Pvalue <- matrix(-1, n_datasets, n_tests) # -1 denotes test wasn't performed
test_time <- matrix(-1, n_datasets, n_tests) # -1 denotes test wasn't performed

for(d in 4:n_datasets) # loop on datasets
{
  switch(datasets[d],  # read dataset 
         'huji'={
           load('../data/APage_APtime.RData')
           input_data <- HUJI.dat
         },
         'AIDS'={
           library(DTDA)
           data("AIDS")
           input_data <- AIDS[,2:3]
         },
         'ICU'={
           input_data <- read.table("../data/ICU_data.txt", header = TRUE)
           input_data[input_data$delta==0,] # censoring is always at 30
         }, 
         'Infection'={
           load('../data/ICU_INF.Rdata')
           input_data <- cbind(X,Y)
           W_max[d] <- max(X)+max(Y)
         }
  ) # end switch 
  
  prms <- c()
  prms$W_max <- W_max[d] 
  for(t in 1:n_tests) # run all tests 
  {
    if(!(bias_method[d] %in% c('truncation', 'Hyperplane_Truncation')) & (test_type[t] %in% c("tsai", 'minP2')))
      next  # these tests run only for truncation 
    if((test_type[t] == 'bootstrap') & (bias_method %in% c('truncation', 'Hyperplane_Truncation', 'huji')))
      next # can't run bootstrap because w can be zero 
    
    set.seed(1)
    tic(paste(datasets[d], ", ", test_type[t], ":"))
    results_test<-runner_single_iteration(input_data, bias_method[d], bias_params[[d]], num_of_statistics, 
                                          path, test_type[t], prms)
    test_time[d,t] <- toc()
    test_Pvalue[d,t] <- length(which(results_test$statistics_under_null>=results_test$True_T))/num_of_statistics
    cat(datasets[d], ', ', test_type[t], ', Pvalue:', test_Pvalue[d,t], '\n')
    if(plot_flag & (t==2)) # permutations
    {
      #      pdf(paste(path, '/../Figures/real_data/', datasets[t], '.pdf', sep=''), 
      #          units="in", width=5, height=5, res=300)      # save image 
      #      tiff(paste(path, '/../Figures/real_data/', datasets[t], '.tiff', sep=''),
      #           units="in", width=5, height=5, res=300)
      bmp(paste(path, '/../Figures/real_data/', datasets[d], '.bmp', sep=''),
          units="in", width=5, height=5, res=300)
      plot(input_data[,1], input_data[,2], main=datasets[d], xlab='x', ylab='y')
      points(results_test$permuted_data[,1], results_test$permuted_data[,2], col='red', pch=3, cex=0.75)
      #      xy_lim <- par("usr")
      #      alpha<-0.4
      #      legend(xy_lim[1]*alpha+xy_lim[2]*(1-alpha), xy_lim[3]*alpha+xy_lim[4]*(1-alpha), legend=c('orig.', 'permuted'),
      #             col=c("black", "red"),  pch=c(1,3)) # No need for legend
      dev.off()      
    }
  } # end loop on tests 
} # end loop on datasets   

# save table with p-values and running times  
results_table <- cbind(test_Pvalue, test_time)
rownames(test_Pvalue) <- datasets
rownames(test_time) <- datasets
colnames(test_Pvalue) <- test_type
colnames(test_time) <- test_type

rownames(results_table) <- datasets
colnames(results_table) <- test_type
#save(test_Pvalue, test_time, test_type, file=paste0(path ,'/../docs/Tables/real_datasets.Rdata'))
#print(xtable(results_table, type = "latex"), 
#      file = paste0(path ,'/../docs/Tables/real_datasets.tex')) # new! save in latex format 

