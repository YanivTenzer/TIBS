# Temporary file to test functions
source('marginal_estimation.R')
source('utilities.R')
xy <- t(rbind(c(1, 2, 4, 5), c(1, 3, 7, 2)))

xy.PDF <- cbind(runif(4), runif(4))
xy.PDF <- t(t(xy.PDF) / colSums(xy.PDF)) # normalize 
colSums(xy.PDF)

xy.CDF <- PDFToCDFMarginals(xy, xy.PDF)

xy.CDF2 <- PDFToCDFMarginals(xy, xy.PDF)
xy.CDF2 - xy.CDF

xy.PDF2 <- CDFToPDFMarginals(xy.CDF)

xy.PDF2 - xy.PDF


