#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <numeric>      // std::iota
//#include <algorithm>    // std::sort, std::stable_sort
//

//#include <time.h>
#include <math.h>

#include <RcppArmadillo.h> // for including R with Rcpp
// [[Rcpp::depends(RcppArmadillo)]]

#include "utilities_ToR.h"


using namespace std;
using namespace arma;
using namespace Rcpp;  // for including R with Rcpp 


// Marginal estimation
//long EstimateMarginals_rcpp(double* data[2], string w_fun, double* params, long n, // inputs 
//	double *xy[2], double* PDFs[2], double* CDFs[2])

// Need the following functions: empirical_cdf, accumulate, w_fun_eval





/** R Code:
##############################################################################
# Estimate marginals Fx, Fy from biased sample with density proportional to Fxy* W
# Parameters:
# data - 2 * n sample(X, Y)
# w.fun - biased sampling w function
# params - new (optional) : what estimation algorithm to use
#
# Output :
	# xy - data
	# CDFs - cdf of the estimated marginals
	# PDFs - pdfs of the estimated marginals
##############################################################################
	EstimateMarginals < -function(data, w.fun, params = c())
{
	if (!missing(params))
	{
		PDF.table < -iterative_marginal_estimation(data, w.fun)
			# get marginal CDFs from PDFs
	}

	if (w.fun % in % c('sum', 'sum_coordinates', 'exponent_minus_sum_abs')) # for w(x, y) > 0 cases
	{ #case 1: strictly positive W, use ML estimator
	  w.inv < -1 / w_fun_eval(data[,1], data[,2], w.fun)
	  Fx < -Fy < -w.inv / sum(w.inv) # normalize
	  PDF.table = as.data.frame(cbind(Fx, Fy))
	  CDF.table = NULL  # why null ?
	}
	else {
		if (w.fun % in % c('survival'))  # what is the definition of w here ?
		{
			# Estimate marginals using Kaplan - Meier estimator
				n < -dim(data)[1]
				require(survival)
				y.srv < -Surv(time = data[, 1], time2 = data[, 2], event = rep(1, n))
				x.srv < -Surv(time = -data[, 2], time2 = -data[, 1], event = rep(1, n))
				Fx < -survfit(x.srv~1)
				Fy < -survfit(y.srv~1)
				CDF.table < -cbind(rev(Fx$surv), 1 - Fy$surv)
				data < -cbind(rev(-Fx$time), Fy$time)
		}
		else {
			if (w.fun % in % c('naive', 'no_bias')) # no bias(assume W(x, y) = 1)
			{
				Fx < -ecdf(data[, 1])
					Fy < -ecdf(data[, 2])
					CDF.table < -cbind(Fx(data[, 1]), Fy(data[, 2]))
			}
			else {
				# general biased sampling function(w.fun can be a function not a string) with exchangable distributions
	#case 2: left truncation, use the estimator of Proposition 1 in the paper for exchangable distributions
					# Augment data to include both xand y values for each axis(due to symmetry)
					augment.data < -1 # new: add xy values
					if (augment.data) # duplicate xand y values
						data < -cbind(union(data[, 1], data[, 2]), union(data[, 1], data[, 2]))
						F1 < -ecdf(data[, 1])
						F2 < -ecdf(data[, 2])
						Fx < -(F1(data[, 1]) + F2(data[, 1])) / 2  # Fx, Fy are the same CDFs evaluated at different data x, y
						Fy < -(F1(data[, 2]) + F2(data[, 2])) / 2
						CDF.table < -cbind(Fx, Fy)
			}
		} # end if
	}  # else on w.fun type
		print('Estimating PDF marginal! returning!')
			PDF.table < -CDFToPDFMarginals(CDF.table)
			save(data, PDF.table, CDF.table, file = 'cdfpdf.Rdata')
			return(list(xy = data, CDFs = CDF.table, PDFs = PDF.table)) # new: return also x, y(might be different than original)
}


##############################################################################
# Iterative Algorithm for marginal estimation from samples of F_{ XY }^ {(w)}
# Parameters:
# data - 2 * n array of(X, Y)
# CDF.table - vector of CDF of X and Y
# Output:
#
##############################################################################
iterative_marginal_estimation < -function(data, w.fun)
{
	epsilon < -0.00000001
		iters < -1000
		f_x < -rep(0, n)
		f_y < -rep(0, n)
		w.mat < -w_fun_eval(data[, 1], data[, 2], w.fun)

		for (t in 1 : iters)
		{
			f_x_prev < -f_x
				f_y_prev < -f_y
				f_x < -1 / (w.mat * f_y)
				f_y < -1 / (t(w.mat) * f_x)
				if (sum(abs(f_x - f_x_prev)) + sum(abs(f_x - f_x_prev)) < epsilon)
					break
		}
	return(c(f_x, f_y)) # PDF.table
}

**/


