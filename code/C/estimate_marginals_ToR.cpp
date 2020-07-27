#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <numeric>      // std::iota
//#include <algorithm>    // std::sort, std::stable_sort
//

//#include <time.h>
#include <math.h>
// #include "utilities.h"
#include <Rcpp.h> // for including R with Rcpp


using namespace std;
using namespace Rcpp;  // for including R with Rcpp 

// Marginal estimation
//long EstimateMarginals_rcpp(double* data[2], string w_fun, double* params, long n, // inputs 
//	double *xy[2], double* PDFs[2], double* CDFs[2])


// Need the following functions: empirical_cdf, accumulate, w_fun_eval





// [[Rcpp::export]]
long EstimateMarginals_rcpp(long n, NumericMatrix data, string w_fun, NumericVector params)  // inputs 
//	double* xy[2], double* PDFs[2], double* CDFs[2])  // outputs 

{
	long i, j;

	long naive_flag = FALSE, pos_flag = FALSE;
	double* w_inv = new double[n];
	double w_inv_sum;

	double* Fx = new double[n];
	double* Fy = new double[n];
	double* F0 = new double[n];
	double* F1 = new double[n];

	double* new_data[2];

	NumericMatrix xy(n, 2);
	NumericMatrix PDFs(n, 2);
	NumericMatrix CDFs(n, 2);

	string pos_w[3] = { "sum", "sum_coordinates", "exponent_minus_sum_abs" };
	string naive_w[2] = { "naive", "no_bias" };

	for (i = 0; i < 3; i++)
		if (w_fun == pos_w[i])
			pos_flag = TRUE;
	for (i = 0; i < 2; i++)
		if (w_fun == naive_w[i])
			naive_flag = TRUE;

	cout << "POS FLAG=" << pos_flag << " NAIVE FLAG=" << naive_flag << endl;

	if (pos_flag) // w_fun % in % c('sum', 'sum_coordinates', 'exponent_minus_sum_abs')) // for w(x, y) > 0 cases
	{ // case 1: strictly positive W, use ML estimator
		for(i=0; i<n; i++)
			w_inv[i] = 1.0 / w_fun_eval(data[0][i], data[1][i], w_fun); // NEED TO IMPLEMENT
		w_inv_sum = accumulate(w_inv, w_inv + n, 0.0);
		for (i = 0; i < n; i++)   // normalize
			PDFs(i,0) = PDFs(i,1) = w_inv[i] / w_inv_sum;
	}
	else {
		// skip survival 
		if (naive_flag) // no bias(assume W(x, y) = 1)
		{
			empirical_cdf(data[0], n, CDFs[0]);
			empirical_cdf(data[1], n, CDFs[1]);
		}
		else { // use W 
				// general biased sampling function(w.fun can be a function not a string) with exchangable distributions
			// case 2: left truncation, use the estimator of Proposition 1 in the paper for exchangable distributions
					// Augment data to include both xand y values for each axis(due to symmetry)
			double augment_data = TRUE;   // new: add xy values
			if (augment_data) // duplicate xand y values
			{
				for (i = 0; i < 2; i++)
					new_data[i] = new double[2 * n];
				for (i = 0; i < 2; i++)
				{
					new_data[i] = new double[2 * n];
					for (j = 0; j < n; j++)
					{
						new_data[i][j] = data[0][j];
						new_data[i][j + n] = data[1][j];
					}
				}
			}
			else // just copy data 
			{
				for (i = 0; i < 2; i++)
				{
					new_data[i] = new double[n];
					for (j = 0; j < n; j++)
						new_data[i][j] = data[i][j];
				}
			}
/*
				F1 <- ecdf(data[, 1])
				F2 <- ecdf(data[, 2])
				Fx <- (F1(data[, 1]) + F2(data[, 1])) / 2  # Fx, Fy are the same CDFs evaluated at different data x, y
				Fy <- (F1(data[, 2]) + F2(data[, 2])) / 2
				CDF.table < -cbind(Fx, Fy)
*/

			// need here permutations !! 


			long *Px = new long[2*n];
			long *Py = new long[2*n];

			double* new_data_sorted[2];
			for (i = 0; i < 2; i++)
				new_data_sorted[i] = new double[2*n];

			sort_with_indexes(new_data[0], n, Px);
			sort_with_indexes(new_data[1], n, Py);


			for (i = 0; i < n; i++)
			{
				new_data_sorted[0][i] = new_data[0][Px[i]];
				new_data_sorted[1][i] = new_data[1][Py[i]];
			}

			empirical_cdf(new_data[0], n, F0);
			empirical_cdf(new_data[1], n, F1);

			long F01, F10; 
			for (i = 0; i < n; i++)
			{
				F01 = binary_search(new_data_sorted[0], new_data_sorted[0] + n, new_data[1][i]);
				F10 = binary_search(new_data_sorted[1], new_data_sorted[1] + n, new_data[0][i]);

				CDFs(i,0) = (F0[i] + F1[i]) / 2;
				CDFs(i,1) = (F0[i] + F1[i]) / 2;  // need to change !!! 
			}
		}
	} // end if

	exit(99);
	return(TRUE);
}  // else on w.fun type





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