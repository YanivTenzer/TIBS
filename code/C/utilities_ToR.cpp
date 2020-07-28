#include <iostream>
#include <math.h>
// #include <Rcpp.h> // for including R with Rcpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "utilities_ToR.h"


using namespace std;
using namespace arma;
using namespace Rcpp;  // for including R with Rcpp 

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). ...
// [[Rcpp::export]]
int timesEight(int x) {
	return x * 8;
}

// [[Rcpp::export]]
int plusFive(int x) {
	return x + 5;
}

// [[Rcpp::export]]
Rcpp::NumericVector add_one_sqrt(NumericVector x) {
	NumericVector y(x.length()); //	y(x);
	y = sqrt(x + 1.0);
	return y;
}

// [[Rcpp::export]]
arma::vec arma_sort(arma::vec x, arma::vec y) {
	return x(arma::sort_index(y));
}

// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
	NumericVector y = clone(x);
	std::sort(y.begin(), y.end());
	return y;
}



// [[Rcpp::export(".mm")]]
arma::mat mm_mult(const arma::mat& lhs,
	const arma::mat& rhs)
{
	return lhs * rhs;
}

// [[Rcpp::export(".vm")]]
arma::mat vm_mult(const arma::vec& lhs,
	const arma::mat& rhs)
{
	return lhs.t() * rhs;
}

// [[Rcpp::export(".mv")]]
arma::mat mv_mult(const arma::mat& lhs,
	const arma::vec& rhs)
{
	return lhs * rhs;
}

// [[Rcpp::export(".vv")]]
arma::mat vv_mult(const arma::vec& lhs,
	const arma::vec& rhs)
{
	return lhs.t() * rhs;
}

/*
* Copied from here: https://figshare.com/articles/Algorithm_of_binary_search_Rcpp_code/3394741
* A binary search divides a range of values(sorted array) into halves by the median point(noted k), and continues
 * tonarrow down the field of search until the pattern(noted x) is found.
 * It is the classic example of a "divide and conquer" algorithm.
 */
 // [[Rcpp::export]]
int binary_search_rcpp(NumericVector array, double pattern) {
	int array_length = array.size();
	int i = 0, j = array_length - 1;
	while (i <= j) {
		int k = (i + j) / 2;
		if (array[k] == pattern) {
			return k; // return index 
		}
		else if (array[k] < pattern) {
			i = k + 1;
		}
		else {
			j = k - 1;
		}
	}
	return -1;
}


// compute empirical cdf 
// [[Rcpp::export]]
NumericVector empirical_cdf_rcpp(NumericVector x)
{
	long i;
	long n = x.length();
	IntegerVector Px(n); 
	NumericVector ecdf(n);

	Px = sort_indexes_rcpp(x); 
	for (i = 0; i < n; i++)
		ecdf[Px[i]] = double(i+1) / double(n); // use double division here
	return(ecdf);
}



// [[Rcpp::export]]
double w_fun_eval_rcpp(double x, double y, string w_fun)
{
	double r = 0.0;

	// allow only string for now 
	// different bias functions (instead of switch)
	if (w_fun == "truncation")
		r = x < y;
	if (w_fun == "Hyperplane_Truncation")
		r = x < y;
	if (w_fun == "exp")
		r = exp((-abs(x) - abs(y)) / 4);
	if (w_fun == "exponent_minus_sum_abs")
		r = exp((-abs(x) - abs(y)) / 4);
	if (w_fun == "stritcly_positive")
		r = exp((-abs(x) - abs(y)) / 4);
	if (w_fun == "huji")
		r = fmax(fmin(65.0 - x - y, 18.0), 0.0);
	if (w_fun == "sum")
		r = x + y;
	if (w_fun == "naive")
		r = 1;

	return(r);
}


// [[Rcpp::export]]
NumericMatrix w_fun_to_mat_rcpp(NumericMatrix data, string w_fun)
{
	long n = data.nrow();  // get sample size
	NumericMatrix w_mat(n, n);
	long i, j; 
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			w_mat(i, j) = w_fun_eval_rcpp(data(i, 0), data(i, 1), w_fun);
	
	return (w_mat);
}



//##################################################################################################
// Compute the modified Hoeffding's test statistic, for the permutation test
// Parameters:
// n - number of samples
// data - n * 2 matrix with(x, y) sample
// grid_points - all possible(x_i, y_j) points.size
// null.expectations.table - 4 * n mass table with pre - computed mass estimates
// 
//  Quardants convension :
//   4 | 1
//   ------
//   3 | 2
//##################################################################################################
// [[Rcpp::export]]
double ComputeStatistic_rcpp(NumericMatrix data, NumericMatrix grid_points, NumericMatrix null_expectations_table)
{
	long n = data.nrow();
	double Obs[4] = { 0 };
	double Exp[4] = { 0 };
	long i, j;
	double Statistic = 0.0;
	long* Rx = new long[n];
	long* Ry = new long[n];

	for (i = 0; i < n; i++) // Slow loop on grid points
	{
		for (j = 0; j < 4; j++) // loop on quadrants
		{
			Exp[j] = null_expectations_table(i, j);
			Obs[j] = 0;
		}
		for (j = 0; j < n; j++)  // loop on data points  
		{
			Rx[j] = data(j, 0) > grid_points(i, 0);
			Ry[j] = data(j, 1) > grid_points(i, 1);
			Obs[0] += double(Rx[j] * Ry[j]);
			Obs[1] += Rx[j];
			Obs[3] += Ry[j];
		}
		Obs[1] -= Obs[0];
		Obs[3] -= Obs[0];
		Obs[2] = n - Obs[0] - Obs[1] - Obs[3];

		if ((Exp[0] > 1) && (Exp[1] > 1) && (Exp[2] > 1) && (Exp[3] > 1))
			for (j = 0; j < 4; j++)
				Statistic += pow((Obs[j] - Exp[j]), 2) / Exp[j];  // set valid statistic when expected is 0 or very small
	} // end loop on grid points

//	cout << "Return Stat: " << Statistic << endl;
	return(Statistic);
}


// Weighted statistic - alternative version (works only for positive w)
// [[Rcpp::export]]
double ComputeStatistic_w_rcpp(NumericMatrix data, NumericMatrix grid_points, string w_fun) 
{
	long n = data.nrow();
	long i, j;
	NumericMatrix w_vec(n); 
	double n_w=0.0;
	double Statistic = 0.0;
	IntegerVector Rx(n);
	IntegerVector Ry(n);




	for (i = 0; i < n; i++)
	{
		w_vec[i] = w_fun_eval_rcpp(data(i, 0), data(i, 1), w_fun);
		n_w += 1.0 / w_vec[i];
	}
	double Obs[4] = { 0 };
	double Exp[4] = { 0 };

	double Rx_sum, Ry_sum, Rx_not_sum, Ry_not_sum;

	for (i = 0; i < n; i++) // Slow loop on grid points
	{
		Rx_sum = Ry_sum = Rx_not_sum = Ry_not_sum = 0.0;
		for(j = 0; j < 4; j++)
			Obs[j] = 0;
		for (j = 0; j < n; j++)  // loop on data points  
		{
			Rx[j] = data(j, 0) > grid_points(i, 0);
			Ry[j] = data(j, 1) > grid_points(i, 1);

			Rx_sum += Rx[j] / w_vec[j];
			Ry_sum += Ry[j] / w_vec[j];
			Rx_not_sum += (1-Rx[j]) / w_vec[j];
			Ry_not_sum += (1-Ry[j]) / w_vec[j];

			Obs[0] += (Rx[j] * Ry[j] / (w_vec[j] * n_w));
			Obs[1] += (Rx[j] * (1-Ry[j]) / (w_vec[j] * n_w));
			Obs[2] += ((1-Rx[j]) * (1-Ry[j]) / (w_vec[j] * n_w));
			Obs[3] += ((1-Rx[j]) * Ry[j] / (w_vec[j] * n_w));
		}
		Exp[0] = Rx_sum * Ry_sum / (n_w * n_w);
		Exp[1] = Rx_sum * Ry_not_sum / (n_w * n_w);
		Exp[2] = Rx_not_sum * Ry_not_sum / (n_w * n_w);
		Exp[3] = Rx_not_sum * Ry_sum / (n_w * n_w);

		if ((Exp[0] > 1) && (Exp[1] > 1) && (Exp[2] > 1) && (Exp[3] > 1))
			for (j = 0; j < 4; j++)
				Statistic += pow((Obs[j] - Exp[j]), 2) / Exp[j];  // set valid statistic when expected is 0 or very small

	}

	return(Statistic);
}


// [[Rcpp::export]]
List GetNullDistribution_rcpp(Rcpp::NumericMatrix pdfs, Rcpp::NumericMatrix w_mat) // why null distribution is used?
{
	long n = pdfs.nrow();
	// Compute the normalizing factor under the null :
	long i, j;
	double z = 0.0;
	NumericMatrix null_distribution(n, n);  // set size 

	if (w_mat.nrow() == 1)  // same w for all 
	{
//		cout << "One WMAT" << endl;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				null_distribution(i, j) = w_mat(1, 1) * pdfs(i, 0) * pdfs(j, 1);
				z += null_distribution(i, j);
			}
	}
	else
	{
//		cout << "Matrix WMAT" << endl;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				null_distribution(i, j) = w_mat(i, j) * pdfs(i, 0) * pdfs(j, 1);
				z += null_distribution(i, j);
			}
	}
//	cout << "Z inside: " << z << endl;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			null_distribution(i, j) /= z;  // normalize 

	List ret;
	ret["Z"] = z; // uppercase Z
	ret["distribution"] = null_distribution;
	return(ret); 
}


// Hoeffding's test in Rcpp

// Get indices of vector sorted 
// [[Rcpp::export]]
IntegerVector sort_indexes_rcpp(NumericVector data) 
{

	IntegerVector index = seq_along(data) - 1;

	sort(index.begin(), index.end(),
		[&](const int& a, const int& b) {
			return (data[a] < data[b]);
		});
	return(index);
}



// double PDFToCDF2d(double** pdf_2d, double* data[2], long n, double** cdf_2d)
// [[Rcpp::export]]
NumericMatrix PDFToCDF2d_rcpp(NumericMatrix pdf_2d, NumericMatrix data)
{
	long i, j;

	long n = pdf_2d.nrow();
	IntegerVector Px(n);
	IntegerVector Py(n);

	NumericVector data0 = data(_, 0);
	NumericVector data1 = data(_, 1);

	Px = sort_indexes_rcpp(data0); //  data(_, 0));
	Py = sort_indexes_rcpp(data1); //  data(_, 1));

	NumericMatrix cdf_2d(n, n);
	NumericMatrix cdf_2d_ret(n, n);

	// cumulative sum on rows 
	for (i = 0; i < n; i++)
	{
		cdf_2d(Px[i], Py[0]) = pdf_2d(Px[i], Py[0]);
		for (j = 1; j < n; j++)
			cdf_2d(Px[i], Py[j]) = cdf_2d(Px[i], Py[j - 1]) + pdf_2d(Px[i], Py[j]);
	}
	// cumulative sum on columns 
	for (i = 1; i < n; i++)
		for (j = 0; j < n; j++)
			cdf_2d(Px[i], Py[j]) = cdf_2d(Px[i - 1], Py[j]) + cdf_2d(Px[i], Py[j]);

	return(cdf_2d); // don't permute back
}



// double GetQuarterExpectedProb(double Point[2], long QId, double* data[2], double** null_distribution_CDF, long n)
// [[Rcpp::export]]
double GetQuarterExpectedProb_rcpp(NumericVector Point, long QId, NumericMatrix data, NumericMatrix null_distribution_CDF)
{
	// convert below R to C++:
	long idx_x = -1, idx_y = -1, i;
	double S = 0.0;
	long n = data.nrow();

	if ((QId == 0) || (QId == 1)) // work with 0-3 coordinates for quadrants
	{
		for (i = 0; i < n; i++)
			if (data(i, 0) > Point[0])
				if ((idx_x == -1) || (data(idx_x, 0) > data(i, 0)))  // take minimal value larger than Point[0]
					idx_x = i;
	}
	else
	{
		for (i = 0; i < n; i++)
			if (data(i, 0) <= Point[0])
				if ((idx_x == -1) || (data(idx_x, 0) < data(i, 0)))
					idx_x = i;
	}
	if ((QId == 0) || (QId == 3))
	{
		for (i = 0; i < n; i++)
			if (data(i, 1) > Point[1])
				if ((idx_y == -1) || (data(idx_y, 1) > data(i, 1)))
					idx_y = i;
	}
	else
	{
		for (i = 0; i < n; i++)
			if (data(i, 1) <= Point[1])
				if ((idx_y == -1) || (data(idx_y, 1) < data(i, 1)))
					idx_y = i;
	}

	if ((idx_x == -1) || (idx_y == -1))
		return(0);

	long idx_x_max = which_max(data(_, 0));
	long idx_y_max = which_max(data(_, 1));

	switch (QId) { // different quardants
	case 0: {S = null_distribution_CDF(idx_x_max, idx_y_max) + null_distribution_CDF(idx_x, idx_y) -
		null_distribution_CDF(idx_x, idx_y_max) - null_distribution_CDF(idx_x_max, idx_y); break; } // "huji" prints "1",
	case 1: {S = null_distribution_CDF(idx_x_max, idx_y) - null_distribution_CDF(idx_x, idx_y); break; }
	case 2: {S = null_distribution_CDF(idx_x, idx_y); break; }
	case 3: {S = null_distribution_CDF(idx_x, idx_y_max) - null_distribution_CDF(idx_x, idx_y); break; }
	} // end switch 
	
	return(S);
}


//double QuarterProbFromBootstrap(double* data[2], double** null_distribution, double* grid_points[2], long n,
//	double* mass_table[4])
// [[Rcpp::export]]
NumericMatrix QuarterProbFromBootstrap_rcpp(NumericMatrix data, NumericMatrix null_distribution, NumericMatrix grid_points)
{
	long i, j;
	long n = data.nrow();
	long n_grid = grid_points.nrow();
	NumericMatrix null_distribution_CDF(n, n);
	NumericVector cur_grid_points(2);
	NumericMatrix mass_table(n_grid, 4);

	null_distribution_CDF = PDFToCDF2d_rcpp(null_distribution, data);  // convert PDF to CDF 

//	cout << "n: " << n << " n_grid: " << n_grid << " null-dist: " << null_distribution.nrow() << ", " << null_distribution.ncol() << endl;
	for (i = 0; i < n_grid; i++)  // loop on grid-points
	{
		cur_grid_points = grid_points(i,_); // 		cur_grid_points[1] = grid_points(i,1);
		for (j = 0; j < 3; j++)
			mass_table(i, j) = n * GetQuarterExpectedProb_rcpp(cur_grid_points, j, data, null_distribution_CDF);
		mass_table(i,3) = n - mass_table(i,2) - mass_table(i,1) - mass_table(i,0);
	}

	return(mass_table);
}  // end function QuarterProbFromBootstrap




// [[Rcpp::export]]
NumericMatrix QuarterProbFromPermutations_rcpp(NumericMatrix data, NumericMatrix P, NumericMatrix grid_points)
{
	long i, j, k, a, b;
	long n = data.nrow();
	NumericMatrix mass_table(n, 4);
	for (i = 0; i < grid_points.nrow(); i++) // loop over grid points 
	{
		for (j = 0; j < n; j++)
		{
			a = data(j, 0) > grid_points(i, 0);
			for (k = 0; k < n; k++)
			{
				b = data(k, 1) > grid_points(i, 1);
				mass_table(i, a + 2 * b)++;
			}
		}
	}
	return(mass_table);
}



/****************************************************************************************************************************/
// Here we add ALL function files from multiple files into a single cpp file (to enable easier compilation with sourcecpp).
// This should be replaced later with a pakcage. 
/****************************************************************************************************************************/

/**
##############################################################################
# Convert marginal PDF to CDF from data
# Parameters:
# data - k * n array of variables
# PDF.table - matrix with every column a vector of PDF of each variable
# Output:
# CDF.table - matrix with every column a vector of CDF of each variable
##############################################################################
**/
// [[Rcpp::export]]
NumericMatrix PDFToCDFMarginals_rcpp(NumericMatrix data, NumericMatrix PDFs)
{
	long c, i;
	long n = PDFs.nrow(); 
	NumericMatrix CDFs(PDFs.nrow(), PDFs.ncol()); //	CDF.table < -array(0L, dim(PDF.table))  # matrix(0, num.samples, num.variables)
	IntegerVector Px(n);

	for (c = 0; c < PDFs.ncol(); c++)
	{
//		cout << "Run Column: " << c << endl;
		Px = sort_indexes_rcpp(data(_, c)); // get indexes permute to order x_i, y_i
		CDFs(Px[0], c) = PDFs(Px[0], c);
		for(i = 1; i < n; i++)
			CDFs(Px[i], c) = CDFs(Px[i-1], c) + PDFs(Px[i], c); //  cumsum(PDFs(Px, c)); // need to change 
	}
	return(CDFs);
}


// [[Rcpp::export]]
NumericMatrix CDFToPDFMarginals_rcpp(NumericMatrix CDFs)
{
	long c, i;
	long n = CDFs.nrow();
	NumericMatrix PDFs(CDFs.nrow(), CDFs.ncol());
	IntegerVector Px(n);
	NumericVector CDF_sorted(n);

	for (c = 0; c < CDFs.ncol(); c++)
	{
		Px = sort_indexes_rcpp(CDFs(_, c));
		CDF_sorted = CDFs(_, c);
		CDF_sorted = CDF_sorted.sort(); //  sort(CDFs(_, c));
		
		PDFs(Px[0], c) = CDF_sorted[0]; // cumsum
		for (i = 1; i < n; i++)
			PDFs(Px[i], c) = CDF_sorted[i] - CDF_sorted[i - 1];
	}
	return(PDFs);
}




// [[Rcpp::export]]
List EstimateMarginals_rcpp(NumericMatrix data, string w_fun)  // inputs  //	double* xy[2], double* PDFs[2], double* CDFs[2])  // outputs 
{
	List ret;
	
	long i, j;
	long n = data.nrow();

	NumericVector w_inv(n); //	double* w_inv = new double[n];
	double w_inv_sum = 0.0;

	NumericMatrix PDFs(2*n, 2);
	NumericMatrix CDFs(2*n, 2);

//	NumericMatrix CDFs_alt(2 * n, 2);

	long naive_flag = FALSE, pos_flag = FALSE;
	string pos_w[3] = { "sum", "sum_coordinates", "exponent_minus_sum_abs" };
	string naive_w[2] = { "naive", "no_bias" };

	for (i = 0; i < 3; i++)
		if (w_fun == pos_w[i])
			pos_flag = TRUE;
	for (i = 0; i < 2; i++)
		if (w_fun == naive_w[i])
			naive_flag = TRUE;

	if (pos_flag) // w_fun % in % c('sum', 'sum_coordinates', 'exponent_minus_sum_abs')) // for w(x, y) > 0 cases
	{ // case 1: strictly positive W, use ML estimator
		for (i = 0; i < n; i++)
		{
			w_inv[i] = 1.0 / w_fun_eval_rcpp(data(i, 0), data(i, 1), w_fun); 
			w_inv_sum += w_inv[i];
		}
		for (i = 0; i < n; i++)   // normalize
			PDFs(i, 0) = PDFs(i, 1) = w_inv[i] / w_inv_sum;
		PDFs = PDFs(Range(0, n-1), _);
		CDFs = PDFToCDFMarginals_rcpp(data, PDFs);  // why null ?
		ret["xy"] = data;
	}
	else {
		// skip survival 
		if (naive_flag) // no bias(assume W(x, y) = 1)
		{
			CDFs = CDFs(Range(0, n - 1), _);
			CDFs(_, 0) = empirical_cdf_rcpp(data(_, 0));
			CDFs(_, 1) = empirical_cdf_rcpp(data(_, 1));
			ret["xy"] = data;
		}
		else { // use W 
				// general biased sampling function(w.fun can be a function not a string) with exchangable distributions
			// case 2: left truncation, use the estimator of Proposition 1 in the paper for exchangable distributions
					// Augment data to include both xand y values for each axis(due to symmetry)
			double augment_data = TRUE;   // new: add xy values
			NumericMatrix new_data(2*n, 2); // largest 
			NumericMatrix new_data_sorted(2 * n, 2);
			if (augment_data) // duplicate x and y values
			{
				for (i = 0; i < 2; i++)
					for (j = 0; j < n; j++)
					{
						new_data(j, i) = data(j, 0);
						new_data(j+n, i) = data(j, 1);
					}
			}
			else // just copy data 
			{
				for (i = 0; i < 2; i++)
					for (j = 0; j < n; j++)
						new_data(j, i) = data(j, i);
			}

			// need here permutations !! 
			NumericVector Px(2 * n);
			NumericVector Py(2 * n);
			Px = sort_indexes_rcpp(new_data(_, 0));
			Py = sort_indexes_rcpp(new_data(_, 1));

	//		Rcout << "Generated Px, Py" << endl; 

			for (i = 0; i < 2*n; i++) // copy sorted 
			{
				new_data_sorted(i, 0) = new_data(Px[i], 0);
				new_data_sorted(i, 1) = new_data(Py[i], 1);
			}
			NumericVector F0 = empirical_cdf_rcpp(new_data(_, 0)); // can be used instead of binary search 
			NumericVector F1 = empirical_cdf_rcpp(new_data(_, 1));
			// alternative: take average of F0 and F1
			CDFs(_, 0) = (F0 + F1) / 2;
			CDFs(_, 1) = CDFs(_, 0); //  (F0 + F1) / 2;

			/** binary search not needed
			long F01, F10;
			for (i = 0; i < 2*n; i++) // loop to 2*n
			{
				F01 = binary_search_rcpp(new_data_sorted(_, 0), new_data(i, 1));  // binary_search_rcpp returns the index 
				F10 = binary_search_rcpp(new_data_sorted(_, 1), new_data(i, 0));

				CDFs(i, 0) = ((F01 + F10) / 2.0 + 1.0) / (2*n);
				CDFs(i, 1) = ((F01 + F10) / 2.0 + 1.0) / (2*n);  // need to change !!! 

//				Rcout << "Set i=" << i << endl; 
			}
			**/

			ret["xy"] = new_data;
		} // end if naive w

//		Rcout << " Now CDF-PDF Marginals" << endl; 
		PDFs = CDFToPDFMarginals_rcpp(CDFs); // convert 
	} // end if positive w

	ret["PDFs"] = PDFs;
	ret["CDFs"] = CDFs;	
//	ret["CDFs_alt"] = CDFs_alt;

	return(ret);

}  // end function 


/**
#################################################################
# sample permutations, using MCMC, over the set of valid permutations,
# with respect to the distribution appears in Eq 8
# Parameters:
# W - matrix with weights
# B - number of permutations to draw
# N - sample size(can be read from data or W ? )
# 
# Output:
# PermutationsTable - A matrix representing the sampled permutations
# P - An n * n matrix with P(i, j) = Pr(pi(i) = j) over the sampled permutations
#########################################################################################
**/

// [[Rcpp::export]]
List PermutationsMCMC_rcpp(NumericMatrix w_mat, List prms) // burn.in = NA, Cycle = NA)  # New: allow non - default burn - in
{
	long n = w_mat.nrow();
	NumericMatrix P(n, n);  // New!matrix with P[i] = j estimate
	long i, j, temp;
	double ratio;

	// Set mcmc default sampling parameters
	if (!prms.containsElementNamed("B")) //   ('B' % in % names(prms)))
		prms["B"] = 1000;
	if (!prms.containsElementNamed("burn_in")) //  ('burn.in' % in % names(prms)))
		prms["burn_in"] = 2 * n;
	if (!prms.containsElementNamed("Cycle")) // (!('Cycle' % in % names(prms)))
		prms["Cycle"] = n;

	long Idx = 0, ctr = 1;
	NumericMatrix PermutationsTable(n, long(prms["B"]));
	IntegerVector Perm = seq(0, n);  // 1 to N-1 
	NumericVector switchIdx(2);
	while (Idx < long(prms["B"]))
	{

	// A Metropolis Hastings algorithm with target stationary distribution \pi
	// Choose the two indices to be switched
		switchIdx[0] = rand() % n; 
		switchIdx[1] = switchIdx[0];
		while (switchIdx[1] == switchIdx[0])
			switchIdx[1] = rand() % n; // sample without replacement form 0 to n-1
		i = switchIdx[1];
		j = switchIdx[2];
		ratio = w_mat(i, Perm[j]) * w_mat(j, Perm[i]) / (w_mat(i, Perm[i]) * w_mat(j, Perm[j]));


		if(	double(rand()) / RAND_MAX < fmin(1.0, ratio) ) // we accept the transition with probability min(1, ratio)
		{
			temp = Perm[i];  // SWAP
			Perm[i] = Perm[j];
			Perm[j] = temp;
			if (ctr == long(prms["burn_in"]) || ((ctr % long(prms["Cycle"]) == 0) && (ctr > long(prms["burn_in"]))))
			{
				for(i=0; i<n; i++)
					PermutationsTable(i,Idx) = Perm[i] + 1;  // add 1 to move from C zero-based indexing to R-one based indexing
				Idx++; 
//				if (Idx%100 == 0)
//					Rcout << "Sample Perm=" << Idx << endl;
			}
			for (i = 0; i < n; i++)
				P(i, Perm[i])++; // Update table counts P
			ctr++;   // update counter only after swap
		}
	}  // end while
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			P(i, j) = P(i, j) / (ctr - 1); //  normalize

	List ret;
	ret["PermutationsTable"] = PermutationsTable;
	ret["P"] = P;
	return(ret); // return list with also P
}


/**
##############################################################################
# Iterative Algorithm for marginal estimation from samples of F_{ XY }^ {(w)}
# Parameters:
# data - 2 * n array of(X, Y)
# CDF.table - vector of CDF of X and Y
# Output:
# 
##############################################################################
**/
// [[Rcpp::export]]
NumericMatrix iterative_marginal_estimation_rcpp(NumericMatrix data, string w_fun)
{
	long n = data.nrow();
	long i; 
	double epsilon = 0.0000001;
	long iters = 1000;
	NumericVector f_x(n);
	NumericVector f_y(n);
	NumericVector f_x_prev(n);
	NumericVector f_y_prev(n);
	NumericVector temp_f_x(n);
	NumericVector temp_f_y(n);

	double change;


	NumericMatrix w_mat = w_fun_to_mat_rcpp(data, w_fun);

	arma::mat w_mat_arma = as<arma::mat>(w_mat);

	long t;
	for (t = 0 ; t < iters; t++)
	{
		
		f_x_prev = f_x;
		f_y_prev = f_y;
		
		temp_f_x = mv_mult(w_mat_arma, f_y);
		for (i = 0; i < n; i++)
			f_x[i] = 1.0 / temp_f_x[i];
		temp_f_y = mv_mult(w_mat_arma.t(), f_x);
		for (i = 0; i < n; i++)
			f_y[i] = 1.0 / temp_f_y[i];

			
//			f_x = 1.0 / (w_mat * f_y);
//		f_y = 1.0 / (w_mat.t() * f_x); // multiply by transpose 

		change = 0.0; 
		for (i = 0; i < n; i++)
			change += (abs(f_x[i] - f_x_prev[i]) + abs(f_y[i] - f_y_prev[i]));
		if (change < epsilon)
			break;
	}
	NumericMatrix f_xy(n, 2);
	f_xy(_, 0) = f_x;
	f_xy(_, 1) = f_y;
	return(f_xy);  // PDF.table
}

/** R Version Below
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


// List of needed functions here: 
// accumulate - NOT NEEDED !
// Completed (?) 
// (V) PermutationsMCMC <- function(W, N, prms) # burn.in=NA, Cycle=NA)  # New: allow non-default burn-in 
// (V) empirical_cdf
// (V) w_fun_eval
// (V) binary_search 
// (V) double PDFToCDF2d(double** pdf_2d, double* data[2], long n, double** cdf_2d)
// (V) QuarterProbFromBootstrap <- function(data, null.distribution, grid.points)
// (V) QuarterProbFromPermutations <- function(data, P, grid.points) #Permutations
// (V) CDFToPDFMarginals <- function(CDF.table)
// (V) PDFToCDFMarginals <- function(data, PDF.table)
// (V) EstimateMarginals <- function(data, w.fun, prms=c())  - Tough (complicated) Depends on many functions 




