#include <iostream>
#include <math.h>
// #include <map>
#include <random>

#include <Rcpp.h> // for including R with Rcpp
// #include "utilities_ToR.h"


using namespace std;
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

/**
// [[Rcpp::export]]
arma::vec arma_sort(arma::vec x, arma::vec y) {
	return x(arma::sort_index(y));
}
**/

// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
	NumericVector y = clone(x);
	std::sort(y.begin(), y.end());
	return y;
}


/**
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
**/



// [[Rcpp::export]]
NumericVector my_mv_mult(NumericMatrix A, NumericVector V)
{
	NumericVector r(A.nrow());
	long i, j;

	for (i = 0; i < A.nrow(); i++)
	{
		r[i] = 0.0;
		for (j = 0; j < A.ncol(); j++)
			r[i] += A(i, j) * V[j];
	}
	return(r);
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


/**
** Copied from here: https://stackoverflow.com/questions/37143283/finding-unique-rows-in-armamat **
template <typename T>
inline bool rows_equal(const T& lhs, const T& rhs, double tol = 0.00000001) {
	return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

// [[Rcpp::export]]
arma::mat unique_rows(const arma::mat& x) {
	unsigned int count = 1, i = 1, j = 1, nr = x.n_rows, nc = x.n_cols;
	arma::mat result(nr, nc);
	result.row(0) = x.row(0);

	for (; i < nr; i++) {
		bool matched = false;
		if (rows_equal(x.row(i), result.row(0))) continue;

		for (j = i + 1; j < nr; j++) {
			if (rows_equal(x.row(i), x.row(j))) {
				matched = true;
				break;
			}
		}

		if (!matched) result.row(count++) = x.row(i);
	}

	return result.rows(0, count - 1);
}
**/

// [[ Rcpp :: export ()]]
//NumericMatrix a4(arma::mat x) {
//	NumericMatrix y = wrap(x);
//	return(y);
//}


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


// should have a similar function for a vector !! 
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


// not needed ? 
/**
double ComputeStatistic_inverse_weighting_rcpp(NumericMatrix data, NumericMatrix grid_points, NumericMatrix w_mat)
{
	long n = data.nrow();
	long i, j, k;
	NumericMatrix w_vec(n);
	double n_w = 0.0;
	double Statistic = 0.0;
	IntegerVector Rx(n);
	IntegerVector Ry(n);
	double Obs[4] = { 0 };
	double Exp[4] = { 0 };
	double Rx_sum, Ry_sum, Rx_not_sum, Ry_not_sum;

	double n_tilde = sum(1 / diag(w_mat));
	double min_Exp = 1.0 / n;
	for (i = 0; i < n; i++) // Slow loop on grid points
		w_vec[i] = w_mat(i, i); // take diagonal 

		for (i = 0; i < n; i++) // Slow loop on grid points
		{
			Rx_sum = Ry_sum = Rx_not_sum = Ry_not_sum = 0.0;
			for (j = 0; j < 4; j++)
				Obs[j] = 0;

			for (j = 0; j < n; j++)  // loop on data points  
			{
				Rx[j] = data(j, 0) > grid_points(i, 0);
				Ry[j] = data(j, 1) > grid_points(i, 1);

				Obs[0] += (Rx[j] * Ry[j] / (w_vec[j] * n_w)); // ComputeStatistic_w_rcpp
				Obs[1] += (Rx[j] * (1 - Ry[j]) / (w_vec[j] * n_w));
				Obs[2] += ((1 - Rx[j]) * (1 - Ry[j]) / (w_vec[j] * n_w));
				Obs[3] += ((1 - Rx[j]) * Ry[j] / (w_vec[j] * n_w));

			}

		  idx1 < -which(Rx * Ry == 1)
		  Obs[1] < -sum(1 / w_vec[idx1])
		  idx2 < -which(Rx * (!Ry) == 1)
		  Obs[2] < -sum(1 / w_vec[idx2])
		  idx3 < -which((!Rx) * (!Ry) == 1)
		  Obs[3] < -sum(1 / w_vec[idx3])
		  idx4 < -which(Ry * (!Rx) == 1)
		  Obs[4] < -sum(1 / w_vec[idx4])

		  Exp[1] < -sum(1 / w_vec[Rx]) * sum(1 / w_vec[Ry])
		  Exp[2] < -sum(1 / w_vec[Rx]) * sum(1 / w_vec[which(!Ry == 1)])
		  Exp[3] < -sum(1 / w_vec[which(!Rx == 1)]) * sum(1 / w_vec[which(!Ry == 1)])
		  Exp[4] < -sum(1 / w_vec[which(!Rx == 1)]) * sum(1 / w_vec[Ry])

			  for (k = 0; k < 4; k++)
			  {
				  Exp[k] /= (n_tilde * n_tilde);
				  Exp[k] /= n_tilde;
			}

		  if (min(Exp) > min_Exp) {
			  Statistic += sum((Obs - Exp) ^ 2 / Exp); // set valid statistic when expected is 0 or very small
		  }
		} // end loop on grid points

	return(Statistic); // return also observed table for diagnostics
}
**/

// [[Rcpp::export]]
List GetNullDistribution_rcpp(NumericMatrix pdfs, NumericMatrix w_mat) // why null distribution is used?
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
	long B = 1000;
	if (prms.containsElementNamed("B")) //   ('B' % in % names(prms)))
		B = prms["B"];
	long burn_in = 2 * n;
	if (prms.containsElementNamed("burn_in")) //  ('burn.in' % in % names(prms)))
		burn_in = prms["burn_in"];
	long Cycle = n; 
	if (prms.containsElementNamed("Cycle")) // (!('Cycle' % in % names(prms)))
		Cycle = prms["Cycle"];

	long Idx = 0, ctr = 1;

//	Rcout << "Perms dimensions: n=" << n << " B = " << B << endl;
	NumericMatrix Permutations(n, long(B));
	IntegerVector Perm = seq(0, n);  // 1 to N-1 
	NumericVector switchIdx(2);
	while (Idx < B)
	{
//		Rcout << "Draw Perm=" << Idx << endl; 
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
			if (ctr == burn_in || ((ctr % Cycle == 0) && (ctr > burn_in)))
			{
				for(i=0; i<n; i++)
					Permutations(i,Idx) = Perm[i] + 1;  // add 1 to move from C zero-based indexing to R-one based indexing
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
	ret["Permutations"] = Permutations;
	ret["P"] = P;
	return(ret); // return list with also P
}

// [[Rcpp::export]]
NumericVector rand_perm(long n)
{
	NumericVector perm(n); 
	long i, j, ctr, r;
	NumericVector is_full(n); 
	for (i = 0; i < n; i++)
		is_full[i] = 0;

	for (i = 0; i < n; i++)
	{
		r = rand() % (n - i);
		ctr = 0;
		for (j = 0; j < r; j++)
			ctr += (1 + is_full[ctr]);

		perm[i] = ctr;
		is_full[perm[i]] = 1;
	}
	return(perm);
}

/**
#############################################################
# Calculate p - value using importance sampling
# (If product of weights is very small/large, we can
	#  multiply the weight function by a constant so the weights
	#  are in some sense centered at 1)
#############################################################
**/

List IS_permute_rcpp(NumericMatrix data, double B, string w_fun) // w = function(x) { 1 }) {
{
	long  n = data.nrow();
	double Tobs = ComputeStatistic_w_rcpp(data, data, w_fun); 
	double reject = 0, sum_p = 0;
	long i, j;
	double pw; 
	double Tb; 
	IntegerVector perm(n);
	NumericMatrix permuted_data(n, 2);
	NumericMatrix w_mat(n, n);
	for (i = 0; i < B; i++) 
	{
		perm = rand_perm(n); //  sample(n); // get a random permutation

		permuted_data(_, 0) = data(_, 0);
		for (j = 0; j < n; j++)
			permuted_data(j, 1) = data(perm[j], 1); // save one example
//		< -data.frame(x = data[1:n, 1], y = data[perm, 2]) // permuted data
		Tb = ComputeStatistic_w_rcpp(permuted_data, permuted_data, w_fun); // grid depends on permuted data
		pw = 1.0;
		for(j=0; j<n; j++)
			pw *= w_fun_eval_rcpp(data(j,0), data(j,1), w_fun);
		reject += (Tb >= Tobs) / pw;
		sum_p += 1 / pw;
	}
	List ret; 
	ret["p_val"] = reject / sum_p;
	ret["Tobs"] = Tobs; 
	return(ret);
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

	// arma::mat w_mat_arma = as<arma::mat>(w_mat);

	long t;
	for (t = 0 ; t < iters; t++)
	{
		
		f_x_prev = f_x;
		f_y_prev = f_y;
		
		temp_f_x = my_mv_mult(w_mat, f_y);
		for (i = 0; i < n; i++)
			f_x[i] = 1.0 / temp_f_x[i];
		temp_f_y = my_mv_mult(transpose(w_mat), f_x); // take transpose 
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


/**
############################################################################################
# Draw a bootstrap sample, given the estimated null distribution Fx, Fy, W
# Use rejection sampling.No need for computing n* n table
# (Problem: what if prob.(rejection) close to 1?)
# Parameters:
# data - n * 2 array with(X, Y) samples
# cdfs - fx and fy - cumulatives (not pdfs!)
# w.fun - Biased sampling function w
# prms - for w max
# n - allow a different sample size
# Output: 
# an n*2 matrix of bootstrap sample 
############################################################################################**/
// [[Rcpp::export]]
NumericMatrix Bootstrap_rcpp(NumericMatrix data, NumericMatrix cdfs, string w_fun, List prms, long n=-1)  // n is optional 
{
	if (n == -1) // is.null(n))
		n = data.nrow();
	NumericMatrix boot_sample = NumericMatrix(n, 2);
//	NumericVector x(n), y(n);
//	NumericVector keep(n);
	double x = 0.0, y = 0.0; //  , r; //  keep;
//	long n_keep = 0;
	long k = 0, ctr=0;
//	long i; //  , j;
	double w_max = prms["W.max"];

//	Rcout << "w_max = " << w_max << endl; 

//	random_device rd; // takes a long time?
//	mt19937 gen(rd());
//	NumericVector v0 = pdfs(_, 0);
//	NumericVector v1 = pdfs(_, 1);

//	vector<double> p0 = as<std::vector<double> >(v0);
//	vector<double> p1 = as<std::vector<double> >(v1);

//	discrete_distribution<> d0(p0.begin(), p0.end()); //  pdfs(_, 0));
//	discrete_distribution<> d1(p1.begin(), p1.end()); //  pdfs(_, 1)); // weighted sampling 
	while (k < n)
	{/**
		r = double(rand()) / RAND_MAX;
		for(i=0; i<cdfs.nrow(); i++)
			if (r <= cdfs(i, 0))
			{
				x = data(i,0);
				break;
			}
		r = double(rand()) / RAND_MAX;
		for (i = 0; i < cdfs.nrow(); i++)
			if (r <= cdfs(i, 1))
			{
				y = data(i,1);
				break;
			} **/
		x = data(binary_search_rcpp(cdfs(_,0), double(rand()) / RAND_MAX), 0); //		rand() data(d0(gen), 0); //  rand() % (n - k); // sample with replacement
		y = data(binary_search_rcpp(cdfs(_, 0), double(rand()) / RAND_MAX), 1); //data(d1(gen), 1); //  rand() % (n - k);
		if (double(rand()) / RAND_MAX < w_fun_eval_rcpp(x, y, w_fun) / w_max) //   exp((-abs(x) - abs(y)) / 4) / w_max) // w_fun_eval_rcpp(x, y, w_fun) / w_max)
		{
			boot_sample(k, 0) = x;
			boot_sample(k++, 1) = y;
		}
		ctr++;
	}

//	cout << "Generated " << ctr << " pairs to get eventually " << n << " samples" << endl; 
	/**
		n_keep = 0;
		// sampling n - k together
		for (j = 0; j < n - k; j++)
		{
			x[j] = data(d0(gen), 0); //  rand() % (n - k); // sample with replacement
			y[j] = data(d1(gen), 1); //  rand() % (n - k);
			keep[j] = rand() < w_fun_eval_rcpp(x[j], y[j], w_fun) / prms["W.max"]);
			n_keep += keep[j];
			if (keep[j])
			{
				boot_sample(k2, 0) = x[j];
				boot_sample(k2, 0) = x[j];
			}		
		}
//		x = data(sample(pdfs.nrow(), n - k, prob = pdfs(_, 0), replace = TRUE), 0); // Sample X ~Fx
//		y = data(sample(pdfs.nrow(), n - k, prob = pdfs(_, 1), replace = TRUE), 1); // Sample Y ~Fy
//		keep = which(as.logical(rbinom(n - k, 1, w_fun_eval_rcpp(x, y, w_fun) / prms["W.max"])));

		if (isempty(keep))
			next;
		boot_sample(1:length(keep) + k, 0) = x[keep];
		boot_sample(1:length(keep) + k, 1) = y[keep];
		k += length(keep);
		**/
	return(boot_sample);
}



// List of needed functions here: 
// unique_matrix - how to perform unique on rows of matrix in cpp?


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







// New: try to write the entire TIBS function in cpp: 

/**
########################################################################
# Perform Test for Independence under general Biased Sampling(TIBS)
# 
# Parameters:
# data - n * 2 array of(x, y) samples(NOT LIST!)
# w.fun - biased sampling function W
# test_type - test to perform
# prms - additional parameters(needed for bootstrap) including B - number of bootstrap / permutation samples to perform
# 
# Output:
# TrueT - test statistic for the data
# statistics_under_null - vector of statistics under the null
########################################################################
**/ /**/
// [[Rcpp::export]]
List TIBS_rcpp(NumericMatrix data, string w_fun, string test_type, List prms)
{

	long ctr;
	double TrueT, NullT;
	List output; 


	// Set defaults
//	long use_cpp = 0;
//	if (prms.containsElementNamed("use_cpp"))
//		use_cpp = prms["use_cpp"];
	long fast_bootstrap = FALSE;
	if (prms.containsElementNamed("fast_bootstrap"))
		fast_bootstrap = prms["fast_bootstrap"];
//	if (!prms.containsElementNamed("minp_eps"))
//		prms["minp_eps"] = NULL; // permDep params. Irrelevant here 
	long PL_expectation = FALSE;
	if (prms.containsElementNamed("PL_expectation"))
		PL_expectation = prms["PL_expectation"];
	long naive_expectation = 0;
	if (prms.containsElementNamed("naive_expectation"))
		naive_expectation = prms["naive_expectation"];
//	if (!prms.containsElementNamed("delta"))
//		prms["delta"] = NA; minP2 params. Irrelevant here

	long B = prms["B"];
	long n = data.nrow();
	long i, j; 

//	Rcout << " Read Input TIBS_RCPP. TEST-TYPE: " << test_type << endl; 


	// 1.Compute weights matrix W : (not needed here, just for permutations test)
	// 2.Create a grid of points, based on the data :
//	arma::mat grid_points_arma = unique_rows(as<arma::mat>(data));  // set unique for ties ? for discrete data

	// NumericMatrix grid_points(grid_points_arma.n_rows, 2); 
	NumericMatrix grid_points(n, 2);
	for (i = 0; i < n; i++) // long(grid_points_arma.n_rows); i++)
		for (j = 0; j < 2; j++)
			grid_points(i, j) = 5.5; //  grid_points_arma(i, j);
//	= as<NumericMatrix>(wrap(grid_points_arma));  // set unique for ties ? for discrete data
	List null_distribution;

	/**/
	// no switch (test_type) for strings in cpp
	if(test_type == "bootstrap_inverse_weighting") 
	{
		List marginals = EstimateMarginals_rcpp(data, w_fun);

//		for (i = 0; i < n; i++)
	//		cout << "i: " << i << "PDF: " << as<NumericVector>(marginals["PDFs"])(i, 0) << ", " << as<NumericVector>(marginals["PDFs"])(i, 1) << 
		//	" CDF: " << as<NumericVector>(marginals["CDFs"])(i, 0) << ", " << as<NumericVector>(marginals["CDFs"])(i, 1) << endl;
		NumericMatrix w_mat = w_fun_to_mat_rcpp(marginals["xy"], w_fun);
		null_distribution = GetNullDistribution_rcpp(marginals["PDFs"], w_mat);
		TrueT = ComputeStatistic_w_rcpp(data, grid_points, w_fun); //  $Statistic
		NumericVector statistics_under_null(B);
		NumericMatrix bootstrap_sample(n, 2);
		for (ctr = 0; ctr < B; ctr++)
		{
			if (ctr % 100 == 0)
				Rcout << "Run Boots=" << ctr << endl;
			bootstrap_sample = Bootstrap_rcpp(marginals["xy"], marginals["CDFs"], w_fun, prms, n); // draw new sample.Problem: which pdf and data ?
			NumericMatrix w_mat_bootstrap = w_fun_to_mat_rcpp(bootstrap_sample, w_fun);
			NullT = ComputeStatistic_w_rcpp(bootstrap_sample, grid_points, w_fun);
			statistics_under_null[ctr] = NullT; //  ["Statistic"] ;
		}
		output["TrueT"] = TrueT;
		output["statistics_under_null"] = statistics_under_null;
	}
	/**/
	if(test_type == "bootstrap") 
	{
//		Rcout << "Start TIBS RCPP Bootstrap" << endl;
		// 3. Estimate the marginals
		List marginals = EstimateMarginals_rcpp(data, w_fun);
//		Rcout << "Estimated Marginals" << endl; 
		NumericMatrix w_mat = w_fun_to_mat_rcpp(marginals["xy"], w_fun); // compute W again for augmented data
		NumericMatrix expectations_table;
		NumericMatrix bootstrap_sample(n, 2);

		// 4. Estimate W(x,y) * Fx * Fy / normalizing.factor
		if (naive_expectation) // here we ignore W(using statistic for unbiased sampling)
		{
			List marginals_naive = EstimateMarginals_rcpp(data, "naive");
			null_distribution = GetNullDistribution_rcpp(marginals_naive["PDFs"], NumericMatrix(1.0)); // scalar matrix of 1
			expectations_table = QuarterProbFromBootstrap_rcpp(marginals_naive["xy"], null_distribution["distribution"], grid_points);
		}
		else
		{
//			Rcout << "Get Null Distribution" << endl;
//			Rcout << "w_mat dim=" << w_mat.nrow() << "Marginals PDFs Dim=" << as<NumericMatrix>(marginals["PDFs"]).nrow() << endl;
			null_distribution = GetNullDistribution_rcpp(as<NumericMatrix>(marginals["PDFs"]), w_mat);
//			Rcout << "Get Expectations Table" << endl;
			expectations_table = QuarterProbFromBootstrap_rcpp(marginals["xy"], null_distribution["distribution"], grid_points);
		}
//		Rcout << "Compute Statistic" << endl;
		// 1. First compute the statistic based on the original data set :
		TrueT = ComputeStatistic_rcpp(data, grid_points, expectations_table);

		// 2. Compute statistic for bootstrap sample :
		NumericVector statistics_under_null(B);
		List null_distribution_bootstrap = null_distribution;
		NumericMatrix w_mat_bootstrap(1,1);
		for (ctr = 0; ctr < B; ctr++) // heavy loop : run on bootstrap
			{
	//			for (i = 0; i < n; i++)
	//				cout << "i: " << i << "PDF: " << as<NumericVector>(marginals["PDFs"])(i, 0) << ", " << as<NumericVector>(marginals["PDFs"])(i, 1) <<
	//				" CDF: " << as<NumericVector>(marginals["CDFs"])(i, 0) << ", " << as<NumericVector>(marginals["CDFs"])(i, 1) << endl;

//				Rcout << "Draw Bootstrap Sample Under Null " << ctr << endl;

				bootstrap_sample = Bootstrap_rcpp(marginals["xy"], marginals["CDFs"], w_fun, prms, n); // draw new sample.Problem: which pdf and data ?

				if (!fast_bootstrap) // re - estimate marginals for null expectation for each bootstrap sample
				{
	//				Rcout << "Estimate Marginals Under Null " << ctr << endl;

					List marginals_bootstrap = EstimateMarginals_rcpp(bootstrap_sample, w_fun);   // Why are the marginals estimated each time ?
				  // 3. Compute weights matrix W :
					if (naive_expectation)
					{
						w_mat_bootstrap(0, 0) = 1.0;   // here we ignore W(using statistic for unbiased sampling)
					}
					else
					{
		//				Rcout << "Compute w_mat Under Null " << ctr << endl;

						w_mat_bootstrap = w_fun_to_mat_rcpp(marginals_bootstrap["xy"], w_fun);
					}
					// 4. Estimate W(x,y) * Fx * FY / normalizing.factor
		//			Rcout << "Get Null Distribution Under Null " << ctr << endl;

					null_distribution_bootstrap = GetNullDistribution_rcpp(marginals_bootstrap["PDFs"], w_mat_bootstrap);
			//		Rcout << "Compute Quarter Under Null " << ctr << endl;

					expectations_table = QuarterProbFromBootstrap_rcpp(
						marginals_bootstrap["xy"], null_distribution_bootstrap["distribution"], grid_points);
				} // if fast bootstrap
		//		Rcout << "Compute Bootstrap Statistic Under Null " << ctr << endl;
				NullT = ComputeStatistic_rcpp(bootstrap_sample, grid_points, expectations_table);

					// null.obs.table  = NullT$obs.table
				statistics_under_null[ctr] = NullT;
			}

		output["TrueT"] = TrueT;
		output["statistics_under_null"] = statistics_under_null;
	}
	/**/
	if(test_type == "permutations")
	{
		NumericMatrix w_mat = w_fun_to_mat_rcpp(data, w_fun);

//		Rcout << "Run Perm MCMC RCPP" << endl;
		List PermutationsList = PermutationsMCMC_rcpp(w_mat, prms);
//		Rcout << "Finished Run Perm MCMC RCPP" << endl;
		NumericMatrix expectations_table;

					
		NumericMatrix P = PermutationsList["P"];
//		Rcout << "Copied Perm P, dim=" << P.nrow() << endl;
		NumericMatrix Permutations = PermutationsList["Permutations"];
//		Rcout << "Copied Perm" << endl;
		NumericMatrix permuted_data(n, 2);


//		Rcout << "Fill Permuted Data " << endl;
		permuted_data(_, 0) = data(_, 0);
		for (i = 0; i < n; i++)
			permuted_data(i, 1) = data(Permutations(i, 0), 1); // save one example

															   //		Rcout << "Finished Permuted Data " << endl; 
		if (naive_expectation) // here we ignore W(using statistic for unbiased sampling)
		{
			List marginals = EstimateMarginals_rcpp(data, "naive");
			NumericMatrix w_mat_scalar(1, 1); w_mat_scalar(1, 1) = 1.0;
			null_distribution = GetNullDistribution_rcpp(marginals["PDFs"], w_mat_scalar);
			expectations_table = QuarterProbFromBootstrap_rcpp(marginals["xy"], null_distribution["distribution"], grid_points);
		}
		else
		{
//			Rcout << "COmpute Expectations Table" << endl; 
			if (PL_expectation)  // get expectations from the bootstrap estimator
			{
				List marginals = EstimateMarginals_rcpp(data, w_fun);
				null_distribution = GetNullDistribution_rcpp(marginals["PDFs"], w_mat);
				expectations_table = QuarterProbFromBootstrap_rcpp(marginals["xy"], null_distribution["distribution"], grid_points);
				w_mat = w_fun_to_mat_rcpp(marginals["xy"], w_fun);
			}
			else
			{
				expectations_table = QuarterProbFromPermutations_rcpp(data, P, grid_points);  // Permutations
			}
		}
//		Rcout << "Compute Statistic" << endl;
		TrueT = ComputeStatistic_rcpp(data, grid_points, expectations_table);

		// Compute the statistics value for each permutation:
		NumericVector statistics_under_null(B);
		for (ctr = 0; ctr < B; ctr++)
		{
//			Rcout << "Compute Permutation Statistic Under Null " << ctr << endl;
			permuted_data(_, 0) = data(_, 0);
			for (i = 0; i < n; i++)
				permuted_data(i, 1) = data(Permutations(i, ctr), 1);
			statistics_under_null[ctr] = ComputeStatistic_rcpp(permuted_data, grid_points, expectations_table);
		}
		output["TrueT"] = TrueT;
		output["statistics_under_null"] = statistics_under_null;
		output["Permutations"] = Permutations;
		output["permuted_data"] = permuted_data;
	} // end permutations test
	/**/
	if(test_type == "permutations_inverse_weighting")
	{
		NumericMatrix w_mat = w_fun_to_mat_rcpp(data, w_fun);
		TrueT = ComputeStatistic_w_rcpp(data, grid_points, w_fun); //  $Statistic
				
		List PermutationsList = PermutationsMCMC_rcpp(w_mat, prms); // burn.in = prms$burn.in, Cycle = prms$Cycle)
		NumericMatrix Permutations = PermutationsList["Permutations"];
		NumericMatrix permuted_data(n, 2);

		// Compute the statistics value for each permutation:
		NumericVector statistics_under_null(B);
		NumericMatrix w_mat_permutation(n, n);
		for (ctr = 0; ctr < B; ctr++)
		{
			permuted_data(_, 0) = data(_, 0);
			for (i = 0; i < n; i++)
				permuted_data(i, 1) = data(Permutations(i, ctr), 1);
//			permuted_sample = cbind(data(_, 0), data(Permutations(_, ctr), 1));
			w_mat_permutation = w_fun_to_mat_rcpp(permuted_data, w_fun);
			NullT = ComputeStatistic_w_rcpp(permuted_data, grid_points, w_fun);
			statistics_under_null[ctr] = NullT;
		}
		output["TrueT"] = TrueT;
		output["statistics_under_null"] = statistics_under_null;
		output["Permutations"] = Permutations;
	} // end permutations with inverse weighting test
	/**/
	if(test_type == "tsai") 
	{ cout << "Can't run Tsao from cpp" << endl; } //   Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}
		
	if(test_type == "minP2") { cout << "Can't run minP2 from cpp" << endl;}
	if(test_type == "importance.sampling") { // new importance sampling permutations test(not working yet)
		cout << "Can't run Importance Sampling from cpp" << endl;} // results = IS.permute(dat, prms$B, w.fun) # W);  // ComputeStatistic.W(dat, grid_points, w.fun)		
// } // end switch 
/**/
	if (!(output.containsElementNamed("permuted_data")))
		output["permuted_data"] = NULL;
/**/
	if (!(output.containsElementNamed("Pvalue"))) //    ("Pvalue" % in % names(output))) // Compute empirical P - value
	{
//		output["Pvalue"] = length(which(output["statistics_under_null"] >= output["TrueT"])) / B;
		double Pvalue = 0.0; 
		NumericVector statistics_under_null = as<NumericVector>(output["statistics_under_null"]);
		for (i = 0; i < B; i++)
			Pvalue += (statistics_under_null[i] >= as<double>(output["TrueT"]));
		output["Pvalue"] = Pvalue / double(B);
	}
	/**/
	return(output);
}