//#include <stdio.h>
//#include <stdlib.h> 
#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <numeric>      // std::iota
//#include <algorithm>    // std::sort, std::stable_sort

//#include <time.h>
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
			return k;
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
		ecdf[Px[i]] = i / n;
	return(ecdf);
}



// [[Rcpp::export]]
double w_fun_eval(double x, double y, string w_fun)
{
	double r = 0.0;
	// const char* w_fun_c = w_fun.c_str();

	// allow only string for now 
	// different bias functions (instead of switch)
	if (w_fun == "truncation")
		r = x < y;
	if (w_fun == "Hyperplane_Truncation")
		r = x < y;
	if (w_fun == "exp")
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
double ComputeStatistic_rcpp(long n, Rcpp::NumericMatrix data, Rcpp::NumericMatrix grid_points, Rcpp::NumericMatrix null_expectations_table)
{
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


//double ComputeStatistic_rcpp(long n, Rcpp::NumericMatrix data, Rcpp::NumericMatrix grid_points, Rcpp::NumericMatrix null_expectations_table)
//double GetNullDistribution(double* pdfs[2], double** w_mat, long n, double** null_distribution)

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
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			null_distribution(i, j) /= z;  // normalize 

	List ret;
	ret["z"] = z;
	ret["distribution"] = null_distribution;
	return(ret); //	return(z);
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
NumericMatrix PDFToCDF2d_rcpp(long n, NumericMatrix pdf_2d, NumericMatrix data)
{
	long i, j;

	IntegerVector Px(n);
	IntegerVector Py(n);

	NumericVector data0 = data(_, 0);
	NumericVector data1 = data(_, 1);

	Px = sort_indexes_rcpp(data0); //  data(_, 0));
	Py = sort_indexes_rcpp(data1); //  data(_, 1));

	NumericMatrix cdf_2d(n, n);

	// cumulative sum on rows 
	for (i = 0; i < n; i++)
	{
		cdf_2d(Py[0], Px[i]) = pdf_2d(Py[0], Px[i]);
		for (j = 1; j < n; j++)
			cdf_2d(Py[j], Px[i]) = cdf_2d(Py[j - 1], Px[i]) + pdf_2d(Py[j], Px[i]);
	}
	// cumulative sum on columns 
	for (j = 0; j < n; j++)
	{
		cdf_2d(Py[j], Px[0]) = pdf_2d(Py[j], Px[0]);
		for (i = 1; i < n; i++)
			cdf_2d(Py[j], Px[i]) = cdf_2d(Py[j], Px[i - 1]) + pdf_2d(Py[j], Px[i]);
	}

	// copy results permuted
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			cdf_2d(Py[j], Px[i]) = cdf_2d(j, i);
	return(cdf_2d);

}



// double GetQuarterExpectedProb(double Point[2], long QId, double* data[2], double** null_distribution_CDF, long n)
// [[Rcpp::export]]
double GetQuarterExpectedProb_rcpp(long n, NumericVector Point, long QId, NumericMatrix data, NumericMatrix null_distribution_CDF)
{
	// convert below R to C++:
	long idx_x = -1, idx_y = -1, i;
	double S;

	if ((QId == 0) || (QId == 1)) // change to 0-3 coordinates for quadrants !!
	{
		for (i = 0; i < n; i++)
			if (data(i, 0) > Point[0])
				if ((idx_x == -1) || (data(idx_x, 0) > data(i, 0)))
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

//	cout << "Data for max-index:";
//	for (i = 0; i < 10; i++)
//		cout << data[0][i] << " , " << data[1][i] << endl;

	long idx_x_max = which_max(data(_, 0));
	long idx_y_max = which_max(data(_, 1));

	switch (QId) { // different quardants
	case 0: {S = null_distribution_CDF(idx_y_max, idx_x_max) + null_distribution_CDF(idx_y, idx_x) -
		null_distribution_CDF(idx_y_max, idx_x) - null_distribution_CDF(idx_y, idx_x_max); break; } // "huji" prints "1",
	case 1: {S = -null_distribution_CDF(idx_y, idx_x_max) - null_distribution_CDF(idx_y, idx_x); break; }
	case 2: {S = -null_distribution_CDF(idx_y, idx_x); break; }
	case 3: {S = -null_distribution_CDF(idx_y_max, idx_x) - null_distribution_CDF(idx_y, idx_x); break; }
	} // end switch 

	cout << "Mass Table Indexes: idx_x=" << idx_x << ", idx_y=" << idx_y << ", idx_x_max=" << idx_x_max << ", idx_y_max=" << idx_y_max << " S=" << S << endl;
	return(S);
}


//double QuarterProbFromBootstrap(double* data[2], double** null_distribution, double* grid_points[2], long n,
//	double* mass_table[4])
// [[Rcpp::export]]
double QuarterProbFromBootstrap_rcpp(long n, NumericMatrix data, NumericMatrix null_distribution, NumericMatrix grid_points, NumericMatrix mass_table)
{
	long i, j;
	NumericMatrix null_distribution_CDF(n, n);
	NumericVector cur_grid_points(2);

	//	cout << "Inside DO PDFToCDF2d\n";
	null_distribution_CDF = PDFToCDF2d_rcpp(n, null_distribution, data);  // convert PDF to CDF 

//	cout << "Inside DO GetQuarterExpectedProb\n";
	for (i = 0; i < n; i++)  // loop on grid-points
	{
		cur_grid_points[0] = grid_points(i,0);
		cur_grid_points[1] = grid_points(i,1);
//		cout << "cur-grid-points i=" << i << " " << cur_grid_points[0] << ", " << cur_grid_points[1] << endl;
		for (j = 0; j < 3; j++)
			mass_table(i,j) = GetQuarterExpectedProb_rcpp(n, cur_grid_points, j, data, null_distribution_CDF);
		mass_table(i,3) = 1 - mass_table(i,2) - mass_table(i,1) - mass_table(i,0);
		//		cout << "Inside Finished i=" << i << endl;
	}

	for (j = 0; j < 4; j++)
		for (i = 0; i < n; i++)  //  normalize to counts
			mass_table(i,j) *= n;

	return(TRUE);
}  // end function QuarterProbFromBootstrap



/** R-version 
###################################################################################################
# compute Expect[Qi(p_j)] for 1 <= i <= 4, and all j, given a grid of points using permutations
# Parameters:
# data - 2 * n array(X, Y)
# P - n * n table representing permutations probabilities Pr(pi(i) = j) (# Permutations - set of permutations)
# grid.points - centers of partitions 2 * k array
#
# Output:
# mass.table - a 4 * #grid - points table with quadrants expectation
###################################################################################################
QuarterProbFromPermutations < -function(data, P, grid.points) #Permutations
{

  #next each point has its probability of being selected we evaluate Expect[Qi(p_j)] by summation
  mass.table < -matrix(0, dim(grid.points)[1], 4)
  for (i in seq(1, dim(grid.points)[1],1))
  {
	x < -grid.points[i,]
	mass.table[i,1] < -sum(P[data[,1] > x[1], data[,2] > x[2]])
	mass.table[i,2] < -sum(P[data[,1] > x[1], data[,2] <= x[2]])
	mass.table[i,3] < -sum(P[data[,1] <= x[1], data[,2] <= x[2]])
	mass.table[i,4] < -sum(P[data[,1] <= x[1], data[,2] > x[2]])
  }
  return(mass.table)
}
**/


NumericMatrix QuarterProbFromPermutations(NumericMatrix data, NumericMatrix P, NumericMatrix grid_points)
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


/**
  // next each point has its probability of being selected we evaluate Expect[Qi(p_j)] by summation
  mass.table < -matrix(0, dim(grid.points)[1], 4)

	
  for (i in seq(1, dim(grid.points)[1],1))
  {
	x < -grid.points[i,]
	mass.table[i,1] < -sum(P[data[,1] > x[1], data[,2] > x[2]])
	mass.table[i,2] < -sum(P[data[,1] > x[1], data[,2] <= x[2]])
	mass.table[i,3] < -sum(P[data[,1] <= x[1], data[,2] <= x[2]])
	mass.table[i,4] < -sum(P[data[,1] <= x[1], data[,2] > x[2]])
  }
  return(mass.table)
}
**/

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
		cout << "Run Column: " << c << endl;
		Px = sort_indexes_rcpp(data(_, c)); // get indexes permute to order x_i, y_i
		CDFs(Px[0], c) = PDFs(Px[0], c);
		for(i = 1; i < n; i++)
			CDFs(Px[i], c) = CDFs(Px[i-1], c) + PDFs(Px[i], c); //  cumsum(PDFs(Px, c)); // need to change 
	}
	return(CDFs);
}

/**
PDFToCDFMarginals < -function(data, PDF.table)
{
	n < -dim(PDF.table)[1]  # number of samples
		CDF.table < -array(0L, dim(PDF.table))  # matrix(0, num.samples, num.variables)
		for (i in 1 : dim(PDF.table)[2])  # loop on variables
		{
		  Px < -sort(data[,i], index.return = TRUE)  # Permute to order x_i, y_i
		  CDF.table[Px$ix,i] < -cumsum(PDF.table[Px$ix,i])
		}
			return(CDF.table)
} **/


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



/**
##############################################################################
# Convert marginal PDF to PDF from data
# Parameters:
# CDF.table - matrix with every column a vector of CDF of each variable
# Output:
# PDF.table - matrix with every column a vector of PDF of each variable
##############################################################################
CDFToPDFMarginals < -function(CDF.table)
{
	n < -dim(CDF.table)[1]  # number of samples
		PDF.table < -array(0L, dim(CDF.table))  # matrix(0, num.samples, num.variables)
		print("DIM CDF -> PDF:")
		print(dim(PDF.table))
		for (i in 1 : dim(CDF.table)[2])  # loop on variables
		{
		  sorted.CDF < -sort(CDF.table[,i], index.return = TRUE)
		  print("sorted.CDF:")
		  print(sorted.CDF)
		  PDF.table[sorted.CDF$ix,i] < -c(sorted.CDF$x[1], sorted.CDF$x[-1] - sorted.CDF$x[-n])
		}
			return(PDF.table)
}
**/


// [[Rcpp::export]]
List EstimateMarginals_rcpp(NumericMatrix data, string w_fun, NumericVector params)  // inputs  //	double* xy[2], double* PDFs[2], double* CDFs[2])  // outputs 
{
	List ret;
	
	long i, j;
	long n = data.nrow();

	NumericVector w_inv(n); //	double* w_inv = new double[n];
	double w_inv_sum = 0.0;

	NumericVector Fx(n);
	NumericVector Fy(n);
	NumericMatrix xy(n, 2);
	NumericMatrix PDFs(n, 2);
	NumericMatrix CDFs(n, 2);

	long naive_flag = FALSE, pos_flag = FALSE;
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
		for (i = 0; i < n; i++)
		{
			w_inv[i] = 1.0 / w_fun_eval(data(i, 0), data(i, 1), w_fun); // NEED TO IMPLEMENT
			w_inv_sum += w_inv[i];
		}
		for (i = 0; i < n; i++)   // normalize
			PDFs(i, 0) = PDFs(i, 1) = w_inv[i] / w_inv_sum;
		CDFs = PDFToCDFMarginals_rcpp(data, PDFs);  // why null ?
		ret["xy"] = data;
	}
	else {
		// skip survival 
		if (naive_flag) // no bias(assume W(x, y) = 1)
		{
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

			for (i = 0; i < 2*n; i++) // loop only up to n? 
			{
				new_data_sorted(i, 0) = new_data(Px[i], 0);
				new_data_sorted(i, 1) = new_data(Py[i], 1);
			}
			NumericVector F0 = empirical_cdf_rcpp(new_data(_, 0));
			NumericVector F1 = empirical_cdf_rcpp(new_data(_, 1));

			long F01, F10;
			for (i = 0; i < 2*n; i++) // loop to n or 2*n? 
			{
				F01 = binary_search_rcpp(new_data_sorted(_, 0), new_data(i, 1));
				F10 = binary_search_rcpp(new_data_sorted(_, 1), new_data(i, 0));
//				F01 = binary_search_rcpp(new_data_sorted(Range(0, n), 0), new_data(i, 1));
//				F10 = binary_search_rcpp(new_data_sorted(Range(0, n), 1), new_data(i, 0));

				CDFs(i, 0) = (F01 + F10) / 2;
				CDFs(i, 1) = (F01 + F10) / 2;  // need to change !!! 
			}
			ret["xy"] = new_data;
		} // end if naive w
		PDFs = CDFToPDFMarginals_rcpp(CDFs); // convert 
	} // end if positive w

	
	ret["PDFs"] = PDFs;
	ret["CDFs"] = CDFs;
	
	return(ret);

}  // end function 






// Next Estimate Marginals: 
// List of needed functions here: 
// PermutationsMCMC <- function(W, N, prms) # burn.in=NA, Cycle=NA)  # New: allow non-default burn-in 
// accumulate - NOT NEEDED !

// Completed (?) 
// (V) empirical_cdf
// (V) w_fun_eval
// (V) binary_search 
// (V) double PDFToCDF2d(double** pdf_2d, double* data[2], long n, double** cdf_2d)
// (V) QuarterProbFromBootstrap <- function(data, null.distribution, grid.points)
// (V) QuarterProbFromPermutations <- function(data, P, grid.points) #Permutations
// (V) CDFToPDFMarginals <- function(CDF.table)
// (V) PDFToCDFMarginals <- function(data, PDF.table)
// (V) EstimateMarginals <- function(data, w.fun, prms=c())  - Tough (complicated) Depends on many functions 




