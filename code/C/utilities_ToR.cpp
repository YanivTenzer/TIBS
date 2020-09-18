#include <iostream>
#include <math.h>
// #include <map>
#include <random>

// #include <Rcpp.h> // for including R with Rcpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
//	cout << "sort perm:" << endl; 
//	for (long i = 0; i < index.length(); i++)
//		cout << index[i] << ", "; 
//	cout << endl;
	return(index);
}

// Gaussian density
double normalPDF(double value)
{
   return exp(-value*value*0.5) /  2.506628; // divide by sqrt(2*pi) 
}


// Gaussian cumulative
double normalCDF(double value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}


/** Copied from here: https://stackoverflow.com/questions/37143283/finding-unique-rows-in-armamat **/
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
/**/


// [[Rcpp::export]]
double mean_rcpp(NumericVector x)
{
	return sum(x)/x.size();
}

// [[Rcpp::export]]
double var_rcpp(NumericVector x, bool bias = true)
{
	// Calculate the mean using C++ function
	long n = x.size();
	double mean = sum(x)/n;
	double sum = 0.0;
	for(int i = 0; i < n; i++)
		sum += pow(x[i] - mean, 2.0); // Square

	return sum/(n-bias); // Return variance
}

NumericVector rowVars_rcpp(NumericMatrix W)
{
	long n = W.nrow();
	NumericVector V(n);
	for(long i=0; i<n; i++)
		V[i] = var_rcpp(W(i,_));
	return V;
}


// [[ Rcpp :: export ()]]
//NumericMatrix a4(arma::mat x) {
//	NumericMatrix y = wrap(x);
//	return(y);
//}


/*
* Copied from here: https://figshare.com/articles/Algorithm_of_binary_search_Rcpp_code/3394741
* A binary search divides a range of values (a sorted array) into halves by the median point(noted k), and continues
 * tonarrow down the field of search until the pattern(noted x) is found.
 * It is the classic example of a "divide and conquer" algorithm.
 */
 // [[Rcpp::export]]
long binary_search_rcpp(NumericVector array, double pattern, long exact_match=0) 
{
	long array_length = array.size();
	long a = 0, b = array_length - 1;
	long mid;
	while (a < b) {
		mid = (a + b) / 2;
		if (array[mid] < pattern) {
			a = mid+1;
		}
		else {
			b = mid;
		}
	}
	if(exact_match && (array[a] != pattern))
		return -1; // didn't find value 
	else
		return a; // return closest index 
}

// Sample from weighted distribution
 // [[Rcpp::export]]
NumericVector weighted_rand_rcpp(long n_samples, NumericVector cdf)
{
	long n = cdf.size();
	NumericVector ret(n);
	double r = double(rand()+1.0) / (RAND_MAX+2.0);
	for(long i=0; i<n; i++)
	{
		ret[i] = binary_search_rcpp(cdf, r, 0);
		if((cdf[ret[i]] < r) || (r == 0))
			Rcout << "!!!! Error sampling: r=" << r << " i=" << ret[i] << " cdf=" << cdf[ret[i]] << " RAND MAX: " << RAND_MAX << endl; 
	}
	return(ret);
}



// Compute empirical cdf 
// [[Rcpp::export]]
NumericVector ecdf_rcpp(NumericVector x)  
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



// Compute one-dimensional empirical cdf
// [[Rcpp::export]]
double ecdf1_rcpp(double x, NumericVector cdf, NumericVector data)
{
	long i, i_max=-1; 
	for(i=0; i<data.length(); i++)
		if(data[i] <= x)
			if((i_max == -1) || (data[i] >= data[i_max]))
				i_max = i;
	if(i_max == -1)
    	return(0.0);
	else
	  return(cdf[i_max]);  
}


// Compute empirical cdf 2d
// [[Rcpp::export]]
double ecdf2_rcpp(NumericVector xy, NumericMatrix cdf_2d, NumericMatrix data)
{
	long i, i_max=-1, j_max=-1;

///	Rcout << "Start loop ecdf2 len=" << data.nrow() << endl; 
	for(i=0; i<data.nrow(); i++)
	{
		if(data(i,0) <= xy[0])
			if((i_max == -1) || (data(i,0) > data(i_max,0)))
				i_max = i;
		if(data(i,1) <= xy[1])
			if((j_max == -1) || (data(i,1) > data(j_max,1)))
				j_max = i;
	}

///	Rcout << "End loop ecdf2 i_max=" << i_max << " j_max=" << j_max << endl; 
///	Rcout << "cdf_2d dims=" << cdf_2d.nrow() << ", " << cdf_2d.ncol() << endl; 
	if( (i_max==-1) || (j_max==-1))
    	return(0.0);
	else
	{
//		Rcout << "INDEX MAX CPP: i,j: " << i_max << ", " << j_max << endl; 
  		return(cdf_2d(i_max,j_max));
	}
}


// Compute matrix empirical cdf 2d from just values
/**
NumericVector ecdf2_rcpp(NumericVector x, NumericVector y)
{
	long i, j;
	long n = x.length();
	IntegerVector Px(n), Py(n); 
	NumericMatrix ecdf2(n, n);

	Px = sort_indexes_rcpp(x); 
	for (i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			ecdf[Px[i]][Py[j]] = double(i+1) / double(n); // use double division here
	return(ecdf2);
}
**/


// Sample a random permutation from the uniform distribution over S_n 
// [[Rcpp::export]]
NumericVector rand_perm(long n)
{
	NumericVector perm(n); 
	long i, temp, r;

	for (i = 0; i < n; i++)
		perm[i] = i; // start with the identity 
	for (i = 0; i < n; i++)
	{
		r = rand() % (n - i) + i;
		temp = perm[i]; // swap 
		perm[i] = perm[r];
		perm[r] = temp; 
	}
	return(perm);
}



// should have a similar function for a vector !! 
// [[Rcpp::export]]
double w_fun_eval_rcpp(double x, double y, string w_fun)
{
	double r = 0.0;

	// allow only string for now 
	// different bias functions (instead of switch)
	if ((w_fun == "truncation") || (w_fun == "Hyperplane_Truncation"))
		r = x <= y; 
	if ((w_fun == "exp") || (w_fun == "exponent_minus_sum_abs") || (w_fun == "stritcly_positive"))
		r = exp((-abs(x) - abs(y)) / 4);
	if (w_fun == "gaussian")
		r = exp((-x*x - y*y) / 2);
	if (w_fun == "huji")
		r = fmax(fmin(65.0 - x - y, 18.0), 0.0);
	if (w_fun == "sum")
		r = x + y; 
	if ((w_fun == "naive") || (w_fun == "const"))
		r = 1;
	if (w_fun == "step")
		r = double(x <= y)+0.1; 

	return(r);
}

// [[Rcpp::export]]
NumericVector w_fun_eval_vec_rcpp(NumericVector x, double y, string w_fun)
{
	long n = x.length();
	long i = 0;
	NumericVector r(n);

	for (i = 0; i < n; i++)
	{
		// allow only string for now 
		// different bias functions (instead of switch)
		if ((w_fun == "truncation") || (w_fun == "Hyperplane_Truncation"))
		{
			r[i] = x[i] <= y; continue;
		}
		if ((w_fun == "exp") || (w_fun == "exponent_minus_sum_abs") || (w_fun == "stritcly_positive"))
		{
			r[i] = exp((-abs(x[i]) - abs(y)) / 4); continue;
		}
		if (w_fun == "gaussian")
			r[i] = exp((-x[i]*x[i] - y*y) / 2);
		if (w_fun == "huji")
		{
			r[i] = fmax(fmin(65.0 - x[i] - y, 18.0), 0.0); continue;
		}
		if (w_fun == "sum")
		{
			r[i] = x[i] + y; continue;
		}
		if (w_fun == "step")
		{
			r[i] = double(x[i] <= y)+0.1; continue;
		}
		if ((w_fun == "naive") || (w_fun == "const"))
		{
			r[i] = 1.0; continue;
		}
	}
	return(r);
}



// [[Rcpp::export]]
NumericMatrix w_fun_to_mat_rcpp(NumericMatrix data, string w_fun)
{
	long n = data.nrow();  // get sample size
	NumericMatrix w_mat(n, n);
	long i; //  , j;
	for (i = 0; i < n; i++)
	{
//		Rcout << "Run VEC!" << endl; 
//		for (j = 0; j < n; j++)
//			w_mat(i, j) = w_fun_eval_rcpp(data(i, 0), data(j, 1), w_fun);
		w_mat(_, i) = w_fun_eval_vec_rcpp(data(_, 0), data(i, 1), w_fun); // fix transform bug 
	}
	return (w_mat);
}


// Determine empirically  w_max for functions where we don't know it (this is approximate and may fail)
/**
double set_w_max_rcpp(double n = 1000, string dependence_type, string w_fun)
{
	w_fun_	
	NumericMatrix xy = SimulateSample_rcpp(n, dependence_type); // need to implement simulate sample 
	return(5 * (max(w_fun_eval_vec_rcpp(xy(_, 1), xy(_, 2), w_fun)) + 5 * std(w_fun_eval_vec_rcpp(xy(_, 1), xy(_, 2), w_fun)) + 1.0));
}
**/


// Set the maximum w value we see in the sample (useful for bootstrap)
// [[Rcpp::export]]
double set_w_max_rcpp_sample(NumericMatrix data, string w_fun)
{
	double epsilon = 0.000001;
	NumericMatrix w_mat = w_fun_to_mat_rcpp(data, w_fun);
	return max(w_mat)*(1+epsilon); 
}



// Determine if weighing function is positive (doesn't support functions)
// [[Rcpp::export]]
long is_pos_w_rcpp(string w_fun, NumericMatrix data, long mat_flag) // , data, mat.flag)
{
	long i; 
	string pos_w[7] = {"sum", "sum_coordinates", "exponent_minus_sum_abs", "const", "naive", "gaussian", "step"}; // all positive w -> should replace by a function checking for positivity 
	string zero_w[2] = {"truncation", "Hyperplane_Truncation"}; // all positive w -> should replace by a function checking for positivity 

	for(long i=0; i<7; i++)
		if(w_fun == pos_w[i])
			return(TRUE);
	for(long i=0; i<2; i++)
		if(w_fun == zero_w[i])
			return(FALSE);

 	if(!mat_flag)  // run and compute values 
	{
		for(i=0; i < data.nrow(); i++)
		    if(w_fun_eval_rcpp(data(i,0), data(i,1), w_fun ) == 0) 
				return(FALSE);
		return(TRUE);
	} else
	{
		return(min(w_fun_to_mat_rcpp(data, w_fun)) > 0); //  test all pairs x_i, y_j 
	}	
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
	NumericVector Rx(n), Ry(n), Eqx(n), Eqy(n);

//	long* Rx = new long[n]; // need to free? 
//	long* Ry = new long[n];
//	long* Eqx = new long[n];
//	long* Eqy = new long[n];


	// New: enable not to give a grid_points as input: set it as data 
	if(grid_points.nrow()<= 1)	
		grid_points = data; 
//	Rcout << "Inside Comp. Stat." << endl; 

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
			Eqx[j] = data(j, 0) == grid_points(i, 0);
			Eqy[j] = data(j, 1) == grid_points(i, 1);
			Obs[0] += double(Rx[j] * Ry[j]);
			Obs[1] += (Rx[j]-Eqx[j]);
			Obs[3] += (Ry[j]-Eqy[j]);
			Obs[2] += (-Eqx[j] - Eqy[j] + Eqx[j] * Eqy[j]);

			if( (Eqx[j]>0) || (Eqy[j]>0) )
				Rcout << "Found data == grid point" << endl; 
		}
		Obs[1] -= Obs[0];
		Obs[3] -= Obs[0];
		Obs[2] += (n - Obs[0] - Obs[1] - Obs[3]);  // here need to check the total number of observed!  why n? 

///		double sum_exp = Exp[0] + Exp[1] + Exp[2] + Exp[3];
///		double sum_obs = Obs[0] + Obs[1] + Obs[2] + Obs[3];
///		Rcout << "SUM-EXP=" << sum_exp << " SUM-OBS=" << sum_obs << endl; 

		if ((Exp[0] > 1) && (Exp[1] > 1) && (Exp[2] > 1) && (Exp[3] > 1))
			for (j = 0; j < 4; j++)
				Statistic += pow((Obs[j] - Exp[j]), 2) / Exp[j];  // set valid statistic when expected is 0 or very small
	} // end loop on grid points

//	Rcout << "Inside Comp. Stat. Finished Loops Stat=" << Statistic << endl; 


//	cout << "Return Stat: " << Statistic << endl;
	return(Statistic);
}


// Weighted statistic - alternative version (works only for positive w)
// [[Rcpp::export]]
double ComputeStatistic_w_rcpp(NumericMatrix data, NumericMatrix grid_points, string w_fun, long counts_flag) 
{
	long n = data.nrow();
	long i, j;
	NumericMatrix w_vec(n); 
	double n_w=0.0;
	double Statistic = 0.0;
	IntegerVector Rx(n);
	IntegerVector Ry(n);
	IntegerVector Lx(n);
	IntegerVector Ly(n);

	// long counts = 0; // 0: use probabilities  for expected, observed / 1: use counts 

	// New: enable not to give a grid_points as input: set it as data 
	if (grid_points.nrow() <= 1)
	{
		Rcout << "SETTING GRID=DATA FOR WEIGHTED STATISTIC" << endl;
		grid_points = data;
	}

	for (i = 0; i < n; i++)
	{
		w_vec[i] = w_fun_eval_rcpp(data(i, 0), data(i, 1), w_fun);
		n_w += 1.0 / w_vec[i];
	}
	double Obs[4] = { 0 };
	double Exp[4] = { 0 };

	double Rx_sum, Ry_sum, Lx_sum, Ly_sum; //  Rx_not_sum, Ry_not_sum;

	for (i = 0; i < n; i++) // Slow loop on grid points
	{
		Rx_sum = Ry_sum = Lx_sum = Ly_sum = 0.0; //  Rx_not_sum = Ry_not_sum = 0.0;
		for(j = 0; j < 4; j++)
			Obs[j] = 0;
		for (j = 0; j < n; j++)  // loop on data points  
		{
			Rx[j] = data(j, 0) > grid_points(i, 0);
			Ry[j] = data(j, 1) > grid_points(i, 1);
			Lx[j] = data(j, 0) < grid_points(i, 0);
			Ly[j] = data(j, 1) < grid_points(i, 1);

///			if( (data(j,0) == grid_points(i,0)) || (data(j,1) == grid_points(i,1)) )
///				Rcout << "found point data == grid!! " << endl;

			Rx_sum += Rx[j] / w_vec[j];
			Ry_sum += Ry[j] / w_vec[j];
			Lx_sum += Lx[j] / w_vec[j];
			Ly_sum += Ly[j] / w_vec[j];


			Obs[0] += (Rx[j] * Ry[j] / (w_vec[j])); //  * n_w));
			Obs[1] += (Rx[j] * Ly[j] / (w_vec[j])); //  * n_w));
			Obs[2] += (Lx[j] * Ly[j] / (w_vec[j])); // * n_w));
			Obs[3] += (Lx[j] * Ry[j] / (w_vec[j])); // * n_w));
		}
		Exp[0] = Rx_sum * Ry_sum / n_w; // (n_w * n_w);
		Exp[1] = Rx_sum * Ly_sum / n_w; // (n_w * n_w);
		Exp[2] = Lx_sum * Ly_sum / n_w; // (n_w * n_w);
		Exp[3] = Lx_sum * Ry_sum / n_w; // (n_w * n_w);

		if(!counts_flag)
			for(j=0; j<4; j++)
			{
				Exp[j] /= n_w;
				Obs[j] /= n_w;
			}


		if ((Exp[0] > 1.0/n) && (Exp[1] > 1.0 / n) && (Exp[2] > 1.0 / n) && (Exp[3] > 1.0 / n))
			for (j = 0; j < 4; j++)
				Statistic += pow((Obs[j] - Exp[j]), 2) / Exp[j];  // set valid statistic when expected is 0 or very small
	}

	return(Statistic);
}


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
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				null_distribution(i, j) = w_mat(1, 1) * pdfs(i, 0) * pdfs(j, 1);
				z += null_distribution(i, j);
			}
	}
	else
	{
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

	// Fix ties 
	for (i = n - 1; i > 0; i--) 
		if (data(Px[i - 1], 0) == data(Px[i], 0))
			for(j = 0; j < n; j++)
				cdf_2d(Px[i - 1], j) = cdf_2d(Px[i], j);
	for (j = n - 1; j > 0; j--)
		if (data(Py[j - 1], 1) == data(Py[j], 1))
			for (i = 0; i < n; i++)
				cdf_2d(i, Py[j - 1]) = cdf_2d(i, Py[j]);


	return(cdf_2d); // don't permute back
}



// double GetQuarterExpectedProb(double Point[2], long QId, double* data[2], double** null_distribution_CDF, long n)
// [[Rcpp::export]]
double GetQuarterExpectedProb_rcpp(NumericVector Point, long QId, NumericMatrix data, NumericMatrix null_distribution_CDF)
{
	// convert below R to C++:
	long idx_x = -1, idx_y = -1, i;
	double s = 0.0;
	long n = data.nrow();

	double epsilon = 0.00000000000001; // should be much smaller than grid perturbation 
	double point_x_minus = Point[0];
	double point_y_minus = Point[1];
	if ((QId == 2) || (QId == 3)) // work with 0-3 coordinates for quadrants
		point_x_minus -= epsilon;
	if ((QId == 1) || (QId == 2))
		point_y_minus -= epsilon; 

	for (i = 0; i < n; i++)
		if (data(i, 0) <= point_x_minus)
			if ((idx_x == -1) || (data(idx_x, 0) < data(i, 0)))  // take minimal value larger than Point[0]
				idx_x = i;
	for (i = 0; i < n; i++)
		if (data(i, 1) <= point_y_minus)
			if ((idx_y == -1) || (data(idx_y, 1) < data(i, 1)))
				idx_y = i;

	double cdf_point = 0.0, cdf_x = 0.0, cdf_y = 0.0;
	if ((idx_x > -1) && (idx_y > -1))
		cdf_point =  null_distribution_CDF(idx_x, idx_y);

	long idx_x_max, idx_y_max;
	if ((QId == 0) || (QId == 1)) 
	{
		if(idx_y > -1)
		{
			idx_x_max = which_max(data(_, 0));
			cdf_y = null_distribution_CDF(idx_x_max, idx_y);
		}
	}
	if ((QId == 0) || (QId == 3)) 
	{
		if(idx_x > -1)
		{
			idx_y_max = which_max(data(_, 1));
			cdf_x = null_distribution_CDF(idx_x, idx_y_max);
		}
	}

	switch (QId) { // different quardants
	case 0: {s = 1 + cdf_point - cdf_x - cdf_y; break; } // "huji" prints "1",
	case 1: {s = cdf_y - cdf_point; break; }
	case 2: {s = cdf_point; break; }
	case 3: {s = cdf_x - cdf_point; break; }
	} // end switch 
	
	return(s);
}

// [[Rcpp::export]]
double GetQuarterExpectedProb_rcpp2(NumericVector Point, long QId, NumericMatrix data, NumericMatrix null_distribution_CDF)
{
	double s = 0.0; 
	double max_x = max(data(_,0)); 
	double max_y = max(data(_,1));
	NumericVector Point_minus(2); //  = Point;
	double epsilon = 0.0000000000001;
	Point_minus[0] = Point[0] - epsilon;
	Point_minus[1] = Point[1] - epsilon;
	NumericVector Point_max(2); //  = Point;

///	Rcout << " Quarter: " << QId << " max_x = " << max_x << " max_y = " << max_y << " Point: " << Point[0] << ", " << Point[1] << endl; 

	switch (QId) { // different quardants
	case 0:	{
		s = 1.0 +  ecdf2_rcpp(Point, null_distribution_CDF, data); ///Rcout << "Did 1st ";
		Point_max[0] = Point[0]; Point_max[1] = max_y; 
		s -= ecdf2_rcpp(Point_max, null_distribution_CDF, data); ///Rcout << "Did 2nd ";
		Point_max[0] = max_x; Point_max[1] = Point[1]; 
		s -= ecdf2_rcpp(Point_max, null_distribution_CDF, data);  break;	}
	case 1: {
		Point_max[0] = max_x; Point_max[1] = Point_minus[1]; 
///		Rcout << "PointMax1: " <<  Point_max[0] << ", " << Point_max[1] << endl; 
		s = ecdf2_rcpp(Point_max, null_distribution_CDF, data);
		Point_max[0] = Point[0]; Point_max[1] = Point_minus[1];
///		Rcout << "PointMax2: " <<  Point_max[0] << ", " << Point_max[1] << endl; 
		s -= ecdf2_rcpp(Point_max, null_distribution_CDF, data); break;}
	case 2: {
		s = ecdf2_rcpp(Point_minus, null_distribution_CDF, data); break;}
	case 3: {
		Point_max[0] = Point_minus[0]; Point_max[1] = max_y;
///		Rcout << "PointMax1: " <<  Point_max[0] << ", " << Point_max[1] << endl; 
		s = ecdf2_rcpp(Point_max, null_distribution_CDF, data);
		Point_max[0] = Point_minus[0]; Point_max[1] = Point[1];
///		Rcout << "PointMax2: " <<  Point_max[0] << ", " << Point_max[1] << endl; 
		s -= ecdf2_rcpp(Point_max, null_distribution_CDF, data);}
	} // end switch 
	return(s);
}


// [[Rcpp::export]]
NumericMatrix QuarterProbFromBootstrap_rcpp(NumericMatrix data, NumericMatrix null_distribution, NumericMatrix grid_points)
{
	long i, j;
	long n = data.nrow();
	long n_grid = grid_points.nrow();
	NumericMatrix null_distribution_CDF(n, n);
	NumericVector cur_grid_points(2);
	NumericMatrix mass_table(n_grid, 4);

///	Rcout << "Running PDF->CDF2d" << endl;
	null_distribution_CDF = PDFToCDF2d_rcpp(null_distribution, data);  // convert PDF to CDF 
///	Rcout << "Did PDF->CDF2d" << endl;

//	cout << "n: " << n << " n_grid: " << n_grid << " null-dist: " << null_distribution.nrow() << ", " << null_distribution.ncol() << endl;
	for (i = 0; i < n_grid; i++)  // loop on grid-points
	{
///		Rcout << "Copy i=" << endl; 
		cur_grid_points = grid_points(i, _); // 		cur_grid_points[1] = grid_points(i,1);
///		Rcout << "Copied i=" << endl; 
		for (j = 0; j < 4; j++)  // compute probabilities (normalize outside the function)
		{
///			Rcout << "Try2 exp (" << i << ", " << j << ")";	
			mass_table(i, j) = GetQuarterExpectedProb_rcpp2(cur_grid_points, j, data, null_distribution_CDF); // try new version 
///			Rcout << " mass2=" << mass_table(i, j); 
		}
///		Rcout << endl;
//		mass_table(i,3) = n - mass_table(i,2) - mass_table(i,1) - mass_table(i,0);  // we can't complete to one 
	}
///	Rcout << "Did loop" << endl; 
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
			a = data(j, 0) <= grid_points(i, 0);
			for (k = 0; k < n; k++)
			{
				b = data(k, 1) <= grid_points(i, 1);
				mass_table(i, 2*a + a*(1-b) + (1-a)*b) += P(j, k);
			}
		}
	}
	return(mass_table);
}


/*
###################################################################################################
# Compute quarter probabilities for the standard bivariate Gaussians 
# with truncation y>x. We assume grid.points satisfy y>x
#
# grid.points - where to evaluate probabilitites 
###################################################################################################*/
// [[Rcpp::export]]
NumericMatrix QuarterProbGaussianAnalytic_rcpp(NumericMatrix grid_points, string w_fun)
{
  	long n = grid_points.nrow();	
  	NumericMatrix mass_table(n, 4); 
	NumericVector x = grid_points(_,0); 
	NumericVector y = grid_points(_,1); 

	for(long i = 0; i < n; i++)
	{
		if((w_fun == "truncation") || (w_fun == "Hyperplane_truncation"))
		{
			mass_table(i,0) =  (1-normalCDF(y[i]))*(1+normalCDF(y[i])-2*normalCDF(y[i]));
			mass_table(i,1) =  pow((normalCDF(x[i])-normalCDF(y[i])), 2);
			mass_table(i,2) =  normalCDF(x[i])*(2*normalCDF(y[i])-normalCDF(x[i]));
			mass_table(i,3) =  2*normalCDF(x[i])*(1-normalCDF(y[i]));
		}
		if((w_fun == "const") || (w_fun == "naive"))
		{
			mass_table(i,0) =  (1-normalCDF(x[i]))*(1-normalCDF(y[i]));
			mass_table(i,1) =  (1-normalCDF(x[i]))*normalCDF(y[i]);
			mass_table(i,2) =  normalCDF(x[i])*normalCDF(y[i]);
			mass_table(i,3) =  normalCDF(x[i])*(1-normalCDF(y[i]));
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

	for (c = 0; c < PDFs.ncol(); c++) // loop on variables 
	{
//		Rcout << "PDF -> CDF COL=" << c << endl; 
		Px = sort_indexes_rcpp(data(_, c)); // get indexes permute to order x_i, y_i

//		cout << "sort perm outside:" << endl;
//		for (i = 0; i < Px.length(); i++)
//			Rcout << Px[i] << ", ";
//		Rcout << endl << "Sorted Indexes 0: " << Px[0] << ", n=" << n << endl; 
		CDFs(Px[0], c) = PDFs(Px[0], c);
		for (i = 1; i < n; i++)
		{
//			Rcout << "Copy Index:" << i << " P[]=" << Px[i-1] << ", " << Px[i] << endl;
			CDFs(Px[i], c) = CDFs(Px[i - 1], c) + PDFs(Px[i], c); //  cumsum(PDFs(Px, c)); // need to change 
		}
//		Rcout << "Computed CDFs" << endl;

		for (i = n - 1; i > 0; i--) // fix ties 
			if (data(Px[i-1], c) == data(Px[i], c))
				CDFs(Px[i-1], c) = CDFs(Px[i], c);
//		Rcout << "Copy Ties" << endl;

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


// Use to keep marginals always sorted by x and y
// [[Rcpp::export]]
List sort_marginals_rcpp(List marginals)
{
	long n = as<NumericMatrix>(marginals["xy"]).nrow();

	IntegerVector perm_inds(n); 
	NumericMatrix sorted_xy(n,2);
  	NumericMatrix sorted_PDFs(n,2);
  	NumericMatrix sorted_CDFs(n,2);
	NumericVector temp(n); // = CDFs_sorted(_, 0);

	// copy values 
	for(long j=0; j<2; j++)
	{
		temp = as<NumericMatrix>(marginals["xy"])(_,j);
		sorted_xy(_,j) = temp.sort();
		temp = as<NumericMatrix>(marginals["CDFs"])(_,j);
		sorted_CDFs(_,j) = temp.sort();
		temp = as<NumericMatrix>(marginals["PDFs"])(_,j);
		perm_inds = sort_indexes_rcpp(temp);
		for(long i=0; i < n; i++)
			sorted_PDFs(i,j) = temp[perm_inds[i]];
	}

  	List marginals_sorted;
	marginals_sorted["xy"] = sorted_xy;
	marginals_sorted["PDFs"] = sorted_PDFs;
	marginals_sorted["CDFs"] = sorted_CDFs;
 	return(marginals_sorted); 
}    


// [[Rcpp::export]]
List EstimateMarginals_rcpp(NumericMatrix data, string w_fun)  
{
	List ret;
	
	long i, j;
	long n = data.nrow();
	NumericVector w_inv(n); 
	double w_inv_sum = 0.0;
	NumericMatrix PDFs(2*n, 2);
	NumericMatrix CDFs(2*n, 2);

//	NumericMatrix CDFs_alt(2 * n, 2);

	long naive_flag = FALSE;
//	string pos_w[4] = {"sum", "sum_coordinates", "exponent_minus_sum_abs", "const"}; // all positive w -> should replace by a function checking for positivity 
	string naive_w[2] = {"naive", "no_bias"};
	long pos_flag = is_pos_w_rcpp(w_fun, data, 1); // take mat_flag = 1

	for (i = 0; i < 2; i++)
		if (w_fun == naive_w[i])
			naive_flag = TRUE;

	if (pos_flag && (!naive_flag)) // w_fun % in % c('sum', 'sum_coordinates', 'exponent_minus_sum_abs')) // for w(x, y) > 0 cases
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
	else { // not positive 
		// skip survival 
		if (naive_flag) // no bias(assume W(x, y) = 1)
		{
			CDFs = CDFs(Range(0, n - 1), _);
			CDFs(_, 0) = ecdf_rcpp(data(_, 0));
			CDFs(_, 1) = ecdf_rcpp(data(_, 1));
			ret["xy"] = data;
		}
		else { // use W , symmetric , exchangable 
				// general biased sampling function(w.fun can be a function not a string) with exchangable distributions
			// case 2: left truncation, use the estimator of Proposition 1 in the paper for exchangable distributions
					// Augment data to include both xand y values for each axis(due to symmetry)
			double augment_data = TRUE;   // new: add xy values
			NumericMatrix new_data(2*n, 2); // largest 
			NumericMatrix new_data_sorted(2*n, 2);
			if(augment_data) // duplicate x and y values
			{
				for (j = 0; j < n; j++)
					for (i = 0; i < 2; i++)
					{
						new_data(j, i) = data(j, 0);  // here 0, 1 are the same 
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

			for (i = 0; i < 2*n; i++) // copy sorted 
			{
				new_data_sorted(i, 0) = new_data(Px[i], 0);
				new_data_sorted(i, 1) = new_data(Py[i], 1);
			}
			NumericVector F0 = ecdf_rcpp(new_data_sorted(_, 0)); // can be used instead of binary search 
			NumericVector F1 = ecdf_rcpp(new_data_sorted(_, 1));
			// alternative: take average of F0 and F1
			CDFs(_, 0) = (F0 + F1) / 2;
			CDFs(_, 1) = CDFs(_, 0); //  (F0 + F1) / 2;

			ret["xy"] = new_data;
		} // end if naive w

		PDFs = CDFToPDFMarginals_rcpp(CDFs); // convert (this part is fast)
	} // end if positive w

	ret["PDFs"] = PDFs;
	ret["CDFs"] = CDFs;	

//	Rcout << " Returning Marginals" << endl;
	return(ret);

}  // end function 


/**
#################################################################
# sample permutations, using MCMC, over the set of valid permutations,
# with respect to the distribution appears in Eq 8
# Parameters:
# w_mat - matrix with weights
# prms: B - number of permutations to draw
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
	long burn_in = 0; //  2 * n;
	if (prms.containsElementNamed("burn_in")) //  ('burn.in' % in % names(prms)))
		burn_in = prms["burn_in"];
	long Cycle = n; 
	if (prms.containsElementNamed("Cycle")) // (!('Cycle' % in % names(prms)))
		Cycle = prms["Cycle"];
	long Idx = 0, ctr = 1;

//	Rcout << "burn.in=" << burn_in << " Cycle=" << Cycle << endl; 

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
		i = switchIdx[0];
		j = switchIdx[1];
		ratio = w_mat(i, Perm[j]) * w_mat(j, Perm[i]) / (w_mat(i, Perm[i]) * w_mat(j, Perm[j]));

		if(	double(rand()) / RAND_MAX < fmin(1.0, ratio) ) // we accept the transition with probability min(1, ratio)
		{
			temp = Perm[i];  // SWAP
			Perm[i] = Perm[j];
			Perm[j] = temp;
			if (ctr == burn_in || ((ctr % Cycle == 0) && (ctr > burn_in)))
			{
				for(i=0; i<n; i++)
					Permutations(i,Idx) = Perm[i] + 1;  // add 1 to move from C zero-based indexing to R-one based indexing. A problem if we use it INSIDE C!!!! 
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


/**
#################################################################
# sample permutations, using Importance Sampling, over the set of valid permutations,
# Parameters:
# w_mat - matrix with weights
# prms, with: B - number of permutations to draw
# 
# Output:
# PermutationsTable - A matrix representing the sampled permutations
# P - An n * n matrix with P(i, j) = Pr(pi(i) = j) over the sampled permutations
#########################################################################################
**/

// [[Rcpp::export]]
List PermutationsIS_rcpp(NumericMatrix w_mat, List prms) // burn.in = NA, Cycle = NA)  # New: allow non - default burn - in
{
	long n = w_mat.nrow();
	long i, j, k, b;

	// Set IS default sampling parameters
	long B = 1000;
	if (prms.containsElementNamed("B")) //   ('B' % in % names(prms)))
		B = prms["B"];
	string importance_sampling_dist = "uniform";
    if(prms.containsElementNamed("importance.sampling.dist"))  // set default uniform distribution 
      importance_sampling_dist = as<string>(prms["importance.sampling.dist"]);

	NumericMatrix Permutations(n, long(B));

//	Rcout << "burn.in=" << burn_in << " Cycle=" << Cycle << endl; 

//	IntegerVector Perm = seq(0, n);  // 1 to N-1 
	NumericVector log_P_IS(B); // probabilities for importance sampling (up to a constant)
	double log_P_IS0 = 0.0;
	
	NumericMatrix log_w_mat(n,n); // = log(w_mat); // (n, n);
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			log_w_mat(i,j) = log(w_mat(i,j));

	// choose IS distribution: 
	if(importance_sampling_dist == "uniform")
	{
		for(b=0; b<B; b++)
			Permutations(_,b) = rand_perm(n);
		for(i=0; i<B; i++)
			log_P_IS[i] = 0.0; // save probabilities on a log-scale
		log_P_IS0 = 0.0;
	}
	if(importance_sampling_dist == "monotone.w") // sort values by average of w
	{
		NumericVector w_mat_row_sums = rowSums(w_mat); 
		IntegerVector w_order = sort_indexes_rcpp(w_mat_row_sums);	// sort in increasing order! (i.e. start with the small ones)

		NumericVector weights(n), weights_cdf(n);
		for(j=0; j<n; j++)
		{				
			weights = w_mat(w_order[j],_);
			for(k=0; k<j; k++)
                weights[w_order[k]] = 0.0;     
//			weights = weights * (1.0/sum(weights)); // normalize
			log_P_IS0 += log(weights[w_order[j]] / sum(weights) );
		}
		for(b=0; b<B; b++)
		{
			for(j=0; j<n; j++)
			{
				weights = w_mat(w_order[j],_);
				for(k=0; k<j; k++)
                   weights[Permutations(w_order[k], b)] = 0.0;     
				weights = weights * (1.0/sum(weights)); // normalize
				weights_cdf[0] = weights[0]; // compute cumsum
				for(k=1; k<n; k++)
					weights_cdf[k] = weights_cdf[k-1] + weights[k]; 	
				Permutations(w_order[j],b) = weighted_rand_rcpp(1, weights_cdf)[0]; // IS this fucked up? 
				if(weights[Permutations(w_order[j],b)] == 0)
					Rcout << "Error! Chose Perm with PRob=0 !!! " << "Weights: " << weights_cdf << " ind: " << Permutations(w_order[j],b) << endl; 
				log_P_IS[b] += log(weights[Permutations(w_order[j],b)]); // Need to update also log_P_IS here
			}
		}
	}

	if(importance_sampling_dist == "match.w") // new: try to match P_W of the original permutation 
	{
		NumericVector w_mat_row_vars = rowVars_rcpp(log_w_mat); // Should be rowVars!!! order variables by variance of rho 
		IntegerVector w_order = sort_indexes_rcpp(-w_mat_row_vars);  // sort in DECREASING order!

		NumericVector log_w_sum_id_vec(n); 
		log_w_sum_id_vec[0] = log_w_mat(w_order[0],w_order[0]); log_P_IS0 = -log(1);
		for(j=1; j<n; j++)
		{
			log_w_sum_id_vec[j] = log_w_sum_id_vec[j-1] + log_w_mat(w_order[j],w_order[j]); 	//	  cumsum(diag(log.w.mat)[w.order]) # we want to get this sum 
		//	log_P_IS0 += (-log(j+1)); // TEMP! set as uniform distirbution 
		}

		double log_w_sum_rand = 0.0;
		double sigma = 2.0; // set variance of importance sampling around true permutation 
		NumericVector weights(n), weights_cdf(n);
		for(j=0; j<n; j++) // update also the log-prob for the identity permutation 
		{
			for(k=0; k<n; k++)
	            weights[k] = exp(-pow(log_w_mat(w_order[j],k) + log_w_sum_rand - log_w_sum_id_vec[j], 2 ) / (2*sigma*sigma)); //  penalize deviations from probability
			for(k=0; k<j; k++) // remove the ones we already occupied
	            weights[w_order[k]] = 0.0;     	
            weights = weights * (1.0/ sum(weights)); 
			log_P_IS0 += log(weights[w_order[j]]); // Need to update also log_P_IS here
            log_w_sum_rand += log_w_mat(w_order[j],w_order[j]);
		}
		for(b=0; b<B; b++)
		{
			log_w_sum_rand = 0.0;
			for(j=0; j<n; j++)
			{
				for(k=0; k<n; k++)
	                weights[k] = exp(-pow(log_w_mat(w_order[j],k) + log_w_sum_rand - log_w_sum_id_vec[j], 2 ) / (2*sigma*sigma)); //  penalize deviations from probability
				for(k=0; k<j; k++) // remove the ones we already occupied
	                weights[w_order[k]] = 0.0;     	
                weights = weights * (1.0/ sum(weights)); 
				weights_cdf[0] = weights[0]; // compute cumsum
				for(k=1; k<n; k++)
					weights_cdf[k] = weights_cdf[k-1] + weights[k]; 
                Permutations(w_order[j],b) = weighted_rand_rcpp(1, weights_cdf)[0];
                log_P_IS[b] += log(weights[Permutations(w_order[j],b)]); // Need to update also log_P_IS here
                log_w_sum_rand += log_w_mat(w_order[j],Permutations(w_order[j],b));
			}
		}

	}


//	Rcout << "log_P_IS, log_P_IS0: " << log_P_IS << endl << log_P_IS0 << endl; 

	// Compute P matrix 

	double P_W_IS0 = -log_P_IS0; // log_P_IS0 on a log-scale 
	for(i=0; i<n; i++)
		P_W_IS0 += log_w_mat(i,i);
	NumericVector P_W_IS = -log_P_IS;  // log_P_IS on a log scale 
	for(b=0; b<B; b++)
		for(i=0; i<n; i++)
			P_W_IS[b] += log_w_mat(i, Permutations(i,b));
//	Rcout << "P_W_IS_LOG: " << P_W_IS << endl; 

	double max_log = max(P_W_IS);
	max_log = max(max_log, P_W_IS0);

//	Rcout << "max log is: " << max_log << endl; 

	P_W_IS0 = exp(P_W_IS0 - max_log);
	for(b=0; b<B; b++)
		P_W_IS[b] = exp(P_W_IS[b]-max_log); // compute P_W/P_I (up to a constant)
//	Rcout << "P_W_IS: " << P_W_IS << endl; 

	// Compute also P_W on log-scale 
 	double log_P_W0 = 0;
	 for(i=0; i<n; i++)
	 	log_P_W0 += log_w_mat(i, i);
    NumericVector log_P_W(B); 
    for(b=0; b<B; b++)
	 	for(i=0; i<n; i++)
    		log_P_W[b] += log_w_mat(i, Permutations(i,b));
//  	double max_w = max(max(log_P_W), log_P_W0);


	// New: compute coefficient of variation: 
	double CV2 = var_rcpp(P_W_IS) / pow(mean_rcpp(P_W_IS), 2); 

	NumericMatrix P(n, n);  // New!matrix with P[i] = j estimate
    for(b=0; b<B; b++)  // next, compute expectations P[i,j]
		for(i=0; i<n; i++)
			P(i, Permutations(i,b)) += P_W_IS[b];		
    P = P * (0*P_W_IS0 + 1.0/sum(P_W_IS)); // normalize   (ignore contribution from the identity permuation)

	NumericVector P_sum = rowSums(P);

/**
	Rcout << "P[i,j] matrix:" << endl;
	for(i=0; i<5; i++)
	{
		for(j=0; j<5; j++)
			Rcout << P(i,j) << " ";
		Rcout << endl;
	}
	Rcout << "Sum rows: ";
	for(i=0; i<5; i++)
		Rcout << P_sum[i] << endl;
	P_sum = colSums(P);
	Rcout << "Sum cols: ";
		for(i=0; i<5; i++)
			Rcout << P_sum[i] << endl;
**/
	List ret; 
	ret["Permutations"] = Permutations;
	ret["P"] = P;
	ret["P.W.IS"] = P_W_IS;
	ret["P.W.IS0"] = P_W_IS0;
	ret["log.P.W"] = log_P_W;
	ret["log.P.W0"] = log_P_W0;
	ret["log.P.IS"] = log_P_IS;
	ret["log.P.IS0"] = log_P_IS0;
	ret["CV2"] = CV2;

//	log.P.W=log.P.W, log.P.W0=log.P.W0, max.w=max.w, 
//              log.P.IS=log.P.IS, log.P.IS0=log.P.IS0)
	return(ret); // return list with also P
}


/**
#############################################################
# Calculate p - value using importance sampling
# (If product of weights is very small/large, we can
	#  multiply the weight function by a constant so the weights
	#  are in some sense centered at 1)
#############################################################
**/
// [[Rcpp::export]]
List IS_permute_rcpp(NumericMatrix data, NumericMatrix grid_points, string w_fun, List prms, string test_type) // NumericMatrix expectations_table, long counts_flag) 
{
	long  n = data.nrow();
	long B = prms["B"];
	long counts_flag = TRUE; 
	if (prms.containsElementNamed("counts.flag"))
		counts_flag = prms["counts.flag"];

	List PermutationsList = PermutationsIS_rcpp(w_fun_to_mat_rcpp(data, w_fun), prms); //  sample permutations 	
	NumericMatrix Permutations = as<NumericMatrix>(PermutationsList["Permutations"]);	
	NumericMatrix P = as<NumericMatrix>(PermutationsList["P"]);	
	NumericVector P_W_IS = as<NumericVector>(PermutationsList["P.W.IS"]); // importance weights of permutations
	double P_W_IS0 = as<double>(PermutationsList["P.W.IS0"]); // importance weight of ID permutation
	NumericMatrix expectations_table(grid_points.nrow(), 4);

	double TrueT; // = ComputeStatistic_w_rcpp(data, data, w_fun); 
	long j, b;
	NumericVector pw(B); 
	NumericVector Tb(B); 
	IntegerVector perm(n);
	NumericMatrix permuted_data(n, 2);

	if(grid_points.ncol() == 1) // no input grid
		grid_points = data;

	long inverse_weight = test_type.find("inverse_weight");
//	Rcout << "Start IS Permute Inverse Weight: " << inverse_weight << " Dim(Exp-Table)=" << expectations_table.nrow() << ", " << expectations_table.ncol() << endl; 
  	if(inverse_weight)
    	TrueT = ComputeStatistic_w_rcpp(data, grid_points, w_fun, counts_flag); // weights. no unique in grid-points 
	else
	{
		NumericMatrix expectations_table =  QuarterProbFromPermutations_rcpp(data, P, grid_points); 
    	TrueT = ComputeStatistic_rcpp(data, grid_points, expectations_table); // no weights
	}

	for (b = 0; b < B; b++) 
	{
		perm = Permutations(_,b); //		perm = rand_perm(n);  // get a random permutation from the uniform disitribution
		permuted_data(_, 0) = data(_, 0);
		for (j = 0; j < n; j++)
			permuted_data(j, 1) = data(perm[j], 1); // permute data 
		if(inverse_weight)
			Tb[b] = ComputeStatistic_w_rcpp(permuted_data, grid_points, w_fun, counts_flag); // grid depends on permuted data. Compute weighted statistic! 
		else 
			Tb[b] = ComputeStatistic_rcpp(permuted_data, grid_points, expectations_table); // grid depends on permuted data. Compute weighted statistic! 
	}

	double Pvalue = P_W_IS0; //1; // P_W_IS0;  weight identity permutation by 1 
	double NormalizingFactor = P_W_IS0; // 1; // P_W_IS0;
	for(b=0; b<B; b++)
	{
		Pvalue += (Tb[b]>=TrueT) * P_W_IS[b];	
		NormalizingFactor += P_W_IS[b];
	}
	Pvalue /= NormalizingFactor;

	double CV2 = var_rcpp(P_W_IS) / pow(mean_rcpp(P_W_IS), 2); // new: add coefficient of variation for importance sampling diagnostics. 

//	Rcout << "reject=" << reject << " denominator.weights=" <<  (pw0 + sum(pw)) << " pval=" << reject / (pw0 + sum(pw)) << endl; 
	List ret; 
	ret["Pvalue"] = Pvalue;  // New: Take into account also the idenity permutation in the pvalue!! 
	ret["TrueT"] = TrueT; 
	ret["statistics.under.null"] = Tb;  // get null statistics (they're not weighted equally)
	ret["perm.weights"] = P_W_IS; // New! return also weight for each permutation !! 
//	ret["pw.max"] = pw_max;
	ret["id.perm.weight"] = P_W_IS0;
	ret["CV2"] = CV2;
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


/**
############################################################################################
# Draw a bootstrap sample, given the estimated null distribution Fx, Fy, W
# Use rejection sampling.No need for computing n* n table
# (Problem: what if prob.(rejection) close to 1?)
# Parameters:
# data - n * 2 array with(X, Y) samples
# cdfs - fx and fy - cumulatives (not pdfs!) SORTED!
# w.fun - Biased sampling function w
# prms - for w max
# n - allow a different sample size
# Output: 
# an n*2 matrix of bootstrap sample 
############################################################################################**/
// [[Rcpp::export]]
List Bootstrap_rcpp(NumericMatrix data, NumericMatrix cdfs, NumericMatrix w_mat,  List prms, long n=-1)  // n is optional .
{
	if (n == -1) // is.null(n))
		n = data.nrow();
	IntegerMatrix boot_indices(n, 2);
	NumericMatrix boot_sample(n, 2);

//	double x = 0.0, y = 0.0; //  , r; //  keep;
	long k = 0, ctr=0; 
	long i=0, j=0;
	double w_max = prms["w.max"];


//	Rcout << "Start Bootstrap inside, w_max=" << w_max << endl; 
//	long m;
	long a, b, t, mid=(n-1)/2;
	double r;
	while (k < n)
	{	
		for(t=0; t<2; t++)
		{
			r = double(rand()) / RAND_MAX;
			a = 0;  b = cdfs.nrow() - 1;  // do binary search 
			while (a < b) 
			{
				mid = (a + b) / 2;

///				Rcout << "Bin search: (a, mid, b) = " << a << ", " << mid << ", " << b << ")" << endl; 
				if (cdfs(mid,t) < r) 
					a = mid+1;
			    else 
					b = mid;
			}
			if(t == 0)
				i = a;
			else
				j = a;
///			Rcout << "r=" << r << " a=" << a << " cdf=" << cdfs(a,t) << endl; 
		}
//		x = data(i, 0); y = data(j, 1); //		rand() data(d0(gen), 0); //  rand() % (n - k); // sample with replacement 		 //data(d1(gen), 1); //  rand() % (n - k);
//		if ((double(rand()) / RAND_MAX) < (w_fun_eval_rcpp(x, y, w_fun) / w_max)) // compare to probability (rejection - can be costly when w_max is large!!)
//		Rcout << "Try k=" << k << " i=" << i << " j=" << j << " w(x,y)=" << w_mat(i, j) << endl; 
		if ((double(rand()) / RAND_MAX) < (w_mat(i, j) / w_max)) // compare to probability (rejection - can be costly when w_max is large!!)
		{
			boot_indices(k, 0) = i; boot_indices(k, 1) = j;
			boot_sample(k, 0) = data(i, 0); boot_sample(k++, 1) = data(j, 1);
		}
		if((i < 0) || (i >=cdfs.nrow()) || (j < 0) || (j>=cdfs.nrow()))
			Rcout << "ERROR! i,j: " << i << ", " << j << " k=" << k << " n=" << n << " data.n=" << data.nrow() << " cdf.n=" << cdfs.nrow() << " w.mat.n=" << w_mat.nrow() << endl; 
		ctr++;
	}
//	Rcout << "Generated " << ctr << " pairs to get eventually " << n << " samples" << endl; 
	List ret;
	ret["sample"] = boot_sample; 
	ret["indices"] = boot_indices;

	return(ret);
}


// A special function for copying bootstrap data 
List BootstrapOrganize_rcpp(List data_marginals, List bootstrap, string w_fun, List prms)
{
	long i, j;
	long n = prms["sample.size"];
	long marginals_n = n; 
	if((w_fun == "truncation") || (w_fun == "Hyperplane_Truncation"))
		marginals_n *= 2; 

	NumericMatrix bootstrap_indices = as<NumericMatrix>(bootstrap["indices"]); 
	List marginals_bootstrap = EstimateMarginals_rcpp(as<NumericMatrix>(bootstrap["sample"]), w_fun);

	NumericMatrix marginals_bootstrap_PDFs = marginals_bootstrap["PDFs"];
  	List marginals_bootstrap_organized; 
	marginals_bootstrap_organized["xy"] = as<NumericMatrix>(data_marginals["xy"]);
	NumericMatrix marginals_bootstrap_organized_PDFs(marginals_n, 2);
	if((w_fun == "truncation") || (w_fun == "Hyperplane_Truncation")) // here data was duplicated 
	{
		for(i=0; i<n; i++)
			for(j=0; j<2; j++)
			{
				marginals_bootstrap_organized_PDFs(bootstrap_indices(i, j), 0) += marginals_bootstrap_PDFs(i+j*n, j); 
				marginals_bootstrap_organized_PDFs(bootstrap_indices(i, j), 1) += marginals_bootstrap_PDFs(i+j*n, j); 
			}
	} else
	{
		for(i=0; i<n; i++)
			for(j=0; j<2; j++)
				marginals_bootstrap_organized_PDFs(bootstrap_indices(i, j), j) += marginals_bootstrap_PDFs(i, j); 
	}
	/**/
	marginals_bootstrap_organized["PDFs"] = marginals_bootstrap_organized_PDFs;

	return(marginals_bootstrap_organized);
}


/**########################################################################
# 
# All steps of TIBS test :
# 1. estimate marginals
# 2. compute null distribution
# 3. compute expectation table
# 4. compute statistic
# 
########################################################################**/
// [[Rcpp::export]]
List TIBS_steps_rcpp(NumericMatrix data, string w_fun, NumericMatrix w_mat, 
	NumericMatrix grid_points, NumericMatrix expectations_table, List prms)
{
	List marginals;
	List null_distribution;
	string use_w;
	long n = data.nrow(); 
//	Rcout << "Start TIBS STEPS" << endl; 
	long naive_expectation = 0;
	if (prms.containsElementNamed("naive.expectation"))
		naive_expectation = prms["naive.expectation"];

	if (naive_expectation) // here we ignore W(using statistic for unbiased sampling)
	{
		use_w = "naive";
		w_mat(1, 1) = 1;
	}
	else
		use_w = w_fun;
		
	if(expectations_table.nrow()==1) // insert empty table // 	if (missing(expectations_table) | isempty(expectations_table)) // modify to c++ code!
	{
		marginals = EstimateMarginals_rcpp(data, use_w);
		if ((!naive_expectation) & (w_mat.nrow() == 1))
			w_mat = w_fun_to_mat_rcpp(marginals["xy"], w_fun);  // compute W again for augmented data
		null_distribution = GetNullDistribution_rcpp(marginals["PDFs"], w_mat);
		expectations_table = double(n) * QuarterProbFromBootstrap_rcpp(marginals["xy"], null_distribution["distribution"], grid_points);
	}
	double T = ComputeStatistic_rcpp(data, grid_points, expectations_table);  // keep same grid points for bootstrap sample ?

	List ret; ret["Statistic"] = T; ret["expectations_table"] = expectations_table; 
	ret["marginals"] = marginals; ret["w_mat"] = w_mat; ret["null_distirbution"] = null_distribution; 
	return(ret);
}



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
**/
// [[Rcpp::export]]
List TIBS_rcpp(NumericMatrix data, string w_fun, string test_type, List prms)
{
	double epsilon = 0.00000000001;
	long ctr;
	double TrueT, NullT;
	List output; // contains all output information we need 
//	Rcout << "Start TIBS " << endl; 

	// Set defaults
	long fast_bootstrap = FALSE;
	if (prms.containsElementNamed("fast.bootstrap"))
		fast_bootstrap = prms["fast.bootstrap"];
//	if (!prms.containsElementNamed("minp_eps"))
//		prms["minp_eps"] = NULL; // permDep params. Irrelevant here 
	long PL_expectation = FALSE;
	if (prms.containsElementNamed("PL.expectation"))
		PL_expectation = prms["PL.expectation"];
	long naive_expectation = 0;
	if (prms.containsElementNamed("naive.expectation"))
		naive_expectation = prms["naive.expectation"];
	long new_bootstrap = TRUE; 
	if (prms.containsElementNamed("new.bootstrap"))
		new_bootstrap = prms["new.bootstrap"];
	long perturb_grid = TRUE; // add a small perturbation to avoid ties (default: TRUE)
	if (prms.containsElementNamed("perturb.grid"))
		perturb_grid = prms["perturb.grid"];
	long counts_flag = TRUE; // use counts and not probabilities in inverse-weighting chi-square statsitic (default: TRUE)
	if (prms.containsElementNamed("counts.flag"))
		counts_flag = prms["counts.flag"];


//	if (!prms.containsElementNamed("delta"))
//		prms["delta"] = NA; minP2 params. Irrelevant here

	long B = prms["B"];
	long n = data.nrow();
	long i, j; 

//	Rcout << "Start TIBS n= " << n << endl; 
 	NumericMatrix new_data(n, 2); // copy data 
	 NumericMatrix grid_points(n, 2); // copy data 
  	for (i = 0; i < n; i++) // long(grid_points_arma.n_rows); i++)
    	for (j = 0; j < 2; j++)
		{
        	new_data(i, j) = data(i, j);
			grid_points(i, j) = data(i, j);
		}
	
//	Rcout << " Read Input TIBS_RCPP. TEST-TYPE: " << test_type << endl; 
	// Create a grid of points, based on the data :
//	NumericMatrix grid_points = new_data;
	if (prms.containsElementNamed("grid.points"))		
	{
		perturb_grid = FALSE;
		grid_points = as<NumericMatrix>(prms["grid.points"]); // copy input grid ! 
	}

	if(perturb_grid)
	{
		NumericMatrix perturb(n, 2);
		perturb(_, 0) = epsilon * ( max(grid_points(_,0))-min(grid_points(_,0)) ) * rnorm(n);
		perturb(_, 1) = epsilon * ( max(grid_points(_,1))-min(grid_points(_,1)) ) * rnorm(n);
		grid_points += perturb; 
	}
/*	arma::mat grid_points_arma = unique_rows(as<arma::mat>(data));  // set unique for ties ? for discrete data
	NumericMatrix grid_points(n, 2);
	for (i = 0; i < n; i++) // long(grid_points_arma.n_rows); i++)
		for (j = 0; j < 2; j++)
			grid_points(i, j) = grid_points_arma(i, j); **/
//	= as<NumericMatrix>(wrap(grid_points_arma));  // set unique for ties ? for discrete data
	List null_distribution;


//	Rcout << "Set TIBS grid points " << endl; 
	// no switch (test_type) for strings in cpp
	if(test_type == "bootstrap_inverse_weighting") 
	{
		TrueT = ComputeStatistic_w_rcpp(data, grid_points, w_fun, counts_flag); //  compute true statistic before data is sorted 
		List marginals = EstimateMarginals_rcpp(data, w_fun);

		NumericMatrix CDFs_sorted = as<NumericMatrix>(marginals["CDFs"]); // get sorted CDFs 
		NumericVector CDF_temp = CDFs_sorted(_, 0);
		CDFs_sorted(_, 0) = CDF_temp.sort(); //  marginals["CDFs"];
		CDF_temp = CDFs_sorted(_, 1);
		CDFs_sorted(_, 1) = CDF_temp.sort(); //  marginals["CDFs"];
		NumericMatrix xy_sorted = as<NumericMatrix>(marginals["xy"]); // get sorted CDFs 
		NumericVector xy_temp = xy_sorted(_, 0);
		xy_sorted(_, 0) = xy_temp.sort(); //  marginals["CDFs"];
		xy_temp = xy_sorted(_, 1);
		xy_sorted(_, 1) = xy_temp.sort(); //  marginals["CDFs"];

		NumericMatrix w_mat = w_fun_to_mat_rcpp(xy_sorted, w_fun);
//		null_distribution = GetNullDistribution_rcpp(marginals["PDFs"], w_mat);
		NumericVector statistics_under_null(B);
		List bootstrap_sample; //		NumericMatrix bootstrap_sample(n, 2);
		for (ctr = 0; ctr < B; ctr++)
		{
			bootstrap_sample = Bootstrap_rcpp(xy_sorted, CDFs_sorted, w_mat,  prms, n); // draw new sample.Problem: which pdf and data ?
			statistics_under_null[ctr] = ComputeStatistic_w_rcpp(bootstrap_sample["sample"], grid_points, w_fun, counts_flag); // should sample grid points again! 
		}
		output["TrueT"] = TrueT;
		output["statistics.under.null"] = statistics_under_null;
	}

	if(test_type == "bootstrap") 
	{
		NumericMatrix w_mat(1, 1);
		NumericMatrix expectations_table(1, 1);
  
		List TrueTList = TIBS_steps_rcpp(data, w_fun, w_mat, grid_points, expectations_table, prms); // new: replace multiple steps by one function
		TrueT = TrueTList["Statistic"];
//		NumericMatrix true_expectations_table = TrueTList["expectations_table"];

	    // sort data: :
    	NumericVector data_temp = data(_, 0);
    	data(_, 0) = data_temp.sort(); //  marginals["CDFs"];
    	data_temp = data(_, 1);
    	data(_, 1) = data_temp.sort(); //  marginals["CDFs"];        
//		Rcout << "Computed TrueT TIBS RCPP Bootstrap" << endl;
    // Run tibs again
		TrueTList = TIBS_steps_rcpp(data, w_fun, w_mat, grid_points, expectations_table, prms); // new: replace multiple steps by one function
		List marginals = TrueTList["marginals"]; // marginals for true data (true data needs to be sorted for this)
		List bootstrap_sample; //  (n, 2);

		// 2. Compute statistic for bootstrap sample :
		NumericVector statistics_under_null(B);
//		List null_distribution_bootstrap = null_distribution;
//		NumericMatrix w_mat_bootstrap(1,1);
		IntegerMatrix bootstrap_indices(n, 2);
		List null_distribution_bootstrap_new;
		NumericMatrix expectations_table_new(n, 4);
		List marginals_sorted;

		NumericMatrix CDFs_sorted = as<NumericMatrix>(marginals["CDFs"]); // get sorted CDFs 
		NumericVector CDF_temp = CDFs_sorted(_, 0);
		CDFs_sorted(_, 0) = CDF_temp.sort(); //  marginals["CDFs"];
		CDF_temp = CDFs_sorted(_, 1);
		CDFs_sorted(_, 1) = CDF_temp.sort(); //  marginals["CDFs"];
		NumericMatrix xy_sorted = as<NumericMatrix>(marginals["xy"]); // get sorted CDFs 
		NumericVector xy_temp = xy_sorted(_, 0);
		xy_sorted(_, 0) = xy_temp.sort(); //  marginals["CDFs"];
		xy_temp = xy_sorted(_, 1);
		xy_sorted(_, 1) = xy_temp.sort(); //  marginals["CDFs"];
		marginals_sorted["xy"] = xy_sorted;
		marginals_sorted["CDFs"] = CDFs_sorted; // we don't need pdfs 	

		List marginals_bootstrap_new;
		w_mat = w_fun_to_mat_rcpp(xy_sorted, w_fun); // compute matrix once (for sorted data)
		if (!(prms.containsElementNamed("w.max")))
			prms["w.max"] = set_w_max_rcpp_sample(data, w_fun); // set w.max for bootstrap sampling 
		else
			prms["w.max"] = max(as<double>(prms["w.max"]), set_w_max_rcpp_sample(data, w_fun)); // set w.max for bootstrap sampling 
		
		for (ctr = 0; ctr < B; ctr++) // heavy loop : run on bootstrap
			{
//				Rcout << "Run Bootstrap B = " << ctr << endl;

				// statistics_under_null[ctr] = TIBS_steps_rcpp(bootstrap_sample, w_fun, w_mat, grid_points, expectations_table, prms); // replace multiple steps by one function
	////			if(ctr==0) // TEMP DEBUG CRUSH !!! 
					bootstrap_sample = Bootstrap_rcpp(xy_sorted, CDFs_sorted /*marginals["CDFs"]*/, w_mat , prms, n); // draw new sample. Problem: which pdf and data ?
				if (!fast_bootstrap) // re - estimate marginals for null expectation for each bootstrap sample
				{
				  if (!new_bootstrap)
					{
				    List NullT = TIBS_steps_rcpp(bootstrap_sample["sample"], w_fun, NumericMatrix(1, 1), grid_points, NumericMatrix(1, 1), prms);
				    List NullT_marginals = NullT["marginals"];
						statistics_under_null[ctr] = NullT["Statistic"];
					} else
					{
//						Rcout << "Organize Marginals B = " << ctr << endl;  
						marginals_bootstrap_new = BootstrapOrganize_rcpp(marginals_sorted, bootstrap_sample, w_fun, prms); // here marginals AREN'T SORTED !!! 
//						Rcout << "Get Null Distribution B = " << ctr << " xy : PDFs : CDFs " << endl;

						null_distribution_bootstrap_new = GetNullDistribution_rcpp(marginals_bootstrap_new["PDFs"], w_mat); /// TrueTList["w_mat"]); // keep w_mat of ORIGINAL DATA! (NOT SORTED!)
//						Rcout << "Compute Expected B = " << ctr << " xy : PDFs : CDFs " << endl;
						expectations_table_new = double(n) * QuarterProbFromBootstrap_rcpp(
							marginals_bootstrap_new["xy"], null_distribution_bootstrap_new["distribution"], grid_points);
//						Rcout << "Compute Statistic Bootstrap B = " << ctr << " nrows= " << grid_points.nrow() << ", " << expectations_table_new.nrow() << endl;
						statistics_under_null[ctr] = ComputeStatistic_rcpp(bootstrap_sample["sample"], grid_points, expectations_table_new); //  NEW!Compute null statistic without recomputing the entire matrix 
//						Rcout << "Finished Statistic Bootstrap B = " << ctr << " Out of " << B << endl;
					}
				} else // if fast bootstrap
		//		Rcout << "Compute Bootstrap Statistic Under Null " << ctr << endl;
					statistics_under_null[ctr] = ComputeStatistic_rcpp(bootstrap_sample["sample"], grid_points, expectations_table);
			} // counter for bootstrap iterations
		
//		Rcout << "Finished all bootstraps " << B << endl;	
//		Rcout << "End True T " << TrueT << endl;	
		output["TrueT"] = TrueT;
//		Rcout << "Now copy null stats " << endl;	
//		Rcout << "Stat. null " << statistics_under_null << endl;
		output["statistics.under.null"] = statistics_under_null;
//		Rcout << "Ended Bootstrap " << endl;
	} // end if bootstrap 
	
	if(test_type == "permutations")
	{
//		Rcout << "Start permutations w_fun to mat" << endl; 
		NumericMatrix w_mat = w_fun_to_mat_rcpp(data, w_fun);
//		Rcout << "Start permutations MCMC" << endl; 
		List PermutationsList = PermutationsMCMC_rcpp(w_mat, prms);
		NumericMatrix expectations_table(1,1);
					
		NumericMatrix P = PermutationsList["P"];
		NumericMatrix Permutations = PermutationsList["Permutations"];
		NumericMatrix permuted_data(n, 2);

//		Rcout << "Start permutations QuarterPRob" << endl; 
		if(!(PL_expectation || naive_expectation))
			expectations_table = QuarterProbFromPermutations_rcpp(data, P, grid_points);
//		Rcout << "Start permutations TIBS-Steps" << endl; 
		List TrueTList = TIBS_steps_rcpp(data, w_fun, w_mat, grid_points, expectations_table, prms); // new: replace multiple steps by one function
		TrueT = TrueTList["Statistic"];
		permuted_data(_, 0) = data(_, 0);
		for (i = 0; i < n; i++)
			permuted_data(i, 1) = data(Permutations(i, 0)-1, 1); // save one example

		// Compute the statistics value for each permutation:
		NumericVector statistics_under_null(B);
		permuted_data(_, 0) = data(_, 0);
//		Rcout << "Start permutations Compute Stat" << endl; 
		for (ctr = 0; ctr < B; ctr++)
		{

			for (i = 0; i < n; i++)
				permuted_data(i, 1) = data(Permutations(i, ctr)-1, 1); // subtract 1 (from R back to C indexing!)
			statistics_under_null[ctr] = ComputeStatistic_rcpp(permuted_data, grid_points, expectations_table);
		}
		output["TrueT"] = TrueT;
		output["statistics.under.null"] = statistics_under_null;
		output["Permutations"] = Permutations;
		output["permuted.data"] = permuted_data;
	} // end permutations test
	/**/
	if(test_type == "permutations_inverse_weighting") // Use inverse-weighing Hoeffding's statistic, with MCMC to give pvalue
	{
		NumericMatrix w_mat = w_fun_to_mat_rcpp(data, w_fun);
		TrueT = ComputeStatistic_w_rcpp(data, grid_points, w_fun, counts_flag); //  $Statistic
				
		List PermutationsList = PermutationsMCMC_rcpp(w_mat, prms); // burn.in = prms$burn.in, Cycle = prms$Cycle)
		NumericMatrix Permutations = PermutationsList["Permutations"];
		NumericMatrix permuted_data(n, 2);

		// Compute the statistics value for each permutation:
		NumericVector statistics_under_null(B);
//		NumericMatrix w_mat_permutation(n, n);
		for (ctr = 0; ctr < B; ctr++)
		{
			permuted_data(_, 0) = data(_, 0);
			for (i = 0; i < n; i++)
				permuted_data(i, 1) = data(Permutations(i, ctr)-1, 1); // subtract 1 (from R back to C indexing!)
//			permuted_sample = cbind(data(_, 0), data(Permutations(_, ctr), 1));
//			w_mat_permutation = w_fun_to_mat_rcpp(permuted_data, w_fun);
			NullT = ComputeStatistic_w_rcpp(permuted_data, grid_points, w_fun, counts_flag);
			statistics_under_null[ctr] = NullT;
		}
		output["TrueT"] = TrueT;
		output["statistics.under.null"] = statistics_under_null;
		output["Permutations"] = Permutations;
	} // end permutations with inverse weighting test

	/**/
	if (test_type == "uniform_importance_sampling") { // Our modified Hoeffding's statistic with importance sampling uniform permutations test (not working yet)
		output = IS_permute_rcpp(data, grid_points, w_fun, prms, test_type);} // instead of w_fun we should give expectation table 

	if (test_type == "monotone_importance_sampling") { // Our modified Hoeffding's statistic with importance sampling uniform permutations test (not working yet)
		prms["importance.sampling.dist"] = "monotone.w";
		output = IS_permute_rcpp(data, grid_points, w_fun, prms, test_type);} // instead of w_fun we should give expectation table 
	if (test_type == "monotone_importance_sampling_inverse_weighting") { // Our modified Hoeffding's statistic with importance sampling uniform permutations test (not working yet)
		prms["importance.sampling.dist"] = "monotone.w";
		output = IS_permute_rcpp(data, grid_points, w_fun, prms, test_type);} // instead of w_fun we should give expectation table 

	if (test_type == "match_importance_sampling") { // Our modified Hoeffding's statistic with importance sampling uniform permutations test (not working yet)
		prms["importance.sampling.dist"] = "match.w";
		output = IS_permute_rcpp(data, grid_points, w_fun, prms, test_type);} // instead of w_fun we should give expectation table 


	if (test_type == "uniform_importance_sampling_inverse_weighting") { // inverse-weighting Hoeffding's statistic with importance sampling uniform permutations  
//		NumericMatrix expectations_table(1, 1); // set empty (0) table 
		output = IS_permute_rcpp(data, grid_points, w_fun, prms, test_type); // need to change!
	} // w = function(x) { 1 }) {
	if(test_type == "tsai") { cout << "Can't run Tsai's test from cpp" << endl; } //   Tsai's test, relevant only for truncation W(x,y)=1_{x<=y}		
	if(test_type == "minP2") { cout << "Can't run minP2 from cpp" << endl;} // minP2 test, relevant only for truncation W(x,y)=1_{x<=y}

	if (!(output.containsElementNamed("permuted.data")))
		output["permuted.data"] = NULL;
	if (!(output.containsElementNamed("Pvalue"))) //    ("Pvalue" % in % names(output))) // Compute empirical P - value
	{
//		output["Pvalue"] = length(which(output["statistics_under_null"] >= output["TrueT"])) / B;
		double Pvalue = 1.0; // include id permutation! 
		NumericVector statistics_under_null = as<NumericVector>(output["statistics.under.null"]);
		for (i = 0; i < B; i++)
			Pvalue += (statistics_under_null[i] >= as<double>(output["TrueT"]));
		output["Pvalue"] = Pvalue / double(B+1); // include id permutation
//		Rcout << "Output Pvalue = " << Pvalue / double(B+1) << endl; 
	}

	for (i = 0; i < n; i++) // long(grid_points_arma.n_rows); i++)
	  for (j = 0; j < 2; j++)
	  {
//		  if(i > n-5)
//			  Rcout << "Copy i,j " << i << ", " << j << " data, new_data " << data(i, j) << ", " << new_data(i, j) << endl; 
	      data(i, j) = new_data(i, j); // copy back !! 
	  }
//	Rcout << "Ended Copy back " << endl; 
		
	return(output);
}
