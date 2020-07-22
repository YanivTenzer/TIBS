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
// #include "utilities.h"
#include <Rcpp.h> // for including R with Rcpp


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

	// First print input to see that it's copied well: 
	/**
	for (i = 0; i < 5; i++) // Slow loop on grid points
		cout << "Data: " << data(i,0) << " " << data(i,1) << endl;

	for (i = 0; i < 5; i++) // Slow loop on grid points
		cout << "grid.points " << grid_points(i, 0) << " " << grid_points(i, 1) << endl;

	for (i = 0; i < 5; i++) // Slow loop on grid points
	{
		for (j = 0; j < 4; j++) // loop on quadrants
			cout << "grid.points " << null_expectations_table(i, j) << " " ;
		cout << endl;
	}
	**/

//	cout << "Starting Loops" << endl;
	for (i = 0; i < n; i++) // Slow loop on grid points
	{
		for (j = 0; j < 4; j++) // loop on quadrants
		{
			Exp[j] = null_expectations_table(i, j);
			Obs[j] = 0;
//			cout << "Expected: " << Exp[j] << endl;
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
//		cout << "Obs: " << Obs[0] << " " << Obs[1] << " " << Obs[2] << " " << Obs[3] << " " << endl;

		if ((Exp[0] > 1) && (Exp[1] > 1) && (Exp[2] > 1) && (Exp[3] > 1))
		{
			for (j = 0; j < 4; j++)
				Statistic += pow((Obs[j] - Exp[j]), 2) / Exp[j];  // set valid statistic when expected is 0 or very small
		}
	} // end loop on grid points

	return(Statistic);
}


// Hoeffding's test in Rcpp

