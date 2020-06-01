// TIBS.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// #include "stdafx.h"
#include "windows.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "utilities.h"

using namespace std;

// Main function 
double TIBS(double* data[2], string w_fun, string test_type, double *params, long n)
{
	double* mass_table[4];
	long i, j;
	for (j = 0; j < 4; j++)
	{
		mass_table[j] = new double[n];
		for (i = 0; i < n; i++)
			mass_table[j][i] = 0.0;
	}

	// compute grid_points 
	double* grid_points[2];

	for (j = 0; j < n; j++)
	{
		grid_points[j] = new double[n];
		for (i = 0; i < n; i++)
			grid_points[j][i] = data[j][i];
	}
	//	grid.points < -unique.matrix(grid.points)  # set unique for ties ? for discrete data

	double** null_distribution = new double* [n];
	for (j = 0; j < n; j++)
		null_distribution[j] = new double[n];

	double** null_expectation_table = new double* [n];
	for (j = 0; j < n; j++)
		null_expectation_table[j] = new double[n];



	/** R Code: 
		marginals < -EstimateMarginals(data, w.fun)
		w.mat = w_fun_to_mat(marginals$xy, w.fun)
		null.distribution < -GetNullDistribution(marginals$PDF, W)
		expectations.table < -QuarterProbFromBootstrap(marginals$xy, null.distribution$null.distribution, grid.points)
	**/

	//	double QuarterProbFromBootstrap(double* data[2], double** null_distribution, double* grid_points[2], long n,
	//		double* mass_table[4])

	double** w_mat = new double* [n];
	for (i = 0; i < n; i++)
		w_mat[i] = new double[n];



	double* marginals_xy[2], *marginals_pdfs[2], *marginals_cdfs[2];
	for (i = 0; i < 2; i++)
	{
		marginals_xy[i] = new double[2*n];
		marginals_pdfs[i] = new double[2*n];
		marginals_cdfs[i] = new double[2*n];
	}

//	stimateMarginals(double* data[2], string w_fun, double* params, long n, // inputs 
//		double* xy[2], double* PDFs[2], double* CDFs[2])
	EstimateMarginals(data, w_fun, params, n, marginals_xy, marginals_pdfs, marginals_cdfs); // estimate marginals 
	w_fun_to_mat(marginals_xy, w_fun, n, w_mat); // 
	double z = GetNullDistribution(marginals_pdfs, w_mat, n, null_distribution);


	// get null expectation
	QuarterProbFromBootstrap(data, null_distribution, grid_points, n, null_expectation_table);
	// R: expectations.table < -QuarterProbFromBootstrap(marginals$xy, null.distribution$null.distribution, grid.points)

	return(ComputeStatistic(n, data, grid_points, null_expectation_table)); // call function and return 
}
