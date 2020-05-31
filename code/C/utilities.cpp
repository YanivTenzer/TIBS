#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

#include <time.h>
#include <math.h>
#include "utilities.h"


using namespace std;

// sort with indexes. From here: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
vector<size_t> sort_indexes(const vector<T>& v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}


// Hoeffding's test in c


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
double ComputeStatistic(long n, double** data, double** grid_points, double **null_expectations_table)
{
	double Obs[4];
	double Exp[4];
	long i, j;
	double  Statistic = 0.0;
	long *Rx = new long[n];
	long* Ry = new long[n];

//	double* obs_table[4];
//	for (j = 0; j < 4; j++)  // loop on quadrants
//		obs_table[j] = new double[n];

	for (i = 0; i < n; i++) // Slow loop on grid points
	{
		for (j = 0; j < 4; j++)  // loop on quadrants 
		{
			Exp[j] = null_expectations_table[i][j];
			Rx[j] = data[j][0] > grid_points[i][0];
			Ry[j] = data[j][1] > grid_points[i][1];
			Obs[0] += Rx[j] * Ry[j];
			Obs[1] += Rx[j];
			Obs[3] += Ry[j];
		}
		Obs[1] -= Obs[0];
		Obs[3] -= Obs[0];
		Obs[2] = n - Obs[1] - Obs[2] - Obs[4];

//		for (j = 0; j < 4; j++)
//			obs_table[j][i] = -Obs[j];
		if ((Exp[0] > 1) && (Exp[1] > 1) && (Exp[2] > 1) && (Exp[3] > 1))
			for (j = 0; j < 4; j++)
				Statistic += pow((Obs[j] - Exp[j]), 2) / Exp[j];  // set valid statistic when expected is 0 or very small
	} // end loop on grid points

	return(Statistic);
}


// Compute expectation table
//###################################################################################################
//# Compute Expect[Qi(p_j)] for 1 <= i <= 4, and all j, given a grid of pointsand bootstrap null distribution
//# Parameters:
//# data - 2 * n array(X, Y)
//# Permutations - set of permutations
//# grid.points - centers of partitions
//#
//# Output:
//# mass.table - a 4 * #grid - points table with quadrants expectation
//###################################################################################################
//QuarterProbFromBootstrap < -function(data, null.distribution, grid.points)
//{
//	mass.table < -matrix(0, dim(grid.points)[1], 4)
//		null.distribution.CDF < -PDFToCDF2d(null.distribution, data)
//
//		for (i in seq(1, dim(grid.points)[1], 1))
//		{
//			for (j in 1 : 3)
//				mass.table[i, j] < -GetQuarterExpectedProb(grid.points[i, ], j, data, null.distribution.CDF)
//				mass.table[i, 4] = 1 - sum(mass.table[i, 1:3]) # , epsilon)
//		}
//	mass.table < -dim(data)[1] * mass.table # normalize to counts
//}


double QuarterProbFromBootstrap( double *data[2], double null_distribution, double *grid_points[2], long n, 
	double *mass_table[4])
{
	long i, j;
	
	double* null_distribution_CDF; 
	null_distribution_CDF = PDFToCDF2d(null_distribution, data);  // convert PDF to CDF 

	for (i = 0; i < n; i++)  // loop on grid-points
		for (j = 0; j < 3; j++)
			mass_table[j][i] = GetQuarterExpectedProb(grid_points[i, ], j, data, null_distribution_CDF);
		mass_table[3][i] = 1 - mass_table[2][i] - mass_table[1][i] - mass_table[0][i]; 
	}
	
	for(j = 0; j < 4; j++)
		for(i = 0; i < n; i++)  //  normalize to counts
			mass_table[j][i] *= n;

	return(TRUE); 
}  // end function QuarterProbFromBootstrap



/**
###################################################################################################
# Compute 2d cumulative distribution.When we have ties we need to correct this
#
# pdf.2d - a two - dim array of probabilities
# data - xy points with probabilities(used for sorting)
###################################################################################################
PDFToCDF2d < -function(pdf.2d, data)
{
	Px < -sort(data[, 1], index.return = TRUE)  # Permute to order x_i, y_i
		Py < -sort(data[, 2], index.return = TRUE)
		cdf.2d < -apply(apply(pdf.2d[Px$ix, Py$ix], 2, cumsum), 1, cumsum)  # cumsum on rows and columns

		# Use data to deal with ties(not working yet)
		#  ties.x < -which(data[-1, 1] == head(data[, 1], -1))
		#  ties.y < -which(data[-1, 2] == head(data[, 2], -1))
		#  for (i in rev(ties.y))
		#    cdf.2d[, i] < -cdf.2d[, i + 1]
		#  for (i in rev(ties.x))
		#    cdf.2d[i, ] < -cdf.2d[i + 1, ]

		return(t(cdf.2d[invPerm(Py$ix), invPerm(Px$ix)]))
}
**/


double PDFToCDF2d(double *pdf_2d, double *data[2])
{
	// sort 
	
	sort_indexes(data[0])

	sort(data[0], data[0] + n); // why +n? 

	Px < -sort(data[, 1], index.return = TRUE)  # Permute to order x_i, y_i
		Py < -sort(data[, 2], index.return = TRUE)
		cdf.2d < -apply(apply(pdf.2d[Px$ix, Py$ix], 2, cumsum), 1, cumsum)  # cumsum on rows and columns

		return(t(cdf.2d[invPerm(Py$ix), invPerm(Px$ix)]))
}


/*
###################################################################################################
# given a quadrant, evaluate the mass function within it
# New version - uses the two - dimensional CDF(not PDF)
# 
# Point - (x, y) value defining quadrants
# QId - which quadrant 1, 2, 3, 4
# data - 2 * n array of(x, y)
# null.distribution.CDF - n * n table of 2d cumulative of Fx * Fy * w(problem with ties)
###################################################################################################
GetQuarterExpectedProb < -function(Point, QId, data, null.distribution.CDF)
{
	if (QId % in % c(1, 2))
	{
		idx.x < -which(data[, 1] > Point[1])
			idx.x < -idx.x[which.min(data[idx.x, 1])]
	}
	else
	{
		idx.x < -which(data[, 1] <= Point[1])
			idx.x < -idx.x[which.max(data[idx.x, 1])]
	}
	if (QId % in % c(1, 4))
	{
		idx.y < -which(data[, 2] > Point[2])
			idx.y < -idx.y[which.min(data[idx.y, 2])]
	}
	else
	{
		idx.y < -which(data[, 2] <= Point[2])
			idx.y < -idx.y[which.max(data[idx.y, 2])]
	}

	if (isempty(idx.x) | isempty(idx.y))
		return(0)
		m < -which.max(data[, 1])
		n < -which.max(data[, 2])

		switch (QId, # First sample from Fxy
			{ S < -null.distribution.CDF[m, n] + null.distribution.CDF[idx.x, idx.y] -
			  null.distribution.CDF[idx.x, n] - null.distribution.CDF[m, idx.y] }, # 1
		{S < -null.distribution.CDF[m, idx.y] - null.distribution.CDF[idx.x, idx.y]}, # 1
		{S < -null.distribution.CDF[idx.x, idx.y]}, # 3
		{S < -null.distribution.CDF[idx.x, n] - null.distribution.CDF[idx.x, idx.y]}) # 4
			return(S)
}
*/

double GetQuarterExpectedProb(double Point[2], long QId, double *data[2], double **null_distribution_CDF, long n)
{
	// convert below R to C++:

	if (QId % in % c(1, 2))
	{
		idx.x < -which(data[, 1] > Point[1])
			idx.x < -idx.x[which.min(data[idx.x, 1])]
	}
	else
	{
		idx.x < -which(data[, 1] <= Point[1])
			idx.x < -idx.x[which.max(data[idx.x, 1])]
	}
	if (QId % in % c(1, 4))
	{
		idx.y < -which(data[, 2] > Point[2])
			idx.y < -idx.y[which.min(data[idx.y, 2])]
	}
	else
	{
		idx.y < -which(data[, 2] <= Point[2])
			idx.y < -idx.y[which.max(data[idx.y, 2])]
	}

	if (isempty(idx.x) | isempty(idx.y))
		return(0)
		m < -which.max(data[, 1])
		n < -which.max(data[, 2])

		switch (QId, # First sample from Fxy
			{ S < -null.distribution.CDF[m, n] + null.distribution.CDF[idx.x, idx.y] -
			  null.distribution.CDF[idx.x, n] - null.distribution.CDF[m, idx.y] }, # 1
		{S < -null.distribution.CDF[m, idx.y] - null.distribution.CDF[idx.x, idx.y]}, # 1
		{S < -null.distribution.CDF[idx.x, idx.y]}, # 3
		{S < -null.distribution.CDF[idx.x, n] - null.distribution.CDF[idx.x, idx.y]}) # 4
			return(S)
}



double TIBS(double *data[2], string w_fun, string test_type, string prms, long n)
{
	double* mass_table[4];
	for (j = 0; j < 4; j++)
	{
		mass_table[j] = new double[n];
		for (i = 0; i < n; i++)
			mass_table[j][i] = 0.0;
	}


	return(ComputeStatistic(n, data, grid_points, null_expectation_table)) // call function and return 
}


// Count total lines in file 
long count_lines_in_file(string file_name)
{
	long count = 0;
	string line;
	ifstream file_s(file_name);
	if (file_s.is_open())
	{
		while (!file_s.eof())
		{
			getline(file_s, line);
			cout << line << endl;
			count++;
		}
		file_s.close();
	}
	return(count);
}

