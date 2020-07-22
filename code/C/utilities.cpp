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
#include "Rcpp.h" // for including R with Rcpp


using namespace std;
using namespace Rcpp;  // for including R with Rcpp 

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
			Exp[j] = null_expectations_table[j][i];
//			if (i < 10)
//				cout << "Ind: " << i << " " << j << " Exp: " << Exp[j] << endl;
		}

		for (j = 0; j < n; j++)  // loop on data points  
		{
			Rx[j] = data[0][j] > grid_points[0][i];
			Ry[j] = data[1][j] > grid_points[1][i];
			Obs[0] += double(Rx[j] * Ry[j]);
			Obs[1] += Rx[j];
			Obs[3] += Ry[j];
		}
		Obs[1] -= Obs[0];
		Obs[3] -= Obs[0];
		Obs[2] = n - Obs[0] - Obs[1] - Obs[3];

		if ((Exp[0] > 1) && (Exp[1] > 1) && (Exp[2] > 1) && (Exp[3] > 1))
		{
//			cout << "Update Statistic" << pow((Obs[j] - Exp[j]), 2) / Exp[j];
			for (j = 0; j < 4; j++)
				Statistic += pow((Obs[j] - Exp[j]), 2) / Exp[j];  // set valid statistic when expected is 0 or very small
		}
	} // end loop on grid points

	return(Statistic);
}


// Hoeffding's test in Rcpp

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
double ComputeStatistic(long n, double** data, double** grid_points, double *null_expectations_table[4])
{
	double Obs[4] = { 0 };
	double Exp[4] = { 0 };
	long i, j;
	double Statistic = 0.0;
	long *Rx = new long[n];
	long *Ry = new long[n];

	for (i = 0; i < n; i++) // Slow loop on grid points
	{
		for (j = 0; j < 4; j++) // loop on quadrants
		{
			Exp[j] = null_expectations_table[j][i];
			if(i < 10)
				cout << "Ind: " << i << " " << j << " Exp: " << Exp[j] << endl;
		}

		for (j = 0; j < n; j++)  // loop on data points  
		{
			Rx[j] = data[0][j] > grid_points[0][i];
			Ry[j] = data[1][j] > grid_points[1][i];
			Obs[0] += double(Rx[j] * Ry[j]);
			Obs[1] += Rx[j];
			Obs[3] += Ry[j];
		}		
		Obs[1] -= Obs[0];
		Obs[3] -= Obs[0];
		Obs[2] = n - Obs[0] - Obs[1] - Obs[3];

		if ((Exp[0] > 1) && (Exp[1] > 1) && (Exp[2] > 1) && (Exp[3] > 1))
		{
			cout << "Update Statistic" << pow((Obs[j] - Exp[j]), 2) / Exp[j];
			for (j = 0; j < 4; j++)
				Statistic += pow((Obs[j] - Exp[j]), 2) / Exp[j];  // set valid statistic when expected is 0 or very small
		}
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
//		null_distribution_CDF < -PDFToCDF2d(null.distribution, data)
//
//		for (i in seq(1, dim(grid.points)[1], 1))
//		{
//			for (j in 1 : 3)
//				mass.table[i, j] < -GetQuarterExpectedProb(grid.points[i, ], j, data, null_distribution_CDF)
//				mass.table[i, 4] = 1 - sum(mass.table[i, 1:3]) # , epsilon)
//		}
//	mass.table < -dim(data)[1] * mass.table # normalize to counts
//}


double QuarterProbFromBootstrap( double *data[2], double **null_distribution, double *grid_points[2], long n, 
	double *mass_table[4])
{
	long i, j;
	
	double** null_distribution_CDF = new double*[n];

	double cur_grid_points[2];
	for (i = 0; i < n; i++)  // loop on grid-points
		null_distribution_CDF[i] = new double[n];
//	cout << "Inside DO PDFToCDF2d\n";
	PDFToCDF2d(null_distribution, data, n, null_distribution_CDF);  // convert PDF to CDF 

//	cout << "Inside DO GetQuarterExpectedProb\n";
	for (i = 0; i < n; i++)  // loop on grid-points
	{
		cur_grid_points[0] = grid_points[0][i];
		cur_grid_points[1] = grid_points[1][i];
		cout << "cur-grid-points i=" << i << " " << cur_grid_points[0] << ", " << cur_grid_points[1] << endl;
		for (j = 0; j < 3; j++)
			mass_table[j][i] = GetQuarterExpectedProb(cur_grid_points, j, data, null_distribution_CDF, n);
		mass_table[3][i] = 1 - mass_table[2][i] - mass_table[1][i] - mass_table[0][i]; 
//		cout << "Inside Finished i=" << i << endl;
	}

	for(j = 0; j < 4; j++)
		for (i = 0; i < n; i++)  //  normalize to counts
		{
			mass_table[j][i] *= n;
			if(i < 10)
				cout << "Index: " << j << " " << i << " Mass-Table: " << mass_table[j][i] << endl;
		}

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

////////////////////////////////////////////////////////////////////////////////////
// Compute 2d cumulative distribution.When we have ties we need to correct this
// Input:
// pdf_2d - a two-dim array n * n of probabilities
// data - xy points with probabilities, array size n*2 (used for sorting)
// Output:
// cdf_2d - a two-dim array n * n of cumulative probabilities
////////////////////////////////////////////////////////////////////////////////////

double PDFToCDF2d(double **pdf_2d, double *data[2], long n, double **cdf_2d)
{
	long i, j;
	// sort 
	long *Px = new long[n];
	long *Py = new long[n];


	sort_with_indexes(data[0], n, Px);
	sort_with_indexes(data[1], n, Py); //	sort(data[0], data[0] + n); // why +n? 

	double** cdf_2d_temp = new double*[n];
	for (i = 0; i < n; i++)
		cdf_2d_temp[i] = new double[n];

	// cumulative sum on rows 
	for(i = 0; i < n; i++)
	{
		cdf_2d_temp[Px[i]][Py[0]] = pdf_2d[Px[i]][Py[0]];
		for (j = 1; j < n; j++)
			cdf_2d_temp[Px[i]][Py[j]] = cdf_2d_temp[Px[i]][Py[j-1]] + pdf_2d[Px[i]][Py[j]];
	}
	// cumulative sum on columns 
	for (j = 0; j < n; j++)
	{
		cdf_2d_temp[Px[0]][Py[j]] = pdf_2d[Px[0]][Py[j]];
		for (i = 1; i < n; i++)
			cdf_2d_temp[Px[i]][Py[j]] = cdf_2d_temp[Px[i-1]][Py[j]] + pdf_2d[Px[i]][Py[j]];
	}
	
	// apply inverse permutation to rows and columns
//	long* Px_inv = new long[n];
//	long* Py_inv = new long[n];
//	inv_perm(Px, n, Px_inv);
//	inv_perm(Py, n, Py_inv);

	// copy results permuted
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			cdf_2d[Px[i]][Py[j]] = cdf_2d[i][j];
	return(0);

}


/*
###################################################################################################
# given a quadrant, evaluate the mass function within it
# New version - uses the two - dimensional CDF(not PDF)
# 
# Point - (x, y) value defining quadrants
# QId - which quadrant 1, 2, 3, 4
# data - 2 * n array of(x, y)
# null_distribution_CDF - n * n table of 2d cumulative of Fx * Fy * w(problem with ties)
###################################################################################################
GetQuarterExpectedProb < -function(Point, QId, data, null_distribution_CDF)
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
			{ S < -null_distribution_CDF[m, n] + null_distribution_CDF[idx.x, idx.y] -
			  null_distribution_CDF[idx.x, n] - null_distribution_CDF[m, idx.y] }, # 1
		{S < -null_distribution_CDF[m, idx.y] - null_distribution_CDF[idx.x, idx.y]}, # 1
		{S < -null_distribution_CDF[idx.x, idx.y]}, # 3
		{S < -null_distribution_CDF[idx.x, n] - null_distribution_CDF[idx.x, idx.y]}) # 4
			return(S)
}
*/

double GetQuarterExpectedProb(double Point[2], long QId, double *data[2], double **null_distribution_CDF, long n)
{
	// convert below R to C++:
	long idx_x = -1, idx_y = -1, i;
	double S; 

	if ((QId == 0) || (QId == 1)) // change to 0-3 coordinates for quadrants !!
	{
		for (i = 0; i < n; i++)
			if (data[0][i] > Point[0])
				if ((idx_x == -1) || (data[0][idx_x] > data[0][i]))
					idx_x = i;
	}
	else
	{
		for (i = 0; i < n; i++)
			if (data[0][i] <= Point[0])
				if ((idx_x == -1) || (data[0][idx_x] < data[0][i]))
					idx_x = i;
	}
	if ((QId == 0) || (QId ==  3))
	{
		for (i = 0; i < n; i++)
			if (data[1][i] > Point[1])
				if ((idx_y == -1) || (data[1][idx_y] > data[1][i]))
					idx_y = i;
	}
	else
	{
		for (i = 0; i < n; i++)
			if (data[1][i] <= Point[1])
				if ((idx_y == -1) || (data[1][idx_y] < data[1][i]))
					idx_y = i;
	}

	if ((idx_x == -1) || (idx_y == -1))
		return(0);

	cout << "Data for max-index:";
	for (i = 0; i < 10; i++)
		cout << data[0][i] << " , " << data[1][i] << endl;

	long idx_x_max = max_index(data[0], n);
	long idx_y_max = max_index(data[1], n);


	switch (QId) { // different quardants
	case 0: {S = null_distribution_CDF[idx_x_max][idx_y_max] + null_distribution_CDF[idx_x][idx_y] -
		null_distribution_CDF[idx_x][idx_y_max] - null_distribution_CDF[idx_x_max][idx_y]; break; } // "huji" prints "1",
	case 1: {S = -null_distribution_CDF[idx_x_max][idx_y] - null_distribution_CDF[idx_x][idx_y]; break; }
	case 2: {S = -null_distribution_CDF[idx_x][idx_y]; break; }
	case 3: {S = -null_distribution_CDF[idx_x][idx_y_max] - null_distribution_CDF[idx_x][idx_y]; break; }
	} // end switch 

	cout << "Mass Table Indexes: idx_x=" << idx_x << ", idx_y=" << idx_y << ", idx_x_max=" << idx_x_max << ", idx_y_max=" << idx_y_max << " S=" << S << endl;

	return(S);
}



/*
###################################################################################
# Estimate the null distribution fx* fy* W(given the estimated PDFs f_x, f_y)
# 
# Parameters:
# pdfs - 2 * n table with marginal distributions fx, fy probabilities
# w - n * n matrix with weights W[i, j] = w(x_i, y_j)
###################################################################################
GetNullDistribution < -function(pdfs, W)
{
	# Compute the normalizing factor under the null :
	null.distribution < -W * (pdfs[, 1] % *%t(pdfs[, 2]))
		Z < -sum(null.distribution)
		null.distribution < -null.distribution / Z
		return(list(null.distribution = null.distribution, normalizing.factors = Z))
}
*/


double GetNullDistribution(double *pdfs[2], double **w_mat, long n, double **null_distribution)
{
	// Compute the normalizing factor under the null :
	long i, j;
	double z = 0;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			null_distribution[i][j] = w_mat[i][j] * pdfs[0][i] * pdfs[1][j];
			if (i < 10)
				cout << i << ", " << j << " W: " << w_mat[i][j] << " pdf0= " << pdfs[0][i] << " pdf1=" << pdfs[1][i] << endl; // printing: 

			z += null_distribution[i][j];
		}
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			null_distribution[i][j] /= z;
//			if (i < 10)
	//			cout << i << ", " << j << " null_distribution_pdf: " << null_distribution[i][j] << endl; // printing: 
		}
	return(z);
}


// General utilities below

// Compute inverse permutation
long inv_perm(long* p, long n, double* inv_p)
{
	long i;
	for (i = 0; i < n; i++)
		inv_p[p[i]] = i; 
	return(TRUE);
}


// compute empirical cdf 
long empirical_cdf(double* x, long n, double* ecdf)
{
	long i; 
	long* Px = new long[n];
	
	sort_with_indexes(x, n, Px); // sort index
	for (i = 0; i < n; i++)
		ecdf[Px[i]] = i / n; 
	return(TRUE);
}


// cumsum 
long cumsum(double* x, long n, double* x_cumsum)
{
	x_cumsum[0] = x[0]; 
	for (long i = 1; i < n; i++)
		x_cumsum[i] = x_cumsum[i - 1] + x[i];
	return(TRUE);
}


// find maximum element of an array
long max_index(double* x, long n)
{
	long i;
	long max_ind = 0;
	double max_val = x[0];

	for (i = 1; i < n; i++)
	{
		if (x[i] > max_val)
		{
			max_val = x[i];
			max_ind = i;
		}
	}

	return(max_ind);
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

//constexpr unsigned int str2int(const char* str, long h = 0)
//{
//	return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
//}

long sort_with_indexes(double* data, long n, long* sort_perm)
{
	unsigned long i;
	vector<int> index(n, 0);
	for (i = 0; i != index.size(); i++)
		index[i] = i;
	sort(index.begin(), index.end(),
		[&](const int& a, const int& b) {
			return (data[a] < data[b]);
		}
	);
	for (i = 0; i != index.size(); i++)
		sort_perm[i] = index[i];
	return(TRUE);

}
