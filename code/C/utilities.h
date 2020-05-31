#ifndef _UTILITIES_
#define _UTILITIES_

double ComputeStatistic(long n, double** data, double** grid_points, double** null_expectations_table);
long count_lines_in_file(std::string file_name);
double QuarterProbFromBootstrap(double* data[2], double null_distribution, double* grid_points[2], long n,
	double* mass_table[4]);
double PDFToCDF2d(double* pdf_2d, double* data[2]);
double TIBS(double* data[2], string w_fun, string test_type, string prms, long n);
double GetQuarterExpectedProb(double Point[2], long QId, double* data[2], double** null_distribution_CDF, long n);
vector<size_t> sort_indexes(const vector<T>& v);
#endif
