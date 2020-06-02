# Utilities from cpp

cppFunction('double ComputeStatistic_cpp(long n, double** data, double** grid_points, double *null_expectations_table[4])
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
}')
