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


int main()
{
    long i, j, d;
    cout << "Hello World!\n";

    const char* datasets[4] = { "huji", "AIDS", "ICU", "Infection" }; // , "Dementia"


    string data_dir = "C:/Users/Or Zuk/Dropbox/BiasedSampling/data";
    string data_file;
//  ofstream data_str;
    string data_arr[1000];


    long n = 100; // Number of data points 
    double *data[2]; // data array 
    double* grid_points[2];
    double** null_expectation_table;

    long num_lines;

    // read data file: 
    for (d = 2; d < 3; d++)  // loop on datasets. Currently only ICU
    {
        switch (d) {
        case 0: {cout << datasets[d] << '\n'; break; } // "huji" prints "1",
        case 1: {cout << "2\n";  break; }// "AIDS"  then prints "2"
        case 2: {
            data_file.assign(data_dir + "/ICU_data.txt");
            ifstream data_str(data_file);

            num_lines = count_lines_in_file(data_file);
            n = num_lines - 1;
            cout << "Total lines:" << num_lines;
            for (j = 0; j < 2; j++)
            {
                data[j] = new double[num_lines];
            }
            
            //            data.open(data_file, ios::in); // | ios::app | ios::binary);
            i = 0;
            //           while (i < 10) {
            //               data >> data_arr[i++];
             //          }

//            while (!data.eof())
            cout << "Input File Name: " << data_file << "\n";
            if (!data_str.is_open()) {
                cout << " Failed to open\n";
            }
            else {
                for (int k = 0; k < num_lines; k++) // loop on lines 
                {
                    getline(data_str, data_arr[i++], '\n');
                    cout << " READ: " << data_arr[i - 1] << '\n';

                    // position = data_arr[i - 1].find(' ', position);
                    // first_section = str1.substr(0, position);

                    stringstream stream(data_arr[i - 1]);

                    if(k>0) // skip first line
                        for (int m = 0; m < 2; m++)
                        {
                            stream >> data[m][k-1];
                            if (m == 1)
                                data[m][k - 1] -= 0.02; // correction 
                            cout << " In reals: ";
                            cout << data[m][k-1];
                        }
                    cout << endl;
                }
//                cout << "FILE IS OPEN:\n";
//                for (int k = 0; k < num_lines-1; k++)
//                    cout << data[0][k] << endl;


            }
            break;
        }
              //            data.close()
        case 3: {cout << "4\n";  break; } // "Infection" then prints "2"

        }
    }

    // Run Hoeffding's test on it: 

    for (j = 0; j < 2; j++)
    {
        grid_points[j] = new double[n];
        for (i = 0; i < n; i++)
            grid_points[j][i] = data[j][i];
    }
    cout << "Hoefdings Test Statistics on Data: " << ComputeStatistic(n, data, grid_points, null_expectation_table) << "\n";


}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
