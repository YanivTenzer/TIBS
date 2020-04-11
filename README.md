# TIBS
TIBS (Testing Independence with Biased Sampling) is a repository of R functions which implement tests for independence under general  known baised sampling: 

More formally, for (x<sub>1</sub>,y<sub>1</sub>), ..., (x<sub>n</sub>, y<sub>n</sub>) ~ F<sub>XY</sub> W with a known weighting function W, TIBS tests the null hypothesis
H<sub>0</sub>: F<sub>XY</sub> = F<sub>X</sub> F<sub>Y</sub>

## Tests: 
The repository implements TIBS - this is a new test based on multiple data-dependent partitions of the (x,y) plane, 
with modified expectation computed based on W.  <br>
In addition, the repository also enables the user to run Tsai's test and the minP2 test - both tests are 
applicable for the special case of truncation: W(x,y) = 1<sub>{y>x}</sub>. 

## Getting started: 
Clone the repository into a directory of your choice. 
Change the path in the function 'run_simulations.R' into the directory you've used. 

The file 'TIBS.R' contains the main function 'TIBS' which produces a test for your dataset. 
For example, to run a TIBS permutation test with B=100 permutations for biased sample W(x,y)=x+y, use: <br>
test_results<-TIBS(biased_data, 'sum', '', 100, 'permutations', prms) <br>
where biased_data is a n*2 array containing your (x<sub>i</sub>,y<sub>i</sub>) sample from F<sub>XY</sub> W


The script 'run_simulations.R' simulates data under multiple dependency distributions Fsub>XY</sub> and biased sampling functions W.
It then applies the TIBS test and other tests to compute significance and type-1 and type-2-errors.
Warning: running time is quite high for this script. to get started, you may want to change the script to 
simulate and test only one dataset, with a small number of iterations and permutations/bootstrap samples (B). 


## Datasets: 
There are currently 4 real-life datasets analyzed by TIBS: 
The huji and AIDS datasets are part of the repository. 
Additional datasets (ICU and Infection) are available upon request. 


The script 'real_life_datasets.R'  runs the TIBS tests on the real datasets.

## Contacts: 
For any questions, bug-reports etc. please contact: <br>
Yaniv Tenzer   yaniv.tenzer@gmail.com <br>
Or Zuk  or.zuk@mail.huji.ac.il
