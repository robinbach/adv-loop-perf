#include <vector>
#include <iostream>
#include <math.h>
#include <Eigen/IterativeLinearSolvers>
#include <iomanip>
#include "dataImport.h"

// parallel computing library: OpenOMP library
// #include <omp.h>

// generate increaseing sequence
#include <numeric>

#include <deque>


using namespace std;
using namespace Eigen;

//Sigmoid function
double sigmoid(double x){
	return 1 / (1 + exp(x));
}

//Posterior probability
double eta(vector<double>& A, vector<double>& B){
	double temp = 0;
	for (int i = 0; i < A.size(); i++){
		temp =temp - A[i] * B[i];
	}
	return sigmoid(temp);
}

//Objective function
double obj(vector<vector<double> >& X, vector<double>& theta, vector<double>& y, double lambda){
	double sum = 0;
	for (int i = 0; i < y.size(); i++){
		double temp = 0;
		for (int j = 0; j < theta.size(); j++){
			temp += theta[j] * X[i][j];
		}
		sum += (1 - y[i])*temp + log(1 + exp(-temp));
	}
	double reg = 0;
	for (int k = 0; k < theta.size(); k++){
		reg += theta[k] * theta[k];
	}
	sum = sum + lambda*reg;
	return sum;
}

//Gradient
void grad(vector<vector<double> >& X, vector<double>& theta, vector<double>& sum, 
	vector<double>& y, double lambda){

	for (int i = 0; i < y.size(); i++){
		double temp = 0;
		for (int j = 0; j < theta.size(); j++){
			temp += theta[j] * X[i][j];
		}
		for (int k = 0; k < theta.size(); k++){
			sum[k] += X[i][k] * (1 - y[i] - sigmoid(temp));
		}
	}

	for (int i = 0; i < sum.size(); i++){
		sum[i] += 2 * lambda*theta[i];
	}
}

//Hessian
void hess(vector<vector<double> >& X, vector<vector<double> >& sum, vector<double>& theta,
	vector<double>& y, double lambda, float startTime){
	sum.clear();
	vector<double> oneRow(theta.size(), 0);

	for (int len = 0; len < theta.size(); len++){
		sum.push_back(oneRow);
	}

	float currentTime = omp_get_wtime( ) - startTime; // (clock() - t) / CLOCKS_PER_SEC
    	cout << "Hess 1.1  finished with " << currentTime << " seconds elapsed." << endl;


	int nthreads = omp_get_num_threads();

	deque<int> y_ids(y.size());

	iota(y_ids.begin(), y_ids.end(), 0);

	int th_id;
	int one_y;

	int theta_size = theta.size();
	// #pragma omp parallel private(th_id, one_y) shared(X, sum, theta, theta_size)
    {
    	// init mapper's local sum pool
    	vector<vector<double> > local_sum;
		vector<double> emptyRow(theta_size, 0);

		for (int len = 0; len < theta_size; len++)
		{
			local_sum.push_back(emptyRow);
		}

    	// th_id = omp_get_thread_num();

    	// #pragma omp master
    	// {
		while(!y_ids.empty())
    	{
    		vector<double> X_row_i;

    		one_y = -1;
            // any access to shared memory should be critical
            // #pragma omp critical
            {
                if(!y_ids.empty())
                {
                    one_y = y_ids.front();
                    y_ids.pop_front();
                    X_row_i = X[one_y];
                }
            }

    		double temp = 0;
    		if(one_y != -1)
    		{
    			for (int j = 0; j < theta_size; j++)
				{	 
					temp += theta[j] * X_row_i[j];
				}

				for (int k = 0; k < theta_size; k++)
				{
					for (int l = 0; l < theta_size; l++)
					{
						local_sum[l][k] += X_row_i[l] * X_row_i[k] * sigmoid(temp)*sigmoid(temp)*exp(temp);
					}
				}
    		}
			
    	}
    	
		for (int k = 0; k < theta_size; k++)
		{
			for (int l = 0; l < theta_size; l++)
			{
				// #pragma omp critical
				{
					sum[l][k] += local_sum[l][k];
				}
			}
		}
    			
		// }
    }

    // #pragma omp barrier
	

	currentTime = omp_get_wtime( ) - startTime; // (clock() - t) / CLOCKS_PER_SEC
    	cout << "Hess 1.2 finished with " << currentTime << " seconds elapsed." << endl;

	for (int k = 0; k < theta.size(); k++){
		for (int l = 0; l < theta.size(); l++){
			if (k==l){
				sum[k][l] += 2 * lambda;
			}
		}
	}
}

int main(int argc, char* argv[]){

	// prepare for parallel computing for Map reduce

	int th_id, nthreads;
	// th_id = omp_get_thread_num();
	// nthreads = omp_get_num_threads();

	// start timing
    float startTime = omp_get_wtime( ); //clock();
    float currentTime;

    // shared variables;
    int n, d;

    vector<vector<double> > X;
    vector<double> y;

    //Initialization
	double lambda = 5;

	// #pragma omp parallel
    {
		nthreads = omp_get_num_threads();
		// #pragma omp master
		{
			cout << "start running with " << nthreads << " cores." << endl;
		}
	}
	// Read in data
	string dataSetFileName = argv[1]; //"covtype_blank.data";
	RawDataSet rawDataSet;

	// only read in the first 30000 instance
	
	rawDataSet.readDataFromFile(dataSetFileName, 30000);
	cout << "Read file done" << endl;
	rawDataSet.printDataMatrixSize();

	
	
	for (int len = 0; len < rawDataSet.rawDataTable.size(); len++){
		y.push_back(rawDataSet.rawDataTable[len].back());



		/////////////////////
		//// NOTICE: code modified to fit the specific test case input: covertype
		////

		// choose covertype 2 as the flag variable
		if (y[len] == 2){
			y[len] = 1;
		}
		else
		{
			y[len] = 0;
		}
		///////////////////////


		
		rawDataSet.rawDataTable[len].pop_back();
		X.push_back(rawDataSet.rawDataTable[len]);
	}

	n = y.size();
	d = X[1].size();
	

	for (int i = 0; i < n; i++){
		X[i].push_back(1);
	}
	

	// #pragma omp parallel private(th_id) shared(nthreads, n, d, X, y, lambda)
	// {
	// 	#pragma omp master
	// 	{
	//Tolerance
	double epsilon = 1.0e-5;
	int iteration = 1000;

	//Gradient descent
	vector<double> theta(d + 1, 0);
	vector<double> theta2(d + 1, 0);

	for (int itr = 0; itr < iteration; itr++)
	{
		vector<double> G(d + 1, 0);
		grad(X, theta, G, y, lambda);

		currentTime = omp_get_wtime( ) - startTime; // (clock() - t) / CLOCKS_PER_SEC
    	cout << "Grad of " <<  itr << " finished with " << currentTime << " seconds elapsed." << endl;

		vector<vector<double> > H;
		hess(X, H, theta, y, lambda, startTime);

		currentTime = omp_get_wtime( ) - startTime; // (clock() - t) / CLOCKS_PER_SEC
    	cout << "Hess of " <<  itr << " finished with " << currentTime << " seconds elapsed." << endl;

		VectorXd x(d + 1), b(d + 1);
		MatrixXd A(d + 1, d + 1);

		for (int i = 0; i < d+1; ++i) {
			A.row(i) = VectorXd::Map(&H[i][0], H[i].size());
		}

		b = VectorXd::Map(&G[0], G.size());

		ConjugateGradient<MatrixXd> cg;
		cg.compute(A);
		x = cg.solve(b);

		//vector<vector<double> > H_inv;
		//MatrixInversion(H, d + 1, H_inv);

		for (int j = 0; j < d + 1; j++){
			double temp = x(j);
			theta2[j] = theta[j];
			theta[j] = theta[j] - temp;
		}

		double obj1 = obj(X, theta, y, lambda);
		double obj2 = obj(X, theta2, y, lambda);
		double delta = abs(obj1 - obj2);
		currentTime = omp_get_wtime( ) - startTime; // (clock() - t) / CLOCKS_PER_SEC
    	cout << "round " <<  itr << " finished with " << currentTime << " seconds elapsed." << endl;
		cout << std::scientific << std::setprecision(10) << obj1 << endl;
		if (delta < epsilon) break;
	}

	currentTime = omp_get_wtime( ) - startTime; // (clock() - t) / CLOCKS_PER_SEC
    cout << "The program finished with "<< currentTime << " seconds elapsed." << endl;

	return 0;
}
