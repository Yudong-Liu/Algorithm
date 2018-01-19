//
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include<time.h>

#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"

using namespace std;


class Repeat_Squaring {
	Matrix A;	
	int no_rows;
	int exponent;

	//matrix pow function
	Matrix mat_pow(Matrix A,int exponent) {
		if (exponent == 0) {
			Matrix One(no_rows,no_rows);//Generate a unit matrix
			for (int i = 1; i <= no_rows; i++) {
				for (int j = 1; j <= no_rows; j++) {
					if(i==j)
						One(i, j) = 1;
					else
						One(i, j) = 0;
				}
			}
			return One;
		}
			
		else if (exponent % 2 == 1)
			return A*mat_pow(A*A, (exponent - 1) / 2);
		else
			return mat_pow(A*A, exponent / 2);
	}

	//Repeated Squaring 
	Matrix repeated_squaring(Matrix A, int no_rows, int exponent) {
		Matrix C;		
		C = mat_pow(A, exponent);
		//cout << setw(10) << C << endl;
		return C;
	}

	//Direct Multiplication
	Matrix direct_mul(Matrix A, int no_rows, int exponent) {
		if (exponent == 0) {
			Matrix One(no_rows, no_rows);//Generate a unit matrix
			for (int i = 1; i <= no_rows; i++) {
				for (int j = 1; j <= no_rows; j++) {
					if (i == j)
						One(i, j) = 1;
					else
						One(i, j) = 0;
				}
			}
			return One;
		}
		
		for (int i = 0; i < exponent; i++) {
			return A*direct_mul(A,no_rows,exponent-1);
		}

	}

	double get_uniform()
	{
		return (((double)rand()) / (RAND_MAX));
	}
	

public:
	//generate a random matrix
	void initialize(int no_rows, Matrix &A) {
		for (int i = 1; i <= no_rows; i++) {
			for (int j = 1; j <= no_rows; j++) {
				A(i, j) = 5 * (get_uniform() - 0.5);
			}
		}
	}
	void calculate(int argc, char * const argv[]) {
		sscanf_s(argv[1], "%d", &no_rows);
		sscanf_s(argv[2], "%d", &exponent);
		cout << "The number of rows/columns in the square matrix is :" << no_rows << endl;
		cout << "The exponent is:" << exponent << endl;
		Matrix A(no_rows,no_rows);
		Matrix c1(no_rows, no_rows);
		Matrix c2(no_rows, no_rows);
		initialize(no_rows, A);
		double diff1,diff2;
		clock_t t1,t2,t3;
		t1= clock();		
		c1=repeated_squaring(A, no_rows, exponent);
		t2 = clock();
		c2=direct_mul(A, no_rows, exponent);
		t3 = clock();
		diff1 = ((double)t2 - (double)t1);
		diff2= ((double)t3 - (double)t2);
		cout << "Repeated Squaring Result:" << endl;
		cout << "It took " << diff1 / CLOCKS_PER_SEC << " seconds to complete" << endl;
		cout << "Direct Multiplication Result:" << endl;
		cout << "It took " << diff2 / CLOCKS_PER_SEC << " seconds to complete" << endl;

		//Test the result of matrix calculation
		cout << "initial matrix:" << endl;
		cout<< setw(10) << A << endl;
		cout << endl;
		cout << "RS:" << endl;
		cout << setw(10) << c1 << endl;
		cout << endl;
		cout << "DM:" << endl;
		cout << setw(10) << c2 << endl;
	
	}
	
};


	int main(int argc, char* argv[])
	{
		Repeat_Squaring A;
		
		A.calculate(argc, argv);
		
		return 0;
	}




