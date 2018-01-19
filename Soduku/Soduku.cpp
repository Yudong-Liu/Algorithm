//Assignment2, finished by Yudong Liu

#include <iostream>
#include "Soduku.h"
#include<cmath>
using namespace std;

int main (int argc, char * const argv[]) 
{
	Sudoku x;
	int choice;

	x.read_puzzle(argc, argv);
	x.print_puzzle();

	cout << "If you want to find one solution, input 1; if you want to find all solutions, input 2" << endl;
	cout << "pleased input your choice:";
	cin >> choice;
	if (choice == 1) {	//find one solution
		x.Solve(0, 0);
		x.print_puzzle();
	}
	else if (choice == 2)	//find all solutions
		x.all_solutions();
	else
		cout << "There is an error!!!";
	
	system("pause");
    return 0;
}
