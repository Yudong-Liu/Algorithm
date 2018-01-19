
#ifndef sudoku
#define sudoku

#include <iostream>
#include <vector>
#include <fstream>
using std::vector;
using namespace std;
class Sudoku 
{ 
	// Private
	int puzzle[9][9];
	int count = 0;
	// Private member function that checks if the named row is valid
	bool row_valid(int row, int col)
	{
		int j;
		for (j = 0; j < 9; j++) 
			if (puzzle[row][j] == puzzle[row][col] && j != col)
				return false;			
		return true;
	}
	
	// Private member function that checks if the named column is valid
	bool col_valid(int row, int col)
	{
		int i;
		for ( i = 0; i < 9; i++) {
			if (puzzle[i][col] == puzzle[row][col] && i != row)
				return false;
		}			
		return true;
	}
	
	// Private member function that checks if the named 3x3 block is valid
	bool block_valid(int row, int col)
	{
		// check 3 x 3 block validity 
		int block_row=0, block_col=0;
		block_row =3*( row / 3 );
		block_col = 3* (col / 3) ;
		int i, j ;
		for (int i = block_row; i < block_row + 3; i++) {
			for (int j = block_col; j < block_col + 3; j++) {
				if (puzzle[i][j] == puzzle[row][col] && i != row && j != col)
					return false;
			}				
		}						
		return true;
	}

	bool valid(int row, int col) {
		return row_valid(row, col) && col_valid(row, col) && block_valid(row, col);
	}

	bool check_zero() {
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) {
				if (0==puzzle[i][j])
					return false;
			}
		}
	}
public:
	void read_puzzle(int argc, char * const argv[])
	{

		fstream readfile;
		int num;
		//char filename[10];
		//cout << "Input the filename:" << endl;
		//cin >> filename;
		readfile.open(argv[1], ios::in);
		if (!readfile.eof()) {
			for (int i = 0; i < 9; i++)
				for (int j = 0; j < 9; j++) {
					readfile >> num;
					puzzle[i][j] = num;
				}					
		}			
		readfile.close();
	}
	
	// Public member function that prints the puzzle when called
	void print_puzzle()
	{	

		if (count == 0) {		
			std::cout << std::endl << "Initial Puzzle"<< std::endl;			
		}		
		else {			
			std::cout << std::endl << "Solution " << count << std::endl;
		}
		count++;
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				// check if we have a legitimate integer between 1 and 9
				if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
				{
					// printing initial value of the puzzle with some formatting
					std::cout << puzzle[i][j] << " ";
				}
				else {
					// printing initial value of the puzzle with some formatting
					std::cout << "X ";
				}
			}
			std::cout << std::endl;
		}
	}

	bool Solve(int row, int col)
	{
		int  k;
		if (row > 8 ) {  //the recursion is finished
			for (int i = row; i < 9; i++) {
				for (int j = col; j < 9; j++) {
					if (0 == puzzle[i][j])
						return true;
				}
			}
		}
		else if (puzzle[row][col] != 0) {//recurse
			if (8==col)
				Solve(row + 1, 0);
			else
				Solve(row, col + 1);
		}
		else if (puzzle[row][col] == 0) {  //try different value
			for (k = 1; k < 10; k++) {
				puzzle[row][col] = k;
				if (valid(row, col) && Solve(row, col)) {
					return true;
				}
			}			
			puzzle[row][col] = 0;
			return false;						
		}
	}
	void all_solutions()
	{
		int  i, j, k;
	 if (check_zero()) { //if there is no zero, end recursion
			print_puzzle();
			return;
		}
		for (i = 0; i<9; i++)
			for (j = 0; j<9; j++)
				if (puzzle[i][j] == 0) {
					for (k = 1; k <10; k++) {
						puzzle[i][j] = k;
						if (valid(i, j)) {
							all_solutions();
						}
					}
					puzzle[i][j] = 0;
					return;
				}
	}

};


#endif