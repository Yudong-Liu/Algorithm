// N Queens Problem via (Backtracking, which is implemented by) Recursion 
// Yudong Liu's Assignment 1

#include <iostream>
//#include "N_Queens_part1_lyd.h"
#include "N Queens.h"
using namespace std;

int main (int argc, char * const argv[]) 
{
    Board x;   
    int board_size;
	sscanf(argv[1], "%d", &board_size);  
	//cin >> board_size;
	x.nQueens(board_size);
    return 0;
}
