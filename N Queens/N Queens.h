//Find all solutions of N Queens
#ifndef N_queens
#include <cmath>
#include<iostream>
#include<stdio.h>
using namespace std;
class Board
{
    int size;   //size is the numbers of queens
    int **chess_board; 

    bool is_this_position_safe(int row, int col) //judge whether the (row, col) position is safe
    {
		int i, j;
		for (i = 0; i<size; i++)
			for (j = 0; j < size; j++) 
			{			
				if (i == row || j == col||abs(i - row) == abs(j - col)) // if there is a chess at the same row or column or diagonal
				{
					if (chess_board[i][j] == 1)
					{
						return false;
					}
				}
			}
		return true;		
    }
    

    void initialize(int n)    // initializes the (n x n) chessboard
    {
        size = n;
		int i, j;
        chess_board=new int*[size];//allocate space to chessboard
		for (i = 0; i < size; i++)
			chess_board[i] = new int[size];
		for (i = 0; i < size; i++)
			for (j = 0; j < size; j++)
				chess_board[i][j] = 0;
    }
    
    
    void print_board() //print the chessboard
    {
    	static int count = 1;
     	printf("\nCase %d:\n", count++);
		for (int i = 0; i < size; i++) 
		{
			for (int j = 0; j < size; j++) 
			{
				if(chess_board[i][j] == 1)
				 cout <<"Q "; 
				else 
					cout <<"_ "; 			
			}
			cout << endl;
		}
	}
    
 
    
    // private member function: recursive backtracking
	bool solve(int col)
	{
		if (col >= size)
			print_board();
			//return true;   //if use this line, we can find only one solution
		else
		{
			for (int i = 0; i < size; i++)
			{
				if (is_this_position_safe(i, col)) 
				{
					chess_board[i][col] = 1;
					if (solve(col + 1))
					{
						return true;
					}
					else
						chess_board[i][col] = 0;	
				}
			}
		}
		return false;
	}
    

public:
    // Solves the n-Queens problem by (recursive) backtracking
    void nQueens(int n)
    {

        initialize(n);
		solve(0);	
        if(n<4)
            cout << "There is no solution to the " << n << "-Queens Problem" << endl;
		if(n>=4)
		cout <<"\n----------\nThis is "<< n << "-Queens Problem Solution" << endl;
    }
};
#endif
