// #include <iostream>
#include <stdint.h>
using namespace std;

void XOR(uint8_t **M, int r1, int r2, int COLSIZE)
{
	for (int c = 0; c < COLSIZE; c++) M[r1][c] ^= M[r2][c];
}

void swapRow(uint8_t **M, int r1, int r2)
{
	uint8_t *t = M[r1];
	M[r1] = M[r2];
	M[r2] = t;
}

int findPivotRow(uint8_t **M, int currentRow, int currentCol, int &pivotCol, int ROWSIZE, int COLSIZE)
{
	for (int c = currentCol; c < COLSIZE; c++)
	{
		for (int r = currentRow; r < ROWSIZE; r++)
		{
			if (M[r][c] == 1)
			{
				pivotCol = c;
				return r;
			}
		}
	}
	return -1; // return -1 when there is nothing else to find
}


void GaussianElimination(uint8_t **M,int ROWSIZE,int COLSIZE)
{
	int currentRow = 0;
	int currentCol = 0;
	while (true) // perform REF
	{
		int pivotCol; 
		int pivotRow = findPivotRow(M,currentRow,currentCol,pivotCol,ROWSIZE,COLSIZE);
		if (pivotRow == -1) break;
		swapRow(M,currentRow,pivotRow);
		for (int r = currentRow+1; r < ROWSIZE; r++)
		{
			if (M[r][pivotCol] == 1) XOR(M,r,currentRow,COLSIZE); // r <-- r XOR currentRow
		}
		currentRow += 1;
		currentCol = pivotCol + 1;
	}
	return;
}



// int main()
// {
// 	uint8_t **M; M = new uint8_t*[4];
// 	for (int i = 0; i < 4; i++) M[i] = new uint8_t[6]();
// 	M[0][0] = 0; M[0][1] = 1; M[0][2] = 1; M[0][3] = 1; M[0][4] = 1; ; M[0][5] = 0;
// 	M[1][0] = 1; M[1][1] = 0; M[1][2] = 1; M[1][3] = 0; M[1][4] = 1; ; M[1][5] = 0;
// 	M[2][0] = 1; M[2][1] = 1; M[2][2] = 1; M[2][3] = 1; M[2][4] = 1; ; M[2][5] = 1;
// 	M[3][0] = 1; M[3][1] = 1; M[3][2] = 0; M[3][3] = 1; M[3][4] = 1; ; M[3][5] = 1;
// 	for (int r = 0; r < 4; r++)
// 	{
// 		for (int c = 0; c < 6; c++) cout << int(M[r][c]) << " ";
// 		cout << endl;
// 	}
// 	cout << endl;
// 	GaussianElimination(M,4,6);

// 	for (int r = 0; r < 4; r++)
// 	{
// 		for (int c = 0; c < 6; c++) cout << int(M[r][c]) << " ";
// 		cout << endl;
// 	}
// }