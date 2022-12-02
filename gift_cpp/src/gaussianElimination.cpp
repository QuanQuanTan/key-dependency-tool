#include "gaussianElimination.h"

void XOR(uint8_t **M, int r1, int r2, int COLSIZE){
	for (int c = 0; c < COLSIZE; c++) M[r1][c] ^= M[r2][c];
}

void swapRow(uint8_t **M, int r1, int r2){
	uint8_t *t = M[r1];
	M[r1] = M[r2];
	M[r2] = t;
}

int findPivotRow(uint8_t **M, int currentRow, int currentCol, int &pivotCol, int ROWSIZE, int COLSIZE){
	for (int c = currentCol; c < COLSIZE; c++){
		for (int r = currentRow; r < ROWSIZE; r++){
			if (M[r][c] == 1){
				pivotCol = c;
				return r;
			}
		}
	}
	return -1; // return -1 when there is nothing else to find
}

void GaussianElimination(uint8_t **M,int ROWSIZE,int COLSIZE){
	int currentRow = 0;
	int currentCol = 0;
	while (true){
		// perform REF
		int pivotCol; 
		int pivotRow = findPivotRow(M,currentRow,currentCol,pivotCol,ROWSIZE,COLSIZE);
		if (pivotRow == -1) break;
		swapRow(M,currentRow,pivotRow);
		for (int r = currentRow+1; r < ROWSIZE; r++){
			if (M[r][pivotCol] == 1) XOR(M,r,currentRow,COLSIZE); // r <-- r XOR currentRow
		}
		currentRow += 1;
		currentCol = pivotCol + 1;
	}
	return;
}

