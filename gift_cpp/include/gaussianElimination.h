#ifndef GELIMINATE_H
#define GELIMINATE_H
#include <stdint.h>
using namespace std;

void XOR(uint8_t **M, int r1, int r2, int COLSIZE);
void swapRow(uint8_t **M, int r1, int r2);
int findPivotRow(uint8_t **M, int currentRow, int currentCol, int &pivotCol, int ROWSIZE, int COLSIZE);
void GaussianElimination(uint8_t **M,int ROWSIZE,int COLSIZE);

#endif