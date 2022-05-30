#include <iostream>
#include "aes.h"
#include "aes_trail.h"
#include <vector>
#include <iomanip>
#include <numeric>
#include <string>

using namespace std;

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl

static void xorArray(uint8_t a[4][4], uint8_t b[4][4])
{
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) a[r][c] ^= b[r][c];
	}
}


void fixedValuesToPath(uint8_t alpha[20][4][4], uint8_t key[20][4][4], uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr)
{
	for (int n = 0; n < nr; n++)
	{
		uint8_t tmp[4][4];
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				diffStates[n][0][4*r+c] = alpha[n][r][c];
				diffStates[n][4][4*r+c] = alpha[n+1][r][c];
				tmp[r][c] = alpha[n+1][r][c];
			}
		}
		xorArray(tmp,key[n]);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][3][4*r+c] = tmp[r][c];}}
		invMC(tmp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][2][4*r+c] = tmp[r][c];}}
		invSR(tmp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][1][4*r+c] = tmp[r][c];}}
	}

	// Settling the fixedStates now
	// '1' corresponds to a value fixed to a specfic subset of values
	for (int n = 0; n < nr; n++)
	{
		// Substitution
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				if (diffStates[n][0][4*r+c] > 0 || fixedStates[n][0][4*r+c] > 0) {fixedStates[n][1][4*r+c] = 1;}
			}
		}
		// ShiftRows
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				fixedStates[n][2][4*r+c] = fixedStates[n][1][4*r+((4+c+r)%4)];
			}
		}
		// MixColumns
		for (int c = 0; c < 4; c++)
		{
			if (fixedStates[n][2][c] + fixedStates[n][2][4+c] + fixedStates[n][2][8+c] + fixedStates[n][2][12+c] == 4)
			{
				fixedStates[n][3][c] = 1;
				fixedStates[n][3][4+c] = 1;
				fixedStates[n][3][8+c] = 1;
				fixedStates[n][3][12+c] = 1;
			}
		}
		if (n < nr -1)
		{
			for (int r = 0; r < 4; r++)
			{
				for (int c = 0; c < 4; c++)
				{
					fixedStates[n][4][4*r+c] = fixedStates[n][3][4*r+c];
					fixedStates[n+1][0][4*r+c] = fixedStates[n][4][4*r+c];
				}
			}
		}
	}
}

void printState(uint8_t diffStates[20][5][16],uint8_t fixedStates[20][5][16], int nr)
{
	for (int n = 0; n < nr; n++)
	{
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(diffStates[n][0][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(diffStates[n][1][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(diffStates[n][2][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(diffStates[n][3][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(diffStates[n][4][4*r+c]) << ",";
			cout << endl;
		}
		cout << endl;
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(fixedStates[n][0][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(fixedStates[n][1][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(fixedStates[n][2][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(fixedStates[n][3][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(2) << hex << int(fixedStates[n][4][4*r+c]) << ",";
			cout << endl;
		}
		cout << endl;
		cout << "-------------------------------------------------------" << endl;
	}
}

bool contains(uint32_t arr[1000], uint32_t arr_length, uint32_t element)
{
	for (int i = 0; i < arr_length; i++)
	{
		if (element == arr[i]) return true;
	}
	return false;
}

vector<vector<uint32_t>> findConstraints(uint8_t diffStates[20][5][16],uint8_t fixedStates[20][5][16],int nr)
{
	uint32_t tup[16][1000] = {0}; //	0000|0000|0000(0000)|0000(0000) // round, pos, input, output
	uint32_t tupLength[16] = {0}; // record the length for each key

	uint32_t tupTemp[16][1000] = {0};
	uint32_t tupTempLength[16] = {0};
	vector<vector<uint32_t>> tupsConstraints;
	for (int n = 0; n < nr; n++)
	{
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				// output a constraint with current tuple
				if (fixedStates[n][0][4*r+c] > 0 && diffStates[n][0][4*r+c] > 0)
				{
					// put this into tupsConstraints
					vector<uint32_t> t;
					for (uint32_t i = 0; i < tupLength[4*r+c]; i++) t.push_back(tup[4*r+c][i]);
					t.push_back((n << (4+2*8)) + ((4*r+c) << (2*8)) + (int(diffStates[n][0][4*r+c]) << 8) + (diffStates[n][1][4*r+c]));
					tupsConstraints.push_back(t);
					tupLength[4*r+c] = 0;
				}
			}
		}
		// adding the tup
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				tup[4*r+c][tupLength[4*r+c]] = (n << (4+2*8)) + ((4*r+c) << (2*8)) + 
					(diffStates[n][0][4*r+c] << 8) + (diffStates[n][1][4*r+c]);
				tupLength[4*r+c]++;
			}
		}



		// applying SR MC to the key
		// row 0
		uint8_t row[4][4] = {{0,5,10,15},{1,6,11,12},{2,7,8,13},{3,4,9,14}};
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				for (uint32_t i = 0; i < tupLength[row[c][0]]; i++) 
				{
					if (contains(tupTemp[4*r+c],tupTempLength[4*r+c],tup[row[c][0]][i])) continue;
					tupTemp[4*r+c][tupTempLength[4*r+c]++] = tup[row[c][0]][i];
				}
				for (uint32_t i = 0; i < tupLength[row[c][1]]; i++) 
				{
					if (contains(tupTemp[4*r+c],tupTempLength[4*r+c],tup[row[c][1]][i])) continue;
					tupTemp[4*r+c][tupTempLength[4*r+c]++] = tup[row[c][1]][i];
				}
				for (uint32_t i = 0; i < tupLength[row[c][2]]; i++) 
				{
					if (contains(tupTemp[4*r+c],tupTempLength[4*r+c],tup[row[c][2]][i])) continue;
					tupTemp[4*r+c][tupTempLength[4*r+c]++] = tup[row[c][2]][i];
				}
				for (uint32_t i = 0; i < tupLength[row[c][3]]; i++) 
				{
					if (contains(tupTemp[4*r+c],tupTempLength[4*r+c],tup[row[c][3]][i])) continue;
					tupTemp[4*r+c][tupTempLength[4*r+c]++] = tup[row[c][3]][i];					
				}
			}
		}
		// put it back to tup from tupTemp for those that are 1
		// erase those with 0
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				if (fixedStates[n][4][4*r+c] == 0) tupLength[4*r+c] = 0;
				else 
				{
					for (int i = 0; i < tupTempLength[4*r+c]; i++) tup[4*r+c][i] = tupTemp[4*r+c][i];
					tupLength[4*r+c] = tupTempLength[4*r+c];
				}
				tupTempLength[4*r+c] = 0;
			}
		}
	}
	return tupsConstraints;
}

void XOR(uint32_t val1[256],uint32_t val2[256])
{
	int sum1 = 0, sum2 = 0;
	for (int i = 0; i < 256; i++)
	{
		sum1 += val1[i];
		sum2 += val2[i];
	}
	if (sum2 == 0) return;
	if (sum1 == 0)
	{
		for (int i = 0; i < 256; i++) val1[i] = val2[i];
		return;
	}
	uint32_t res[256] = {0};
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < 256; j++) 
		{
			if (val1[i] * val2[j] > 0) res[i^j] += 1; // boolean only!
			else res[i^j] += 0;
			// res[i^j] += val1[i] * val2[j];
		}
	}
	for (int i = 0; i < 256; i++)
	{
		val1[i] = res[i];
	}
}

void reduceKey(uint32_t *val, uint32_t valLength)
{
	int counter = 0;
	int a[2]; a[0] = 1; a[1] = 1;
	for (int i = 0; i < (1<<(valLength*8)); i++)
	{
		if (val[i] > 0)
		{
			if (counter == 0) a[counter++] = val[i];
			else 
			{	
				a[1] = val[i];
				a[0] = gcd(a[0],a[1]);
			}
		}
	}
	for (int i = 0; i < (1<<(valLength*8)); i++) val[i] /= a[0];
}

bool isLinear(vector<uint32_t> constraint)
{
	uint32_t ANDval = (1 << 8) - 1;
	for (uint32_t i = 0; i < constraint.size();i++) if ((constraint[i] & ANDval) == 0) return false;
	return true;
}

void mul2Array(uint32_t arr[256])
{
	uint32_t tmp[256] = {0};
	for (int i = 0; i < 256; i++) tmp[mul2(i)] = arr[i];
	for (int i = 0; i < 256; i++) arr[i] = tmp[i];
}

void mul3Array(uint32_t arr[256])
{
	uint32_t tmp[256] = {0};
	for (int i = 0; i < 256; i++) tmp[mul3(i)] = arr[i];
	for (int i = 0; i < 256; i++) arr[i] = tmp[i];
}


vector<int> resolveLinear(vector<uint32_t> constraint, uint32_t* &rightPartValues)
{
	uint8_t MCArray[4][4] = {{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};
	vector<int> keys;
	// the last int in the constraint is the one with the right part (use xddt)
	int t = constraint[constraint.size()-1];
	uint32_t ANDval = 255;
	uint32_t xddt[256] = {0};
	computeXDDT((t >> 8) & ANDval,t & ANDval,xddt);

	for (int i = 0; i < 256; i++) {rightPartValues[i] = xddt[i];}
	keys.push_back(constraint[constraint.size()-1]);
	for (uint32_t i = 0; i < constraint.size()-1; i++)
	{
		uint32_t yddt[256] = {0};
		computeYDDT((constraint[i] >> 8) & ANDval,constraint[i] & ANDval,yddt);
		if (MCArray[((t >> 16) & 0xf) / 4][((constraint[i] >> 16) & 0xf) / 4] == 2) mul2Array(yddt);
		else if (MCArray[((t >> 16) & 0xf) / 4][((constraint[i] >> 16) & 0xf) / 4] == 3) mul3Array(yddt);
		XOR(rightPartValues,yddt);
		reduceKey(rightPartValues,1);
	}
	return keys;
}



bool isZero(uint32_t *array, uint32_t length)
{
	for (uint32_t i = 0; i < length; i++)
	{
		if (array[i] != 0) return false;
	}
	return true;
}

bool isIn(vector<uint32_t> vec, int x)
{
	for (int i = 0; i < vec.size(); i++)
	{
		if (vec[i] == x) return true;
	}
	return false;
}


void Overlap(uint32_t arr1[256], uint32_t arr2[256])
{
	for (int i = 0; i < 256; i++) arr1[i] = arr1[i] * arr2[i];
}

int findIndex(vector<uint32_t> roundPos, int element)
{
	for (uint32_t i = 0; i < roundPos.size(); i++)
	{
		if (roundPos[i] == element) return i;
	}
	return -1;
}

void substArray(uint32_t arr[256])
{
	uint32_t tmp[256] = {0};
	for (int i = 0; i < 256; i++) tmp[getSbox(i)] = arr[i];
	for (int i = 0; i < 256; i++) arr[i] = tmp[i];
}

void propagateLinear(vector<vector<int>> &linearKeys, uint32_t **linearValues, uint32_t valuesLength, int TK_NUM, uint32_t lastLinearKeyRound, uint32_t keySize = 128)
{
	// need to check 1. is it all uniformed?
	// 2. is there a collision already?
	int required;
	int nextCol = 0;
	uint32_t nr, col;
	int N = keySize / 32;
	vector<uint32_t> roundPos;


	for (int i = 0; i < valuesLength; i++) linearValues[i] = linearValues[i];
	for (int i = 0; i < linearKeys.size(); i++) roundPos.push_back(linearKeys[i][0]>>16);
	if (keySize == 128) required = 11 * 4;
	if (keySize == 192) required = 13 * 4;
	if (keySize == 256) required = 15 * 4;
	
	bool zeroComputed[32] = {0};
	uint32_t *valuesInComputation; valuesInComputation = new uint32_t[2 * 16 * 256](); // 2 rounds * 16 pos * 256 possibilities // this stores the keys that are still computing
	uint32_t *valuesCompleted; valuesCompleted = new uint32_t[2 * 16 * 256]; // 2 rounds * 16 pos * 256 possibilities // this stores the keys that are still computing
	bool zeroComputing[32] = {0};
	// only for AES-128
	while (nextCol < (lastLinearKeyRound *4))
	{
		nr = nextCol / 4; col = nextCol % 4;
		
		for (int n = 0; n < 2; n++) 
		{
			for (int r = 0; r < 4; r++)
			{
				for (int c = 0; c < 4; c++)
				{
					if (isZero(&valuesInComputation[n * (16*256) + (4*r+c) * 256],256)) zeroComputing[n*16+4*r+c] = 1;
					else zeroComputing[n*16+4*r+c] = 0;
				}
			}
		}
		// this is the key schedule to propagate the values
		if (nextCol < N) // insert key into roundValues when it's initializing 
		{

		}
		else if ((nextCol >= N) && ((nextCol % N) == 0)) // Rotate nr-1, Sbox nr-1, XOR RCON[nextCol/N] XOR nr-N
		{
			// start retrieving values from the previous valuesCompleted instead of the linearValues
			// for first row
			
			if ((zeroComputed[((nextCol % N) / 4) * 16 + 4 * 0 + col] == 0) && (zeroComputed[((nextCol % N) / 4) * 16 + 4 * (0+1) + (col-1) % N] == 0))
			{
				XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*0+col)*256],&valuesCompleted[((nextCol % N)/ 4) * (256*16) + (4*1 + ((col-1) % N)) * 256]);
				substArray(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*0+col)*256]);
				XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*0+col)*256],&valuesCompleted[((nextCol % N)/ 4) * (256*16) + (4*0 + col) * 256]);
				uint32_t tmp[256] = {0}; tmp[getConstants(nr-1)] = 1;
				XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*0+col)*256],tmp);
			}
			for (int r = 1; r < 4; r++)
			{

				if ((zeroComputed[((nextCol % N) / 4) * 16 + 4 * r + col] == 0) && (zeroComputed[((nextCol % N) / 4) * 16 + 4 * ((r+1)%4) + ((col-1)%N)] == 0))
				{
					XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256],&valuesCompleted[((nextCol % N)/ 4) * (256*16) + (4* ((r+1)%4) + ((col-1)%N)) * 256]);
					substArray(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256]);
					XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256],&valuesCompleted[((nextCol % N)/ 4) * (256*16) + (4* r + col) * 256]);
				}
			}
		}
		else if ((nextCol > N) && (N > 6) && ((nextCol % N) == 4)) // only for AES-256 Sbox XOR for mod 4
		{
			for (int r = 0; r < 4; r++)
			{
				if ((zeroComputed[((nextCol % N) / 4) * 16 + 4 * r + col] == 0) && (zeroComputing[((nextCol % N) / 4) * 16 + 4 * r + col-1] == 0))
				{
					XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256],&valuesInComputation[((nextCol % N)/ 4) * (256*16) + (4* r + col-1) * 256]);
					substArray(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256]);
					XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256],&valuesCompleted[((nextCol % N)/ 4) * (256*16) + (4* r + col) * 256]);
				}
			}
		}
		else // Just XOR only
		{
			for (int r = 0; r < 4; r++)
			{
				if ((zeroComputed[((nextCol % N) / 4) * 16 + 4 * r + col] == 0) && (zeroComputing[((nextCol % N) / 4) * 16 + 4 * r + col-1] == 0))
				{
					XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256],&valuesInComputation[((nextCol % N)/ 4) * (256*16) + (4* r + col-1) * 256]);
					XOR(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256],&valuesCompleted[((nextCol % N)/ 4) * (256*16) + (4* r + col) * 256]);
				}

			}
		}
		// now, if all have been computed, we can update and shift it to valuesCompleted
		if ((nextCol % N) == N-1)
		{
			
			for (int n = 0; n < 2; n++) 
			{
				for (int r = 0; r < 4; r++)
				{
					for (int c = 0; c < 4; c++)
					{
						if (isZero(&valuesInComputation[n * (16*256) + (4*r+c) * 256],256)) zeroComputing[n*16+4*r+c] = 1;
						else zeroComputing[n*16+4*r+c] = 0;
					}
				}
			}

			// combining propagated values with new restrictions from trail
			for (int n = N - 1; n >= 0 ; n--)
			{
				nr = (nextCol - n) / 4;
				col = (nextCol - n) % 4;
				for (int r = 0; r < 4; r++)
				{
					if (isIn(roundPos,(nr << 4) + 4*r + col))
					{
						int index = findIndex(roundPos,(nr << 4) + 4*r + col);
						if (zeroComputing[((nextCol % N) / 4) * (16) + (4*r+col)] == 1){ for (int i = 0; i < 256; i++) valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256 + i] = linearValues[index][i];}
						else Overlap(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+col)*256],linearValues[index]);
					}
				}
			}
			// check if there is already an incompatiblility
			for (int n = 0; n < 2; n++)
			{
				for (int r = 0; r < 4; r++)
				{
					for (int c = 0; c < 4; c++)
					{
						if (zeroComputing[n * 16 + 4*r + c] == 0)
						{
							if (isZero(&valuesInComputation[((nextCol % N) / 4) * (16*256) + (4*r+c)*256],256) == 1)
							{
								cout << "incompatability occurs" << endl;
								printout(n);
								printout(r);
								printout(c);
								printout(nextCol);
							}
						}
					}
				}
			}

			// Update from computation to completed
			// refreshing computation
			for (int i = 0; i < 256 * 16 * 2; i++) 
			{
				valuesCompleted[i] = valuesInComputation[i];
				valuesInComputation[i] = 0;
			}
			for (int n = 0; n < 2; n++) 
			{
				for (int r = 0; r < 4; r++)
				{
					for (int c = 0; c < 4; c++)
					{
						reduceKey(&valuesCompleted[n * (16*256) + (4*r+c) * 256],1);
					}
				}
			}

			// update the zeroComputed
			for (int n = 0; n < 2; n++) 
			{
				for (int r = 0; r < 4; r++)
				{
					for (int c = 0; c < 4; c++)
					{
						if (isZero(&valuesCompleted[n * (16*256) + (4*r+c) * 256],256)) zeroComputed[n*16+4*r+c] = 1;
						else zeroComputed[n*16+4*r+c] = 0;
					}
				}
			}

		}
		nextCol++;
		printout(nextCol);
		
	}
}


int main()
{
	rand_generator.seed(time(0));
	uint8_t alpha[20][4][4] = {0};
	uint8_t key_diff[3][4][4] = {0};
	uint8_t diffStates[20][5][16];
	uint8_t fixedStates[20][5][16] = {0};
	uint8_t actualKeyDiff[20][4][4] = {0};
	uint32_t nr;
	int TK_NUM;
	for (int q = 0; q < 10000; q++)
	{
	// AES trails
	// generateBestTrail(alpha,actualKeyDiff,0x10,3,nr,TK_NUM);
	// TK1_2009_241(alpha,actualKeyDiff,nr,TK_NUM); // AES trail
	justACrazyTrail(alpha,actualKeyDiff,nr,TK_NUM);



	fixedValuesToPath(alpha,actualKeyDiff,diffStates,fixedStates,nr);
	printState(diffStates,fixedStates,nr);

	vector<vector<uint32_t>> constraints = findConstraints(diffStates,fixedStates,nr);
	
	// just for printing
	uint32_t ANDval = (1 << 8) - 1;
	for (uint32_t i = 0; i < constraints.size(); i++)
	{
		for (uint32_t j = 0; j < constraints[i].size(); j++)
		{
			cout << hex << "(" << (constraints[i][j] >> (4+2*8)) << "," << ((constraints[i][j] >> (2*8)) & 0xf) << ","
			<< ((constraints[i][j] >> 8) & ANDval) << "," << (constraints[i][j] & ANDval) << ") ";
		}
		cout << isLinear(constraints[i]);
		cout << endl;
	}


	vector<vector<uint32_t>> linearConstraints;
	vector<vector<uint32_t>> nonLinearConstraints;
	for (uint32_t i = 0; i < constraints.size(); i++)
	{
		if (isLinear(constraints[i])) linearConstraints.push_back(constraints[i]);
		else if (!isLinear(constraints[i])) nonLinearConstraints.push_back(constraints[i]);
		else {cout << "error!" << endl;}
	}
	uint32_t **linearValues; linearValues = new uint32_t*[linearConstraints.size()]();
	uint32_t *linearValuesLength; linearValuesLength = new uint32_t[linearConstraints.size()]();
	vector<vector<int>> linearKeys;
	cout << "Allocated " << linearConstraints.size() << " for linearValues" << endl;
	// check linear propagation 
	// check 1. perm schedule
	// check 2. multiple keys how to propagate?
	for (uint32_t i = 0; i < linearConstraints.size(); i++)
	{
		uint32_t *values; values = new uint32_t[(1<<8)]();
		uint32_t valuesLength = 1;
		cout << "Resolving linear constraints[" << i << "]...";
		vector<int> keys = resolveLinear(linearConstraints[i],values);
		linearValues[i] = values;
		linearKeys.push_back(keys);
		cout << "Completed" << endl;
	}
	uint32_t lastLinearKeyRound = linearKeys[linearKeys.size()-1][linearKeys[linearKeys.size()-1].size()-1] >> 20;
	propagateLinear(linearKeys,linearValues,linearConstraints.size(),TK_NUM,lastLinearKeyRound,128);
	}
	return 0;
}