// including 8 bit skinny in here
#include "skinny.h"
#include "trails.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <bits/stdc++.h>
#include "SAT.h"
#include <cryptominisat5/cryptominisat.h>
#include <iomanip>
#include <cmath>
using namespace std;

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl

#ifndef SIZE
#define SIZE 4;
#endif


void fixedValuesToPath(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr)
{
	uint8_t key_diff_temp[4][4][4] = {0};
	for (int i = 0; i < 4; i++) { for (int r = 0; r < 4; r++) { for (int c = 0; c < 4; c++) key_diff_temp[i][r][c] = key_diff[i][r][c];}}
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
		invMC(tmp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][3][4*r+c] = tmp[r][c];}}
		invSR(tmp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][2][4*r+c] = tmp[r][c];}}
		AddRoundTweakey(tmp,key_diff_temp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][1][4*r+c] = tmp[r][c];}}
		if (SIZE == 4) key_schedule_round_function(key_diff_temp);
		else key_schedule_round_function8(key_diff_temp);
	}

	// Settling the fixedStates now
	// '1' corresponds to a value fixed to a specfic subset of values
	for (int n = 0; n < nr; n++)
	{
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				if (diffStates[n][0][4*r+c] > 0 || fixedStates[n][0][4*r+c] > 0) {fixedStates[n][1][4*r+c] = 1; fixedStates[n][2][4*r+c] = 1;}
			}
		}

		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				fixedStates[n][3][4*r+c] = fixedStates[n][2][4*r+((4+c-r)%4)];
			}
		}
		for (int c = 0; c < 4; c++)
		{
			if (fixedStates[n][3][c] > 0 && fixedStates[n][3][8+c] > 0 && fixedStates[n][3][12+c] > 0) fixedStates[n][4][c] = 1;
			if (fixedStates[n][3][c] > 0) fixedStates[n][4][4+c] = 1;
			if (fixedStates[n][3][4+c] > 0 && fixedStates[n][3][8+c] > 0) fixedStates[n][4][8+c] = 1;
			if (fixedStates[n][3][c] > 0 && fixedStates[n][3][8+c] > 0) fixedStates[n][4][12+c] = 1;
		}
		if (n < nr -1)
		{
			for (int r = 0; r < 4; r++)
			{
				for (int c = 0; c < 4; c++)
				{
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
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][0][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][1][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][2][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][3][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][4][4*r+c]) << ",";
			cout << endl;
		}
		cout << endl;
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][0][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][1][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][2][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][3][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][4][4*r+c]) << ",";
			cout << endl;
		}
		cout << endl;
		cout << "-------------------------------------------------------" << endl;
	}
}

vector<vector<uint32_t>> findConstraints(uint8_t diffStates[20][5][16],uint8_t fixedStates[20][5][16],int nr)
{
	uint32_t tup[16][50] = {0}; //	0000|0000|0000(0000)|0000(0000) // round, pos, input, output
	uint32_t tupLength[16] = {0}; // record the length for each key

	uint32_t tupTemp[16][50] = {0};
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
					t.push_back((n << (4+2*SIZE)) + ((4*r+c) << (2*SIZE)) + (int(diffStates[n][0][4*r+c]) << SIZE) + (diffStates[n][1][4*r+c]));
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
				tup[4*r+c][tupLength[4*r+c]] = (n << (4+2*SIZE)) + ((4*r+c) << (2*SIZE)) + 
					(diffStates[n][0][4*r+c] << SIZE) + (diffStates[n][1][4*r+c]);
				tupLength[4*r+c]++;
			}
		}



		// applying SR MC to the key
		// row 0
		for (uint32_t i = 0; i < tupLength[0]; i++) tupTemp[0][tupTempLength[0]++] = tup[0][i];
		for (uint32_t i = 0; i < tupLength[10]; i++) tupTemp[0][tupTempLength[0]++] = tup[10][i];
		for (uint32_t i = 0; i < tupLength[13]; i++) tupTemp[0][tupTempLength[0]++] = tup[13][i];

		for (uint32_t i = 0; i < tupLength[1]; i++) tupTemp[1][tupTempLength[1]++] = tup[1][i];
		for (uint32_t i = 0; i < tupLength[11]; i++) tupTemp[1][tupTempLength[1]++] = tup[11][i];
		for (uint32_t i = 0; i < tupLength[14]; i++) tupTemp[1][tupTempLength[1]++] = tup[14][i];

		for (uint32_t i = 0; i < tupLength[2]; i++) tupTemp[2][tupTempLength[2]++] = tup[2][i];
		for (uint32_t i = 0; i < tupLength[8]; i++) tupTemp[2][tupTempLength[2]++] = tup[8][i];
		for (uint32_t i = 0; i < tupLength[15]; i++) tupTemp[2][tupTempLength[2]++] = tup[15][i];

		for (uint32_t i = 0; i < tupLength[3]; i++) tupTemp[3][tupTempLength[3]++] = tup[3][i];
		for (uint32_t i = 0; i < tupLength[9]; i++) tupTemp[3][tupTempLength[3]++] = tup[9][i];
		for (uint32_t i = 0; i < tupLength[12]; i++) tupTemp[3][tupTempLength[3]++] = tup[12][i];

		// row 1
 		for (uint32_t i = 0; i < tupLength[0]; i++) tupTemp[4][tupTempLength[4]++] = tup[0][i];	
		for (uint32_t i = 0; i < tupLength[1]; i++) tupTemp[5][tupTempLength[5]++] = tup[1][i];	
		for (uint32_t i = 0; i < tupLength[2]; i++) tupTemp[6][tupTempLength[6]++] = tup[2][i];	
		for (uint32_t i = 0; i < tupLength[3]; i++) tupTemp[7][tupTempLength[7]++] = tup[3][i];	
		
		// row 2
		for (uint32_t i = 0; i < tupLength[7]; i++) tupTemp[8][tupTempLength[8]++] = tup[7][i];
		for (uint32_t i = 0; i < tupLength[10]; i++) tupTemp[8][tupTempLength[8]++] = tup[10][i];

		for (uint32_t i = 0; i < tupLength[4]; i++) tupTemp[9][tupTempLength[9]++] = tup[4][i];
		for (uint32_t i = 0; i < tupLength[11]; i++) tupTemp[9][tupTempLength[9]++] = tup[11][i];

		for (uint32_t i = 0; i < tupLength[5]; i++) tupTemp[10][tupTempLength[10]++] = tup[5][i];
		for (uint32_t i = 0; i < tupLength[8]; i++) tupTemp[10][tupTempLength[10]++] = tup[8][i];

		for (uint32_t i = 0; i < tupLength[6]; i++) tupTemp[11][tupTempLength[11]++] = tup[6][i];
		for (uint32_t i = 0; i < tupLength[9]; i++) tupTemp[11][tupTempLength[11]++] = tup[9][i];

		// row 3
		for (uint32_t i = 0; i < tupLength[0]; i++) tupTemp[12][tupTempLength[12]++] = tup[0][i];
		for (uint32_t i = 0; i < tupLength[10]; i++) tupTemp[12][tupTempLength[12]++] = tup[10][i];

		for (uint32_t i = 0; i < tupLength[1]; i++) tupTemp[13][tupTempLength[13]++] = tup[1][i];
		for (uint32_t i = 0; i < tupLength[11]; i++) tupTemp[13][tupTempLength[13]++] = tup[11][i];

		for (uint32_t i = 0; i < tupLength[2]; i++) tupTemp[14][tupTempLength[14]++] = tup[2][i];
		for (uint32_t i = 0; i < tupLength[8]; i++) tupTemp[14][tupTempLength[14]++] = tup[8][i];

		for (uint32_t i = 0; i < tupLength[3]; i++) tupTemp[15][tupTempLength[15]++] = tup[3][i];
		for (uint32_t i = 0; i < tupLength[9]; i++) tupTemp[15][tupTempLength[15]++] = tup[9][i];

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

void XDDT(uint32_t first, uint32_t second, uint32_t output[1<<SIZE])
{
	if (SIZE == 4)
	{
		for (int v1 = 0; v1 < (1<<SIZE); v1++)
		{
			for (int v2 = 0; v2 < (1<<SIZE); v2++)
			{
				if (((v1^v2) == first) && ((getSbox(v1)^getSbox(v2)) == second)) output[v1]++;
			}
		}
	}
	else if (SIZE == 8)
	{
		for (int v1 = 0; v1 < (1<<SIZE); v1++)
		{
			for (int v2 = 0; v2 < (1<<SIZE); v2++)
			{
				if (((v1^v2) == first) && ((getSbox8(v1)^getSbox8(v2)) == second)) output[v1]++;
			}
		}
	}
}

void YDDT(uint32_t first, uint32_t second, uint32_t output[1<<SIZE])
{
	if (SIZE == 4)
	{
		for (int v1 = 0; v1 < (1<<SIZE); v1++)
		{
			for (int v2 = 0; v2 < (1<<SIZE); v2++)
			{
				if (((v1^v2) == first) && ((getSbox(v1)^getSbox(v2)) == second)) output[getSbox(v1)]++;
			}
		}
	}
	else if (SIZE == 8)
	{
		for (int v1 = 0; v1 < (1<<SIZE); v1++)
		{
			for (int v2 = 0; v2 < (1<<SIZE); v2++)
			{
				if (((v1^v2) == first) && ((getSbox8(v1)^getSbox8(v2)) == second)) output[getSbox8(v1)]++;
			}
		}
	}
}

void XOR(uint32_t val1[(1<<SIZE)],uint32_t val2[(1<<SIZE)])
{
	int sum1 = 0, sum2 = 0;
	for (int i = 0; i < (1<<SIZE); i++)
	{
		sum1 += val1[i];
		sum2 += val2[i];
	}
	if (sum2 == 0) return;
	if (sum1 == 0)
	{
		for (int i = 0; i < (1<<SIZE); i++) val1[i] = val2[i];
		return;
	}
	uint32_t res[(1<<SIZE)] = {0};
	for (int i = 0; i < (1<<SIZE); i++)
	{
		for (int j = 0; j < (1<<SIZE); j++) 
		{
			res[i^j] += val1[i] * val2[j];
		}
	}
	for (int i = 0; i < (1<<SIZE); i++)
	{
		val1[i] = res[i];
	}
}

void reduceKey(uint32_t *val, uint32_t valLength)
{
	int counter = 0;
	int a[2]; a[0] = 1; a[1] = 1;
	for (int i = 0; i < (1<<(valLength*SIZE)); i++)
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
	for (int i = 0; i < (1<<(valLength*SIZE)); i++) val[i] /= a[0];
}

bool isLinear(vector<uint32_t> constraint)
{
	uint32_t ANDval = (1 << SIZE) - 1;
	for (uint32_t i = 0; i < constraint.size();i++) if ((constraint[i] & ANDval) == 0) return false;
	return true;
}

void addRoundConstant(uint32_t tup, uint32_t val[(1<<SIZE)])
{
	uint32_t constantArray[(1<<SIZE)] = {0};
	if (((tup >> (SIZE*2)) & 0xf) == 0) constantArray[getConstants(tup >> (4+SIZE*2)) & 0xf] = 1;
	else if (((tup >> (SIZE*2)) & 0xf) == 4) constantArray[(getConstants(tup >> (4+SIZE*2)) >> 4) & 0b11] = 1;
	else if (((tup >> (SIZE*2)) & 0xf) == 8) constantArray[2] = 1;
	else constantArray[0] = 1;
	XOR(val,constantArray);
	reduceKey(val,1);
}

vector<int> resolveLinear(vector<uint32_t> constraint, uint32_t* &rightPartValues)
{
	vector<int> keys;
	// the last int in the constraint is the one with the right part (use xddt)
	int t = constraint[constraint.size()-1];
	uint32_t ANDval = (1 << SIZE) - 1;
	uint32_t xddt[1<<SIZE] = {0};
	XDDT((t >> SIZE) & ANDval,t & ANDval,xddt);
	for (int i = 0; i < (1<<SIZE); i++) {rightPartValues[i] = xddt[i];}
	for (uint32_t i = 0; i < constraint.size()-1; i++)
	{
		if (((constraint[i] >> (2*SIZE)) & 0xf) < 8) keys.push_back(constraint[i]);
		addRoundConstant(constraint[i],rightPartValues);
		uint32_t yddt[1<<SIZE] = {0};
		YDDT((constraint[i] >> SIZE) & ANDval,constraint[i] & ANDval,yddt);
		XOR(rightPartValues,yddt);
		reduceKey(rightPartValues,1);
	}
	return keys;
}



// here comes the nonlinear parts!

vector<int> resolveLeft(uint32_t zeroPoint, vector<uint32_t> constraint, uint32_t values[(1<<SIZE)], vector<uint32_t> &zeroDependencies)
{
	uint64_t ROUNDval = (1 << 4) - 1;
	uint32_t ANDval = (1 << SIZE) - 1;
	vector<int> keys;
	int inverseMC[16][3] = {{0,10,13},{1,11,14},{2,8,15},{3,9,12},{0,-1,-1},{1,-1,-1},{2,-1,-1},{3,-1,-1},
							{7,10,-1},{4,11,-1},{5,8,-1},{6,9,-1},{0,10,-1},{1,11,-1},{2,8,-1},{3,9,-1}};
	uint32_t round = (zeroPoint >> (4+2*SIZE));
	uint32_t pos = (zeroPoint >> (2*SIZE)) & 0xf;
	int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
	values[0] = 1;
	for (uint32_t i = 0; i < constraint.size(); i++)
	{
		if (!((((constraint[i] >> (4+2*SIZE)) == (round-1)) && (((constraint[i] >> (2*SIZE)) & ROUNDval) == inversePos[0])) ||
			 (((constraint[i] >> (4+2*SIZE)) == (round-1)) && (((constraint[i] >> (2*SIZE)) & ROUNDval) == inversePos[1])) ||
			 (((constraint[i] >> (4+2*SIZE)) == (round-1)) && (((constraint[i] >> (2*SIZE)) & ROUNDval) == inversePos[2])))) continue;
		if (((constraint[i] >> (2*SIZE)) & ANDval) < 8) keys.push_back(constraint[i]);
		if ((constraint[i] & ANDval) == 0) zeroDependencies.push_back(constraint[i]); 
		else
		{
			uint32_t yddt[1<<SIZE] = {0};
			YDDT((constraint[i] >> SIZE) & ANDval,constraint[i] & ANDval,yddt);
			addRoundConstant(constraint[i],values);
			XOR(values,yddt);
			reduceKey(values,1);
		}
	}
	return keys;
}

vector<int> resolveLast(uint32_t last, vector<uint32_t> constraint, uint32_t values[1<<SIZE], vector<uint32_t> &zeroDependencies)
{
	// this is almost similar to resolveLeft, with the additional xddt required
	vector<int> keys;
	uint32_t ROUNDval = (1 << 4) - 1;
	uint32_t ANDval = (1 << SIZE) - 1;
	int inverseMC[16][3] = {{0,10,13},{1,11,14},{2,8,15},{3,9,12},{0,-1,-1},{1,-1,-1},{2,-1,-1},{3,-1,-1},
							{7,10,-1},{4,11,-1},{5,8,-1},{6,9,-1},{0,10,-1},{1,11,-1},{2,8,-1},{3,9,-1}};
	uint32_t round = (last >> (4+2*SIZE));
	uint32_t pos = (last >> (2*SIZE)) & 0xf;
	int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
	// put the xddt into it
	values[0] = 1;
	uint32_t xddt[1<<SIZE] = {0};
	XDDT((last >> SIZE) & ANDval,last & ANDval,xddt);
	XOR(values,xddt); // this is the only difference from resolveLeft (with xddt)
	for (uint32_t i = 0; i < constraint.size(); i++)
	{
		if (!((((constraint[i] >> (4+2*SIZE)) == (round-1)) && (((constraint[i] >> (2*SIZE)) & ROUNDval) == inversePos[0])) ||
			 (((constraint[i] >> (4+2*SIZE)) == (round-1)) && (((constraint[i] >> (2*SIZE)) & ROUNDval) == inversePos[1])) ||
			 (((constraint[i] >> (4+2*SIZE)) == (round-1)) && (((constraint[i] >> (2*SIZE)) & ROUNDval) == inversePos[2])))) continue;
		if (((constraint[i] >> (2*SIZE)) & ANDval) < 8) keys.push_back(constraint[i]);
		if ((constraint[i] & ANDval) == 0) zeroDependencies.push_back(constraint[i]); 
		else
		{
			uint32_t yddt[1<<SIZE] = {0};
			YDDT((constraint[i] >> SIZE) & ANDval,constraint[i] & ANDval,yddt);
			addRoundConstant(constraint[i],values);
			XOR(values,yddt);
			reduceKey(values,1);
		}
	}
	return keys;
}

void AddKeyAndSubst(uint32_t* &value,uint32_t &valueLength)
{
	uint32_t valueSize = (1 << ((valueLength-1) * SIZE));
	uint32_t destinationSize = (1 << ((valueLength+1)*SIZE));
	uint32_t* destination;
	destination = new uint32_t [destinationSize]();
	if (SIZE == 4)
	{
		for (uint32_t i = 0; i < valueSize; i++)
		{
			for (int k = 0; k < 16; k++)
			{
				for (int v = 0; v < 16; v++)
				{
					destination[i*(16*16)+k*16+getSbox(v^k)] += value[i * 16 + v];
				}
			}
		}
	}
	else if (SIZE == 8)
	{
		for (uint32_t i = 0; i < valueSize; i++)
		{
			for (int k = 0; k < 256; k++)
			{
				for (int v = 0; v < 256; v++)
				{
					destination[i*(256*256)+k*256+getSbox8(v^k)] += value[i * 256 + v];
				}
			}
		}
	}

	delete[] value;
	value = destination;
	valueLength++;
}


void resolveZeroDependencies(uint32_t* valuesToAdd, uint32_t valuesToAddLength, uint32_t* &destination, uint32_t &destinationLength)
{
	uint32_t valuesKeySize = 1 << ((valuesToAddLength-1) * SIZE);
	uint32_t destinationKeySize = 1 << ((destinationLength-1) * SIZE);
	uint32_t *dest;
 	dest = new uint32_t [(1ULL<<((valuesToAddLength+destinationLength-1)*SIZE))]();
	for (uint32_t i = 0; i < valuesKeySize; i++)
	{
		for (uint32_t j = 0; j < destinationKeySize; j++)
		{	
			XOR(&dest[i*(destinationKeySize*(1 << SIZE))+j*(1 << SIZE)],&valuesToAdd[i*(1 << SIZE)]);
			XOR(&dest[i*(destinationKeySize*(1 << SIZE))+j*(1 << SIZE)],&destination[j*(1 << SIZE)]);
		}
	}
	reduceKey(dest,valuesToAddLength + destinationLength - 1);
	delete[] destination;
	destination = dest;
	destinationLength = valuesToAddLength + destinationLength - 1;
}

vector<int> resolveNonLinear(vector<uint32_t> constraint, uint32_t* &values, uint32_t &valuesLength)
{
	// first, we need to locate the zeros (which determine how many rounds/ connection points)
	sort(constraint.begin(),constraint.end());
	vector<uint32_t> zeroPoints; 
	
	uint32_t ANDval = (1 << SIZE) - 1;
	for (uint32_t i = 0; i < constraint.size(); i++) 
	{
		if ((constraint[i] & ANDval) == 0) zeroPoints.push_back(constraint[i]);
	}
	uint32_t* zeroValues[zeroPoints.size()] = {0};
	uint32_t zeroValuesLength[zeroPoints.size()] = {0}; // this will tell me how many keys are already involved
	vector<vector<int>> zeroKeys;
	for (uint32_t i = 0; i < zeroPoints.size(); i++)
	{
		vector<uint32_t> zeroDependencies;
		uint32_t* val;
		val = new uint32_t [(1<<SIZE)]();
		uint32_t valLength = 1; // by default, the value is 1. Unless there are more than 1 zeroDependencies adding to the same zeroPoint
		vector<int> keys = resolveLeft(zeroPoints[i],constraint, val, zeroDependencies); // zeroDependencies are the ones that require us to add the dimension of key
		// settle if there is any zeroDependencies
		for (uint32_t j = 0; j < zeroDependencies.size(); j++)
		{
			uint32_t index;
			// find what index of zeroPoint is the zeroDependices[j]
			for (uint32_t k = 0; k < zeroPoints.size(); k++) { if (zeroPoints[k] == zeroDependencies[j]) index = k;}
			resolveZeroDependencies(zeroValues[index],zeroValuesLength[index],val,valLength);
			keys.insert(keys.begin(),int(-1)); // spliting 
			keys.insert(keys.begin(),zeroKeys[index].begin(),zeroKeys[index].end());
		}
		AddKeyAndSubst(val,valLength);
		zeroKeys.push_back(keys);
		zeroValues[i] = val;
		zeroValuesLength[i] = valLength;
	}
	// Settle the last tuple!
	uint32_t* val;
	val = new uint32_t [1<<SIZE]();
	uint32_t valLength = 1;
	vector<uint32_t> zeroDependencies;
	vector<int> keys = resolveLast(constraint[constraint.size()-1],constraint, val, zeroDependencies);
	if (zeroDependencies.size() == 0) // sanity check
	{
		cout << "zeroDependencies cannot be zero at this point" << endl;
		exit(0);
	}
	for (uint32_t j = 0; j < zeroDependencies.size(); j++)
	{
		uint32_t index;
		// find what index of zeroPoint is the zeroDependices[j]
		for (uint32_t k = 0; k < zeroPoints.size(); k++) { if (zeroPoints[k] == zeroDependencies[j]) index = k;}
		resolveZeroDependencies(zeroValues[index],zeroValuesLength[index],val,valLength);
		keys.insert(keys.begin(),int(-1)); // spliting 
		keys.insert(keys.begin(),zeroKeys[index].begin(),zeroKeys[index].end());
	}
	values = val;
	valuesLength = valLength;
	return keys;
}

// Propagating to the last round
void propagateLinear(vector<int> &keys, uint32_t* &values, uint32_t nr, uint32_t &valuesLength, int TK_NUM)
{
	// round number should be the same for all the keys here
	uint32_t n = keys[0] >> (4+2*SIZE);
	uint32_t callsToLFSR = ((nr - n)/2) * keys.size();
	uint32_t TK2_columns[1<<(2*SIZE)];
	uint32_t TK3_columns[1<<(2*SIZE)];
	if (SIZE == 4)
	{
		for (int col = 0; col < 16; col++)
		{
			int tmp2 = col;
			int tmp3 = col;
			for (uint32_t i = 0; i < callsToLFSR; i++)
			{
				tmp2 = LFSR(tmp2,2); tmp3 = LFSR(tmp3,3);
			}
			TK2_columns[col] = tmp2; // index --> TK2_column[index]
			TK3_columns[col] = tmp3; // index --> TK2_column[index]
		}
	}
	else if (SIZE == 8)
	{
		for (int col = 0; col < 256; col++)
		{
			int tmp2 = col;
			int tmp3 = col;
			for (uint32_t i = 0; i < callsToLFSR; i++)
			{
				tmp2 = LFSR8(tmp2,2); tmp3 = LFSR8(tmp3,3);
			}
			TK2_columns[col] = tmp2; // index --> TK2_column[index]
			TK3_columns[col] = tmp3; // index --> TK2_column[index]
		}
	}
	uint32_t ANDval2 = (1 << (2*SIZE)) - 1;
	// permutate the keys first
	for (uint32_t i = 0; i < keys.size(); i++)
	{
		while (uint32_t(keys[i] >> (4+2*SIZE)) < nr) 
		{
			for (int j = 0; j < 16; j++)
			{
				if (getPermSchedule(j) == ((keys[i] >> (2*SIZE)) & 0xf))
				{
					keys[i] = (((keys[i] >> (4+2*SIZE)) + 1) << (4+2*SIZE)) ^ (j << (2*SIZE)) ^ (keys[i] & ANDval2);
					break;
				}
			}
		}
	}
	uint32_t* newKey;
	valuesLength = TK_NUM;
	newKey = new uint32_t[1<<(SIZE*TK_NUM)]();
	if (TK_NUM == 1) return;
	else if (TK_NUM == 2)
	{
		valuesLength = 2;
		for (int tk1 = 0; tk1 < (1<<SIZE); tk1++)
		{
			for (int diff = 0; diff < (1<<SIZE); diff++)
			{
				if (values[diff] == 0) continue;
				int tk2 = tk1 ^ diff;
				newKey[tk1*(1<<SIZE)+TK2_columns[tk2]] += values[diff];
			}
		}
	}
	else if (TK_NUM == 3)
	{
		valuesLength = 3;
		for (int tk1 = 0; tk1 < (1<<SIZE); tk1++)
		{
			for (int tk2 = 0; tk2 < (1<<SIZE); tk2++)
			{
				for (int diff = 0; diff < (1<<SIZE); diff++)
				{
					if (values[diff] == 0) continue;
					int tk3 = tk1 ^ tk2 ^ diff;
					newKey[tk1*(1<<(2*SIZE))+TK2_columns[tk2]*(1<<SIZE)+TK3_columns[tk3]] += values[diff];
				}
			}
		}
	}
	values = newKey;
	return;
}


// need to revisit nonlinear once again
bool propagateNonLinear(vector<int> &keys, uint32_t* &values, uint32_t nr, uint32_t &valuesLength, int TK_NUM)
{
	if ((SIZE*valuesLength*TK_NUM) > 39)
	{
		cout << "Expected size of more than 2^39. Skipping..." << endl;
		return false;
	}
	
	uint32_t keyCounter = 0;
	int *TK2_columns; 
	int *TK3_columns; 
	int *callsToLFSR;
	TK2_columns = new int [valuesLength*(1<<SIZE)](); // to record the column transformation due to LFSR
	TK3_columns = new int [valuesLength*(1<<SIZE)](); // to record the column transformation due to LFSR
	callsToLFSR = new int [valuesLength](); // to record the number of LFSR required
	int loopCounter = 0;

	while (keyCounter < keys.size())
	{
		uint32_t n = keys[keyCounter] >> (4+2*SIZE);
		uint32_t noKeys = 0;
		while (keys[keyCounter] != -1) 
		{
			noKeys++;
			keyCounter++;
			if (keyCounter >= keys.size()) break;
		}
		callsToLFSR[loopCounter] = ((nr - n)/2) * noKeys;

		if (SIZE == 4)
		{
			for (int col = 0; col < 16; col++)
			{
				int tmp2 = col;
				int tmp3 = col;
				for (int i = 0; i < callsToLFSR[loopCounter]; i++)
				{
					tmp2 = LFSR(tmp2,2); tmp3 = LFSR(tmp3,3);
				}
				TK2_columns[loopCounter*16+col] = tmp2; // index --> TK2_column[index]
				TK3_columns[loopCounter*16+col] = tmp3; // index --> TK2_column[index]
			}
		}
		else if (SIZE == 8)
		{
			for (int col = 0; col < 256; col++)
			{
				int tmp2 = col;
				int tmp3 = col;
				for (int i = 0; i < callsToLFSR[loopCounter]; i++)
				{
					tmp2 = LFSR8(tmp2,2); tmp3 = LFSR8(tmp3,3);
				}
				TK2_columns[loopCounter*256+col] = tmp2; // index --> TK2_column[index]
				TK3_columns[loopCounter*256+col] = tmp3; // index --> TK2_column[index]
			}
		}
		loopCounter++;
		keyCounter++;
	}
	uint32_t ANDval2 = (1 << (2*SIZE)) - 1;
	// permutate the keys now
	for (uint32_t i = 0; i < keys.size(); i++)
	{
		if (keys[i] == -1) continue;
		while (uint32_t(keys[i] >> (4+2*SIZE)) < nr) 
		{
			for (int j = 0; j < (1<<SIZE); j++)
			{
				if (getPermSchedule(j) == ((keys[i] >> (2*SIZE)) & 0xf))
				{
					keys[i] = (((keys[i] >> (4+2*SIZE)) + 1) << (4+2*SIZE)) ^ (j << (2*SIZE)) ^ (keys[i] & ANDval2);
					break;
				}
			}
		}
	}
	uint32_t* newKey;
	cout << "Allocating an array of size: " << dec << (SIZE*valuesLength*TK_NUM) << endl;
	newKey = new uint32_t [1ULL<<(SIZE*valuesLength*TK_NUM)](); // size of the resulting array
	cout << "Successfully constructed the array" << endl;
	uint32_t ANDval = (1 << SIZE) - 1;

	if (TK_NUM == 1) return true;
	else if (TK_NUM == 2)
	{
		for (uint64_t tk1 = 0; tk1 < (1<<(SIZE*valuesLength)); tk1++)
		{
			for (uint64_t diff = 0; diff < (1<<(SIZE*valuesLength)); diff++)
			{
				uint64_t index = ((tk1 >> (SIZE*(valuesLength-1))) & ANDval) * (1<<SIZE) + TK2_columns[0*(1<<SIZE)+(((tk1 >> (SIZE*(valuesLength-1))) & ANDval) ^ ((diff >> (SIZE*(valuesLength-1))) & ANDval))];
				for (uint32_t i = 1; i < valuesLength; i++)
				{
					index = index << (TK_NUM*SIZE);
					index += ((tk1 >> (SIZE*(valuesLength-i-1))) & ANDval) * (1<<SIZE) + TK2_columns[i*(1<<SIZE)+(((tk1 >> (SIZE*(valuesLength-i-1))) & ANDval) ^ ((diff >> (SIZE*(valuesLength-i-1))) & (1<<SIZE)))];
				}
				newKey[index] += values[diff];
			}
		}
	}
	else if (TK_NUM == 3)
	{
		for (uint64_t tk1 = 0; tk1 < (1<<(SIZE*valuesLength)); tk1++)
		{
			for (uint64_t tk2 = 0; tk2 < (1<<(SIZE*valuesLength)); tk2++)
			{
				for (uint64_t diff = 0; diff < (1<<(SIZE*valuesLength)); diff++)
				{
					uint64_t index = ((tk1 >> (SIZE*(valuesLength-1))) & ANDval) * (1<<SIZE) * (1<<SIZE) + TK2_columns[0*(1<<SIZE)+((tk2 >> (SIZE*(valuesLength-1))) & ANDval)] * (1<<SIZE) + 
					TK3_columns[0*(1<<SIZE)+((((tk1^tk2) >> (SIZE*(valuesLength-1))) & ANDval) ^ ((diff >> (SIZE*(valuesLength-1))) & ANDval))];
					for (uint32_t i = 1; i < valuesLength; i++)
					{
						index = index << (TK_NUM*SIZE);
						index += ((tk1 >> (SIZE*(valuesLength-i-1))) & ANDval) * (1<<SIZE) * (1<<SIZE) + TK2_columns[i*(1<<SIZE)+((tk2 >> (SIZE*(valuesLength-i-1))) & ANDval)] * (1<<SIZE) +
						TK3_columns[i*(1<<SIZE)+((((tk1^tk2) >> ((1<<SIZE)*(valuesLength-i-1))) & ANDval) ^ ((diff >> (SIZE*(valuesLength-i-1))) & ANDval))];
					}
					newKey[index] += values[diff];
				}
			}
		}
	}
	valuesLength = valuesLength * TK_NUM;


	values = newKey;
	return true;
}


bool isZero(uint32_t *array, uint32_t length)
{
	for (uint32_t i = 0; i < length; i++)
	{
		if (array[i] != 0) return false;
	}
	return true;
}



void compareConstraints(uint32_t **linearValues, uint32_t *linearValuesLength, vector<vector<int>> linearKeys, uint32_t **nonLinearValues, uint32_t *nonLinearValuesLength, vector<vector<int>> nonLinearKeys, int TK_NUM)
{
	cout << "Comparing Constraints..." << endl;
	// linearValues is split into (TK1,TK2)
	// comparing linear with linear
	uint32_t *valuesByCells; valuesByCells = new uint32_t[(1ULL << (SIZE*TK_NUM))*16]; // 16 cells, with 2**(4 * TK_NUM) keys for each cell 
	for (uint64_t i = 0; i < (16*(1 << (SIZE*TK_NUM))); i++) valuesByCells[i] = 1;
	
	// linear
	for (uint32_t i = 0; i < linearKeys.size(); i++)
	{
		for (int j = 0; j < (1 << (SIZE*TK_NUM)); j++) valuesByCells[((linearKeys[i][0] >> (2*SIZE)) & 0xf) * (1 << (SIZE*TK_NUM)) + j] *= linearValues[i][j];
	}


	for (int cell = 0; cell < 16; cell++)
	{
		if (isZero(&valuesByCells[cell*(1<<(SIZE*TK_NUM))],(1<<(SIZE*TK_NUM))))
		{
			cout << "impossible due to cell " << cell << endl;
		}
	}

	cout << "Postprocessing..." << endl;
	// post processing for nonlinear constraints to reduce the number of clauses

	for (uint32_t i = 0; i < nonLinearKeys.size(); i++)
	{
		vector<int> nonLinearKeysTmp;
		for (uint32_t j = 0; j < nonLinearKeys[i].size(); j++)
		{
			if (nonLinearKeys[i][j] != -1) nonLinearKeysTmp.push_back(nonLinearKeys[i][j]);
		}
	
		int dim = nonLinearValuesLength[i] / TK_NUM;
		for (int j = 0; j < dim; j++)
		{
			uint32_t cell = ((nonLinearKeys[i][j] >> (2*SIZE)) & 0xf);
			// check if any of the values in this cell that is not accepted from the linear constraints. We remove them
			for (uint32_t k = 0; k < (1 << (TK_NUM*SIZE)); k++)
			{
				if (valuesByCells[(1ULL<< (SIZE*TK_NUM))*cell + k] == 0) // Eliminate all possible values from this
				{
					for (uint64_t l = 0; l < (1ULL<<((dim-1)*TK_NUM*SIZE)); l++)
					{
						uint64_t upperVal = (l >> ((dim-1-j)*TK_NUM*SIZE)) << ((dim-j)*TK_NUM*SIZE);
						uint64_t lowerVal = l & ((1ULL << ((dim-1-j)*TK_NUM*SIZE)) - 1);
						nonLinearValues[i][upperVal + k*(1ULL << ((dim-1-j)*TK_NUM*SIZE)) + lowerVal] = 0;
					}
				}
			}
		}
	}




	delete valuesByCells;

	cout << "Using the SAT solver..." << endl;
	// using the SAT solver
	SATsolver(linearKeys,linearValues,linearValuesLength,nonLinearKeys,nonLinearValues,nonLinearValuesLength,TK_NUM);
}

bool propagateNonLinearBool(vector<int> &keys, uint32_t* &values, uint32_t nr, uint32_t &valuesLength, int TK_NUM, uint16_t * &boolVals)
{
	// only keeping the bool information
	if ((SIZE*valuesLength*TK_NUM) > 39)
	{
		cout << "Expected size of more than 2^39. Skipping..." << endl;
		return false;
	}
	
	uint32_t keyCounter = 0;
	int *TK2_columns; 
	int *TK3_columns; 
	int *callsToLFSR;
	TK2_columns = new int [valuesLength*(1<<SIZE)](); // to record the column transformation due to LFSR
	TK3_columns = new int [valuesLength*(1<<SIZE)](); // to record the column transformation due to LFSR
	callsToLFSR = new int [valuesLength](); // to record the number of LFSR required
	int loopCounter = 0;

	while (keyCounter < keys.size())
	{
		int n = keys[keyCounter] >> (4+2*SIZE);
		int noKeys = 0;
		while (keys[keyCounter] != -1) 
		{
			noKeys++;
			keyCounter++;
			if (keyCounter >= keys.size()) break;
		}
		callsToLFSR[loopCounter] = ((nr - n)/2) * noKeys;

		if (SIZE == 4)
		{
			for (int col = 0; col < 16; col++)
			{
				int tmp2 = col;
				int tmp3 = col;
				for (int i = 0; i < callsToLFSR[loopCounter]; i++)
				{
					tmp2 = LFSR(tmp2,2); tmp3 = LFSR(tmp3,3);
				}
				TK2_columns[loopCounter*16+col] = tmp2; // index --> TK2_column[index]
				TK3_columns[loopCounter*16+col] = tmp3; // index --> TK2_column[index]
			}
		}
		else if (SIZE == 8)
		{
			for (int col = 0; col < 256; col++)
			{
				int tmp2 = col;
				int tmp3 = col;
				for (int i = 0; i < callsToLFSR[loopCounter]; i++)
				{
					tmp2 = LFSR8(tmp2,2); tmp3 = LFSR8(tmp3,3);
				}
				TK2_columns[loopCounter*256+col] = tmp2; // index --> TK2_column[index]
				TK3_columns[loopCounter*256+col] = tmp3; // index --> TK2_column[index]
			}
		}
		loopCounter++;
		keyCounter++;
	}
	// permutate the keys now
	uint32_t ANDval2 = (1 << (2*SIZE)) - 1;
	for (uint32_t i = 0; i < keys.size(); i++)
	{
		if (keys[i] == -1) continue;
		while ((keys[i] >> (4+2*SIZE)) < nr) 
		{
			for (int j = 0; j < 16; j++)
			{
				if (getPermSchedule(j) == ((keys[i] >> (2*SIZE)) & 0xf))
				{
					keys[i] = (((keys[i] >> (4+2*SIZE)) + 1) << (4+2*SIZE)) ^ (j << (2*SIZE)) ^ (keys[i] & ANDval2);
					break;
				}
			}
		}
	}
	uint16_t* newKey;
	cout << "Allocating an array of size: " << dec << (SIZE*valuesLength*TK_NUM-4) << endl;
	newKey = new uint16_t [1ULL<<(SIZE*valuesLength*TK_NUM-4)](); // size of the resulting array
	cout << "Successfully constructed the array" << endl;
	uint32_t ANDval = (1 << SIZE) - 1;

	if (TK_NUM == 1)
	{
		for (int i = 0; i < (1<<(SIZE)); i++) newKey[i] = values[i];
	}
	else if (TK_NUM == 2)
	{
		for (int tk1 = 0; tk1 < (1<<(SIZE*valuesLength)); tk1++)
		{
			for (int diff = 0; diff < (1<<(SIZE*valuesLength)); diff++)
			{
				uint64_t index = ((tk1 >> (SIZE*(valuesLength-1))) & ANDval) * (1<<SIZE) + TK2_columns[0*(1<<SIZE)+(((tk1 >> (SIZE*(valuesLength-1))) & ANDval) ^ ((diff >> (SIZE*(valuesLength-1))) & ANDval))];
				for (int i = 1; i < valuesLength; i++)
				{
					index = index << (TK_NUM*SIZE);
					index += ((tk1 >> (SIZE*(valuesLength-i-1))) & ANDval) * (1<<SIZE) + TK2_columns[i*(1<<SIZE)+(((tk1 >> (SIZE*(valuesLength-i-1))) & ANDval) ^ ((diff >> (SIZE*(valuesLength-i-1))) & ANDval))];
				}
				// cout << index << endl;
				if (values[diff] > 0) newKey[index/16] = newKey[index/16] | (1ULL << (15-(index%16)));
				else newKey[index/16] = (newKey[index/16] & ((1ULL<<(index%16))-1) << (16-(index%16))) + (newKey[index/16] & ((1ULL << (15-(index%15)))-1));
			}
		}
	}
	else if (TK_NUM == 3)
	{
		for (int tk1 = 0; tk1 < (1<<(SIZE*valuesLength)); tk1++)
		{
			for (int tk2 = 0; tk2 < (1<<(SIZE*valuesLength)); tk2++)
			{
				for (int diff = 0; diff < (1<<(SIZE*valuesLength)); diff++)
				{
					uint64_t index = ((tk1 >> (SIZE*(valuesLength-1))) & ANDval) * (1<<SIZE) * (1<<SIZE) + TK2_columns[0*(1<<SIZE)+((tk2 >> (SIZE*(valuesLength-1))) & ANDval)] * (1<<SIZE) + 
					TK3_columns[0*(1<<SIZE)+((((tk1^tk2) >> (SIZE*(valuesLength-1))) & ANDval) ^ ((diff >> (SIZE*(valuesLength-1))) & ANDval))];
					for (int i = 1; i < valuesLength; i++)
					{
						index = index << (TK_NUM*SIZE);
						index += ((tk1 >> (SIZE*(valuesLength-i-1))) & ANDval) * (1<<SIZE) * (1<<SIZE) + TK2_columns[i*(1<<SIZE)+((tk2 >> (SIZE*(valuesLength-i-1))) & ANDval)] * (1<<SIZE) +
						TK3_columns[i*(1<<SIZE)+((((tk1^tk2) >> (SIZE*(valuesLength-i-1))) & ANDval) ^ ((diff >> (SIZE*(valuesLength-i-1))) & ANDval))];
					}
					if (values[diff] > 0)
					{
						newKey[index/16] = newKey[index/16] | (1ULL << (15-(index%16)));
					}
					else
					{
						newKey[index/16] = (newKey[index/16] & ((1ULL<<(index%16))-1) << (16-(index%16))) + (newKey[index/16] & ((1ULL << (15-(index%16)))-1));
					}
				}
			}
		}
	}
	valuesLength = valuesLength * TK_NUM;
	boolVals = newKey;
	return true;
}



// what if there are more than 1 key?!
bool compareLinearConstraintsOnly(uint32_t **linearValues, uint32_t *linearValuesLength, vector<vector<int>> linearKeys, int TK_NUM)
{
	cout << "Comparing Constraints..." << endl;
	// linearValues is split into (TK1,TK2)
	// comparing linear with linear
	uint32_t *valuesByCells; valuesByCells = new uint32_t[(1ULL << (SIZE*TK_NUM))*16]; // 16 cells, with 2**(4 * TK_NUM) keys for each cell 
	for (int i = 0; i < (16*(1 << (SIZE*TK_NUM))); i++) valuesByCells[i] = 1;
	
	// linear
	for (uint32_t i = 0; i < linearKeys.size(); i++)
	{
		for (int j = 0; j < (1 << (SIZE*TK_NUM)); j++) valuesByCells[((linearKeys[i][0] >> (2*SIZE)) & 0xf) * (1 << (SIZE*TK_NUM)) + j] *= linearValues[i][j];
	}
	
	bool flag = true;
	for (int cell = 0; cell < 16; cell++)
	{
		if (isZero(&valuesByCells[cell*(1<<(SIZE*TK_NUM))],(1<<(SIZE*TK_NUM))))
		{
			cout << "impossible due to cell " << cell << endl;
			flag = false;
		}
	}
	return flag;
}

void compareConstraintsBool(uint32_t **linearValues, uint32_t *linearValuesLength, vector<vector<int>> linearKeys, uint16_t **nonLinearValuesBool, uint32_t *nonLinearValuesLength, vector<vector<int>> nonLinearKeys, int TK_NUM)
{
	cout << "Comparing Constraints..." << endl;
	// linearValues is split into (TK1,TK2)
	// comparing linear with linear
	uint32_t *valuesByCells; valuesByCells = new uint32_t[(1ULL << (SIZE*TK_NUM))*16]; // 16 cells, with 2**(4 * TK_NUM) keys for each cell 
	for (int i = 0; i < (16*(1 << (SIZE*TK_NUM))); i++) valuesByCells[i] = 1;
	
	// linear
	for (uint32_t i = 0; i < linearKeys.size(); i++)
	{
		for (int j = 0; j < (1 << (SIZE*TK_NUM)); j++) valuesByCells[((linearKeys[i][0] >> (2*SIZE)) & 0xf) * (1 << (SIZE*TK_NUM)) + j] *= linearValues[i][j];
	}
	

	for (int cell = 0; cell < 16; cell++)
	{
		if (isZero(&valuesByCells[cell*(1<<(SIZE*TK_NUM))],(1<<(SIZE*TK_NUM))))
		{
			cout << "impossible due to cell " << cell << endl;
		}
	}

	cout << "Postprocessing..." << endl;
	// post processing for nonlinear constraints to reduce the number of clauses

	for (uint32_t i = 0; i < nonLinearKeys.size(); i++)
	{
		vector<int> nonLinearKeysTmp;
		for (uint32_t j = 0; j < nonLinearKeys[i].size(); j++)
		{
			if (nonLinearKeys[i][j] != -1) nonLinearKeysTmp.push_back(nonLinearKeys[i][j]);
		}
	
		int dim = nonLinearValuesLength[i] / TK_NUM;
		for (int j = 0; j < dim; j++)
		{
			uint32_t cell = ((nonLinearKeys[i][j] >> (2*SIZE)) & 0xf);
			// check if any of the values in this cell that is not accepted from the linear constraints. We remove them
			for (int k = 0; k < (1 << (TK_NUM*SIZE)); k++)
			{
				if (valuesByCells[(1ULL<< (SIZE*TK_NUM))*cell + k] == 0) // Eliminate all possible values from this
				{
					for (uint64_t l = 0; l < (1ULL<<((dim-1)*TK_NUM*SIZE)); l++)
					{
						uint64_t upperVal = (l >> ((dim-1-j)*TK_NUM*SIZE)) << ((dim-j)*TK_NUM*SIZE);
						uint64_t lowerVal = l & ((1ULL << ((dim-1-j)*TK_NUM*SIZE)) - 1);
						uint64_t index = upperVal + k*(1ULL << ((dim-1-j)*TK_NUM*SIZE)) + lowerVal;
						uint64_t remainder = index % 16;
						nonLinearValuesBool[i][index/16] = (nonLinearValuesBool[i][index/16] & (((1ULL << remainder)-1) << (16-remainder)))
															+ (nonLinearValuesBool[i][index/16] & ((1ULL << (15-remainder)) - 1));
					}
				}
			}
		}
	}

	delete valuesByCells;

	cout << "Using the SAT solver..." << endl;
	// using the SAT solver
	SATsolverBool(linearKeys,linearValues,linearValuesLength,nonLinearKeys,nonLinearValuesBool,nonLinearValuesLength,TK_NUM);
}


int computeProbTrail(uint8_t alpha[20][4][4], uint8_t key_diff[4][4][4], uint32_t nr, int TK_NUM, vector<vector<uint32_t>> linearConstraints, uint32_t **linearValues, uint32_t *linearValuesLength, vector<vector<int>> linearKeys)
{
	int prob = 0;
	int DDT4[16][16] = {0};
	computeDDT(DDT4);
	int DDT8[256][256] = {0};
	computeDDT8(DDT8);
	uint32_t ANDval = (1 << SIZE) - 1;
	uint8_t key_diff_temp[4][4][4];
	for (int i = 0; i < 4; i++){ for (int j = 0; j < 4; j++) { for (int k = 0; k < 4; k++) key_diff_temp[i][j][k] = key_diff[i][j][k];}}
	uint8_t before[4][4], after[4][4];
	for (int n = 0; n < nr; n++)
	{
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				before[r][c] = alpha[n][r][c];
				after[r][c] = alpha[n+1][r][c];
			}
		}
		invMC(after); invSR(after); AddRoundTweakey(after,key_diff_temp);
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				if (before[r][c] == 0) continue;
				if (SIZE == 4) prob += int(-log2(DDT4[before[r][c]][after[r][c]]/16.0));
				if (SIZE == 8) prob += int(-log2(DDT8[before[r][c]][after[r][c]]/256.0));
			}
		}

		if (SIZE == 4) key_schedule_round_function(key_diff_temp);
		else key_schedule_round_function8(key_diff_temp);
	}
	cout << "Original probability: " << dec << prob << endl;
	for (int i = 0; i < linearConstraints.size(); i++)
	{
		uint8_t inp = ((linearConstraints[i][linearConstraints[i].size()-1] >> SIZE) & ANDval);
		uint8_t outp = ((linearConstraints[i][linearConstraints[i].size()-1]) & ANDval);
		if (SIZE == 4) prob -= int(-log2(DDT4[inp][outp]/16.0));
		if (SIZE == 8) prob -= int(-log2(DDT8[inp][outp]/256.0));
	}
	cout << "Probability that does not depend on key: " << prob << endl;
	int noKeys = 0;
	uint32_t *valuesByCells; valuesByCells = new uint32_t[(1ULL << (SIZE*TK_NUM))*16]; // 16 cells, with 2**(4 * TK_NUM) keys for each cell 
	for (int i = 0; i < (16*(1 << (SIZE*TK_NUM))); i++) valuesByCells[i] = 1;
	
	// linear
	for (uint32_t i = 0; i < linearKeys.size(); i++)
	{
		for (int j = 0; j < (1 << (SIZE*TK_NUM)); j++) valuesByCells[((linearKeys[i][0] >> (2*SIZE)) & 0xf) * (1 << (SIZE*TK_NUM)) + j] *= linearValues[i][j];
	}

	for (int cell = 0; cell < 16; cell++)
	{
		int s = 0;
		for (int i = 0; i < (1<<(SIZE*TK_NUM)); i++)
		{
			// cout << int(valuesByCells[cell*(1<<(SIZE*TK_NUM))+i]) << " ";
			if (int(valuesByCells[cell*(1<<(SIZE*TK_NUM))+i]) > 0) s += 1;
		}
		// cout << endl;
		if (s == 0) 
		{
			return 0;
		}
		noKeys += int(-log2((s+0.0)/(1<<(SIZE*TK_NUM))));
	}
	delete valuesByCells;
	return noKeys;
}

int main()
{
	// for (int q = 1; q < 12; q++)
	for (int u_index = 0; u_index < 5; u_index++)
	{
		for (int v_index = 0; v_index < 7; v_index++)
	{
	// 	for (int w_index = 0; w_index < 18; w_index++)
	// {
	uint8_t alpha[20][4][4];
	uint8_t key_diff[4][4][4];
	uint32_t nr;
	int TK_NUM;
	// u = 0x08,0x0c,0x0d,0x0e,0x0f
	// v = 0x28, 0x29, 0x2c, 0x2d, 0x2e, 0x2f, 0xb9	
	int u_list[5] = {0x08,0x0c,0x0d,0x0e,0x0f};
	int v_list[7] = {0x28, 0x29, 0x2c, 0x2d, 0x2e, 0x2f, 0xb9};
	// int w_list[18] = {0x2a,0x2d,0x2e,0x2f,0x6b,0x7c,0xb8,0xb9,0xba,0xbd,0xbe,0xbf,0xeb,0xec,0xef,0xfb,0xfc,0xfe};
	int u = u_list[u_index];
	int v = v_list[v_index];
	// int w = w_list[w_index];
	cout << dec << "u:" << u << " v:" << v << endl;
	// int u = u_list[u_index];
	// int v = v_list[v_index];
	// cout << "q : " << q << endl;

	// 9/21 trails failed
 	// if (q == 1) SK_2020_1402(alpha,key_diff,nr,TK_NUM);
	// if (q == 2) TK1_2020_1402(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 3) TK2_2020_1402(alpha,key_diff,nr,TK_NUM); 
	// if (q == 4) TK3_2020_1402(alpha,key_diff,nr,TK_NUM);
	// if (q == 5) get_4TK2_1_17_L_2020_1317(alpha,key_diff,nr,TK_NUM); // this is used as a first test
	// if (q == 6) get_4TK2_1_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // without nonLinearConstraint[0], it's possible
	// if (q == 7) get_4TK2_1_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// if (q == 8) get_4TK2_1_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // without nonLinearConstraint[0], it's possible
	// if (q == 9) get_4TK2_2_17_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 10) get_4TK2_2_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// if (q == 11) get_4TK2_2_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 12) get_4TK2_2_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 13) get_4TK3_1_22_L_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// if (q == 14) get_4TK3_1_22_U_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 15) get_4TK3_1_22_L_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// if (q == 16) get_4TK3_1_22_U_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 17) get_4TK3_1_22_U_BMD3_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 18) get_4TK3_1_23_U_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 19) get_4TK3_1_23_U_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible linear
	
	
	// get_4TK2_18_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // u = 6, v = 2,10,12,13 // all ok
	// get_4TK3_22_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // u = 8,9, v = 5 // all ok
	// 8 bit trails
	// 8/13 trails failed
	// SIZE = 8;


	// if (q == 1) get_8TK2_1_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // possible
	// if (q == 2) get_8TK2_1_18_U_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 4) get_8TK2_1_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 5) get_8TK2_1_20_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 6) get_8TK2_1_20_U_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// if (q == 7) get_8TK2_1_21_L_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// if (q == 9) SK_tosc_2017_i4_99_129(alpha,key_diff,nr,TK_NUM); // imposible
	// if (q == 10) TK1_8_2020_1402(alpha,key_diff,nr,TK_NUM); // ok  
	// if (q == 11) TK2_8_2020_1402(alpha,key_diff,nr,TK_NUM); // ok

	// u = 0x40,0x50,0xc0,0xd0  v =  0x2a,0x2e,0x2f,0x6a,0xba,0xbe,0xbf,0xea,0xef,0xfa,0xfe
	// w = 0x2a,0x2d,0x2e,0x2f,0x6b,0x7c,0xb8,0xb9,0xba,0xbd,0xbe,0xbf,0xeb,0xec,0xef,0xfb,0xfc,0xfe
 	// get_8TK2_19_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v,w); // all ok

	// u = 0x08,0x0c,0x0d,0x0e,0x0f
	// v = 0x28, 0x29, 0x2c, 0x2d, 0x2e, 0x2f, 0xb9	
	get_8TK3_22_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // all ok
	

	// getTK1_8_5Rounds(alpha,key_diff,nr,TK_NUM); // this is used as a test case
	uint8_t diffStates[20][5][16] = {0}; // records the path 
	uint8_t fixedStates[20][5][16] = {0}; // '1' represents a state that is fixed in value
	fixedValuesToPath(alpha,key_diff,diffStates,fixedStates,nr);
	printState(diffStates,fixedStates,nr);
	vector<vector<uint32_t>> constraints = findConstraints(diffStates,fixedStates,nr);
	
	// just for printing
	uint32_t ANDval = (1 << SIZE) - 1;
	for (uint32_t i = 0; i < constraints.size(); i++)
	{
		for (uint32_t j = 0; j < constraints[i].size(); j++)
		{
			cout << hex << "(" << (constraints[i][j] >> (4+2*SIZE)) << "," << ((constraints[i][j] >> (2*SIZE)) & 0xf) << ","
			<< ((constraints[i][j] >> SIZE) & ANDval) << "," << (constraints[i][j] & ANDval) << ") ";
		}
		cout << isLinear(constraints[i]);
		cout << endl;
	}
	// end of just for printing

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


	for (uint32_t i = 0; i < linearConstraints.size(); i++)
	{
		uint32_t *values; values = new uint32_t[(1<<SIZE)]();
		uint32_t valuesLength = 1;
		cout << "Resolving linear constraints[" << i << "]...";
		vector<int> keys = resolveLinear(linearConstraints[i],values);
		cout << "Completed" << endl;
		cout << "Propagating linear constraints[" << i << "]...";
		propagateLinear(keys,values,nr,valuesLength,TK_NUM);
		cout << "Completed" << endl;
		linearValues[i] = values;
		linearValuesLength[i] = valuesLength;
		linearKeys.push_back(keys);
	}
	cout << "After propagating, the keys are:" << endl;
	for (int i = 0; i < linearKeys.size(); i++) cout << hex << linearKeys[i][0] << endl;
	// add the linear one here first
	if (compareLinearConstraintsOnly(linearValues,linearValuesLength,linearKeys,TK_NUM))
	{
		cout << "creating nonLinear arrays..."; 
		uint32_t **nonLinearValues; nonLinearValues = new uint32_t*[nonLinearConstraints.size()]();
		uint32_t *nonLinearValuesLength; nonLinearValuesLength = new uint32_t[nonLinearConstraints.size()]();
		uint16_t **nonLinearValuesBool; nonLinearValuesBool = new uint16_t*[nonLinearConstraints.size()]();
		cout << "Allocated " << nonLinearConstraints.size() << " for nonLinearValues" << endl;
		vector<vector<int>> nonLinearKeys;
		cout << "Completed" << endl;
		// nonlinear constraints (boolean)
		if (true)
		{
			for (uint32_t i = 0; i < nonLinearConstraints.size(); i++)
			{
				uint32_t *values;
				uint32_t valuesLength;
				uint16_t *boolVals;
				cout << "Resolving nonLinearConstraint[" << i << "]...";
				vector<int> keys = resolveNonLinear(nonLinearConstraints[i],values,valuesLength); // note that a -1 implies a split 
				cout << "Completed" << endl;
				cout << "Propagating nonLinearConstraint[" << i << "]...";
				bool t = propagateNonLinearBool(keys,values,nr,valuesLength,TK_NUM,boolVals);
				cout << "Completed" << endl;
				if (t)
				{
					nonLinearValuesBool[i] = boolVals;
					nonLinearValuesLength[i] = valuesLength;
					nonLinearKeys.push_back(keys);
				}
				else
				{
					nonLinearValuesLength[i] = 0;
				}
			}
			cout << "Successfully propagating all keys" << endl;
				// continue from here. Construct one that can accept boolean
			compareConstraintsBool(linearValues,linearValuesLength,linearKeys,nonLinearValuesBool,nonLinearValuesLength,nonLinearKeys,TK_NUM);
		}

		// nonlinear constraints (non-boolean)
		if (false)
		{
			for (uint32_t i = 0; i < nonLinearConstraints.size(); i++)
			{
				uint32_t *values;
				uint32_t valuesLength;
				cout << "resolving nonLinearConstraint[" << i << "]...";
				vector<int> keys = resolveNonLinear(nonLinearConstraints[i],values,valuesLength); // note that a -1 implies a split 
				cout << "Completed" << endl;

				cout << "propagating nonLinearConstraint[" << i << "]...";
				bool t = propagateNonLinear(keys,values,nr,valuesLength,TK_NUM);
				cout << "Completed" << endl;
				if (t)
				{
					nonLinearValues[i] = values;
					nonLinearValuesLength[i] = valuesLength;
					nonLinearKeys.push_back(keys);
				}
				else
				{
					nonLinearValuesLength[i] = 0;
				}
			}
			cout << "Successfully propagating all keys" << endl;
			compareConstraints(linearValues,linearValuesLength,linearKeys,nonLinearValues,nonLinearValuesLength,nonLinearKeys,TK_NUM);
		}
	}
	int noKeys = computeProbTrail(alpha,key_diff,nr,TK_NUM,linearConstraints,linearValues,linearValuesLength,linearKeys);
	cout << "no. of keys that allow the characteristic: " << noKeys << endl;
	}
}
	return 0;
}
