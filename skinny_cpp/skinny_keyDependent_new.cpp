#include "skinny.h"
#include "trails.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <bits/stdc++.h>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <unordered_set>
using namespace std;

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl

#ifndef SIZE
#define SIZE 4;
#endif
// general requirements
int TK_NUM;
uint32_t nr;
uint8_t alpha[20][4][4];
uint8_t key_diff[4][4][4];
uint8_t diffStates[20][5][16] = {0};
uint8_t fixedStates[20][5][16] = {0};
uint32_t ANDval = (1<<SIZE)-1;
uint64_t ROUNDval = (1 << 4) - 1; 
int inverseMC[16][3] = {{0,10,13},{1,11,14},{2,8,15},{3,9,12},{0,-1,-1},{1,-1,-1},{2,-1,-1},{3,-1,-1},{7,10,-1},{4,11,-1},{5,8,-1},{6,9,-1},{0,10,-1},{1,11,-1},{2,8,-1},{3,9,-1}};
uint32_t keyProb = 0;
// linear constraints requirements
vector<vector<uint32_t>> linearConstraints;
vector<uint32_t> linearKeys;
vector<uint32_t> propagatedLinearKeys;
vector<int> linearProbReduce; // contain the probability (in log2) for a linear constraint
vector<int> linearProb; // contain the probability (in log2) for a linear constraint
uint32_t **linearValues;
uint32_t *linearValuesLength;
uint64_t *linearPossibleKeys;

// non-linear constraints requirements
vector<vector<uint32_t>> nonLinearConstraints;
vector<vector<uint32_t>> nonLinearKeys;
vector<uint32_t> propagatedNonLinearKeys;
uint32_t **nonLinearValues;
uint32_t *nonLinearValuesLength;
uint32_t **nonLinearLast; // contain the count for the last Sbox transition (XDDT)
uint64_t *nonLinearPossibleKeys; // keeps track if a key is possible due to the constraint. (For counting possible keys)
uint16_t *combinedTimes; // number of "values" at the end of each tuple
uint32_t **nonLinearDistribution; // 2*jth entry keep track of the count, 2*j+1 keep track of the total. Divide to get prob
uint32_t *nonLinearDistributionLength; 


bool isLinear(vector<uint32_t> constraint)
{
	for (uint32_t i = 0; i < constraint.size();i++) if ((constraint[i] & ANDval) == 0) return false;
	return true;
}


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
	#if SIZE == 4
	int DDT4[16][16] = {0};
	computeDDT(DDT4);

	#elif SIZE == 8
	int DDT8[256][256] = {0};
	computeDDT8(DDT8);
	#endif

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
					#if SIZE == 4
						keyProb += -log2(DDT4[diffStates[n][0][4*r+c]][diffStates[n][1][4*r+c]]/16.0);
					#elif SIZE == 8
						keyProb += -log2(DDT8[diffStates[n][0][4*r+c]][diffStates[n][1][4*r+c]]/256.0);
					#endif
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
					for (uint32_t i = 0; i < tupTempLength[4*r+c]; i++) tup[4*r+c][i] = tupTemp[4*r+c][i];
					tupLength[4*r+c] = tupTempLength[4*r+c];
				}
				tupTempLength[4*r+c] = 0;
			}
		}
	}
	return tupsConstraints;
}

void printConstraint(vector<uint32_t> constraint)
{
	for (uint32_t j = 0; j < constraint.size(); j++)
	{
		cout << hex << "(" << (constraint[j] >> (4+2*SIZE)) << "," << ((constraint[j] >> (2*SIZE)) & 0xf) << ","
		<< ((constraint[j] >> SIZE) & ANDval) << "," << (constraint[j] & ANDval) << ") ";
	}
}

void printConstraints(vector<vector<uint32_t>> constraints)
{
	cout << "The constraints are (in the form of (round, position, input_diff, output_diff)" << endl;
	for (uint32_t i = 0; i < constraints.size(); i++)
	{
		printConstraint(constraints[i]);
		cout << isLinear(constraints[i]);
		cout << endl;
	}
}

void reduceKey(uint32_t* &val, uint32_t valLength) // this reduces the exact values, but according to proportions
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
				if (a[0] == 1) return; // terminate if it's 1
			}
		}
	}
 	for (int i = 0; i < (1<<(valLength*SIZE)); i++) val[i] /= a[0];
}

void XDDT(uint32_t first, uint32_t second, uint32_t output[1<<SIZE])
{
	#if SIZE == 4
		for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++)
		{
			for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++)
			{
				if (((v1^v2) == first) && ((getSbox(v1)^getSbox(v2)) == second)) output[v1]++;
			}
		}
	#elif SIZE == 8
		for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++)
		{
			for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++)
			{
				if (((v1^v2) == first) && ((getSbox8(v1)^getSbox8(v2)) == second)) output[v1]++;
			}
		}
	#else 
		#error Unsupported choice setting
	#endif
}

void YDDT(uint32_t first, uint32_t second, uint32_t output[1<<SIZE])
{
	#if SIZE == 4
		for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++)
		{
			for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++)
			{
				if (((v1^v2) == first) && ((getSbox(v1)^getSbox(v2)) == second)) output[getSbox(v1)]++;
			}
		}
	#elif SIZE == 8
		for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++)
		{
			for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++)
			{
				if (((v1^v2) == first) && ((getSbox8(v1)^getSbox8(v2)) == second)) output[getSbox8(v1)]++;
			}
		}
	#else 
		#error Unsupported choice setting
	#endif
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
		for (int j = 0; j < (1<<SIZE); j++) res[i^j] += val1[i] * val2[j];
	}
	for (int i = 0; i < (1<<SIZE); i++) val1[i] = res[i];
}

void addRoundConstant(uint32_t tup, uint32_t val[(1<<SIZE)])
{
	// return;
	uint32_t constantArray[(1<<SIZE)] = {0};
	if (((tup >> (SIZE*2)) & 0xf) == 0) constantArray[getConstants(tup >> (4+SIZE*2)) & 0xf] = 1;
	else if (((tup >> (SIZE*2)) & 0xf) == 4) constantArray[(getConstants(tup >> (4+SIZE*2)) >> 4) & 0b11] = 1;
	else if (((tup >> (SIZE*2)) & 0xf) == 8) constantArray[2] = 1;
	XOR(val,constantArray);
	reduceKey(val,1);
}

uint32_t resolveLinear(vector<uint32_t> constraint, uint32_t* &rightPartValues)
{
	#if SIZE == 4
	int DDT4[16][16] = {0};
	computeDDT(DDT4);

	#elif SIZE == 8
	int DDT8[256][256] = {0};
	computeDDT8(DDT8);
	#endif

	uint32_t key = 0;
	int t = constraint[constraint.size()-1];
	uint32_t xddt[1<<SIZE] = {0};
	XDDT((t >> SIZE) & ANDval,t & ANDval,xddt);
	for (int i = 0; i < (1<<SIZE); i++) {rightPartValues[i] = xddt[i];}
	for (uint32_t i = 0; i < constraint.size()-1; i++)
	{
		if (((constraint[i] >> (2*SIZE)) & 0xf) < 8) key = constraint[i];
		addRoundConstant(constraint[i],rightPartValues);
		uint32_t yddt[1<<SIZE] = {0};
		YDDT((constraint[i] >> SIZE) & ANDval,constraint[i] & ANDval,yddt);
		XOR(rightPartValues,yddt);
	}	
	int totalCount = 0;
	int singleValue = 0;
	for (int i = 0; i < 1<<(SIZE); i++)
	{
		totalCount += rightPartValues[i];
		if (rightPartValues[i] > 0) singleValue = rightPartValues[i];
	}
	#if SIZE == 4
	linearProb.push_back(-log2((DDT4[(constraint[constraint.size()-1] >> SIZE) & ANDval][constraint[constraint.size()-1] & ANDval]+0.0)/(1<<SIZE)));
	#elif SIZE == 8
	linearProb.push_back(-log2((DDT8[(constraint[constraint.size()-1] >> SIZE) & ANDval][constraint[constraint.size()-1] & ANDval]+0.0)/(1<<SIZE)));
	#endif
	linearProbReduce.push_back(SIZE-int(-log2((singleValue+0.0)/(totalCount))));
	reduceKey(rightPartValues,1);
	return key;
}

vector<uint32_t> resolveLeft(uint32_t zeroPoint, vector<uint32_t> constraint, uint32_t values[(1<<SIZE)], vector<uint32_t> &zeroDependencies, bool last = false)
{
	vector<uint32_t> keys;
	uint32_t round = (zeroPoint >> (4+2*SIZE));
	uint32_t pos = (zeroPoint >> (2*SIZE)) & 0xf;
	int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
	for (uint32_t i = 0; i < constraint.size(); i++)
	{
		uint32_t constraintRound = constraint[i] >> (4+2*SIZE);
		int constraintPos = (constraint[i] >> (2*SIZE)) & ROUNDval;
		// check if the value is involved in this zeroPoint
		if (constraintRound != round-1) continue;
		if (constraintPos != inversePos[0] && constraintPos != inversePos[1] && constraintPos != inversePos[2]) continue;
		if (constraintPos < 8) keys.push_back(constraint[i]);
		if ((constraint[i] & ANDval) == 0) zeroDependencies.push_back(constraint[i]);
		else if (last == false) // if last is true. follows resolveLast (i.e. ignore yddt here, but settling outside)
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

bool AddKeyNonLinear(uint32_t* &value,uint32_t &valueLength)
{
	// this function extends valueLength by 1
	uint32_t keySize = (1ULL << ((valueLength-1) * SIZE)); // the last "SIZE" is dedicated for values, not key
	uint64_t destinationSize = (1ULL << ((valueLength+1)*SIZE));
	if ((valueLength+1)*SIZE >= 36) 
	{
		cout << "this constraint exceed 2^36 space. Skipped this constraint." << endl;
		return false;
	}
	uint32_t* destination; destination = new uint32_t [destinationSize]();
	#if SIZE == 4
		for (uint64_t i = 0; i < keySize; i++)
		{
			for (int k = 0; k < 16; k++)
			{
				for (int v = 0; v < 16; v++)
				{
					destination[i*(16*16)+k*16+v] += value[i * 16 + (v^k)];
				}
			}
		}
	#elif SIZE == 8
		for (uint64_t i = 0; i < keySize; i++)
		{
			for (int k = 0; k < 256; k++)
			{
				for (int v = 0; v < 256; v++)
				{
					destination[i*(256*256)+k*256+v] += value[i * 256 + (v^k)];
				}
			}
		}
	#endif
		delete[] value;
		value = destination;
		valueLength++;
	return true;
}

void SubstNonLinear(uint32_t* &value,uint32_t &valueLength)
{
	uint32_t destinationSize = (1ULL << ((valueLength)*SIZE));
	uint32_t* destination; destination = new uint32_t [destinationSize]();
	#if SIZE == 4
		for (uint64_t i = 0; i < (1ULL<<((valueLength-1)*SIZE)); i++)
		{
			for (int v = 0; v < 16; v++)
			{
				destination[i*16+getSbox(v)] += value[i*16 + v];
			}
		}

	#elif SIZE == 8
		for (uint64_t i = 0; i < (1ULL<<((valueLength-1)*SIZE)); i++)
		{
			for (int v = 0; v < 256; v++)
			{
				destination[i*256+getSbox8(v)] += value[i*256 + v];
			}
		}
	#endif
	delete[] value;
	value = destination;
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


void AddFinalYDDT(vector<uint32_t> constraint, uint32_t* &val, uint32_t valLength)
{
	uint32_t last = constraint[constraint.size()-1];
	uint32_t *valLast; valLast = new uint32_t[1<<SIZE]();
	uint32_t round = (last >> (4+2*SIZE));
	uint32_t pos = (last >> (2*SIZE)) & 0xf;
	int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
	for (uint32_t i = 0; i < constraint.size(); i++)
	{
		uint32_t constraintRound = constraint[i] >> (4+2*SIZE);
		int constraintPos = (constraint[i] >> (2*SIZE)) & ROUNDval;
		// check if the value is involved in this zeroPoint
		if (constraintRound != round-1) continue;
		if (constraintPos != inversePos[0] && constraintPos != inversePos[1] && constraintPos != inversePos[2]) continue;
		if ((constraint[i] & ANDval) == 0) continue;
		uint32_t yddt[1<<SIZE] = {0};
		uint32_t tmp[1<<SIZE] = {0};
		YDDT((constraint[i] >> SIZE) & ANDval,constraint[i] & ANDval,yddt);
		addRoundConstant(constraint[i],valLast);
		XOR(valLast,yddt);
		XOR(valLast,tmp);
		reduceKey(valLast,1);
	}
	uint32_t *tmp; tmp = new uint32_t[1ULL<<(SIZE*valLength)]();
	// uint32_t s = 0;
	// for (uint64_t i = 0; i < (1ULL<<(valLength*SIZE)); i++) s += val[i];
	for (int j = 0; j < (1<<SIZE); j++)
	{
		if (valLast[j] == 0) continue;
		for (uint64_t i = 0; i < (1ULL<<(valLength*SIZE)); i++)
		{
			tmp[uint64_t(i^j)] += val[i];
		}
	}
	val = tmp;
}

bool resolveNonLinear(vector<uint32_t> constraint, uint32_t* &values, uint32_t &valuesLength, uint32_t lastValues[1<<SIZE], vector<uint32_t> &output)
{
	// first, we need to locate the zeros (which determine how many rounds/ connection points)
	sort(constraint.begin(),constraint.end());
	vector<uint32_t> zeroPoints; 
	for (uint32_t i = 0; i < constraint.size(); i++) 
	{
		if ((constraint[i] & ANDval) == 0) zeroPoints.push_back(constraint[i]);
	}
	uint32_t* zeroValues[zeroPoints.size()] = {0};
	uint32_t zeroValuesLength[zeroPoints.size()] = {0}; // this will tell me how many keys are already involved
	vector<vector<uint32_t>> zeroKeys;
	for (uint32_t i = 0; i < zeroPoints.size(); i++)
	{
		vector<uint32_t> zeroDependencies;
		uint32_t* val;
		val = new uint32_t[(1<<SIZE)]();
		uint32_t valLength = 1; // by default, the value is 1. Unless there are more than 1 zeroDependencies adding to the same zeroPoint
		vector<uint32_t> keys = resolveLeft(zeroPoints[i],constraint, val, zeroDependencies); // zeroDependencies are the ones that require us to add the dimension of key
		// settle if there is any zeroDependencies
		for (uint32_t j = 0; j < zeroDependencies.size(); j++)
		{
			uint32_t index;
			for (uint32_t k = 0; k < zeroPoints.size(); k++) { if (zeroPoints[k] == zeroDependencies[j]) index = k;}
			resolveZeroDependencies(zeroValues[index],zeroValuesLength[index],val,valLength);
			keys.insert(keys.begin(),zeroKeys[index].begin(),zeroKeys[index].end());
		}
		/* for AddKeyAnsSubst 
		what this does is to put the fixed key at all the positions except for the last. Which
		corresponds to the value that it will produce after the Sbox
		for example, if val has nonzero entries at 1,3,8,a, then the output will be in the form of
		(k,S[k^1] for all k),...,(k,S[k^a] for all k)
		*/
		if (!AddKeyNonLinear(val,valLength)) {return false;}
		SubstNonLinear(val,valLength);
		zeroKeys.push_back(keys);
		zeroValues[i] = val;
		zeroValuesLength[i] = valLength;
	}
	// Settle the last tuple!
	uint32_t* val; val = new uint32_t[1<<SIZE]();
	uint32_t valLength = 1;
	vector<uint32_t> zeroDependencies;
	vector<uint32_t> keys = resolveLeft(constraint[constraint.size()-1], constraint, lastValues, zeroDependencies, true); // last argument changed it from resolveLeft to resolveLast
	// putting lastValues using XDDT
	uint32_t xddt[1<<SIZE] = {0};
	XDDT((constraint[constraint.size()-1] >> SIZE) & ANDval,constraint[constraint.size()-1] & ANDval,xddt);
	XOR(lastValues,xddt);

	if (zeroDependencies.size() == 0) // sanity check
	{
		cout << "zeroDependencies cannot be zero at this point" << endl;
		exit(0);
	}
	for (uint32_t j = 0; j < zeroDependencies.size(); j++)
	{
		uint32_t index;
		for (uint32_t k = 0; k < zeroPoints.size(); k++) { if (zeroPoints[k] == zeroDependencies[j]) index = k;}
		// at this point, we should capture the info: val gives the keys that we want, zeroValues has the entire distribution.
		resolveZeroDependencies(zeroValues[index],zeroValuesLength[index],val,valLength);
		keys.insert(keys.begin(),zeroKeys[index].begin(),zeroKeys[index].end());
	}
	if (!AddKeyNonLinear(val,valLength)) {return false;}
	AddFinalYDDT(constraint,val,valLength);
	values = val;
	valuesLength = valLength;
	output = keys;
	return true;
}

bool combineNonLinearConstraints(uint32_t index_i, uint32_t index_j)
{
	vector<uint32_t> combinedKeys;
	vector<uint32_t> index2; // index2 contains the index of the key in tmpValue that belongs to nonLinearKeys[index_j]
	for (uint32_t i = 0; i < nonLinearKeys[index_i].size(); i++) combinedKeys.push_back(nonLinearKeys[index_i][i]);
	for (uint32_t j = 0; j < nonLinearKeys[index_j].size(); j++)
	{
		auto it = find(combinedKeys.begin(), combinedKeys.end(), nonLinearKeys[index_j][j]);
		if (it == combinedKeys.end()) combinedKeys.push_back(nonLinearKeys[index_j][j]);
		else index2.push_back(it - combinedKeys.begin());
	}

	if ((combinedKeys.size()+combinedTimes[index_i]+combinedTimes[index_j])*SIZE > 36) 
	{
		cout << "Warning: space required exceeded 2^36 entries!" << endl;
		cout << "Please do an experiment instead!" << endl;
		exit(0);
		return false;
	}
	// allocating the necessary array
	uint32_t tmpValuesLength = combinedKeys.size()+combinedTimes[index_i]+combinedTimes[index_j];
	uint32_t *tmpValues = new uint32_t[(1UL<<(tmpValuesLength*SIZE))]();


	// the main combining loop
	for (uint64_t i = 0; i < (1ULL << (tmpValuesLength*SIZE)); i++)
	{
		uint64_t firstIndex = (i >> ((tmpValuesLength-nonLinearValuesLength[index_i]+1)*SIZE));
		firstIndex = (firstIndex << (SIZE*combinedTimes[index_i])) + ((i >> (SIZE*combinedTimes[index_j])) & ((1<<(SIZE*combinedTimes[index_i]))-1));
		uint64_t secondIndex = 0;
		for (uint32_t j = 0; j < index2.size(); j++)
		{
			secondIndex = secondIndex << SIZE;
			secondIndex += (i >> ((tmpValuesLength-index2[j]-1)*SIZE)) & ((1<<SIZE)-1);
		}
		secondIndex = (secondIndex << (SIZE*combinedTimes[index_j])) + (i & ((1<<(SIZE*combinedTimes[index_j]))-1));
		tmpValues[i] = nonLinearValues[index_i][firstIndex] * nonLinearValues[index_j][secondIndex];
	}
	nonLinearValues[index_i] = tmpValues;
	nonLinearValuesLength[index_i] = tmpValuesLength;
	// changing the nonLinearLast
	uint32_t *newLastValues = new uint32_t[1UL<<(SIZE*(combinedTimes[index_i]+combinedTimes[index_j]))]();
	for (int i = 0; i < (1 << (SIZE*combinedTimes[index_i])); i++)
	{
		for (int j = 0; j < (1 << (SIZE*combinedTimes[index_j])); j++)
		{
			newLastValues[(i<<(SIZE*combinedTimes[index_j]))+j] = nonLinearLast[index_i][i] * nonLinearLast[index_j][j];
		}
	}
	nonLinearLast[index_i] = newLastValues;

	// update combinedTimes
	combinedTimes[index_i] += combinedTimes[index_j];

	// deleting arrays
	for (uint32_t j = index_j; j < nonLinearKeys.size()-1; j++)
	{
		nonLinearValues[j] = nonLinearValues[j+1];
		nonLinearValuesLength[j] = nonLinearValuesLength[j+1];
		nonLinearLast[j] = nonLinearLast[j+1];
		combinedTimes[j] = combinedTimes[j+1];
	}
	delete[] nonLinearValues[nonLinearKeys.size()-1];
	nonLinearValuesLength[nonLinearKeys.size()-1] = 0;
	delete[] nonLinearLast[nonLinearKeys.size()-1];
	combinedTimes[nonLinearKeys.size()-1] = 0;


	nonLinearKeys.erase(nonLinearKeys.begin()+index_i);
	nonLinearKeys.insert(nonLinearKeys.begin()+index_i,combinedKeys);
	nonLinearKeys.erase(nonLinearKeys.begin()+index_j);
	nonLinearConstraints.erase(nonLinearConstraints.begin()+index_j);
	return true;
}

// bool intersect(vector<uint32_t> constraint1, vector<uint32_t> constraint2, bool constraint1Linear, bool constraint2Linear)
// {
// 	// if this intersection function can tell which kind of intersection
// 	for (uint32_t i = 0; i < constraint1.size(); i++)
// 	{
// 		if ((constraint1Linear == true) && (i == constraint1.size()-1)) continue;
// 		for (uint32_t j = 0; j < constraint2.size(); j++)
// 		{
// 			if ((constraint1[i] >> (2*SIZE)) == (constraint2[i] >> (2*SIZE))) 
// 			{
// 				if ((i == 0 && j == constraint2.size()-1) || (j == 0 && i == constraint1.size()-1))
// 				{
// 					cout << "intersection is of type 2 or 3" << endl;
// 				}
// 				return true;
// 			}
// 		}
// 	}
// 	return false;
// }

void combineConstraints()
{
	// use linear values to eliminate nonlinear values
	for (uint16_t i = 0; i < nonLinearKeys.size(); i++)
	{
		uint16_t j = 0;
		while (j < linearKeys.size())
		{
			bool involved = false;
			for (uint16_t k = 0; k < nonLinearKeys[i].size(); k++)
			{
				if ((linearKeys[j] >> (2*SIZE)) == (nonLinearKeys[i][k] >> (2*SIZE))) // found an intersection, we remove it from nonLinearValues
				{
					involved = true;
					for (uint64_t l1 = 0; l1 < (1<<SIZE); l1++) // loop through all possible key values
					{
						if (linearValues[j][l1] == 0)
						{
							for (uint64_t l2 = 0; l2 < (1ULL << (SIZE*(nonLinearValuesLength[i]-1))); l2++)
							{
								uint64_t msn = (l2 >> (SIZE*(nonLinearValuesLength[i] - k - 1))) << (SIZE*(nonLinearValuesLength[i] - k));
								uint64_t lsn = l2 & ((1ULL << (SIZE*(nonLinearValuesLength[i] - k - 1))) -1);
								nonLinearValues[i][msn + (uint64_t(l1) << (SIZE*(nonLinearValuesLength[i] - k - 1))) + lsn] = 0;
							}
						}
						else
						{
							bool zeroFlag = true;
							for (uint64_t l2 = 0; l2 < (1ULL << (SIZE*(nonLinearValuesLength[i]-1))); l2++)
							{
								uint64_t msn = (uint64_t(l2) >> (SIZE*(nonLinearValuesLength[i] - k - 1))) << (SIZE*(nonLinearValuesLength[i] - k));
								uint64_t lsn = l2 & ((1ULL << (SIZE*(nonLinearValuesLength[i] - k - 1))) -1);
								if (nonLinearValues[i][msn+(uint64_t(l1) << (SIZE*(nonLinearValuesLength[i] - k - 1)))+lsn] > 0) 
								{
									zeroFlag = false;
									break;
								}
							}
							if (zeroFlag) linearValues[j][l1] = 0;
						}
					}
				}
			}
			if (involved)
			{
				for (uint16_t l = j; l < linearKeys.size()-1; l++)
				{
					linearValues[l] = linearValues[l+1];
					linearValuesLength[l] = linearValuesLength[l+1];
				}
				linearValuesLength[linearKeys.size()-1] = 0;
				linearKeys.erase(linearKeys.begin()+j);
				linearConstraints.erase(linearConstraints.begin()+j);
				linearProbReduce.erase(linearProbReduce.begin()+j);
				linearProb.erase(linearProb.begin()+j);
				j--;
			}
			// else if (intersect(linearConstraints[j],nonLinearConstraints[i],true,false))
			// {
			// 	cout << "linear constraint " << j << " and nonlinear constraint " << i << " are not independent, but treating them as independent" << endl;
			// }
			j++;
		}
	}

	// linear with linear
	for (uint32_t i = 0; i < linearKeys.size(); i++)
	{
		uint32_t j = i + 1;
		while (j < linearKeys.size())
		{
			if ((linearKeys[i] >> (2*SIZE)) == (linearKeys[j] >> (2*SIZE)))
			{
				for (uint32_t k = 0; k < (1 << SIZE); k++)
				{
					if (linearValues[j][k] == 0) linearValues[i][k] = 0;
				}
				for (uint16_t l = j; l < linearKeys.size()-1; l++)
				{
					linearValues[l] = linearValues[l+1];
					linearValuesLength[l] = linearValuesLength[l+1];
				}
				linearValuesLength[linearKeys.size()-1] = 0;
				linearKeys.erase(linearKeys.begin()+j);
				linearConstraints.erase(linearConstraints.begin()+j);
				linearProbReduce.erase(linearProbReduce.begin()+j);
				linearProb.erase(linearProb.begin()+j);
				j--;
			}
			j++;
		}
	}
	// nonlinear with nonlinear
	// work here

	uint32_t i = 0;
	uint32_t j = 0;
	while (i < nonLinearKeys.size())
	{
		j = i + 1;
		while (j < nonLinearKeys.size())
		{
			bool combined = false;
			for (uint32_t k0 = 0; k0 < nonLinearKeys[i].size(); k0++)
			{
				for (uint32_t k1 = 0; k1 < nonLinearKeys[j].size(); k1++)
				{
					if ((nonLinearKeys[i][k0] >> (2*SIZE)) == ((nonLinearKeys[j][k1]) >> 2*SIZE)) // if the round and pos are the same. Combine them
					{
						combined = combineNonLinearConstraints(i,j);
						if (combined)
						{
							j -= 1;	
							break;
						}
					}
				}
				if (combined) break;
			}
			// if (combined == false && intersect(nonLinearConstraints[i],nonLinearConstraints[j],false,false))
			// {
			// 	cout << "nonlinear constraint " << i << " and nonlinear constraint " << j << " are not independent, but treating them as independent" << endl;
			// }
			j++;
		}
		i++;
	}

}

void computeNonLinearDistribution()
{
	for (uint32_t i = 0; i < nonLinearConstraints.size(); i++)
	{
		nonLinearDistribution[i] = new uint32_t[1ULL<<(((nonLinearValuesLength[i]-combinedTimes[i])*SIZE)+1)]();
		nonLinearDistributionLength[i] = nonLinearValuesLength[i]-combinedTimes[i];
		for (uint64_t j = 0; j < (1ULL<<(nonLinearDistributionLength[i]*SIZE)); j++)
		{
			for (uint32_t k = 0; k < (1UL<<(combinedTimes[i]*SIZE)); k++)
			{
				if (nonLinearLast[i][k] > 0) nonLinearDistribution[i][2*j] += nonLinearValues[i][(j<<(SIZE*combinedTimes[i]))+k]; // the count
				nonLinearDistribution[i][2*j+1] += nonLinearValues[i][(j<<(SIZE*combinedTimes[i]))+k]; // the total
			}
			if (nonLinearDistribution[i][2*j] > 0) nonLinearPossibleKeys[i]++; // only add in when there is a possible value due to the keys
		}
	}

}

void computeLinearPossibleKeys()
{
	for (uint32_t i = 0; i < linearConstraints.size(); i++)
	{
		for (uint32_t j = 0; j < (1<<SIZE); j++)
		{
			if (linearValues[i][j] > 0) linearPossibleKeys[i]++;
		}
		if (linearPossibleKeys[i] == 0)
		{
			cout << "there are no keys satisfying the linear constraints!" << endl;
			exit(0);
		}
	}
}

double getOriginalProbability()
{
	double prob = 0;
	uint8_t key_diff_temp[4][4][4];
	for (int i = 0; i < 4; i++){ for (int j = 0; j < 4; j++) { for (int k = 0; k < 4; k++) key_diff_temp[i][j][k] = key_diff[i][j][k];}}
	uint8_t before[4][4], after[4][4];

	#if SIZE == 4
	int DDT4[16][16] = {0};
	computeDDT(DDT4);

	#elif SIZE == 8
	int DDT8[256][256] = {0};
	computeDDT8(DDT8);
	#endif

	for (uint32_t n = 0; n < nr; n++) // this loop computes the original prob
	{
		for (uint8_t r = 0; r < 4; r++)
		{
			for (uint8_t c = 0; c < 4; c++)
			{
				before[r][c] = alpha[n][r][c];
				after[r][c] = alpha[n+1][r][c];
			}
		}
		invMC(after); invSR(after); AddRoundTweakey(after,key_diff_temp);
		for (uint8_t r = 0; r < 4; r++)
		{
			for (uint8_t c = 0; c < 4; c++)
			{
				if (before[r][c] == 0) continue;
				#if SIZE == 4
					prob += -log2(DDT4[before[r][c]][after[r][c]]/16.0);
				#elif SIZE == 8
					prob += -log2(DDT8[before[r][c]][after[r][c]]/256.0);
				#endif
			}
		}
		#if SIZE == 4
			key_schedule_round_function(key_diff_temp);
		#elif SIZE == 8
			key_schedule_round_function8(key_diff_temp);
		#endif
	}
	return prob;
}


void remove(vector<int> &v) // only used for getSubgraph
{
    auto itr = v.begin();
    unordered_set<int> s;
 
    for (auto curr = v.begin(); curr != v.end(); ++curr)
    {
        if (s.insert(*curr).second) {
            *itr++ = *curr;
        }
    }
    v.erase(itr, v.end());
}

vector<vector<int>> getSubgraph(vector<vector<int>> G, bool keysDetermined[16], vector<int> &keyPositions)
{
	// find the start point
	uint8_t start = 0;
	while (keysDetermined[start] != 0 && start < 16) start++;
	vector<int> allConstraints = G[start];
	keyPositions.push_back(start);
	keysDetermined[start] = 1;
    vector<vector<int>> constraintIndex; constraintIndex.push_back(allConstraints);
    
	uint8_t i = 0;
	while (i < allConstraints.size())
	{
		for (uint8_t j = start+1; j < 16; j++)
		{
			if (keysDetermined[j] == 1) continue; // visited
			auto it = find(G[j].begin(),G[j].end(),allConstraints[i]);
			if (it == G[j].end()) continue;
			allConstraints.insert(allConstraints.end(),G[j].begin(),G[j].end());
			constraintIndex.push_back(G[j]);
			keyPositions.push_back(j);
			keysDetermined[j] = 1;
			remove(allConstraints); // remove duplicates
		}
		i++;
	}
	return constraintIndex;
}

double subGraphKey(vector<int> keyPositions, vector<vector<int>> constraintsIndex)
{
	// first, find how many 
	double currentKeys = 0;
	vector<int> accountedFor;
	for (uint32_t i = 0; i < keyPositions.size(); i++)
	{
		if (constraintsIndex[i].size() > TK_NUM)
		{
			cout << "Intersection exceeded more than TK size. Please go for experiment" << endl;
			exit(0);
			// if SAT works, then there exists exactly one key
			currentKeys += 0;
		}
		if (constraintsIndex[i].size() < TK_NUM) currentKeys += (TK_NUM-constraintsIndex[i].size()) * SIZE; // we have to expand if there is no multiplication. implicit
		for (uint32_t j = 0; j < constraintsIndex[i].size(); j++)
		{
			auto it = find(accountedFor.begin(),accountedFor.end(),constraintsIndex[i][j]);
			if (it != accountedFor.end()) continue;
			else if (constraintsIndex[i][j] >= 0) currentKeys += log2(linearPossibleKeys[constraintsIndex[i][j]]);
			else currentKeys += log2(nonLinearPossibleKeys[-constraintsIndex[i][j]-1]); 
		}
		// currentKeys += (TK_NUM-1) * SIZE;
		accountedFor.insert(accountedFor.end(),constraintsIndex[i].begin(), constraintsIndex[i].end());
		remove(accountedFor);
	}
	currentKeys -= (TK_NUM-1) * SIZE * keyPositions.size(); // note that we minus the TK, so we can compare the size in terms of the "XORed" keys
	return currentKeys;
}



double computeKeys()
{
	bool keysDetermined[16] = {0}; // this determines if the keys have been accounted for (use permSchedule to get to round nr)
	double tmpKeys = 0;
	uint32_t n; uint32_t pos;

	// we first compute the number of keys
	vector<vector<int>> propagatedCell; // this keep track of what cell is affected by which constraint
	for (uint32_t i = 0; i < 16; i++)
	{
		vector<int> tmp;
		propagatedCell.push_back(tmp);
	}
	for (uint32_t i = 0; i < linearKeys.size(); i++)
	{
		n = linearKeys[i] >> (2*SIZE+4);
		pos = (linearKeys[i] >> (2*SIZE)) & 0xf;
		for (uint32_t k = n; k < nr; k++) pos = getPermSchedule2(pos);
		propagatedCell[pos].push_back(i);
	}
	for (uint32_t i = 0; i < nonLinearKeys.size(); i++)
	{
		for (uint32_t j = 0; j < nonLinearKeys[i].size(); j++)
		{
			n = nonLinearKeys[i][j] >> (2*SIZE+4);
			pos = (nonLinearKeys[i][j] >> (2*SIZE)) & 0xf;
			for (uint32_t k = n; k < nr; k++) pos = getPermSchedule2(pos);
			propagatedCell[pos].push_back(-i-1);
		}
	}
	// propagatedCell is filled at this point
	// Finding a subgraph using propagatedCell
	// 2 vectors in total. The first one will record the set of constraints index. nonlinear constraints start at negatives. -1 --> 0, -2 --> 1, ...
	// The second vector will record the set of key cells involved
	bool flag = true;
	while (flag)
	{
		vector<int> keyPositions;
		vector<vector<int>> constraintsIndex = getSubgraph(propagatedCell, keysDetermined, keyPositions);  
		if ((keyPositions.size() == 1) && (constraintsIndex[0].size() == 0)) tmpKeys += SIZE;
		else tmpKeys += subGraphKey(keyPositions,constraintsIndex);
		// to check if we continue
		flag = false;
		for (int i = 0; i < 16; i++) {if (keysDetermined[i] == 0) flag = true;} // continue until all the key positions are exhausted
	}
	
	return tmpKeys;
}

void computeProbability()
{
	double prob = getOriginalProbability();
	cout << "original prob: " << prob << endl;
	#if SIZE == 4
	int DDT4[16][16] = {0};
	computeDDT(DDT4);

	#elif SIZE == 8
	int DDT8[256][256] = {0};
	computeDDT8(DDT8);
	#endif

	prob -= keyProb; // independent of the ones involved in key constraints (nonlinear constraints only)

	for (uint32_t i = 0; i < linearKeys.size(); i++)
	{
		prob += linearProb[i] - linearProbReduce[i];
	}
	if (nonLinearKeys.size() == 0)
	{
		cout << "probability: " << prob << endl;
		return;
	}
	uint64_t **probArray; probArray = new uint64_t*[nonLinearKeys.size()]();
	uint64_t *probLength; probLength = new uint64_t[nonLinearKeys.size()]();
	uint64_t *probBase; probBase = new uint64_t[nonLinearKeys.size()]();
	for (uint32_t i = 0; i < nonLinearKeys.size(); i++)
	{
		uint64_t max = 0;
		uint64_t base = 1;
		vector<uint64_t> nonZerosProb;
		for (uint64_t j = 0; j < (1ULL<<(nonLinearDistributionLength[i]*SIZE)); j++)	
		{
			if (nonLinearDistribution[i][2*j] > 0) 
			{
				auto it = find(nonZerosProb.begin(),nonZerosProb.end(),nonLinearDistribution[i][2*j]);
				if (it == nonZerosProb.end()) nonZerosProb.push_back(nonLinearDistribution[i][2*j]);
				if (nonLinearDistribution[i][2*j] > max) max = nonLinearDistribution[i][2*j];
				base = nonLinearDistribution[i][2*j+1];
			}
		}
		if (nonZerosProb.size() == 0) 
		{
			cout << "probability: " << 0 << endl;
		return;
		}
		uint64_t divisor = nonZerosProb[0];
		if (divisor == 0)
		{
			cout << "divisor is zero!" << endl;
			exit(0);
		}
		for (uint32_t j = 1; j < nonZerosProb.size(); j++) divisor = gcd(divisor, nonZerosProb[j]);
		probLength[i] = max/divisor+1;
		probArray[i] = new uint64_t[probLength[i]]();
		base = base/divisor;
		for (uint64_t j = 0; j < (1ULL<<(nonLinearDistributionLength[i]*SIZE)); j++)	
		{
			if (nonLinearDistribution[i][2*j] == 0) continue;
			probArray[i][nonLinearDistribution[i][2*j]/divisor]++;
		}
		vector<uint64_t> nonZerosNumbers;
		for (uint64_t j = 0; j < probLength[i]; j++)
		{
			if (probArray[i][j] > 0) nonZerosNumbers.push_back(probArray[i][j]);
		}
		divisor = nonZerosNumbers[0];
		for (uint32_t j = 1; j < nonZerosNumbers.size(); j++) divisor = gcd(divisor, nonZerosNumbers[j]);
		for (uint64_t j = 0; j < probLength[i]; j++) probArray[i][j] /= divisor;
		probBase[i] = base;
	}

	// combine the nonlinear
	uint64_t combinedProbBase = probBase[0];
	uint64_t combinedMax = probLength[0];
	for (uint32_t i = 1; i < nonLinearKeys.size(); i++) 
	{
		combinedProbBase *= probBase[i]; 
		combinedMax *= probLength[i];
	}
	uint64_t *combinedProbArray; combinedProbArray = new uint64_t[combinedMax]();
	uint64_t *indices = new uint64_t[nonLinearKeys.size()]();
	while (true)
    {
        bool flag = true;
        for (uint32_t i = 0; i < nonLinearKeys.size(); i++) { if (probArray[i][indices[i]] == 0) flag = false;}
        if (flag)
        {
        	int index = indices[0];
        	uint64_t value = probArray[0][indices[0]];
        	for (uint32_t i = 1; i < nonLinearKeys.size(); i++)
        	{
        		index *= indices[i];
        		value *= probArray[i][indices[i]];
        	}
            combinedProbArray[index] += value;
        }
        int indexToMove = nonLinearKeys.size()-1;
        while ((indices[indexToMove] >= probLength[indexToMove]-1) && (indexToMove >= 0)) indexToMove--;
        if (indexToMove == -1) break;
        indices[indexToMove]++;
        for (uint32_t i = indexToMove+1; i < nonLinearKeys.size(); i++) indices[i] = 0;
    }
	double sum = 0;
	for (uint32_t i = 1; i < combinedMax; i++) sum += combinedProbArray[i];
	cout << "The distribution: " << endl;
	for (uint32_t i = 1; i < combinedMax; i++)
	{
		if (combinedProbArray[i] > 0)
		{
			cout << prob - log2((i+0.0)/combinedProbBase) << " : " << combinedProbArray[i]/(sum)*100 << endl;
		}
	}
}

void propagateConstraints()
{

	for (int i = 0; i < linearKeys.size(); i++)
	{
		uint32_t key = linearKeys[i];
		while ((key >> (2*SIZE + 4)) < nr) key = ((key >> (2*SIZE + 4)) + 1) << (2*SIZE+4) ^ (getPermSchedule2((key >> (2*SIZE))&0xf) << (2*SIZE)) ^ (key & ((1<<(2*SIZE))-1));
		linearKeys[i] = key;
	}

	for (int i = 0; i < nonLinearKeys.size(); i++)
	{
		for (int j = 0; j < nonLinearKeys[i].size(); j++)
		{
			uint32_t key = nonLinearKeys[i][j];
			while ((key >> (2*SIZE + 4)) < nr) key = ((key >> (2*SIZE + 4)) + 1) << (2*SIZE+4) ^ (getPermSchedule2((key >> (2*SIZE))&0xf) << (2*SIZE)) ^ (key & ((1<<(2*SIZE))-1));
			nonLinearKeys[i][j] = key;
		}
	}
}

bool sizeCheck(vector<uint32_t> constraint)
{
	int key = 0;
	for (uint32_t i = 0; i < constraint.size()-1; i++)
	{
		if (((constraint[i] >> (2*SIZE)) & 0xf) < 8) key++;
	}
	if ((key+1) * SIZE > 36) return false;
	return true;
}



int main(int argc, char *argv[])
{
	// SKINNY64
	SK_2020_1402(alpha,key_diff,nr,TK_NUM);
	// TK1_2020_1402(alpha,key_diff,nr,TK_NUM); // impossible in liFnear
	// TK2_2020_1402(alpha,key_diff,nr,TK_NUM); 
	// TK3_2020_1402(alpha,key_diff,nr,TK_NUM);
	// get_4TK2_1_17_L_2020_1317(alpha,key_diff,nr,TK_NUM); // this is used as a first test
	// get_4TK2_1_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // without nonLinearConstraint[0], it's possible
	// get_4TK2_1_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// get_4TK2_1_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // without nonLinearConstraint[0], it's possible
	// get_4TK2_2_17_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// get_4TK2_2_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// get_4TK2_2_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// get_4TK2_2_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// get_4TK3_1_22_L_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// get_4TK3_1_22_U_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// get_4TK3_1_22_L_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	// get_4TK3_1_22_U_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// get_4TK3_1_22_U_BMD3_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// get_4TK3_1_23_U_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	// get_4TK3_1_23_U_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible linear

	// get_4TK2_18_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // u = 6, v = 2,10,12,13 // all ok
	// get_4TK3_22_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // u = 8,9, v = 5 // all ok


	// 8 bit trails
	// SKINNY128
	// get_8TK2_1_18_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_1_18_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_1_19_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_1_20_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_1_20_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_1_21_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// SK_tosc_2017_i4_99_129(alpha,key_diff,nr,TK_NUM);
	// TK1_8_2020_1402(alpha,key_diff,nr,TK_NUM);
	// TK2_8_2020_1402(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_2_21_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_22_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_22_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_23_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_23_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_24_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_24_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_25_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_25_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_2_18_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_2_18_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_2_19_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_2_19_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_2_20_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_25_L_2020_1317(alpha,key_diff,nr,TK_NUM);

	// uint8_t u1[4]= {0x40,0x50,0xc0,0xd0};  uint8_t v1[11] =  {0x2a,0x2e,0x2f,0x6a,0xba,0xbe,0xbf,0xea,0xef,0xfa,0xfe};
	// uint8_t w1[18] = {0x2a,0x2d,0x2e,0x2f,0x6b,0x7c,0xb8,0xb9,0xba,0xbd,0xbe,0xbf,0xeb,0xec,0xef,0xfb,0xfc,0xfe};
	// u = 0x40,0x50,0xc0,0xd0  v =  0x2a,0x2e,0x2f,0x6a,0xba,0xbe,0xbf,0xea,0xef,0xfa,0xfe
	// w = 0x2a,0x2d,0x2e,0x2f,0x6b,0x7c,0xb8,0xb9,0xba,0xbd,0xbe,0xbf,0xeb,0xec,0xef,0xfb,0xfc,0xfe
	// uint8_t u = u1[q % 4];
	// uint8_t v = v1[(q / 4) % 11];
	// uint8_t w = w1[(q / 44)];
 	// get_8TK2_19_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v,w); // all ok
	// uint8_t u1[5] = {0x08,0x0c,0x0d,0x0e,0x0f};
	// uint8_t v1[7] = {0x28, 0x29, 0x2c, 0x2d, 0x2e, 0x2f, 0xb9};
	// uint8_t u = u1[q % 5];
	// uint8_t v = v1[(q/5)];
	
	// u = 0x08,0x0c,0x0d,0x0e,0x0f
	// v = 0x28, 0x29, 0x2c, 0x2d, 0x2e, 0x2f, 0xb9	
	// get_8TK3_22_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // all ok
	
	fixedValuesToPath(alpha,key_diff,diffStates,fixedStates,nr);
	vector<vector<uint32_t>> constraints = findConstraints(diffStates,fixedStates,nr);
	// printState(diffStates,fixedStates,nr);
	printConstraints(constraints);
	for (uint32_t i = 0; i < constraints.size(); i++) // split linear and nonlinear constraints
	{
		if (isLinear(constraints[i])) linearConstraints.push_back(constraints[i]); 
		else if (!isLinear(constraints[i]))
		{
			if (sizeCheck(constraints[i])) nonLinearConstraints.push_back(constraints[i]);
		}
		else {cout << "error!" << endl;}
	}
	// work on the linear constraints
	linearValues = new uint32_t*[linearConstraints.size()]();
	linearValuesLength = new uint32_t[linearConstraints.size()]();
	linearPossibleKeys = new uint64_t[linearConstraints.size()]();
	cout << "Allocated " << linearConstraints.size() << " for linear constraints" << endl;

	// work on the nonlinear constraints
	nonLinearValues = new uint32_t*[nonLinearConstraints.size()]();
	nonLinearValuesLength = new uint32_t[nonLinearConstraints.size()]();
	nonLinearLast = new uint32_t*[nonLinearConstraints.size()]();
	combinedTimes = new uint16_t[nonLinearConstraints.size()]();
	for (uint32_t i = 0; i < nonLinearConstraints.size(); i++) combinedTimes[i] = 1; 
	nonLinearDistribution = new uint32_t*[nonLinearConstraints.size()]();
	nonLinearDistributionLength = new uint32_t[nonLinearConstraints.size()]();
	nonLinearPossibleKeys = new uint64_t[nonLinearConstraints.size()]();
	cout << "Allocated " << nonLinearConstraints.size() << " for nonlinear constraints" << endl;

	for (uint32_t i = 0; i < linearConstraints.size(); i++)
	{
		uint32_t *values; values = new uint32_t[1<<SIZE]();
		uint32_t valuesLength = 1;
		cout << "Resolving linear constraints[" << i << "]...";
		uint32_t keys = resolveLinear(linearConstraints[i],values);
		cout << "Completed!" << endl;
		linearValues[i] = values;
		linearValuesLength[i] = valuesLength;
		linearKeys.push_back(keys);
	}
	uint32_t i = 0;
	while (i < nonLinearConstraints.size())
	{
		uint32_t *values;
		uint32_t valuesLength;
		cout << "Resolving nonlinear constraints[" << i << "]...";
		nonLinearLast[i] = new uint32_t[1<<SIZE]();
		vector<uint32_t> keys;
		if (resolveNonLinear(nonLinearConstraints[i],values,valuesLength,nonLinearLast[i],keys))
		{
			cout << "Completed!" << endl;
			nonLinearValues[i] = values;
			nonLinearValuesLength[i] = valuesLength;
			nonLinearKeys.push_back(keys);
		}
		else
		{
			int c = nonLinearConstraints[i][nonLinearConstraints[i].size()-1];
			#if SIZE == 4
			int DDT4[16][16] = {0};
			computeDDT(DDT4);
			keyProb -= log2((DDT4[(c>>SIZE)&ANDval][c&ANDval]+0.0)/(1<<SIZE)); // adding back the keyProb since we are ignoring the key constraint
			#elif SIZE == 8
			int DDT8[256][256] = {0};
			computeDDT8(DDT8);
			keyProb -= log2((DDT8[(c>>SIZE)&ANDval][c&ANDval]+0.0)/(1<<SIZE));
			#endif
			nonLinearConstraints.erase(nonLinearConstraints.begin()+i);
			i--;
			delete[] nonLinearLast[i];
		}
		i++;
	}
	if (TK_NUM >= 2)
	{	
		cout << "Combining constraints...";
		combineConstraints();
		cout << "Completed" << endl;
		computeLinearPossibleKeys();
		computeNonLinearDistribution();
		// next step
		double keys = computeKeys();
		cout << "number of keys: " << keys << endl;
		computeProbability();
	}
	else
	{
		cout << "Propagating the constraints...";
		propagateConstraints();
		cout << "Completed" << endl;
		cout << "Combining constraints..." << endl;;
		combineConstraints();
		cout << "Completed" << endl;
		computeLinearPossibleKeys();
		computeNonLinearDistribution();
		double keys = computeKeys();
		cout << "number of keys: " << keys << endl;
		computeProbability();
	}

	return 0;
}