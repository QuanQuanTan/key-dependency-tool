#include <iostream>
#include "deoxys.h"
#include "deoxys_trail.h"
#include <vector>
#include <iomanip>
#include <numeric>
#include <string>
#include "SAT.h"
#include <cryptominisat5/cryptominisat.h>

using namespace std;

void fixedValuesToPath(uint8_t alpha[20][4][4], uint8_t key_diff[3][4][4], uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr)
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
		invMC(tmp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][3][4*r+c] = tmp[r][c];}}
		invSR(tmp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][2][4*r+c] = tmp[r][c];}}
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){tmp[r][c] = alpha[n][r][c];}}
		AddRoundTweakey(tmp,key_diff);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][1][4*r+c] = tmp[r][c];}}
		keyScheduleRoundFunction(key_diff);
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
				if (diffStates[n][1][4*r+c] > 0 || fixedStates[n][1][4*r+c] > 0) {fixedStates[n][2][4*r+c] = 1;}
			}
		}
		// ShiftRows
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				fixedStates[n][3][4*r+c] = fixedStates[n][2][4*r+((4+c+r)%4)];
			}
		}
		// MixColumns
		for (int c = 0; c < 4; c++)
		{
			if (fixedStates[n][3][c] + fixedStates[n][3][4+c] + fixedStates[n][3][8+c] + fixedStates[n][3][12+c] == 4)
			{
				fixedStates[n][4][c] = 1;
				fixedStates[n][4][4+c] = 1;
				fixedStates[n][4][8+c] = 1;
				fixedStates[n][4][12+c] = 1;
			}
		}
		if (n < nr -1)
		{
			for (int r = 0; r < 4; r++)
			{
				for (int c = 0; c < 4; c++)
				{
					fixedStates[n+1][0][4*r+c] = fixedStates[n][4][4*r+c];
					fixedStates[n+1][1][4*r+c] = fixedStates[n][4][4*r+c];
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
				if (fixedStates[n][1][4*r+c] > 0 && diffStates[n][1][4*r+c] > 0)
				{
					// put this into tupsConstraints
					vector<uint32_t> t;
					for (uint32_t i = 0; i < tupLength[4*r+c]; i++) t.push_back(tup[4*r+c][i]);
					t.push_back((n << (4+2*8)) + ((4*r+c) << (2*8)) + (int(diffStates[n][1][4*r+c]) << 8) + (diffStates[n][2][4*r+c]));
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
					(diffStates[n][1][4*r+c] << 8) + (diffStates[n][2][4*r+c]);
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
			res[i^j] += val1[i] * val2[j];
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

// check here onwards
void addRoundConstant(uint32_t tup, uint32_t val[256])
{
	uint32_t constantArray[256] = {0};
	if (((tup >> 16) & 0xf) == 0) constantArray[1] = 1;
	else if (((tup >> 16) & 0xf) == 4) constantArray[2] = 1;
	else if (((tup >> 16) & 0xf) == 8) constantArray[4] = 1;
	else if (((tup >> 16) & 0xf) == 12) constantArray[8] = 1;
	else if (((tup >> 16) & 0xf) == 1) constantArray[getConstants((tup >> 16) & 0xf)] = 1;
	else if (((tup >> 16) & 0xf) == 5) constantArray[getConstants((tup >> 16) & 0xf)] = 1;
	else if (((tup >> 16) & 0xf) == 9) constantArray[getConstants((tup >> 16) & 0xf)] = 1;
	else if (((tup >> 16) & 0xf) == 3) constantArray[getConstants((tup >> 16) & 0xf)] = 1;
	else constantArray[0] = 1;
	XOR(val,constantArray);
	reduceKey(val,1);
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
	addRoundConstant(t,rightPartValues);

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

// Propagating to the last round
void propagateLinear(vector<int> &keys, uint32_t* &values, uint32_t nr, uint32_t &valuesLength, int TK_NUM)
{
	// round number should be the same for all the keys here
	uint32_t n = keys[0] >> (4+2*8);
	uint32_t callsToLFSR = nr - n - 1;
	uint32_t TK2_columns[1<<(2*8)];
	uint32_t TK3_columns[1<<(2*8)];
	for (int col = 0; col < 256; col++)
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
	uint32_t ANDval2 = (1 << (2*8)) - 1;
	// permutate the keys first
	for (uint32_t i = 0; i < keys.size(); i++)
	{
		while (uint32_t(keys[i] >> (4+2*8)) < nr) 
		{
			for (int j = 0; j < 16; j++)
			{
				if (getPermSchedule(j) == ((keys[i] >> (2*8)) & 0xf))
				{
					keys[i] = (((keys[i] >> (4+2*8)) + 1) << (4+2*8)) ^ (j << (2*8)) ^ (keys[i] & ANDval2);
					break;
				}
			}
		}
	}
	uint32_t* newKey;
	valuesLength = TK_NUM;
	newKey = new uint32_t[1<<(8*TK_NUM)]();
	if (TK_NUM == 1) return;
	else if (TK_NUM == 2)
	{
		valuesLength = 2;
		for (int tk1 = 0; tk1 < 256; tk1++)
		{
			for (int diff = 0; diff < 256; diff++)
			{
				if (values[diff] == 0) continue;
				int tk2 = tk1 ^ diff;
				newKey[tk1*256+TK2_columns[tk2]] += values[diff];
			}
		}
	}
	else if (TK_NUM == 3)
	{
		valuesLength = 3;
		for (int tk1 = 0; tk1 < 256; tk1++)
		{
			for (int tk2 = 0; tk2 < 256; tk2++)
			{
				for (int diff = 0; diff < 256; diff++)
				{
					if (values[diff] == 0) continue;
					int tk3 = tk1 ^ tk2 ^ diff;
					newKey[tk1*(1<<(2*8))+TK2_columns[tk2]*256+TK3_columns[tk3]] += values[diff];
				}
			}
		}
	}
	values = newKey;
	return;
}

bool isZero(uint32_t *array, uint32_t length)
{
	for (uint32_t i = 0; i < length; i++)
	{
		if (array[i] != 0) return false;
	}
	return true;
}

int computeProbTrail(uint8_t alpha[20][4][4], uint8_t key_diff[3][4][4], uint32_t nr, int TK_NUM, vector<vector<uint32_t>> linearConstraints, uint32_t **linearValues, uint32_t *linearValuesLength, vector<vector<int>> linearKeys)
{
	int prob = 0;
	uint32_t DDT[256][256] = {0};
	computeDDT(DDT);
	uint32_t ANDval = (1 << 8) - 1;
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
				prob += int(-log2(DDT[before[r][c]][after[r][c]]/256.0));
			}
		}

		keyScheduleRoundFunction(key_diff_temp);
	}
	cout << "Original probability: " << dec << prob << endl;
	for (int i = 0; i < linearConstraints.size(); i++)
	{
		uint8_t inp = ((linearConstraints[i][linearConstraints[i].size()-1] >> 8) & ANDval);
		uint8_t outp = ((linearConstraints[i][linearConstraints[i].size()-1]) & ANDval);
		prob -= int(-log2(DDT[inp][outp]/256.0));
		// cout << log2(DDT[inp][outp]/256.0) << endl;
	}
	cout << "Probability that does not depend on key: " << dec << prob << endl;
	int noKeys = 0;
	uint32_t *valuesByCells; valuesByCells = new uint32_t[(1ULL << (8*TK_NUM))*16]; // 16 cells, with 2**(4 * TK_NUM) keys for each cell 
	for (int i = 0; i < (16*(1 << (8*TK_NUM))); i++) valuesByCells[i] = 1;
	
	// linear
	for (uint32_t i = 0; i < linearKeys.size(); i++)
	{
		for (int j = 0; j < (1 << (8*TK_NUM)); j++) valuesByCells[((linearKeys[i][0] >> (2*8)) & 0xf) * (1 << (8*TK_NUM)) + j] *= linearValues[i][j];
	}

	for (int cell = 0; cell < 16; cell++)
	{
		int s = 0;
		for (int i = 0; i < (1<<(8*TK_NUM)); i++)
		{
			// cout << int(valuesByCells[cell*(1<<(SIZE*TK_NUM))+i]) << " ";
			if (int(valuesByCells[cell*(1<<(8*TK_NUM))+i]) > 0) s += 1;
		}
		if (s == 0) 
		{
			return 0;
		}
		noKeys += int(-log2((s+0.0)/(1<<(8*TK_NUM))));
	}
	delete valuesByCells;
	return noKeys;
}


int main()
{
	rand_generator.seed(time(0));
	uint8_t alpha[20][4][4] = {0};
	uint8_t key_diff[3][4][4] = {0};
	uint8_t diffStates[20][5][16];
	uint8_t fixedStates[20][5][16] = {0};
	uint8_t actualKey[20][4][4] = {0};
	uint32_t nr;
	int TK_NUM;

	// Deoxys trails
	// TK2_8U_2017_693(alpha,key_diff,nr,TK_NUM);
	generateBestTrail(alpha,key_diff,0x10,3,nr,TK_NUM);
	// create2Fail_trail(alpha,key_diff,nr,TK_NUM);
	// justACrazyTrail(alpha,key_diff,nr,TK_NUM);


	fixedValuesToPath(alpha,key_diff,diffStates,fixedStates,nr);
	printState(diffStates,fixedStates,nr);


	vector<vector<uint32_t>> constraints = findConstraints(diffStates,fixedStates,nr);
	
	// just for printing
	cout << "The constraints are:" << endl;
	uint32_t ANDval = (1 << 8) - 1;
	for (uint32_t i = 0; i < constraints.size(); i++)
	{
		for (uint32_t j = 0; j < constraints[i].size(); j++)
		{
			cout << hex << "(" << (constraints[i][j] >> (4+2*8)) << "," << ((constraints[i][j] >> (2*8)) & 0xf) << ","
			<< ((constraints[i][j] >> 8) & ANDval) << "," << (constraints[i][j] & ANDval) << ") ";
		}
		cout << "Linear?: " << isLinear(constraints[i]);
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
	cout << "Allocated " << linearConstraints.size() << " slots for linear constraints" << endl;

	// check linear propagation 
	// check 1. perm schedule
	// check 2. multiple keys how to propagate?
	for (uint32_t i = 0; i < linearConstraints.size(); i++)
	{
		uint32_t *values; values = new uint32_t[(1<<8)]();
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
	// check this
	SATsolverLinear(linearKeys,linearValues,linearValuesLength,TK_NUM);
	int noKeys = computeProbTrail(alpha,key_diff,nr,TK_NUM,linearConstraints,linearValues,linearValuesLength,linearKeys);
	cout << "no. of keys that allow the characteristic: " << noKeys << endl;
	return 0;
}