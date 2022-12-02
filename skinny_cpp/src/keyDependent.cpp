
#include "skinny.h"
#include "trails.h"
#include "keyDependent.h"
using namespace std;

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl
#define DEBUG 1
// general requirements
int TK_NUM;
uint32_t nr;
uint8_t alpha[20][4][4]; // differential characteristic
uint8_t key_diff[4][4][4]; // related key 
uint8_t diffStates[20][5][16] = {0}; // ??
uint8_t fixedStates[20][5][16] = {0}; // ??
uint32_t ANDval = (1<<SIZE)-1;
uint64_t ROUNDval = (1 << 4) - 1;
int inverseMC[16][3] = {{0,10,13},{1,11,14},{2,8,15},{3,9,12},{0,-1,-1},{1,-1,-1},{2,-1,-1},{3,-1,-1},{7,10,-1},{4,11,-1},{5,8,-1},{6,9,-1},{0,10,-1},{1,11,-1},{2,8,-1},{3,9,-1}};
double keyProb = 0;
vector<uint32_t> keyProbCells;


// linear constraints requirements
vector<vector<uint32_t>> linearConstraints;
vector<uint32_t> linearKeys;
vector<uint32_t> propagatedLinearKeys;
vector<int> linearProbReduce; // contain the probability (in log2) for a linear constraint due to the constraint
vector<int> linearProb; // contain the probability (in log2) for a linear constraint according to DDT
uint32_t **linearLast;
uint32_t **linearDistribution;
uint32_t *linearDistributionLength;
uint64_t *linearPossibleKeys;
uint32_t *linearDistributionBase;

// non-linear constraints requirements
vector<vector<uint32_t>> nonLinearConstraints;
vector<vector<uint32_t>> nonLinearKeys;
vector<uint32_t> propagatedNonLinearKeys;
uint32_t **nonLinearValues;
uint32_t *nonLinearValuesLength;
uint32_t **nonLinearLast; // contain the count for the last Sbox transition (getXDDT)
uint64_t *nonLinearPossibleKeys; // keeps track if a key is possible due to the constraint. (For counting possible keys)
uint16_t *combinedTimes; // number of "values" at the end of each tuple
uint32_t **nonLinearDistribution; // 2*jth entry keep track of the count, 2*j+1 keep track of the total. Divide to get prob
uint32_t *nonLinearDistributionLength; 
uint32_t *nonLinearDistributionBase;


// special constraints
vector<vector<uint32_t>> dualKeys; // keys in the form of k0 + k4 = {...}
vector<uint32_t> singleKeys; // keys in the form of k4 = {...}
vector<vector<uint32_t>> sameRoundConstraints;
uint32_t** DValues;
uint32_t** SValues;
uint32_t DValuesLength;
uint32_t SValuesLength;

// experimental constraints
vector<vector<uint32_t>> expNonLinearKeys;
vector<uint32_t> expLinearKeys;

bool isLinear(vector<uint32_t> constraint){
	// a constraint is in the form of {(A_1|B_1|C_1|D_1), ... , (A_n|B_n|C_n|D_n)}
	// C,D (4/8 bits depending on SIZE), records the input,output difference
	// B (4 bits), records the position of the difference
	// A (flexible, up to 32 bits in total), records the round number
	// a constraint is nonlinear if there exist a C_j|D_j = 0.
	for (uint32_t i = 0; i < constraint.size();i++) if ((constraint[i] & ANDval) == 0) return false;
	return true;
}

void fixedValuesToPath(uint8_t characteristic[20][4][4], uint8_t key_diff[4][4][4], uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr){
	// the inputs are characteristic and key_diff
	// diffStates records all the differences progressing from before SB to after MC
	// fixedStates (boolean) records if the value of the particular state is fixed to a certain subset or not (is the output of the Sbox XOR with another element with
	// a uniform distribution?)

	uint8_t key_diff_temp[4][4][4] = {0};
	for (int i = 0; i < 4; i++) { for (int r = 0; r < 4; r++) { for (int c = 0; c < 4; c++) key_diff_temp[i][r][c] = key_diff[i][r][c];}}
	// this loop just putting all of the characteristic into diffStates for visual check.
	for (int n = 0; n < nr; n++)
	{
		uint8_t tmp[4][4];
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				diffStates[n][0][4*r+c] = characteristic[n][r][c];
				diffStates[n][4][4*r+c] = characteristic[n+1][r][c];
				tmp[r][c] = characteristic[n+1][r][c];
			}
		}
		invMC(tmp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][3][4*r+c] = tmp[r][c];}}
		invSR(tmp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][2][4*r+c] = tmp[r][c];}}
		AddRoundTweakey(tmp,key_diff_temp);
		for (int r=0;r<4;r++){for(int c=0;c<4;c++){diffStates[n][1][4*r+c] = tmp[r][c];}}
		#if SIZE == 4
		key_schedule_round_function(key_diff_temp);
		#elif SIZE == 8
		key_schedule_round_function8(key_diff_temp);
		#endif
	}

	// Settling the fixedStates now
	// '1' corresponds to a value fixed to a specfic subset of values
	for (int n = 0; n < nr; n++)
	{
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				// for the state after subst, if sbox is active, OR if the state is already fixed, then the state after subst/ add key is fixed
				if (diffStates[n][0][4*r+c] > 0 || fixedStates[n][0][4*r+c] > 0) {fixedStates[n][1][4*r+c] = 1; fixedStates[n][2][4*r+c] = 1;}
			}
		}

		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				// performing shiftrows for the fixed state
				fixedStates[n][3][4*r+c] = fixedStates[n][2][4*r+((4+c-r)%4)];
			}
		}
		// performing mixcolumn for the fixed state
		for (int c = 0; c < 4; c++)
		{
			if (fixedStates[n][3][c] > 0 && fixedStates[n][3][8+c] > 0 && fixedStates[n][3][12+c] > 0) fixedStates[n][4][c] = 1;
			if (fixedStates[n][3][c] > 0) fixedStates[n][4][4+c] = 1;
			if (fixedStates[n][3][4+c] > 0 && fixedStates[n][3][8+c] > 0) fixedStates[n][4][8+c] = 1;
			if (fixedStates[n][3][c] > 0 && fixedStates[n][3][8+c] > 0) fixedStates[n][4][12+c] = 1;
		}
		// transfer the values from after MC to before SB
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



vector<vector<uint32_t>> findConstraints(uint8_t diffStates[20][5][16],uint8_t fixedStates[20][5][16],int nr){
	// function to obtain constraints out from diffStates and fixedStates
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
	for (int n = 0; n < nr; n++){
		for (int r = 0; r < 4; r++){
			for (int c = 0; c < 4; c++){
				// if a state is fixed and the difference is active, we have a constraint! Put it into tupsConstraints
				if (fixedStates[n][0][4*r+c] > 0 && diffStates[n][0][4*r+c] > 0){
					vector<uint32_t> t;
					// translating the array to vector
					for (uint32_t i = 0; i < tupLength[4*r+c]; i++) t.push_back(tup[4*r+c][i]);
					uint32_t v = (n << (4+2*SIZE)) + ((4*r+c) << (2*SIZE)) + (int(diffStates[n][0][4*r+c]) << SIZE) + (diffStates[n][1][4*r+c]);
					t.push_back(v);
					tupsConstraints.push_back(t);
					tupLength[4*r+c] = 0; // erase the data since it has been taken into account

					auto it = find(keyProbCells.begin(),keyProbCells.end(),v);
					if (it == keyProbCells.end()) 
					{
						keyProbCells.push_back(v);
						#if SIZE == 4
							keyProb += -log2(DDT4[diffStates[n][0][4*r+c]][diffStates[n][1][4*r+c]]/16.0);
						#elif SIZE == 8
							keyProb += -log2(DDT8[diffStates[n][0][4*r+c]][diffStates[n][1][4*r+c]]/256.0);
						#endif
					}
				}
			}
		}
		// add the difference into tup
		for (int r = 0; r < 4; r++){
			for (int c = 0; c < 4; c++){
				tup[4*r+c][tupLength[4*r+c]] = (n << (4+2*SIZE)) + ((4*r+c) << (2*SIZE)) + 
					(diffStates[n][0][4*r+c] << SIZE) + (diffStates[n][1][4*r+c]);
				tupLength[4*r+c]++;
			}
		}

		// applying SR MC to the key
		// row 0
		for (uint32_t i = 0; i < 4; i++){
			for (uint32_t j = 0; j < tupLength[i]; j++){
				tupTemp[i][tupTempLength[i]++] = tup[i][j];
			}	
			for (uint32_t j = 0; j < tupLength[8+((2+i)%4)]; j++){
				tupTemp[i][tupTempLength[i]++] = tup[8+((2+i)%4)][j];
			}	
			for (uint32_t j = 0; j < tupLength[12+((1+i)%4)]; j++){
				tupTemp[i][tupTempLength[i]++] = tup[12+((1+i)%4)][j];
			}	
		}
		// row 1
		for (uint32_t i = 0; i < 4; i++){
			for (uint32_t j = 0; j < tupLength[i]; j++){
				tupTemp[4+i][tupTempLength[4+i]++] = tup[i][j];
			}	
		}
		// row 2
		for (uint32_t i = 0; i < 4; i++){
			for (uint32_t j = 0; j < tupLength[4+((3+i)%4)]; j++){
				tupTemp[8+i][tupTempLength[8+i]++] = tup[4+((3+i)%4)][j];
			}	
			for (uint32_t j = 0; j < tupLength[8+((2+i)%4)]; j++){
				tupTemp[8+i][tupTempLength[8+i]++] = tup[8+((2+i)%4)][j];
			}	
		}
		// row 3
		for (uint32_t i = 0; i < 4; i++){
			for (uint32_t j = 0; j < tupLength[i]; j++){
				tupTemp[12+i][tupTempLength[12+i]++] = tup[i][j];
			}	
			for (uint32_t j = 0; j < tupLength[8+((2+i)%4)]; j++){
				tupTemp[12+i][tupTempLength[12+i]++] = tup[8+((2+i)%4)][j];
			}	
		}

		// put it back to tup from tupTemp for those that are 1
		// erase those with 0
		for (int r = 0; r < 4; r++){
			for (int c = 0; c < 4; c++){
				if (fixedStates[n][4][4*r+c] == 0) tupLength[4*r+c] = 0;
				else {
					for (uint32_t i = 0; i < tupTempLength[4*r+c]; i++) tup[4*r+c][i] = tupTemp[4*r+c][i];
					tupLength[4*r+c] = tupTempLength[4*r+c];
				}
				tupTempLength[4*r+c] = 0;
			}
		}
	}
	return tupsConstraints;
}


void reduceKey(uint32_t* &val, uint32_t valLength) {
	// this reduces the exact values, but according to proportions
	// this helps to make the code run more efficiently
	int counter = 0;
	int a[2]; a[0] = 1; a[1] = 1;
	for (int i = 0; i < (1<<(valLength*SIZE)); i++){
		// we ignre values that are 0
		if (val[i] > 0) {
			if (counter == 0) a[counter++] = val[i];
			else {	
				a[1] = val[i];
				a[0] = gcd(a[0],a[1]);
				if (a[0] == 1) return; // terminate if it's 1
			}
		}
	}
 	for (int i = 0; i < (1<<(valLength*SIZE)); i++) val[i] /= a[0];
}

void getXDDT(uint32_t input_diff, uint32_t output_diff, uint32_t output[1<<SIZE]){
	// this gets the X_{DDT} (input values that is valid for the differential transition)
	#if SIZE == 4
	for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++){
		for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++){
			if (((v1^v2) == input_diff) && ((getSbox(v1)^getSbox(v2)) == output_diff)) output[v1]++;
		}
	}
	#elif SIZE == 8
	for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++){
		for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++){
			if (((v1^v2) == input_diff) && ((getSbox8(v1)^getSbox8(v2)) == output_diff)) output[v1]++;
		}
	}
	#else 
		#error Unsupported choice setting
	#endif
}

void getYDDT(uint32_t input_diff, uint32_t output_diff, uint32_t output[1<<SIZE]){
	// this gets the Y_{DDT} (output values that is valid for the differential transition)
	#if SIZE == 4
	for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++){
		for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++){
			if (((v1^v2) == input_diff) && ((getSbox(v1)^getSbox(v2)) == output_diff)) output[getSbox(v1)]++;
		}
	}
	#elif SIZE == 8
	for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++){
		for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++){
			if (((v1^v2) == input_diff) && ((getSbox8(v1)^getSbox8(v2)) == output_diff)) output[getSbox8(v1)]++;
		}
	}
	#else 
		#error Unsupported choice setting
	#endif
}

void XOR(uint32_t val1[(1<<SIZE)],uint32_t val2[(1<<SIZE)]){
	// input --> val1 and val2 with the INDICES representing the value we want to XOR
	// output --> val1
	// if val2 is [0,0,...,0], return val1 unchanged
	// if val1 is [0,0,...,0], replace val1 values with val2
	uint32_t sum1 = 0, sum2 = 0;
	for (uint32_t i = 0; i < (1<<SIZE); i++){
		sum1 += val1[i];
		sum2 += val2[i];
	}
	if (sum2 == 0) return;
	if (sum1 == 0){
		for (uint32_t i = 0; i < (1<<SIZE); i++) val1[i] = val2[i];
		return;
	}
	uint32_t res[(1<<SIZE)] = {0};
	for (int i = 0; i < (1<<SIZE); i++){
		for (int j = 0; j < (1<<SIZE); j++) res[i^j] += val1[i] * val2[j];
	}
	for (int i = 0; i < (1<<SIZE); i++) val1[i] = res[i];
}

void addRoundConstant(uint32_t tup, uint32_t val[(1<<SIZE)]){
	// this function add the round constant in the context of tuples (1,2,4,4) etc
	uint32_t constantArray[(1<<SIZE)] = {0};
	if (((tup >> (SIZE*2)) & 0xf) == 0) constantArray[getConstants(tup >> (4+SIZE*2)) & 0xf] = 1;
	else if (((tup >> (SIZE*2)) & 0xf) == 4) constantArray[(getConstants(tup >> (4+SIZE*2)) >> 4) & 0b11] = 1;
	else if (((tup >> (SIZE*2)) & 0xf) == 8) constantArray[2] = 1;
	XOR(val,constantArray);
	reduceKey(val,1);
}

uint32_t resolveLinear(vector<uint32_t> constraint, uint32_t *distribution, uint32_t &distributionLength, uint32_t &distributionBase){
	// input: constraint
	// output: key: the tuple with the key involved in the constraint (round number and position recorded)
	uint32_t last[1<<SIZE] = {0};
	uint32_t values[1<<SIZE] = {0};

	uint32_t key = 0;
	int t = constraint[constraint.size()-1]; // the last tuple in constraint is the RHS 
	getXDDT((t >> SIZE) & ANDval,t & ANDval,last); // get the X_{DDT} for the RHS

	for (uint32_t i = 0; i < constraint.size()-1; i++) {
		if (((constraint[i] >> (2*SIZE)) & 0xf) < 8) key = constraint[i]; // in Skinny, exactly 1 key is involved in each linear constraint
		addRoundConstant(constraint[i],values); 
		uint32_t yddt[1<<SIZE] = {0};
		getYDDT((constraint[i] >> SIZE) & ANDval,constraint[i] & ANDval,yddt);
		XOR(values,yddt); // XOR all the values together since they are connected via linear layers only
	}	

	for (int k = 0; k < (1<<SIZE); k++){
		for (int i = 0; i < (1<<SIZE); i++){
			for (int j = 0; j < (1<<SIZE); j++) {
				if ((i ^ k) == j){
					distribution[k] += values[i] * last[j];
					if (k==0) distributionBase += values[i];
				}
			}
		}
	}
	// reduce the value by the distribution base
	uint32_t a = distributionBase;
	for (int k = 1; k < (1<<SIZE); k++)
	{
		if (distribution[k] != 0) a = gcd(a,distribution[k]);
	}
	for (int k = 0; k < (1<<SIZE); k++) distribution[k] /= a;
	distributionBase /= a;

	distributionLength = 1;
	return key;
}

vector<uint32_t> resolveLeft(uint32_t zeroPoint, vector<uint32_t> constraint, uint32_t values[(1<<SIZE)], vector<uint32_t> &zeroDependencies, bool last = false){
	// part of resolveNonLinear to deal with the left side of an inactive Sbox
	// this part is very similar to how we do with resolveLinear, just without involving X_{DDT}
	// input: zeroPoint, the tuple showing where the inactive Sbox lies. What we want are tuples that lie on the left of this
	// input: constraint, all the tuples involved in this constraint
	// output: all the keys involved here as a vector
	// values: recording down the XOR of all the values here
	// zeroDependencies: record down any inactive Sboxes on the left of what we detected
	vector<uint32_t> keys;
	uint32_t round = (zeroPoint >> (4+2*SIZE)); 
	uint32_t pos = (zeroPoint >> (2*SIZE)) & 0xf;
	int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
	for (uint32_t i = 0; i < constraint.size(); i++)
	{
		uint32_t constraintRound = constraint[i] >> (4+2*SIZE);
		int constraintPos = (constraint[i] >> (2*SIZE)) & ROUNDval;
		// check if the value is involved in this zeroPoint
		// if it's not on the left of the zeroPoint, ignore
		// if it's on the left, but not connected to zeroPoint, ignore
		if (constraintRound != round-1) continue;
		if (constraintPos != inversePos[0] && constraintPos != inversePos[1] && constraintPos != inversePos[2]) continue;
		if (constraintPos < 8) keys.push_back(constraint[i]); // record the key
		if ((constraint[i] & ANDval) == 0) zeroDependencies.push_back(constraint[i]); // if the constraint is to the left of zeroPoint and it's another inactiveSbox, save it
		else if (last == false) // record down the XOR values of all these constraints
		{
			uint32_t yddt[1<<SIZE] = {0};
			getYDDT((constraint[i] >> SIZE) & ANDval,constraint[i] & ANDval,yddt);
			addRoundConstant(constraint[i],values);
			XOR(values,yddt);
			reduceKey(values,1);
		}
	}
	return keys;
}

bool AddKeyNonLinear(uint32_t* &value,uint32_t &valueLength){
	// we have to extend the dimension to accomodate the values for each key
	// thus, valueLength will be increased by 1
	// input: value
	// output: value, but a new dimension has been added such that the keys are taking into consideration
	uint32_t keySize = (1ULL << ((valueLength-1) * SIZE)); // the last "SIZE" is dedicated for values, not key
	uint64_t destinationSize = (1ULL << ((valueLength+1)*SIZE));
	if ((valueLength+1)*SIZE >= 36) {
		cout << "this constraint exceed 2^36 space. Skipping this constraint." << endl;
		return false;
	}
	uint32_t* destination; destination = new uint32_t [destinationSize]();
	#if SIZE == 4
		for (uint64_t i = 0; i < keySize; i++){
			for (int k = 0; k < 16; k++){
				for (int v = 0; v < 16; v++){
					destination[i*(16*16)+k*16+v] += value[i * 16 + (v^k)];
				}
			}
		}
	#elif SIZE == 8
		for (uint64_t i = 0; i < keySize; i++){
			for (int k = 0; k < 256; k++){
				for (int v = 0; v < 256; v++){
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

void SubstNonLinear(uint32_t* &value,uint32_t &valueLength){
	// This function applies subst to the values. (Note only the last SIZE bits are for values, the rest are for keys)
	uint32_t destinationSize = (1ULL << ((valueLength)*SIZE));
	uint32_t* destination; destination = new uint32_t [destinationSize]();
	#if SIZE == 4
		for (uint64_t i = 0; i < (1ULL<<((valueLength-1)*SIZE)); i++){
			for (int v = 0; v < 16; v++){
				destination[i*16+getSbox(v)] += value[i*16 + v];
			}
		}

	#elif SIZE == 8
		for (uint64_t i = 0; i < (1ULL<<((valueLength-1)*SIZE)); i++){
			for (int v = 0; v < 256; v++){
				destination[i*256+getSbox8(v)] += value[i*256 + v];
			}
		}
	#endif
	delete[] value;
	value = destination;
}

void resolveZeroDependencies(uint32_t* valuesToAdd, uint32_t valuesToAddLength, uint32_t* &destination, uint32_t &destinationLength){
	//  this function combines the keys from valuesToAdd and destination
	// for example valuesToAdd = (k2,v1), destination = (k0,k1,v0)
	// the output is (k0,k1,k2,v0^v1)
	uint32_t valuesKeySize = 1 << ((valuesToAddLength-1) * SIZE);
	uint32_t destinationKeySize = 1 << ((destinationLength-1) * SIZE);
	uint32_t *dest;
 	dest = new uint32_t [(1ULL<<((valuesToAddLength+destinationLength-1)*SIZE))]();
	for (uint32_t i = 0; i < valuesKeySize; i++){
		for (uint32_t j = 0; j < destinationKeySize; j++){	
			XOR(&dest[i*(destinationKeySize*(1 << SIZE))+j*(1 << SIZE)],&valuesToAdd[i*(1 << SIZE)]);
			XOR(&dest[i*(destinationKeySize*(1 << SIZE))+j*(1 << SIZE)],&destination[j*(1 << SIZE)]);
		}
	}
	reduceKey(dest,valuesToAddLength + destinationLength - 1);
	delete[] destination;
	destination = dest;
	destinationLength = valuesToAddLength + destinationLength - 1;
}


void AddFinalgetYDDT(vector<uint32_t> constraint, uint32_t* &val, uint32_t valLength){
	// This is to deal with last YDDT before the RHS (if it exists)
	uint32_t last = constraint[constraint.size()-1];
	uint32_t *valLast; valLast = new uint32_t[1<<SIZE]();
	uint32_t round = (last >> (4+2*SIZE));
	uint32_t pos = (last >> (2*SIZE)) & 0xf;
	int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
	for (uint32_t i = 0; i < constraint.size(); i++){
		uint32_t constraintRound = constraint[i] >> (4+2*SIZE);
		int constraintPos = (constraint[i] >> (2*SIZE)) & ROUNDval;
		// check if the value is involved in this zeroPoint
		if (constraintRound != round-1) continue;
		if (constraintPos != inversePos[0] && constraintPos != inversePos[1] && constraintPos != inversePos[2]) continue;
		if ((constraint[i] & ANDval) == 0) continue;
		uint32_t yddt[1<<SIZE] = {0};
		uint32_t tmp[1<<SIZE] = {0};
		getYDDT((constraint[i] >> SIZE) & ANDval,constraint[i] & ANDval,yddt);
		addRoundConstant(constraint[i],valLast);
		XOR(valLast,yddt);
		XOR(valLast,tmp);
		reduceKey(valLast,1);
	}
	uint32_t *tmp; tmp = new uint32_t[1ULL<<(SIZE*valLength)]();
	// uint32_t s = 0;
	// for (uint64_t i = 0; i < (1ULL<<(valLength*SIZE)); i++) s += val[i];
	for (int j = 0; j < (1<<SIZE); j++){
		if (valLast[j] == 0) continue;
		for (uint64_t i = 0; i < (1ULL<<(valLength*SIZE)); i++)
		{
			tmp[uint64_t(i^j)] += val[i];
		}
	}
	val = tmp;
}

bool resolveNonLinear(vector<uint32_t> constraint, uint32_t* &values, uint32_t &valuesLength, uint32_t lastValues[1<<SIZE], vector<uint32_t> &output){
	// this is the main function for solving the nonlinear functions
	// first, we need to locate the zeros (which determine how many rounds/ connection points)
	sort(constraint.begin(),constraint.end());
	vector<uint32_t> zeroPoints; 
	for (uint32_t i = 0; i < constraint.size(); i++) {
		if ((constraint[i] & ANDval) == 0) zeroPoints.push_back(constraint[i]);
	}
	uint32_t* zeroValues[zeroPoints.size()] = {0};
	uint32_t zeroValuesLength[zeroPoints.size()] = {0}; // this will tell me how many keys are already involved


	vector<vector<uint32_t>> zeroKeys;
	for (uint32_t i = 0; i < zeroPoints.size(); i++) {
		// the zeroPoint has been sorted from left to right. 
		vector<uint32_t> zeroDependencies; // for each zeroPoint, we want to find out what are the tuples connect to the left of this
		uint32_t* val; val = new uint32_t[(1<<SIZE)]();
		uint32_t valLength = 1; // by default, the value is 1. Unless there are more than 1 zeroDependencies adding to the same zeroPoint
		vector<uint32_t> keys = resolveLeft(zeroPoints[i],constraint, val, zeroDependencies); // zeroDependencies are the ones that require us to add the dimension of key
		// settle if there is any zeroDependencies
		for (uint32_t j = 0; j < zeroDependencies.size(); j++){
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
		if (!AddKeyNonLinear(val,valLength)) {return false;} // if the space required is too huge, this will terminate nonLinear search further
		SubstNonLinear(val,valLength); // Substitution layer is to be added after every zeroPoint
		zeroKeys.push_back(keys);
		zeroValues[i] = val;
		zeroValuesLength[i] = valLength;
	}
	// Settle the last tuple!
	uint32_t* val; val = new uint32_t[1<<SIZE]();
	uint32_t valLength = 1;
	vector<uint32_t> zeroDependencies;
	vector<uint32_t> keys = resolveLeft(constraint[constraint.size()-1], constraint, lastValues, zeroDependencies, true); // last argument changed it from resolveLeft to resolveLast
	// putting lastValues using getXDDT
	uint32_t xddt[1<<SIZE] = {0};
	getXDDT((constraint[constraint.size()-1] >> SIZE) & ANDval,constraint[constraint.size()-1] & ANDval,xddt);
	XOR(lastValues,xddt);

	if (zeroDependencies.size() == 0){
		// sanity check
		cout << "zeroDependencies cannot be zero at this point" << endl;
		exit(0);
	}
	for (uint32_t j = 0; j < zeroDependencies.size(); j++){
		uint32_t index;
		for (uint32_t k = 0; k < zeroPoints.size(); k++) { if (zeroPoints[k] == zeroDependencies[j]) index = k;}
		// at this point, we should capture the info: val gives the keys that we want, zeroValues has the entire distribution.
		resolveZeroDependencies(zeroValues[index],zeroValuesLength[index],val,valLength);
		keys.insert(keys.begin(),zeroKeys[index].begin(),zeroKeys[index].end());
	}
	if (!AddKeyNonLinear(val,valLength)) {return false;}
	AddFinalgetYDDT(constraint,val,valLength);
	values = val;
	valuesLength = valLength;
	output = keys;
	return true;
}

void computeNonLinearDistribution(){
	for (uint32_t i = 0; i < nonLinearKeys.size(); i++){
		nonLinearDistribution[i] = new uint32_t[1ULL<<((nonLinearValuesLength[i]-1)*SIZE)]();
		nonLinearDistributionLength[i] = nonLinearValuesLength[i]-1;
		// printout(nonLinearDistributionLength[i]);
		// printout(nonLinearKeys[i].size());
		for (uint64_t j = 0; j < (1ULL<<((nonLinearValuesLength[i]-1)*SIZE)); j++){
			for (uint32_t k = 0; k < (1<<SIZE); k++)
			{
				if (nonLinearLast[i][k] > 0) nonLinearDistribution[i][j] += nonLinearValues[i][(j<<SIZE)+k]; // the count
				if (j == 0) nonLinearDistributionBase[i] += nonLinearValues[i][(j<<SIZE)+k]; // the total
			}
		}
	}
}


void computeNonLinearPossibleKeys(){
	for (uint32_t i = 0; i < nonLinearKeys.size(); i++){
		for (uint64_t j = 0; j < (1ULL<<(nonLinearDistributionLength[i]*SIZE)); j++){
			if (nonLinearDistribution[i][j] > 0) nonLinearPossibleKeys[i]++;
		}
		if (nonLinearPossibleKeys[i] == 0){
			cout << "nonLinearConstraints[" << i << "] is not satisified!" << endl;
			cout << "there are no keys satisfying the nonlinear constraints!" << endl;
			exit(0);
		}
	}
	
}


void computeLinearPossibleKeys(){
	for (uint32_t i = 0; i < linearKeys.size(); i++){
		for (uint32_t j = 0; j < (1<<SIZE); j++){
			if (linearDistribution[i][j] > 0) linearPossibleKeys[i]++;
		}
		if (linearPossibleKeys[i] == 0){
			cout << "there are no keys satisfying the linear constraints!" << endl;
			exit(0);
		}
	}
}

double getOriginalProbability(){
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

	for (uint32_t n = 0; n < nr; n++) {
		// this loop computes the original prob
		for (uint8_t r = 0; r < 4; r++){
			for (uint8_t c = 0; c < 4; c++){
				before[r][c] = alpha[n][r][c];
				after[r][c] = alpha[n+1][r][c];
			}
		}
		invMC(after); invSR(after); AddRoundTweakey(after,key_diff_temp);
		for (uint8_t r = 0; r < 4; r++){
			for (uint8_t c = 0; c < 4; c++){
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

void remove(vector<int> &v) {
	// only used for getSubgraph
    auto itr = v.begin();
    unordered_set<int> s;
 
    for (auto curr = v.begin(); curr != v.end(); ++curr){
        if (s.insert(*curr).second) {
            *itr++ = *curr;
        }
    }
    v.erase(itr, v.end());
}

vector<vector<int>> getSubgraph(vector<vector<int>> G, bool Determined[16], vector<int> &keyPositions){
	// find the start point
	uint8_t start = 0;
	while (Determined[start] != 0 && start < 16) start++;
	vector<int> allConstraints = G[start];
	keyPositions.push_back(start);
	Determined[start] = 1;
    vector<vector<int>> constraintIndex; constraintIndex.push_back(allConstraints);
    
	uint8_t i = 0;
	while (i < allConstraints.size()){
		for (uint8_t j = start+1; j < 16; j++){
			if (Determined[j] == 1) continue; // visited
			auto it = find(G[j].begin(),G[j].end(),allConstraints[i]);
			if (it == G[j].end()) continue;
			allConstraints.insert(allConstraints.end(),G[j].begin(),G[j].end());
			constraintIndex.push_back(G[j]);
			keyPositions.push_back(j);
			Determined[j] = 1;
			remove(allConstraints); // remove duplicates
		}
		i++;
	}
	return constraintIndex;
}

double subGraphKey(vector<int> keyPositions, vector<vector<int>> constraintsIndex){
	// first, find how many 
	double currentKeys = 0;
	// for (int i = 0; i < keyPositions.size(); i++)
	// {
	// 	cout << keyPositions[i] << " ";
	// }
	// cout << endl;
	vector<int> accountedFor;
	for (uint32_t i = 0; i < keyPositions.size(); i++){
		if (constraintsIndex[i].size() > TK_NUM){
			cout << "Intersection exceeded more than TK size. Please go for experiment" << endl;
			exit(0);
			// if SAT works, then there exists exactly one key
			currentKeys += 0;
		}
		if (constraintsIndex[i].size() < TK_NUM) currentKeys += (TK_NUM-constraintsIndex[i].size()) * SIZE; // we have to expand if there is no multiplication. implicit
		for (uint32_t j = 0; j < constraintsIndex[i].size(); j++){
			auto it = find(accountedFor.begin(),accountedFor.end(),constraintsIndex[i][j]);
			if (it != accountedFor.end()) continue;
			else if (constraintsIndex[i][j] >= 0) currentKeys += log2(linearPossibleKeys[constraintsIndex[i][j]]); // linear keys
			else currentKeys += log2(nonLinearPossibleKeys[-constraintsIndex[i][j]-1]);  // nonlinearkeys
		}
		// currentKeys += (TK_NUM-1) * SIZE;
		accountedFor.insert(accountedFor.end(),constraintsIndex[i].begin(), constraintsIndex[i].end());
		remove(accountedFor);
	}
	currentKeys -= (TK_NUM-1) * SIZE * keyPositions.size(); // note that we minus the TK, so we can compare the size in terms of the "XORed" keys
	// printout(currentKeys);
	return currentKeys;
}

double computeKeys(){
	bool Determined[16] = {0}; // this determines if the keys have been accounted for (use permSchedule to get to round nr)
	double tmpKeys = 0;
	uint32_t n; uint32_t pos;

	// we first compute the number of keys
	vector<vector<int>> propagatedCell; // this keep track of what cell is affected by which constraint
	for (uint32_t i = 0; i < 16; i++){
		vector<int> tmp;
		propagatedCell.push_back(tmp);
	}
	#if DEBUG == 1
	for (uint32_t i = 0; i < linearKeys.size(); i++){
		n = linearKeys[i] >> (2*SIZE+4);
		pos = (linearKeys[i] >> (2*SIZE)) & 0xf;
		for (uint32_t k = n; k < 4; k++) pos = getInvPermSchedule(pos); // propagate to round nr
		propagatedCell[pos].push_back(i);
	}
	for (uint32_t i = 0; i < nonLinearKeys.size(); i++){
		for (uint32_t j = 0; j < nonLinearKeys[i].size(); j++){
			n = nonLinearKeys[i][j] >> (2*SIZE+4);
			pos = (nonLinearKeys[i][j] >> (2*SIZE)) & 0xf;
			for (uint32_t k = n; k < 4; k++) pos = getInvPermSchedule(pos);  
			propagatedCell[pos].push_back(-i-1);
		}
	}
	#else
	for (uint32_t i = 0; i < linearKeys.size(); i++){
		n = linearKeys[i] >> (2*SIZE+4);
		pos = (linearKeys[i] >> (2*SIZE)) & 0xf;
		for (uint32_t k = n; k < nr; k++) pos = getInvPermSchedule(pos); // propagate to round nr
		propagatedCell[pos].push_back(i);
	}
	for (uint32_t i = 0; i < nonLinearKeys.size(); i++){
		for (uint32_t j = 0; j < nonLinearKeys[i].size(); j++){
			n = nonLinearKeys[i][j] >> (2*SIZE+4);
			pos = (nonLinearKeys[i][j] >> (2*SIZE)) & 0xf;
			for (uint32_t k = n; k < nr; k++) pos = getInvPermSchedule(pos);  
			propagatedCell[pos].push_back(-i-1);
		}
	}
	#endif
	// for (int i = 0; i < 16; i++)
	// {
	// 	printout(i);
	// 	for (int j = 0; j < propagatedCell[i].size(); j++)
	// 	{
	// 		cout << propagatedCell[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }
	// propagatedCell is filled at this point. propgatedCell[4] = {1,0,-1} --> key 4 is affected by linearConstraint[0,1] and nonLinearConstraints[0]
	// Finding a subgraph using propagatedCell
	// 2 vectors in total. The first one will record the set of constraints index. nonlinear constraints start at negatives. -1 --> 0, -2 --> 1, ...
	// The second vector will record the set of key cells involved
	bool flag = true;
	while (flag){
		vector<int> keyPositions;
		vector<vector<int>> constraintsIndex = getSubgraph(propagatedCell, Determined, keyPositions); // constraintsindex break them down into each individual key
		if ((keyPositions.size() == 1) && (constraintsIndex[0].size() == 0)) 
		{
			// cout << "no restrictions: ";
			// cout << keyPositions[0] << endl;
			tmpKeys += SIZE; // when there is no constraints, all keys are possible
		}
		else tmpKeys += subGraphKey(keyPositions,constraintsIndex); // otherwise, count the keys
		// to check if we continue
		flag = false;
		for (int i = 0; i < 16; i++) {if (Determined[i] == 0) flag = true;} // continue until all the key positions are exhausted
	}
	
	return tmpKeys;
}

void computeProbability() // TBC
{
	double prob = getOriginalProbability();
	cout << "original prob: " << prob << endl;
	// #if SIZE == 4
	// int DDT4[16][16] = {0};
	// computeDDT(DDT4);	

	// #elif SIZE == 8
	// int DDT8[256][256] = {0};
	// computeDDT8(DDT8);
	// #endif

	prob -= keyProb; // independent of the ones involved in key constraints (nonlinear constraints only)

	// printout(prob);
	// printout(keyProb);
	// for (uint32_t i = 0; i < linearKeys.size(); i++){ // obselete
	// 	prob += linearProb[i] - linearProbReduce[i];
	// }
	uint64_t **probArray_Linear; probArray_Linear = new uint64_t*[linearKeys.size()]();
	uint64_t *probLength_Linear; probLength_Linear = new uint64_t[linearKeys.size()]();
	uint64_t *probBase_Linear; probBase_Linear = new uint64_t[linearKeys.size()]();

	for (uint32_t i = 0; i < linearKeys.size(); i++){
		uint64_t max = 0;
		uint64_t base = 0;
		vector<uint64_t> nonZerosProb;
		for (uint64_t j = 0; j < (1ULL<<SIZE); j++)	{
			if (linearDistribution[i][j] > 0) {
				auto it = find(nonZerosProb.begin(),nonZerosProb.end(),linearDistribution[i][j]);
				if (it == nonZerosProb.end()) nonZerosProb.push_back(linearDistribution[i][j]); // find out all the diferent non-zero probabilities
				if (linearDistribution[i][j] > max) max = linearDistribution[i][j]; // find out the max of the these non-zero probabilities
			}
		}
		base = linearDistributionBase[i];
		if (nonZerosProb.size() == 0) {
			cout << "probability: " << 0 << "as linearConstraint[" << i << "] is not satisfied!" << endl;
			return;
		}
		probLength_Linear[i] = max+1; 
		probArray_Linear[i] = new uint64_t[probLength_Linear[i]]();
		for (uint64_t j = 0; j < (1ULL<<SIZE); j++)	{
			if (linearDistribution[i][j] == 0) continue;
			probArray_Linear[i][linearDistribution[i][j]]++; // we have a probArray which gives us the distribution 
		}
		probBase_Linear[i] = base; 
		vector<uint64_t> nonZerosNumbers;
		for (uint64_t j = 0; j < probLength_Linear[i]; j++){
			if (probArray_Linear[i][j] > 0) nonZerosNumbers.push_back(probArray_Linear[i][j]); // only get the nonzero numbers in nonZerosNumbers
		}
		uint64_t divisor = nonZerosNumbers[0];
		for (uint32_t j = 1; j < nonZerosNumbers.size(); j++) divisor = gcd(divisor, nonZerosNumbers[j]);
		for (uint64_t j = 0; j < probLength_Linear[i]; j++) probArray_Linear[i][j] /= divisor;
		
	}
	uint64_t **probArray_NonLinear; probArray_NonLinear = new uint64_t*[nonLinearKeys.size()]();
	uint64_t *probLength_NonLinear; probLength_NonLinear = new uint64_t[nonLinearKeys.size()]();
	uint64_t *probBase_NonLinear; probBase_NonLinear = new uint64_t[nonLinearKeys.size()]();
	for (uint32_t i = 0; i < nonLinearKeys.size(); i++){
		uint64_t max = 0;
		uint64_t base = 1;
		vector<uint64_t> nonZerosProb;
		for (uint64_t j = 0; j < (1ULL<<(nonLinearDistributionLength[i]*SIZE)); j++)	{
			if (nonLinearDistribution[i][j] > 0) {
				auto it = find(nonZerosProb.begin(),nonZerosProb.end(),nonLinearDistribution[i][j]);
				if (it == nonZerosProb.end()) nonZerosProb.push_back(nonLinearDistribution[i][j]); // find out all the diferent non-zero probabilities
				if (nonLinearDistribution[i][j] > max) max = nonLinearDistribution[i][j]; // find out the max of the these non-zero probabilities
			}
		}
		base = nonLinearDistributionBase[i]; // this is the base of the nonLinearDistribution
		if (nonZerosProb.size() == 0) {
			cout << "probability: " << 0 << "as nonLinearConstraint[" << i << "] is not satisfied!" << endl;
			return;
		}
		uint64_t divisor = base;
		if (divisor == 0){
			cout << "divisor is zero!" << endl;
			exit(0);
		}
		for (uint32_t j = 0; j < nonZerosProb.size(); j++) divisor = gcd(divisor, nonZerosProb[j]); // dividing them to get the smallest divisor
		probLength_NonLinear[i] = max/divisor+1; 
		probArray_NonLinear[i] = new uint64_t[probLength_NonLinear[i]]();
		base = base/divisor; // re-adjusting the base according to the gcd
		for (uint64_t j = 0; j < (1ULL<<(nonLinearDistributionLength[i]*SIZE)); j++)	{
			if (nonLinearDistribution[i][j] == 0) continue;
			probArray_NonLinear[i][nonLinearDistribution[i][j]/divisor]++; // we have a probArray which gives us the distribution 
		}
		probBase_NonLinear[i] = base;
		vector<uint64_t> nonZerosNumbers;
		for (uint64_t j = 0; j < probLength_NonLinear[i]; j++){
			if (probArray_NonLinear[i][j] > 0) nonZerosNumbers.push_back(probArray_NonLinear[i][j]); // only get the nonzero numbers in nonZerosNumbers
		}
		divisor = nonZerosNumbers[0];
		for (uint32_t j = 1; j < nonZerosNumbers.size(); j++) divisor = gcd(divisor, nonZerosNumbers[j]);
		for (uint64_t j = 0; j < probLength_NonLinear[i]; j++) probArray_NonLinear[i][j] /= divisor;
	}
	// combine all the probabilities
	uint64_t **probArray; probArray = new uint64_t*[linearKeys.size() + nonLinearKeys.size()];
	uint64_t *probLength; probLength = new uint64_t[linearKeys.size() + nonLinearKeys.size()]();
	uint64_t *probBase; probBase = new uint64_t[linearKeys.size() + nonLinearKeys.size()]();

	for (int i = 0; i < linearKeys.size(); i++) probArray[i] = probArray_Linear[i];
	for (int i = 0; i < nonLinearKeys.size(); i++) probArray[linearKeys.size() + i] = probArray_NonLinear[i];

	for (int i = 0; i < linearKeys.size(); i++) probLength[i] = probLength_Linear[i];
	for (int i = 0; i < nonLinearKeys.size(); i++) probLength[linearKeys.size() + i] = probLength_NonLinear[i];

	for (int i = 0; i < linearKeys.size(); i++) probBase[i] = probBase_Linear[i];
	for (int i = 0; i < nonLinearKeys.size(); i++) probBase[linearKeys.size() + i] = probBase_NonLinear[i];
	// combine the nonlinear
	uint64_t combinedProbBase = probBase[0];
	uint64_t combinedMax = probLength[0];
	for (uint32_t i = 1; i < linearKeys.size() + nonLinearKeys.size(); i++) {
		combinedProbBase *= probBase[i]; 
		combinedMax *= probLength[i];
	}
	// printout(combinedMax);
	// printout(combinedProbBase);
	uint64_t *combinedProbArray; combinedProbArray = new uint64_t[combinedMax]();
	uint64_t *indices = new uint64_t[linearKeys.size() + nonLinearKeys.size()]();

	// for (int i = linearKeys.size(); i < linearKeys.size()+nonLinearKeys.size(); i++)
	// {
	// 	printout(i);
	// 	for (int j = 0; j < probLength[i]; j++)
	// 	{
	// 		cout << probArray[i][j] << " " << probBase[i] << endl;
	// 	}
	// 	cout << endl;
	// }

	// for (int i = 0; i < linearKeys.size()+nonLinearKeys.size(); i++)
	// {
	// 	printout(i);
	// 	for (int j = 0; j < probLength[i]; j++)
	// 	{
	// 		cout << probArray[i][j] << " " << probBase[i] << endl;
	// 	}
	// 	cout << endl;
	// }
	while (true){
        int flag = 1;
        for (uint32_t i = 0; i < linearKeys.size() + nonLinearKeys.size(); i++) { if (probArray[i][indices[i]] == 0) flag = 0;}
        if (flag) {
        	int index = indices[0];
        	uint64_t value = probArray[0][indices[0]];
        	for (uint32_t i = 1; i < linearKeys.size() + nonLinearKeys.size(); i++){
        		index *= indices[i];
        		value *= probArray[i][indices[i]];
        	}
            combinedProbArray[index] += value;
        }
        int indexToMove = linearKeys.size() + nonLinearKeys.size()-1;
        while ((indices[indexToMove] >= probLength[indexToMove]-1) && (indexToMove >= 0)) indexToMove--;
        if (indexToMove == -1) break;
        indices[indexToMove]++;
        for (uint32_t i = indexToMove+1; i < linearKeys.size() + nonLinearKeys.size(); i++) indices[i] = 0;
    }
	
	double sum = 0;
	for (uint32_t i = 1; i < combinedMax; i++) sum += combinedProbArray[i];
	cout << "The distribution: " << endl;
	for (uint64_t i = 1; i < combinedMax; i++){
		if (combinedProbArray[i] > 0){
			cout << prob - log2((i+0.0)/combinedProbBase) << " : " << combinedProbArray[i]/(sum)*100 << endl;
		}
	}
}


void propagateKeys()
{
	#if DEBUG == 1
	for (int i = 0; i < linearKeys.size(); i++){
		uint32_t key = linearKeys[i];
		while ((key >> (2*SIZE + 4)) < 4) key = ((key >> (2*SIZE + 4)) + 1) << (2*SIZE+4) ^ (getInvPermSchedule((key >> (2*SIZE))&0xf) << (2*SIZE)) ^ (key & ((1<<(2*SIZE))-1));
		linearKeys[i] = key;
	}

	for (int i = 0; i < nonLinearKeys.size(); i++){
		for (int j = 0; j < nonLinearKeys[i].size(); j++){
			uint32_t key = nonLinearKeys[i][j];
			while ((key >> (2*SIZE + 4)) < 4) key = ((key >> (2*SIZE + 4)) + 1) << (2*SIZE+4) ^ (getInvPermSchedule((key >> (2*SIZE))&0xf) << (2*SIZE)) ^ (key & ((1<<(2*SIZE))-1));
			nonLinearKeys[i][j] = key;
		}
	}
	#else
	for (int i = 0; i < linearKeys.size(); i++){
		uint32_t key = linearKeys[i];
		while ((key >> (2*SIZE + 4)) < nr) key = ((key >> (2*SIZE + 4)) + 1) << (2*SIZE+4) ^ (getInvPermSchedule((key >> (2*SIZE))&0xf) << (2*SIZE)) ^ (key & ((1<<(2*SIZE))-1));
		linearKeys[i] = key;
	}

	for (int i = 0; i < nonLinearKeys.size(); i++){
		for (int j = 0; j < nonLinearKeys[i].size(); j++){
			uint32_t key = nonLinearKeys[i][j];
			while ((key >> (2*SIZE + 4)) < nr) key = ((key >> (2*SIZE + 4)) + 1) << (2*SIZE+4) ^ (getInvPermSchedule((key >> (2*SIZE))&0xf) << (2*SIZE)) ^ (key & ((1<<(2*SIZE))-1));
			nonLinearKeys[i][j] = key;
		}
	}
	#endif
}


bool sizeCheck(vector<uint32_t> constraint)
{
	int key = 0;
	for (uint32_t i = 0; i < constraint.size()-1; i++){
		if (((constraint[i] >> (2*SIZE)) & 0xf) < 8) key++;
	}
	if ((key+1) * SIZE > 36) return false;
	return true;
}

// from here on, change to accomodate SIZE 8
vector<vector<uint32_t>> getEquations(){
	// 0 --> beforeSbox
	// 1 --> afterSbox
	// 2 --> afterKey
	// 3 --> afterSR
	// 4 --> afterMC
	// this computes equations within the same round
	// in M, each element is (y0,y1,y2,y3,x0,x1,x2,x3,k0,k2)
	// the first element is (y0 xor y1) intersect (x2 xor x3)
	vector<vector<uint32_t>> v;
	for (int n = 0; n < nr; n++)
	{
		for (int col = 0; col < 4; col++)
		{
			// uint16_t M[10] = {0b1100001100,0b1010110111,0b1001001000,0b0110111011,0b0101001000,
			// 			  	  0b0011110011,0b1110011111,0b1101100110,0b1011011111,0b1111111111};

			uint16_t M[11] = {0b1100001100,0b1010110111,0b1001000100,0b0110111011,0b0101001000,0b0011110011,
			0b1110011101,0b1101100101,0b1011011111,0b0111010001,0b1111111111};
			for (int i = 0; i < 11; i++)
			{
				vector<uint32_t> tmp;
				for (int j = 0; j < 4; j++) { if (((M[i] >> (9-j)) & 1) == 1) tmp.push_back(((n+1)<<5)+((4*j+col) << 1) + 0);}
				tmp.push_back(-1);
				for (int j = 4; j < 8; j++) { if (((M[i] >> (9-j)) & 1) == 1) tmp.push_back((n<<5)+(  (4*(j-4) + ((12-j+col) % 4))<< 1) + 0);}
				tmp.push_back(-1);
				if ((M[i] >> 1) & 1) tmp.push_back((n<<5)+(col << 1) + 1);
				if ((M[i] >> 0) & 1) tmp.push_back((n<<5)+(4+((3+col) % 4) << 1) + 1);
				auto it = find(v.begin(), v.end(), tmp);
				if (it != v.end()) continue;
				v.push_back(tmp);
			}
		}
			
	}
	return v;
}


// ADD SAME ROUND CONSTANTS
void getSameRoundConstraints(vector<vector<uint32_t>> equations, uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr){
	uint32_t n0,n1,pos0,pos1;
	vector<uint32_t> currentDKey;
	uint32_t currentSKey;
	#if SIZE == 4
	int DDT4[16][16] = {0};
	computeDDT(DDT4);

	#elif SIZE == 8
	int DDT8[256][256] = {0};
	computeDDT8(DDT8);
	#endif
	// first, we detect if all the 'y's are fixed AND all the 'x's are fixed. This will give us a constraint on them or the key involved
	for (int i = 0; i < equations.size(); i++){
		// output side (lhs)
		uint32_t lhs_values[1<<SIZE] = {0}, rhs_values[1<<SIZE] = {0};

		int j = 0;
		bool skip = false;
		while (j < equations[i].size()) // output side
		{
			if (equations[i][j] == -1) break;
			n0 = equations[i][j] >> 5;
			pos0 = (equations[i][j] >> 1) & 0xf;
			if (diffStates[n0][0][pos0] == 0) {
				// if there is no restrictions, all values are allowed
				for (int k = 0; k < (1<<SIZE); k++) lhs_values[k] = 1;
				skip = true;
				break;
			} 
			// if there is a difference, then perhaps there is a subset of values it allows
			uint32_t xddt[1<<SIZE] = {0};
			getXDDT(diffStates[n0][0][pos0],diffStates[n0][1][pos0],xddt);
			XOR(lhs_values,xddt);
			j++;
		}
		if (skip) continue;
		j++;
		while (j < equations[i].size()) {
			// input side
			if (equations[i][j] == -1) break;
			n0 = equations[i][j] >> 5;
			pos0 = (equations[i][j] >> 1) & 0xf;
			if (diffStates[n0][1][pos0] == 0) {
				// if there is no restrictions, all values are allowed
				for (int k = 0; k < (1<<SIZE); k++) rhs_values[k] = 1; 
					skip = true;
				break;
			} 
			uint32_t yddt[1<<SIZE] = {0};
			getYDDT(diffStates[n0][0][pos0],diffStates[n0][1][pos0],yddt);
			XOR(rhs_values,yddt);
			addRoundConstant((n0 << (4+2*SIZE)) + (pos0 << (2*SIZE)) + (diffStates[n0][0][pos0] << SIZE) + diffStates[n0][1][pos0], rhs_values);
			j++;
		}
		if (skip) continue;
		// find the keys that allow this to pass
		uint32_t keys[1<<SIZE] = {0};
		j++;
		// it's either we have no keys at all, or both keys are present
		if (j < equations[i].size()) // there are keys involved
		{
			if (j < equations[i].size() - 1)  // both keys are present
			{
				n0 = equations[i][j] >> 5;
				pos0 = (equations[i][j] >> 1) & 0xf;
				n1 = equations[i][j+1] >> 5;
				pos1 = (equations[i][j+1] >> 1) & 0xf;
				vector<uint32_t> tmp = {(n0 << (4+2*SIZE)) + (pos0 << 2*SIZE) + (diffStates[n0][0][pos0] << SIZE) + (diffStates[n0][1][pos0] << 0),
										 (n1 << (4+2*SIZE)) + (pos1 << 2*SIZE) + (diffStates[n1][0][pos1] << SIZE) + (diffStates[n1][1][pos1] << 0)};
				if (tmp == currentDKey){
					// we XOR the keys together
					XOR(lhs_values,rhs_values);
					for (int k = 0; k < (1<<SIZE); k++) DValues[DValuesLength-1][k] = DValues[DValuesLength-1][k] * lhs_values[k];
					reduceKey(DValues[DValuesLength-1],1);
				}
				else{
					vector<uint32_t> constraint;
					bool swap = false;
					for (int k = 0; k < equations[i].size(); k++)
					{
						if (equations[i][k] == -1) {
							swap = true;
							continue;
						}
						if (equations[i][k] & 1 == 1) continue;
						int n = equations[i][k] >> 5;
						int pos = (equations[i][k] >> 1) & 0xf;
						uint32_t v = (n << (4+2*SIZE)) + (pos << (2*SIZE)) + (diffStates[n][0][pos] << SIZE) + (diffStates[n][1][pos] << 0);
						if (!swap) constraint.push_back(v);
						else constraint.insert(constraint.begin(),v);
					}
					// pushback when there exists a constraint not taken into account yet
					// when no key is involved, does it reduces the probability? 
					bool included = true;
					int n = constraint[constraint.size()-1] >> (4+2*SIZE);
					for (int k = 0; k < constraint.size(); k++)
					{
						if ((constraint[k] >> (4+2*SIZE)) != n) continue;
						auto it = find(keyProbCells.begin(),keyProbCells.end(),constraint[k]);
						if (it == keyProbCells.end()) 
						{
							included = false;
							int r = (constraint[k] >> (2*SIZE)) & 0xf; 
							#if SIZE == 4
							keyProb += -log2(DDT4[diffStates[n][0][r]][diffStates[n][1][r]]/16.0);
							#elif SIZE == 8
							keyProb += -log2(DDT8[diffStates[n][0][r]][diffStates[n][1][r]]/256.0);
							#endif
							keyProbCells.push_back(constraint[k]);
						}
					}
					if (included == false) sameRoundConstraints.push_back(constraint);

					XOR(lhs_values,rhs_values);
					DValuesLength++;
					dualKeys.push_back(tmp);
					currentDKey = tmp;
					for (int k = 0; k < (1<<SIZE); k++) DValues[DValuesLength-1][k] = lhs_values[k];
					reduceKey(DValues[DValuesLength-1],1);
				}
			}
			else {
				n0 = equations[i][j] >> 5;
				pos0 = (equations[i][j] >> 1) & 0xf;
				uint32_t tmp = (n0 << (4+2*SIZE)) + (pos0 << 2*SIZE) + (diffStates[n0][0][pos0] << SIZE) + (diffStates[n0][1][pos0] << 0);
				if (tmp == currentSKey){
					// we XOR the keys together
					XOR(lhs_values,rhs_values);
					for (int k = 0; k < (1<<SIZE); k++) SValues[SValuesLength-1][k] = SValues[SValuesLength-1][k] * lhs_values[k];
					reduceKey(SValues[SValuesLength-1],1);
				}
				else{
					vector<uint32_t> constraint;
					bool swap = false;

					for (int k = 0; k < equations[i].size(); k++)
					{
						if (equations[i][k] == -1){
							swap = true;
							continue;
						}
						if (equations[i][k] & 1 == 1) continue;
						int n = equations[i][k] >> 5;
						int pos = (equations[i][k] >> 1) & 0xf;
						if (!swap) constraint.push_back((n << (4+2*SIZE)) + (pos << (2*SIZE)) + (diffStates[n][0][pos] << SIZE) + (diffStates[n][1][pos] << 0));
						else constraint.insert(constraint.begin(),(n << (4+2*SIZE)) + (pos << (2*SIZE)) + (diffStates[n][0][pos] << SIZE) + (diffStates[n][1][pos] << 0));
					}
					bool included = true;
					int n = constraint[constraint.size()-1] >> (4+2*SIZE);
					for (int k = 0; k < constraint.size(); k++)
					{
						if ((constraint[k] >> (4+2*SIZE)) != n) continue;
						auto it = find(keyProbCells.begin(),keyProbCells.end(),constraint[k]);
						if (it == keyProbCells.end()) 
						{
							included = false;
							int r = (constraint[k] >> (2*SIZE)) & 0xf; 
							#if SIZE == 4
							keyProb += -log2(DDT4[diffStates[n][0][r]][diffStates[n][1][r]]/16.0);
							#elif SIZE == 8
							keyProb += -log2(DDT8[diffStates[n][0][r]][diffStates[n][1][r]]/256.0);
							#endif
							keyProbCells.push_back(constraint[k]);
						}
					}
					if (included == false) sameRoundConstraints.push_back(constraint);
					XOR(lhs_values,rhs_values);
					SValuesLength++;
					singleKeys.push_back(tmp);
					currentSKey = tmp;
					for (int k = 0; k < (1<<SIZE); k++) SValues[SValuesLength-1][k] = lhs_values[k];
					reduceKey(SValues[SValuesLength-1],1);
				}
			}
			
		}
		else // no key: evaluate immediately to see if it's compatible or not
		{
			int sum = 0;
			for (int j = 0; j < (1<<SIZE); j++) sum += lhs_values[j] * rhs_values[j];
			
			cout << endl;
			if (sum == 0) {
				cout << "This trail is invalid due to incompatibilities in round " << n0 << endl;
				cout << "The following equation must hold true but it cannot be satisfied:" << endl;
				string s = "";
				for (int j = 0; j < equations[i].size()-1; j++)
				{
					if (equations[i][j] == -1) s = s.substr(0,s.length()-3) + " = ";
					else s = s + "s_" + to_string((equations[i][j] >> 1) & 0xf) + "^" + to_string(equations[i][j] >> 5) + " + ";
				}
				s = s.substr(0,s.length()-3);
				cout << s << endl;
				cout << "Exiting..." << endl;
				exit(0);
			}
		}
	}
	return;
}

void addLinearConstraint(uint32_t key, uint32_t* val){
	// this function focuses on adding a same round constraint into the linear constraints
	linearKeys.push_back(key);
	// linearValues[linearKeys.size()-1] = val;
	linearDistribution[linearKeys.size()-1] = new uint32_t[1<<SIZE]();
	linearDistributionLength[linearKeys.size()-1] = 1;
	int total = 0;
	for (int i = 0; i < (1<<SIZE); i++){
		linearDistribution[linearKeys.size()-1][i] = val[i];
		
	}
	linearDistributionBase[linearKeys.size()-1] = 1;
	linearDistributionLength[linearKeys.size()-1] = 1;
}

void addNonLinearConstraint(vector<uint32_t> keys, uint32_t* val){
	// this function focuses on adding a same round constraint into the linear constraints
	nonLinearKeys.push_back(keys);
	// expand val into 2D
	uint32_t* newVal; newVal = new uint32_t[1<<(2*SIZE)]();
	for (int i = 0; i < (1<<SIZE); i++){
		for (int j = 0; j < (1<<SIZE); j++){
			newVal[(i << SIZE) + j] = val[i^j];
		}
	}
	nonLinearDistribution[nonLinearKeys.size()-1] = newVal;
	nonLinearDistributionLength[nonLinearKeys.size()-1] = 2;
	nonLinearDistributionBase[nonLinearKeys.size()-1] = 1;
}


void combineConstraints(){
	// we combine constraints that act on the same keys

	// combining nonLinear with linear
	for (uint16_t i = 0; i < nonLinearKeys.size(); i++)
	{
		uint16_t j = 0;
		while (j < linearKeys.size()){
			bool involved = false;
			for (uint16_t k = 0; k < nonLinearKeys[i].size(); k++){
				if ((linearKeys[j] >> (2*SIZE)) == (nonLinearKeys[i][k] >> (2*SIZE))){ 
					// found an intersection. The involved key is in the k^th index
					// cout << hex << "combining linear " << j << " and nonlinear " << i << endl;
					involved = true;
					nonLinearDistributionBase[i] *= linearDistributionBase[j];

					for (uint64_t l1 = 0; l1 < (1<<SIZE); l1++) {
						// loop through all possible key values
						for (uint64_t l2 = 0; l2 < (1ULL << (SIZE*(nonLinearDistributionLength[i]-1))); l2++){
							uint64_t msn = (l2 >> (SIZE*(nonLinearDistributionLength[i] - k - 1))) << (SIZE*(nonLinearDistributionLength[i] - k));
							uint64_t lsn = l2 & ((1ULL << (SIZE*(nonLinearDistributionLength[i] - k - 1))) -1);
							nonLinearDistribution[i][msn+(l1 << (SIZE*(nonLinearDistributionLength[i] - k - 1)))+lsn] *= linearDistribution[j][l1]; // the probabilites here should have a larger base		
						}
					}
				}
			}
			if (involved){
				// if (j == 5 && i == 2)
				// {
				// 	for (int k = 0; k < 256; k++)
				// 	{
				// 		cout << linearDistribution[5][k] << " ";
				// 	}
				// 	cout << endl;
				// 	printout(linearDistributionBase[5]);
				// 	for (int k = 0; k < (1<<(SIZE*nonLinearDistributionLength[2])); k++)
				// 	{
				// 		cout << nonLinearDistribution[2][k] << " ";
				// 	}
				// 	cout << endl;
				// 	cout << nonLinearDistributionBase[2] << endl;
				// }
				// we need to reduce the nonLinearDistribution
				uint64_t divisor = nonLinearDistributionBase[i];
				for (int k = 0; k < (1ULL<<(SIZE*nonLinearDistributionLength[i])); k++)  {
					if (nonLinearDistribution[i][k] > 0) divisor = gcd(divisor,nonLinearDistribution[i][k]);
				}
				for (int k = 0; k < (1ULL<<(SIZE*nonLinearDistributionLength[i])); k++) nonLinearDistribution[i][k] = nonLinearDistribution[i][k]/divisor;
				nonLinearDistributionBase[i] = nonLinearDistributionBase[i]/divisor;

				// printout(nonLinearDistributionBase[2]);


				for (uint16_t l = j; l < linearKeys.size()-1; l++){
					linearDistribution[l] = linearDistribution[l+1];
					linearDistributionLength[l] = linearDistributionLength[l+1];
					linearDistributionBase[l] = linearDistributionBase[l+1];
				}
				linearDistributionLength[linearKeys.size()-1] = 0;
				linearDistributionBase[linearKeys.size()-1] = 0;
				linearKeys.erase(linearKeys.begin()+j);
				j--;
			}
			j++;
		}
	}

	// combining linear with linear
	for (uint32_t i = 0; i < linearKeys.size(); i++){
		uint32_t j = i + 1;
		while (j < linearKeys.size()){
			if ((linearKeys[i] >> (2*SIZE)) == (linearKeys[j] >> (2*SIZE))){
				// cout << "combining linear " << j << " and linear " << i << endl;
				for (uint32_t k = 0; k < (1 << SIZE); k++){				
					linearDistribution[i][k] = linearDistribution[i][k] * linearDistribution[j][k];
				}
				linearDistributionBase[i] *= linearDistributionBase[j];
				for (uint16_t l = j; l < linearKeys.size()-1; l++){
					linearDistribution[l] = linearDistribution[l+1];
					linearDistributionLength[l] = linearDistributionLength[l+1];
					linearDistributionBase[l] = linearDistributionBase[l+1];
				}

				uint64_t divisor = linearDistributionBase[i];
				for (int k = 0; k < (1ULL<<(SIZE*linearDistributionLength[i])); k++)  {
					if (linearDistribution[i][k] > 0) divisor = gcd(divisor,linearDistribution[i][k]);
				}
				for (int k = 0; k < (1ULL<<(SIZE*linearDistributionLength[i])); k++) linearDistribution[i][k] = linearDistribution[i][k]/divisor;
				linearDistributionBase[i] = linearDistributionBase[i]/divisor;

				linearDistributionLength[linearKeys.size()-1] = 0;
				linearDistributionBase[linearKeys.size()-1] = 0;
				linearKeys.erase(linearKeys.begin()+j);
				j--;
			}
			j++;
		}
	}

	uint32_t i = 0; // TBC
	uint32_t j = 0;
	while (i < nonLinearKeys.size()){
		j = i + 1;
		while (j < nonLinearKeys.size()){
			bool combined = false;
			for (uint32_t k0 = 0; k0 < nonLinearKeys[i].size(); k0++){
				for (uint32_t k1 = 0; k1 < nonLinearKeys[j].size(); k1++){
					if ((nonLinearKeys[i][k0] >> (2*SIZE)) == ((nonLinearKeys[j][k1]) >> 2*SIZE)) {
						// if the round and pos are the same. Combine them
						// cout << "combining nonlinear " << j << " and nonlinear " << i << endl;
						combined = combineNonLinearConstraints(i,j);
						if (combined){

							uint64_t divisor = nonLinearDistributionBase[i];
							for (int k = 0; k < (1ULL<<(SIZE*nonLinearDistributionLength[i])); k++)  {
								if (nonLinearDistribution[i][k] > 0) divisor = gcd(divisor,nonLinearDistribution[i][k]);
							}
							for (int k = 0; k < (1ULL<<(SIZE*nonLinearDistributionLength[i])); k++) nonLinearDistribution[i][k] = nonLinearDistribution[i][k]/divisor;
							nonLinearDistributionBase[i] = nonLinearDistributionBase[i]/divisor;
							j -= 1;	
							break;
						}
					}
				}
				if (combined) break;
			}
			j++;
		}
		i++;
	}
}


bool combineNonLinearConstraints(uint32_t index_i, uint32_t index_j){
	vector<uint32_t> combinedKeys;
	vector<uint32_t> overlap_index;
	int unique2 = 0;
	for (uint32_t i = 0; i < nonLinearKeys[index_i].size(); i++) combinedKeys.push_back(nonLinearKeys[index_i][i]);
	for (uint32_t j = 0; j < nonLinearKeys[index_j].size(); j++) {
		// we only combine keys that are diffSet(i,j) but we need to remember the position of the keys of index_j in the combinedKeys vector.
		// This is stored in index2
		bool flag = false;
		for (uint32_t i = 0; i < combinedKeys.size(); i++)
		{
			if ((combinedKeys[i] >> (2*SIZE)) == (nonLinearKeys[index_j][j] >> (2*SIZE))) 
			{
				overlap_index.push_back(i);
				flag = true;
				break;
			}
		}
		if (!flag) {
			combinedKeys.push_back(nonLinearKeys[index_j][j]);
			unique2++;
		}
	}
	if (combinedKeys.size()*SIZE > 36){
		cout << "Memory required by combining nonlinear & nonlinear constraints exceeding 2^36. Recommend to do an experiment instead. Exiting..." << endl;
		exit(0);
	}

	uint32_t tmpDistributionLength = combinedKeys.size();
	uint32_t* tmpDistribution = new uint32_t[1ULL<<(tmpDistributionLength*SIZE)]();

	// main loop to combine the distributions
	// e.g. index_i has key k0,k1,k2 index_j has key k1,k4
	// combined key will have k0,k1,k2,k4
	// firstIndex is the index of index_i
	// secondIndex is the index of index_j
	for (uint64_t i = 0; i < (1ULL << (tmpDistributionLength*SIZE)); i++){
		uint64_t firstIndex = i >> (unique2*SIZE); 
		uint64_t secondIndex = 0; // this is a little harder to compute
		for (uint16_t j = 0; j < overlap_index.size(); j++){
			secondIndex = secondIndex << SIZE;
			secondIndex += (i >> ((tmpDistributionLength-overlap_index[j]-1)*SIZE)) & ((1<<SIZE)-1);
		}
		// we have to add the remaining secondIndex 
		secondIndex = secondIndex << (unique2*SIZE);
		secondIndex ^= i & (1<<(unique2*SIZE));

		tmpDistribution[i] = nonLinearDistribution[index_i][firstIndex] * nonLinearDistribution[index_j][secondIndex];
	}
	nonLinearDistribution[index_i] = tmpDistribution;
	nonLinearDistributionLength[index_i] = tmpDistributionLength;
	nonLinearDistributionBase[index_i] *= nonLinearDistributionBase[index_j];

	for (uint32_t j = index_j; j < nonLinearKeys.size()-1; j++)
	{
		nonLinearDistribution[j] = nonLinearDistribution[j+1];
		nonLinearDistributionLength[j] = nonLinearDistributionLength[j+1];
		nonLinearDistributionBase[j] = nonLinearDistributionBase[j+1];
	}
	delete[] nonLinearDistribution[nonLinearKeys.size()-1];
	nonLinearDistributionLength[nonLinearKeys.size()-1] = 0;
	nonLinearDistributionBase[nonLinearKeys.size()-1] = 0;

	nonLinearKeys.erase(nonLinearKeys.begin()+index_i);
	nonLinearKeys.insert(nonLinearKeys.begin()+index_i,combinedKeys);
	nonLinearKeys.erase(nonLinearKeys.begin()+index_j);
	return true;
}
	