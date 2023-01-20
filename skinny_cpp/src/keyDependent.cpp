
#include "skinny.h"
#include "trails.h"
#include "keyDependent.h"
#include "printDetails.h"
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
vector<uint32_t> keyProbCells; // this cell records the outputs that are already accounted for
vector<uint32_t> output_tuples;
// special constraints
vector<vector<uint32_t>> HOLinearConstraints;
vector<vector<uint32_t>> HOLinearOutput;
vector<vector<uint32_t>> HOLinearKeys;


vector<vector<uint32_t>> constraintsGroupIndex;

map<uint32_t,uint64_t*> savedValues; 

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
					output_tuples.push_back(t[t.size()-1]);
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
void reduceKey(uint64_t* &val, uint32_t valLength) {
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
void findUniqueKeys(vector<uint32_t> &uniqueKeys, vector<vector<uint32_t>> constraints){
	for (int i = 0; i < constraints.size(); i++){
		int maxn = findMaxN(constraints[i]);
		for (int j = 0; j < constraints[i].size(); j++){
			if ((constraints[i][j] >> (4 + 2*SIZE)) == maxn) continue;
			if ((((constraints[i][j] >> (2*SIZE))) & 0xf) >= 8) continue;
			auto it = find(uniqueKeys.begin(),uniqueKeys.end(),constraints[i][j]);
			if (it != uniqueKeys.end()) continue;
			uniqueKeys.push_back(constraints[i][j]);
		}
	}
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
void getHOLinearConstraints(vector<vector<uint32_t>> equations, uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr){
	// equations are all the possible HO linear constraints that one can ever get. 
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
	
	for (int i = 0; i < equations.size(); i++){
		bool fixed = true;
		int n, pos;
		int max_n = 0;
		vector<uint32_t> constrained;
		int count = 0;
		for (int j = 0; j < equations[i].size(); j++){
			if (count >= 2) continue;
			if (equations[i][j] == -1) {
				count++;
				continue;
			}
			n = equations[i][j] >> 5;
			pos = (equations[i][j] >> 1) & 0xf;
			if (diffStates[n][0][pos] == 0) fixed = false;
			constrained.insert(constrained.begin(),(n << (4 + 2*SIZE)) + (pos << (2*SIZE)) + (diffStates[n][0][pos] << SIZE) + diffStates[n][1][pos]);
			if (max_n < n) max_n = n;
		}
		if (fixed){
			// we check if it has been accounted for. For the outputs
			// first, we locate the output terms of this constraints
			vector<uint32_t> output_terms;
			for (int j = 0; j < constrained.size(); j++){
				if (max_n == (constrained[j] >> (4 + 2*SIZE))) output_terms.push_back(constrained[j]);
			}
			// now, we check it against all the known output
			bool allIncluded = true;
			vector<uint32_t> extra_vec; // contains the ones that are not yet included
			for (int k = 0; k < output_terms.size(); k++) {
				auto it = find(output_tuples.begin(),output_tuples.end(), output_terms[k]);
				if (it == output_tuples.end()){
					extra_vec.push_back(output_terms[k]);
				}
				else{
					allIncluded = false;
				}
			}
			// if it's false, we include it in
			// if (allIncluded == true) continue;
			HOLinearConstraints.push_back(constrained);
			HOLinearOutput.push_back(output_terms);
			output_tuples.insert(output_tuples.end(),extra_vec.begin(),extra_vec.end());
			vector<uint32_t> tmpKeys;
			for (int j = 0; j < constrained.size(); j++)
			{
				if (((constrained[j] >> (4+2*SIZE)) < max_n) && (((constrained[j] >> (2*SIZE)) & 0xf) < 8))
				{
					tmpKeys.push_back(constrained[j]);
				}
			}
			HOLinearKeys.push_back(tmpKeys);
			for (int j = 0; j < extra_vec.size(); j++)
			{
				n = extra_vec[j] >> (4+2*SIZE);
				pos = (extra_vec[j] >> (2*SIZE)) & 0xf;
				#if SIZE == 4
				keyProb += -log2(DDT4[diffStates[n][0][pos]][diffStates[n][1][pos]]/16.0);
				#elif SIZE == 8
				keyProb += -log2(DDT8[diffStates[n][0][pos]][diffStates[n][1][pos]]/256.0);
				#endif
			}
		}
	}
	return;
}


void checkHOLinearCompatbility()
{
	for (int i = 0; i < HOLinearKeys.size(); i++)
	{
		// we evaluate directly if there is no key
		// note that only equations involving exactly 2 outputs can have no keys
		// which means there is exactly just 1 uniform variable for 2 constraints to share
		if (HOLinearKeys[i].size() > 0) continue;
		uint32_t uniform[1<<SIZE] = {0}, inp[2][1<<SIZE]={0}, outp[2][1<<SIZE] = {0};
		for (int j = 0; j < (1<<SIZE); j++) uniform[j] = 1;
		for (int j = 0; j < HOLinearOutput[i].size(); j++)
		{
			getXDDT((HOLinearOutput[i][j] >> SIZE) & ANDval,HOLinearOutput[i][j] & ANDval,outp[j]);
			int n = (HOLinearOutput[i][j] >> (4+2*SIZE));
			int pos = ((HOLinearOutput[i][j] >> (2*SIZE)) & 0xf);
			int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
			for (int k = 0; k < 3; k++)
			{
				if (inversePos[k] == -1) continue;
				for (int l = 0; l < HOLinearConstraints[i].size(); l++)
				{
					if ((((n-1) << 4) + inversePos[k]) == (HOLinearConstraints[i][l] >> 2*SIZE)) 
					{
						uint32_t yddt[1<<SIZE] = {0};
						getYDDT((HOLinearConstraints[i][l] >> SIZE) & ANDval,HOLinearConstraints[i][l] & ANDval,yddt);
						XOR(inp[j],yddt);
						addRoundConstant(HOLinearConstraints[i][l],inp[j]);
					}
				}
			}
		}

		// the next step is to check the requirements on the uniform part
		int sum[2] = {0};
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < (1<<SIZE); k++) sum[j] += inp[j][k];
			if (sum[j] == 0)
			{
				// the uniform must take up the output values directly
				for (int k = 0; k < (1<<SIZE); k++) {
					uniform[k] *= outp[j][k];
				}
			}
			else
			{
				// we keep values of uniform in which the xor can reach output[0]
				for (int k = 0; k < (1<<SIZE); k++)
				{
					if (uniform[k] == 0) continue; // if it is not working anymore, don't bother checking
					int flag = 0;
					for (int l = 0; l < (1<<SIZE); l++)
					{
						if (inp[j][l] == 0) continue;
						if (outp[j][l ^ k] >= 1) flag++;
					}
					uniform[k] *= flag;
				}
			}
		}

		// if uniform = 0, there is no solution for this equation
		int sum_uniform = 0;
		int keyCount = 0;
		for (int j = 0; j < (1<<SIZE); j++) 
		{
			sum_uniform += uniform[j];
			if (uniform[j] > 0) keyCount++;
		}
		if (sum_uniform == 0)
		{
			cout << endl; 
			cout << endl;
			cout << "This characteristic is not possible due to the following constraint." << endl;
			printConstraint(HOLinearConstraints[i]);
			// cout << "We will continue to search for more impossibilities..." << endl; 
			cout << endl;
			exit(0);
		}
	}
}

uint32_t findMaxN(vector<uint32_t> constraint)
{
	uint32_t maxN = 0;
	for (int i = 0; i < constraint.size(); i++)
	{
		maxN = max(maxN,constraint[i] >> (4+2*SIZE));
	}
	return maxN;
}

vector<vector<vector<uint32_t>>> splitConstraints(vector<vector<uint32_t>> constraints, vector<vector<uint32_t>> HOLinearConstraints)
{
	vector<vector<vector<uint32_t>>> constraintsGroup;
	for (int i = 0; i < constraints.size() + HOLinearConstraints.size(); i++) {
		vector<uint32_t> tmp; constraintsGroupIndex.push_back(tmp);
	}
	// find out the maxN of all the constraints
	int maxN = 0;
	for (int i = 0; i < constraints.size(); i++)
	{
		int N = findMaxN(constraints[i]);
		if (N > maxN) maxN = N;
	}
	for (int i = 0; i < HOLinearConstraints.size(); i++)
	{
		int N = findMaxN(HOLinearConstraints[i]);
		if (N > maxN) maxN = N;
	}


	// find the constraintsGroupIndex for the constraints
	for (uint32_t i = 0; i < constraints.size(); i++)
	{
		for (uint32_t j = 0; j < constraints[i].size()-1; j++)
		{
			int n = constraints[i][j] >> (4 + 2*SIZE);
			int n_fixed = n;
			int pos = (constraints[i][j] >> (2*SIZE)) & 0xf;
			if (pos >= 8) continue;
			while (n < maxN) {
				n++;
				pos = getInvPermSchedule(pos);
			}
			constraintsGroupIndex[i].push_back((n_fixed << 4) ^ pos);
		}
	}

	// for this, we have to first add in the hidden variables in HOLinearConstraints
	vector<vector<uint32_t>> tmpHOLinearConstraints;
	for (int i = 0; i < HOLinearConstraints.size(); i++)
	{
		int maxN = findMaxN(HOLinearConstraints[i]);
		vector<uint32_t> tmp;
		for (int j = 0; j < HOLinearConstraints[i].size(); j++){
			if ((HOLinearConstraints[i][j] >> (4 + 2 * SIZE)) == maxN){
				tmp.push_back(HOLinearConstraints[i][j]);
				int pos = (HOLinearConstraints[i][j] >> (2*SIZE)) & 0xf;
				int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
				for (int k = 0; k < 3; k++){
					if (inversePos[k] != -1){
						tmp.insert(tmp.begin(),((maxN-1) << 4 + 2*SIZE) + (inversePos[k] << (2*SIZE)) + 
							(diffStates[maxN-1][0][inversePos[k]] << SIZE) + diffStates[maxN-1][1][inversePos[k]]);
					}
				}
			}
		}
		sort( tmp.begin(), tmp.end() );
		tmp.erase( unique( tmp.begin(), tmp.end() ), tmp.end() );
		tmpHOLinearConstraints.push_back(tmp);
	}



	// find the constraintsGroupIndex for the HO constraints
	for (uint32_t i = 0; i < tmpHOLinearConstraints.size(); i++)
	{
		int thisMaxN = findMaxN(tmpHOLinearConstraints[i]);
		for (uint32_t j = 0; j < tmpHOLinearConstraints[i].size(); j++)
		{
			int fixed_n = tmpHOLinearConstraints[i][j] >> (4 + 2*SIZE);
			int n = fixed_n;
			int pos = (tmpHOLinearConstraints[i][j] >> (2*SIZE)) & 0xf;
			if (pos >= 8) continue;
			while (n < maxN) {
				n++;
				pos = getInvPermSchedule(pos);
			}
			constraintsGroupIndex[i+constraints.size()].push_back((fixed_n << 4) ^ pos); // we keep the original round
		}
	}
	// combine constraints and HOLinearConstraints
	constraints.insert(constraints.end(),tmpHOLinearConstraints.begin(),tmpHOLinearConstraints.end());
	for (int i = 0; i < constraints.size(); i++) {
		vector<vector<uint32_t>> tmp2 = {constraints[i]};
		constraintsGroup.push_back(tmp2);
	}

	// check if two entries in constraintIndex have any keys in common. In these cases, we combine them together
	int i = 0;
	while (i < constraintsGroupIndex.size())
	{
		int j = i + 1;
		while (j < constraintsGroupIndex.size())
		{
			bool together = false;
			for (int k = 0; k < constraintsGroupIndex[i].size(); k++){
				for (int l = 0; l < constraintsGroupIndex[j].size(); l++){
					if (constraintsGroupIndex[i][k] == constraintsGroupIndex[j][l]) together = true;
				}
			}
			if (together)
			{
				constraintsGroup[i].insert(constraintsGroup[i].end(),constraintsGroup[j].begin(),constraintsGroup[j].end());
				constraintsGroupIndex[i].insert(constraintsGroupIndex[i].end(),constraintsGroupIndex[j].begin(),constraintsGroupIndex[j].end());
				constraintsGroup.erase(constraintsGroup.begin()+j);
				constraintsGroupIndex.erase(constraintsGroupIndex.begin()+j);
				j--;
			}
			j++;
		}
		i++;
	}

	// next, we look at the key schedule. If the TK value is actually < number of different rounds for a single key, we have to combine them
	int A[16][30] = {0};
	for (int i = 0; i < constraintsGroup.size(); i++)
	{
		for (int j = 0; j < constraintsGroupIndex[i].size(); j++)
		{
			A[constraintsGroupIndex[i][j]&0xf][constraintsGroupIndex[i][j]>>4] = 1;
		}
	}
	for (int i = 0; i < 16; i++)
	{
		int sum = 0;
		for (int j = 0; j < 30; j++){
			sum += A[i][j];
		}
		if (sum > TK_NUM)
		{
			int k = 0;
			// first, we find the first constraintGroup that has this key. Then we combine any other groups with this key into this constraintGroup
			bool breakFlag = false;
			for (k = 0; k < constraintsGroupIndex.size(); k++)
			{
				for (int l = 0; l < constraintsGroupIndex[k].size(); l++){
					if ((constraintsGroupIndex[k][l] & 0xf) == i) {
						breakFlag = true;
						break;
					}
				}
				if (breakFlag) break;
			}
			int l = k + 1;
			while (l < constraintsGroupIndex.size())
			{
				bool together = false;
				for (int m = 0; m < constraintsGroupIndex[k].size(); m++){
					for (int n = 0; n < constraintsGroupIndex[l].size(); n++){
						if ((constraintsGroupIndex[k][m] & 0xf) == (constraintsGroupIndex[l][n] & 0xf)) {
							together = true;
						}
					}
				}
				if (together)
				{
					constraintsGroup[k].insert(constraintsGroup[k].end(),constraintsGroup[l].begin(),constraintsGroup[l].end());
					constraintsGroupIndex[k].insert(constraintsGroupIndex[k].end(),constraintsGroupIndex[l].begin(),constraintsGroupIndex[l].end());
					constraintsGroup.erase(constraintsGroup.begin()+l);
					constraintsGroupIndex.erase(constraintsGroupIndex.begin()+l);
					l = k;
				}
				l++;
			}
		}
	}


	// next, we check if how many overlaps of non-key tuples. If there are more than 1 overlap, we have to combine them
	vector<vector<uint32_t>> elements;
	for (int i = 0; i < constraintsGroup.size(); i++){
		vector<uint32_t> tmp;
		for (int j = 0; j < constraintsGroup[i].size(); j++){
			for (int k = 0; k < constraintsGroup[i][j].size(); k++){
				tmp.push_back(constraintsGroup[i][j][k]);
			}
		}
		elements.push_back(tmp);
	}
	
	i = 0;
	cout << endl;
	while (i < elements.size())
	{
		int j = i + 1;
		while (j < elements.size())
		{
			vector<uint32_t> together;
			for (int k = 0; k < elements[i].size(); k++){
				for (int l = 0; l < elements[j].size(); l++){
					if (elements[i][k] == elements[j][l]) 
					{
						auto it = find(together.begin(),together.end(),elements[i][k]);
						if (it == together.end()) together.push_back(elements[i][k]);
					}
				}
			}
			if (together.size() > 1)
			{
				constraintsGroup[i].insert(constraintsGroup[i].end(),constraintsGroup[j].begin(),constraintsGroup[j].end());
				constraintsGroupIndex[i].insert(constraintsGroupIndex[i].end(),constraintsGroupIndex[j].begin(),constraintsGroupIndex[j].end());
				constraintsGroup.erase(constraintsGroup.begin()+j);
				constraintsGroupIndex.erase(constraintsGroupIndex.begin()+j);
				elements[i].insert(elements[i].end(),elements[j].begin(),elements[j].end());
				elements.erase(elements.begin()+j);
				j = i;
			}
			j++;
		}
		i++;
	}

	for (int i = 0; i < constraintsGroupIndex.size(); i++)
	{
		sort( constraintsGroupIndex[i].begin(), constraintsGroupIndex[i].end() );
		constraintsGroupIndex[i].erase( unique( constraintsGroupIndex[i].begin(), constraintsGroupIndex[i].end() ), constraintsGroupIndex[i].end() );
	}
	return constraintsGroup;
}

bool computable(vector<vector<uint32_t>> constraints)
{
	// if computationally too intensive, return false
	#if SIZE == 4
	int DDT4[16][16] = {0};
	computeDDT(DDT4);

	#elif SIZE == 8
	int DDT8[256][256] = {0};
	computeDDT8(DDT8);
	#endif
	vector<uint32_t> inputs, others, keys;
	double runTime = 0;
	for (int i = 0; i < constraints.size(); i++){
		int maxN = findMaxN(constraints[i]);
		for (int j = 0; j < constraints[i].size(); j++){
			if ((constraints[i][j] >> (4+2*SIZE)) == maxN) others.push_back(constraints[i][j]);
			else if ((constraints[i][j] & ANDval) == 0) others.push_back(constraints[i][j]);
		}
	}

	for (int i = 0; i < constraints.size(); i++){
		for (int j = 0; j < constraints[i].size(); j++){
			auto it = find(others.begin(),others.end(),constraints[i][j]);
			auto it2 = find(inputs.begin(),inputs.end(),constraints[i][j]);
			if (it == others.end() && it2 == inputs.end()){
				if (((constraints[i][j] >> (2*SIZE)) & 0xf) < 8) keys.push_back(constraints[i][j]);
				inputs.push_back(constraints[i][j]);
				uint8_t inp = (constraints[i][j] >> SIZE) & ANDval;
				uint8_t outp = constraints[i][j] & ANDval;
				#if SIZE == 4
				runTime += log2(DDT4[inp][outp]);
				#elif SIZE == 8
				runTime += log2(DDT8[inp][outp]);
				#endif
			}
			if (((constraints[i][j] & ANDval) == 0) && (((constraints[i][j] >> 2*SIZE) & 0xf) < 8)) keys.push_back(constraints[i][j]);
		}
	}

	// check what is the size of the keys
	int maxN = findMaxN(keys);
	int index[16] = {0};
	vector<vector<uint32_t>> KEYS;
	for (int i = 0; i < 16; i++) index[i] = -1;
	for (int i = 0; i < keys.size(); i++)
	{
		int n = keys[i] >> (4+2*SIZE);
		int pos = (keys[i] >> (2*SIZE)) & 0xf;
		while (n < maxN) {
			n++;
			pos = getInvPermSchedule(pos);
		}
		
		if (index[pos] == -1)
		{
			index[pos] = KEYS.size();
			vector<uint32_t> key1 = {keys[i]};
			KEYS.push_back(key1);
		}
		else
		{
			int start_n = KEYS[index[pos]][0] >> (4+2*SIZE);
			if (start_n < n)
			{ 
				while (start_n < n-2)
				{
					KEYS[index[pos]].push_back(-1);
					start_n += 2;
				}
				KEYS[index[pos]].push_back(keys[i]);
			}
		}
	}
	for (int i = 0; i < KEYS.size(); i++){
		if (KEYS[i].size() == 1) runTime += SIZE;
		else runTime += TK_NUM * SIZE;
	}
	if (runTime > 28) return false;
	return true;
}