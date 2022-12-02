#include "gift.h"
#include "trails.h"
#include "gaussianElimination.h"
#include "SAT.h"
#include "keyDependent.h"


#if !defined SIZE
#error SIZE is not defined.
#endif

using namespace std;

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl

vector<vector<int>> keyConstraints;
int saved_prob = 0;
void printState(uint8_t alpha[30][16], uint32_t nr){
	/*
	This function prints out the before and after substitution values for each round of the characreristic
	*/
	uint8_t beforeSbox[16], afterSbox[16];
	for (int n = 0; n < nr; n++)
	{
		for (int i = 0; i < 16; i++)
		{
			beforeSbox[i] = alpha[n][i];
			afterSbox[i] = alpha[n+1][i];
		}
		#if SIZE == 4
		invPermutation(afterSbox);
		#elif SIZE == 8
		invPermutation8(afterSbox);
		#endif
		for (int i = 0; i < 16; i++) cout << setw(SIZE/4) << hex << int(beforeSbox[i]) << " ";
		cout << endl;
		for (int i = 0; i < 16; i++) cout << setw(SIZE/4) << hex << int(afterSbox[i]) << " ";
		cout << endl;
		cout << "----------------------------------------------------------------" << endl;
	}
}

// void findEquations(uint8_t xEquations[16][16][32], uint8_t yEquations[16][16][32]){
// 	/*
// 	This function helps to find all the possible linear equations associating with each sbox difference transition
// 	For instance, if an Sbox difference is assumed to be 1 --> 3 and a linear equation arose from it is x0 + x2 = 1
// 	Then we will have xEquations[1][3][0b10100] = 1 
// 	note the breakdown of the last array index of xEquation/yEquation is x0|x1|x2|x3|b
// 	*/
// 	uint32_t XDDT[16] = {0};
// 	uint32_t YDDT[16] = {0};
// 	bool flagY, flagX;
// 	for (int i = 0; i < 16; i++)
// 	{
// 		for (int j = 0; j < 16; j++)
// 		{
// 			computeYDDT(i,j,YDDT);
// 			computeXDDT(i,j,XDDT);
// 			vector<int> tmpX,tmpY// 			for (int k = 0; k < 16; k++)
// 			{
// 				if (YDDT[k] > 0) tmpY.push_back(k);
// 				if (XDDT[k] > 0) tmpX.push_back(k);
// 			}
// 			if (tmpY.size() > 0)
// 			{
// 				for (int b = 1; b < 16; b++)
// 				{
// 					flagY = true;
// 					for (uint32_t z = 0; z < tmpY.size(); z++) {if ((__builtin_popcount(tmpY[z] & b) & 1) != 0) flagY = false;}
// 					if (flagY) yEquations[i][j][(b<<1)+0] = 1;
// 					flagY = true;
// 					for (uint32_t z = 0; z < tmpY.size(); z++) {if ((__builtin_popcount(tmpY[z] & b) & 1) != 1) flagY = false;}
// 					if (flagY) yEquations[i][j][(b<<1)+1] = 1;
// 				}
// 			}
// 			if (tmpX.size() > 0)
// 			{
// 				for (int b = 1; b < 16; b++)
// 				{
// 					flagX = true;
// 					for (uint32_t z = 0; z < tmpX.size(); z++) { if ((__builtin_popcount(tmpX[z] & b) & 1) != 0) flagX = false;}
// 					if (flagX) xEquations[i][j][(b<<1)+0] = 1;
// 					flagX = true;
// 					for (uint32_t z = 0; z < tmpY.size(); z++) { if ((__builtin_popcount(tmpX[z] & b) & 1) != 1) flagX = false;}
// 					if (flagX) xEquations[i][j][(b<<1)+1] = 1;
// 				}
// 			}
// 		}
// 	}
// }


void findEquations(uint8_t xEquations[16][16][32], uint8_t yEquations[16][16][32]){
	/*
	This function helps to find all the possible linear equations associating with each sbox difference transition
	For instance, if an Sbox difference is assumed to be 1 --> 3 and a linear equation arose from it is x0 + x2 = 1
	Then we will have xEquations[1][3][0b10100] = 1 
	note the breakdown of the last array index of xEquation/yEquation is x0|x1|x2|x3|b
	*/
	uint32_t XDDT[16] = {0};
	uint32_t YDDT[16] = {0};
	bool flagY, flagX;
	for (int i = 1; i < 16; i++)
	{
		for (int j = 1; j < 16; j++)
		{
			computeYDDT(i,j,YDDT);
			computeXDDT(i,j,XDDT);
			vector<int> tmpX,tmpY;
			for (int k = 0; k < 16; k++)
			{
				if (YDDT[k] > 0) tmpY.push_back(k);
				if (XDDT[k] > 0) tmpX.push_back(k);
			}
			if (tmpY.size() > 0)
			{
				uint8_t **M; M = new uint8_t* [32];
				for (int k = 0; k < 32; k++) M[k] = new uint8_t [5]();
				int count = 0;
				for (int b = 1; b < 16; b++)
				{
					flagY = true;
					for (uint32_t z = 0; z < tmpY.size(); z++) {if ((__builtin_popcount(tmpY[z] & b) & 1) != 0) flagY = false;}
					if (flagY)
					{
						M[count][0] = (b>>3)&0b1;
						M[count][1] = (b>>2)&0b1;
						M[count][2] = (b>>1)&0b1;
						M[count][3] = (b>>0)&0b1;
						M[count][4] = 0;
						count++;
					}
					
					flagY = true;
					for (uint32_t z = 0; z < tmpY.size(); z++) {if ((__builtin_popcount(tmpY[z] & b) & 1) != 1) flagY = false;}
					if (flagY)
					{
						M[count][0] = (b>>3)&0b1;
						M[count][1] = (b>>2)&0b1;
						M[count][2] = (b>>1)&0b1;
						M[count][3] = (b>>0)&0b1;
						M[count][4] = 1;
						count++;
					}
					
				}
				GaussianElimination(M,32,5);


				for (int r = 0; r < 32; r++)
				{
					uint8_t b = 0;
					for (int c = 0; c < 5; c++)
					{
						b = b << 1;
						b += M[r][c];
					}
					if (b > 0) yEquations[i][j][b] = 1;
				}
				for (int k = 0; k < 32; k++) delete[] M[k];
				delete[] M;
			}
			if (tmpX.size() > 0)
			{
				uint8_t **M; M = new uint8_t* [32];
				for (int k = 0; k < 32; k++) M[k] = new uint8_t [5]();
				int count = 0;
				for (int b = 1; b < 16; b++)
				{
					flagX = true;
					for (uint32_t z = 0; z < tmpX.size(); z++) {if ((__builtin_popcount(tmpX[z] & b) & 1) != 0) flagX = false;}
					if (flagX)
					{
						M[count][0] = (b>>3)&0b1;
						M[count][1] = (b>>2)&0b1;
						M[count][2] = (b>>1)&0b1;
						M[count][3] = (b>>0)&0b1;
						M[count][4] = 0;
						count++;
					}
					flagX = true;
					for (uint32_t z = 0; z < tmpX.size(); z++) {if ((__builtin_popcount(tmpX[z] & b) & 1) != 1) flagX = false;}
					if (flagX)
					{
						M[count][0] = (b>>3)&0b1;
						M[count][1] = (b>>2)&0b1;
						M[count][2] = (b>>1)&0b1;
						M[count][3] = (b>>0)&0b1;
						M[count][4] = 1;
						count++;
					}
					
				}
				GaussianElimination(M,32,5);

				for (int r = 0; r < 32; r++)
				{
					uint8_t b = 0;
					for (int c = 0; c < 5; c++)
					{
						b = b << 1;
						b += M[r][c];
					}
					if (b > 0) xEquations[i][j][b] = 1;
				}
				for (int k = 0; k < 32; k++) delete[] M[k];
				delete[] M;
			}
			
		}
	}
}

vector<vector<int>> findYConstraints(uint8_t alpha[30][16], int nr, vector<vector<vector<int>>> &positions){
	// finding all the y (output) constraints via YDDT
	// for y constraints, we apply the permutation directly on them.
	uint8_t beforeSbox[16] = {0};
	uint8_t afterSbox[16] = {0};
	uint8_t yEquations[16][16][32] = {0}, xEquations[16][16][32] = {0};
	findEquations(xEquations,yEquations); // find the equations associated with each differential transition
	vector<vector<int>> parity;
	for (int n = 0; n < (nr-1) * 32; n++) {
		// initializing a vector of vectors
		vector<vector<int>> tmp;
		vector<int> tmp2;
		positions.push_back(tmp);
		parity.push_back(tmp2);
	}

	for (int n = 0; n < nr-1; n++){
		for (int c = 0; c < 16; c++){
			beforeSbox[c] = alpha[n][c]; afterSbox[c] = alpha[n+1][c]; // need to be inversed
		}
		#if SIZE == 4
		invPermutation(afterSbox);
		for (int c = 0; c < 16; c++){// for each Sbox
			if (beforeSbox[c] == 0) continue; // skip if inactive
			for (int j = 0; j < 32; j++){
				if (yEquations[beforeSbox[c]][afterSbox[c]][j] > 0) // if there is an equation, store it in positions and parity
				{
					vector<int> tmp;
 					for (int k = 4; k >= 1; k--) // k = 0 refers to the the parity (ignored here)
					{
						// position[round*32+sbox_position] = round * 128 + bit_position (63.......0)
						if ((j >> k) & 0b1) tmp.push_back(int(n * 128) + int(getPerm(63 - 4*c - (4-k))));
					}
					positions[n*32 + 15-c].push_back(tmp);
					parity[n*32 + 15-c].push_back(j & 0b1);
				}
			}
		}
		#elif SIZE == 8
		invPermutation8(afterSbox);
		for (int c = 0; c < 16; c++) // for each pair of Sboxes
		{
			uint8_t topBeforeSbox = beforeSbox[c] >> 4;
			uint8_t topAfterSbox = afterSbox[c] >> 4;
			uint8_t botBeforeSbox = beforeSbox[c] & 0xf;
			uint8_t botAfterSbox = afterSbox[c] & 0xf;
			if (topBeforeSbox != 0)
			{
				for (int j = 0; j < 32; j++)
				{
					if (yEquations[topBeforeSbox][topAfterSbox][j] > 0)
					{
						vector<int> tmp;
	 					for (int k = 4; k >= 1; k--)
						{
							if ((j >> k) & 0b1) tmp.push_back(int(n * 128) + int(getPerm8(127 - 8*c - (4-k) )));
						}
						positions[n*32 + 31-2*c].push_back(tmp);
						parity[n*32 + 31-2*c].push_back(j & 0b1);
					}
				}
			}
			if (botBeforeSbox != 0)
			{
				for (int j = 0; j < 32; j++)
				{
					if (yEquations[botBeforeSbox][botAfterSbox][j] > 0)
					{
						vector<int> tmp;
	 					for (int k = 4; k >= 1; k--)
						{
							if ((j >> k) & 0b1) tmp.push_back(int(n * 128) + int(getPerm8(127 - 8*c - (8-k))));
						}
						positions[n*32 + 31-2*c-1].push_back(tmp);
						parity[n*32 + 31-2*c-1].push_back(j & 0b1);
					}
				}
			}
		}
		#endif
	}
	return parity;
}

vector<vector<int>> findXConstraints(uint8_t alpha[30][16], int nr, vector<vector<vector<int>>> &positions){
	// finding all the constraints via YDDT
	uint8_t beforeSbox[16] = {0};
	uint8_t afterSbox[16] = {0};
	uint8_t xEquations[16][16][32] = {0}, yEquations[16][16][32] = {0};
	findEquations(xEquations,yEquations); // find the equations associated with each differential transition
	vector<vector<int>> parity;
	for (int n = 0; n < (nr-1) * 32; n++) // initialization
	{
		vector<vector<int>> tmp;
		vector<int> tmp2;
		positions.push_back(tmp);
		parity.push_back(tmp2);
	}
	for (int n = 1; n < nr; n++)
	{
		for (int i = 0; i < 16; i++)
		{
			beforeSbox[i] = alpha[n][i];
			afterSbox[i] = alpha[n+1][i]; // needs to be inversed 
		}
		#if SIZE == 4
		invPermutation(afterSbox);
		for (int i = 0; i < 16; i++)
		{
			if (beforeSbox[i] == 0) continue;
			for (int j = 0; j < 32; j++)
			{
				if (xEquations[beforeSbox[i]][afterSbox[i]][j] > 0)
				{
					vector<int> tmp;
 					for (int k = 4; k >= 1; k--)
					{
						// the round numbers here are stored such that we match the previous round's y constraints
						if ((j >> k) & 0b1) tmp.push_back((n-1) * 128 + int(63 - 4*i - (4-k))); 
					}
					positions[(n-1)*32 + 15-i].push_back(tmp);
					parity[(n-1)*32 + 15-i].push_back(j & 0b1); 
				}
			}
		}
		#elif SIZE == 8
		invPermutation8(afterSbox);
		for (int i = 0; i < 16; i++)
		{
			uint8_t topBeforeSbox = beforeSbox[i] >> 4;
			uint8_t topAfterSbox = afterSbox[i] >> 4;
			uint8_t botBeforeSbox = beforeSbox[i] & 0xf;
			uint8_t botAfterSbox = afterSbox[i] & 0xf;
			if (topBeforeSbox != 0)
			{
				for (int j = 0; j < 32; j++)
				{
					if (xEquations[topBeforeSbox][topAfterSbox][j] > 0)
					{
						vector<int> tmp;
	 					for (int k = 4; k >= 1; k--)
						{
							// the round numbers here are stored such that we match the previous round's y constraints
							if ((j >> k) & 0b1) tmp.push_back((n-1) * 128 + int(127 - 8*i - (4-k)));
						}
						positions[(n-1)*32 + 31-2*i].push_back(tmp);
						parity[(n-1)*32 + 31-2*i].push_back(j & 0b1);
					}
				}
			}
			if (botBeforeSbox != 0)
			{
				for (int j = 0; j < 32; j++)
				{
					if (xEquations[botBeforeSbox][botAfterSbox][j] > 0)
					{
						vector<int> tmp;
	 					for (int k = 4; k >= 1; k--)
						{
							// the round numbers here are stored such that we match the previous round's y constraints
							if ((j >> k) & 0b1) tmp.push_back((n-1) * 128 + int(127 - 8*i - (8-k)));
						}
						positions[(n-1)*32 + 31-2*i-1].push_back(tmp);
						parity[(n-1)*32 + 31-2*i-1].push_back(j & 0b1);
					}
				}
			}
		}
		#endif
	}
	return parity;
}

vector<int> checkXPos(vector<int> constraint){
	// find the Sbox and round number from each position (uniquness)
	// constraint elements are in the form of round * 128 + bit_position
	// return a vector showing round * 32 + sbox_position 15...0
	vector<int> pos;
	for (uint32_t i = 0; i < constraint.size(); i++)
	{
		int sbox = (constraint[i] & 0b1111111) / 4; 
		int n = constraint[i] >> 7;
		 if (std::find(pos.begin(), pos.end(), n*32+sbox) == pos.end())
		 {
		 	pos.push_back(n*32+sbox);
		 }
	}
	return pos;
}


vector<int> vectorXOR(vector<int> a, vector<int> b) {
	// this has a maximum of 128, preserving the round number
	// XOR-ing two vectors. If they have the same variable even number of times, then they cancel them off
	int M[64 * SIZE/4] = {0};
	for (uint32_t i = 0; i < a.size(); i++) M[a[i] & 0b1111111] ^= 1;
	for (uint32_t i = 0; i < b.size(); i++) M[b[i] & 0b1111111] ^= 1;
	// M[i] = 1 if ith element should be in the output
	int n;
	if (a.size() > 0) n = a[0] >> 7;
	else if (b.size() > 0) n = b[0] >> 7; // get the round number
	vector<int> tmp;
	for (int i = 0; i < (64 * SIZE/4); i++) {if (M[i] == 1) {tmp.push_back(n*128+i);}}
	return tmp;
}


// here 


void insertKeyConstraint(vector<int> vec, int parity)
{
	vector<int> keys; 
	vector<int> cons;
	for (uint32_t i = 0; i < vec.size(); i++)
	{
		#if SIZE == 4
		if ((((vec[i] & 0b1111111) % 4) == 0) || (((vec[i] & 0b1111111) % 4) == 1)){
			if (std::find(keys.begin(), keys.end(), vec[i]) == keys.end()){
				keys.push_back(vec[i]);
			}
		} 
		if ((vec[i] & 0b1111111) == 63) cons.push_back(vec[i]);
		#elif SIZE == 8
		if ((((vec[i] & 0b1111111) % 4) == 1) || (((vec[i] & 0b1111111) % 4) == 2)) {
			if (std::find(keys.begin(), keys.end(), vec[i]) == keys.end()){
				keys.push_back(vec[i]);
			}
		}
		if ((vec[i] & 0b1111111) == 127) cons.push_back(vec[i]);
		#endif

		if ((vec[i] & 0b1111111) == 23) cons.push_back(vec[i]);
		if ((vec[i] & 0b1111111) == 19) cons.push_back(vec[i]);
		if ((vec[i] & 0b1111111) == 15) cons.push_back(vec[i]);
		if ((vec[i] & 0b1111111) == 11) cons.push_back(vec[i]);
		if ((vec[i] & 0b1111111) ==  7) cons.push_back(vec[i]);
		if ((vec[i] & 0b1111111) ==  3) cons.push_back(vec[i]);
	}
	for (int i = 0; i < vec.size(); i++)
	{
		cout << dec << (vec[i]>>7) << "," << (vec[i] & 0b1111111) << " ";
	}
	cout << endl;
	int con_parity = 0; // compute the parity for the constants (if any)
	for (int i = 0; i < cons.size(); i++)
	{
		int n = cons[i] >> 7;
		int pos = (cons[i] & 0b1111111);
		if (pos == (1<<SIZE)-1) con_parity ^= 1;
		if (pos == 23) con_parity ^= (getConstants(n) >> 5) & 1;
		if (pos == 19) con_parity ^= (getConstants(n) >> 4) & 1;
		if (pos == 15) con_parity ^= (getConstants(n) >> 3) & 1;
		if (pos == 11) con_parity ^= (getConstants(n) >> 2) & 1;
		if (pos ==  7) con_parity ^= (getConstants(n) >> 1) & 1;
		if (pos ==  3) con_parity ^= (getConstants(n) >> 0) & 1;
	}
	if (keys.size() == 0)
	{
		if (cons.size() == 0) cout << "no key/constants are found for this constraint" << endl;
		else cout << "a constant is found for this constraint" << endl;
		// immediately evaluate if the constraint is satisfied immediately
		if (con_parity != parity) cout << "constraint does not satisfy in any case!" << endl;
		else {
			cout << "constraint is directly satisfied. Save 2^-1 prob here!" << endl;
			saved_prob += 1;
		}
	}
	else if (keys.size() > 0 && cons.size() == 0)
	{
		for (int i = 0; i < keys.size(); i++) cout << (keys[i] >> 7) << "," << (keys[i] & 0b1111111) << " ";
		cout << endl;
		if (cons.size() == 0) cout << "no constants are found for this constraint but a key is involved" << endl;
		else cout << "constants and keys are involved for this constraint" << endl;
		keys.push_back(parity ^ con_parity);
		// auto it = find(keyConstraints.begin(), keyConstraints.end(), tmp);
		// if (it != keyConstraints.end()) return; // add in only if it doesn't exist
        keyConstraints.push_back(keys);
	}

}

// void insertKeyConstraint(vector<int> vec, int parity)
// {
// 	/*
// 	This function will check if the elements in the vector will pass through a key/constraint/nothing
// 	It will print out the appropriate message
// 	*/
// 	vector<int> tmp; 
// 	// this will detect if a key is involved and put it into tmp
// 	for (uint32_t i = 0; i < vec.size(); i++)
// 	{
// 		#if SIZE == 4
// 		if ((((vec[i] & 0b1111111) % 4) == 0) || (((vec[i] & 0b1111111) % 4) == 1)) 
// 		{
// 			if (std::find(tmp.begin(), tmp.end(), vec[i]) == tmp.end()){
// 				tmp.push_back(vec[i]);
// 			}
// 		}
// 		#elif SIZE == 8
// 		if ((((vec[i] & 0b1111111) % 4) == 1) || (((vec[i] & 0b1111111) % 4) == 2))
// 		{
// 			if (std::find(tmp.begin(), tmp.end(), vec[i]) == tmp.end()){
// 				tmp.push_back(vec[i]);
// 			}
// 		}
// 		#endif
// 	}
// 	printout("tmp");
// 	for (int i = 0; i < tmp.size(); i++)
// 	{
// 		cout << (tmp[i] >> 7) << "," << (tmp[i] & 0b1111111) << " ";
// 	}
// 	cout << endl;

// 	// check from here onwards. But first, check what's in vec


// 	if ((tmp.size() == 0) && (parity == 1))
// 	{
// 		bool constantFlag = false;
// 		for (uint32_t i = 0; i < tmp.size(); i++)
// 		{
// 			if ((tmp[i] == 64*(SIZE/4)-1) || (tmp[i] == 23) || (tmp[i] == 19) || (tmp[i] == 15) || 
// 				(tmp[i] == 11) || (tmp[i] == 7) ||  (tmp[i] == 3)) constantFlag = true;
// 		}
// 		cout << "Round constant is involved, but ok (a different round constant may not work)" << endl;
// 	}
// 	if ((tmp.size() == 0) && (parity == 1)) // if it's incompatible 
// 	{
// 		for (uint32_t i = 0; i < vec.size(); i++) cout << (vec[i] >> 7) << "," << (vec[i] & 0b1111111) << " ";
// 		cout << endl;
// 		bool constantFlag = false;
// 		for (uint32_t i = 0; i < tmp.size(); i++)
// 		{
// 			if ((tmp[i] == 64*(SIZE/4)-1) || (tmp[i] == 23) || (tmp[i] == 19) || (tmp[i] == 15) || 
// 				(tmp[i] == 11) || (tmp[i] == 7) ||  (tmp[i] == 3)) constantFlag = true;
// 		}
// 		if (constantFlag) cout << "Round constant incompatibility occurs (a different round constant may help)" << endl;
// 		else cout << "Incompatible due to Sbox restrictions (no key nor round constant can help)" << endl;
		
// 	}
// 	else if (tmp.size() > 0) // key dependent constraints
// 	{
// 		tmp.push_back(parity);
// 		auto it = find(keyConstraints.begin(), keyConstraints.end(), tmp);
// 		if (it != keyConstraints.end()) return; // add in only if it doesn't exist
//         keyConstraints.push_back(tmp);
// 	}
// }

void searchOne(vector<int> pos, vector<int> constraint, int parity, vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity)
{
	for (uint32_t i = 0; i < xPositions[pos[0]].size(); i++)
	{
		// cout << "comparing: ";
		// for (uint32_t j = 0; j < constraint.size(); j++) cout << (constraint[j] >> 7) << "," << (constraint[j] & 0b1111111) << " ";
		// cout << "p: " << parity << " ";
		// cout << endl << "xpos: ";
		// for (uint32_t j = 0; j < xPositions[pos[0]][i].size(); j++) cout << (xPositions[pos[0]][i][j] >> 7) << "," << (xPositions[pos[0]][i][j] & 0b1111111) << " ";
		// cout << "p: " << xParity[pos[0]][i] << endl;
		vector<int> tmp = vectorXOR(constraint,xPositions[pos[0]][i]);
		if (tmp.size() == 0) 
		{
			// cout << "got in here --------------------------" << endl;
			parity ^= xParity[pos[0]][i];
			vector<int> tmp2 = constraint;
			tmp2.insert(tmp2.end(),xPositions[pos[0]][i].begin(),xPositions[pos[0]][i].end());
			insertKeyConstraint(tmp2,parity);
		}
	}
	return;
}

void searchTwo(vector<int> pos, vector<int> constraint, int parity, vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity)
{
	for (uint32_t i0 = 0; i0 < xPositions[pos[0]].size(); i0++)
	{
		for (uint32_t i1 = 0; i1 < xPositions[pos[1]].size(); i1++)
		{
			// cout << "comparing: ";
			// for (uint32_t j = 0; j < constraint.size(); j++) cout << (constraint[j] >> 6) << "," << (constraint[j] & 0b111111) << " ";
			// cout << "p: " << parity << " " << endl;
			// cout << "xpos0: ";
			// for (uint32_t j = 0; j < xPositions[pos[0]][i0].size(); j++) cout << (xPositions[pos[0]][i0][j] >> 6) << "," << (xPositions[pos[0]][i0][j] & 0b111111) << " ";
			// cout << "p: " << xParity[pos[0]][i0] << " " << endl;
			// cout << "xpos1: ";
			// for (uint32_t j = 0; j < xPositions[pos[1]][i1].size(); j++) cout << (xPositions[pos[1]][i1][j] >> 6) << "," << (xPositions[pos[1]][i1][j] & 0b111111) << " ";
			// cout << "p: " << xParity[pos[1]][i1] << endl;
			vector<int> tmp = vectorXOR(constraint,xPositions[pos[0]][i0]);
			tmp = vectorXOR(tmp,xPositions[pos[1]][i1]);

			if (tmp.size() == 0) // only if they XOR to zero --> the equation will be formed 
			{
				// cout << "got in here --------------------------" << endl;
				parity ^= xParity[pos[0]][i0] ^ xParity[pos[1]][i1];
				vector<int> tmp2 = constraint;
				tmp2.insert(tmp2.end(),xPositions[pos[0]][i0].begin(),xPositions[pos[0]][i0].end());
				tmp2.insert(tmp2.end(),xPositions[pos[1]][i1].begin(),xPositions[pos[1]][i1].end());
				insertKeyConstraint(tmp2,parity);
			}
		}
	}
	return;
}

void searchThree(vector<int> pos, vector<int> constraint, int parity, vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity)
{
	for (uint32_t i0 = 0; i0 < xPositions[pos[0]].size(); i0++)
	{
		for (uint32_t i1 = 0; i1 < xPositions[pos[1]].size(); i1++)
		{
			for (uint32_t i2 = 0; i2 < xPositions[pos[2]].size(); i2++)
			{
				// cout << "comparing: ";
				// for (uint32_t j = 0; j < constraint.size(); j++) cout << (constraint[j] >> 6) << "," << (constraint[j] & 0b111111) << " ";
				// cout << "p: " << parity << " ";
				// cout << endl << "xpos: ";
				// for (uint32_t j = 0; j < xPositions[pos[0]][i0].size(); j++) cout << (xPositions[pos[0]][i0][j] >> 6) << "," << (xPositions[pos[0]][i0][j] & 0b111111) << " ";
				// cout << "p: " << xParity[pos[0]][i0] << " " << endl;
				// for (uint32_t j = 0; j < xPositions[pos[1]][i1].size(); j++) cout << (xPositions[pos[1]][i1][j] >> 6) << "," << (xPositions[pos[1]][i1][j] & 0b111111) << " ";
				// cout << "p: " << xParity[pos[1]][i1] << " " << endl;
				// for (uint32_t j = 0; j < xPositions[pos[2]][i2].size(); j++) cout << (xPositions[pos[2]][i2][j] >> 6) << "," << (xPositions[pos[2]][i2][j] & 0b111111) << " ";
				// cout << "p: " << xParity[pos[2]][i2] << endl;
				vector<int> tmp = vectorXOR(constraint,xPositions[pos[0]][i0]);
				tmp = vectorXOR(tmp,xPositions[pos[1]][i1]);
				tmp = vectorXOR(tmp,xPositions[pos[2]][i2]);
				if (tmp.size() == 0) 
				{
					// cout << "got in here --------------------------" << endl;
					parity ^= xParity[pos[0]][i0] ^ xParity[pos[1]][i1] ^ xParity[pos[2]][i2];
					vector<int> tmp2 = constraint;
					tmp2.insert(tmp2.end(),xPositions[pos[0]][i0].begin(),xPositions[pos[0]][i0].end());
					tmp2.insert(tmp2.end(),xPositions[pos[1]][i1].begin(),xPositions[pos[1]][i1].end());
					tmp2.insert(tmp2.end(),xPositions[pos[2]][i2].begin(),xPositions[pos[2]][i2].end());
					insertKeyConstraint(tmp2,parity);
				}
			}
		}
	}
	return;
}

void searchFour(vector<int> pos, vector<int> constraint, int parity, vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity)
{
	for (uint32_t i0 = 0; i0 < xPositions[pos[0]].size(); i0++)
	{
		for (uint32_t i1 = 0; i1 < xPositions[pos[1]].size(); i1++)
		{
			for (uint32_t i2 = 0; i2 < xPositions[pos[2]].size(); i2++)
			{
				for (uint32_t i3 = 0; i3 < xPositions[pos[3]].size(); i3++)
				{
					// cout << "comparing: ";
					// for (uint32_t j = 0; j < constraint.size(); j++) cout << (constraint[j] >> 6) << "," << (constraint[j] & 0b111111) << " ";
					// cout << "p: " << parity << " ";
					// cout << endl << "xpos: ";
					// for (uint32_t j = 0; j < xPositions[pos[0]][i0].size(); j++) cout << (xPositions[pos[0]][i0][j] >> 6) << "," << (xPositions[pos[0]][i0][j] & 0b111111) << " ";
					// cout << "p: " << xParity[pos[0]][i0] << " " << endl;
					// for (uint32_t j = 0; j < xPositions[pos[1]][i1].size(); j++) cout << (xPositions[pos[1]][i1][j] >> 6) << "," << (xPositions[pos[1]][i1][j] & 0b111111) << " ";
					// cout << "p: " << xParity[pos[1]][i1] << " " << endl;
					// for (uint32_t j = 0; j < xPositions[pos[2]][i2].size(); j++) cout << (xPositions[pos[2]][i2][j] >> 6) << "," << (xPositions[pos[2]][i2][j] & 0b111111) << " ";
					// cout << "p: " << xParity[pos[2]][i2] << " " << endl;
					// for (uint32_t j = 0; j < xPositions[pos[3]][i3].size(); j++) cout << (xPositions[pos[3]][i3][j] >> 6) << "," << (xPositions[pos[3]][i3][j] & 0b111111) << " ";
					// cout << "p: " << xParity[pos[3]][i3] << endl;
					vector<int> tmp = vectorXOR(constraint,xPositions[pos[0]][i0]);
					tmp = vectorXOR(tmp,xPositions[pos[1]][i1]);
					tmp = vectorXOR(tmp,xPositions[pos[2]][i2]);
					tmp = vectorXOR(tmp,xPositions[pos[3]][i3]);
					if (tmp.size() == 0) 
					{
						// cout << "got in here --------------------------" << endl;
						parity ^= xParity[pos[0]][i0] ^ xParity[pos[1]][i1] ^ xParity[pos[2]][i2] ^ xParity[pos[3]][i3];
						vector<int> tmp2 = constraint;
						tmp2.insert(tmp2.end(),xPositions[pos[0]][i0].begin(),xPositions[pos[0]][i0].end());
						tmp2.insert(tmp2.end(),xPositions[pos[1]][i1].begin(),xPositions[pos[1]][i1].end());
						tmp2.insert(tmp2.end(),xPositions[pos[2]][i2].begin(),xPositions[pos[2]][i2].end());
						tmp2.insert(tmp2.end(),xPositions[pos[3]][i3].begin(),xPositions[pos[3]][i3].end());
						insertKeyConstraint(tmp2,parity);
					}
				}
			}
		}
	}
	return;
}
// here
void findLinearConstraints(vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity, 
								  vector<vector<vector<int>>> yPositions, vector<vector<int>> yParity, uint32_t nr)
{
	// cout << "The inputs: " << endl;
	// for (int i = 0; i < xPositions.size(); i++)
	// {
	// 	for (int j = 0; j < xPositions[i].size(); j++)
	// 	{
	// 		if (xPositions[i].size() == 0) continue;
	// 		for (int k = 0; k < xPositions[i][j].size(); k++) cout << (xPositions[i][j][k] >> 7) << "," << (xPositions[i][j][k] & 0b1111111) << " ";
	// 		cout << "p: " << xParity[i][j] << endl;
	// 	}
	// }
	// // cout << "The outputs: " << endl;
	// for (int i = 0; i < yPositions.size(); i++)
	// {
	// 	for (int j = 0; j < yPositions[i].size(); j++)
	// 	{
	// 		if (yPositions[i].size() == 0) continue;
	// 		for (int k = 0; k < yPositions[i][j].size(); k++) cout << (yPositions[i][j][k] >> 7) << "," << (yPositions[i][j][k] & 0b1111111) << " ";
	// 		cout << "p: " << yParity[i][j] << endl;
	// 	}
	// }
	// exhaustively searching Ycells 0,1,2,3 of round 0
	// for each 4 Sboxes
	for (int n = 0; n < nr-1; n++)
	{
		for (int jumpFourCells = 0; jumpFourCells < (4 * SIZE/4); jumpFourCells++)
		{
			// singular cell
			// only 1 Y cell is involved
			for (int cell0 = 4*jumpFourCells+0; cell0 < 4*jumpFourCells+4; cell0++){
				for (uint32_t i = 0; i < yPositions[n*32+cell0].size(); i++){
					vector<int> constraint = yPositions[n*32+cell0][i];
					int parity = yParity[n*32+cell0][i];
					vector<int> pos = checkXPos(constraint); // for each Sbox, there are at most 4 output bits, with the particular equation, find the Sboxes in the next round that is linked
					if (pos.size() == 0) { cout << "Something is wrong here!. Exiting..." << endl; exit(0);}
					else if (pos.size() == 1) searchOne(pos,constraint,parity,xPositions,xParity); // if there is only 1 Sbox involved for X and Y
					else if (pos.size() == 2) searchTwo(pos,constraint,parity,xPositions,xParity);  // if there is only 2 Sboxes involved for X and 1 Sbox for Y
					else if (pos.size() == 3) searchThree(pos,constraint,parity,xPositions,xParity);
					else if (pos.size() == 4) searchFour(pos,constraint,parity,xPositions,xParity);
				}
			}
			// double cell
			// 2 Y cells are involved
			for (int cell0 = 4*jumpFourCells+0; cell0 < 4*jumpFourCells+4; cell0++)
			{
				for (int cell1 = cell0+1; cell1 < 4*jumpFourCells+4; cell1++)
				{
					for (uint32_t i0 = 0; i0 < yPositions[n*32+cell0].size(); i0++)
					{
						for (uint32_t i1 = 0; i1 < yPositions[n*32+cell1].size(); i1++)
						{
							vector<int> constraint = yPositions[n*32+cell0][i0];
							constraint = vectorXOR(constraint,yPositions[n*32+cell1][i1]);
							int parity = yParity[n*32+cell0][i0] ^ yParity[n*32+cell1][i1];
							vector<int> pos = checkXPos(constraint);
							if (pos.size() == 0) { cout << "Something is wrong here!. Exiting..." << endl; exit(0);}
							else if (pos.size() == 1) searchOne(pos,constraint,parity,xPositions,xParity);
							else if (pos.size() == 2) searchTwo(pos,constraint,parity,xPositions,xParity);
							else if (pos.size() == 3) searchThree(pos,constraint,parity,xPositions,xParity);
							else if (pos.size() == 4) searchFour(pos,constraint,parity,xPositions,xParity);
						}
					}
				}
			}

			// triple

			for (int cell0 = 4*jumpFourCells+0; cell0 < 4*jumpFourCells+4; cell0++)
			{
				for (int cell1 = cell0+1; cell1 < 4*jumpFourCells+4; cell1++)
				{
					for (int cell2 = cell1+1; cell2 < 4*jumpFourCells+4; cell2++)
					{
						for (uint32_t i0 = 0; i0 < yPositions[n*32+cell0].size(); i0++)
						{
							for (uint32_t i1 = 0; i1 < yPositions[n*32+cell1].size(); i1++)
							{
								for (uint32_t i2 = 0; i2 < yPositions[n*32+cell2].size(); i2++)
								{
									vector<int> constraint = yPositions[n*32+cell0][i0];
									constraint = vectorXOR(constraint,yPositions[n*32+cell1][i1]);
									constraint = vectorXOR(constraint,yPositions[n*32+cell2][i2]);
									int parity = yParity[n*32+cell0][i0] ^ yParity[n*32+cell1][i1] ^ yParity[n*32+cell2][i2];
									vector<int> pos = checkXPos(constraint);
									if (pos.size() == 0) { cout << "Something is wrong here!. Exiting..." << endl; exit(0);}
									else if (pos.size() == 1) searchOne(pos,constraint,parity,xPositions,xParity);
									else if (pos.size() == 2) searchTwo(pos,constraint,parity,xPositions,xParity);
									else if (pos.size() == 3) searchThree(pos,constraint,parity,xPositions,xParity);
									else if (pos.size() == 4) searchFour(pos,constraint,parity,xPositions,xParity);
								}
							}
						}
					}
				}
			}
			// QUADRUPLE
			for (int cell0 = 4*jumpFourCells+0; cell0 < 4*jumpFourCells+4; cell0++)
			{
				for (int cell1 = cell0+1; cell1 < 4*jumpFourCells+4; cell1++)
				{
					for (int cell2 = cell1+1; cell2 < 4*jumpFourCells+4; cell2++)
					{
						for (int cell3 = cell2+1; cell3 < 4*jumpFourCells+4; cell3++)
						{
							for (uint32_t i0 = 0; i0 < yPositions[n*32+cell0].size(); i0++)
							{
								for (uint32_t i1 = 0; i1 < yPositions[n*32+cell1].size(); i1++)
								{
									for (uint32_t i2 = 0; i2 < yPositions[n*32+cell2].size(); i2++)
									{
										for (uint32_t i3 = 0; i3 < yPositions[n*32+cell3].size(); i3++)
										{
											vector<int> constraint = yPositions[n*32+cell0][i0];
											constraint = vectorXOR(constraint,yPositions[n*32+cell1][i1]);
											constraint = vectorXOR(constraint,yPositions[n*32+cell2][i2]);
											constraint = vectorXOR(constraint,yPositions[n*32+cell3][i3]);
											int parity = yParity[n*32+cell0][i0] ^ yParity[n*32+cell1][i1] ^ yParity[n*32+cell2][i2] ^ yParity[n*32+cell3][i3];
											vector<int> pos = checkXPos(constraint);
											if (pos.size() == 0) { cout << "Something is wrong here!. Exiting..." << endl; exit(0);}
											else if (pos.size() == 1) searchOne(pos,constraint,parity,xPositions,xParity);
											else if (pos.size() == 2) searchTwo(pos,constraint,parity,xPositions,xParity);
											else if (pos.size() == 3) searchThree(pos,constraint,parity,xPositions,xParity);
											else if (pos.size() == 4) searchFour(pos,constraint,parity,xPositions,xParity);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void AddConstants(vector<vector<vector<int>>> &yPositions, vector<vector<int>> &yParity, int nr)
{
	for (int n = 0; n < nr-1; n++)
	{
		for (int i = 0; i < 16; i++)
		{
			for (uint32_t j = 0; j < yPositions[n*32+i].size(); j++)
			{
				for (uint32_t k = 0; k < yPositions[n*32+i][j].size(); k++)
				{
					if ((yPositions[n*32+i][j][k] & 0b1111111) == 63) yParity[n*32+i][j] ^= 1;
					if ((yPositions[n*32+i][j][k] & 0b1111111) == 23) yParity[n*32+i][j] ^= ((getConstants(n) >> 5) & 0b1);
					if ((yPositions[n*32+i][j][k] & 0b1111111) == 19) yParity[n*32+i][j] ^= ((getConstants(n) >> 4) & 0b1);
					if ((yPositions[n*32+i][j][k] & 0b1111111) == 15) yParity[n*32+i][j] ^= ((getConstants(n) >> 3) & 0b1);
					if ((yPositions[n*32+i][j][k] & 0b1111111) == 11) yParity[n*32+i][j] ^= ((getConstants(n) >> 2) & 0b1);
					if ((yPositions[n*32+i][j][k] & 0b1111111) ==  7) yParity[n*32+i][j] ^= ((getConstants(n) >> 1) & 0b1);
					if ((yPositions[n*32+i][j][k] & 0b1111111) ==  3) yParity[n*32+i][j] ^= ((getConstants(n) >> 0) & 0b1);
				}
			}
		}
	}
}

vector<vector<int>> propagateLinear(int nr)
{
	uint64_t masterKey[30][128] = {0};
	for (int i = 0; i < 128; i++) masterKey[0][i] = 127-i; // initialization of the key
	for (int i = 1; i < 30; i++) // key schedule
	{
		for (int j = 0; j < 8; j++)
		{
			for (int k = 0; k < 16; k++)
			{
				if (j == 0) {masterKey[i][k] = masterKey[i-1][127 + (((14+k) % 16) - 31)];}
				if (j == 1) {masterKey[i][j * 16 + k] = masterKey[i-1][127 + (((4+k) % 16) - 15)];}
				if (j > 1) masterKey[i][j*16+k] = masterKey[i-1][((j-2)*16 + k)];
			}
		}
	} 
	vector<vector<int>> masterKeyConstraints;
	// cout << "Keys before keyschedule" << endl;

	// replacing the subkeys by their master key schedule
	#if SIZE == 4
	for (uint32_t i = 0; i < keyConstraints.size(); i++)
	{
		vector<int> tmp;
		for (uint32_t j = 0; j < keyConstraints[i].size()-1; j++)
		{
			int n = (keyConstraints[i][j] >> 7);
			int pos = keyConstraints[i][j] & 0b1111111;
			int keyIndex = pos / 4;
			int uvIndex = pos % 4; // u -->  1, v --> 0
			int index;
			if (uvIndex == 0) index = 127 - keyIndex;
			else if (uvIndex == 1) index = 127 - 16 - keyIndex;
			else { cout << "what?" << endl; exit(0);}
			tmp.push_back(masterKey[n][index]);
		}
		tmp.push_back(keyConstraints[i][keyConstraints[i].size()-1]);
		masterKeyConstraints.push_back(tmp);
	}

	#elif SIZE == 8
	for (uint32_t i = 0; i < keyConstraints.size(); i++)
	{
		vector<int> tmp;
		for (uint32_t j = 0; j < keyConstraints[i].size()-1; j++)
		{
			int n = (keyConstraints[i][j] >> 7);
			int pos = keyConstraints[i][j] & 0b1111111;
			int keyIndex = pos / 4; 
			int uvIndex = pos % 4; // u == 2, v == 1
			int index;
			if (uvIndex == 2) index = 63 - keyIndex;
			else if (uvIndex == 1) index = 127 - keyIndex;
			else {cout << "what?" << endl; exit(0);}
			tmp.push_back(masterKey[n][index]);
		}
		tmp.push_back(keyConstraints[i][keyConstraints[i].size()-1]); // adding the parity
		masterKeyConstraints.push_back(tmp);
	}
	#endif
	return masterKeyConstraints;
}

bool isZeroMatrix(uint8_t **M, int ROWSIZE, int COLSIZE)
{
	for (int r = 0; r < ROWSIZE; r++)
	{
		for (int c = 0; c < COLSIZE; c++)
		{
			if (M[r][c] != 0) {return false;}
		}
	}
	return true;
}


void computeDimension(vector<vector<int>> masterKeyConstraints)
{
	printout("masterKeyConstraints");
	for (int i = 0; i < masterKeyConstraints.size(); i++)
	{
		for (int j = 0; j < masterKeyConstraints[i].size(); j++)
		{
			cout << (masterKeyConstraints[i][j] >> 7) << "," << (masterKeyConstraints[i][j] & 0b1111111) << " ";
		}
		cout << endl;
	}
	cout << endl;
	if (masterKeyConstraints.size() == 0) 
	{
		cout << "Dimension: " << 128 << endl; 
		return;
	}
	
	// there is at most 128 possible key values + 1 RHS
	uint8_t **M; M = new uint8_t*[masterKeyConstraints.size()];
	// setting up the matrix
	for (int i = 0; i < masterKeyConstraints.size(); i++) M[i] = new uint8_t[129]();
	for (int i = 0; i < masterKeyConstraints.size(); i++)
	{
		for (int j = 0; j < masterKeyConstraints[i].size()-1; j++)
		{
			M[i][masterKeyConstraints[i][j] & 0b1111111] = 1;
		}
		M[i][128] = masterKeyConstraints[i][masterKeyConstraints[i].size()-1]; // constant
	}
	// perform Gaussian Elimination
	GaussianElimination(M,masterKeyConstraints.size(),129);
	for (int r = masterKeyConstraints.size()-1; r >= 0; r--)
	{
		for (int c = 0; c < 129; c++)
		{
			if ((M[r][c] == 1) && (c == 128)) {cout << "impossible due to multiple key constraints!" << endl; return;}
			else if ((M[r][c] == 1) && (c < 128)) // compute the dimension
			{
				cout << "Dimension: " << 127 - r << endl;
				return;
			}
		}
	}
	delete[] M;
	return;
}

void computeProb(uint8_t alpha[30][16], uint32_t nr)
{
	int prob = 0;
	uint32_t DDT[16][16] = {0};
	computeDDT(DDT);

	uint8_t state_temp_before[16];
	uint8_t state_temp_after[16];
	for (uint32_t n = 0; n < nr; n++)
	{
		for (int i = 0; i < 16; i++) state_temp_before[i] = alpha[n][i];
		for (int i = 0; i < 16; i++) state_temp_after[i] = alpha[n+1][i];
		#if SIZE == 4
		invPermutation(state_temp_after);
		for (int i = 0; i < 16; i++) prob += -log2(DDT[state_temp_before[i]][state_temp_after[i]]/16.0);
		#elif SIZE == 8
		invPermutation8(state_temp_after);
		for (int i = 0; i < 16; i++)
		{
			prob += -log2(DDT[state_temp_before[i]>>4][state_temp_after[i]>>4]/16.0);
			prob += -log2(DDT[state_temp_before[i]&0xf][state_temp_after[i]&0xf]/16.0);
		}
		#endif
	}
	cout << "Stated prob: " << prob << endl;
}

