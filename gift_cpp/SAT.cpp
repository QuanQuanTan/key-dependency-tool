#include <cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <vector>
#include "SAT.h"
#include "gift.h"


using namespace std;
using namespace CMSat;

uint32_t freeVar;
uint64_t noOfClauses;
int indexLUT[1<<16] = {0};
uint32_t sboxFreeVar[16] = {0};

/*
How the bits are being arranged in indexLUT:
// nr (5 bits) | bit index. Note bit 0 is rightmost bit (7 bits) | x wire or y wire: x = 0, y = 1
// 1 (key) | 00000 = nr (5 bits) | 0000000 = bit index (7 bits) | 0
*/
void addSbox(uint8_t input, uint8_t output, CMSat::SATSolver &solver, int sbox_index, int n, int size = 4)
{
	/*
	This function adds the restrictions given by a particular input/output pair into the solver
	if input, output == 3,5 
	Then we will add the following:
	not(t) \/ not(x0) 
	not(t) \/ not(x1) 
	not(t) \/ x2
	not(t) \/ x3

	not(t) \/ not(y0)
	not(t) \/ y1
	not(t) \/ not(y2)
	not(t) \/ y3
	if t is true then the remaining must be true
	t is stored in sboxFreeVar for a clause later (t0 \/ t1 \/ ... \/ t15) 
	*/
	vector<Lit> clause_Lit;
	int tmpIndex;
	if (size == 4)
	{
		for (int i = 0; i < 4; i++) // adding the effect of constants in for the x wires. Since by default, we include constants in the y wires of the ACTIVE sboxes only
		{
			if (getPerm(sbox_index * 4 + i) == 63) output ^= (1 << i);
			if (getPerm(sbox_index * 4 + i) == 23) output ^= (((getConstants(n) >> 5) & 1) << i);
			if (getPerm(sbox_index * 4 + i) == 19) output ^= (((getConstants(n) >> 4) & 1) << i);
			if (getPerm(sbox_index * 4 + i) == 15) output ^= (((getConstants(n) >> 3) & 1) << i);
			if (getPerm(sbox_index * 4 + i) == 11) output ^= (((getConstants(n) >> 2) & 1) << i);
			if (getPerm(sbox_index * 4 + i) ==  7) output ^= (((getConstants(n) >> 1) & 1) << i);
			if (getPerm(sbox_index * 4 + i) ==  3) output ^= (((getConstants(n) >> 0) & 1) << i);
		}
		sboxFreeVar[input] = freeVar++; // storing the ``t" for a clause later
		// adding the input
		for (int i = 0; i < 4; i++)
		{
			clause_Lit.clear(); 
			clause_Lit.push_back(Lit(sboxFreeVar[input], true)); // true = adding a ``not"
			tmpIndex = ((n-1) << 8) + (((4 * sbox_index + (3-i)) << 1) + 1); // +1 is given here because we are looking at y of the prev round
			// checking if this bit has been referenced in indexLUT. If not, add that in
			if (indexLUT[tmpIndex] == -1) indexLUT[tmpIndex] = freeVar++; 

			// adding the index into the clause
			clause_Lit.push_back(Lit(indexLUT[tmpIndex], (((input >> (3-i)) & 1) + 1) % 2)); 
			solver.add_clause(clause_Lit); noOfClauses++; 
		}
		// adding the output
		for (int i = 0; i < 4; i++)
		{	
			clause_Lit.clear(); 
			clause_Lit.push_back(Lit(sboxFreeVar[input], true));
			tmpIndex = (n << 8) + (getPerm(4 * sbox_index + (3-i)) << 1); // no +1 is given here because we are looking x of the next round
			// checking if this bit has been referenced in indexLUT. If not, add that in
			if (indexLUT[tmpIndex] == -1) indexLUT[tmpIndex] = freeVar++;
			// adding the index into the clause
			clause_Lit.push_back(Lit(indexLUT[tmpIndex], (((output >> (3-i)) & 1) + 1) % 2));
			solver.add_clause(clause_Lit);
			noOfClauses++;
		}
	}
	else if (size == 8)
	{
		for (int i = 0; i < 4; i++) // adding the effect of constants in for the x wires. Since by default, we include constants in the y wires of the ACTIVE sboxes only
		{
			if (getPerm8(sbox_index * 4 + i) == 127) output ^= (1 << i);
			if (getPerm8(sbox_index * 4 + i) == 23) output ^= (((getConstants(n) >> 5) & 1) << i);
			if (getPerm8(sbox_index * 4 + i) == 19) output ^= (((getConstants(n) >> 4) & 1) << i);
			if (getPerm8(sbox_index * 4 + i) == 15) output ^= (((getConstants(n) >> 3) & 1) << i);
			if (getPerm8(sbox_index * 4 + i) == 11) output ^= (((getConstants(n) >> 2) & 1) << i);
			if (getPerm8(sbox_index * 4 + i) ==  7) output ^= (((getConstants(n) >> 1) & 1) << i);
			if (getPerm8(sbox_index * 4 + i) ==  3) output ^= (((getConstants(n) >> 0) & 1) << i);
		}	
		sboxFreeVar[input] = freeVar++; // storing the ``t" for a clause later
		// adding the input
		for (int i = 0; i < 4; i++)
		{
			clause_Lit.clear(); 
			clause_Lit.push_back(Lit(sboxFreeVar[input], true)); // true = adding a ``not"
			tmpIndex = ((n-1) << 8) + (((4 * sbox_index + (3-i)) << 1) + 1); // +1 is given here because we are looking at y of the prev round
			// checking if this bit has been referenced in indexLUT. If not, add that in
			if (indexLUT[tmpIndex] == -1) indexLUT[tmpIndex] = freeVar++; 

			// adding the index into the clause
			clause_Lit.push_back(Lit(indexLUT[tmpIndex], (((input >> (3-i)) & 1) + 1) % 2)); 
			solver.add_clause(clause_Lit); noOfClauses++; 
		}
		// adding the output
		for (int i = 0; i < 4; i++)
		{	
			clause_Lit.clear(); 
			clause_Lit.push_back(Lit(sboxFreeVar[input], true));
			tmpIndex = (n << 8) + (getPerm8(4 * sbox_index + (3-i)) << 1); // no +1 is given here because we are looking x of the next round
			// checking if this bit has been referenced in indexLUT. If not, add that in
			if (indexLUT[tmpIndex] == -1) indexLUT[tmpIndex] = freeVar++;
			// adding the index into the clause
			clause_Lit.push_back(Lit(indexLUT[tmpIndex], (((output >> (3-i)) & 1) + 1) % 2));
			solver.add_clause(clause_Lit);
			noOfClauses++;
		}
	}
}

void SATSolver(uint8_t alpha[30][16], vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity, vector<vector<vector<int>>> yPositions, vector<vector<int>> yParity, int nr, int size)
{
	noOfClauses = 0;
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
	// do it in batches of 4 rounds (GIFT-64 key rotation)
	/*
	For GIFT-64, this is what we do:
	1. split the rounds into 4 different batches: 4n, 4n+1, 4n+2, 4n+3
	2. (Reason) 4n+1 involves inactive Sboxes from round 1,5,9,...
	3. The inactive Sboxes involve keys from k0,k1 (round 0) and k2,k3 (round 1)
	4. According to the key schedule, the same 16 * 2 are used once again.
	5. Once we have retrieved them, we use cryptominisat to find out if they are compatible

	For GIFT-128, doing similar techniques with GIFT-64 will lead to a memory explosion. Instead, we simply analyze each round independently
	*/
	if (size == 4)
	{
		for (int batch = 0; batch < 4; batch++)
		{
			cout << "Performing batch number: " << batch << endl;

			// temp storage to sieve out the shortlisted batch
			vector<vector<vector<int>>> xPositions_temp;
			vector<vector<int>> xParity_temp;
			vector<vector<vector<int>>> yPositions_temp;
			vector<vector<int>> yParity_temp;
			int tmpIndex;
			// solver related declarations
			CMSat::SATSolver solver;
		    vector<uint32_t> clause;
		    vector<Lit> clause_Lit;
		    solver.set_num_threads(8);
		    solver.new_vars((1 << 16));
		    freeVar = 0; // this is the next free variable available to be used
		    for (int i = 0; i < (1 << 16); i++) indexLUT[i] = -1; // by default, we do not have any mapping to freeVar yet

			// shortlisting the batch
			for (int n = 0; n < nr-1; n++)
			{
				if (((n % 4) == ((batch + 0) % 4)) || ((n % 4) == ((batch + 1) % 4)))
				{
					for (int c = 0; c < 32; c++)
					{
						xPositions_temp.push_back(xPositions[n*32+c]);
						yPositions_temp.push_back(yPositions[n*32+c]);
						xParity_temp.push_back(xParity[n*32+c]);
						yParity_temp.push_back(yParity[n*32+c]);
					}
				}
			}
		    
		    // add the x (xor) equations into the solver
		    for (uint32_t i = 0; i < xPositions_temp.size(); i++)
		    {
		    	for (uint32_t j = 0; j < xPositions_temp[i].size(); j++)
		    	{
		    		clause.clear();
		    		for (uint32_t k = 0; k < xPositions_temp[i][j].size(); k++)
		    		{
		    			if (indexLUT[xPositions_temp[i][j][k] << 1] == -1) indexLUT[xPositions_temp[i][j][k] << 1] = freeVar++;
		    			clause.push_back(indexLUT[xPositions_temp[i][j][k] << 1]); // xPositions have a 0 in the LSB
		    		}
		    		solver.add_xor_clause(clause,bool(xParity_temp[i][j]));
		    		noOfClauses++;
		    	}
		    }
		    // add the y (xor) equations into the solver
		    for (uint32_t i = 0; i < yPositions_temp.size(); i++)
		    {
		    	for (uint32_t j = 0; j < yPositions_temp[i].size(); j++)
		    	{
		    		clause.clear();
		    		for (uint32_t k = 0; k < yPositions_temp[i][j].size(); k++)
		    		{
		    			// cout << (yPositions_temp[i][j][k] >> 7) << "," << (yPositions_temp[i][j][k] % 128) << " ";
		    			if (indexLUT[(yPositions_temp[i][j][k] << 1)+1] == -1) indexLUT[(yPositions_temp[i][j][k] << 1)+1] = freeVar++;
		    			clause.push_back(indexLUT[(yPositions_temp[i][j][k] << 1)+1]); // yPositions have a 1 in the LSB
		    		}
		    		// cout << "p: " << yParity_temp[i][j] << endl;
		    		solver.add_xor_clause(clause,bool(yParity_temp[i][j]));
		    		noOfClauses++;
		    	}
		    }

		    // adding in the inactive Sboxes
		    for (int n = 1; n < nr-1; n++)
		    {
		    	if (((n % 4) == ((batch + 0) % 4)) || ((n % 4) == ((batch + 1) % 4))) // they must be from the correct batch
		    	{
			    	for (int c = 0; c < 16; c++)
			    	{
			    		if (alpha[n][c] == 0) // if it's inactive, add them in
			    		{
							int sbox = 15-c;
							for (int i = 0; i < 16; i++) addSbox(i,getSbox(i),solver,sbox,n); // 16 possible sbox transitions
							clause_Lit.clear();
							for (int i = 0; i < 16; i++) clause_Lit.push_back(Lit(sboxFreeVar[i], false)); // at least one freeVar must be 1 to ensure one Sbox transition
							solver.add_clause(clause_Lit);
							noOfClauses++;
			    		}
			    	}
		    	}
		    }

		    // establishing the link from previous round to this round and from this round to the next round with wires/ keys
		    for (int n = 0; n < nr-1; n++)
		    {
		    	if (((n % 4) == ((batch + 0) % 4)) || ((n % 4) == ((batch + 1) % 4)))
		    	{
			    	for (int c = 0; c < 64; c++) // x + y ( + k) = 0 
			    	{
			    		clause.clear();
			    		tmpIndex = (n << 8) + (c << 1);
			    		if (indexLUT[tmpIndex] == -1) { indexLUT[tmpIndex] = freeVar++;}
			    		if (indexLUT[tmpIndex + 1] == -1) { indexLUT[tmpIndex + 1] = freeVar++;}
			    		clause.push_back(indexLUT[tmpIndex]);
			    		clause.push_back(indexLUT[tmpIndex + 1]);
			    		if (c % 4 == 0) // key involved (v)
			    		{
			    			tmpIndex = (1 << 13) + masterKey[n][127-(c/4)];
			    			if (indexLUT[tmpIndex] == -1) {indexLUT[tmpIndex] = freeVar++;}
			    			clause.push_back(indexLUT[tmpIndex]);
			    		}
			    		if (c % 4 == 1) // key involved (u)
			    		{
			    			tmpIndex = (1 << 13) + masterKey[n][127-(16+c/4)];
			    			if (indexLUT[tmpIndex] == -1) {indexLUT[tmpIndex] = freeVar++;}
			    			clause.push_back(indexLUT[tmpIndex]);
			    		}
			    		solver.add_xor_clause(clause,false);
			    		noOfClauses++;
			    	}
			    }
		    }
		    cout << "no of clauses: " << noOfClauses << endl;
			CMSat::lbool ret1 = solver.solve();
			cout << "Is it possible? : ";
		    cout << ret1 << endl;
		}
	}
	else if (size == 8)
	{
		for (int n = 0; n < nr-2; n++)
		{
			// temp storage to sieve out the shortlisted batch
			vector<vector<vector<int>>> xPositions_temp;
			vector<vector<int>> xParity_temp;
			vector<vector<vector<int>>> yPositions_temp;
			vector<vector<int>> yParity_temp;
			int tmpIndex;
			// solver related declarations
			CMSat::SATSolver solver;
		    vector<uint32_t> clause;
		    vector<Lit> clause_Lit;
		    solver.set_num_threads(8);
		    solver.new_vars((1 << 16));
		    freeVar = 0; // this is the next free variable available to be used
		    for (int i = 0; i < (1 << 16); i++) indexLUT[i] = -1; // by default, we do not have any mapping to freeVar yet

		    for (int i = 0; i < 2; i++)
		    {
		    	for (int c = 0; c < 32; c++)
				{
					xPositions_temp.push_back(xPositions[(n+i)*32+c]);
					yPositions_temp.push_back(yPositions[(n+i)*32+c]);
					xParity_temp.push_back(xParity[(n+i)*32+c]);
					yParity_temp.push_back(yParity[(n+i)*32+c]);
				}
		    }

		    // add the x (xor) equations into the solver
		    for (uint32_t i = 0; i < xPositions_temp.size(); i++)
		    {
		    	for (uint32_t j = 0; j < xPositions_temp[i].size(); j++)
		    	{
		    		clause.clear();
		    		for (uint32_t k = 0; k < xPositions_temp[i][j].size(); k++)
		    		{
		    			if (indexLUT[xPositions_temp[i][j][k] << 1] == -1) indexLUT[xPositions_temp[i][j][k] << 1] = freeVar++;
		    			clause.push_back(indexLUT[xPositions_temp[i][j][k] << 1]); // xPositions have a 0 in the LSB
		    		}
		    		solver.add_xor_clause(clause,bool(xParity_temp[i][j]));
		    		noOfClauses++;
		    	}
		    }
		    // add the y (xor) equations into the solver
		    for (uint32_t i = 0; i < yPositions_temp.size(); i++)
		    {
		    	for (uint32_t j = 0; j < yPositions_temp[i].size(); j++)
		    	{
		    		clause.clear();
		    		for (uint32_t k = 0; k < yPositions_temp[i][j].size(); k++)
		    		{
		    			if (indexLUT[(yPositions_temp[i][j][k] << 1)+1] == -1) indexLUT[(yPositions_temp[i][j][k] << 1)+1] = freeVar++;
		    			clause.push_back(indexLUT[(yPositions_temp[i][j][k] << 1)+1]); // yPositions have a 1 in the LSB
		    		}
		    		solver.add_xor_clause(clause,bool(yParity_temp[i][j]));
		    		noOfClauses++;
		    	}
		    }
		    
		    // adding in the inactive Sboxes. In GIFT128, we split into first 4 bits and the last 4 bits
	    	for (int c = 0; c < 16; c++)
	    	{
	    		if ((alpha[n+1][c] >> 4) == 0) // if it's inactive
	    		{
					int sbox = 31-2*c;
					for (int i = 0; i < 16; i++) addSbox(i,getSbox(i),solver,sbox,n+1,8); // 16 possible sbox transitions
					clause_Lit.clear();
					for (int i = 0; i < 16; i++) clause_Lit.push_back(Lit(sboxFreeVar[i], false)); // at least one freeVar must be 1 to ensure one Sbox transition
					solver.add_clause(clause_Lit);
					noOfClauses++;
	    		}
	    		if ((alpha[n+1][c] % 16) == 0) // if it's inactive
	    		{
					int sbox = 31-2*c-1;
					for (int i = 0; i < 16; i++) addSbox(i,getSbox(i),solver,sbox,n+1,8); // 16 possible sbox transitions
					clause_Lit.clear();
					for (int i = 0; i < 16; i++) clause_Lit.push_back(Lit(sboxFreeVar[i], false)); // at least one freeVar must be 1 to ensure one Sbox transition
					solver.add_clause(clause_Lit);
					noOfClauses++;
	    		}
	    	}

	    	// establishing the link from previous round to this round and from this round to the next round with wires/ keys
		    for (int i = 0; i < 2; i++)
		    {
				for (int c = 0; c < 128; c++)
		    	{
		    		clause.clear();
		    		tmpIndex = ((n+i) << 8) + (c << 1);
		    		if (indexLUT[tmpIndex] == -1) { indexLUT[tmpIndex] = freeVar++;}
		    		if (indexLUT[tmpIndex + 1] == -1) { indexLUT[tmpIndex + 1] = freeVar++;}
		    		clause.push_back(indexLUT[tmpIndex]);
		    		clause.push_back(indexLUT[tmpIndex + 1]);
		    		if (c % 4 == 2) // key involved (u --> k5 || k4)
		    		{
		    			tmpIndex = (1 << 13) + masterKey[n+i][127-(64+c/4)];
		    			if (indexLUT[tmpIndex] == -1) {indexLUT[tmpIndex] = freeVar++;}
		    			clause.push_back(indexLUT[(1 << 13) + masterKey[n+i][127-(64+c/4)]]);
		    		}
		    		if (c % 4 == 1) // key involved (u --> k1 || k0)
		    		{
		    			tmpIndex = (1 << 13) + masterKey[n+i][127-(c/4)];
		    			if (indexLUT[tmpIndex] == -1) {indexLUT[tmpIndex] = freeVar++;}
		    			clause.push_back(indexLUT[tmpIndex]);
		    		}
		    		solver.add_xor_clause(clause,false);
		    		noOfClauses++;
		    	}
		    }
		    cout << "no of clauses: " << noOfClauses << endl;
			CMSat::lbool ret1 = solver.solve();
			cout << "Is it possible? : ";
		    cout << ret1 << endl;
		}
	}
}