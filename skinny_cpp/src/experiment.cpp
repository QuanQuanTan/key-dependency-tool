#include <iostream>
#include <vector>
#include <map>
#include <random>
#include <cmath>
#include <ctime>
#include "skinny.h"
#include "keyDependent.h"
#include "trails.h"
#include "printDetails.h"
#include <fstream>
#include <thread>
#include <mutex>

using namespace std;
random_device rd;
mt19937 rand_generator(rd());
uniform_int_distribution<uint16_t> rand16(0,0xffff);
uniform_int_distribution<uint64_t> rand64(0,0xffffffffffffffff);
#if SIZE == 4
uniform_int_distribution<uint8_t> rand4(0,0xf);
#elif SIZE == 8
uniform_int_distribution<uint8_t> rand8(0,0xff);	
#endif

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl

uint64_t possibleKeys = 0;
uint64_t totalKeys = 0;
vector<double> counts;

mutex term_mutex;
mutex setUp_mutex;



// for each node, we use the distribution given over here to compute
// for each successful right pair, we save the values into this



struct Node
{
	uint32_t index, pos;
	uint8_t outputdiff, inputdiff;
	uint8_t value0, value1;
	vector<Node*> prev;
	vector<Node*> next;
	int keyIndex = -1;
	uint8_t constant = 0;
	uint64_t possibleValues[1<<SIZE] = {0};
	uint64_t possibleValuesSum = 0;
	std::discrete_distribution<> dist;
	vector<uint32_t> possibleValuesVec;
};

vector<uint32_t> outputs, inputs, intermediates;
void printNode(Node* n){
	printout(int(n->index)); cout << endl;
	printout(int(n->inputdiff)); cout << endl;
	printout(int(n->outputdiff)); cout << endl;
	printout(int(n->constant)); cout << endl;
	for (int i = 0; i < (1<<SIZE); i++) cout << int(n->possibleValues[i]) << " ";
	cout << endl;		 
}


vector<Node*> getPrev(uint32_t element, Node* elementAddress, map<uint32_t,Node*> &constraintsAddress, vector<uint32_t> &keys, map<uint32_t,uint64_t*> &localSavedValues, bool force = false)
{
	uint32_t n = element >> (4+2*SIZE);
	uint32_t pos = (element >> (2*SIZE)) & 0xf;
	int inversePos[3] = {inverseMC[pos][0],inverseMC[pos][1],inverseMC[pos][2]};
	vector<Node*> prev;
	for (int i = 0; i < 3; i++)
	{
		if (inversePos[i] == -1) break;
		int npos = ((n-1) << 4) + inversePos[i];
		bool contains = false;
		for (auto& elem : constraintsAddress)
		{
		    if ((elem.first >> (2*SIZE)) == npos) 
		    {
		    	
		    	prev.push_back(elem.second);
		    	constraintsAddress[elem.first]->next.push_back(elementAddress);
		    	contains = true;
		    	break;
		    }
		}
		if (!contains && force){
			// this should be triggered when reach a higher order linear constraint
			Node* nn = new Node();
			nn->pos = inversePos[i];
			nn->inputdiff = diffStates[n-1][0][inversePos[i]];
			nn->outputdiff = diffStates[n-1][1][inversePos[i]];
			nn->index = ((n-1) << (4+2*SIZE)) + (nn->pos << (2*SIZE)) + (nn->inputdiff << SIZE) + (nn->outputdiff);
			// cout << hex << "force creating " << nn->index << " by element " << element << endl;
			
			auto it = localSavedValues.find(nn->index);
			if (it == localSavedValues.end()) localSavedValues[nn->index] = new uint64_t[1<<SIZE](); 

			nn->prev = getPrev(nn->index, nn, constraintsAddress,keys,localSavedValues); // recursively add
			nn->next.push_back(elementAddress); 

			// update constants
			uint32_t constant = getConstants(n-1);
			if (pos == 0) nn->constant = constant & 0xf;
			else if (pos == 4) nn->constant = (constant >> 4) & 0b11;
			else if (pos == 8) nn->constant = 0x2;

			// update possible Values
			uint32_t values[1<<SIZE] = {0};
			getYDDT(nn->inputdiff,nn->outputdiff,values);
			for (int k = 0; k < (1<<SIZE); k++)
			{
				nn->possibleValues[k] += values[k];
				nn->possibleValuesSum += values[k];
				for (int l = 0; l < values[k]; l++){
					nn->possibleValuesVec.push_back(k);
				}
			}

			// update key information
			if (pos < 8) {
				keys.push_back(nn->index);
			}
			constraintsAddress[nn->index] = nn;
			prev.push_back(nn);
		}
	}

	return prev;
}

vector<uint32_t> getKeys(vector<uint32_t> constraint)
{
	vector<uint32_t> keys;
	uint32_t maxN = findMaxN(constraint);
	for (int j = 0; j < constraint.size(); j++)
	{
		if ((constraint[j] >> (4+2*SIZE)) == maxN) continue;
		int pos = (constraint[j] >> (2*SIZE)) & 0xf;
		if (pos < 8) keys.push_back(constraint[j]);
	}

	return keys;
}

uint32_t getKeyDiff(uint32_t element, const uint8_t key_diff[4][4][4])
{
 	uint8_t key_diff_temp[4][4][4];
	for (int i = 0; i < 4; i++){for (int j = 0; j < 4; j++){for (int k = 0; k < 4; k++) key_diff_temp[i][j][k] = key_diff[i][j][k];}}
	int nr = element >> (4 + 2*SIZE);
	int n = 0;
	while (n < nr)
	{
		#if SIZE == 4
		key_schedule_round_function(key_diff_temp);
		#elif SIZE == 8
		key_schedule_round_function8(key_diff_temp);
		#endif
		n++;
	}
	int pos = (element >> (2*SIZE)) & 0xf;
	uint32_t keyDiff = (key_diff_temp[0][pos/4][pos%4] << (2*SIZE)) ^ (key_diff_temp[1][pos/4][pos%4] << (SIZE)) ^ (key_diff_temp[2][pos/4][pos%4]);
	return keyDiff;
}

vector<uint32_t> intersection(vector<uint32_t> v1, vector<uint32_t> v2){
    vector<uint32_t> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
    return v3;
}

void getPossibleKeys(vector<vector<uint32_t>> constraints, map<uint32_t,uint32_t> keyAddress, map<uint32_t,vector<uint32_t>> &keyPossible){
	for (int i = 0; i < constraints.size(); i++){
		// we only proceed with this constraint if 1. linear, 2. (TK_NUM = 0/1) OR (TK_NUM = 2/3 AND it is the first of it's key in terms of position)
		if ((constraints[i][constraints[i].size()-2] >> (4 + 2 * SIZE)) == (constraints[i][constraints[i].size()-1] >> (4 + 2 * SIZE))) continue;
		if (!isLinear(constraints[i])) continue;
		bool cont = false;
		uint32_t key = getKeys(constraints[i])[0];
		if (TK_NUM > 1){
			if ((keyAddress[key] & 0xf) == 0) cont = true;
		}
		else cont = true;
		if (!cont) continue;
		uint32_t keyPos = (keyAddress[key] >> 4) << 4;;
		vector<vector<uint32_t>> values;
		for (int j = 0; j < constraints[i].size()-1; j++){
			vector<uint32_t> tmp;
			uint32_t* yddt; yddt = new uint32_t[1<<SIZE]();
			getYDDT((constraints[i][j] >> SIZE) & ANDval, constraints[i][j] & ANDval, yddt);
			addRoundConstant(constraints[i][j],yddt);
			for (int k = 0; k < (1<<SIZE); k++){
				if (yddt[k] > 0) tmp.push_back(k);
			} 
			values.push_back(tmp);
		}
		uint32_t* xddt; xddt = new uint32_t[1<<SIZE]();
		getXDDT((constraints[i][constraints[i].size()-1] >> SIZE) & ANDval, constraints[i][constraints[i].size()-1] & ANDval, xddt);
		// now, we find the possible key values
		uint64_t keys[1<<SIZE] = {0};
		int indicesMax[constraints[i].size()-1] = {0};
		for (int j = 0; j < constraints[i].size()-1; j++) indicesMax[j] = values[j].size();
		for (int k = 0; k < (1<<SIZE); k++){
			int indices[constraints[i].size()-1] = {0};
			
			while (true){
				bool breakFlag = false;
				uint8_t v = k;
				for (int l = 0; l < constraints[i].size()-1; l++) v ^= values[l][indices[l]];
				if (xddt[v] == 1) {
					keys[k]++;
					breakFlag = true;
				}

				// increment of index
				int indexToMove = constraints[i].size()-2;
				while (indices[indexToMove] == indicesMax[indexToMove] - 1){
					indices[indexToMove] = 0;
					indexToMove--;
					if (indexToMove == -1){
						breakFlag = true;
						break;
					}
				}
				if (breakFlag) break;
				indices[indexToMove]++;
			}
		}
		// putting the keys into a vector
		vector<uint32_t> keyVec;
		for (int j = 0; j < (1<<SIZE); j++){
			if (keys[j] > 0) keyVec.push_back(j);
		}
		// if keyPossible is empty, we put it in, otherwise we have to find the intersection
		auto it = keyPossible.find(keyPos);
		if (it == keyPossible.end()){
			keyPossible[keyPos] = keyVec;
		}
		else{
			keyPossible[keyPos] = intersection(keyVec,keyPossible[keyPos]);
		}
	}
}



void organizeKeys(vector<uint32_t> keys, const uint8_t key_diff[4][4][4], vector<uint32_t> &keyDiff, map<uint32_t,uint32_t> &keyAddress, vector<vector<uint8_t>> &KEYS0,vector<vector<uint8_t>> &KEYS1)
{
	vector<vector<uint32_t>> KEYS;
	// get the set of keys involved
	// it should be a data structure that is easy to update using LFSRs
	sort(keys.begin(),keys.end());
	int index[16] = {0};
	for (int i = 0; i < 16; i++) index[i] = -1;
	int maxN = findMaxN(keys);
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
			keyDiff.push_back(getKeyDiff(keys[i],key_diff));
			index[pos] = KEYS.size();
			keyAddress[keys[i]] = (KEYS.size() << 4);
			vector<uint32_t> key = {keys[i]};
			KEYS.push_back(key);
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
				keyAddress[keys[i]] = (index[pos] << 4) + KEYS[index[pos]].size();
				KEYS[index[pos]].push_back(keys[i]);
			}
		}
	}

	// replicating to form KEY0 and KEY1
	for (int i = 0; i < KEYS.size(); i++)
	{
		vector<uint8_t> tmp0;
		for (int j = 0; j < KEYS[i].size(); j++) tmp0.push_back(0);
		KEYS0.push_back(tmp0);
	}
	KEYS1 = KEYS0;
}


void setUp(vector<vector<uint32_t>> constraints, const uint8_t key_diff[4][4][4], map<uint32_t,Node*> &constraintsAddress, map<uint32_t,uint32_t> &keyAddress, vector<uint32_t> &keyDiff, vector<uint32_t> &keys, vector<vector<uint8_t>>&KEYS0, vector<vector<uint8_t>>&KEYS1, map<uint32_t,uint64_t*> &localSavedValues, map<uint32_t,vector<uint32_t>> &keyPossible)
{	
	// insert the input vectors
	for (int i = 0; i < constraints.size(); i++)
	{
		uint32_t maxN = findMaxN(constraints[i]);
		// split the constraint into the various vectors
		// put the inputs in first
		for (int j = 0; j < constraints[i].size(); j++)
		{
			if (((constraints[i][j] >> (4+2*SIZE)) < maxN) && ((constraints[i][j] & ANDval) > 0))
			{
				auto it = constraintsAddress.find(constraints[i][j]);
				if (it != constraintsAddress.end()) continue;
				Node* n = new Node(); 
				// create a global memory of the node (meant for recording the distribution of values)
				// if it hasn't been created before
				
				localSavedValues[constraints[i][j]] = new uint64_t[1<<SIZE](); 
				n->index = constraints[i][j];

				// update constants
				int nr = constraints[i][j] >> (4+2*SIZE);
				int pos = (constraints[i][j] >> (2*SIZE)) & 0xf;
				uint32_t constant = getConstants(nr);
				if (pos == 0) n->constant = constant & 0xf;
				else if (pos == 4) n->constant = (constant >> 4) & 0b11;
				else if (pos == 8) n->constant = 0x2;

				// update pos
				n->pos = pos;
				// update input/output differences
				n->inputdiff = (constraints[i][j] >> SIZE) & ANDval;
				n->outputdiff = constraints[i][j] & ANDval;

				// if the cell is not used before, we use the DDT to determine the possibleValues
				// else, we use the record that we have
				auto it2 = savedValues.find(constraints[i][j]);
				if (it2 == savedValues.end())
				{
					// update possible Values
					uint32_t values[1<<SIZE] = {0};
					getYDDT(n->inputdiff,n->outputdiff,values);

					for (int k = 0; k < (1<<SIZE); k++)
					{
						n->possibleValues[k] += values[k];
						n->possibleValuesSum += values[k];
						for (int l = 0; l < values[k]; l++) n->possibleValuesVec.push_back(k);
					}
				}
				else if (it2 != savedValues.end())
				{
					// use the values from savedValues
					for (int k = 0; k < (1<<SIZE); k++)
					{
						n->possibleValues[k] += savedValues[constraints[i][j]][k];
						n->possibleValuesSum += savedValues[constraints[i][j]][k];
						for (int l = 0; l < savedValues[constraints[i][j]][k]; l++) n->possibleValuesVec.push_back(k);
					}
				}
				// update key information
				if (pos < 8) {
					keys.push_back(constraints[i][j]);
				}

				constraintsAddress[constraints[i][j]] = n;
			}
		}
	}
	// insert the intermediates
	for (int i = 0; i < constraints.size(); i++)
	{
		for (int j = 0; j < constraints[i].size(); j++)
		{
			if ((constraints[i][j] & ANDval) == 0)
			{
				Node* n = new Node();
				n->index = constraints[i][j];
				n->prev = getPrev(constraints[i][j], n, constraintsAddress,keys, localSavedValues);
				auto it = localSavedValues.find(constraints[i][j]);
				if (it == localSavedValues.end()) localSavedValues[constraints[i][j]] = new uint64_t[1<<SIZE](); 

				// update constants
				int nr = constraints[i][j] >> (4+2*SIZE);
				int pos = (constraints[i][j] >> (2*SIZE)) & 0xf;
				uint32_t constant = getConstants(nr);
				if (pos == 0) n->constant = constant & 0xf;
				else if (pos == 4) n->constant = (constant >> 4) & 0b11;
				else if (pos == 8) n->constant = 0x2;
				// update pos
				n->pos = pos;
				// update input/output differences
				n->inputdiff = (constraints[i][j] >> SIZE) & ANDval;
				n->outputdiff = constraints[i][j] & ANDval;

				auto it2 = savedValues.find(constraints[i][j]);
				if (it2 == savedValues.end())
				{
					// update possible Values
					uint32_t values[1<<SIZE] = {0};
					getYDDT(n->inputdiff,n->outputdiff,values);

					for (int k = 0; k < (1<<SIZE); k++)
					{
						n->possibleValues[k] = values[k];
						n->possibleValuesSum += values[k];
						for (int l = 0; l < values[k]; l++) n->possibleValuesVec.push_back(k);
					}
				}
				else if (it2 != savedValues.end())
				{
					// use the values from savedValues
					for (int k = 0; k < (1<<SIZE); k++)
					{
						for (int l = 0; l < savedValues[constraints[i][j]][k]; l++)
						{
							n->possibleValues[k] = savedValues[constraints[i][j]][k];
							n->possibleValuesSum += savedValues[constraints[i][j]][k];
							for (int l = 0; l < savedValues[constraints[i][j]][k]; l++) n->possibleValuesVec.push_back(k);
						}
					}
				}

				// update key information
				if (pos < 8) {
					keys.push_back(constraints[i][j]);
				}
				constraintsAddress[constraints[i][j]] = n;
				
			}
		}
	}
	// insert the outputs
	for (int i = 0; i < constraints.size(); i++)
	{
		for (int j = 0; j < constraints[i].size(); j++)
		{
			uint32_t maxN = findMaxN(constraints[i]);
			if ((constraints[i][j] >> (4+2*SIZE)) == maxN)
			{
				// check if constraintAddress contains it already (it might be that the output of one constraint is involved in the input of the other)
				// when that is the case, we add connections to it without creating new nodes
				auto it = constraintsAddress.find(constraints[i][j]);
				
				if (it == constraintsAddress.end())
				{
					auto it = localSavedValues.find(constraints[i][j]);
					if (it == localSavedValues.end()) localSavedValues[constraints[i][j]] = new uint64_t[1<<SIZE](); 

					Node* n = new Node();
					// update constants
					int nr = constraints[i][j] >> (4+2*SIZE);
					int pos = (constraints[i][j] >> (2*SIZE)) & 0xf;
					uint32_t constant = getConstants(nr);
					if (pos == 0) n->constant = constant & 0xf;
					else if (pos == 4) n->constant = (constant >> 4) & 0b11;
					else if (pos == 8) n->constant = 0x2;
					// update pos
					n->pos = pos;
					// update input/output differences
					n->inputdiff = (constraints[i][j] >> SIZE) & ANDval;
					n->outputdiff = constraints[i][j] & ANDval;
					n->index = constraints[i][j];
					auto it2 = savedValues.find(constraints[i][j]);
					if (it2 == savedValues.end())
					{
						// update possible Values
						uint32_t values[1<<SIZE] = {0};
						getYDDT(n->inputdiff,n->outputdiff,values);

						for (int k = 0; k < (1<<SIZE); k++)
						{
							n->possibleValues[k] = values[k];
							n->possibleValuesSum += values[k];
							for (int l = 0; l < values[k]; l++) n->possibleValuesVec.push_back(k);
						}
					}
					else if (it2 != savedValues.end())
					{
						// use the values from savedValues
						for (int k = 0; k < (1<<SIZE); k++)
						{
							// cout << k << " : " << savedValues[constraints[i][j]][k] << endl;
							n->possibleValues[k] = savedValues[constraints[i][j]][k];
							n->possibleValuesSum += savedValues[constraints[i][j]][k];
							for (int l = 0; l < savedValues[constraints[i][j]][k]; l++) n->possibleValuesVec.push_back(k);
						}
					}
					// update key information
					n->prev = getPrev(constraints[i][j], n, constraintsAddress,keys,localSavedValues,true);
					constraintsAddress[constraints[i][j]] = n;
				}
				else 
				{
					constraintsAddress[constraints[i][j]]->prev = getPrev(constraints[i][j], constraintsAddress[constraints[i][j]], constraintsAddress,keys, localSavedValues,true);
				}
				outputs.push_back(constraints[i][j]);
			}
		}
	}
	organizeKeys(keys,key_diff,keyDiff,keyAddress,KEYS0,KEYS1);
	getPossibleKeys(constraints, keyAddress, keyPossible);
	// add keys into all the nodes
	for (auto it = keyAddress.begin(); it != keyAddress.end(); it++){
		constraintsAddress[it->first]->keyIndex = it->second;
	}
}

bool evaluateNode(Node* n,const vector<vector<uint8_t>> &KEYS0, const vector<vector<uint8_t>> &KEYS1, int ind = -1, uint64_t val = 0)
{
	if (n->prev.size() == 0 && ind == -1)
	{
		int v = int(val % n->possibleValuesSum);
		int current = -1;
		while (v >= 0)
		{
			current++;
			v -= n->possibleValues[current];
		}
		n->value0 = current;
		n->value1 = n->value0 ^ n->outputdiff;
		return true;
	}
	else if (n->prev.size() > 0)
	{
		n->value0 = 0;
		n->value1 = 0;
		for (int i = 0; i < n->prev.size(); i++)
		{
			n->value0 ^= n->prev[i]->value0 ^ n->prev[i]->constant;
			n->value1 ^= n->prev[i]->value1 ^ n->prev[i]->constant;
			if (n->prev[i]->keyIndex != -1)
			{
				 n->value0 ^= KEYS0[n->prev[i]->keyIndex >> 4][n->prev[i]->keyIndex & 0xf];
				 n->value1 ^= KEYS1[n->prev[i]->keyIndex >> 4][n->prev[i]->keyIndex & 0xf];
			}
		}
		#if SIZE==4
		n->value0 = getSbox(n->value0);
		n->value1 = getSbox(n->value1);
		#elif SIZE==8
		n->value0 = getSbox8(n->value0);
		n->value1 = getSbox8(n->value1);
		#endif
	}
	else if ((n->prev.size() == 0) && (ind != -1)) // for theoretical computation
	{
		int current = -1;
		while (ind >= 0)
		{
			current++;
			ind -= n->possibleValues[current];
		}
		n->value0 = current;
		n->value1 = n->value0 ^ n->outputdiff;
	}
	return false;
}

void initializeKeys(vector<vector<uint8_t>> &KEYS0, vector<vector<uint8_t>> &KEYS1, vector<uint32_t> &keyDiff, map<uint32_t,vector<uint32_t>> keyPossible)
{
	// we only take key from keyPossible (if it exists) TBC
	for (int i = 0; i < KEYS0.size(); i++)
	{
		uint8_t tk[3] = {0}, tkp[3] = {0};
		for (int j = 0; j < TK_NUM; j++) {
			#if SIZE==4
			tk[j] = rand4(rand_generator);
			#elif SIZE==8
			tk[j] = rand8(rand_generator);
			#endif
			tkp[j] = tk[j] ^ ((keyDiff[i] >> ((2-j)*SIZE)) & ANDval);
		}
		auto it = keyPossible.find(i<<4);
		if (it != keyPossible.end()){
			uint8_t t = rand64(rand_generator) % keyPossible[i<<4].size();
			tk[TK_NUM-1] = keyPossible[i<<4][t];
			for (int j = 0; j < TK_NUM-1; j++) tk[TK_NUM-1] ^= tk[j];
			tkp[TK_NUM-1] = tk[TK_NUM-1] ^ ((keyDiff[i] >> ((3-TK_NUM)*SIZE)) & ANDval);
		}
		
		KEYS0[i][0] = 0;
		KEYS1[i][0] = 0;
		for (int j = 0; j < TK_NUM; j++) {
			KEYS0[i][0] ^= tk[j];
			KEYS1[i][0] ^= tkp[j];
		}

		for (int j = 1; j < KEYS0[i].size(); j++)
		{
			KEYS0[i][j] = 0;
			KEYS1[i][j] = 0;
			// update 
			for (int k = 1; k < TK_NUM; k++) {
				#if SIZE==4
				tk[k] = LFSR(tk[k],k+1);
				tkp[k] = LFSR(tkp[k],k+1);
				#elif SIZE==8
				tk[k] = LFSR8(tk[k],k+1);
				tkp[k] = LFSR8(tkp[k],k+1);
				#endif
			}
			for (int k = 0; k < TK_NUM; k++){
				KEYS0[i][j] ^= tk[k];
				KEYS1[i][j] ^= tkp[k];	
			}
		}
	}
}

bool runExperiment(uint64_t keys, uint64_t runs, vector<double> &counts_local, uint64_t &pks, uint64_t &tks, vector<vector<uint8_t>> KEYS0, vector<vector<uint8_t>> KEYS1, vector<uint32_t> keyDiff, map<uint32_t,Node*> constraintsAddress, map<uint32_t,uint64_t*> &localSavedValues, map<uint32_t,vector<uint32_t>> keyPossible, bool terminateIfPossible = false, bool randomKeys = false)
{
	uint64_t k = 0;
	tks = 0;
	pks = 0;
	vector<Node*> basicNodes; 
	for (auto it = constraintsAddress.begin(); it != constraintsAddress.end(); it++){
		if (it->second->prev.size() == 0) basicNodes.push_back(it->second);
	}
	uint32_t* random; random = new uint32_t [(1ULL<<(runs-1)) * basicNodes.size()];
	while (k < 128)
	{
		initializeKeys(KEYS0,KEYS1,keyDiff,keyPossible);
		uint64_t count = 0; 
		uint64_t next = 0;
		
		uint64_t random_size = (1ULL<<(runs-1)) * basicNodes.size();
		for (uint64_t i = 0; i < random_size/2; i++){
			random[i] = rand64(rand_generator);
		}

		for (uint64_t i = 0; i < (1ULL<<runs); i++)
		{
			if (next >= random_size/2 - basicNodes.size()){
				for (uint64_t j = 0; j < random_size/2; j++){
					random[j] = rand64(rand_generator);
				}
				next = 0;
			}
			
			if (randomKeys) initializeKeys(KEYS0,KEYS1,keyDiff,keyPossible);
			bool flag = true;
			for (auto it = constraintsAddress.begin(); it != constraintsAddress.end(); it++)
			{
				if (evaluateNode(it->second,KEYS0,KEYS1,-1,random[next])) next++;
				Node* n = it->second;
				if ((n->value0 ^ n->value1) != n->outputdiff) 
				{
					flag = false;
					break;
				}
			}
			if (!flag) continue;
			count++;
			// save only if the counts work
			for (auto it = constraintsAddress.begin(); it != constraintsAddress.end(); it++)
			{
				localSavedValues[it->first][it->second->value0]++;
				localSavedValues[it->first][it->second->value1]++;
			}
		}
		
		if (count > 0) 
		{
			if (terminateIfPossible) {
				delete[] random;
				return true;
			}
			counts_local.push_back(-log2((count+0.0)/(1<<runs)));
			k++;
			cout << -log2((count+0.0)/(1<<runs)) << endl;
		}
		tks++;
		if ((log2(tks) >= keys) && (k == 0)) break;
	}
	delete[] random;
	pks = counts_local.size();
	if (terminateIfPossible) return false;
	return true;
}


double getNumOutputs(vector<vector<uint32_t>> constraints)
{
	#if SIZE == 4
	int DDT4[16][16] = {0};
	computeDDT(DDT4);

	#elif SIZE == 8
	int DDT8[256][256] = {0};
	computeDDT8(DDT8);
	#endif
	vector<uint32_t> unique;
	double outputs = 0;
	for (int i = 0; i < constraints.size(); i++)
	{
		int maxN = findMaxN(constraints[i]);
		for (int j = 0; j < constraints[i].size(); j++)
		{
			if ((constraints[i][j] >> (4+2*SIZE)) == maxN) 
			{
				auto it = find(unique.begin(), unique.end(), constraints[i][j]);
				if (it == unique.end())
				{
					#if SIZE == 4
					outputs += -log2(DDT4[(constraints[i][j] >> SIZE) & ANDval][constraints[i][j] & ANDval]/16.0);
					#elif SIZE == 8
					outputs += -log2(DDT8[(constraints[i][j] >> SIZE) & ANDval][constraints[i][j] & ANDval]/256.0);
					#endif
					unique.push_back(constraints[i][j]);
				}
			}
		}
	}
	return outputs;
}

void removeSetUp(map<uint32_t,Node*> &constraintsAddress)
{
	for (auto it = constraintsAddress.begin(); it != constraintsAddress.end(); it++){
		delete it->second;
	}
}

void parallelWorker(vector<vector<uint32_t>> constraints, const uint8_t key_diff[4][4][4], map<uint32_t,uint64_t*> &enclosingSavedValues)
{

	map<uint32_t, Node*> constraintsAddress;
	map<uint32_t, uint32_t> keyAddress;
	vector<uint32_t> keyDiff; // TK0||TK1||TK2 (how to get this information?)
	vector<uint32_t> keys;
	vector<vector<uint8_t>> KEYS0;
	vector<vector<uint8_t>> KEYS1;
	bool final = true;
	map<uint32_t,uint64_t*> localSavedValues;
	map<uint32_t,vector<uint32_t>> keyPossible;
	setUp_mutex.lock();
	setUp(constraints,key_diff,constraintsAddress,keyAddress,keyDiff,keys,KEYS0,KEYS1,localSavedValues,keyPossible);
	double numOutputs = getNumOutputs(constraints);
	uint64_t numKeys = uint64_t(keys.size());
	vector<double> counts_local;
	uint64_t pks = 0; // possible keys
	uint64_t tks = 0; // total keys
	// should be 4
	cout << "computing with 128 valid keys or 2^" << max(uint64_t(8),numKeys*3) << " possible keys and 2^" << max(uint64_t(16),uint64_t(numOutputs+6)) << " random plaintexts" << endl;
	setUp_mutex.unlock();
	// cout << "computing with 2^" << max(uint64_t(8),numKeys*3) << " and 2^" << max(uint64_t(10),uint64_t(numOutputs*SIZE/2)) << endl;
 	runExperiment(max(uint64_t(8),numKeys*3),max(uint64_t(16),uint64_t(numOutputs+6)),counts_local,pks,tks,KEYS0,KEYS1,keyDiff,constraintsAddress,localSavedValues,keyPossible);	
	lock_guard<mutex> guard( term_mutex );
	counts.insert(counts.end(), counts_local.begin(), counts_local.end());
	
	for (auto it = keyPossible.begin(); it != keyPossible.end(); it++){
		uint64_t divisor = gcd(it->second.size(),1<<SIZE);
		pks *= it->second.size()/divisor;
		tks *= (1<<SIZE)/divisor;
	}
	possibleKeys += pks;
	totalKeys += tks;
	// update the localSavedValues into a enclosingSavedValues and eventually savedValues
	for (auto it = localSavedValues.begin(); it != localSavedValues.end(); it++)
	{
		auto it2 = enclosingSavedValues.find(it->first);	
		if (it2 == enclosingSavedValues.end()) enclosingSavedValues[it->first] = localSavedValues[it->first];
		else {
			for (int i = 0; i < (1<<SIZE); i++) enclosingSavedValues[it->first][i] += localSavedValues[it->first][i];
		}
	}
	removeSetUp(constraintsAddress);
}

double progressiveExperiment(vector<vector<uint32_t>> constraints)
{
	// we try some constraints first, if they are already impossible, don't bother trying further
	for (int i = 2; i < constraints.size(); i++)
	{
		break;
		map<uint32_t, Node*> constraintsAddress;
		map<uint32_t, uint32_t> keyAddress;
		vector<uint32_t> keyDiff;
		vector<uint32_t> keys;
		vector<vector<uint8_t>> KEYS0;
		vector<vector<uint8_t>> KEYS1;
		vector<vector<uint32_t>> tmpConstraints;
		map<uint32_t,vector<uint32_t>> keyPossible;
		for (int j = 0; j < i; j++) tmpConstraints.push_back(constraints[j]);
		map<uint32_t, uint64_t*> localSavedValues;
		setUp(tmpConstraints,key_diff,constraintsAddress,keyAddress,keyDiff,keys,KEYS0,KEYS1,localSavedValues,keyPossible);
		double numOutputs = getNumOutputs(tmpConstraints);
		uint64_t numKeys = KEYS0.size();
		uint64_t pks = 0, tks = 0;
		vector<double> counts_local;
		cout << "trying with just constraint 0 to constraint " << i-1 << endl;
		cout << "runnning with 2^" << numKeys*2 << " keys with 2^" << int(numOutputs+4) << " probability" << endl;
		if (!runExperiment(numKeys*2,uint64_t(numOutputs+4),counts_local,pks,tks,KEYS0,KEYS1,keyDiff,constraintsAddress,localSavedValues,keyPossible,true,true)) 
		{
			cout << "it is impossible with constraint 0 to constraint " << i - 1 << endl; 
			exit(0);
		}

	}
	counts.clear();
	possibleKeys = 0;
	totalKeys = 0;
	vector<thread> threads;
	int THREAD_NUMS = 16;
	map<uint32_t,uint64_t*> enclosingSavedValues;
	for (int i = 0; i < THREAD_NUMS; i++)
	{
		threads.push_back(
			thread( parallelWorker, constraints, key_diff, ref(enclosingSavedValues))
			);
	}
	for ( auto & th : threads ) 
	{
		th.join();
	}

	// putting the enclosingSavedValues into savedValues
	for (auto it = enclosingSavedValues.begin(); it != enclosingSavedValues.end(); it++)
	{
		auto it2 = savedValues.find(it->first);
		if (it2 == savedValues.end()) savedValues[it->first] = enclosingSavedValues[it->first];
		else {
			for (int i = 0; i < (1<<SIZE); i++) savedValues[it->first][i] *= enclosingSavedValues[it->first][i];
		}
	}
	if (possibleKeys == 0){
		cout << "Unable to find any keys that satisfy. Exiting..." << endl; exit(0);
	}
	return -log2((possibleKeys+0.0)/totalKeys);
}

void saveExp(string file)
{
	ofstream fout;
	fout.open(file);
	for (uint32_t i = 0; i < counts.size(); i++)
	{
		fout << counts[i]; fout << "\n";
	}
	fout.close();
}


double addSavedValues(vector<uint32_t> uniqueKeys, map<double,double> &distribution, uint8_t key_diff[4][4][4], vector<vector<uint32_t>> constraints){
	map<uint32_t,Node*> constraintsAddress;
	map<uint32_t,uint32_t> keyAddress;
	vector<uint32_t> keyDiff;
	vector<uint32_t> keys;
	vector<vector<uint8_t>> KEYS0,KEYS1;
	map<uint32_t,uint64_t*> localSavedValues;
	map<uint32_t,vector<uint32_t>> keyPossible;
	setUp(constraints,key_diff,constraintsAddress,keyAddress,keyDiff,keys,KEYS0,KEYS1,localSavedValues,keyPossible);
	for (int i = 0; i < KEYS0.size(); i++) {
		if ((KEYS0[i].size() > 1) && (TK_NUM > 1)) {
			cout << "KEYS0[i].size() > 1" << endl;
			cout << "fell into a case not covered by this algorithm" << endl;
		}
	}

	vector<uint32_t> inputs, others;
	for (auto it = constraintsAddress.begin(); it != constraintsAddress.end(); it++){
		if (it->second->prev.size() == 0) inputs.push_back(it->first);
		else others.push_back(it->first);
	}
	double runTime = 0;
	// for (int i = 0; i < constraints.size(); i++){
	// 	int maxN = findMaxN(constraints[i]);
	// 	for (int j = 0; j < constraints[i].size(); j++){
	// 		if ((constraints[i][j] >> (4+2*SIZE)) == maxN) others.push_back(constraints[i][j]);
	// 		else if ((constraints[i][j] & ANDval) == 0) others.push_back(constraints[i][j]);
	// 	}
	// }

	// for (int i = 0; i < constraints.size(); i++){
	// 	for (int j = 0; j < constraints[i].size(); j++){
	// 		auto it = find(others.begin(),others.end(),constraints[i][j]);
	// 		auto it2 = find(inputs.begin(),inputs.end(),constraints[i][j]);
	// 		if (it == others.end() && it2 == inputs.end())inputs.push_back(constraints[i][j]);
	// 	}
	// }
	int indices[constraintsAddress.size()] = {0};
	int indicesMax[constraintsAddress.size()] = {0};
	int i = 0;
	for (auto it = constraintsAddress.begin(); it != constraintsAddress.end(); it++){
		auto it2 = find(inputs.begin(),inputs.end(),it->first);
		if (it2 != inputs.end()) indicesMax[i] = it->second->possibleValuesSum;
		i++;
	}

	bool flag,breakFlag;
	int indexToMove;

	vector<uint64_t> counts;
	uint64_t total = 1;
	for (int i = 0; i < constraintsAddress.size(); i++){
		if (indicesMax[i] > 0) total *= indicesMax[i];
	}
	int totalKeys = 0;
	totalKeys += KEYS0.size();
	for (int i = 0; i < KEYS0.size(); i++){
		if (KEYS0[i].size() == 1) continue;
		else totalKeys += TK_NUM-1;
	}
	for (uint64_t k = 0; k < (1ULL << (totalKeys*SIZE)); k++){
		// set up the key
		uint64_t tmp_k = k;
		for (int i = 0; i < uniqueKeys.size(); i++){
			if (KEYS0[keyAddress[uniqueKeys[i]]>>4].size() == 1){
				uint8_t u = tmp_k & ANDval;
				tmp_k = tmp_k >> SIZE;
				KEYS0[keyAddress[uniqueKeys[i]]>>4][0] = u;
				KEYS1[keyAddress[uniqueKeys[i]]>>4][0] = KEYS0[keyAddress[uniqueKeys[i]]>>4][0];
				for (int j = 0; j < TK_NUM; j++) {
					KEYS1[keyAddress[uniqueKeys[i]]>>4][0] ^= ((keyDiff[keyAddress[uniqueKeys[i]]>>4] >> ((2-j)*SIZE)) & ANDval);
				}
			}
			else{
				uint8_t u0[TK_NUM], u1[TK_NUM];
				for (int l = 0; l < TK_NUM; l++){
					u0[l] = tmp_k & ANDval;
					u1[l] = u0[l] ^ ((keyDiff[keyAddress[uniqueKeys[i]]>>4] >> ((2-l)*SIZE)) & ANDval);
					tmp_k = tmp_k >> SIZE;
				}
				for (int j = 0; j < KEYS0[keyAddress[uniqueKeys[i]]>>4].size(); j++){
					KEYS0[keyAddress[uniqueKeys[i]]>>4][j] = 0;
					KEYS1[keyAddress[uniqueKeys[i]]>>4][j] = 0;
					for (int l = 0; l < TK_NUM; l++) {
						KEYS0[keyAddress[uniqueKeys[i]]>>4][j] ^= u0[l];
						KEYS1[keyAddress[uniqueKeys[i]]>>4][j] ^= u1[l];
					}
					for (int l = 1; l < TK_NUM; l++){
						#if SIZE==4
						u0[l] = LFSR(u0[l],l+1);
						u1[l] = LFSR(u1[l],l+1);
						#elif SIZE==8
						u0[l] = LFSR8(u0[l],l+1);
						u1[l] = LFSR8(u1[l],l+1);
						#endif
					}
				}
			}

		}
		uint64_t count = 0;
		while (true){
			flag = true;
			int i = 0;
			for (auto it = constraintsAddress.begin(); it != constraintsAddress.end(); it++){
				evaluateNode(it->second,KEYS0,KEYS1,indices[i]);
				i++;
				Node* n = it->second;
				if ((n->value0 ^ n->value1) != n->outputdiff) {
					flag = false;
					break;
				}
			}
			indexToMove = constraintsAddress.size() - 1;
			breakFlag = false;
			if (flag){
				for (auto it = constraintsAddress.begin(); it != constraintsAddress.end(); it++){
					localSavedValues[it->first][it->second->value0] += 1;
				}
				count++;
			}
			while (indices[indexToMove] >= indicesMax[indexToMove] - 1){
				indices[indexToMove] = 0;
				indexToMove--;
				if (indexToMove == -1){
					breakFlag = true;
					break;
				}
			}

			if (breakFlag) break;
			indices[indexToMove]++;
		}
		if (count > 0) counts.push_back(count);
	}
	// add them into savedValues
	for (auto it = localSavedValues.begin(); it != localSavedValues.end(); it++){
		auto it2 = savedValues.find(it->first);
		if (it2 == savedValues.end()){
			savedValues[it->first] = it->second;
		}
		else {
			for (int i = 0; i < (1<<SIZE); i++){
				savedValues[it->first][i] *= it->second[i];
			}
		}
		reduceKey(savedValues[it->first],1);
	}
	for (int i = 0; i < counts.size(); i++){
		double val = log2(counts[i]/(total+0.0));
		if (distribution.find(val) == distribution.end()){
			distribution[val] = 1;
		}
		else{
			distribution[val]++;
		}
	}
	for(auto it = distribution.begin(); it != distribution.end(); it++){
		it->second = it->second/(counts.size()+0.0) * 100;
	}
	if (counts.size() == 0) {cout << "This constraintGroup is impossible... Exiting..." << endl; exit(0);}
	if (counts.size() == 1ULL<<(totalKeys*SIZE)) return 0;
	else return log2(counts.size())-totalKeys*SIZE;
}

// int main()
// {
	
// 	// TODO
// 	// encode the higher order constraints!!
// 	// TK_NUM = 1;
// 	// get_4TK2_1_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // 6
// 	TK2_8_2020_1402(alpha,key_diff,nr,TK_NUM);
// 	fixedValuesToPath(alpha,key_diff,diffStates,fixedStates,nr);
// 	// uint8_t key_diff[4][4][4] = {0};

// 	// vector<uint32_t> constraint1 = {0x41b1,0x5518};
// 	// vector<uint32_t> constraint2 = {0x8284,0x9642};
// 	// vector<vector<uint32_t>> constraints = {constraint1,constraint2};

// 	// vector<uint32_t> constraint3 = {0x2652,0x14a5,0x1ba5,0x2900,0x3b21}; // nonlinear constraints
// 	// vector<vector<uint32_t>> constraints = {constraint3};
	
// 	// vector<uint32_t> constraint4 = {0x0d21, 0x0a21, 0x1418, 0x101b}; // higher order linear constraints
// 	// vector<vector<uint32_t>> constraints = {constraint4}; 

// 	vector<uint32_t> constraint1 = {0x210120, 0x352080};
// 	vector<uint32_t> constraint2 = {0x408003, 0x540320};
// 	vector<uint32_t> constraint3 = {0x408003, 0x4a8002, 0x5c0120};
// 	vector<uint32_t> constraint4 = {0x408003, 0x4a8002, 0x540320, 0x5c0120};
// 	vector<uint32_t> constraint5 = {0x520320, 0x662080};
// 	vector<uint32_t> constraint6 = {0x520320, 0x479002, 0x4a8002, 0x580000, 0x6e2093};
// 	vector<uint32_t> constraint7 = {0x520320, 0x479002, 0x4a8002, 0x580000, 0x5f0320, 0x620000, 0x68b080, 0x7e8003};
// 	vector<uint32_t> constraint8 = {0x520320, 0x550120, 0x580000, 0x6a2080, 0x6e2093};
// 	vector<uint32_t> constraint9 = {0x550120, 0x479002, 0x4a8002, 0x580000, 0x6a2080};
// 	vector<uint32_t> constraint10 = {0x662080, 0x540320, 0x5b0120, 0x690000, 0x7b8003};
// 	vector<vector<uint32_t>> constraints = {constraint1,constraint2,constraint3,constraint4,constraint5,
// 											constraint6,constraint7,constraint8,constraint9,constraint10,}; 
// 	// vector<vector<uint32_t>> constraints = {constraint4,constraint6,constraint7}; // impossible
// 	// vector<vector<uint32_t>> constraints = {constraint1,constraint2,constraint3,constraint4,constraint5,constraint6,constraint7,constraint8};
	
// 	progressiveExperiment(constraints);
// 	saveExp("tmp.txt");
// 	// cout << "Starting second experiment" << endl;
// 	// cout << "savedValues at this point:" << endl;
// 	// for (auto it = savedValues.begin(); it != savedValues.end(); it++)
// 	// {
// 	// 	cout << it->first << ": ";
// 	// 	for (int i = 0; i < (1<<SIZE); i++)
// 	// 	{
// 	// 		cout << it->second[i] << " ";
// 	// 	}
// 	// 	cout << endl;
// 	// }
// 	// vector<vector<uint32_t>> constraints2 = {constraint6};
// 	// progressiveExperiment(constraints2);
// 	return 0;
// }