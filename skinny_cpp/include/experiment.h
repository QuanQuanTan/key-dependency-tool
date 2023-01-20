#ifndef EXP_H
#define EXP_H

#include <vector>
#include <map>
#include <iostream>
using namespace std;

#if !defined SIZE
#error SIZE is not defined.
#endif

struct Node
{
	uint32_t index, pos;
	uint8_t outputdiff, inputdiff;
	uint8_t value0, value1;
	vector<Node*> prev;
	vector<Node*> next;
	int keyIndex = -1;
	uint8_t constant = 0;
	uint8_t possibleValues[1<<SIZE];
	uint64_t possibleValuesLength = 0;
	vector<uint32_t> possibleValuesVec;
};

vector<Node*> getPrev(uint32_t element, Node* elementAddress, map<uint32_t,Node*> constraintsAddress, vector<uint32_t> keys, map<uint32_t,uint64_t*> &localSavedValues, bool force = false);
vector<uint32_t> getKeys(vector<uint32_t> constraint);
uint32_t getKeyDiff(uint32_t element, const uint8_t key_diff[4][4][4]);
bool runExperiment(uint64_t keys, uint64_t runs, vector<double> counts_local, double pks, double totalKeys, vector<vector<uint8_t>> KEYS0, vector<vector<uint8_t>> KEYS1, vector<uint32_t> keyDiff, map<uint32_t,Node*> constraintsAddress ,bool terminateIfPossible = false, bool randomKeys = false);
int getNumOutputs(vector<vector<uint32_t>> constraints);
double progressiveExperiment(vector<vector<uint32_t>> constraints);

void printNode(Node* n);
void getPossibleKeys(vector<vector<uint32_t>> constraints, map<uint32_t,uint32_t> keyAddress, map<uint32_t,vector<uint32_t>> &keyPossible);
void organizeKeys(vector<uint32_t> keys, const uint8_t key_diff[4][4][4], vector<uint32_t> &keyDiff, map<uint32_t,uint32_t> &keyAddress, vector<vector<uint8_t>> &KEYS0,vector<vector<uint8_t>> &KEYS1);
void setUp(vector<vector<uint32_t>> constraints, const uint8_t key_diff[4][4][4], map<uint32_t,Node*> &constraintsAddress, map<uint32_t,uint32_t> &keyAddress, vector<uint32_t> &keyDiff, vector<uint32_t> &keys, vector<vector<uint8_t>>&KEYS0, vector<vector<uint8_t>>&KEYS1, map<uint32_t,uint64_t*> &localSavedValues, map<uint32_t,vector<uint32_t>> &keyPossible);
void evaluateNode(Node* n,vector<vector<uint8_t>> KEYS0, vector<vector<uint8_t>> KEYS1, int ind = -1);
void initializeKeys(vector<vector<uint8_t>> &KEYS0, vector<vector<uint8_t>> &KEYS1, vector<uint32_t> &keyDiff, map<uint32_t,vector<uint32_t>> keyPossible);
void removeSetUp(map<uint32_t,Node*> &constraintsAddress);
void parallelWorker(vector<vector<uint32_t>> constraints, const uint8_t key_diff[4][4][4], map<uint32_t,uint64_t*> &enclosingSavedValues);
void saveExp(string file);
double addSavedValues(vector<uint32_t> uniqueKeys, map<double,double> &distribution, uint8_t key_diff[4][4][4], vector<vector<uint32_t>> constraints);




#endif