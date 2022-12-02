#ifndef KEYDEP_H
#define KEYDEP_H

#include <iostream>
#include <vector>
#include <numeric>
#include <bits/stdc++.h>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <unordered_set>

using namespace std;

#if !defined SIZE
#error SIZE is not defined.
#endif

// general requirements
extern int TK_NUM;
extern uint32_t nr;
extern uint8_t alpha[20][4][4]; // differential characteristic
extern uint8_t key_diff[4][4][4]; // related key 
extern uint8_t diffStates[20][5][16];
extern uint8_t fixedStates[20][5][16];
extern uint32_t ANDval;
extern uint64_t ROUNDval;
extern int inverseMC[16][3];
extern double keyProb;


// linear constraints requirements
extern vector<vector<uint32_t>> linearConstraints;
extern vector<uint32_t> linearKeys;
extern vector<uint32_t> propagatedLinearKeys;
extern vector<int> linearProbReduce; // contain the probability (in log2) for a linear constraint
extern vector<int> linearProb; // contain the probability (in log2) for a linear constraint
extern uint32_t ** linearLast;
extern uint32_t **linearDistribution;
extern uint32_t *linearDistributionLength;
extern uint32_t *linearDistributionBase;
extern uint64_t *linearPossibleKeys;

// non-linear constraints requirements
extern vector<vector<uint32_t>> nonLinearConstraints;
extern vector<vector<uint32_t>> nonLinearKeys;
extern vector<uint32_t> propagatedNonLinearKeys;
extern uint32_t **nonLinearValues;
extern uint32_t *nonLinearValuesLength;
extern uint32_t **nonLinearLast; // contain the count for the last Sbox transition (getXDDT)
extern uint64_t *nonLinearPossibleKeys; // keeps track if a key is possible due to the constraint. (For counting possible keys)
extern uint16_t *combinedTimes; // number of "values" at the end of each tuple
extern uint32_t **nonLinearDistribution; // 2*jth entry keep track of the count, 2*j+1 keep track of the total. Divide to get prob
extern uint32_t *nonLinearDistributionLength; 
extern uint32_t *nonLinearDistributionBase;


// same round constraints
extern vector<vector<uint32_t>> sameRoundConstraints;

// special constraints
extern vector<vector<uint32_t>> dualKeys;
extern vector<uint32_t> singleKeys;
extern uint32_t** DValues;
extern uint32_t** SValues;
extern uint32_t DValuesLength;
extern uint32_t SValuesLength;


// experimental constraints
extern vector<vector<uint32_t>> expNonLinearKeys;
extern vector<uint32_t> expLinearKeys;

bool isLinear(vector<uint32_t> constraint);
void fixedValuesToPath(uint8_t characteristic[20][4][4], uint8_t key_diff[4][4][4], uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr);
vector<vector<uint32_t>> findConstraints(uint8_t diffStates[20][5][16],uint8_t fixedStates[20][5][16],int nr);
void reduceKey(uint32_t* &val, uint32_t valLength) ;
void getXDDT(uint32_t input_diff, uint32_t output_diff, uint32_t output[1<<SIZE]);
void getYDDT(uint32_t input_diff, uint32_t output_diff, uint32_t output[1<<SIZE]);
void XOR(uint32_t val1[(1<<SIZE)],uint32_t val2[(1<<SIZE)]);
void addRoundConstant(uint32_t tup, uint32_t val[(1<<SIZE)]);
uint32_t resolveLinear(vector<uint32_t> constraint, uint32_t* distribution, uint32_t &distributionLength, uint32_t &distributionBase);
vector<uint32_t> resolveLeft(uint32_t zeroPoint, vector<uint32_t> constraint, uint32_t values[(1<<SIZE)], vector<uint32_t> &zeroDependencies, bool last);
bool AddKeyNonLinear(uint32_t* &value,uint32_t &valueLength);
void SubstNonLinear(uint32_t* &value,uint32_t &valueLength);
void resolveZeroDependencies(uint32_t* valuesToAdd, uint32_t valuesToAddLength, uint32_t* &destination, uint32_t &destinationLength);
void AddFinalgetYDDT(vector<uint32_t> constraint, uint32_t* &val, uint32_t valLength);
bool resolveNonLinear(vector<uint32_t> constraint, uint32_t* &values, uint32_t &valuesLength, uint32_t lastValues[1<<SIZE], vector<uint32_t> &output);
bool combineNonLinearConstraints(uint32_t index_i, uint32_t index_j);
void combineConstraints();
void computeNonLinearDistribution();
void computeNonLinearPossibleKeys();
void computeLinearPossibleKeys();
double getOriginalProbability();
void remove(vector<int> &v) ;
vector<vector<int>> getSubgraph(vector<vector<int>> G, bool keysDetermined[16], vector<int> &keyPositions);
double subGraphKey(vector<int> keyPositions, vector<vector<int>> constraintsIndex);
double computeKeys();
void computeProbability();
bool sizeCheck(vector<uint32_t> constraint);
vector<vector<uint32_t>> getEquations();
void getSameRoundConstraints(vector<vector<uint32_t>> equations, uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr);
void addLinearConstraint(uint32_t key, uint32_t* val);
void addNonLinearConstraint(vector<uint32_t> key, uint32_t* val);
void propagateKeys();
void printSameRoundConstraints();
#endif