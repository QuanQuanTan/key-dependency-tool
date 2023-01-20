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

// same round constraints
extern vector<vector<uint32_t>> HOLinearConstraints;
extern vector<vector<uint32_t>> HOLinearKeys;

extern vector<vector<uint32_t>> constraintsGroupIndex;

extern map<uint32_t,uint64_t*> savedValues; 

bool isLinear(vector<uint32_t> constraint);
void fixedValuesToPath(uint8_t characteristic[20][4][4], uint8_t key_diff[4][4][4], uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr);
vector<vector<uint32_t>> findConstraints(uint8_t diffStates[20][5][16],uint8_t fixedStates[20][5][16],int nr);
void reduceKey(uint32_t* &val, uint32_t valLength);
void reduceKey(uint64_t* &val, uint32_t valLength);
void getXDDT(uint32_t input_diff, uint32_t output_diff, uint32_t output[1<<SIZE]);
void getYDDT(uint32_t input_diff, uint32_t output_diff, uint32_t output[1<<SIZE]);
void XOR(uint32_t val1[(1<<SIZE)],uint32_t val2[(1<<SIZE)]);
void addRoundConstant(uint32_t tup, uint32_t val[(1<<SIZE)]);
double getOriginalProbability();
void findUniqueKeys(vector<uint32_t> &uniqueKeys, vector<vector<uint32_t>> constraints);
vector<vector<uint32_t>> getEquations();
void getHOLinearConstraints(vector<vector<uint32_t>> equations, uint8_t diffStates[20][5][16], uint8_t fixedStates[20][5][16], int nr);
void checkHOLinearCompatbility();
uint32_t findMaxN(vector<uint32_t> constraint);
vector<vector<vector<uint32_t>>> splitConstraints(vector<vector<uint32_t>> constraints, vector<vector<uint32_t>> HOLinearConstraints);
bool computable(vector<vector<uint32_t>> constraints);

#endif