#ifndef KEYDEP_H
#define KEYDEP_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <bits/stdc++.h>
#include <cryptominisat5/cryptominisat.h>

extern vector<vector<int>> keyConstraints;
extern int saved_prob;

#if !defined SIZE
#error SIZE is not defined.
#endif


void printState(uint8_t alpha[30][16], uint32_t nr);
void findEquations(uint8_t xEquations[16][16][32], uint8_t yEquations[16][16][32]);
void insertKeyConstraint(vector<int> vec, int parity);
void searchOne(vector<int> pos, vector<int> constraint, int parity, vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity);
void searchTwo(vector<int> pos, vector<int> constraint, int parity, vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity);
void searchThree(vector<int> pos, vector<int> constraint, int parity, vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity);
void searchFour(vector<int> pos, vector<int> constraint, int parity, vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity);
void findLinearConstraints(vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity, 
						   vector<vector<vector<int>>> yPositions, vector<vector<int>> yParity, uint32_t nr);
void AddConstants(vector<vector<vector<int>>> &yPositions, vector<vector<int>> &yParity, int nr);
void computeDimension(vector<vector<int>> masterKeyConstraints);
void computeProb(uint8_t alpha[30][16], uint32_t nr);

vector<vector<int>> findYConstraints(uint8_t alpha[30][16], int nr, vector<vector<vector<int>>> &positions);
vector<vector<int>> findXConstraints(uint8_t alpha[30][16], int nr, vector<vector<vector<int>>> &positions);
vector<int> checkXPos(vector<int> constraint);
vector<int> vectorXOR(vector<int> a, vector<int> b);
bool isZeroMatrix(uint8_t **M, int ROWSIZE, int COLSIZE);
vector<vector<int>> propagateLinear(int nr);
#endif