#ifndef PRINTDET_H
#define PRINTDET_H
#include <vector>
#include <iostream>

using namespace std;

void printState(uint8_t diffStates[20][5][16],uint8_t fixedStates[20][5][16], int nr);
void printConstraint(vector<uint32_t> constraint);
void printAllConstraints(vector<vector<uint32_t>> constraints);
void printConstraints(vector<vector<uint32_t>> constraints);
#endif