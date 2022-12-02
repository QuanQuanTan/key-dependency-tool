#ifndef SAT_H
#define SAT_H

#include <vector>
#include <stdint.h>
#include <cryptominisat5/cryptominisat.h>
#if !defined SIZE
#error SIZE is not defined.
#endif

using namespace std;
// void Tseitin(vector<CMSat::Lit> &clause, CMSat::SATSolver &solver, vector<uint32_t> cells, uint32_t values, int TK_NUM, uint32_t freeVar, uint64_t &noOfClauses);
void addSbox(uint8_t input, uint8_t output, CMSat::SATSolver &solver, int sbox_index, int n, int size);
void solve(uint8_t alpha[30][16], vector<vector<vector<int>>> xPositions, vector<vector<int>> xParity, vector<vector<vector<int>>> yPositions, vector<vector<int>> yParity, int nr, int size);

#endif