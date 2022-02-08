#ifndef SAT_H
#define SAT_H

#include <vector>
#include <cryptominisat5/cryptominisat.h>
using namespace std;

void Tseitin(vector<CMSat::Lit> &clause, CMSat::SATSolver &solver, vector<uint32_t> cells, uint32_t values, int TK_NUM, uint32_t freeVar, uint64_t &noOfClauses);
void SATsolver(vector<vector<int>> linearKeys,uint32_t **linearValues, uint32_t *linearValuesLength, 
	vector<vector<int>> nonLinearKeys, uint32_t **nonLinearValues, uint32_t *nonLinearValuesLength, int TK_NUM);
void SATsolverBool(vector<vector<int>> linearKeys,uint32_t **linearValues, uint32_t *linearValuesLength,
	vector<vector<int>> nonLinearKeys,uint16_t **nonLinearValuesBool, uint32_t *nonLinearValuesLength, int TK_NUM);
void SATsolverLinear(vector<vector<int>> linearKeys, uint32_t **linearValues, uint32_t *linearValuesLength, int TK_NUM);
#endif