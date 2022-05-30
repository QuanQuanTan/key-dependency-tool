#include <cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <vector>
#include "SAT.h"


using namespace std;
using namespace CMSat;

#ifndef SIZE
#define SIZE 4
#endif

void Tseitin(vector<CMSat::Lit> &clause, CMSat::SATSolver &solver, vector<uint32_t> cells, uint32_t values, int TK_NUM, uint32_t freeVar, uint64_t &noOfClauses)
{
	// create a CNF formula using Tseitin transformation
	// except for the last step where we OR all the freeVar together
	for (uint32_t i = 0; i < cells.size(); i++)
	{
		int cellNum = cells[i];
		int cellValue = values >> ((SIZE*TK_NUM) * (cells.size()-i-1)); // note: no elimination of those in higher significant bits
		for (int tk = 0; tk < TK_NUM; tk++)
		{
			for (int s = 0; s < SIZE; s++)
			{
				noOfClauses++;
				clause.clear();
				clause.push_back(CMSat::Lit(cellNum*SIZE*TK_NUM+SIZE*tk+s,bool(((cellValue >> (SIZE*(TK_NUM-tk)-s-1)) + 1) % 2)));
				clause.push_back(CMSat::Lit(freeVar,true));
				solver.add_clause(clause);
			}
		}
	}
}


void SATsolver(vector<vector<int>> linearKeys,uint32_t **linearValues, uint32_t *linearValuesLength,
	vector<vector<int>> nonLinearKeys,uint32_t **nonLinearValues, uint32_t *nonLinearValuesLength, int TK_NUM)
{
	uint64_t noOfClauses = 0;

    CMSat::SATSolver solver;
    vector<CMSat::Lit> clause;

    //Let's use 4 threads
    solver.set_num_threads(8);

    // number of variables = 4 bits * 16 cells * TK2
    // here
    solver.new_vars(SIZE*16*TK_NUM); // 
    uint64_t freeVarPoint = SIZE*16*TK_NUM;
    for (uint32_t i = 0; i < linearKeys.size(); i++)
    {
    	vector<uint32_t> cells;
    	for (uint32_t j = 0; j < linearKeys[i].size(); j++) cells.push_back((linearKeys[i][j] >> (2*SIZE)) & 0xf);
    	uint64_t freeVarIncrement = freeVarPoint;
    	uint64_t incrementCounter = 0;
    	for (uint64_t j = 0; j < (1ULL<<(SIZE*linearValuesLength[i])); j++)
    	{
    		if (linearValues[i][j] > 0)
    		{
    			if (freeVarIncrement - freeVarPoint >= 4096*incrementCounter) 
    			{
    				incrementCounter++;
    				solver.new_vars(4096);
    			}
    			Tseitin(clause,solver,cells,j,TK_NUM,freeVarIncrement,noOfClauses);
    			freeVarIncrement++;
    		}
    	}
    	// combine to add the final clause in (final step of the Tsetin transformation)
    	clause.clear();
    	for (uint64_t j = freeVarPoint; j < freeVarIncrement; j++)
    	{
    		clause.push_back(CMSat::Lit(j,false));
    	}
    	noOfClauses++;
		solver.add_clause(clause);
		freeVarPoint = freeVarIncrement; // increase the counter by this amount
    }
	CMSat::lbool ret1 = solver.solve();
	cout << "Is it possible (linear)? : ";
    cout << ret1 << endl;
    if (ret1 == CMSat::l_False) return;


    // nonLinear case
    for (uint32_t i = 0; i < nonLinearKeys.size(); i++)
    {
    	cout << "Looking at nonLinearKey index " << i << "..." << endl;
    	vector<uint32_t> cells;
    	for (uint32_t j = 0; j < nonLinearKeys[i].size(); j++) 
    	{
    		if (nonLinearKeys[i][j] == -1) continue;
    		cells.push_back((nonLinearKeys[i][j] >> (2*SIZE)) & 0xf);
    	}
    	if (cells.size() >= 4)
    	{
    		cout << "nonLinearKey index " << i << " is too much for cryptominisat. Skipping..." << endl;
    		continue;	
    	} 
    	uint64_t freeVarIncrement = freeVarPoint;
    	uint64_t incrementCounter = 0;
    	for (uint64_t j = 0; j < (1ULL<<(SIZE*nonLinearValuesLength[i])); j++)
    	{
    		if (nonLinearValues[i][j] > 0)
    		{
    			if (freeVarIncrement - freeVarPoint >= 16384*incrementCounter) 
    			{
    				incrementCounter++;
    				solver.new_vars(16384);
    			}
    			Tseitin(clause,solver,cells,j,TK_NUM,freeVarIncrement,noOfClauses);
    			freeVarIncrement++;
    		}
    	}
    	// combine to add the final clause in (final step of the Tsetin transformation)
    	clause.clear();

    	for (uint64_t j = freeVarPoint; j < freeVarIncrement; j++)
    	{
    		clause.push_back(CMSat::Lit(j,false));
    	}
    	noOfClauses++;
		solver.add_clause(clause);
		freeVarPoint = freeVarIncrement; // increase the counter by this amount
    }

	cout << dec << "number of variables: " << solver.nVars() << endl;
	cout << dec << "number of clauses: " << noOfClauses << endl;

    cout << "Solving now..." << endl;
	CMSat::lbool ret = solver.solve();
    // assert(ret == CMSat::l_False);
    cout << "Is it possible (everything)? : ";
    cout << ret << endl;
}


void SATsolverBool(vector<vector<int>> linearKeys,uint32_t **linearValues, uint32_t *linearValuesLength,
	vector<vector<int>> nonLinearKeys,uint16_t **nonLinearValuesBool, uint32_t *nonLinearValuesLength, int TK_NUM)
{
	uint64_t noOfClauses = 0;

    CMSat::SATSolver solver;
    vector<CMSat::Lit> clause;

    //Let's use 4 threads
    solver.set_num_threads(8);

    // number of variables = 4 bits * 16 cells * TK2
    solver.new_vars(SIZE*16*TK_NUM); // 
    uint64_t freeVarPoint = SIZE*16*TK_NUM;
    for (uint32_t i = 0; i < linearKeys.size(); i++)
    {
    	vector<uint32_t> cells;
    	for (uint32_t j = 0; j < linearKeys[i].size(); j++) cells.push_back((linearKeys[i][j] >> (2*SIZE)) & 0xf);
    	uint64_t freeVarIncrement = freeVarPoint;
    	uint64_t incrementCounter = 0;
    	for (uint64_t j = 0; j < (1ULL<<(SIZE*linearValuesLength[i])); j++)
    	{
    		if (linearValues[i][j] > 0)
    		{
    			if (freeVarIncrement - freeVarPoint >= 4096*incrementCounter) 
    			{
    				incrementCounter++;
    				solver.new_vars(4096);
    			}
    			Tseitin(clause,solver,cells,j,TK_NUM,freeVarIncrement,noOfClauses);
    			freeVarIncrement++;
    		}
    	}
    	// combine to add the final clause in (final step of the Tsetin transformation)
    	clause.clear();
    	for (uint64_t j = freeVarPoint; j < freeVarIncrement; j++)
    	{
    		clause.push_back(CMSat::Lit(j,false));
    	}
    	noOfClauses++;
		solver.add_clause(clause);
		freeVarPoint = freeVarIncrement; // increase the counter by this amount
    }
	CMSat::lbool ret1 = solver.solve();
	cout << "Is it possible (linear)? : ";
    cout << ret1 << endl;
    if (ret1 == CMSat::l_False) return;

    // nonLinear case
    for (uint32_t i = 0; i < nonLinearKeys.size(); i++)
    {
    	cout << "Looking at nonLinearKey index " << i << "..." << endl;
    	vector<uint32_t> cells;
    	for (uint32_t j = 0; j < nonLinearKeys[i].size(); j++) 
    	{
    		if (nonLinearKeys[i][j] == -1) continue;
    		cells.push_back((nonLinearKeys[i][j] >> (2*SIZE)) & 0xf);
    	}
    	if (cells.size() >= 4)
    	{
    		cout << "nonLinearKey index " << i << " is too much for cryptominisat. Skipping..." << endl;
    		continue;
    	} 
    	uint64_t freeVarIncrement = freeVarPoint;
    	uint64_t incrementCounter = 0;
    	for (uint64_t j = 0; j < (1ULL<<(SIZE*nonLinearValuesLength[i])); j++)
    	{
    		uint64_t index = j / 16;
    		uint16_t remainder = j % 16;
    		if (((nonLinearValuesBool[i][index] >> (15-remainder)) & 1) > 0)
    		{
    			if (freeVarIncrement - freeVarPoint >= 16384*incrementCounter) 
    			{
    				incrementCounter++;
    				solver.new_vars(16384);
    			}
    			Tseitin(clause,solver,cells,j,TK_NUM,freeVarIncrement,noOfClauses);
    			freeVarIncrement++;
    		}
    	}
    	// combine to add the final clause in (final step of the Tsetin transformation)
    	clause.clear();

    	for (uint64_t j = freeVarPoint; j < freeVarIncrement; j++)
    	{
    		clause.push_back(CMSat::Lit(j,false));
    	}
    	noOfClauses++;
		solver.add_clause(clause);
		freeVarPoint = freeVarIncrement; // increase the counter by this amount
    }

	cout << dec << "number of variables: " << solver.nVars() << endl;
	cout << dec << "number of clauses: " << noOfClauses << endl;

    cout << "Solving now..." << endl;
	CMSat::lbool ret = solver.solve();
    // assert(ret == CMSat::l_False);
    cout << "Is it possible (everything)? : ";
    cout << ret << endl;
}