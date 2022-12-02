#include <iostream>
#include "printDetails.h"
#include <vector>
#include <iomanip>
#include "skinny.h"
#include "keyDependent.h"
using namespace std;

void printState(uint8_t diffStates[20][5][16],uint8_t fixedStates[20][5][16], int nr){
	// for visual check
	for (int n = 0; n < nr; n++)
	{
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][0][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][1][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][2][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][3][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(diffStates[n][4][4*r+c]) << ",";
			cout << endl;
		}
		cout << endl;
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][0][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][1][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][2][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][3][4*r+c]) << ",";
			cout << " ";
			for (int c = 0; c < 4; c++) cout << setfill('0') << setw(SIZE/4) << hex << int(fixedStates[n][4][4*r+c]) << ",";
			cout << endl;
		}
		cout << endl;
		cout << "-------------------------------------------------------" << endl;
	}
}

void printConstraint(vector<uint32_t> constraint){
	// for printing purposes only
	for (uint32_t j = 0; j < constraint.size(); j++){
		cout << hex << "(" << (constraint[j] >> (4+2*SIZE)) << "," << ((constraint[j] >> (2*SIZE)) & 0xf) << ","
		<< ((constraint[j] >> SIZE) & ANDval) << "," << (constraint[j] & ANDval) << ") ";
	}
}

void printConstraints(vector<vector<uint32_t>> constraints){
	// for printing purposes only
	int max_n = 0; // this computes the maximum round of all constraints
	for (int i = 0; i < constraints.size(); i++)
	{
		for (int j = 0; j < constraints[i].size(); j++)
		{
			if ((constraints[i][j] >> (4+2*SIZE)) > max_n)
			{
				max_n = (constraints[i][j] >> (4+2*SIZE));
			}
		}
	}
	cout << "The constraints are (in the form of (round, position, input_diff, output_diff)" << endl;
	cout << "The number infront of each constraint corresponds to the grouping that each constraint falls under." << endl;
	cout << "Each group is treated as independent in this computation" << endl;
	cout << "Some constraints may have more than one group identity. These groups have to be considered together" << endl;
	for (uint32_t i = 0; i < constraints.size(); i++){
		int last_n = constraints[i][constraints[i].size()-1] >> (4 + 2*SIZE); // last_n is the max round of this particular constraint
		vector<uint32_t> key; 
		for (int k = 0; k < constraints[i].size(); k++){
			if (((constraints[i][k] >> (2*SIZE)) & 0xf) < 8) key.push_back(constraints[i][k]);
		}
		vector<int> positions;
		for (int k = 0; k < key.size(); k++){
			int n = key[k] >> (4+2*SIZE);
			if (n == last_n) continue; // key is not involved when it is on the last round
			int pos = (key[k] >> (2*SIZE)) & 0xf;
			while (n < max_n)
			{
				pos = getInvPermSchedule(pos);
				n++;
			}
			positions.push_back(pos);
		}
		for (int k = 0; k < positions.size(); k++) cout << hex << positions[k] << " ";
		printConstraint(constraints[i]);
		if (isLinear(constraints[i])) cout << "linear" << endl;
		else cout << "nonlinear" << endl;
	}
	for (uint32_t i = 0; i < sameRoundConstraints.size(); i++){
		vector<uint32_t> key; 
		for (int k = 0; k < sameRoundConstraints[i].size(); k++){
			if (((sameRoundConstraints[i][k] >> (2*SIZE)) & 0xf) < 8) key.push_back(sameRoundConstraints[i][k]);
		}
		vector<int> positions;
		for (int k = 0; k < key.size(); k++){
			int n = key[k] >> (4+2*SIZE);
			int pos = (key[k] >> (2*SIZE)) & 0xf;
			while (n < max_n)
			{
				pos = getInvPermSchedule(pos);
				n++;
			}
			positions.push_back(pos);
		}
		for (int k = 0; k < positions.size(); k++) cout << positions[k] << " ";
		printConstraint(sameRoundConstraints[i]);
		cout << " higher-order linear" << endl;
	}
}
