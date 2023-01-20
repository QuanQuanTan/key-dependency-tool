#include "skinny.h"
#include "trails.h"
#include "keyDependent.h"
#include "printDetails.h"
#include "experiment.h"
#include <sstream>
#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl

int TRAIL_INDEX;
vector<vector<double>> probabilities;
vector<vector<double>> prob_distribution;
map<int,double> finalDistribution;
double keysReduction = 0;

void testTrail()
{

	double prob = getOriginalProbability();
	cout << "original prob: " << prob << endl;
	prob -= keyProb; // independent of the ones involved in key constraints (nonlinear constraints only)
	cout << dec << "probability if we exclude those that are not involved in any constraints: " << prob << endl;


	fixedValuesToPath(alpha,key_diff,diffStates,fixedStates,nr);
	printState(diffStates,fixedStates,nr);
	// this finds linear/nonlinear constraints
	vector<vector<uint32_t>> constraints = findConstraints(diffStates,fixedStates,nr);
	// this is for higher order constraints
	vector<vector<uint32_t>> equations = getEquations();
	getHOLinearConstraints(equations,diffStates,fixedStates,nr);
	printAllConstraints(constraints);
	// check for any obvious incompatibility by those constraints without key
	checkHOLinearCompatbility();
	
	cout << dec << "Original prob: " << getOriginalProbability() << endl;
	cout << dec << "Independent prob (prob ignoring the ones depend on keys): " << getOriginalProbability()-keyProb << endl;
	if (keyProb == 0) return;

	vector<vector<vector<uint32_t>>> constraintsGroup = splitConstraints(constraints,HOLinearConstraints);
	cout << "These are the splits:" << endl;
	for (int i = 0; i < constraintsGroup.size(); i++)
	{
		for (int j = 0; j < constraintsGroup[i].size(); j++)
		{
			printConstraint(constraintsGroup[i][j]);
			cout << endl;
		}
		cout << "-----------------------------------------------" << endl;
	}

		cout << "===============================================" << endl;
	// we try to solve those without any higher order linear constraints

	vector<vector<vector<uint32_t>>> experimentials;
	for (int l = 0; l < constraintsGroup.size(); l++){
		bool computed = computable(constraintsGroup[l]);
		if (computed){
			cout << "Computing the following constraints theoretically: " << endl;
			printConstraints(constraintsGroup[l]); cout << endl;
			uint32_t* goodKeys;
			vector<uint32_t> uniqueKeys;
			findUniqueKeys(uniqueKeys,constraintsGroup[l]);
			map<double,double> distribution;
			double keys = addSavedValues(uniqueKeys,distribution,key_diff,constraintsGroup[l]);
			cout << "key space reduction: " << keys << endl;
			keysReduction -= keys;

			vector<double> tmp, tmp2;
			probabilities.push_back(tmp);
			prob_distribution.push_back(tmp2);
			for (auto it = distribution.begin(); it != distribution.end(); it++)
			{
				probabilities[probabilities.size()-1].push_back(-it->first);
				prob_distribution[prob_distribution.size()-1].push_back(it->second);
				cout << it->first << ": " << it->second << endl;
			}
			cout << "-------------------------------------" << endl;
		}
		if (!computed) experimentials.push_back(constraintsGroup[l]); // to be processed later
	}

	// conclude the theoretical distributions
	cout << "==============================================" << endl;
	cout << "For the theoretical side: " << endl;
	cout << "The following are the distribution: " << endl;
	int index[probabilities.size()] = {0};
	int indexToMove = probabilities.size()-1;
	while (true)
	{
		double prob = 0;
		double prob_dist = 1;
		for (int i = 0; i < probabilities.size(); i++)
		{
			prob += probabilities[i][index[i]];
			prob_dist *= prob_distribution[i][index[i]]/100.0;
		}
		auto it = finalDistribution.find(int(prob*1000));
		if (it == finalDistribution.end()) finalDistribution[int(prob*1000)] = prob_dist * 100;
		else finalDistribution[int(prob*1000)] += prob_dist * 100;
		index[indexToMove]++;
		while (index[indexToMove] >= probabilities[indexToMove].size())
		{
			index[indexToMove] = 0;
			indexToMove--;
			if (indexToMove == -1) break;
			index[indexToMove]++;
		}
		if (indexToMove == -1) break;
		indexToMove = probabilities.size()-1;
	}
	cout << "Key reduction: " << keysReduction << endl;
	for (auto it = finalDistribution.begin(); it != finalDistribution.end(); it++)
	{
		cout << (it->first/1000.0) << ": " << it->second << endl;
	}
	cout << "==============================================" << endl;
	savedValues.clear();
	for (int l = 0; l < experimentials.size(); l++)
	{
		cout << "----------------------------------------------" << endl;
		cout << "Computing the following experimentally:" << endl;
		printConstraints(experimentials[l]);
		double prob = progressiveExperiment(experimentials[l]);
		cout << "key space reduction: " << prob << endl;
		keysReduction += prob;
		string filename = "exp_" + to_string(SIZE) + "_" + to_string(TRAIL_INDEX) + "_" + to_string(l);
		saveExp(filename + ".txt");;
		string cmd = "python3 parseExperiment.py " + filename + ".txt";
		system(cmd.c_str());

		ifstream fin; fin.open(filename + "_dist.txt");
		string line;
		double proba, dist;
		vector<double> prob_tmp, dist_tmp;
		cout << "prob, dist: " << endl;
		while (std::getline(fin, line)) {
			istringstream ss(line);
			ss >> proba >> dist;
			cout << proba << " " << dist << endl;
			prob_tmp.push_back(proba);
			dist_tmp.push_back(dist);
		}
		fin.close();
		probabilities.push_back(prob_tmp);
		prob_distribution.push_back(dist_tmp);
	}

	// for (int i = 0; i < probabilities.size(); i++){
	// 	for (int j = 0; j < probabilities[i].size(); j++){
	// 		cout << probabilities[i][j] << " ";
	// 	}
	// 	cout << endl;
	// 	for (int j = 0; j < prob_distribution[i].size(); j++){
	// 		cout << prob_distribution[i][j] << " ";
	// 	}
	// 	cout << endl;
	// 	cout << "----------------------" << endl;
	// }

	cout << "==============================================" << endl;
	cout << "Overall: " << endl;
	cout << "The following is the probability distribution: " << endl;
	int index2[probabilities.size()] = {0};
	indexToMove = probabilities.size()-1;
	finalDistribution.clear();
	while (true)
	{
		double prob = getOriginalProbability() - keyProb;
		double prob_dist = 1;
		for (int i = 0; i < probabilities.size(); i++)
		{
			prob += probabilities[i][index2[i]];
			prob_dist *= prob_distribution[i][index2[i]]/100.0;
		}
		auto it = finalDistribution.find(int(prob*1000));
		if (it == finalDistribution.end()) finalDistribution[int(prob*1000)] = prob_dist * 100;
		else finalDistribution[int(prob*1000)] += prob_dist * 100;
		index2[indexToMove]++;
		while (index2[indexToMove] >= probabilities[indexToMove].size())
		{
			index2[indexToMove] = 0;
			indexToMove--;
			if (indexToMove == -1) break;
			index2[indexToMove]++;
		}
		if (indexToMove == -1) break;
		indexToMove = probabilities.size()-1;
	}
	cout << "Key reduction: " << keysReduction << endl;
	for (auto it = finalDistribution.begin(); it != finalDistribution.end(); it++)
	{
		cout << (it->first/1000.0) << ": " << it->second << endl;
	}
	cout << "==============================================" << endl;

}


int main(int argc, char **argv)
{
	if (argc == 1)
	{
		cout << "Please insert the trail number required" << endl;
		exit(0);
	}
	TRAIL_INDEX = atoi(argv[1]); // I am not going to check it
	// SKINNY64
	#if SIZE == 4
	if (TRAIL_INDEX==0) Corrected_TK1_2020_1402(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==1) SK_2020_1402(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==2) TK1_2020_1402(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (TRAIL_INDEX==3) TK2_2020_1402(alpha,key_diff,nr,TK_NUM); 
	if (TRAIL_INDEX==4) TK3_2020_1402(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==5) get_4TK2_1_17_L_2020_1317(alpha,key_diff,nr,TK_NUM); // this is used as a first test
	if (TRAIL_INDEX==6) get_4TK2_1_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // without nonLinearConstraint[0], it's possible
	if (TRAIL_INDEX==7) get_4TK2_1_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	if (TRAIL_INDEX==8) get_4TK2_1_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // without nonLinearConstraint[0], it's possible
	if (TRAIL_INDEX==9) get_4TK2_2_17_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (TRAIL_INDEX==10) get_4TK2_2_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	if (TRAIL_INDEX==11) get_4TK2_2_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (TRAIL_INDEX==12) get_4TK2_2_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (TRAIL_INDEX==13) get_4TK3_1_22_L_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	if (TRAIL_INDEX==14) get_4TK3_1_22_U_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (TRAIL_INDEX==15) get_4TK3_1_22_L_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	if (TRAIL_INDEX==16) get_4TK3_1_22_U_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (TRAIL_INDEX==17) get_4TK3_1_22_U_BMD3_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (TRAIL_INDEX==18) get_4TK3_1_23_U_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (TRAIL_INDEX==19) get_4TK3_1_23_U_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible linear
	// get_4TK2_18_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // u = 6, v = 2,10,12,13 // all ok
	// get_4TK3_22_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // u = 8,9, v = 5 // all ok


	// 8 bit trails
	// SKINNY128
	#elif SIZE == 8
	if (TRAIL_INDEX==1) TK1_8_2020_1402(alpha,key_diff,nr,TK_NUM); 
	if (TRAIL_INDEX==2) SK_tosc_2017_i4_99_129(alpha,key_diff,nr,TK_NUM);
	

	if (TRAIL_INDEX==3) TK2_8_2020_1402(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==4) get_8TK2_1_18_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==5) get_8TK2_1_18_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==6) get_8TK2_1_19_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==7) get_8TK2_1_20_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==8) get_8TK2_1_20_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==9) get_8TK2_1_21_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	
	if (TRAIL_INDEX==10) get_8TK2_2_18_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==11) get_8TK2_2_18_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==12) get_8TK2_2_19_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==13) get_8TK2_2_19_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==14) get_8TK2_2_20_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK2_2_21_U_2020_1317(alpha,key_diff,nr,TK_NUM); // repeated trail
	// uint8_t u1[4]= {0x40,0x50,0xc0,0xd0};  uint8_t v1[11] =  {0x2a,0x2e,0x2f,0x6a,0xba,0xbe,0xbf,0xea,0xef,0xfa,0xfe};
	// uint8_t w1[18] = {0x2a,0x2d,0x2e,0x2f,0x6b,0x7c,0xb8,0xb9,0xba,0xbd,0xbe,0xbf,0xeb,0xec,0xef,0xfb,0xfc,0xfe};
	// u = 0x40,0x50,0xc0,0xd0  v =  0x2a,0x2e,0x2f,0x6a,0xba,0xbe,0xbf,0xea,0xef,0xfa,0xfe
	// w = 0x2a,0x2d,0x2e,0x2f,0x6b,0x7c,0xb8,0xb9,0xba,0xbd,0xbe,0xbf,0xeb,0xec,0xef,0xfb,0xfc,0xfe
	// uint8_t u = u1[q % 4];
	// uint8_t v = v1[(q / 4) % 11];
	// uint8_t w = w1[(q / 44)];
	// int u = 0x40, v = 0x2a, w = 0x2a;
 	// get_8TK2_19_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v,w); // all ok


	if (TRAIL_INDEX==15) get_8TK3_1_22_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==16) get_8TK3_1_22_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==17) get_8TK3_1_23_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==18) get_8TK3_1_23_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_24_L_2020_1317(alpha,key_diff,nr,TK_NUM); // repeated
	if (TRAIL_INDEX==19) get_8TK3_1_24_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (TRAIL_INDEX==20) get_8TK3_1_25_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	#endif
	// get_8TK3_1_25_U_2020_1317(alpha,key_diff,nr,TK_NUM);	// repeated trails

	// uint8_t u1[5] = {0x08,0x0c,0x0d,0x0e,0x0f};
	// uint8_t v1[7] = {0x28, 0x29, 0x2c, 0x2d, 0x2e, 0x2f, 0xb9};
	// uint8_t u = u1[q % 5];
	// uint8_t v = v1[(q/5)];
	
	// u = 0x08,0x0c,0x0d,0x0e,0x0f
	// v = 0x28, 0x29, 0x2c, 0x2d, 0x2e, 0x2f, 0xb9	
	// int u = 0x08, v = 0x28;
	// get_8TK3_22_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // all ok
	
	
	testTrail();
	return 0;
}

