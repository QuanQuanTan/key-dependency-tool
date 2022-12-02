#include "skinny.h"
#include "trails.h"
#include "keyDependent.h"
#include "printDetails.h"

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl
void testTrail()
{
	cout << "Note that for this version of the algorithm, we are able to detect higher order linear constraints, but unable to process it. Hence, these higher-order constraints will be avoided. In the meantime, please perform an experiment." << endl;
	DValues = new uint32_t* [4*nr];
	for (int i = 0; i < 4*nr; i++) DValues[i] = new uint32_t [1<<SIZE]();
	DValuesLength = 0;

	SValues = new uint32_t* [4*nr];
	for (int i = 0; i < 4*nr; i++) SValues[i] = new uint32_t [1<<SIZE]();
	SValuesLength = 0;

	fixedValuesToPath(alpha,key_diff,diffStates,fixedStates,nr);
	vector<vector<uint32_t>> constraints = findConstraints(diffStates,fixedStates,nr);
	vector<vector<uint32_t>> equations = getEquations();
	
	// printState(diffStates,fixedStates,nr);
	getSameRoundConstraints(equations,diffStates,fixedStates,nr);
	printConstraints(constraints);
	cout << dec << "keyProb:" << keyProb << endl;

	cout << "Original prob: " << getOriginalProbability() << endl;
	cout << "Independent prob: " << getOriginalProbability()-keyProb << endl;
	// cout << "We have " << DValuesLength << " same round constraints in the form \"k0 + k4\" = {...}" << endl;
	// cout << "We have " << SValuesLength << " same round constraints in the form \"k0\" = {...} (Equivalent to linear constraints)" << endl;
	for (uint32_t i = 0; i < constraints.size(); i++) // split linear and nonlinear constraints
	{
		if (isLinear(constraints[i])) linearConstraints.push_back(constraints[i]); 
		else if (!isLinear(constraints[i]))
		{
			if (sizeCheck(constraints[i])) nonLinearConstraints.push_back(constraints[i]);
		}
		else {cout << "error!" << endl;}
	}
	// work on the linear constraints

	linearDistribution = new uint32_t*[linearConstraints.size()+SValuesLength]();
	linearDistributionLength = new uint32_t[linearConstraints.size()+SValuesLength]();
	linearDistributionBase = new uint32_t[linearConstraints.size()+SValuesLength]();
	linearPossibleKeys = new uint64_t[linearConstraints.size()+SValuesLength]();
	cout << "Allocated " << linearConstraints.size() << " units for linear constraints" << endl;

	// work on the nonlinear constraints
	nonLinearValues = new uint32_t*[nonLinearConstraints.size()+DValuesLength]();
	nonLinearValuesLength = new uint32_t[nonLinearConstraints.size()+DValuesLength]();
	nonLinearLast = new uint32_t*[nonLinearConstraints.size()+DValuesLength]();
	combinedTimes = new uint16_t[nonLinearConstraints.size()+DValuesLength]();
	for (uint32_t i = 0; i < nonLinearConstraints.size()+DValuesLength; i++) combinedTimes[i] = 1; 
	nonLinearDistribution = new uint32_t*[nonLinearConstraints.size()+DValuesLength]();
	nonLinearDistributionLength = new uint32_t[nonLinearConstraints.size()+DValuesLength]();
	nonLinearDistributionBase = new uint32_t[nonLinearConstraints.size()+DValuesLength]();
	nonLinearPossibleKeys = new uint64_t[nonLinearConstraints.size()+DValuesLength]();
	cout << "Allocated " << nonLinearConstraints.size() << " units for nonlinear constraints" << endl;
	cout << "Total number of constraints: " << DValuesLength + SValuesLength + linearConstraints.size() + nonLinearConstraints.size() << endl;
	if (DValuesLength=SValuesLength+linearConstraints.size()+nonLinearConstraints.size() == 0) exit(0);
	// linearKeys.clear();
	// linearConstraints.clear();
	// printout(linearKeys.size());
	for (uint32_t i = 0; i < linearConstraints.size(); i++)
	{
		uint32_t *distribution; distribution = new uint32_t[1<<(1+SIZE)]();
		uint32_t distributionLength = 0;
		uint32_t distributionBase = 0;
		uint32_t valuesLength = 1;
		cout << "Resolving linear constraints[" << i << "]...";
		uint32_t keys = resolveLinear(linearConstraints[i],distribution,distributionLength,distributionBase);
		cout << "Completed!" << endl;
		linearDistribution[i] = distribution;
		linearDistributionLength[i] = distributionLength;
		linearDistributionBase[i] = distributionBase;
		linearKeys.push_back(keys);
	}

	// cout << "Adding in same round constraints into linearConstraints...";
	// for (uint32_t i = 0; i < SValuesLength; i++){
	// 	addLinearConstraint(singleKeys[i],SValues[i]);
	// }
	// cout << "Completed!" << endl;
	// for (int i = SValuesLength; i < 4*nr; i++) delete[] SValues[i];

	uint32_t i = 0;

	while (i < nonLinearConstraints.size())
	{
		uint32_t *values;
		uint32_t valuesLength;
		cout << "Resolving nonlinear constraints[" << i << "]...";
		nonLinearLast[i] = new uint32_t[1<<SIZE]();
		vector<uint32_t> keys;
		if (resolveNonLinear(nonLinearConstraints[i],values,valuesLength,nonLinearLast[i],keys))
		{
			nonLinearValues[i] = values;
			nonLinearValuesLength[i] = valuesLength;
			nonLinearKeys.push_back(keys);
			cout << "Completed!" << endl;
		}
		else
		{
			// obselete?
			cout << "Obselete code reached" << endl;
			int c = nonLinearConstraints[i][nonLinearConstraints[i].size()-1];
			#if SIZE == 4
			int DDT4[16][16] = {0};
			computeDDT(DDT4);
			keyProb -= log2((DDT4[(c>>SIZE)&ANDval][c&ANDval]+0.0)/(1<<SIZE)); // adding back the keyProb since we are ignoring the key constraint
			#elif SIZE == 8
			int DDT8[256][256] = {0};
			computeDDT8(DDT8);
			keyProb -= log2((DDT8[(c>>SIZE)&ANDval][c&ANDval]+0.0)/(1<<SIZE));
			#endif
			nonLinearConstraints.erase(nonLinearConstraints.begin()+i);
			i--;
			delete[] nonLinearLast[i];
		}
		i++;
	}
	cout << "Computing the nonlinear distribution..." << endl;
	computeNonLinearDistribution();
	
	cout << "Completed!" << endl;
	// Adding sameRoundConstraints into nonlienar
	// cout << "Adding in same round constraints into nonLinearConstraints..." << endl;
	// for (uint32_t i = 0; i < DValuesLength; i++){
	// 	addNonLinearConstraint(dualKeys[i],DValues[i]);
	// }
	// cout << "Completed!" << endl;
	if (TK_NUM >= 2)
	{

		cout << "Combining constraints...";
		combineConstraints(); // a step to minimize the number of constraints
		
		
		cout << "Completed" << endl;
		computeLinearPossibleKeys();
		computeNonLinearPossibleKeys();
		// next step
		cout << "Computing keys...";
		double keys = computeKeys();
		cout << "key space reduction: " << SIZE*16 - keys << endl;
		computeProbability();
	}
	else
	{
		cout << "Propagating the constraints...";
		propagateKeys();
		cout << "Completed" << endl;
		cout << "Combining constraints...";
		combineConstraints();
		cout << "Completed" << endl;
		computeLinearPossibleKeys();
		computeNonLinearPossibleKeys();
		cout << "Computing keys...";
		double keys = computeKeys();
		cout << "Completed!" << endl;
		cout << "key space reduction: " << SIZE*16 - keys << endl;
		computeProbability();
	}
}


int main(int argc, char **argv)
{
	if (argc == 1)
	{
		cout << "Please insert the trail number required" << endl;
		exit(0);
	}
	int Q = atoi(argv[1]); // I am not going to check it
	// SKINNY64
	#if SIZE == 4
	if (Q==0) Corrected_TK1_2020_1402(alpha,key_diff,nr,TK_NUM);
	if (Q==1) SK_2020_1402(alpha,key_diff,nr,TK_NUM);
	if (Q==2) TK1_2020_1402(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (Q==3) TK2_2020_1402(alpha,key_diff,nr,TK_NUM); 
	if (Q==4) TK3_2020_1402(alpha,key_diff,nr,TK_NUM);
	if (Q==5) get_4TK2_1_17_L_2020_1317(alpha,key_diff,nr,TK_NUM); // this is used as a first test
	if (Q==6) get_4TK2_1_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // without nonLinearConstraint[0], it's possible
	if (Q==7) get_4TK2_1_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	if (Q==8) get_4TK2_1_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // without nonLinearConstraint[0], it's possible
	if (Q==9) get_4TK2_2_17_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (Q==10) get_4TK2_2_17_U_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	if (Q==11) get_4TK2_2_18_L_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (Q==12) get_4TK2_2_19_U_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (Q==13) get_4TK3_1_22_L_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	if (Q==14) get_4TK3_1_22_U_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (Q==15) get_4TK3_1_22_L_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // ok
	if (Q==16) get_4TK3_1_22_U_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (Q==17) get_4TK3_1_22_U_BMD3_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (Q==18) get_4TK3_1_23_U_BMD1_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible in linear
	if (Q==19) get_4TK3_1_23_U_BMD2_2020_1317(alpha,key_diff,nr,TK_NUM); // impossible linear
	// get_4TK2_18_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // u = 6, v = 2,10,12,13 // all ok
	// get_4TK3_22_L_2021_656(alpha,key_diff,nr,TK_NUM,u,v); // u = 8,9, v = 5 // all ok


	// 8 bit trails
	// SKINNY128
	#elif SIZE == 8
	if (Q==1) TK1_8_2020_1402(alpha,key_diff,nr,TK_NUM); 
	if (Q==2) SK_tosc_2017_i4_99_129(alpha,key_diff,nr,TK_NUM);
	

	if (Q==3) TK2_8_2020_1402(alpha,key_diff,nr,TK_NUM);
	if (Q==4) get_8TK2_1_18_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==5) get_8TK2_1_18_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==6) get_8TK2_1_19_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==7) get_8TK2_1_20_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==8) get_8TK2_1_20_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==9) get_8TK2_1_21_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	
	if (Q==10) get_8TK2_2_18_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==11) get_8TK2_2_18_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==12) get_8TK2_2_19_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==13) get_8TK2_2_19_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==14) get_8TK2_2_20_U_2020_1317(alpha,key_diff,nr,TK_NUM);
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


	if (Q==15) get_8TK3_1_22_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==16) get_8TK3_1_22_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==17) get_8TK3_1_23_L_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==18) get_8TK3_1_23_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_24_L_2020_1317(alpha,key_diff,nr,TK_NUM); // repeated
	if (Q==19) get_8TK3_1_24_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	if (Q==20) get_8TK3_1_25_L_2020_1317(alpha,key_diff,nr,TK_NUM);
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