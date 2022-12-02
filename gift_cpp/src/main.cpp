#include <iostream>
#include "gift.h"
#include "trails.h"
#include "gaussianElimination.h"
#include "SAT.h"
#include "keyDependent.h"

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl

void testTrail(uint8_t alpha[30][16], int nr)
{
	vector<vector<vector<int>>> xPositions, yPositions;
	vector<vector<int>> newXConstraints, newYConstraints;
	vector<int> newXParity, newYParity;
	keyConstraints.clear();

	printState(alpha,nr);
	computeProb(alpha,nr);
	vector<vector<int>> yParity = findYConstraints(alpha,nr,yPositions);
	AddConstants(yPositions,yParity, nr);
	vector<vector<int>> xParity = findXConstraints(alpha,nr,xPositions);
	// for (int i = 0; i < yPositions.size(); i++)
	// {
	// 	for (int j = 0; j < yPositions[i].size(); j++)
	// 	{
	// 		for (int k = 0; k < yPositions[i][j].size(); k++)
	// 		{
	// 			cout << dec << (yPositions[i][j][k] >> 7) << "," << (yPositions[i][j][k] & 0b1111111) << " ";
	// 		}
	// 		cout << endl;
	// 	}
	// }
	// cout << "-------------------------------" << endl;
	// for (int i = 0; i < xPositions.size(); i++)
	// {
	// 	for (int j = 0; j < xPositions[i].size(); j++)
	// 	{
	// 		for (int k = 0; k < xPositions[i][j].size(); k++)
	// 		{
	// 			cout << dec << (xPositions[i][j][k] >> 7) << "," << (xPositions[i][j][k] & 0b1111111) << " ";
	// 		}
	// 		cout << endl;
	// 	}
	// }
	findLinearConstraints(xPositions,xParity,yPositions,yParity,nr);
	vector<vector<int>> masterKeyConstraints = propagateLinear(nr);
	computeDimension(masterKeyConstraints);
	cout << saved_prob << " number of constraints are satisfied automatically" << endl;
	solve(alpha,xPositions,xParity,yPositions,yParity,nr,SIZE);


}

int main(int argc, char **argv)
{
	if (argc == 1)
	{
		cout << "Please insert the trail number required" << endl;
		exit(0);
	}
	int Q = atoi(argv[1]); // I am not going to check it
	uint8_t alpha[30][16] = {0};
	uint32_t nr;
	saved_prob = 0;
	#if SIZE == 4
	if (Q==0) SK_4_2021_1179_1(alpha,nr);
	if (Q==1) SK_4_2021_1179_2(alpha,nr);
	if (Q==2) SK_4_2021_1179_3(alpha,nr);
	if (Q==3) SK_4_2019_49_9(alpha,nr);
	if (Q==4) SK_4_2019_49_12(alpha,nr);
	if (Q==5) SK_4_2019_49_13(alpha,nr);
	if (Q==6) SK_4_2018_390_Table4(alpha,nr);
	if (Q==7) SK_4_2018_390_Table6(alpha,nr);
	// randomTrail4(alpha,nr);
	// nonLinearTrail4(alpha,nr);
	#elif SIZE == 8
	if (Q==0) SK_8_12_2019_25(alpha,nr);
	if (Q==1) SK_8_13_2019_25(alpha,nr);
	if (Q==2) SK_8_21_2019_25(alpha,nr);
	if (Q==3) SK_8_2019_49_21(alpha,nr);
	if (Q==4) SK_8_2018_390_Table10(alpha,nr); // impossible
	if (Q==5) SK_8_2018_390_Table15(alpha,nr);
	if (Q==6) SK_8_2018_390_Table16(alpha,nr);
	#endif
	testTrail(alpha,nr);
	return 0;
}