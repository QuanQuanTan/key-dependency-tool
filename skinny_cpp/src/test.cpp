#include "trails.h"
#include "skinny.h"
#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

#if !defined SIZE
#error SIZE is not defined.
#endif


mt19937 rand_generator;

void printDiff(uint8_t state0[4][4], uint8_t state1[4][4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << hex << setw(2) << setfill('0') << int(state0[i][j] ^ state1[i][j]) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void printState(uint8_t state0[4][4],uint8_t state1[4][4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << hex << setw(2) << setfill('0') << int(state0[i][j]) << " ";
		}
		cout << " | ";
		for (int j = 0; j < 4; j++)
		{
			cout << hex << setw(2) << setfill('0') << int(state1[i][j]) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void test_trail(uint8_t alpha[20][4][4],uint8_t key_diff[4][4][4], uint32_t nr, int TK_NUM, uint64_t numTrails,uint64_t numKeys, int start)
{

	uint8_t masterKey0[4][4][4]={0}, masterKey1[4][4][4]={0}, currentKey0[4][4][4]={0}, currentKey1[4][4][4]={0};
	uint8_t state0[4][4], state1[4][4];
	uniform_int_distribution<uint64_t> rand_distribution(0,0xffffffffffffffff);	
	for (uint64_t q1 = 0; q1 < numKeys; q1++)
	{
		int count = 0;
		for (int i = 0; i < TK_NUM; i++)
		{
			#if SIZE == 4
			uint64_t rand_num = rand_distribution(rand_generator);
			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					masterKey0[i][j][k] = rand_num & 0xf;
					masterKey1[i][j][k] = masterKey0[i][j][k] ^ key_diff[i][j][k];
					rand_num = rand_num >> 4;
				}
			}
			#elif SIZE == 8
			uint64_t rand_num0 = rand_distribution(rand_generator);
			uint64_t rand_num1 = rand_distribution(rand_generator);
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					masterKey0[i][j][k] = rand_num0 & 0xff;
					masterKey1[i][j][k] = masterKey0[i][j][k] ^ key_diff[i][j][k];
					rand_num0 = rand_num0 >> 8;
				}
			}
			for (int j = 2; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					masterKey0[i][j][k] = rand_num1 & 0xff;
					masterKey1[i][j][k] = masterKey0[i][j][k] ^ key_diff[i][j][k];
					rand_num0 = rand_num1 >> 8;
				}
			}
			#endif
		}
		// int w = 2;
		// cout << "here" << endl;
		for (uint64_t q2 = 0; q2 < numTrails; q2++)
		{
			// if (log2(q2) > w)
			// {
			// 	cout << log2(q2) << endl;
			// 	w++;
			// }
			for (int i = 0; i < TK_NUM; i++){
				for (int j = 0; j < 4; j++){
					for (int k = 0; k < 4; k++){
						currentKey0[i][j][k] = masterKey0[i][j][k];
						currentKey1[i][j][k] = masterKey1[i][j][k];
					}
				}
			}
			for (int i = 0; i < start; i++)
			{
				#if SIZE == 4
				key_schedule_round_function(currentKey0);
				key_schedule_round_function(currentKey1);
				#elif SIZE == 8
				key_schedule_round_function8(currentKey0);
				key_schedule_round_function8(currentKey1);
				#endif
			}
			#if SIZE == 4
			uint64_t rand_num = rand_distribution(rand_generator);
			for (int i = 0; i < 4; i++) // randomizing the plaintext
			{
				for (int j = 0; j < 4; j++) 
				{
					state0[i][j] = rand_num & 0xf;
					state1[i][j] = state0[i][j] ^ alpha[start][i][j];
					rand_num = rand_num >> 4;
				}
			}
			#elif SIZE == 8
			uint64_t rand_num0 = rand_distribution(rand_generator);
			uint64_t rand_num1 = rand_distribution(rand_generator);
			for (int i = 0; i < 2; i++) // randomizing the plaintext
			{
				for (int j = 0; j < 4; j++) 
				{
					state0[i][j] = rand_num0 & 0xff;
					state1[i][j] = state0[i][j] ^ alpha[start][i][j];
					rand_num0 = rand_num0 >> 8;
				}
			}
			for (int i = 2; i < 4; i++) // randomizing the plaintext
			{
				for (int j = 0; j < 4; j++) 
				{
					state0[i][j] = rand_num1 & 0xff;
					state1[i][j] = state0[i][j] ^ alpha[start][i][j];
					rand_num1 = rand_num1 >> 8;
				}
			}
			#endif
			bool flag = false;
			for (int r = start; r < start+nr; r++)
			{
				#if SIZE == 4
				Substitution(state0);
				Substitution(state1);
				#elif SIZE == 8
				Substitution8(state0);
				Substitution8(state1);
				#endif
				// if (r==1) printState(state0,state1);
				// cout << "r: " << r << endl;
				// if (r==2) printDiff(state0,state1);
				AddConstant(state0,r);
				AddConstant(state1,r);

				// if (r==1) printDiff(state0,state1);
				AddRoundTweakey(state0,currentKey0);
				AddRoundTweakey(state1,currentKey1);
				// if (r==2) printDiff(state0,state1);
				// if (r==1) printDiff(state0,state1);
				SR(state0);
				SR(state1);
				// if (r==1) printDiff(state0,state1);
				MC(state0);
				MC(state1);
				#if SIZE == 4
				key_schedule_round_function(currentKey0);
				key_schedule_round_function(currentKey1);
				#elif SIZE == 8
				key_schedule_round_function8(currentKey0);
				key_schedule_round_function8(currentKey1);
				#endif
				// if (r==1) printDiff(state0,state1);
				// cout << "---------------------" << endl;
				// if (r==2) printDiff(state0,state1);
				// cout << "---------------------" << endl;
				flag = true;
				for (int i = 0; i < 4; i++){
					for (int j = 0; j < 4; j++){
						if ((state0[i][j] ^ state1[i][j]) != alpha[r+1][i][j]) flag = false;
					}
				}
				if (flag == false) {
					break;
				}
			}
			if (flag == true) count++;
		}
		cout << dec << "key #" << q1 << ": " << -log2((count+0.0)/numTrails) << endl;
	}
}




int main(){	
	uint8_t alpha[20][4][4] = {0};
	uint8_t key_diff[4][4][4] = {0};
	uint32_t nr;
	uint32_t trail_rounds;
	int TK_NUM;
	uint64_t numKeys = 1<<6;
	uint64_t numTrails = 1ULL<<23;

	int start = 0x0;
	// SK_2020_1402(alpha,key_diff,nr,TK_NUM);
	// TK1_2020_1402(alpha, key_diff ,nr, TK_NUM);
	// get_4TK2_2_17_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// get_8TK3_1_23_U_2020_1317(alpha,key_diff,nr,TK_NUM);
	// Corrected_TK1_2020_1402(alpha, key_diff,trail_rounds, TK_NUM);
	get_4TK3_22_L_2021_656(alpha,key_diff,nr,TK_NUM,8,5);
	// nr = 3;
	test_trail(alpha,key_diff,nr,TK_NUM,numTrails,numKeys,start);
	return 0;
}
