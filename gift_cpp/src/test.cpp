#include <iostream>
#include <random>
#include "trails.h"
#include "gift.h"

#define printout(x) cout << #x << ": " << x << " at LINE " << __LINE__ << endl

mt19937 rand_generator;

void test(uint8_t alpha[30][16], int nr, uint64_t numTrails, uint64_t numKeys, int start)
{
	uniform_int_distribution<uint64_t> rand_distribution(0,0xffffffffffffffff);	
	uint16_t masterKey[8];
	uint16_t copyKey[8];
	uint8_t state0[16], state1[16];
	for (uint64_t q1 = 0; q1 < numKeys; q1++)
	{
		int count = 0;
		// establish the key
		uint64_t rand_num = rand_distribution(rand_generator);	
		for (int i = 0; i < 4; i++){
			masterKey[i] = rand_num & 0xffff;
			rand_num = rand_num >> 16;
		}
		rand_num = rand_distribution(rand_generator);
		for (int i = 4; i < 8; i++){
			masterKey[i] = rand_num & 0xffff;
			rand_num = rand_num >> 16;
		}
		for (uint64_t q2 = 0; q2 < numTrails; q2++)
		{
			// randomize the state
			#if SIZE == 4
			rand_num = rand_distribution(rand_generator);	
			for (int i = 0; i < 16; i++)
			{
				state0[i] = rand_num & 0xf;
				state1[i] = state0[i] ^ alpha[start][i];
				rand_num = rand_num >> 4;
			}
			#elif SIZE == 8
			rand_num = rand_distribution(rand_generator);	
			for (int i = 0; i < 8; i++){
				state0[i] = rand_num & 0xff;
				state1[i] = state0[i] ^ alpha[start][i];
				rand_num = rand_num >> 8;
			}
			rand_num = rand_distribution(rand_generator);	
			for (int i = 8; i < 16; i++){
				state0[i] = rand_num & 0xff;
				state1[i] = state0[i] ^ alpha[start][i];
				rand_num = rand_num >> 8;
			}
			#endif
			// copy the masterKey out
			for (int i = 0; i < 8; i++) copyKey[i] = masterKey[i];
			bool flag = false;
			for (int n = start; n < start + nr; n++)
			{
				RoundFunction(state0,&copyKey[6],nr);
				RoundFunction(state1,&copyKey[6],nr);
				keyScheduleRoundFunction(copyKey);
				flag = true;
				for (int i = 0; i < 16; i++){
					if ((state0[i] ^ state1[i]) != alpha[n+1][i]) flag = false;
				}
				if (!flag) break;
			}
			if (flag) count++;
		}
		cout << dec << "key #" << q1 << ": " << -log2((count+0.0)/numTrails) << endl;
	}	
}




int main()
{
	uint8_t alpha[30][16] = {0};
	uint32_t nr;
	uint64_t numKeys = 1<<9;
	uint64_t numTrails = 1ULL<<16;
	// SK_4_2021_1179_1(alpha,nr);
	// SK_4_2021_1179_2(alpha,nr);
	// SK_4_2021_1179_3(alpha,nr);
	// SK_4_2019_49_9(alpha,nr);
	SK_4_2019_49_12(alpha,nr);
	// SK_4_2019_49_13(alpha,nr);
	// SK_4_2018_390_Table4(alpha,nr);
	// SK_4_2018_390_Table6(alpha,nr);
	// randomTrail4(alpha,nr);
	// nonLinearTrail4(alpha,nr);
	// SK_8_12_2019_25(alpha,nr);
	// SK_8_13_2019_25(alpha,nr);
	// SK_8_21_2019_25(alpha,nr);
	// SK_8_2019_49_21(alpha,nr);
	// SK_8_2018_390_Table10(alpha,nr); // impossible
	// SK_8_2018_390_Table15(alpha,nr);
	// SK_8_2018_390_Table16(alpha,nr);
	uint32_t n = 3;
	uint32_t start = 9;

	test(alpha,n,numTrails,numKeys,start);

	return 0;
}