#include "deoxys.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <string>

using namespace std;

static mt19937 rand_generator;

static void conversion(uint8_t alpha[20][4][4],uint32_t trail_rounds,uint8_t key_diff[3][4][4])
{
	uint8_t key_diff_temp[3][4][4];
	for (int i = 0; i < 3; i++){ for (int j = 0; j < 4; j++) { for (int k = 0; k < 4; k++) key_diff_temp[i][j][k] = key_diff[i][j][k];}}
	for (int i = 0; i < trail_rounds-1; i++){ keyScheduleRoundFunction(key_diff_temp);}
	SR(alpha[trail_rounds]);
	MC(alpha[trail_rounds]);
}

static bool sanityCheck(uint8_t alpha[20][4][4],uint32_t trail_rounds,uint8_t key_diff[3][4][4])
{
	uint32_t DDT[256][256] = {0};
	computeDDT(DDT);

	uint8_t key_diff_temp[4][4][4];
	uint8_t state_temp_before[4][4];
	uint8_t state_temp_after[4][4];
	for (int i = 0; i < 4; i++){ for (int j = 0; j < 4; j++) { for (int k = 0; k < 4; k++) key_diff_temp[i][j][k] = key_diff[i][j][k];}}
	for (int n = 0; n < trail_rounds; n++)
	{
		for (int r = 0; r < 4; r++){ for (int c = 0; c < 4; c++){ state_temp_before[r][c] = alpha[n][r][c]; }}
		for (int r = 0; r < 4; r++){ for (int c = 0; c < 4; c++){ state_temp_after[r][c] = alpha[n+1][r][c]; }}
		AddRoundTweakey(state_temp_before,key_diff_temp);
		invMC(state_temp_after); invSR(state_temp_after); 
		// for (int r = 0; r < 4; r++)
		// {
		// 	for (int c = 0; c < 4; c++)
		// 	{
		// 		std::cout << std::setw(2) << std::hex << int(state_temp_before[r][c]) << " ";
		// 	}
		// 	std::cout << "  ";
		// 	for (int c = 0; c < 4; c++)
		// 	{
		// 		std::cout << std::setw(2) << std::hex << int(state_temp_after[r][c]) << " ";
		// 	}
		// 	std::cout << "  ";
		// 	for (int c = 0; c < 4; c++)
		// 	{
		// 		std::cout << std::setw(2) << std::hex << int(key_diff_temp[0][r][c]) << " ";
		// 	}
		// 	std::cout << "  ";
		// 	for (int c = 0; c < 4; c++)
		// 	{
		// 		std::cout << std::setw(2) << std::hex << int(key_diff_temp[1][r][c]) << " ";
		// 	}
		// 	std::cout << std::endl;
		// }
		// std::cout << std::endl;
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) 
			{
				if (DDT[state_temp_before[r][c]][state_temp_after[r][c]] == 0) 
				{
					std::cout << "Error in round: " << n << ", pos: " << 4*r+c << ", input: " << int(state_temp_before[r][c]) << ", output: " << int(state_temp_after[r][c]) << std::endl;
					return false;
				}
			}
		}
		keyScheduleRoundFunction(key_diff_temp);
	}
	return true;
}

void generateBestTrail(uint8_t alpha[20][4][4], uint8_t key_diff[3][4][4], uint8_t a, uint8_t pos, uint32_t &trail_rounds, int &TK_NUM)
{
	TK_NUM = 1;
	trail_rounds = 4;
	uniform_int_distribution<uint8_t> rng8(0, 0xff);
	uint32_t DDT[256][256] = {0};
	computeDDT(DDT);
	uint8_t rand_num;
	// a is the value of the active byte at the end of the second round
	// pos is the position of the active cell (0 - 15)
	uint8_t tmp[4][4] = {0}; tmp[pos/4][pos%4] = a;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) alpha[2][r][c] = tmp[r][c];
	}
	// backwards
	for (int n = 1; n >= 0; n--)
	{
		invMC(tmp);
		invSR(tmp);
		// invSubstitution (random for now?)
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				if (tmp[r][c] != 0)
				{
					while (true)
					{
						rand_num = rng8(rand_generator);
						if (DDT[rand_num][tmp[r][c]] != 0)
						{
							tmp[r][c] = rand_num;
							break;
						}
					}
				}
			}
		}
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) alpha[n][r][c] = tmp[r][c];
		}
	}
	// forward 

	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) tmp[r][c] = alpha[2][r][c];
	}
	for (int n = 3; n < 5; n++)
	{
		// Substitution
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				if (tmp[r][c] != 0)
				{
					while (true)
					{
						rand_num = rng8(rand_generator);
						if (DDT[tmp[r][c]][rand_num] != 0)
						{
							tmp[r][c] = rand_num;
							break;
						}
					}
				}
			}
		}
		// SR & MC
		SR(tmp);
		MC(tmp);
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) alpha[n][r][c] = tmp[r][c];
		}
	}
	
	// empty key diff
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++) key_diff[i][j][k] = 0;
		}
	}
	if (!sanityCheck(alpha,trail_rounds,key_diff))
	{
		std::cout << "Insanity!" << std::endl;
		exit(0);
	}
	
}

void TK2_8U_2017_693(uint8_t alpha[20][4][4], uint8_t key_diff[3][4][4], uint32_t& trail_rounds, int &TK_NUM)
{
	TK_NUM = 2;
	trail_rounds = 4;
	uint8_t alpha_before_temp[20][4][4] = {\
	{{0x00,0xb9,0x00,0x00},{0x00,0x00,0xd1,0x00},{0x00,0x00,0x00,0xab},{0x61,0x00,0x00,0x97}},
	{{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0xe5,0x00,0x00},{0x00,0x00,0x00,0x00}},
	{{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00}},
	{{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00}}};

	uint8_t alpha_after_temp[20][4][4] = {\
	{{0x00,0xb9,0x00,0x00},{0x00,0x00,0xd1,0x00},{0x00,0x00,0x00,0xab},{0x61,0x00,0x00,0x00}}, // these are trash lines
	{{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00}}, // these are trash lines
	{{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00}}, // these are trash lines
	{{0x02,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00}}};
	uint8_t key_diff_temp[3][4][4] = {
	{{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x46}},
	{{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0xd1}},
	{{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00},{0x00,0x00,0x00,0x00}}};
	for (int r = 0; r < 20; r++) {for (int i = 0; i < 4; i++) { for (int j = 0; j < 4; j++) {alpha[r][i][j] = alpha_before_temp[r][i][j];}}}
	for (int i = 0; i < 4; i++) { for (int j = 0; j < 4; j++) {alpha[trail_rounds][i][j] = alpha_after_temp[trail_rounds-1][i][j];}}
	for (int i = 0; i < 3; i++){ for (int j = 0; j < 4; j++) { for (int k = 0; k < 4; k++) key_diff[i][j][k] = key_diff_temp[i][j][k];}}
	conversion(alpha,trail_rounds,key_diff);
	if (!sanityCheck(alpha,trail_rounds,key_diff))
	{
		std::cout << "Insanity!" << std::endl;
		exit(0);
	}
}

void create2Fail_trail(uint8_t alpha[20][4][4], uint8_t key_diff[3][4][4], uint32_t& trail_rounds, int &TK_NUM)
{
	uint8_t rand_num;
	uniform_int_distribution<uint8_t> rng8(0, 0xff);
	uint8_t tmp[4][4] = {0};
	uint32_t DDT[256][256] = {0};
	computeDDT(DDT);
	
	tmp[0][0] = rng8(rand_generator);
	MC(tmp);
	tmp[3][1] = rng8(rand_generator);
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) alpha[2][r][c] = tmp[r][c];
	}
	
	// forward
	// Substitution
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (tmp[r][c] != 0)
			{
				while (true)
				{
					rand_num = rng8(rand_generator);
					if (DDT[tmp[r][c]][rand_num] != 0)
					{
						tmp[r][c] = rand_num;
						break;
					}
				}
			}
		}
	}
	SR(tmp);
	MC(tmp);
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) alpha[3][r][c] = tmp[r][c];
	}

	// backward
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) tmp[r][c] = alpha[2][r][c];
	}
	for (int n = 1; n >= 0; n--)
	{
		invMC(tmp);
		invSR(tmp);
		// inverseSub
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				if (tmp[r][c] != 0)
				{
					while (true)
					{
						rand_num = rng8(rand_generator);
						if (DDT[rand_num][tmp[r][c]] != 0)
						{
							tmp[r][c] = rand_num;
							break;
						}
					}
				}
			}
		}
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) alpha[n][r][c] = tmp[r][c];
		}
	}
	trail_rounds = 3;
	TK_NUM = 2;
	if (!sanityCheck(alpha,trail_rounds,key_diff))
	{
		std::cout << "Insanity!" << std::endl;
		exit(0);
	}
}