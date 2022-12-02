#include "gift.h"

using namespace std;



uint8_t Sbox[16] = {0x1, 0xa, 0x4, 0xc, 0x6, 0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe};
uint8_t perm[64] = {0,17,34,51,48,1,18,35,32,49,2,19,16,33,50,3,4,21,38,55,52,5,22,39,36,53,6,23,20,37,54,
				7,8,25,42,59,56,9,26,43,40,57,10,27,24,41,58,11,12,29,46,63,60,13,30,47,44,61,14,31,28,45,62,15};
uint8_t perm8[128] = {0,33,66,99,96,1,34,67,64,97,2,35,32,65,98,3,4,37,70,103,100,5,38,71,68,101,6,39,36,69,102,
					7,8,41,74,107,104,9,42,75,72,105,10,43,40,73,106,11,12,45,78,111,108,13,46,79,76,109,14,47,
					44,77,110,15,16,49,82,115,112,17,50,83,80,113,18,51,48,81,114,19,20,53,86,119,116,21,54,87,
					84,117,22,55,52,85,118,23,24,57,90,123,120,25,58,91,88,121,26,59,56,89,122,27,28,61,94,127,
					124,29,62,95,92,125,30,63,60,93,126,31};
uint8_t Constants[48] = {0x01,0x03,0x07,0x0F,0x1F,0x3E,0x3D,0x3B,0x37,0x2F,0x1E,0x3C,0x39,0x33,0x27,0x0E,
	                     0x1D,0x3A,0x35,0x2B,0x16,0x2C,0x18,0x30,0x21,0x02,0x05,0x0B,0x17,0x2E,0x1C,0x38,
	                     0x31,0x23,0x06,0x0D,0x1B,0x36,0x2D,0x1A,0x34,0x29,0x12,0x24,0x08,0x11,0x22,0x04};
uint8_t invSbox[16] = {0xd,0x0,0x8,0x6,0x2,0xc,0x4,0xb,0xe,0x7,0x1,0xa,0x3,0x9,0xf,0x5};
uint8_t invPerm[64] = {0,5,10,15,16,21,26,31,32,37,42,47,48,53,58,63,12,1,6,11,28,17,22,27,44,33,38,43,60,49,54,59,
						8,13,2,7,24,29,18,23,40,45,34,39,56,61,50,55,4,9,14,3,20,25,30,19,36,41,46,35,52,57,62,51};
uint8_t invPerm8[128] = {0,5,10,15,16,21,26,31,32,37,42,47,48,53,58,63,64,69,74,79,80,85,90,95,96,101,106,111,112,117,
						122,127,12,1,6,11,28,17,22,27,44,33,38,43,60,49,54,59,76,65,70,75,92,81,86,91,108,97,102,107,
						124,113,118,123,8,13,2,7,24,29,18,23,40,45,34,39,56,61,50,55,72,77,66,71,88,93,82,87,104,109,
						98,103,120,125,114,119,4,9,14,3,20,25,30,19,36,41,46,35,52,57,62,51,68,73,78,67,84,89,94,83,
						100,105,110,99,116,121,126,115};

uint8_t getSbox(uint8_t val) {return Sbox[val];}
uint8_t getInvSbox(uint8_t val) {return invSbox[val];}
uint8_t getConstants(uint8_t val) {return Constants[val];}
uint8_t getPerm(uint8_t val) {return perm[val];}
uint8_t getPerm8(uint8_t val) {return perm8[val];}
uint8_t getinvPerm(uint8_t val) {return invPerm[val];}
uint8_t getinvPerm8(uint8_t val) {return invPerm8[val];}

void Substitution(uint8_t state[16])
{
	for (int i = 0; i < 16; i++) state[i] = Sbox[state[i]];
}


void invSubstitution(uint8_t state[16])
{
	for (int i = 0; i < 16; i++) state[i] = invSbox[state[i]];
}


void Permutation(uint8_t state[16])
{
	uint64_t s1 = 0, s2 = 0;
	for (int i = 0; i < 16; i++)
	{
		s1 = s1 << 4;
		s1 ^= state[i];
	}
	for (int pos = 0; pos < 64; pos++) s2 ^= ((s1 >> pos) & 0b1) << perm[pos];
	for (int i = 15; i >= 0; i--)
	{
		state[i] = s2 & 0xf;
		s2 = s2 >> 4;
	}
}

void invPermutation(uint8_t state[16])
{
	uint64_t s1 = 0, s2 = 0;
	for (int i = 0; i < 16; i++)
	{
		s1 = s1 << 4;
		s1 ^= state[i];
	}
	for (int pos = 0; pos < 64; pos++) s2 ^= ((s1 >> pos) & 0b1) << invPerm[pos];
	for (int i = 15; i >= 0; i--)
	{
		state[i] = s2 & 0xf;
		s2 = s2 >> 4;
	}
}

void AddRoundKey(uint8_t state[16], uint16_t key[2])
{
	uint16_t ukey = key[0];
	uint16_t vkey = key[1];
	for (int i = 15; i >= 0; i--)
	{
		state[i] = state[i] ^ (((ukey >> (15-i)) & 0b1) << 1);
		state[i] = state[i] ^ ((vkey >> (15-i)) & 0b1);
	}
}

void AddConstant(uint8_t state[16], uint8_t constant)
{
	state[0] ^= (1 << 3);
	state[10] ^= ((constant >> 5) & 0b1) << 3;
	state[11] ^= ((constant >> 4) & 0b1) << 3;
	state[12] ^= ((constant >> 3) & 0b1) << 3;
	state[13] ^= ((constant >> 2) & 0b1) << 3;
	state[14] ^= ((constant >> 1) & 0b1) << 3;
	state[15] ^= ((constant >> 0) & 0b1) << 3;
}

void keyScheduleRoundFunction(uint16_t masterKey[8])
{
	uint16_t tmp[8] = {0};
	tmp[0] = (masterKey[6] >> 2) | ((masterKey[6] << (16-2)) & 0xffff);
	tmp[1] = (masterKey[7] >> 12) | ((masterKey[7] << (16-12)) & 0xffff);
	for (int i = 2; i < 8; i++) tmp[i] = masterKey[i-2];
	for (int i = 0; i < 8; i++) masterKey[i] = tmp[i];
}

void RoundFunction(uint8_t state[16], uint16_t key[2], int nr)
{
	Substitution(state);
	Permutation(state);
	AddRoundKey(state,key);
	AddConstant(state,Constants[nr]);
}

void invRoundFunction(uint8_t state[16], uint16_t key[2], int nr)
{
	AddConstant(state,Constants[nr]);
	AddRoundKey(state,key);
	invPermutation(state);
	invSubstitution(state);
}

void gift_encrypt(uint8_t state[16], uint16_t masterKey[8])
{
	for (int nr = 0; nr < 28; nr++)
	{
		RoundFunction(state,&masterKey[6],nr);
		keyScheduleRoundFunction(masterKey);
	}
}

void gift_decrypt(uint8_t state[16], uint16_t masterKey[8])
{
	uint16_t expandedKey[40][8] = {0};
	for (int nr = 0; nr < 28; nr++)
	{
		for (int i = 0; i < 8; i++) expandedKey[nr][i] = masterKey[i];
		keyScheduleRoundFunction(masterKey);
	}
	for (int nr = 27; nr >= 0; nr--)
	{
		invRoundFunction(state,&expandedKey[nr][6],nr);
	}
}

void test_vectors()
{
	uint8_t state11[16] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
	uint8_t state12[16] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
	uint16_t masterKey11[8] = {0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000};
	uint16_t masterKey12[8] = {0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000};
	gift_encrypt(state11,masterKey11);
	uint8_t expected1[16] = {0xf,0x6,0x2,0xb,0xc,0x3,0xe,0xf,0x3,0x4,0xf,0x7,0x7,0x5,0xa,0xc};
	for (int i = 0; i < 16; i++)
	{
		if (state11[i] != expected1[i])
		{
			cout << "GIFT-64 Test vector 1 (enc) failed" << endl;
			exit(0);
		}
	}
	gift_decrypt(state11,masterKey12);
	for (int i = 0; i < 16; i++)
	{
		if (state11[i] != state12[i])
		{
			cout << "GIFT-64 Test vector 1 (dec) failed" << endl;
			exit(0);
		}
	}

	uint8_t state21[16] = {0xf,0xe,0xd,0xc,0xb,0xa,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0};
	uint8_t state22[16] = {0xf,0xe,0xd,0xc,0xb,0xa,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0};
	uint16_t masterKey21[8] = {0xfedc,0xba98,0x7654,0x3210,0xfedc,0xba98,0x7654,0x3210};
	uint16_t masterKey22[8] = {0xfedc,0xba98,0x7654,0x3210,0xfedc,0xba98,0x7654,0x3210};
	uint8_t expected2[16] = {0xc,0x1,0xb,0x7,0x1,0xf,0x6,0x6,0x1,0x6,0x0,0xf,0xf,0x5,0x8,0x7};
	gift_encrypt(state21,masterKey21);
	for (int i = 0; i < 16; i++)
	{
		if (state21[i] != expected2[i])
		{
			cout << "GIFT-64 Test vector 2 (enc) failed" << endl;
			exit(0);
		}
	}
	gift_decrypt(state21,masterKey22);
	for (int i = 0; i < 16; i++)
	{
		if (state21[i] != state22[i])
		{
			cout << "GIFT-64 Test vector 2 (dec) failed" << endl;
			exit(0);
		}
	}

	uint8_t state31[16] = {0xc,0x4,0x5,0x0,0xc,0x7,0x7,0x2,0x7,0xa,0x9,0xb,0x8,0xa,0x7,0xd};
	uint8_t state32[16] = {0xc,0x4,0x5,0x0,0xc,0x7,0x7,0x2,0x7,0xa,0x9,0xb,0x8,0xa,0x7,0xd};
	uint16_t masterKey31[8] = {0xbd91,0x731e,0xb6bc,0x2713,0xa1f9,0xf6ff,0xc750,0x44e7};
	uint16_t masterKey32[8] = {0xbd91,0x731e,0xb6bc,0x2713,0xa1f9,0xf6ff,0xc750,0x44e7};
	uint8_t expected3[16] = {0xe,0x3,0x2,0x7,0x2,0x8,0x8,0x5,0xf,0xa,0x9,0x4,0xb,0xa,0x8,0xb};
	gift_encrypt(state31,masterKey31);
	for (int i = 0; i < 16; i++)
	{
		if (state31[i] != expected3[i])
		{
			cout << "GIFT-64 Test vector 3 (enc) failed" << endl;
			exit(0);
		}
	}
	gift_decrypt(state31,masterKey32);
	for (int i = 0; i < 16; i++)
	{
		if (state31[i] != state32[i])
		{
			cout << "GIFT-64 Test vector 3 (dec) failed" << endl;
			exit(0);
		}
	}
}


void Substitution8(uint8_t state[16])
{
	for (int i = 0; i < 16; i++)
	{
		state[i] = (Sbox[state[i] >> 4] << 4) ^ Sbox[state[i] & 0xf];
	}
}
void invSubstitution8(uint8_t state[16])
{
	for (int i = 0; i < 16; i++)
	{
		state[i] = (invSbox[state[i] >> 4] << 4) ^ invSbox[state[i] & 0xf];
	}
}

void Permutation8(uint8_t state[16])
{
	uint64_t s1[2] = {0,0}, s2[2] = {0,0};
	for (int i = 0; i < 8; i++)
	{
		s1[0] = s1[0] << 8;
		s1[0] ^= state[i];
		s1[1] = s1[1] << 8;
		s1[1] ^= state[8+i];
	}
	// settle the position here!
	for (int pos = 64; pos < 128; pos++)
	{
		if (perm8[pos] < 64) s2[1] ^= ((s1[0] >> (pos-64)) & 0b1) << perm8[pos];
		else s2[0] ^= ((s1[0] >> (pos-64)) & 0b1) << (perm8[pos]-64);
	}
	for (int pos = 0; pos < 64; pos++)
	{
		if (perm8[pos] < 64) s2[1] ^= ((s1[1] >> pos) & 0b1) << perm8[pos];
		else s2[0] ^= ((s1[1] >> pos) & 0b1) << (perm8[pos]-64);
	}
	for (int i = 15; i >= 8; i--)
	{
		state[i] = s2[1] & 0xff;
		s2[1] = s2[1] >> 8;
	}
	for (int i = 7; i >= 0; i--)
	{
		state[i] = s2[0] & 0xff;
		s2[0] = s2[0] >> 8;
	}
}

void invPermutation8(uint8_t state[16])
{
	uint64_t s1[2] = {0,0}, s2[2] = {0,0};
	for (int i = 0; i < 8; i++)
	{
		s1[0] = s1[0] << 8;
		s1[0] ^= state[i];
		s1[1] = s1[1] << 8;
		s1[1] ^= state[8+i];
	}
	// settle the position here!
	for (int pos = 64; pos < 128; pos++)
	{
		if (invPerm8[pos] < 64) s2[1] ^= ((s1[0] >> (pos-64)) & 0b1) << invPerm8[pos];
		else s2[0] ^= ((s1[0] >> (pos-64)) & 0b1) << (invPerm8[pos]-64);
	}
	for (int pos = 0; pos < 64; pos++)
	{
		if (invPerm8[pos] < 64) s2[1] ^= ((s1[1] >> pos) & 0b1) << invPerm8[pos];
		else s2[0] ^= ((s1[1] >> pos) & 0b1) << (invPerm8[pos]-64);
	}
	for (int i = 15; i >= 8; i--)
	{
		state[i] = s2[1] & 0xff;
		s2[1] = s2[1] >> 8;
	}
	for (int i = 7; i >= 0; i--)
	{
		state[i] = s2[0] & 0xff;
		s2[0] = s2[0] >> 8;
	}
}

void AddRoundKey8(uint8_t state[16], uint16_t key[4])
{
	uint32_t ukey = (key[0] << 16) ^ key[1];
	uint32_t vkey = (key[2] << 16) ^ key[3];
	for (int i = 15; i >= 0; i--)
	{
		state[i] = state[i] ^ (((ukey >> (31-2*i)) & 0b1) << 6);
		state[i] = state[i] ^ (((vkey >> (31-2*i)) & 0b1) << 5);
		state[i] = state[i] ^ (((ukey >> (31-2*i-1)) & 0b1) << 2);
		state[i] = state[i] ^ (((vkey >> (31-2*i-1)) & 0b1) << 1);
	}
}

void AddConstant8(uint8_t state[16], uint8_t constant)
{
	state[0] ^= (1 << 7);
	state[13] ^= ((constant >> 5) & 0b1) << 7;
	state[13] ^= ((constant >> 4) & 0b1) << 3;
	state[14] ^= ((constant >> 3) & 0b1) << 7;
	state[14] ^= ((constant >> 2) & 0b1) << 3;
	state[15] ^= ((constant >> 1) & 0b1) << 7;
	state[15] ^= ((constant >> 0) & 0b1) << 3;
}

void RoundFunction8(uint8_t state[16], uint16_t key[4], int nr)
{
	Substitution8(state);
	Permutation8(state);
	AddRoundKey8(state,key);
	AddConstant8(state,Constants[nr]);
}

void invRoundFunction8(uint8_t state[16], uint16_t key[4], int nr)
{
	AddConstant8(state,Constants[nr]);
	AddRoundKey8(state,key);
	invPermutation8(state);
	invSubstitution8(state);
}

void gift_encrypt8(uint8_t state[16], uint16_t masterKey[8])
{
	uint16_t requiredKey[4] = {0};
	for (int nr = 0; nr < 40; nr++)
	{
		requiredKey[0] = masterKey[2];
		requiredKey[1] = masterKey[3];
		requiredKey[2] = masterKey[6]; 
		requiredKey[3] = masterKey[7];
		RoundFunction8(state,requiredKey,nr);
		keyScheduleRoundFunction(masterKey);
	}
}

void gift_decrypt8(uint8_t state[16], uint16_t masterKey[8])
{
	uint16_t expandedKey[40][16] = {0};
	uint16_t requiredKey[4] = {0};
	for (int nr = 0; nr < 40; nr++)
	{
		for (int i = 0; i < 16; i++) expandedKey[nr][i] = masterKey[i];
		keyScheduleRoundFunction(masterKey);
	}
	for (int nr = 39; nr >= 0; nr--)
	{
		requiredKey[0] = expandedKey[nr][2];
		requiredKey[1] = expandedKey[nr][3];
		requiredKey[2] = expandedKey[nr][6];
		requiredKey[3] = expandedKey[nr][7];
		invRoundFunction8(state,requiredKey,nr);
	}
}

void test_vectors8()
{
	uint8_t state11[16] = {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};
	uint8_t state12[16] = {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};
	uint16_t masterKey11[8] = {0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000};
	uint16_t masterKey12[8] = {0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000};
	gift_encrypt8(state11,masterKey11);
	uint8_t expected1[16] = {0xcd,0x0b,0xd7,0x38,0x38,0x8a,0xd3,0xf6,0x68,0xb1,0x5a,0x36,0xce,0xb6,0xff,0x92};
	for (int i = 0; i < 16; i++)
	{
		if (state11[i] != expected1[i])
		{
			cout << "GIFT128 - Test vector 1 (enc) failed" << endl;
			exit(0);
		}
	}
	gift_decrypt8(state11,masterKey12);
	for (int i = 0; i < 16; i++)
	{
		if (state11[i] != state12[i])
		{
			cout << "GIFT128 - Test vector 1 (dec) failed" << endl;
			exit(0);
		}
	}

	uint8_t state21[16] = {0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10}; 
	uint8_t state22[16] = {0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10,0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10}; 
	uint16_t masterKey21[8] = {0xfedc,0xba98,0x7654,0x3210,0xfedc,0xba98,0x7654,0x3210}; 
	uint16_t masterKey22[8] = {0xfedc,0xba98,0x7654,0x3210,0xfedc,0xba98,0x7654,0x3210}; 
	uint8_t expected2[16] = {0x84,0x22,0x24,0x1a,0x6d,0xbf,0x5a,0x93,0x46,0xaf,0x46,0x84,0x09,0xee,0x01,0x52};
	gift_encrypt8(state21,masterKey21);
	for (int i = 0; i < 16; i++)
	{
		if (state21[i] != expected2[i])
		{
			cout << "GIFT128 - Test vector 2 (enc) failed" << endl;
			exit(0);
		}
	}
	gift_decrypt8(state21,masterKey22);
	for (int i = 0; i < 16; i++)
	{
		if (state21[i] != state22[i])
		{
			cout << "GIFT128 - Test vector 2 (dec) failed" << endl;
			exit(0);
		}
	}
	uint8_t state31[16] = {0xe3,0x9c,0x14,0x1f,0xa5,0x7d,0xba,0x43,0xf0,0x8a,0x85,0xb6,0xa9,0x1f,0x86,0xc1};
	uint8_t state32[16] = {0xe3,0x9c,0x14,0x1f,0xa5,0x7d,0xba,0x43,0xf0,0x8a,0x85,0xb6,0xa9,0x1f,0x86,0xc1};
	uint16_t masterKey31[8] = {0xd0f5,0xc59a,0x7700,0xd3e7,0x9902,0x8fa9,0xf90a,0xd837};
	uint16_t masterKey32[8] = {0xd0f5,0xc59a,0x7700,0xd3e7,0x9902,0x8fa9,0xf90a,0xd837};
	uint8_t expected3[16] = {0x13,0xed,0xe6,0x7c,0xbd,0xcc,0x3d,0xbf,0x40,0x0a,0x62,0xd6,0x97,0x72,0x65,0xea};
	gift_encrypt8(state31,masterKey31);
	for (int i = 0; i < 16; i++)
	{
		if (state31[i] != expected3[i])
		{
			cout << "GIFT128 - Test vector 3 (enc) failed" << endl;
			exit(0);
		}
	}
	gift_decrypt8(state31,masterKey32);
	for (int i = 0; i < 16; i++)
	{
		if (state31[i] != state32[i])
		{
			cout << "GIFT128 - Test vector 3 (dec) failed" << endl;
			exit(0);
		}
	}
}
void computeDDT(uint32_t DDT[16][16])
{
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++) DDT[i][j] = 0;
	}
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++) DDT[i^j][Sbox[i]^Sbox[j]]++;
	}
}

void computeXDDT(uint8_t x, uint8_t y, uint32_t xddt[16])
{
	for (int i = 0; i < 16; i++) xddt[i] = 0;
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			if (((i^j) == x) && ((Sbox[i]^Sbox[j]) == y)) xddt[i]++;
		}
	}
}

void computeYDDT(uint8_t x, uint8_t y, uint32_t yddt[16])
{
	for (int i = 0; i < 16; i++) yddt[i] = 0;
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			if (((i^j) == x) && ((Sbox[i]^Sbox[j]) == y)) yddt[Sbox[i]]++;
		}
	}
}


// int main()
// {
// 	mt19937 rand_generator;
// 	uniform_int_distribution<uint16_t> rng8(0, 0xff);
// 	rand_generator.seed(time(0));
// 	// uint8_t before[16] = {0x00,0x00,0x00,0x00,0x70,0x60,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};
// 	// uint8_t middle[16] = {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0xa0,0x00,0x00};
// 	// uint8_t after[16] = {0x00,0x00,0x00,0x10,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00};
// 	uint8_t before[16] = {0x00,0x10,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x80,0x20,0x00,0x00};
// 	uint8_t middle[16] = {0x00,0x00,0x00,0x00,0x80,0x00,0x00,0x20,0x00,0x00,0x00,0x50,0x00,0x00,0x00,0x20};
// 	uint8_t after[16] = {0x00,0x00,0x01,0x00,0x00,0x20,0x08,0x00,0x00,0x14,0x04,0x04,0x00,0x02,0x02,0x02};
// 	uint8_t after_tmp[16];
// 	for (int i = 0; i < 16; i++) after_tmp[i] = after[i];
// 	invPermutation8(after_tmp);
// 	for (int i = 0; i < 16; i++) cout << hex << int(after_tmp[i]) << " ";
// 	cout << endl;
// // exit(0);
// 	int count = 0;
// 	int k_log2 = 1;
// 	for (int k = 0; k < (1<<24); k++)
// 	{
// 		if (log2(k) > k_log2)
// 		{
// 			cout << dec << k_log2 << endl;
// 			k_log2++;
// 		}
// 		uint8_t state1[16], state2[16];
// 		for (int i = 0; i < 16; i++)
// 		{
// 			uint8_t rand_num = rng8(rand_generator);
// 			state1[i] = rand_num;
// 			state2[i] = rand_num ^ before[i];
// 		}
// 		uint16_t key1[4];
// 		uint16_t key2[4];
// 		for (int i = 0; i < 4; i++)
// 		{
// 			uint8_t rand_num = rng8(rand_generator);
// 			uint8_t rand_num2 = rng8(rand_generator);
// 			key1[i] = (uint16_t(rand_num) << 8) ^ uint16_t(rand_num2);
// 		}
// 		for (int i = 0; i < 4; i++)
// 		{
// 			uint8_t rand_num = rng8(rand_generator);
// 			uint8_t rand_num2 = rng8(rand_generator);
// 			key2[i] = (uint16_t(rand_num) << 8) ^ uint16_t(rand_num2);
// 		}
// 		Substitution8(state1);
// 		Substitution8(state2);
// 		Permutation8(state1);
// 		Permutation8(state2);
// 		AddRoundKey8(state1,key1);
// 		AddRoundKey8(state2,key1);
// 		bool t1 = true;
// 		for (int i = 0; i < 16; i++)
// 		{
// 			if ((state1[i]^state2[i]) != middle[i]) t1 = false;
// 			// cout << (state1[i] ^ state2[i]) << " ";
// 		}
// 		if (!t1) continue;
// 		Substitution8(state1);
// 		Substitution8(state2);
// 		Permutation8(state1);
// 		Permutation8(state2);
// 		AddRoundKey8(state1,key2);
// 		AddRoundKey8(state2,key2);
// 		bool t = true;
// 		for (int i = 0; i < 16; i++)
// 		{
// 			if ((state1[i]^state2[i]) != after[i]) t = false;
// 			// cout << (state1[i] ^ state2[i]) << " ";
// 		}
// 		// cout << endl;
// 		if (t) count++;
// 		if (t) cout << "count: " << count/(k+1.0) << endl;
// 	}
	



// }