#include <iostream>
#include <cstring>
#include "skinny.h"

using namespace std;
static int perm_schedule[4][4] = {{9,15,8,13},{10,14,12,11},{0,1,2,3},{4,5,6,7}};
static int inv_perm_schedule[4][4] = {{8,9,10,11},{12,13,14,15},{2,0,4,7},{6,3,5,1}};
static uint8_t Sbox[16] = {0xc,0x6,0x9,0x0,0x1,0xa,0x2,0xb,0x3,0x8,0x5,0xd,0x4,0xe,0x7,0xf};
static uint8_t invSbox[16] = {0x3,0x4,0x6,0x8,0xc,0xa,0x1,0xe,0x9,0x2,0x5,0x7,0x0,0xb,0xd,0xf};

static uint8_t Constants[62] = {0x01,0x03,0x07,0x0f,0x1f,0x3e,0x3d,0x3b,0x37,0x2f,0x1e,0x3c,0x39,0x33,0x27,0x0e,
						 0x1d,0x3a,0x35,0x2b,0x16,0x2c,0x18,0x30,0x21,0x02,0x05,0x0b,0x17,0x2e,0x1c,0x38,
						 0x31,0x23,0x06,0x0d,0x1b,0x36,0x2d,0x1a,0x34,0x29,0x12,0x24,0x08,0x11,0x22,0x04,
						 0x09,0x13,0x26,0x0c,0x19,0x32,0x25,0x0a,0x15,0x2a,0x14,0x28,0x10,0x20};

// #define printout(x) cout << #x << ": " << x << endl;
// #define printstate(s) for (int r=0;r<4;r++) {for(int c=0;c<4;c++) {cout<<hex<<int(s[r][c]);} cout<<" ";}	
// #define printstatediff(s1,s2) for (int r=0;r<4;r++) {for(int c=0;c<4;c++) {cout<<hex<<int(s1[r][c]^s2[r][c]);} cout<<" ";} cout << endl;	
// #define printkey(k) for(int tk=0;tk<TK_NUM;tk++) {for(int r=0;r<4;r++) {for(int c=0;c<4;c++) {cout<<hex<<int(k[tk][r][c]);}cout<<" ";}cout<<endl;}	
// #define printkeydiff(k1,k2) for(int tk=0;tk<TK_NUM;tk++) {for(int r=0;r<4;r++) {for(int c=0;c<4;c++) {cout<<hex<<int(k1[tk][r][c]^k2[tk][r][c]);}cout<<" ";}cout<<endl;}	

uint8_t getSbox(uint8_t val) {return Sbox[val];}
uint8_t getInvSbox(uint8_t val) {return invSbox[val];}
uint8_t getConstants(uint8_t val) {return Constants[val];}
uint8_t getPermSchedule(uint8_t val) {return perm_schedule[val/4][val%4];}

void computeDDT(int DDT[16][16])
{
	for (int i = 0; i < 16; i++) { for (int j = 0; j < 16; j++){ DDT[i][j] = 0;}}
	for (int i = 0 ; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			DDT[i^j][Sbox[i]^Sbox[j]]++;
		}
	}
}

uint8_t LFSR(uint8_t n, int v)
{
	uint8_t x3 = (n >> 3) % 2;
	uint8_t x2 = (n >> 2) % 2;
	uint8_t x1 = (n >> 1) % 2;
	uint8_t x0 = (n >> 0) % 2;
	if (v == 2) return (x2 << 3) ^ (x1 << 2) ^ (x0 << 1) ^ ((x2 ^ x3) << 0);
	if (v == 3) return ((x0 ^ x3) << 3) ^ (x3 << 2) ^ (x2 << 1) ^ (x1 << 0);
	if (v == 4) return (x2 << 3) ^ (x1 << 2) ^ ((x0 ^ x2) << 1) ^ ((x1 ^ x2 ^ x3) << 0);
}

uint8_t invLFSR(uint8_t n, int v)
{
	uint8_t x3 = (n >> 3) % 2;
	uint8_t x2 = (n >> 2) % 2;
	uint8_t x1 = (n >> 1) % 2;
	uint8_t x0 = (n >> 0) % 2;
	if (v == 2) return ((x0 ^ x3) << 3) ^ (x3 << 2) ^ (x2 << 1) ^ (x1 << 0);
	if (v == 3) return (x2 << 3) ^ (x1 << 2) ^ (x0 << 1) ^ ((x2 ^ x3) << 0);
	if (v == 4) return ((x0 ^ x2 ^ x3) << 3) ^ (x3 << 2) ^ (x2 << 1) ^ ((x1 ^ x3) << 0);
}

void perm(uint8_t key_state[4][4])
{
	uint8_t tmp[4][4];
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			tmp[r][c] = key_state[perm_schedule[r][c]/4][perm_schedule[r][c]%4];
		}
	}
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) key_state[r][c] = tmp[r][c];
	}

}

void invPerm(uint8_t key_state[4][4])
{
	uint8_t tmp[4][4];
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) tmp[r][c] = key_state[inv_perm_schedule[r][c]/4][inv_perm_schedule[r][c]%4];
	}
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) key_state[r][c] = tmp[r][c];
	}
}

void key_schedule_round_function(uint8_t key_states[4][4][4])
{
	for (int i = 0; i < 4; i++) perm(key_states[i]);
	for (int i = 1; i < 4; i++)
	{
		for (int r = 0; r < 2; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				key_states[i][r][c] = LFSR(key_states[i][r][c],i+1);
			}
		}
	}
}

void inv_key_schedule_round_function(uint8_t key_states[4][4][4])
{
	for (int i = 1; i < 4; i++)
	{
		for (int r = 0; r < 2; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				key_states[i][r][c] = invLFSR(key_states[i][r][c],i+1);
			}
		}
	}
	for (int i = 0; i < 4; i++) invPerm(key_states[i]);
}

// change this
void SR(uint8_t state[4][4])
{
	uint8_t tmp[4][4];
	tmp[0][0] = state[0][0]; tmp[0][1] = state[0][1]; tmp[0][2] = state[0][2]; tmp[0][3] = state[0][3];
	tmp[1][0] = state[1][3]; tmp[1][1] = state[1][0]; tmp[1][2] = state[1][1]; tmp[1][3] = state[1][2];
	tmp[2][0] = state[2][2]; tmp[2][1] = state[2][3]; tmp[2][2] = state[2][0]; tmp[2][3] = state[2][1];
	tmp[3][0] = state[3][1]; tmp[3][1] = state[3][2]; tmp[3][2] = state[3][3]; tmp[3][3] = state[3][0];
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			state[r][c] = tmp[r][c];
		}
	}
}

void invSR(uint8_t state[4][4])
{
	uint8_t tmp[4][4];
	tmp[0][0] = state[0][0]; tmp[0][1] = state[0][1]; tmp[0][2] = state[0][2]; tmp[0][3] = state[0][3];
	tmp[1][0] = state[1][1]; tmp[1][1] = state[1][2]; tmp[1][2] = state[1][3]; tmp[1][3] = state[1][0];
	tmp[2][0] = state[2][2]; tmp[2][1] = state[2][3]; tmp[2][2] = state[2][0]; tmp[2][3] = state[2][1];
	tmp[3][0] = state[3][3]; tmp[3][1] = state[3][0]; tmp[3][2] = state[3][1]; tmp[3][3] = state[3][2];
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			state[r][c] = tmp[r][c];
		}
	}
}

void MC(uint8_t state[4][4])
{
	uint8_t tmp[4][4];
	for (int c = 0; c < 4; c++)
	{
		tmp[1][c] = state[0][c];
		tmp[2][c] = state[1][c] ^ state[2][c];
		tmp[3][c] = state[0][c] ^ state[2][c];
		tmp[0][c] = tmp[3][c] ^ state[3][c];
	}
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			state[r][c] = tmp[r][c];
		}
	}
}

void invMC(uint8_t state[4][4])
{
	uint8_t tmp[4][4];
	for (int c = 0; c < 4; c++)
	{
		tmp[0][c] = state[1][c];
		tmp[3][c] = state[0][c] ^ state[3][c];
		tmp[2][c] = tmp[0][c] ^ state[3][c];
		tmp[1][c] = tmp[2][c] ^ state[2][c];
	}
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			state[r][c] = tmp[r][c];
		}
	}
}

void Substitution(uint8_t state[4][4])
{
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			state[r][c] = Sbox[state[r][c]];
		}
	}
}

void invSubstitution(uint8_t state[4][4])
{
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			state[r][c] = invSbox[state[r][c]];
		}
	}
}

void AddRoundTweakey(uint8_t state[4][4],uint8_t TK[4][4][4])
{
	for (int r = 0; r < 2; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			for (int k = 0; k < 4; k++)
			{
				state[r][c] ^= TK[k][r][c];
			}
		}
	}
}

void AddConstant(uint8_t state[4][4], int r)
{
	uint8_t cons = Constants[r];
	uint8_t c0 = cons & 0xf;
	uint8_t c1 = (cons >> 4) & 0x3;
	uint8_t c2 = 0x2;
	state[0][0] ^= c0;
	state[1][0] ^= c1;
	state[2][0] ^= c2;
}

void RoundFunction(uint8_t state[4][4], uint8_t TK[4][4][4], int r)
{
	Substitution(state);
 	AddConstant(state,r);
	AddRoundTweakey(state,TK);
	SR(state);
	MC(state);
}

void InvRoundFunction(uint8_t state[4][4], uint8_t TK[4][4][4], int r)
{
	invMC(state);
	invSR(state);
	AddRoundTweakey(state,TK);
	AddConstant(state,r);
	invSubstitution(state);
}

void Skinny_Enc(uint8_t state[4][4], uint8_t TK[4][4][4], int nr)
{
	for (int r = 0; r < nr; r++)
	{
		RoundFunction(state,TK,r);
		key_schedule_round_function(TK);
	}
}

void Skinny_Dec(uint8_t state[4][4], uint8_t TK[4][4][4], int nr)
{
	uint8_t TK_ALL[44][4][4][4];
	for (int rd = 0; rd < nr; rd++)
	{
		for (int k = 0; k < 4; k++)
		{
			for (int r = 0; r < 4; r++)
			{
				for (int c = 0; c < 4; c++)
				{
					TK_ALL[rd][k][r][c] = TK[k][r][c];
				}
			}
		}
		key_schedule_round_function(TK);
	}

	for (int r = nr-1; r >= 0; r--) InvRoundFunction(state,TK_ALL[r],r);
}

void reset_key(uint8_t key_tmp[4][4][4], uint8_t key[4][4][4])
{
	for (int k = 0; k < 4; k++)
	{
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++) key_tmp[k][r][c] = key[k][r][c];
		}
	}
}

void reset_plaintext(uint8_t p_tmp[4][4], uint8_t p[4][4])
{
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++) 
		{
			p_tmp[r][c] = p[r][c];
		}
	}
}

void test_Vectors()
{
	int TK_NUM;
	uint8_t TK[4][4][4] = {0};
	cout << "Testing for encryption..." << endl;
	cout << "test vector 1: ";
	int rounds = 32;
	TK_NUM = 1;
	TK[0][0][0] = 0xf; TK[0][0][1] = 0x5; TK[0][0][2] = 0x2; TK[0][0][3] = 0x6; 
	TK[0][1][0] = 0x9; TK[0][1][1] = 0x8; TK[0][1][2] = 0x2; TK[0][1][3] = 0x6; 
	TK[0][2][0] = 0xf; TK[0][2][1] = 0xc; TK[0][2][2] = 0x6; TK[0][2][3] = 0x8; 
	TK[0][3][0] = 0x1; TK[0][3][1] = 0x2; TK[0][3][2] = 0x3; TK[0][3][3] = 0x8;  
	
	uint8_t state1[4][4] = {{0x0,0x6,0x0,0x3},{0x4,0xf,0x9,0x5},{0x7,0x7,0x2,0x4},{0xd,0x1,0x9,0xd}};
	Skinny_Enc(state1,TK,rounds);
	uint8_t expected1[4][4] = {{0xb,0xb,0x3,0x9},{0xd,0xf,0xb,0x2},{0x4,0x2,0x9,0xb},{0x8,0xa,0xc,0x7}};
	bool flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state1[r][c] != expected1[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;
	cout << "test vector 2: ";
	rounds = 36;
	TK_NUM = 2;

	TK[0][0][0] = 0x9; TK[0][0][1] = 0xe; TK[0][0][2] = 0xb; TK[0][0][3] = 0x9; 
	TK[0][1][0] = 0x3; TK[0][1][1] = 0x6; TK[0][1][2] = 0x4; TK[0][1][3] = 0x0; 
	TK[0][2][0] = 0xd; TK[0][2][1] = 0x0; TK[0][2][2] = 0x8; TK[0][2][3] = 0x8; 
	TK[0][3][0] = 0xd; TK[0][3][1] = 0xa; TK[0][3][2] = 0x6; TK[0][3][3] = 0x3; 

	TK[1][0][0] = 0x7; TK[1][0][1] = 0x6; TK[1][0][2] = 0xa; TK[1][0][3] = 0x3; 
	TK[1][1][0] = 0x9; TK[1][1][1] = 0xd; TK[1][1][2] = 0x1; TK[1][1][3] = 0xc; 
	TK[1][2][0] = 0x8; TK[1][2][1] = 0xb; TK[1][2][2] = 0xe; TK[1][2][3] = 0xa; 
	TK[1][3][0] = 0x7; TK[1][3][1] = 0x1; TK[1][3][2] = 0xe; TK[1][3][3] = 0x1; 

	uint8_t state2[4][4] = {{0xc,0xf,0x1,0x6},{0xc,0xf,0xe,0x8},{0xf,0xd,0x0,0xf},{0x9,0x8,0xa,0xa}};
	Skinny_Enc(state2,TK,rounds);
	uint8_t expected2[4][4] = {{0x6,0xc,0xe,0xd},{0xa,0x1,0xf,0x4},{0x3,0xd,0xe,0x9},{0x2,0xb,0x9,0xe}};
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state2[r][c] != expected2[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

	cout << "test vector 3: ";
	rounds = 40;
	TK_NUM = 3;
	TK[0][0][0] = 0xe; TK[0][0][1] = 0xd; TK[0][0][2] = 0x0; TK[0][0][3] = 0x0; 
	TK[0][1][0] = 0xc; TK[0][1][1] = 0x8; TK[0][1][2] = 0x5; TK[0][1][3] = 0xb; 
	TK[0][2][0] = 0x1; TK[0][2][1] = 0x2; TK[0][2][2] = 0x0; TK[0][2][3] = 0xd; 
	TK[0][3][0] = 0x6; TK[0][3][1] = 0x8; TK[0][3][2] = 0x6; TK[0][3][3] = 0x1;

	TK[1][0][0] = 0x8; TK[1][0][1] = 0x7; TK[1][0][2] = 0x5; TK[1][0][3] = 0x3; 
	TK[1][1][0] = 0xe; TK[1][1][1] = 0x2; TK[1][1][2] = 0x4; TK[1][1][3] = 0xb; 
	TK[1][2][0] = 0xf; TK[1][2][1] = 0xd; TK[1][2][2] = 0x9; TK[1][2][3] = 0x0; 
	TK[1][3][0] = 0x8; TK[1][3][1] = 0xf; TK[1][3][2] = 0x6; TK[1][3][3] = 0x0;

	TK[2][0][0] = 0xb; TK[2][0][1] = 0x2; TK[2][0][2] = 0xd; TK[2][0][3] = 0xb; 
	TK[2][1][0] = 0xb; TK[2][1][1] = 0x4; TK[2][1][2] = 0x1; TK[2][1][3] = 0xb; 
	TK[2][2][0] = 0x4; TK[2][2][1] = 0x2; TK[2][2][2] = 0x2; TK[2][2][3] = 0xd; 
	TK[2][3][0] = 0xf; TK[2][3][1] = 0xc; TK[2][3][2] = 0xd; TK[2][3][3] = 0x0;

	uint8_t state3[4][4] = {{0x5,0x3,0x0,0xc},{0x6,0x1,0xd,0x3},{0x5,0xe,0x8,0x6},{0x6,0x3,0xc,0x3}};
	uint8_t expected3[4][4] = {{0xd,0xd,0x2,0xc},{0xf,0x1,0xa,0x8},{0xf,0x3,0x3,0x0},{0x3,0x0,0x3,0xc}};
	Skinny_Enc(state3,TK,rounds);
	
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state3[r][c] != expected3[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

	cout << "Testing for decryption..." << endl;
	cout << "test vector 4: ";
	rounds = 32;
	TK_NUM = 1;
	TK[0][0][0] = 0xf; TK[0][0][1] = 0x5; TK[0][0][2] = 0x2; TK[0][0][3] = 0x6; 
	TK[0][1][0] = 0x9; TK[0][1][1] = 0x8; TK[0][1][2] = 0x2; TK[0][1][3] = 0x6; 
	TK[0][2][0] = 0xf; TK[0][2][1] = 0xc; TK[0][2][2] = 0x6; TK[0][2][3] = 0x8; 
	TK[0][3][0] = 0x1; TK[0][3][1] = 0x2; TK[0][3][2] = 0x3; TK[0][3][3] = 0x8;  
	uint8_t state4[4][4] = {{0xb,0xb,0x3,0x9},{0xd,0xf,0xb,0x2},{0x4,0x2,0x9,0xb},{0x8,0xa,0xc,0x7}};
	uint8_t expected4[4][4] = {{0x0,0x6,0x0,0x3},{0x4,0xf,0x9,0x5},{0x7,0x7,0x2,0x4},{0xd,0x1,0x9,0xd}};
	Skinny_Dec(state4,TK,rounds);
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state4[r][c] != expected4[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

	cout << "test vector 5: ";
	rounds = 36;
	TK_NUM = 2;

	TK[0][0][0] = 0x9; TK[0][0][1] = 0xe; TK[0][0][2] = 0xb; TK[0][0][3] = 0x9; 
	TK[0][1][0] = 0x3; TK[0][1][1] = 0x6; TK[0][1][2] = 0x4; TK[0][1][3] = 0x0; 
	TK[0][2][0] = 0xd; TK[0][2][1] = 0x0; TK[0][2][2] = 0x8; TK[0][2][3] = 0x8; 
	TK[0][3][0] = 0xd; TK[0][3][1] = 0xa; TK[0][3][2] = 0x6; TK[0][3][3] = 0x3; 

	TK[1][0][0] = 0x7; TK[1][0][1] = 0x6; TK[1][0][2] = 0xa; TK[1][0][3] = 0x3; 
	TK[1][1][0] = 0x9; TK[1][1][1] = 0xd; TK[1][1][2] = 0x1; TK[1][1][3] = 0xc; 
	TK[1][2][0] = 0x8; TK[1][2][1] = 0xb; TK[1][2][2] = 0xe; TK[1][2][3] = 0xa; 
	TK[1][3][0] = 0x7; TK[1][3][1] = 0x1; TK[1][3][2] = 0xe; TK[1][3][3] = 0x1; 

	uint8_t state5[4][4] = {{0x6,0xc,0xe,0xd},{0xa,0x1,0xf,0x4},{0x3,0xd,0xe,0x9},{0x2,0xb,0x9,0xe}};
	uint8_t expected5[4][4] = {{0xc,0xf,0x1,0x6},{0xc,0xf,0xe,0x8},{0xf,0xd,0x0,0xf},{0x9,0x8,0xa,0xa}};
	Skinny_Dec(state5,TK,rounds);
	
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state5[r][c] != expected5[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

	cout << "test vector 6: ";
	rounds = 40;
	TK_NUM = 3;
	TK[0][0][0] = 0xe; TK[0][0][1] = 0xd; TK[0][0][2] = 0x0; TK[0][0][3] = 0x0; 
	TK[0][1][0] = 0xc; TK[0][1][1] = 0x8; TK[0][1][2] = 0x5; TK[0][1][3] = 0xb; 
	TK[0][2][0] = 0x1; TK[0][2][1] = 0x2; TK[0][2][2] = 0x0; TK[0][2][3] = 0xd; 
	TK[0][3][0] = 0x6; TK[0][3][1] = 0x8; TK[0][3][2] = 0x6; TK[0][3][3] = 0x1;

	TK[1][0][0] = 0x8; TK[1][0][1] = 0x7; TK[1][0][2] = 0x5; TK[1][0][3] = 0x3; 
	TK[1][1][0] = 0xe; TK[1][1][1] = 0x2; TK[1][1][2] = 0x4; TK[1][1][3] = 0xb; 
	TK[1][2][0] = 0xf; TK[1][2][1] = 0xd; TK[1][2][2] = 0x9; TK[1][2][3] = 0x0; 
	TK[1][3][0] = 0x8; TK[1][3][1] = 0xf; TK[1][3][2] = 0x6; TK[1][3][3] = 0x0;

	TK[2][0][0] = 0xb; TK[2][0][1] = 0x2; TK[2][0][2] = 0xd; TK[2][0][3] = 0xb; 
	TK[2][1][0] = 0xb; TK[2][1][1] = 0x4; TK[2][1][2] = 0x1; TK[2][1][3] = 0xb; 
	TK[2][2][0] = 0x4; TK[2][2][1] = 0x2; TK[2][2][2] = 0x2; TK[2][2][3] = 0xd; 
	TK[2][3][0] = 0xf; TK[2][3][1] = 0xc; TK[2][3][2] = 0xd; TK[2][3][3] = 0x0;

	uint8_t state6[4][4] = {{0xd,0xd,0x2,0xc},{0xf,0x1,0xa,0x8},{0xf,0x3,0x3,0x0},{0x3,0x0,0x3,0xc}};
	uint8_t expected6[4][4] = {{0x5,0x3,0x0,0xc},{0x6,0x1,0xd,0x3},{0x5,0xe,0x8,0x6},{0x6,0x3,0xc,0x3}};
	Skinny_Dec(state6,TK,rounds);
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state6[r][c] != expected6[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

}



// SKINNY8
static uint8_t Sbox8[256] = {0x65, 0x4c, 0x6a, 0x42, 0x4b, 0x63, 0x43, 0x6b, 0x55, 0x75, 0x5a, 0x7a, 0x53, 0x73, 0x5b, 0x7b,
					    	0x35, 0x8c, 0x3a, 0x81, 0x89, 0x33, 0x80, 0x3b, 0x95, 0x25, 0x98, 0x2a, 0x90, 0x23, 0x99, 0x2b,
					    	0xe5, 0xcc, 0xe8, 0xc1, 0xc9, 0xe0, 0xc0, 0xe9, 0xd5, 0xf5, 0xd8, 0xf8, 0xd0, 0xf0, 0xd9, 0xf9,
					    	0xa5, 0x1c, 0xa8, 0x12, 0x1b, 0xa0, 0x13, 0xa9, 0x05, 0xb5, 0x0a, 0xb8, 0x03, 0xb0, 0x0b, 0xb9,
					    	0x32, 0x88, 0x3c, 0x85, 0x8d, 0x34, 0x84, 0x3d, 0x91, 0x22, 0x9c, 0x2c, 0x94, 0x24, 0x9d, 0x2d,
					    	0x62, 0x4a, 0x6c, 0x45, 0x4d, 0x64, 0x44, 0x6d, 0x52, 0x72, 0x5c, 0x7c, 0x54, 0x74, 0x5d, 0x7d,
					    	0xa1, 0x1a, 0xac, 0x15, 0x1d, 0xa4, 0x14, 0xad, 0x02, 0xb1, 0x0c, 0xbc, 0x04, 0xb4, 0x0d, 0xbd,
					    	0xe1, 0xc8, 0xec, 0xc5, 0xcd, 0xe4, 0xc4, 0xed, 0xd1, 0xf1, 0xdc, 0xfc, 0xd4, 0xf4, 0xdd, 0xfd,
					    	0x36, 0x8e, 0x38, 0x82, 0x8b, 0x30, 0x83, 0x39, 0x96, 0x26, 0x9a, 0x28, 0x93, 0x20, 0x9b, 0x29,
					    	0x66, 0x4e, 0x68, 0x41, 0x49, 0x60, 0x40, 0x69, 0x56, 0x76, 0x58, 0x78, 0x50, 0x70, 0x59, 0x79,
					    	0xa6, 0x1e, 0xaa, 0x11, 0x19, 0xa3, 0x10, 0xab, 0x06, 0xb6, 0x08, 0xba, 0x00, 0xb3, 0x09, 0xbb,
					    	0xe6, 0xce, 0xea, 0xc2, 0xcb, 0xe3, 0xc3, 0xeb, 0xd6, 0xf6, 0xda, 0xfa, 0xd3, 0xf3, 0xdb, 0xfb,
					    	0x31, 0x8a, 0x3e, 0x86, 0x8f, 0x37, 0x87, 0x3f, 0x92, 0x21, 0x9e, 0x2e, 0x97, 0x27, 0x9f, 0x2f,
					    	0x61, 0x48, 0x6e, 0x46, 0x4f, 0x67, 0x47, 0x6f, 0x51, 0x71, 0x5e, 0x7e, 0x57, 0x77, 0x5f, 0x7f,
					    	0xa2, 0x18, 0xae, 0x16, 0x1f, 0xa7, 0x17, 0xaf, 0x01, 0xb2, 0x0e, 0xbe, 0x07, 0xb7, 0x0f, 0xbf,
					    	0xe2, 0xca, 0xee, 0xc6, 0xcf, 0xe7, 0xc7, 0xef, 0xd2, 0xf2, 0xde, 0xfe, 0xd7, 0xf7, 0xdf, 0xff};
static uint8_t invSbox8[256] = {0xac, 0xe8, 0x68, 0x3c, 0x6c, 0x38, 0xa8, 0xec, 0xaa, 0xae, 0x3a, 0x3e, 0x6a, 0x6e, 0xea, 0xee, 
							   0xa6, 0xa3, 0x33, 0x36, 0x66, 0x63, 0xe3, 0xe6, 0xe1, 0xa4, 0x61, 0x34, 0x31, 0x64, 0xa1, 0xe4, 
							   0x8d, 0xc9, 0x49, 0x1d, 0x4d, 0x19, 0x89, 0xcd, 0x8b, 0x8f, 0x1b, 0x1f, 0x4b, 0x4f, 0xcb, 0xcf, 
							   0x85, 0xc0, 0x40, 0x15, 0x45, 0x10, 0x80, 0xc5, 0x82, 0x87, 0x12, 0x17, 0x42, 0x47, 0xc2, 0xc7, 
							   0x96, 0x93, 0x03, 0x06, 0x56, 0x53, 0xd3, 0xd6, 0xd1, 0x94, 0x51, 0x04, 0x01, 0x54, 0x91, 0xd4, 
							   0x9c, 0xd8, 0x58, 0x0c, 0x5c, 0x08, 0x98, 0xdc, 0x9a, 0x9e, 0x0a, 0x0e, 0x5a, 0x5e, 0xda, 0xde, 
							   0x95, 0xd0, 0x50, 0x05, 0x55, 0x00, 0x90, 0xd5, 0x92, 0x97, 0x02, 0x07, 0x52, 0x57, 0xd2, 0xd7, 
							   0x9d, 0xd9, 0x59, 0x0d, 0x5d, 0x09, 0x99, 0xdd, 0x9b, 0x9f, 0x0b, 0x0f, 0x5b, 0x5f, 0xdb, 0xdf, 
							   0x16, 0x13, 0x83, 0x86, 0x46, 0x43, 0xc3, 0xc6, 0x41, 0x14, 0xc1, 0x84, 0x11, 0x44, 0x81, 0xc4, 
							   0x1c, 0x48, 0xc8, 0x8c, 0x4c, 0x18, 0x88, 0xcc, 0x1a, 0x1e, 0x8a, 0x8e, 0x4a, 0x4e, 0xca, 0xce, 
							   0x35, 0x60, 0xe0, 0xa5, 0x65, 0x30, 0xa0, 0xe5, 0x32, 0x37, 0xa2, 0xa7, 0x62, 0x67, 0xe2, 0xe7, 
							   0x3d, 0x69, 0xe9, 0xad, 0x6d, 0x39, 0xa9, 0xed, 0x3b, 0x3f, 0xab, 0xaf, 0x6b, 0x6f, 0xeb, 0xef, 
							   0x26, 0x23, 0xb3, 0xb6, 0x76, 0x73, 0xf3, 0xf6, 0x71, 0x24, 0xf1, 0xb4, 0x21, 0x74, 0xb1, 0xf4, 
							   0x2c, 0x78, 0xf8, 0xbc, 0x7c, 0x28, 0xb8, 0xfc, 0x2a, 0x2e, 0xba, 0xbe, 0x7a, 0x7e, 0xfa, 0xfe, 
							   0x25, 0x70, 0xf0, 0xb5, 0x75, 0x20, 0xb0, 0xf5, 0x22, 0x27, 0xb2, 0xb7, 0x72, 0x77, 0xf2, 0xf7, 
							   0x2d, 0x79, 0xf9, 0xbd, 0x7d, 0x29, 0xb9, 0xfd, 0x2b, 0x2f, 0xbb, 0xbf, 0x7b, 0x7f, 0xfb, 0xff};

uint8_t getSbox8(uint8_t val) {return Sbox8[val];}
uint8_t getInvSbox8(uint8_t val) {return invSbox8[val];}

void computeDDT8(int DDT[256][256])
{
	for (int i = 0; i < 256; i++) { for (int j = 0; j < 256; j++){ DDT[i][j] = 0;}}
	for (int i = 0 ; i < 256; i++)
	{
		for (int j = 0; j < 256; j++)
		{
			DDT[i^j][Sbox8[i]^Sbox8[j]]++;
		}
	}
}

uint8_t LFSR8(uint8_t n, int v)
{
	uint8_t x7 = (n >> 7) % 2;
	uint8_t x6 = (n >> 6) % 2;
	uint8_t x5 = (n >> 5) % 2;
	uint8_t x4 = (n >> 4) % 2;
	uint8_t x3 = (n >> 3) % 2;
	uint8_t x2 = (n >> 2) % 2;
	uint8_t x1 = (n >> 1) % 2;
	uint8_t x0 = (n >> 0) % 2;
	if (v == 2) return (x6 << 7) ^ (x5 << 6) ^ (x4 << 5) ^ (x3 << 4) ^ (x2 << 3) ^ (x1 << 2) ^ (x0 << 1) ^ ((x7 ^ x5) << 0);
	else if (v == 3) return ((x0 ^ x6) << 7) ^ (x7 << 6) ^ (x6 << 5) ^ (x5 << 4) ^ (x4 << 3) ^ (x3 << 2) ^ (x2 << 1) ^ (x1 << 0);
	else if (v == 4) return 0;
}

uint8_t invLFSR8(uint8_t n, int v)
{
	uint8_t x7 = (n >> 7) % 2;
	uint8_t x6 = (n >> 6) % 2;
	uint8_t x5 = (n >> 5) % 2;
	uint8_t x4 = (n >> 4) % 2;
	uint8_t x3 = (n >> 3) % 2;
	uint8_t x2 = (n >> 2) % 2;
	uint8_t x1 = (n >> 1) % 2;
	uint8_t x0 = (n >> 0) % 2;
	if (v == 2) return ((x6 ^ x0) << 7) ^ (x7 << 6) ^ (x6 << 5) ^ (x5 << 4) ^ (x4 << 3) ^ (x3 << 2) ^ (x2 << 1) ^ (x1 << 0);
	else if (v == 3) return (x6 << 7) ^ (x5 << 6) ^ (x4 << 5) ^ (x3 << 4) ^ (x2 << 3)^ (x1 << 2)^ (x0 << 1) ^ (x7 ^ x5);
	else if (v == 4) return 0;
}

void key_schedule_round_function8(uint8_t key_states[4][4][4])
{
	for (int i = 0; i < 4; i++) perm(key_states[i]);
	for (int i = 1; i < 4; i++)
	{
		for (int r = 0; r < 2; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				key_states[i][r][c] = LFSR8(key_states[i][r][c],i+1);
			}
		}
	}
}

void inv_key_schedule_round_function8(uint8_t key_states[4][4][4])
{
	for (int i = 1; i < 4; i++)
	{
		for (int r = 0; r < 2; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				key_states[i][r][c] = invLFSR8(key_states[i][r][c],i+1);
			}
		}
	}
	for (int i = 0; i < 4; i++) invPerm(key_states[i]);
}


void Substitution8(uint8_t state[4][4])
{
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			state[r][c] = Sbox8[state[r][c]];
		}
	}
}

void invSubstitution8(uint8_t state[4][4])
{
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			state[r][c] = invSbox8[state[r][c]];
		}
	}
}
void RoundFunction8(uint8_t state[4][4], uint8_t TK[4][4][4], int r)
{
	Substitution8(state);
 	AddConstant(state,r);
	AddRoundTweakey(state,TK);
	SR(state);
	MC(state);

}
void InvRoundFunction8(uint8_t state[4][4], uint8_t TK[4][4][4], int r)
{
	invMC(state);
	invSR(state);
	AddRoundTweakey(state,TK);
	AddConstant(state,r);
	invSubstitution8(state);
}
void Skinny_Enc8(uint8_t state[4][4], uint8_t TK[4][4][4], int nr)
{
	for (int r = 0; r < nr; r++)
	{
		RoundFunction8(state,TK,r);
		key_schedule_round_function8(TK);
	}
}
void Skinny_Dec8(uint8_t state[4][4], uint8_t TK[4][4][4], int nr)
{
	uint8_t TK_ALL[100][4][4][4];
	for (int rd = 0; rd < nr; rd++)
	{
		for (int k = 0; k < 4; k++)
		{
			for (int r = 0; r < 4; r++)
			{
				for (int c = 0; c < 4; c++)
				{
					TK_ALL[rd][k][r][c] = TK[k][r][c];
				}
			}
		}
		key_schedule_round_function8(TK);
	}

	for (int r = nr-1; r >= 0; r--) InvRoundFunction8(state,TK_ALL[r],r);
}


void test_Vectors8()
{
	int TK_NUM;
	uint8_t TK[4][4][4] = {0};
	cout << "Testing for encryption..." << endl;
	cout << "test vector 1: ";
	int rounds = 40;
	TK_NUM = 1;
	TK[0][0][0] = 0x4f; TK[0][0][1] = 0x55; TK[0][0][2] = 0xcf; TK[0][0][3] = 0xb0; 
	TK[0][1][0] = 0x52; TK[0][1][1] = 0x0c; TK[0][1][2] = 0xac; TK[0][1][3] = 0x52; 
	TK[0][2][0] = 0xfd; TK[0][2][1] = 0x92; TK[0][2][2] = 0xc1; TK[0][2][3] = 0x5f; 
	TK[0][3][0] = 0x37; TK[0][3][1] = 0x07; TK[0][3][2] = 0x3e; TK[0][3][3] = 0x93;  
	
	uint8_t state1[4][4] = {{0xf2,0x0a,0xdb,0x0e},{0xb0,0x8b,0x64,0x8a},{0x3b,0x2e,0xee,0xd1},{0xf0,0xad,0xda,0x14}};
	Skinny_Enc8(state1,TK,rounds);
	uint8_t expected1[4][4] = {{0x22,0xff,0x30,0xd4},{0x98,0xea,0x62,0xd7},{0xe4,0x5b,0x47,0x6e},{0x33,0x67,0x5b,0x74}};
	bool flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state1[r][c] != expected1[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;
	cout << "test vector 2: ";
	rounds = 48;
	TK_NUM = 2;
	TK[0][0][0] = 0x00; TK[0][0][1] = 0x9c; TK[0][0][2] = 0xec; TK[0][0][3] = 0x81; 
	TK[0][1][0] = 0x60; TK[0][1][1] = 0x5d; TK[0][1][2] = 0x4a; TK[0][1][3] = 0xc1; 
	TK[0][2][0] = 0xd2; TK[0][2][1] = 0xae; TK[0][2][2] = 0x9e; TK[0][2][3] = 0x30; 
	TK[0][3][0] = 0x85; TK[0][3][1] = 0xd7; TK[0][3][2] = 0xa1; TK[0][3][3] = 0xf3; 

	TK[1][0][0] = 0x1a; TK[1][0][1] = 0xc1; TK[1][0][2] = 0x23; TK[1][0][3] = 0xeb; 
	TK[1][1][0] = 0xfc; TK[1][1][1] = 0x00; TK[1][1][2] = 0xfd; TK[1][1][3] = 0xdc; 
	TK[1][2][0] = 0xf0; TK[1][2][1] = 0x10; TK[1][2][2] = 0x46; TK[1][2][3] = 0xce; 
	TK[1][3][0] = 0xed; TK[1][3][1] = 0xdf; TK[1][3][2] = 0xca; TK[1][3][3] = 0xb3; 

	uint8_t state2[4][4] = {{0x3a,0x0c,0x47,0x76},{0x7a,0x26,0xa6,0x8d},{0xd3,0x82,0xa6,0x95},{0xe7,0x02,0x2e,0x25}};
	Skinny_Enc8(state2,TK,rounds);
	uint8_t expected2[4][4] = {{0xb7,0x31,0xd9,0x8a},{0x4b,0xde,0x14,0x7a},{0x7e,0xd4,0xa6,0xf1},{0x6b,0x9b,0x58,0x7f}};
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state2[r][c] != expected2[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

	cout << "test vector 3: ";
	rounds = 56;
	TK_NUM = 3;
	
	TK[0][0][0] = 0xdf; TK[0][0][1] = 0x88; TK[0][0][2] = 0x95; TK[0][0][3] = 0x48; 
	TK[0][1][0] = 0xcf; TK[0][1][1] = 0xc7; TK[0][1][2] = 0xea; TK[0][1][3] = 0x52; 
	TK[0][2][0] = 0xd2; TK[0][2][1] = 0x96; TK[0][2][2] = 0x33; TK[0][2][3] = 0x93; 
	TK[0][3][0] = 0x01; TK[0][3][1] = 0x79; TK[0][3][2] = 0x74; TK[0][3][3] = 0x49;
	
	TK[1][0][0] = 0xab; TK[1][0][1] = 0x58; TK[1][0][2] = 0x8a; TK[1][0][3] = 0x34; 
	TK[1][1][0] = 0xa4; TK[1][1][1] = 0x7f; TK[1][1][2] = 0x1a; TK[1][1][3] = 0xb2; 
	TK[1][2][0] = 0xdf; TK[1][2][1] = 0xe9; TK[1][2][2] = 0xc8; TK[1][2][3] = 0x29; 
	TK[1][3][0] = 0x3f; TK[1][3][1] = 0xbe; TK[1][3][2] = 0xa9; TK[1][3][3] = 0xa5;
	
	TK[2][0][0] = 0xab; TK[2][0][1] = 0x1a; TK[2][0][2] = 0xfa; TK[2][0][3] = 0xc2; 
	TK[2][1][0] = 0x61; TK[2][1][1] = 0x10; TK[2][1][2] = 0x12; TK[2][1][3] = 0xcd; 
	TK[2][2][0] = 0x8c; TK[2][2][1] = 0xef; TK[2][2][2] = 0x95; TK[2][2][3] = 0x26; 
	TK[2][3][0] = 0x18; TK[2][3][1] = 0xc3; TK[2][3][2] = 0xeb; TK[2][3][3] = 0xe8;

	uint8_t state3[4][4] = {{0xa3,0x99,0x4b,0x66},{0xad,0x85,0xa3,0x45},{0x9f,0x44,0xe9,0x2b},{0x08,0xf5,0x50,0xcb}};
	uint8_t expected3[4][4] = {{0x94,0xec,0xf5,0x89},{0xe2,0x01,0x7c,0x60},{0x1b,0x38,0xc6,0x34},{0x6a,0x10,0xdc,0xfa}};
	Skinny_Enc8(state3,TK,rounds);
	
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state3[r][c] != expected3[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

	cout << "Testing for decryption..." << endl;
	cout << "test vector 4: ";
	rounds = 40;
	TK_NUM = 1;
	TK[0][0][0] = 0x4f; TK[0][0][1] = 0x55; TK[0][0][2] = 0xcf; TK[0][0][3] = 0xb0; 
	TK[0][1][0] = 0x52; TK[0][1][1] = 0x0c; TK[0][1][2] = 0xac; TK[0][1][3] = 0x52; 
	TK[0][2][0] = 0xfd; TK[0][2][1] = 0x92; TK[0][2][2] = 0xc1; TK[0][2][3] = 0x5f; 
	TK[0][3][0] = 0x37; TK[0][3][1] = 0x07; TK[0][3][2] = 0x3e; TK[0][3][3] = 0x93;

	TK[1][0][0] = 0x00; TK[1][0][1] = 0x00; TK[1][0][2] = 0x00; TK[1][0][3] = 0x00; 
	TK[1][1][0] = 0x00; TK[1][1][1] = 0x00; TK[1][1][2] = 0x00; TK[1][1][3] = 0x00; 
	TK[1][2][0] = 0x00; TK[1][2][1] = 0x00; TK[1][2][2] = 0x00; TK[1][2][3] = 0x00; 
	TK[1][3][0] = 0x00; TK[1][3][1] = 0x00; TK[1][3][2] = 0x00; TK[1][3][3] = 0x00;  

	TK[2][0][0] = 0x00; TK[2][0][1] = 0x00; TK[2][0][2] = 0x00; TK[2][0][3] = 0x00; 
	TK[2][1][0] = 0x00; TK[2][1][1] = 0x00; TK[2][1][2] = 0x00; TK[2][1][3] = 0x00; 
	TK[2][2][0] = 0x00; TK[2][2][1] = 0x00; TK[2][2][2] = 0x00; TK[2][2][3] = 0x00; 
	TK[2][3][0] = 0x00; TK[2][3][1] = 0x00; TK[2][3][2] = 0x00; TK[2][3][3] = 0x00;  


	uint8_t state4[4][4] = {{0x22,0xff,0x30,0xd4},{0x98,0xea,0x62,0xd7},{0xe4,0x5b,0x47,0x6e},{0x33,0x67,0x5b,0x74}};
	uint8_t expected4[4][4] = {{0xf2,0x0a,0xdb,0x0e},{0xb0,0x8b,0x64,0x8a},{0x3b,0x2e,0xee,0xd1},{0xf0,0xad,0xda,0x14}};
	Skinny_Dec8(state4,TK,rounds);
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state4[r][c] != expected4[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}

		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

	cout << "test vector 5: ";
	rounds = 48;
	TK_NUM = 2;

	TK[0][0][0] = 0x00; TK[0][0][1] = 0x9c; TK[0][0][2] = 0xec; TK[0][0][3] = 0x81; 
	TK[0][1][0] = 0x60; TK[0][1][1] = 0x5d; TK[0][1][2] = 0x4a; TK[0][1][3] = 0xc1; 
	TK[0][2][0] = 0xd2; TK[0][2][1] = 0xae; TK[0][2][2] = 0x9e; TK[0][2][3] = 0x30; 
	TK[0][3][0] = 0x85; TK[0][3][1] = 0xd7; TK[0][3][2] = 0xa1; TK[0][3][3] = 0xf3; 

	TK[1][0][0] = 0x1a; TK[1][0][1] = 0xc1; TK[1][0][2] = 0x23; TK[1][0][3] = 0xeb; 
	TK[1][1][0] = 0xfc; TK[1][1][1] = 0x00; TK[1][1][2] = 0xfd; TK[1][1][3] = 0xdc; 
	TK[1][2][0] = 0xf0; TK[1][2][1] = 0x10; TK[1][2][2] = 0x46; TK[1][2][3] = 0xce; 
	TK[1][3][0] = 0xed; TK[1][3][1] = 0xdf; TK[1][3][2] = 0xca; TK[1][3][3] = 0xb3; 

	TK[2][0][0] = 0x00; TK[2][0][1] = 0x00; TK[2][0][2] = 0x00; TK[2][0][3] = 0x00; 
	TK[2][1][0] = 0x00; TK[2][1][1] = 0x00; TK[2][1][2] = 0x00; TK[2][1][3] = 0x00; 
	TK[2][2][0] = 0x00; TK[2][2][1] = 0x00; TK[2][2][2] = 0x00; TK[2][2][3] = 0x00; 
	TK[2][3][0] = 0x00; TK[2][3][1] = 0x00; TK[2][3][2] = 0x00; TK[2][3][3] = 0x00;


	uint8_t state5[4][4] = {{0xb7,0x31,0xd9,0x8a},{0x4b,0xde,0x14,0x7a},{0x7e,0xd4,0xa6,0xf1},{0x6b,0x9b,0x58,0x7f}};
	uint8_t expected5[4][4] = {{0x3a,0x0c,0x47,0x76},{0x7a,0x26,0xa6,0x8d},{0xd3,0x82,0xa6,0x95},{0xe7,0x02,0x2e,0x25}};
	Skinny_Dec8(state5,TK,rounds);
	
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state5[r][c] != expected5[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

	cout << "test vector 6: ";
	rounds = 56;
	TK_NUM = 3;
	TK[0][0][0] = 0xdf; TK[0][0][1] = 0x88; TK[0][0][2] = 0x95; TK[0][0][3] = 0x48; 
	TK[0][1][0] = 0xcf; TK[0][1][1] = 0xc7; TK[0][1][2] = 0xea; TK[0][1][3] = 0x52; 
	TK[0][2][0] = 0xd2; TK[0][2][1] = 0x96; TK[0][2][2] = 0x33; TK[0][2][3] = 0x93; 
	TK[0][3][0] = 0x01; TK[0][3][1] = 0x79; TK[0][3][2] = 0x74; TK[0][3][3] = 0x49;
	
	TK[1][0][0] = 0xab; TK[1][0][1] = 0x58; TK[1][0][2] = 0x8a; TK[1][0][3] = 0x34; 
	TK[1][1][0] = 0xa4; TK[1][1][1] = 0x7f; TK[1][1][2] = 0x1a; TK[1][1][3] = 0xb2; 
	TK[1][2][0] = 0xdf; TK[1][2][1] = 0xe9; TK[1][2][2] = 0xc8; TK[1][2][3] = 0x29; 
	TK[1][3][0] = 0x3f; TK[1][3][1] = 0xbe; TK[1][3][2] = 0xa9; TK[1][3][3] = 0xa5;
	
	TK[2][0][0] = 0xab; TK[2][0][1] = 0x1a; TK[2][0][2] = 0xfa; TK[2][0][3] = 0xc2; 
	TK[2][1][0] = 0x61; TK[2][1][1] = 0x10; TK[2][1][2] = 0x12; TK[2][1][3] = 0xcd; 
	TK[2][2][0] = 0x8c; TK[2][2][1] = 0xef; TK[2][2][2] = 0x95; TK[2][2][3] = 0x26; 
	TK[2][3][0] = 0x18; TK[2][3][1] = 0xc3; TK[2][3][2] = 0xeb; TK[2][3][3] = 0xe8;

	uint8_t state6[4][4] = {{0x94,0xec,0xf5,0x89},{0xe2,0x01,0x7c,0x60},{0x1b,0x38,0xc6,0x34},{0x6a,0x10,0xdc,0xfa}};
	uint8_t expected6[4][4] = {{0xa3,0x99,0x4b,0x66},{0xad,0x85,0xa3,0x45},{0x9f,0x44,0xe9,0x2b},{0x08,0xf5,0x50,0xcb}};
	Skinny_Dec8(state6,TK,rounds);
	flag = true;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			if (state6[r][c] != expected6[r][c])
			{
				cout << "failed!" << endl; flag = false; break;
			}
		}
		if (flag == false) break;
	}
	if (flag) cout << "passed!" << endl;

}