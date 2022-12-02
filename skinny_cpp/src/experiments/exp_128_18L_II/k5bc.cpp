#include <iostream>
#include <random>
#include <vector>
#include <cmath>

using namespace std;

mt19937 rand_generator(time(nullptr));

#define SIZE 8

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
uniform_int_distribution<uint16_t> rand_distribution0(0,0xffff);
vector<uint8_t> getYDDT(uint32_t input_diff, uint32_t output_diff){
	uint32_t output[1<<SIZE] = {0};
	// this gets the Y_{DDT} (output values that is valid for the differential transition)
	#if SIZE == 4
	for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++){
		for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++){
			if (((v1^v2) == input_diff) && ((Sbox[v1]^Sbox[v2]) == output_diff)) output[Sbox[v1]]++;
		}
	}
	#elif SIZE == 8
	for (uint16_t v1 = 0; v1 < (1<<SIZE); v1++){
		for (uint16_t v2 = 0; v2 < (1<<SIZE); v2++){
			if (((v1^v2) == input_diff) && ((Sbox8[v1]^Sbox8[v2]) == output_diff)) output[Sbox8[v1]]++;
		}
	}
	#else 
		#error Unsupported choice setting
	#endif
	vector<uint8_t> res;
	for (int i = 0; i < (1<<SIZE); i++)
	{
		if (output[i] > 0) res.push_back(i);
	}
	return res;
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



uint8_t select(vector<uint8_t> v)
{
	uint16_t rand_num = rand_distribution0(rand_generator);
	return v[rand_num % v.size()];
}


void k5bc(uint64_t numKeys, uint64_t numTrails)
{
	// # c (0,2,20,90) (1,6,90,2) 1
	// # cn 5 (2,4,2,8) (1,6,90,2) (1,9,80,2) (2,b,0,0) (3,9,8,10) 0 (2,4 constant 0)
	// # b c (0,5,20,80) (0,2,20,90) (1,a,80,2) (1,e,90,2)
	// # b (0,5,20,80) (1,6,90,2) (1,a,80,2) (1,e,90,2) # a combination of 2 constraints above
	uniform_int_distribution<uint8_t> rand_distribution1(0,0xff);	
	uint8_t k5,kb,kc_0,kc_1,kc_1_next;
	uint8_t k5p,kbp,kcp_0,kcp_1,kcp_1_next;
	uint8_t rand_num;

	uint8_t s02,s24,s19,s05,s08;
	uint8_t s02p,s24p,s19p,s05p,s08p;
	uint8_t s16,s16p,s1a,s1ap,s2b,s2bp,s1e,s1ep,s39,s39p;
	vector<uint8_t> y2090 = getYDDT(0x20,0x90);
	vector<uint8_t> y0208 = getYDDT(0x02,0x08);
	vector<uint8_t> y9002 = getYDDT(0x90,0x02);
	vector<uint8_t> y8002 = getYDDT(0x80,0x02);
	vector<uint8_t> y2080 = getYDDT(0x20,0x80);
	vector<uint8_t> keyc = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191};
	uint32_t badKeys = 0;
	for (uint64_t i = 0; i < numKeys; i++)
	{
		rand_num = rand_distribution1(rand_generator); k5 = rand_num; k5p = k5;
		rand_num = rand_distribution1(rand_generator); kb = rand_num; kbp = kb;
		rand_num = rand_distribution1(rand_generator); kc_0 = rand_num; kcp_0 = kc_0;
		rand_num = rand_distribution1(rand_generator); kc_1 = kc_0 ^ select(keyc); kcp_1 = kc_1;
		kc_1_next = LFSR8(kc_1,2);
		kcp_1_next = LFSR8(kcp_1,2);
		uint64_t count = 0;

		for (uint64_t j = 0; j < numTrails; j++)
		{
			s02 = select(y2090); s02p = s02 ^ 0x90;
			s24 = select(y0208); s24p = s24 ^ 0x08;
			s19 = select(y8002); s19p = s19 ^ 0x02;
			s05 = select(y2080); s05p = s05 ^ 0x80;
			s08 = rand_distribution1(rand_generator); s08p = s08;
			// # c (0,2,20,90) (1,6,90,2) 1
			// # cn 5 (2,4,2,8) (1,6,90,2) (1,9,80,2) (2,b,0,0) (3,9,8,10) 0 (2,4 constant 0)
			// # b c (0,5,20,80) (0,2,20,90) (1,a,80,2) (1,e,90,2)
			// # b (0,5,20,80) (1,6,90,2) (1,a,80,2) (1,e,90,2) # a combination of 2 constraints above
			s16 = Sbox8[s02 ^ kc_0 ^ kc_1];
			s16p = Sbox8[s02p ^ kcp_0 ^ kcp_1];
			s1a = Sbox8[s05 ^ s08 ^ kb];
			s1ap = Sbox8[s05p ^ s08p ^ kbp];
			s2b = Sbox8[s16 ^ s19 ^ k5];
			s2bp = Sbox8[s16p ^ s19p ^ k5p];
			s1e = Sbox8[s02 ^ s08 ^ kc_0 ^ kc_1];
			s1ep = Sbox8[s02p ^ s08p ^ kcp_0 ^ kcp_1];
			s39 = Sbox8[s24 ^ 0 ^ s2b ^ kc_0 ^ kc_1_next];
			s39p = Sbox8[s24p ^ 0 ^ s2bp ^ kcp_0 ^ kcp_1_next];

			if ((s16 ^ s16p) != 0x02) continue;
			if ((s1a ^ s1ap) != 0x02) continue;
			if ((s1e ^ s1ep) != 0x02) continue;
			if ((s39 ^ s39p) != 0x10) continue;

			count++;
		}
		if (count > 0)
		{
			cout << log2((count+0.0)/numTrails) << endl;
		}
		else
		{
			badKeys++;
		}
	}
	cout << "badKeys: " << badKeys << endl;
}

int main()
{
	uint64_t numKeys = 1<<10;
	uint64_t numTrails = 1<<12;
	k5bc(numKeys,numTrails);
	return 0;
}