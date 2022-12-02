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


void k28(uint64_t numKeys, uint64_t numTrails)
{
	// # 8 (0,1,20,80) (0,b,20,80) (0,e,40,4) (1,1,4,5) 1
	// # 8 (0,1,20,80) (1,5,80,2) 1
	// # 2 (1,0,80,2) (2,4,2,8) 1 (1,0 constant 3)
	// # 2 8 8n (1,0,80,2) (1,a,80,2) (0,1,20,80) (0,b,20,80) (1,d,0,0) (2,0,0,0) (2,a,2,8) (2,d,2,c) (3,0,4,6) 0 (1,0 constant=3 , 2,0 constant = 7 )
	// # 2 8 8n (1,0,80,2) (1,a,80,2) (0,1,20,80) (0,b,20,80) (1,d,0,0) (2,0,0,0) (2,a,2,8) (3,c,8,10) 0 (1,0 constant=3 , 2,0 constant = 7 )
	uniform_int_distribution<uint8_t> rand_distribution1(0,0xff);	
	uint8_t k8_0,k8_1,k8_1_next;
	uint8_t k8p_0,k8p_1,k8p_1_next;
	uint8_t k2_0,k2_1;
	uint8_t k2p_0,k2p_1;
	uint8_t rand_num;


	uint8_t s01,s0b,s0e,s1a,s10,s2a,s2d;
	uint8_t s01p,s0bp,s0ep,s1ap,s10p,s2ap,s2dp;
	uint8_t s11,s11p,s15,s15p,s24,s24p,s1d,s1dp,s20,s20p,s30,s30p,s3c,s3cp;
	vector<uint8_t> y2080 = getYDDT(0x20,0x80);
	vector<uint8_t> y4004 = getYDDT(0x40,0x04);
	vector<uint8_t> y8002 = getYDDT(0x80,0x02);
	vector<uint8_t> y0208 = getYDDT(0x02,0x08);
	vector<uint8_t> y020c = getYDDT(0x02,0x0c);
	vector<uint8_t> key2 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};
	vector<uint8_t> key8 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191};
	uint32_t badKeys = 0;
	for (uint64_t i = 0; i < numKeys; i++)
	{
		rand_num = rand_distribution1(rand_generator); k8_0 = rand_num; k8p_0 = k8_0;
		k8_1 = k8_0 ^ select(key8); k8p_1 = k8_1;
		rand_num = rand_distribution1(rand_generator); k2_0 = rand_num; k2p_0 = k2_0;
		k2_1 = k2_0 ^ select(key2); k2p_1 = k2_1;
		k8_1_next = LFSR8(k8_1,2);
		k8p_1_next = LFSR8(k8p_1,2);
		uint64_t count = 0;
		
		for (uint64_t j = 0; j < numTrails; j++)
		{
			s01 = select(y2080); s01p = s01 ^ 0x80;
			s0b = select(y2080); s0bp = s0b ^ 0x80;
			s0e = select(y4004); s0ep = s0e ^ 0x04;
			s1a = select(y8002); s1ap = s1a ^ 0x02;
			s10 = select(y8002); s10p = s10 ^ 0x02;
			s2a = select(y0208); s2ap = s2a ^ 0x08;
			s2d = select(y020c); s2dp = s2d ^ 0x0c;
			// # 8 (0,1,20,80) (0,b,20,80) (0,e,40,4) (1,1,4,5) 1
			// # 8 (0,1,20,80) (1,5,80,2) 1
			// # 2 (1,0,80,2) (2,4,2,8) 1 (1,0 constant 3)
			// # 2 8 8n (1,0,80,2) (1,a,80,2) (0,1,20,80) (0,b,20,80) (1,d,0,0) (2,0,0,0) (2,a,2,8) (2,d,2,c) (3,0,4,6) 0 (1,0 constant=3 , 2,0 constant = 7 )
			// # 2 8 8n (1,0,80,2) (1,a,80,2) (0,1,20,80) (0,b,20,80) (1,d,0,0) (2,0,0,0) (2,a,2,8) (3,c,8,10) 0 (1,0 constant=3 , 2,0 constant = 7 )
			s11 = Sbox8[s01 ^ s0b ^ s0e ^ k8_0 ^ k8_1]; s11p = Sbox8[s01p ^ s0bp ^ s0ep ^ k8p_0 ^ k8p_1];
			if ((s11 ^ s11p) != 0x05) continue;
			s15 = Sbox8[s01 ^ k8_0 ^ k8_1]; s15p = Sbox8[s01p ^ k8p_0 ^ k8p_1];
			if ((s15 ^ s15p) != 0x02) continue;
			s24 = Sbox8[s10 ^ 3 ^ k2_0 ^ k2_1]; s24p = Sbox8[s10p ^ 3 ^ k2p_0 ^ k2p_1];
			if ((s24 ^ s24p) != 0x08) continue;
			s1d = Sbox8[s01 ^ s0b ^ k8_0 ^ k8_1]; s1dp = Sbox8[s01p ^ s0bp ^ k8p_0 ^ k8p_1];
			s20 = Sbox8[s10 ^ 3 ^ s1a ^ s1d ^ k2_0 ^ k2_1]; s20p = Sbox8[s10p ^ 3 ^ s1ap ^ s1dp ^ k2p_0 ^ k2p_1];
			s30 = Sbox8[s2d ^ s20 ^ 7 ^ s2a ^ k8_0 ^ k8_1_next]; s30p = Sbox8[s2dp ^ s20p ^ 7 ^ s2ap ^ k8p_0 ^ k8p_1_next];
			if ((s30 ^ s30p) != 0x06) continue;
			s3c = Sbox8[s2a ^ s20 ^ 7 ^ k8_0 ^ k8_1_next]; s3cp = Sbox8[s2ap ^ s20p ^ 7 ^ k8p_0 ^ k8p_1_next];
			if ((s3c ^ s3cp) != 0x10) continue;
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
	uint64_t numKeys = 1<<8; 
	uint64_t numTrails = 1<<24;
	k28(numKeys,numTrails);

	return 0;
}