#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <ctime>

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


void k06c(uint64_t numKeys, uint64_t numTrails)
{
	// # 0 (5,1,5,5) (5,b,4,5) (5,e,5,5) (6,1,5,5) 1
	// # 0 (5,1,5,5) (6,5,5,5) 1
	// # 6 (3,2,40,4) (4,6,4,1) 1
	// # 6 (5,4,5,1),(5,0xb,4,5) (6,9,4,5) 1 *

	// # c (4,2,4,5) (5,6,5,1) 1

	// # c 0 (4,2,4,5) (3,7,40,4) (3,a,40,4) (4,8,0,0) (5,e,5,5) 0 *
	// # 6 0 (5,e,5,5) (5,4,5,1) (5,1,5,5) (6,1,5,5) (6,9,4,5) [This must hold if the 2 linear constraints above hold]
	// # 6 (5,e,5,5) (5,b,4,5) (5,4,5,1) (6,1,5,5) (6,5,5,5) (6,9,4,5) [This must hold if the 3 linear constraints above hold]
	uniform_int_distribution<uint8_t> rand_distribution1(0,0xff);	
	uint8_t kc,k0_0,k0_1,k0_1_next,k6_0,k6_1,k6_1_next;
	uint8_t kcp,k0p_0,k0p_1,k0p_1_next,k6p_0,k6p_1,k6p_1_next;


	uint8_t s51,s5b,s32,s54,s42,s37,s3a;
	uint8_t s51p,s5bp,s32p,s54p,s42p,s37p,s3ap;

	uint8_t s48,s5e,s61,s65,s46,s69,s56;
	uint8_t s48p,s5ep,s61p,s65p,s46p,s69p,s56p;

	vector<uint8_t> y0505 = getYDDT(0x05,0x05);
	vector<uint8_t> y0405 = getYDDT(0x04,0x05);
	vector<uint8_t> y4004 = getYDDT(0x40,0x04);
	vector<uint8_t> y0501 = getYDDT(0x05,0x01);
	vector<uint8_t> key0 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};
	vector<uint8_t> key6 = {0, 1, 4, 5, 8, 9, 12, 13, 16, 17, 20, 21, 24, 25, 28, 29, 32, 33, 36, 37, 40, 41, 44, 45, 48, 49, 52, 53, 56, 57, 60, 61, 64, 65, 68, 69, 72, 73, 76, 77, 80, 81, 84, 85, 88, 89, 92, 93, 96, 97, 100, 101, 104, 105, 108, 109, 112, 113, 116, 117, 120, 121, 124, 125, 128, 129, 132, 133, 136, 137, 140, 141, 144, 145, 148, 149, 152, 153, 156, 157, 160, 161, 164, 165, 168, 169, 172, 173, 176, 177, 180, 181, 184, 185, 188, 189, 192, 193, 196, 197, 200, 201, 204, 205, 208, 209, 212, 213, 216, 217, 220, 221, 224, 225, 228, 229, 232, 233, 236, 237, 240, 241, 244, 245, 248, 249, 252, 253};
	vector<uint8_t> keyc = {0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23, 32, 33, 34, 35, 36, 37, 38, 39, 48, 49, 50, 51, 52, 53, 54, 55, 64, 65, 66, 67, 68, 69, 70, 71, 80, 81, 82, 83, 84, 85, 86, 87, 96, 97, 98, 99, 100, 101, 102, 103, 112, 113, 114, 115, 116, 117, 118, 119, 128, 129, 130, 131, 132, 133, 134, 135, 144, 145, 146, 147, 148, 149, 150, 151, 160, 161, 162, 163, 164, 165, 166, 167, 176, 177, 178, 179, 180, 181, 182, 183, 192, 193, 194, 195, 196, 197, 198, 199, 208, 209, 210, 211, 212, 213, 214, 215, 224, 225, 226, 227, 228, 229, 230, 231, 240, 241, 242, 243, 244, 245, 246, 247};
	uint32_t badKeys = 0;
	uint8_t rand_num;
	for (uint64_t i = 0; i < numKeys; i++)
	{
		rand_num = rand_distribution1(rand_generator); kc = select(keyc); kcp = kc;
		rand_num = rand_distribution1(rand_generator); k0_0 = rand_num; k0p_0 = k0_0;
		rand_num = rand_distribution1(rand_generator); k0_1 = k0_0 ^ select(key0); k0p_1 = k0_1;
		rand_num = rand_distribution1(rand_generator); k6_0 = rand_num; k6p_0 = k6_0;
		rand_num = rand_distribution1(rand_generator); k6_1 = k6_0 ^ select(key6); k6p_1 = k6_1;
		k0_1_next = LFSR8(k0_1,2);
		k0p_1_next = LFSR8(k0p_1,2);
		k6_1_next = LFSR8(k6_1,2);
		k6p_1_next = LFSR8(k6p_1,2);
		uint64_t count = 0;

		for (uint64_t j = 0; j < numTrails; j++)
		{
			s51 = select(y0505); s51p = s51 ^ 0x05;
			s5b = select(y0405); s5bp = s5b ^ 0x05;
			s32 = select(y4004); s32p = s32 ^ 0x04;
			s54 = select(y0501); s54p = s54 ^ 0x01;
			s42 = select(y0405); s42p = s42 ^ 0x05;
			s37 = select(y4004); s37p = s37 ^ 0x04;
			s3a = select(y4004); s3ap = s3a ^ 0x04;

			s48 = Sbox8[k0_0 ^ k0_1 ^ s37 ^ s3a];
			s48p = Sbox8[k0p_0 ^ k0p_1 ^ s37p ^ s3ap];
			s5e = Sbox8[kc ^ s42 ^ s48 ^ 0x2];
			s5ep = Sbox8[kcp ^ s42p ^ s48p ^ 0x2];
			if ((s5e ^ s5ep) != 0x05) continue;

			s61 = Sbox8[k0_0 ^ k0_1_next ^ s51 ^ s5b ^ s5e];
			s61p = Sbox8[k0p_0 ^ k0p_1_next ^ s51p ^ s5bp ^ s5ep];
			if ((s61 ^ s61p) != 0x05) continue;
			s65 = Sbox8[k0_0 ^ k0_1_next ^ s51];
			s65p = Sbox8[k0p_0 ^ k0p_1_next ^ s51p];
			if ((s65 ^ s65p) != 0x05) continue;
			s46 = Sbox8[k6_0 ^ k6_1 ^ s32];
			s46p = Sbox8[k6p_0 ^ k6p_1 ^ s32p];
			if ((s46 ^ s46p) != 0x01) continue;
			s69 = Sbox8[k6_0 ^ k6_1_next ^ s54 ^ 0x1 ^ s5b];
			s69p = Sbox8[k6p_0 ^ k6p_1_next ^ s54p ^ 0x1 ^ s5bp];
			if ((s69 ^ s69p) != 0x05) continue;
			s56 = Sbox8[kc ^ s42];
			s56p = Sbox8[kcp ^ s42p];
			if ((s56 ^ s56p) != 0x01) continue;

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
	uint64_t numKeys = 1<<5;
	uint64_t numTrails = 1<<24;
	k06c(numKeys,numTrails);

	return 0;
}