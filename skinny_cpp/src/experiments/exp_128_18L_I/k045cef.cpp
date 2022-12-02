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


void k045cef(uint64_t numKeys, uint64_t numTrails)
{
	// # 4 (0,0,10,40) (1,4,40,6) 1 (0,0 constant = 1)
	// # c (1,2,4,6) (2,6,6,21) 1
	// # c (1,2,4,6) (1,8,5,26) (2,e,20,90) 1 (1,8 constant = 2)
	// # f (1,3,40,6) (2,7,6,21) 1
	

	// # f 5 (1,3,40,6) (0,4,10,40) (0,b,10,40) (1,9,0,0) (2,f,6,21) 0 (0,4 constant = 0)
	// # 5n e (2,6,6,21) (1,4,40,6) (1,b,4,6) (2,9,0,0) (3,b,21,20) 0 (1,4 constant = 0)
	// # c 4n (1,2,4,6) (1,8,5,26) (1,f,21,20) (2,2,0,0) (2,8,26,21) (3,e,21,20) 0 (1,8 constant = 2,  2,8 constant = 2)
	// # 0 4 (0,7,4,5) (0,0,10,40) (1,8,5,26) (1,c,40,6) (0,0 constant = 1)
	// # 0 (0,7,4,5) (1,4,40,6) (1,8,5,26) (1,c,40,6) (constant 1,8 = 2)
	uniform_int_distribution<uint8_t> rand_distribution1(0,0xff);	
	uint8_t k0,kc,ke,kf;
	uint8_t k0p,kcp,kep,kfp;

	uint8_t k4_0,k4_1,k4_1_next,k5_0,k5_1,k5_1_next;
	uint8_t k4p_0,k4p_1,k4p_1_next,k5p_0,k5p_1,k5p_1_next;

	uint8_t s00,s12,s13,s04,s0b,s1b,s1f,s28,s07,s0a;
	uint8_t s00p,s12p,s13p,s04p,s0bp,s1bp,s1fp,s28p,s07p,s0ap;

	uint8_t s18,s1c,s14,s26,s2e,s27,s19,s29,s22,s2f,s3b,s3e;
	uint8_t s18p,s1cp,s14p,s26p,s2ep,s27p,s19p,s29p,s22p,s2fp,s3bp,s3ep;


	vector<uint8_t> y1040 = getYDDT(0x10,0x40);
	vector<uint8_t> y0406 = getYDDT(0x04,0x06);
	vector<uint8_t> y4006 = getYDDT(0x40,0x06);
	vector<uint8_t> y0621 = getYDDT(0x06,0x21);
	vector<uint8_t> y2120 = getYDDT(0x21,0x20);
	vector<uint8_t> y2621 = getYDDT(0x26,0x21);
	vector<uint8_t> y0405 = getYDDT(0x04,0x05);


	vector<uint8_t> key4 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127};
	vector<uint8_t> keyc = {128, 129, 130, 131, 132, 133, 134, 135, 144, 145, 146, 147, 148, 149, 150, 151, 160, 161, 162, 163, 164, 165, 166, 167, 176, 177, 178, 179, 180, 181, 182, 183, 192, 193, 194, 195, 196, 197, 198, 199, 208, 209, 210, 211, 212, 213, 214, 215, 224, 225, 226, 227, 228, 229, 230, 231, 240, 241, 242, 243, 244, 245, 246, 247};
	vector<uint8_t> keyf = {10, 11, 12, 13, 26, 27, 28, 29, 42, 43, 44, 45, 58, 59, 60, 61, 74, 75, 76, 77, 90, 91, 92, 93, 106, 107, 108, 109, 122, 123, 124, 125, 138, 139, 140, 141, 154, 155, 156, 157, 170, 171, 172, 173, 186, 187, 188, 189, 202, 203, 204, 205, 218, 219, 220, 221, 234, 235, 236, 237, 250, 251, 252, 253};
	uint32_t badKeys = 0;
	uint8_t rand_num;
	for (uint64_t i = 0; i < numKeys; i++)
	{
		rand_num = rand_distribution1(rand_generator); k0 = rand_num; k0p = k0;
		kc = select(keyc); kcp = kc;
		rand_num = rand_distribution1(rand_generator); ke = rand_num; kep = ke;
		kf = select(keyf); kfp = kf;

		rand_num = rand_distribution1(rand_generator); k4_0 = rand_num; k4p_0 = k4_0;
		k4_1 = k4_0 ^ select(key4); k4p_1 = k4_1;
		rand_num = rand_distribution1(rand_generator); k5_0 = rand_num; k5p_0 = k5_0;
		rand_num = rand_distribution1(rand_generator); k5_1 = rand_num; k5p_1 = k5_1;

		k4_1_next = LFSR8(k4_1,2);
		k4p_1_next = LFSR8(k4p_1,2);
		k5_1_next = LFSR8(k5_1,2);
		k5p_1_next = LFSR8(k5p_1,2);
		uint64_t count = 0;

		for (uint64_t j = 0; j < numTrails; j++)
		{
			s00 = select(y1040); s00p = s00 ^ 0x40;
			s12 = select(y0406); s12p = s12 ^ 0x06;
			s13 = select(y4006); s13p = s13 ^ 0x06;
			s04 = select(y1040); s04p = s04 ^ 0x40;
			s0b = select(y1040); s0bp = s0b ^ 0x40;
			s1b = select(y0406); s1bp = s1b ^ 0x06;
			s1f = select(y2120); s1fp = s1f ^ 0x20;
			s28 = select(y2621); s28p = s28 ^ 0x21;
			s07 = select(y0405); s07p = s07 ^ 0x05;			
			s0a = rand_distribution1(rand_generator); s0ap = s0a;
			// # 4 (0,0,10,40) (1,4,40,6) 1 (0,0 constant = 1)
			// # c (1,2,4,6) (2,6,6,21) 1
			// # c (1,2,4,6) (1,8,5,26) (2,e,20,90) 1 (1,8 constant = 2)
			// # f (1,3,40,6) (2,7,6,21) 1
			
			// # f 5 (1,3,40,6) (0,4,10,40) (0,b,10,40) (1,9,0,0) (2,f,6,21) 0 (0,4 constant = 0)
			// # 5n e (2,6,6,21) (1,4,40,6) (1,b,4,6) (2,9,0,0) (3,b,21,20) 0 (1,4 constant = 0)
			// # c 4n (1,2,4,6) (1,8,5,26) (1,f,21,20) (2,2,0,0) (2,8,26,21) (3,e,21,20) 0 (1,8 constant = 2,  2,8 constant = 2)
			// # 0 4 (0,7,4,5) (0,0,10,40) (1,8,5,26) (1,c,40,6) (0,0 constant = 1)
			// # 0 (0,7,4,5) (1,4,40,6) (1,8,5,26) (1,c,40,6) (constant 1,8 = 2)
			s18 = Sbox8[s07 ^ s0a ^ k0];
			s18p = Sbox8[s07p ^ s0ap ^ k0p];
			if ((s18 ^ s18p) != 0x26) continue;
			s1c = Sbox8[s00 ^ 1 ^ s0a ^ k4_0 ^ k4_1];
			s1cp = Sbox8[s00p ^ 1 ^ s0ap ^ k4p_0 ^ k4p_1];
			if ((s1c ^ s1cp) != 0x06) continue;

			s14 = Sbox8[s00 ^ 1 ^ k4_0 ^ k4_1];
			s14p = Sbox8[s00p ^ 1 ^ k4p_0 ^ k4p_1];
			if ((s14 ^ s14p) != 0x06) continue;
			s26 = Sbox8[s12 ^ kc];
			s26p = Sbox8[s12p ^ kcp];
			if ((s26 ^ s26p) != 0x21) continue;
			s2e = Sbox8[s12 ^ s18 ^ 2 ^ kc];
			s2ep = Sbox8[s12p ^ s18p ^ 2 ^ kcp];
			if ((s2e ^ s2ep) != 0x90) continue;
			s27 = Sbox8[s13 ^ kf];
			s27p = Sbox8[s13p ^ kfp];
			if ((s27 ^ s27p) != 0x21) continue;

			// # f 5 (1,3,40,6) (0,4,10,40) (0,b,10,40) (1,9,0,0) (2,f,6,21) 0 (0,4 constant = 0)
			// # 5n e (2,6,6,21) (1,4,40,6) (1,b,4,6) (2,9,0,0) (3,b,21,20) 0 (1,4 constant = 0)
			// # c 4n (1,2,4,6) (1,8,5,26) (1,f,21,20) (2,2,0,0) (2,8,26,21) (3,e,21,20) 0 (1,8 constant = 2,  2,8 constant = 2)
			// # 0 4 (0,7,4,5) (0,0,10,40) (1,8,5,26) (1,c,40,6) (0,0 constant = 1)
			// # 0 (0,7,4,5) (1,4,40,6) (1,8,5,26) (1,c,40,6) (constant 1,8 = 2)
			s19 = Sbox8[s04 ^ 0 ^ s0b ^ k5_0 ^ k5_1];
			s19p = Sbox8[s04p ^ 0 ^ s0bp ^ k5p_0 ^ k5p_1];
			s29 = Sbox8[s14 ^ 0 ^ s1b ^ ke];
			s29p = Sbox8[s14p ^ 0 ^ s1bp ^ kep];
			s22 = Sbox8[s12 ^ s18 ^ 2 ^ s1f ^ kc];
			s22p = Sbox8[s12p ^ s18p ^ 2 ^ s1fp ^ kcp];
			s2f = Sbox8[s19 ^ s13 ^ kf];
			s2fp = Sbox8[s19p ^ s13p ^ kfp];
			if ((s2f ^ s2fp) != 0x21) continue;
			s3b = Sbox8[s26 ^ s29 ^ k5_0 ^ k5_1_next];
			s3bp = Sbox8[s26p ^ s29p ^ k5p_0 ^ k5p_1_next];
			if ((s3b ^ s3bp) != 0x20) continue;
			s3e = Sbox8[s22 ^ s28 ^ 2 ^ k4_0 ^ k4_1_next];
			s3ep = Sbox8[s22p ^ s28p ^ 2 ^ k4p_0 ^ k4p_1_next];
			if ((s3e ^ s3ep) != 0x20) continue;

			count++;
		}
		cout << count << endl;
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
	uint64_t numKeys = 1<<4;
	uint64_t numTrails = 1<<28;
	k045cef(numKeys,numTrails);

	return 0;
}