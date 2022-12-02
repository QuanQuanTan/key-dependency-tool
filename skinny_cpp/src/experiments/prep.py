import numpy as np
from secrets import randbits
from collections import Counter
import itertools

Sbox4 = [0xc,0x6,0x9,0x0,0x1,0xa,0x2,0xb,0x3,0x8,0x5,0xd,0x4,0xe,0x7,0xf]
Sbox8 = [0x65, 0x4c, 0x6a, 0x42, 0x4b, 0x63, 0x43, 0x6b, 0x55, 0x75, 0x5a, 0x7a, 0x53, 0x73, 0x5b, 0x7b,
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
							0xe2, 0xca, 0xee, 0xc6, 0xcf, 0xe7, 0xc7, 0xef, 0xd2, 0xf2, 0xde, 0xfe, 0xd7, 0xf7, 0xdf, 0xff]
invSbox8 = [172, 232, 104, 60, 108, 56, 168, 236, 170, 174, 58, 62, 106, 110, 234, 238, 166, 163, 51, 54, 102, 99, 227, 230, 225, 164, 97, 52, 49, 100, 161, 228, 141, 201, 73, 29, 77, 25, 137, 205, 139, 143, 27, 31, 75, 79, 203, 207, 133, 192, 64, 21, 69, 16, 128, 197, 130, 135, 18, 23, 66, 71, 194, 199, 150, 147, 3, 6, 86, 83, 211, 214, 209, 148, 81, 4, 1, 84, 145, 212, 156, 216, 88, 12, 92, 8, 152, 220, 154, 158, 10, 14, 90, 94, 218, 222, 149, 208, 80, 5, 85, 0, 144, 213, 146, 151, 2, 7, 82, 87, 210, 215, 157, 217, 89, 13, 93, 9, 153, 221, 155, 159, 11, 15, 91, 95, 219, 223, 22, 19, 131, 134, 70, 67, 195, 198, 65, 20, 193, 132, 17, 68, 129, 196, 28, 72, 200, 140, 76, 24, 136, 204, 26, 30, 138, 142, 74, 78, 202, 206, 53, 96, 224, 165, 101, 48, 160, 229, 50, 55, 162, 167, 98, 103, 226, 231, 61, 105, 233, 173, 109, 57, 169, 237, 59, 63, 171, 175, 107, 111, 235, 239, 38, 35, 179, 182, 118, 115, 243, 246, 113, 36, 241, 180, 33, 116, 177, 244, 44, 120, 248, 188, 124, 40, 184, 252, 42, 46, 186, 190, 122, 126, 250, 254, 37, 112, 240, 181, 117, 32, 176, 245, 34, 39, 178, 183, 114, 119, 242, 247, 45, 121, 249, 189, 125, 41, 185, 253, 43, 47, 187, 191, 123, 127, 251, 255]

def getSbox8(x):
	return Sbox8[x]

def getSbox4(x):
	return Sbox4[x]

def getInvSbox8(x):
	return invSbox8[x]

def getDDT4():
	ddt = np.zeros((16,16),dtype=int)
	for i in range(16):
		for j in range(16):
			ddt[i^j][Sbox4[i]^Sbox4[j]] += 1
	return ddt

def getDDT8():
	ddt = np.zeros((256,256),dtype=int)
	for i in range(256):
		for j in range(256):
			ddt[i^j][Sbox8[i]^Sbox8[j]] += 1
	return ddt

def getXDDT8():
	xddt = [[[] for _ in range(256)] for _ in range(256)]
	for i in range(256):
		for j in range(256):
			xddt[i^j][Sbox8[i]^Sbox8[j]].append(i)
	return xddt


def getYDDT8():
	yddt = [[[] for _ in range(256)] for _ in range(256)]
	for i in range(256):
		for j in range(256):
			yddt[i^j][Sbox8[i]^Sbox8[j]].append(Sbox8[i])
	return yddt



def getXDDT4():
	xddt = [[[] for _ in range(16)] for _ in range(16)]
	for i in range(16):
		for j in range(16):
			xddt[i^j][Sbox4[i]^Sbox4[j]].append(i)
	return xddt


def getYDDT4():
	yddt = [[[] for _ in range(16)] for _ in range(16)]
	for i in range(16):
		for j in range(16):
			yddt[i^j][Sbox4[i]^Sbox4[j]].append(Sbox4[i])
	return yddt

def XOR(A,B):
	return [a^b for a in A for b in B]


def LFSR4(n,v):
	x3 = (n >> 3) % 2;
	x2 = (n >> 2) % 2;
	x1 = (n >> 1) % 2;
	x0 = (n >> 0) % 2;
	if (v == 2): return (x2 << 3) ^ (x1 << 2) ^ (x0 << 1) ^ ((x2 ^ x3) << 0);
	if (v == 3): return ((x0 ^ x3) << 3) ^ (x3 << 2) ^ (x2 << 1) ^ (x1 << 0);
	if (v == 4): return (x2 << 3) ^ (x1 << 2) ^ ((x0 ^ x2) << 1) ^ ((x1 ^ x2 ^ x3) << 0);


def LFSR8(n,v):
	x7 = (n >> 7) % 2;
	x6 = (n >> 6) % 2;
	x5 = (n >> 5) % 2;
	x4 = (n >> 4) % 2;
	x3 = (n >> 3) % 2;
	x2 = (n >> 2) % 2;
	x1 = (n >> 1) % 2;
	x0 = (n >> 0) % 2;
	if (v == 2): return (x6 << 7) ^ (x5 << 6) ^ (x4 << 5) ^ (x3 << 4) ^ (x2 << 3) ^ (x1 << 2) ^ (x0 << 1) ^ ((x7 ^ x5) << 0);
	elif (v == 3): return ((x0 ^ x6) << 7) ^ (x7 << 6) ^ (x6 << 5) ^ (x5 << 4) ^ (x4 << 3) ^ (x3 << 2) ^ (x2 << 1) ^ (x1 << 0);
	elif (v == 4): return 0;

def select(lst):
	x = randbits(8)
	return lst[x % len(lst)]

def Sbox8Pair(s0,s1):
	return Sbox8[s0],Sbox8[s1]

def Sbox4Pair(s0,s1):
	return Sbox4[s0],Sbox4[s1]

def computeLinearProb8(left,right,keys=[i for i in range(256)],constant=[0,0,0]):
	XDDT8 = getXDDT8()
	YDDT8 = getYDDT8()
	# find out how many of the them allow left to go
	# to the right
	rvalues = XDDT8[right[2]][right[3]]
	count = [0 for _ in range(len(keys))]
	if len(left) == 1:
		lvalues0 = YDDT8[left[0][2]][left[0][3]]; lvalues0 = [l^constant[0] for l in lvalues0]
		for k in keys:
			for l0 in lvalues0:
				if (k ^ l0) in rvalues: count[k] += 1
		distribution = [-np.log2(c/len(lvalues0)) for c in count if c > 0]
	elif len(left) == 2:
		lvalues0 = YDDT8[left[0][2]][left[0][3]]; lvalues0 = [l^constant[0] for l in lvalues0]
		lvalues1 = YDDT8[left[1][2]][left[1][3]]; lvalues1 = [l^constant[1] for l in lvalues1]
		for k in keys:
			for l0 in lvalues0:
				for l1 in lvalues1:
					if k ^ l0 ^ l1 in rvalues: count[k] += 1
		distribution = [-np.log2(c/len(lvalues0)/len(lvalues1)) for c in count if c > 0]
	elif len(left) == 3:
		lvalues0 = YDDT8[left[0][2]][left[0][3]]; lvalues0 = [l^constant[0] for l in lvalues0]
		lvalues1 = YDDT8[left[1][2]][left[1][3]]; lvalues1 = [l^constant[1] for l in lvalues1]
		lvalues2 = YDDT8[left[2][2]][left[2][3]]; lvalues2 = [l^constant[2] for l in lvalues2]
		for k in keys:
			for l0 in lvalues0:
				for l1 in lvalues1:
					for l2 in lvalues2:
						if k ^ l0 ^ l1 ^ l2 in rvalues: count[k] += 1
		distribution = [-np.log2(c/len(lvalues0)/len(lvalues1)/len(lvalues2)) for c in count if c > 0]
	
	possibleKeys =  [i for i in range(256) if count[i] > 0]
	distribution = Counter(distribution)
	total = sum(distribution.values())
	for k,v,in distribution.items():
		distribution[k] = v/total
	return possibleKeys,distribution

def computeLinearProb4(left,right,keys=[i for i in range(16)],constant=[0,0,0]):
	XDDT4 = getXDDT4()
	YDDT4 = getYDDT4()
	# find out how many of the them allow left to go
	# to the right
	rvalues = XDDT4[right[2]][right[3]]
	count = [0 for _ in range(len(keys))]
	if len(left) == 1:
		lvalues0 = YDDT4[left[0][2]][left[0][3]]; lvalues0 = [l^constant[0] for l in lvalues0]
		for k in keys:
			for l0 in lvalues0:
				if (k ^ l0) in rvalues: count[k] += 1
		distribution = [-np.log2(c/len(lvalues0)) for c in count if c > 0]
	elif len(left) == 2:
		lvalues0 = YDDT4[left[0][2]][left[0][3]]; lvalues0 = [l^constant[0] for l in lvalues0]
		lvalues1 = YDDT4[left[1][2]][left[1][3]]; lvalues1 = [l^constant[1] for l in lvalues1]
		for k in keys:
			for l0 in lvalues0:
				for l1 in lvalues1:
					if k ^ l0 ^ l1 in rvalues: count[k] += 1
		distribution = [-np.log2(c/len(lvalues0)/len(lvalues1)) for c in count if c > 0]
	elif len(left) == 3:
		lvalues0 = YDDT4[left[0][2]][left[0][3]]; lvalues0 = [l^constant[0] for l in lvalues0]
		lvalues1 = YDDT4[left[1][2]][left[1][3]]; lvalues1 = [l^constant[1] for l in lvalues1]
		lvalues2 = YDDT4[left[2][2]][left[2][3]]; lvalues2 = [l^constant[2] for l in lvalues2]
		for k in keys:
			for l0 in lvalues0:
				for l1 in lvalues1:
					for l2 in lvalues2:
						if k ^ l0 ^ l1 ^ l2 in rvalues: count[k] += 1
		distribution = [-np.log2(c/len(lvalues0)/len(lvalues1)/len(lvalues2)) for c in count if c > 0]
	
	possibleKeys =  [i for i in range(16) if count[i] > 0]
	distribution = Counter(distribution)
	total = sum(distribution.values())
	for k,v,in distribution.items():
		distribution[k] = v/total
	return possibleKeys,distribution

def overlapKeys(*keys):
	possibleKeys = []
	for k in range(256):
		flag = True
		for key in keys:
			if k not in key:
				flag = False
				break
		if flag == True:
			possibleKeys.append(k)
	return possibleKeys

def combineDictionary(lst):
	# converting the list of dictionaries to list of (lists of tuples)
	newLst = []
	for i in range(len(lst)):
		if lst[i] == 0: continue
		newD = []
		for k,v in lst[i].items(): newD.append((k,v))
		newLst.append(newD)
	newDictionary = {}
	for t in itertools.product(*newLst):
		prob = 0
		proportion = 1
		for element in t: 
			prob += element[0]
			proportion *= element[1]
		if prob not in newDictionary.keys():
			newDictionary[prob] = proportion
		else:
			newDictionary[prob] += proportion
	return newDictionary


	
