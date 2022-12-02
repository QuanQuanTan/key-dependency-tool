import numpy as np
import secrets
from collections import Counter
import matplotlib.pyplot as plt

SIZE = 8
Sbox = [0xc,0x6,0x9,0x0,0x1,0xa,0x2,0xb,0x3,0x8,0x5,0xd,0x4,0xe,0x7,0xf]
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

def getSbox8(x):
	return Sbox8[x]

def getDDT4():
	ddt = np.zeros((16,16),dtype=int)
	for i in range(16):
		for j in range(16):
			ddt[i^j][Sbox[i]^Sbox[j]] += 1
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



def getXDDT():
	xddt = [[[] for _ in range(16)] for _ in range(16)]
	for i in range(16):
		for j in range(16):
			xddt[i^j][Sbox[i]^Sbox[j]].append(i)
	return xddt


def getYDDT():
	yddt = [[[] for _ in range(16)] for _ in range(16)]
	for i in range(16):
		for j in range(16):
			yddt[i^j][Sbox[i]^Sbox[j]].append(Sbox[i])
	return yddt

def XOR(A,B):
	return [a^b for a in A for b in B]


def LFSR(n,v):
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
	x = secrets.randbits(int(np.log2(len(lst))))
	return lst[x % len(lst)]


def computeLinearKeys(left,right):
	# compute the number of keys
	keys = []
	for k in range(1<<8):
		for l in left:
			if k ^ l in right: 
				keys.append(k)
	return list(set(keys))

def computeLinearProb(*argv,debug=False):
	keys = argv[-1]
	# find the prob
	# we assume that for all keys, the prob should be the same
	length = (len(argv)-1)//2
	counts = {}
	keys = [[0 for _ in range(256)] for _ in range(length)]
	for i in range(length):
		left = argv[2*i]; right = argv[2*i+1]
		for k in range(256):
			for l in left:
				if k ^ l in right:
					keys[i][k] += 1/len(left)

	key = []
	for i in range(256):
		tmp = 1
		for j in range(length):
			tmp *= keys[j][i]
		key.append(tmp)	

	counts = Counter(key)
	new_counts = {}
	s = sum(counts.values())
	for k,v in counts.items():
		if k == 0: continue
		new_counts[-np.log2(k)] = v/s
	return new_counts # prob, distribution




XDDT = getXDDT()
YDDT = getYDDT()
XDDT8 = getXDDT8()
YDDT8 = getYDDT8()
DDT8 = getDDT8()
DDT4 = getDDT4()


# linear
# (0,3,25,25) (0,9,50,4) (0,c,15,61) (1,3,40,6) 14		# (4,1)							
# (0,0,10,40) (1,4,40,6) 14								# (4,4)
# (0,3,25,25) (1,7,25,a5) 14							# (4,1)
# (0,3,25,25) (0,9,50,4) (1,f,21,20) 14					# (4,1)		
# (1,1,40,6) (2,5,6,21) 14								# (4,8)
# (1,2,4,6) (2,6,6,21) 14									# (4,c)
# (1,3,40,6) (2,7,6,21) 14								# (4,f)
# (1,2,4,6) (1,8,5,26) (2,e,20,90) 14						# (4,c)		
# (2,0,6,21) (3,4,21,20) 14								# (4,2)
# (3,1,90,3) (3,b,21,20) (4,d,20,90) 14					# (4,9)			

# nonlinear
# (1,3,40,6) (0,4,10,40) (0,b,10,40) (1,9,0,0) (2,f,6,21) 0 				# (4,f) (4,5)
# (2,6,6,21) (1,4,40,6) (1,b,4,6) (2,9,0,0) (3,b,21,20) 0					# (4,5) (4,e)
# (1,2,4,6) (1,8,5,26) (1,f,21,20) (2,2,0,0) (2,8,26,21) (3,e,21,20) 0 	# (4,c) (4,4)

# same round
# doesn't restrict anything. ignored

probs = [0 for _ in range(16)]
keys = 0
def unrestricted():
	global probs,keys
	# these keys are not restricted at all
	# 0,3,6,7,a,b,d
	print('For this experiment, we move all keys to round 4 (including the key positions)')
	print('The keys at position 0,3,6,7,a,b,d are free')
	probs[0x0] = {0:1}
	probs[0x3] = {0:1}
	probs[0x6] = {0:1}
	probs[0x7] = {0:1}
	probs[0xa] = {0:1}
	probs[0xb] = {0:1}
	probs[0xd] = {0:1}
	keys += 0



def key1():
	global probs,keys
	# # key 1
	# (0,3,25,25) (0,9,50,4) (0,c,15,61) (1,3,40,6) 14		# (4,1)
	# (0,3,25,25) (1,7,25,a5) 14								# (4,1)
	# (0,3,25,25) (0,9,50,4) (1,f,21,20) 14					# (4,1)
	left0 = XOR(XOR(YDDT8[0x25][0x25],YDDT8[0x50][0x04]),YDDT8[0x15][0x61])
	right0 = XDDT8[0x40][0x06]
	left1 = YDDT8[0x25][0x25]
	right1 = XDDT8[0x25][0xa5]
	left2 = XOR(YDDT8[0x25][0x25],YDDT8[0x50][0x04])
	right2 = XDDT8[0x21][0x20]

	k0 = computeLinearKeys(left0,right0)
	k1 = computeLinearKeys(left1,right1)
	k2 = computeLinearKeys(left2,right2)
	k = []
	for i in range(256):
		if i in k0 and i in k1 and i in k2: k.append(i)	
	prob = computeLinearProb(left0,right0,left1,right1,left2,right2,k)
	key = 8 - np.log2(len(k))
	print('key1:',np.log2(key))
	
	probs[1] = prob
	keys += key



def key2():
	global probs,keys
	# # key 2
	# (2,0,6,21) (3,4,21,20) 14								# (4,2)
	left = YDDT8[0x6][0x21]
	right = XDDT8[0x21][0x20]
	k = computeLinearKeys(left,right)
	prob = computeLinearProb(left,right,k)
	key = 8 - np.log2(len(k))

	print('key2:',key)
	probs[2] = prob
	keys += key


def key4c(RUN=False):
	global probs,keys
	# # key 4,c
	# (0,0,10,40) (1,4,40,6) 14									# (4,4)
	# (1,2,4,6) (2,6,6,21) 14									# (4,c)
	# (1,2,4,6) (1,8,5,26) (2,e,20,90) 14						# (4,c)		
	# (1,2,4,6) (1,8,5,26) (1,f,21,20) (2,2,0,0) (2,8,26,21) (3,e,21,20) 0 	# (4,c) (4,4)
	if not RUN:
		
		key = 16 - np.log2(128*64)
		probs[0xc] = {5.83:0.25, 7.5:0.5,9:0.25}
		probs[0x4] = {0:1}
		print('key 4,c:',np.log2(128*64))
		keys += key
	else:
		key44 = computeLinearKeys(YDDT8[0x10][0x40],XDDT8[0x40][0x06])
		key4c_0 = computeLinearKeys(YDDT8[0x04][0x06],XDDT8[0x06][0x21])
		key4c_1 = computeLinearKeys(XOR(YDDT8[0x04][0x06],YDDT8[0x05][0x26]),XDDT8[0x20][0x90])
		key4c = []
		for i in key4c_0:
			for j in key4c_1:
				if i == j:
					key4c.append(i)
		runs = 1 << 14
		k = 1 << 12
		counts = []
		for _ in range(k):
			k0_tk1 = secrets.randbits(8) # (4,c)
			k0_tk2 = k0_tk1 ^ key4c[secrets.randbits(16) % len(key4c)] # (4,c)
			k1_tk1 = secrets.randbits(8) # (4,4)
			k1_tk2 = k1_tk1 ^ key44[secrets.randbits(16) % len(key44)] # (4,4)
			count = 0
			for _ in range(runs):
				state00 = YDDT8[0x10][0x40][secrets.randbits(16)%len(YDDT8[0x10][0x40])]; state00p = state00 ^ 0x40
				state12 = YDDT8[0x04][0x06][secrets.randbits(16)%len(YDDT8[0x04][0x06])]; state12p = state12 ^ 0x06
				state18 = YDDT8[0x05][0x26][secrets.randbits(16)%len(YDDT8[0x05][0x26])] ^ 0x02; state18p = state18 ^ 0x26 # constants
				state1f = YDDT8[0x21][0x20][secrets.randbits(16)%len(YDDT8[0x21][0x20])]; state1fp = state1f ^ 0x20
				state28 = YDDT8[0x26][0x21][secrets.randbits(16)%len(YDDT8[0x26][0x21])] ^ 0x02; state28p = state28 ^ 0x21 # constants

				state14 = Sbox8[state00 ^ k1_tk1 ^ k1_tk2]; state14p = Sbox8[state00p ^ k1_tk1 ^ k1_tk2]; 
				if state14 ^ state14p != 0x06: continue
				k1_tk2new = LFSR8(k1_tk2,2) # convert from (0,0) to (2,2)
				if Sbox8[state12 ^ k0_tk1 ^ k0_tk2] ^ Sbox8[state12p ^ k0_tk1 ^ k0_tk2] != 0x21: continue
				if Sbox8[state12 ^ state18 ^ k0_tk1 ^ k0_tk2] ^ Sbox8[state12p ^ state18p ^ k0_tk1 ^ k0_tk2] != 0x90: continue
				
				state22 = Sbox8[state12 ^ state18 ^ state1f ^ k0_tk1 ^ k0_tk2]; state22p = Sbox8[state12p ^ state18p ^ state1fp ^ k0_tk1 ^ k0_tk2];
				if state22 ^ state22p != 0:
					assert False	
				state3e = Sbox8[state22 ^ state28 ^ k1_tk1 ^ k1_tk2new]; state3ep = Sbox8[state22p ^ state28p ^ k1_tk1 ^ k1_tk2new];
				if state3e ^ state3ep != 0x20: continue
				count += 1
			if count != 0:
				counts.append(-np.log2((count/runs)))
		counts.sort()
		plt.plot(counts)
		plt.savefig('experiment_18L1_k4c.png')
		# print('number of keys:',-np.log2(len(counts)/k))
		totalKeys = len(key44)*len(key4c)*(len(counts)/k)
		print('number of keys for 4 and c:', np.log2(totalKeys))

def key8():
	global probs,keys
	# key 8
	#(1,1,40,6) (2,5,6,21) 14								# (4,8)
	left = YDDT8[0x40][0x06]
	right = XDDT8[0x06][0x21]

	k = computeLinearKeys(left,right)
	prob = computeLinearProb(left,right,k)
	key = 8 - np.log2(len(k))
	print('key8:',key)
	probs[8] = prob
	keys += key

def key9():
	global probs,keys
	# key 9
	# (3,1,90,3) (3,b,21,20) (4,d,20,90) 14					# (4,9)			
	left = XOR(YDDT8[0x90][0x03],YDDT8[0x21][0x20])
	right = XDDT8[0x20][0x90]

	k = computeLinearKeys(left,right)
	prob = computeLinearProb(left,right,k,debug=True)
	key = 8 - np.log2(len(k))
	print('key9:',key)
	probs[9] = prob
	keys += key


def key5ef(RUN=False):
	global probs, keys
	key4f = computeLinearKeys(YDDT8[0x40][0x06],XDDT8[0x06][0x21])
	key45 = [i for i in range(256)]
	key4e = [i for i in range(256)]
	if not RUN:
		key = 24-np.log2(len(key4f)*len(key45)*len(key4e))
		print('key 5,e,f: ',key)
		probs[0x5] = {4.6091389865467285: 6.4453125,5.163159819112441: 17.2607421875,5.725482960145787: 20.2392578125,6.42055063503151: 23.974609375}
		probs[0xe] = {0:1}
		probs[0xf] = {0:1}
		keys += key
	else:
		########(1,3,40,6) (2,7,6,21) 14												# (4,f)
		########(1,3,40,6) (0,4,10,40) (0,b,10,40) (1,9,0,0) (2,f,6,21) 0 				# (4,f) (4,5)
		########(2,6,6,21) (1,4,40,6) (1,b,4,6) (2,9,0,0) (3,b,21,20) 0					# (4,5) (4,e)
		counts = []
		runs = 1 << 14
		k = 1 << 12
		for _ in range(k):
			k0_tk1 = secrets.randbits(8) # (4,f)
			k0_tk2 = k0_tk1 ^ key4f[secrets.randbits(16) % len(key4f)] # (4,f)
			k1_tk1 = secrets.randbits(8) # (4,5)
			k1_tk2 = k1_tk1 ^ key45[secrets.randbits(16) % len(key45)] # (4,5)
			k2_tk1 = secrets.randbits(8) # (4,e)
			k2_tk2 = k2_tk1 ^ key4e[secrets.randbits(16) % len(key4e)] # (4,e)
			
			count = 0
			for _ in range(runs):
				state13 = YDDT8[0x40][0x06][secrets.randbits(16)%len(YDDT8[0x40][0x06])]; state13p = state13 ^ 0x06
				state04 = YDDT8[0x10][0x40][secrets.randbits(16)%len(YDDT8[0x10][0x40])]; state04p = state04 ^ 0x40
				state0b = YDDT8[0x10][0x40][secrets.randbits(16)%len(YDDT8[0x10][0x40])]; state0bp = state0b ^ 0x40
				state26 = YDDT8[0x06][0x21][secrets.randbits(16)%len(YDDT8[0x06][0x21])]; state26p = state26 ^ 0x21
				state14 = YDDT8[0x40][0x06][secrets.randbits(16)%len(YDDT8[0x40][0x06])]; state14p = state14 ^ 0x06
				state1b = YDDT8[0x04][0x06][secrets.randbits(16)%len(YDDT8[0x04][0x06])]; state1bp = state1b ^ 0x06

				state27 = Sbox8[state13 ^ k0_tk1 ^ k0_tk2]; state27p = Sbox8[state13p ^ k0_tk1 ^ k0_tk2]; 
				if state27 ^ state27p != 0x21: continue
				state19 = Sbox8[state04 ^ state0b ^ k1_tk1 ^ k1_tk2]; state19p = Sbox8[state04p ^ state0bp ^ k1_tk1 ^ k1_tk2]
				state2f = Sbox8[state19 ^ state13 ^ k2_tk1 ^ k2_tk2]; state2fp = Sbox8[state19p ^ state13p ^ k2_tk1 ^ k2_tk2]; 
				if state2f ^ state2fp != 0x21: continue
				k1_tk2new = LFSR8(k1_tk2,2) # convert from (0,4) to (2,6)
				state29 = Sbox8[state14 ^ state1b ^ k2_tk1 ^ k2_tk2]; state29p = Sbox8[state14p ^ state1bp ^ k2_tk1  ^ k2_tk2];
				state3b = Sbox8[state29 ^ state26 ^ k1_tk1 ^ k1_tk2new]; state3bp = Sbox8[state29p ^ state26p ^ k1_tk1 ^ k1_tk2new];
				if state3b ^ state3bp != 0x20: continue
				count += 1
			if count != 0:
				counts.append(-np.log2((count/runs)))
		counts.sort()
		plt.plot(counts)
		plt.savefig('experiment_18L1_k5ef.png')
		print('number of keys:',-np.log2(len(counts)/k))
		totalKeys = len(key4f)*len(key45)*len(key4e) * (len(counts)/k)
		print('number of keys for 5ef:', np.log2(totalKeys))

if __name__ == '__main__':
	unrestricted();
	key1();
	key2();
	key4c(RUN=True);
	key8();
	key9();
	key5ef(RUN=True);
	print('size of key:',-keys)
	independent_prob = 63
	print(probs)
