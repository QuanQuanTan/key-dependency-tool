import numpy as np
from secrets import randbits
from collections import Counter
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import prep


XDDT8 = prep.getXDDT8()
YDDT8 = prep.getYDDT8()
DDT8 = prep.getDDT8()

probs_per_cell = [0 for _ in range(16)]
keys_reduction = [0 for _ in range(16)]


def key06c(RUN=True):
	# 0 (5,1,5,5) (5,b,4,5) (5,e,5,5) (6,1,5,5) 1
	# 0 (5,1,5,5) (6,5,5,5) 1
	# 6 (3,2,40,4) (4,6,4,1) 1
	# 6 (5,4,5,1),(5,0xb,4,5) (6,9,4,5) 1

	# c (4,2,4,5) (5,6,5,1) 1

	# c 0 (4,2,4,5) (3,7,40,4) (3,a,40,4) (4,8,0,0) (5,e,5,5) 0
	# 6 0 (5,e,5,5) (5,4,5,1) (5,1,5,5) (6,1,5,5) (6,9,4,5) [This must hold if the 2 linear constraints above hold]
	# 6 (5,e,5,5) (5,b,4,5) (5,4,5,1) (6,1,5,5) (6,5,5,5) (6,9,4,5) [This must hold if the 3 linear constraints above hold]
	left0 = [(5,1,5,5),(5,0xb,4,5),(5,0xe,5,5)]; right0 = (6,1,5,5)
	left1 = [(5,1,5,5)]; right1 = (6,5,5,5)
	left2 = [(3,2,0x40,4)]; right2 = (4,6,4,1)
	left3 = [(5,4,5,1),(5,0xb,4,5)]; right3 = (6,9,4,5)
	left4 = [(4,2,4,5)]; right4 = (5,6,5,1)

	possibleKeys0,_ = prep.computeLinearProb8(left0,right0)
	possibleKeys1,_ = prep.computeLinearProb8(left1,right1)
	possibleKeys0 = prep.overlapKeys(possibleKeys0,possibleKeys1)

	possibleKeys2,_ = prep.computeLinearProb8(left2,right2)
	possibleKeys3,_ = prep.computeLinearProb8(left3,right3,constant=[3,0])
	possibleKeys6 = prep.overlapKeys(possibleKeys2,possibleKeys3)	

	possibleKeysc,_ = prep.computeLinearProb8(left4,right4)
	keys_reduction[0x0] += -np.log2((256-len(possibleKeys0))/256)
	keys_reduction[0x0] += -np.log2((256-len(possibleKeys6))/256)
	keys_reduction[0x0] += -np.log2((256-len(possibleKeysc))/256)

	numTrails = 1 << 20
	numKeys = 1 << 6
	counts = []
	for _ in range(numKeys):
		count = 0
		k0 = prep.select(possibleKeys0); k0p = k0
		k6 = prep.select(possibleKeys6); k6p = k6
		kc = prep.select(possibleKeysc); kcp = kc

		for _ in range(numTrails):
			s51 = prep.select(YDDT8[0x05][0x05]); s51p = s51 ^ 0x05;
			s5b = prep.select(YDDT8[0x04][0x05]); s5bp = s5b ^ 0x05;
			s5e = prep.select(YDDT8[0x05][0x05]); s5ep = s5e ^ 0x05;
			s32 = prep.select(YDDT8[0x40][0x04]); s32p = s32 ^ 0x04;
			s54 = prep.select(YDDT8[0x05][0x01]); s54p = s54 ^ 0x01;
			s42 = prep.select(YDDT8[0x04][0x05]); s42p = s42 ^ 0x05;
			s37 = prep.select(YDDT8[0x40][0x04]); s37p = s37 ^ 0x04;
			s3a = prep.select(YDDT8[0x40][0x04]); s3ap = s3a ^ 0x04;

			s61,s61p = prep.Sbox8Pair(k0 ^ s51 ^ s5b ^ s5e,k0p ^ s51p ^ s5bp ^ s5ep) 
			s65,s65p = prep.Sbox8Pair(k0 ^ s51, k0p ^ s51p)
			s46,s46p = prep.Sbox8Pair(k6 ^ s32,k6p ^ s32p)
			s69,s69p = prep.Sbox8Pair(k6 ^ s54 ^ 0x1 ^ s5b,k6p ^ s54p ^ 0x1 ^ s5bp)
			s56,s56p = prep.Sbox8Pair(kc ^ s42,kcp ^ s42p)
			if s61 ^ s61p != 0x05: continue
			if s65 ^ s65p != 0x05: continue
			if s46 ^ s46p != 0x01: continue
			if s69 ^ s69p != 0x05: continue
			if s56 ^ s56p != 0x01: continue

			# nonlinear cases
			# c 0 (4,2,4,5) (3,7,40,4) (3,a,40,4) (4,8,0,0) (5,e,5,5) 0
			s48,s48p = prep.Sbox8Pair(k0 ^ s37 ^ s3a, k0p ^ s37p ^ s3ap)
			s5e_new,s5ep_new = prep.Sbox8Pair(kc ^ s42 ^ s48 ^ 0x2, kcp ^ s42p ^ s48p ^ 0x2)
			if s5e_new != s5e or s5ep_new != s5ep: continue

			# same round (ignored)
			# 6 0 (5,e,5,5) (5,4,5,1) (5,1,5,5) (6,1,5,5) (6,9,4,5)
			# 6 (5,e,5,5) (5,b,4,5) (5,4,5,1) (6,1,5,5) (6,5,5,5) (6,9,4,5)
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if (numKeys == len(counts)): pass
	else: keys_reduction[0] += -np.log2((numKeys-len(counts))/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_06c.png')
	plt.clf()
	# np.savetxt('k06c_counts.txt',counts)
	print(counts)



def key3f(RUN=True):
	# 3 f (1,4,8,10) (0,6,2,8) (0,9,2,8) (1,b,0,0) (2,9,10,40) 0
	numTrails = 1 << 16
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k3 = prep.select([i for i in range(256)]); k3p = k3
		kf = prep.select([i for i in range(256)]); kfp = kf
		for _ in range(numTrails):
			s14 = prep.select(YDDT8[0x08][0x10]); s14p = s14 ^ 0x10;
			s06 = prep.select(YDDT8[0x02][0x08]); s06p = s06 ^ 0x08;
			s09 = prep.select(YDDT8[0x02][0x08]); s09p = s09 ^ 0x08;
			s1b,s1bp = prep.Sbox8Pair(s06 ^ s09 ^ kf, s06p ^ s09p ^ kfp)
			s29,s29p = prep.Sbox8Pair(s1b ^ s14 ^ 0x3 ^ k3, s1bp ^ s14p ^ 0x3 ^ k3p)
			if s29 ^ s29p != 0x40: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if len(counts) == numKeys: keys_reduction[3] = 0
	else: keys_reduction[3] = -np.log2((numKeys-len(counts))/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_3f.png')
	plt.clf()
	np.savetxt('k3f_counts.txt',counts)



def key4():
	# 4 (5,2,5,5) (6,6,5,5) 1
	left = [(5,2,5,5)]
	right = (6,6,5,5)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[0x4] = 0
	keys_reduction[4] = -np.log2((256-possibleKeys)/256)
	probs_per_cell[4] = distribution
	print(keys_reduction[4])
	print(distribution)

def key5b(RUN=True):
	# b (6,3,5,5) (7,7,5,5) 	
	# 5 b (1,2,4,1) (2,6,0,0) (2,9,10,40) (3,b,40,4) 0
	possibleKeys,_ = prep.computeLinearProb8([(6,3,5,5)],(7,7,5,5))
	keys_reduction[5] = -np.log2((256-len(possibleKeys))/256)
	numTrails = 1 << 16
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k5 = prep.select([i for i in range(256)]); k5p = k5 ^ 1
		kb = prep.select(possibleKeys); kbp = kb
		for _ in range(numTrails):
			s63 = prep.select(YDDT8[0x05][0x05]); s63p = s63 ^ 0x05;
			s12 = prep.select(YDDT8[0x04][0x01]); s12p = s12 ^ 0x01;
			s29 = prep.select(YDDT8[0x10][0x40]); s29p = s29 ^ 0x40;
			s77,s77p = prep.Sbox8Pair(s63 ^ kb, s63p ^ kbp)
			s26,s26p = prep.Sbox8Pair(s12 ^ k5, s12p ^ k5p)
			s3b,s3bp = prep.Sbox8Pair(s26 ^ s29 ^ kb, s26p ^ s29p ^ kbp)
			if s77 ^ s77p != 0x05: continue
			if s3b ^ s3bp != 0x04: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if len(counts) == numKeys: pass
	else: keys_reduction[5] += -np.log2((numKeys-len(counts))/numKeys) 
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_5b.png')
	plt.clf()
	np.savetxt('k5b_counts.txt',counts)

def key8():
	# 8 (0,3,2,8) (1,7,8,10) 1
	left = [(0,3,2,8)]
	right = (1,7,8,0x10)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[0x8] = 0
	keys_reduction[8] = -np.log2((256-possibleKeys)/256)
	probs_per_cell[8] = distribution
	print(keys_reduction[8])
	print(distribution)

def key9(RUN=True):
	# 9 (2,3,10,40) (3,7,40,4) 1
	# 9 (6,1,5,5) (7,5,5,1) 1
	left0 = [(2,3,0x10,0x40)] 
	right0 = (3,7,0x40,4)
	left1 = [(6,1,5,5)]
	right1 = (7,5,5,1)
	possibleKeys0,distribution = prep.computeLinearProb8(left0,right0)
	possibleKeys1,distribution = prep.computeLinearProb8(left1,right1)
	# find the overlap of the keys
	possibleKeys = prep.overlapKeys(possibleKeys0,possibleKeys1)
	keys_reduction[9] += -np.log2((256-len(possibleKeys))/256)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k9 = prep.select(possibleKeys); k9p = k9
		for _ in range(numTrails):
			s23 = prep.select(YDDT8[0x10][0x40]); s23p = s23 ^ 0x40
			s61 = prep.select(YDDT8[0x05][0x05]); s61p = s61 ^ 0x05
			s37,s37p = prep.Sbox8Pair(s23 ^ k9, s23p ^ k9p)
			s75,s75p = prep.Sbox8Pair(s61 ^ k9, s61p ^ k9p)
			if s37 ^ s37p != 0x04: continue
			if s75 ^ s75p != 0x01: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[9] += -np.log2((numKeys-len(counts))/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_9.png')
	plt.clf()
	np.savetxt('k9_counts.txt',counts)


	


def keya():
	# a (4,0,4,5) (5,4,5,1) 1
	left = [(4,0,4,5) ]
	right = (5,4,5,1)
	possibleKeys,distribution = prep.computeLinearProb8(left,right,constant=[0xf])
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[0xa] = 0
	keys_reduction[0xa] = -np.log2((256-possibleKeys)/256)
	probs_per_cell[0xa] = distribution
	print(keys_reduction[0xa])
	print(distribution)

def keyd():
	# d (4,6,4,1) (4,9,5,5) (5,b,4,5) 1
	left = [(4,6,4,1),(4,9,5,5)]
	right = (5,0xb,4,5) 
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[0xd] = 0
	else: keys_reduction[0xd] = -np.log2((256-possibleKeys)/256)
	probs_per_cell[0xd] = distribution
	print(keys_reduction[0xd])
	print(distribution)

def keye():
	# e (0,0,2,8) (1,4,8,10) 1
	left = [(0,0,2,8)]
	right = (1,4,8,0x10)
	possibleKeys,distribution = prep.computeLinearProb8(left,right,constant=[0x1])
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[0xe] = 0
	else: keys_reduction[0xe] = -np.log2((256-possibleKeys)/256)
	probs_per_cell[0xe] = distribution
	print(keys_reduction[0xe])
	print(distribution)


if __name__ == "__main__":
	# key4()
	# key8()
	# keya()
	# keyd()
	# keye()
	# key3f()
	# key9()
	# key5b()
	key06c()