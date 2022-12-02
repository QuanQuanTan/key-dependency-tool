import numpy as np
from secrets import randbits
from collections import Counter
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, '../')
import prep



XDDT8 = prep.getXDDT8()
YDDT8 = prep.getYDDT8()
DDT8 = prep.getDDT8()

probs_per_cell = [0 for _ in range(16)]
keys_reduction = [0 for _ in range(16)]


def key04bc(RUN=False): # CPP
	index = 0
	if RUN == False:
		left1 = [(2,1,1,0x20)]; right1 = (3,5,0x20,0x80)
		left3 = [(5,2,3,0x20)]; right3 = (6,6,0x20,0x80)
		possibleKeys4,_ = prep.computeLinearProb8(left1,right1)
		possibleKeysc,_ = prep.computeLinearProb8(left3,right3)	

		keys_reduction[index] += -np.log2(len(possibleKeys4)/256)
		keys_reduction[index] += -np.log2(len(possibleKeysc)/256)
		probs_per_cell[index] = {14.482261618655691: 1.0} 
		keys_reduction[index] += 1.4902249956730629
		return 

def key1(RUN=False):
	index = 1
	left0 = [(4,3,0x90,0x2)]; right0 = (5,7,2,9)
	left1 = [(4,3,0x90,0x2),(4,9,5,1)]; right1 = (5,0xf,3,0x20)
	tmp1,_ = prep.computeLinearProb8(left0,right0)
	tmp2,_ = prep.computeLinearProb8(left1,right1)
	# find the overlap of the keys
	pk1 = prep.overlapKeys(tmp1,tmp2)
	if len(pk1) < 256: keys_reduction[index] = -np.log2(pk1/256)
	if RUN == False:
		probs_per_cell[index] = {4.19547937356867: 0.4833984375, 5.847296237503679: 0.0322265625, 6.438553916497797: 0.109375, 6.99413411148911: 0.2099609375, 7.978682523152651: 0.1650390625}
		keys_reduction[index] += 0
		return 
	# 1 (4,3,90,2) (5,7,2,9) 1
	# 1 (4,3,90,2) (4,9,5,1) (5,f,3,20) 1

	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k1 = prep.select(pk1); k1p = k1
		for _ in range(numTrails):
			s43 = prep.select(YDDT8[0x90][0x02]); s43p = s43 ^ 0x02
			s49 = prep.select(YDDT8[0x05][0x01]); s49p = s49 ^ 0x01
			s57,s57p = prep.Sbox8Pair(s43 ^ k1, s43p ^ k1p)
			s5f,s5fp = prep.Sbox8Pair(s43 ^ s49 ^ k1, s43p ^ s49p ^ k1p)
			if s57 ^ s57p != 0x09: continue
			if s5f ^ s5fp != 0x20: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_1.png')
	plt.clf()
	np.savetxt('k1_counts.txt',counts)

def key2(RUN=False):
	index = 2
	left0 = [(4,1,5,1)] 
	right0 = (5,5,1,0x20)
	pk2,_ = prep.computeLinearProb8(left0,right0)
	keys_reduction[index] += -np.log2(len(pk2)/256)
	if RUN ==False:	
		probs_per_cell[index] = {2.0001741773162975: 0.26639344262295084, 2.9988188211946345: 0.7336065573770492}
		keys_reduction[index] += 1.0692626624371138
		return
	# 2 (4,1,5,1) (5,5,1,20) 1
	# 2n (6,0,20,80) (7,4,80,3) 1 (6,0 constant = B) 
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k2_tk1 = randbits(8); k2p_tk1 = k2_tk1
		k2_tk2 = prep.select(pk2) ^ k2_tk1; k2p_tk2 = k2_tk2 
		
		k2_tk2_next = prep.LFSR8(k2_tk2,2); k2p_tk2_next = prep.LFSR8(k2p_tk2,2)
		for _ in range(numTrails):
			s41 = prep.select(YDDT8[0x05][0x01]); s41p = s41 ^ 0x01
			s60 = prep.select(YDDT8[0x20][0x80]); s60p = s60 ^ 0x80
			s55,s55p = prep.Sbox8Pair(s41 ^ k2_tk1 ^ k2_tk2, s41p ^ k2p_tk1 ^ k2p_tk2)
			if s55 ^ s55p != 0x20: continue
			s74,s74p = prep.Sbox8Pair(s60 ^ k2_tk1 ^ k2_tk2_next ^ 0xB, s60p ^ k2p_tk1 ^ k2p_tk2_next ^ 0xB)
			if s74 ^ s74p != 0x03: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_2.png')
	plt.clf()
	np.savetxt('k2_counts.txt',counts)

def key3e(RUN=False):
	index = 3
	if RUN ==False:
		probs_per_cell[index] = {0.999671806839783: 0.31108144192256343, 1.9985215388179751: 0.6889185580774366}
		keys_reduction[index] += 0.45117809154124894
		return
	# 5 e (6,6,20,80) (5,4,3,20) (5,b,1,20) (6,9,0,0) (7,b,80,3) 0 (5,4 constant = 3)
	numTrails = 1 << 16
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k5 = randbits(8); k5p = k5
		ke = randbits(8); kep = ke
		for _ in range(numTrails):
			s66 = prep.select(YDDT8[0x20][0x80]); s66p = s66 ^ 0x80;
			s54 = prep.select(YDDT8[0x03][0x20]); s54p = s54 ^ 0x20;
			s5b = prep.select(YDDT8[0x01][0x20]); s5bp = s5b ^ 0x20;
			s69,s69p = prep.Sbox8Pair(s54 ^ 0x3 ^ s5b ^ ke, s54p ^ 0x3 ^ s5bp ^ kep)
			s7b,s7bp = prep.Sbox8Pair(s66 ^ s69 ^ k5, s66p ^ s69p ^ k5p)
			if s7b ^ s7bp != 0x03: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if len(counts) == numKeys: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_3e.png')
	plt.clf()
	np.savetxt('k3e_counts.txt',counts)

def key8():
	index = 8
	# 8 (5,1,1,20) (6,5,20,80) 1 
	left = [(5,1,1,0x20)]
	right = (6,5,0x20,0x80)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	if possibleKeys == 256: keys_reduction[index] = 0
	keys_reduction[index] = -np.log2(len(possibleKeys)/256)
	probs_per_cell[index] = distribution


def keyd(RUN=False):
	index = 0xd
	left0 = [(1,2,4,5)]; right0 = (2,6,5,5)
	left1 = [(1,2,4,5),(1,8,0x40,4)]; right1 = (2,0xe,1,0x20)
	tmp1,_ = prep.computeLinearProb8(left0,right0)
	tmp2,_ = prep.computeLinearProb8(left1,right1)
	# find the overlap of the keys
	pkd = prep.overlapKeys(tmp1,tmp2)
	keys_reduction[index] = -np.log2(len(pkd)/256)
	if RUN ==False:
		probs_per_cell[index] = {4.000456782640217: 1.0}
		keys_reduction[index] += 0
		return
	# d (1,2,4,5) (2,6,5,5) 1
	# d (1,2,4,5) (1,8,40,4) (2,e,1,20) 1 (1,8 constant = 2)

	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		kd = prep.select(possibleKeys); kdp = kd
		for _ in range(numTrails):
			s12 = prep.select(YDDT8[0x04][0x05]); s12p = s12 ^ 0x05
			s18 = prep.select(YDDT8[0x40][0x04]); s18p = s18 ^ 0x04
			s26,s26p = prep.Sbox8Pair(s12 ^ kd, s12p ^ kdp)
			s2e,s2ep = prep.Sbox8Pair(s12 ^ s18 ^ 0x2 ^ kd, s12p ^ s18p ^ 0x2 ^ kdp)
			if s26 ^ s26p != 0x05: continue
			if s2e ^ s2ep != 0x20: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[d] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_d.png')
	plt.clf()
	np.savetxt('kd_counts.txt',counts)

def key9():
	index = 9
	# 9 (7,1,93,ea) (7,b,80,3) (8,d,3,bc) 1 
	left = [(7,1,0x93,0xea),(7,0xb,0x80,3)]
	right = (8,0xd,3,0xbc)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(len(possibleKeys)/256)
	probs_per_cell[index] = distribution

def keyf():
	index = 0xf
	# f (5,3,3,20) (6,7,20,80) 1
	left = [(5,3,3,0x20)]
	right = (6,7,0x20,0x80)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	if possibleKeys == 256: keys_reduction[index] = 0
	keys_reduction[index] = -np.log2(len(possibleKeys)/256)
	probs_per_cell[index] = distribution





if __name__ == "__main__":
	key04bc()
	key1()
	key2()
	key3e()
	key8()
	keyd()
	key9()
	keyf()
	print(keys_reduction)
	print(probs_per_cell)
	print('sum keyreduction:', sum(keys_reduction))
	print(prep.combineDictionary(probs_per_cell))