import numpy as np
from secrets import randbits
from collections import Counter
import matplotlib.pyplot as plt
import sys
sys.path.insert(1, '../')
import prep


XDDT4 = prep.getXDDT4()
YDDT4 = prep.getYDDT4()
DDT4 = prep.getDDT4()

probs_per_cell = [0 for _ in range(16)]
keys_reduction = [0 for _ in range(16)]





def key2(RUN=False):
	index = 0x2
	if RUN == False:
		left0 = [(0,1,0x1,0x8)]; right0 = (1,5,0x8,0x4)
		possibleKeys2, _ = prep.computeLinearProb4(left0,right0)
		if len(possibleKeys2) < 16: keys_reduction[index] += -np.log2(len(possibleKeys2)/16)
		keys_reduction[index] += 1.9019679170394732
		probs_per_cell[index] = {0.0: 1.0}
		return
	# 2 (2,0,5,2) (3,4,2,5) 1 *
	# 2 (0,1,1,8) (1,5,8,4) 1
	left0 = [(0,1,0x1,0x8)]; right0 = (1,5,0x8,0x4)
	possibleKeys2, _ = prep.computeLinearProb4(left0,right0)
	
	if len(possibleKeys2) < 16: keys_reduction[index] = -np.log2(len(possibleKeys2)/16)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k2_0 = randbits(4); k2p_0 = k2_0
		k2_1 = prep.select(possibleKeys2) ^ k2_0; k2p_1 = k2_1 
		for _ in range(numTrails):
			s01 = prep.select(YDDT4[0x1][0x8]); s01p = s01 ^ 0x8
			s20 = prep.select(YDDT4[0x5][0x2]); s20p = s20 ^ 0x2
			
			s15,s15p = prep.Sbox4Pair(s01 ^ k2_0 ^ k2_1, s01p ^ k2p_0 ^ k2p_1)
			k2_1new = prep.LFSR4(k2_1,2); k2p_1new = prep.LFSR4(k2p_1,2)
			s34,s34p = prep.Sbox4Pair(s20 ^ 7 ^ k2_0 ^ k2_1new, s20p ^ 7 ^ k2p_0 ^ k2p_1new); 

			if s15 ^ s15p != 0x4: continue
			if s34 ^ s34p != 0x5: continue

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



def key4():
	index = 0x4
	# 4 (0,0,5,2) (1,4,2,5) 1
	left = [(0,0,5,0x2)]
	right = (1,4,0x2,5)
	possibleKeys,distribution = prep.computeLinearProb4(left,right,constant=[1])
	possibleKeys = len(possibleKeys)
	if possibleKeys == 16: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/16)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)


def key5e(RUN=False):
	index = 0x5
	if RUN == False:
		keys_reduction[index] += 0.18461670418646142
		probs_per_cell[index] = {1.0005689931098047: 0.13762486126526083, 1.4128647272344896: 0.3018867924528302, 1.9997807381686585: 0.2885682574916759, 2.998792658597951: 0.27192008879023305}
		return
	# 5 e (2,6,5,a) (1,4,2,5) (1,b,9,5) (2,9,0,0) (3,b,a,5) 0 * 
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k5 = randbits(4); k5p = k5
		ke = randbits(4); kep = ke
		
		for _ in range(numTrails):
			s26 = prep.select(YDDT4[0x5][0xa]); s26p = s26 ^ 0xa
			s14 = prep.select(YDDT4[0x2][0x5]); s14p = s14 ^ 0x5
			s1b = prep.select(YDDT4[0x9][0x5]); s1bp = s1b ^ 0x5
			
			s29,s29p = prep.Sbox4Pair(s14 ^ 0 ^ s1b ^ ke, s14p ^ 0 ^ s1bp ^ kep)
			s3b,s3bp = prep.Sbox4Pair(s29 ^ s26 ^ k5, s29p ^ s26p ^ k5p)

			if s3b ^ s3bp != 0x5: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_5.png')
	plt.clf()
	np.savetxt('k5_counts.txt',counts)



def key8():
	index = 0x8
	# 8 (1,1,8,5) (2,5,5,2) 1
	left = [(1,1,0x8,5)]
	right = (2,5,5,2)
	possibleKeys,distribution = prep.computeLinearProb4(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 16: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/16)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)


def key9():
	index = 0x9
	# 9 (3,1,a,5) (3,b,a,5) (4,d,5,a) 1
	left = [(3,1,0xa,5),(3,0xb,0xa,5)]
	right = (4,0xd,5,0xa)
	possibleKeys,distribution = prep.computeLinearProb4(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 16: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/16)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)





def keybc(RUN=False):
	index = 0xb
	if RUN == False:
		left = [(1,2,2,5)]; right = (2,6,5,0xa)
		possibleKeysc, _ = prep.computeLinearProb4(left,right)
		if len(possibleKeysc) < 16: keys_reduction[index] += -np.log2(len(possibleKeysc)/16)
		keys_reduction[index] += 1.8707169830550336
		probs_per_cell[index] = {2.0004148042969527: 1.0}
		return
	# c (1,2,2,5) (2,6,5,a) 1
	# b c (1,5,8,4) (1,2,2,5) (2,a,4,2) (2,e,5,a)
	# b (1,5,8,4) (2,6,5,a) (2,a,4,2) (2,e,5,a)
	left = [(1,2,2,5)]; right = (2,6,5,0xa)

	possibleKeysc, _ = prep.computeLinearProb4(left,right)
	if len(possibleKeysc) < 16: keys_reduction[index] = -np.log2(len(possibleKeysc)/16)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		kc = prep.select(possibleKeysc); kcp = kc
		kb = randbits(4); kbp = kb
		for _ in range(numTrails):
			s12 = prep.select(YDDT4[0x2][0x5]); s12p = s12 ^ 0x5
			s15 = prep.select(YDDT4[0x8][0x4]); s15p = s15 ^ 0x5
			s18 = randbits(4); s18p = s18
				
			s26,s26p = prep.Sbox4Pair(s12 ^ kc ,s12p ^ kcp)
			s2a,s2ap = prep.Sbox4Pair(s15 ^ s18 ^ kb, s15p ^ s18p ^ kbp)
			s2e,s2ep = prep.Sbox4Pair(s12 ^ s18 ^ kc, s12p ^ s18p ^ kcp)

			if s26 ^ s26p != 0xa: continue
			if s2a ^ s2ap != 0x2: continue
			if s2e ^ s2ep != 0xa: continue

			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_bc.png')
	plt.clf()
	np.savetxt('kbc_counts.txt',counts)



def keyf():
	index = 0xf
	# f (1,3,9,5) (2,7,5,2) 1
	left = [(1,3,0x9,0x5)]
	right = (2,7,0x5,0x2)
	possibleKeys,distribution = prep.computeLinearProb4(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 16: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/16)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)


if __name__ == "__main__":
	key2()
	key4()
	key5e()
	key8()
	key9()
	keybc()
	keyf()
	print(keys_reduction)
	print(probs_per_cell)
	print('sum keyreduction:', sum(keys_reduction))
	print(prep.combineDictionary(probs_per_cell))
