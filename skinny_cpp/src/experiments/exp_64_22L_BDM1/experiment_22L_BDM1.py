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
	index = 2
	if RUN == False:
		left0 = [(0,1,0x8,0x4)]; right0 = (1,5,0x4,0x2)
		possibleKeys2, _ = prep.computeLinearProb4(left0,right0)
		if len(possibleKeys2) < 16: keys_reduction[index] = -np.log2(len(possibleKeys2)/16)
		keys_reduction[index] += 2.123483053435
		probs_per_cell[index] = {1.001212418614345: 1.0}
		return
	# 2 (2,0,2,3) (3,4,3,c) 1 *
	# 2 (0,1,8,4) (1,5,4,2) 1

	left0 = [(0,1,0x8,0x4)]; right0 = (1,5,0x4,0x2)
	possibleKeys2, _ = prep.computeLinearProb4(left0,right0)
	if len(possibleKeys2) < 16: keys_reduction[index] = -np.log2(len(possibleKeys2)/16)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k2_0 = randbits(4); k2p_0 = k2_0
		k2_1 = randbits(4); k2p_1 = k2_1
		k2_2 = prep.select(possibleKeys2) ^ k2_0 ^ k2_1; k2p_2 = k2_2
		for _ in range(numTrails):
			s01 = prep.select(YDDT4[0x8][0x4]); s01p = s01 ^ 0x4
			s20 = prep.select(YDDT4[0x2][0x3]); s20p = s20 ^ 0x3
			
			s15,s15p = prep.Sbox4Pair(s01 ^ k2_0 ^ k2_1 ^ k2_2, s01p ^ k2p_0 ^ k2p_1 ^ k2p_2)
			k2_1new = prep.LFSR4(k2_1,2); k2p_1new = prep.LFSR4(k2p_1,2); 
			k2_2new = prep.LFSR4(k2_2,3); k2p_2new = prep.LFSR4(k2p_2,3); 
			s34,s34p = prep.Sbox4Pair(s20 ^ 7 ^ k2_0 ^ k2_1new ^ k2_2new, s20p ^ 7 ^ k2p_0 ^ k2p_1new ^ k2p_2new); 

			if s15 ^ s15p != 0x2: continue
			if s34 ^ s34p != 0xc: continue

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
	index = 4
	# 4 (0,0,9,4) (1,4,4,2) 1 *
	left = [(0,0,9,4)]
	right = (1,4,4,2)
	possibleKeys,distribution = prep.computeLinearProb4(left,right,constant=[1])
	possibleKeys = len(possibleKeys)
	if possibleKeys == 16: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/16)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)



def key5e(RUN=False):
	index = 5
	if RUN == False:
		keys_reduction[index] += 1.1607962119030562
		probs_per_cell[index] = {0.9992445889072284: 0.13755458515283842, 1.998008475615602: 0.8624454148471615}
		return
	# 5 e (2,6,2,3) (1,4,4,2) (1,b,4,2) (2,9,0,0) (3,b,3,c) 0 *
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k5 = randbits(4); k5p = k5
		ke = randbits(4); kep = ke
		
		for _ in range(numTrails):
			s26 = prep.select(YDDT4[0x2][0x3]); s26p = s26 ^ 0x3
			s14 = prep.select(YDDT4[0x4][0x2]); s14p = s14 ^ 0x2
			s1b = prep.select(YDDT4[0x4][0x2]); s1bp = s1b ^ 0x2
			
			s29,s29p = prep.Sbox4Pair(s14 ^ 0 ^ s1b ^ ke, s14p ^ 0 ^ s1bp ^ kep)
			s3b,s3bp = prep.Sbox4Pair(s29 ^ s26 ^ k5, s29p ^ s26p ^ k5p)

			if s3b ^ s3bp != 0xc: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_5b.png')
	plt.clf()
	np.savetxt('k5b_counts.txt',counts)

def key8():
	index = 8
	# 8 (1,1,4,2) (2,5,2,3) 1
	left = [(1,1,0x4,2)]
	right = (2,5,2,3)
	possibleKeys,distribution = prep.computeLinearProb4(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 16: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/16)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)


def key9():
	index = 9
	# 9 (3,1,6,4) (3,b,3,c) (4,d,c,2) 1
	left = [(3,1,6,4),(3,0xb,3,0xc)]
	right = (4,0xd,0xc,2)
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
		left0 = [(0,1,0x8,0x4)]; right0 = (1,5,0x4,0x2)
		possibleKeys2, _ = prep.computeLinearProb4(left0,right0)
		if len(possibleKeys2) < 16: keys_reduction[index] = -np.log2(len(possibleKeys2)/16)
		keys_reduction[index] += 1.0516327684153224
		probs_per_cell[index] = {3.003080778706044: 1.0} 
		return
	# c (1,2,4,2) (2,6,2,3) 1
	# b c (1,5,4,2) (1,2,4,2) (2,a,2,3) (2,e,2,6)
	# b (1,5,4,2) (2,6,2,3) (2,a,2,3) (2,e,2,6)


	left = [(1,2,4,2)]; right = (2,6,2,3)
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
			s12 = prep.select(YDDT4[0x4][0x2]); s12p = s12 ^ 0x2
			s15 = prep.select(YDDT4[0x4][0x2]); s15p = s15 ^ 0x2
			s18 = randbits(4); s18p = s18
				
			s26,s26p = prep.Sbox4Pair(s12 ^ kc, s12p ^ kcp)
			s2a,s2ap = prep.Sbox4Pair(s15 ^ s18 ^ kb, s15p ^ s18p ^ kbp)
			s2e,s2ep = prep.Sbox4Pair(s12 ^ s18 ^ kc, s12p ^ s18p ^ kcp)

			if s26 ^ s26p != 0x3: continue
			if s2a ^ s2ap != 0x3: continue
			if s2e ^ s2ep != 0x6: continue

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
	# f (1,3,4,2) (2,7,2,3) 1
	left = [(1,3,4,2)]
	right = (2,7,2,3)
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