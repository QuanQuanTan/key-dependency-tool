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

def key2(RUN=False):
	index = 0x2
	left0 = [(0,1,0x21,0x20)]; right0 = (1,5,0x20,0x80)
	pk2, _ = prep.computeLinearProb8(left0,right0)
	if len(pk2) < 256: keys_reduction[index] += -np.log2(len(pk2)/256)
	if RUN == False:
		probs_per_cell[index] = {2.839493197885927: 0.04734848484848485, 3.09115040268763: 0.20265151515151514, 3.4074832360405494: 0.11553030303030302, 3.6745075197156947: 0.1553030303030303, 3.8579897734267528: 0.2746212121212121, 4.416737770831407: 0.17424242424242425, 4.994192760635604: 0.030303030303030304}
		keys_reduction[index] += 0.9556058806415466
		return
	# 2 (0,1,21,20) (1,5,20,80) 1
	# 2n (2,0,90,3) (3,4,3,20) 1 (2,0 constant = 7)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k2_tk1 = randbits(8); k2p_tk1 = k2_tk1
		k2_tk2 = randbits(8); k2p_tk2 = k2_tk2
		k2_tk3 = prep.select(pk2) ^ k2_tk1 ^ k2_tk2; k2p_tk3 = k2_tk3

		k2_tk2_next = prep.LFSR8(k2_tk2,2); k2p_tk2_next = prep.LFSR8(k2p_tk2,2); 
		k2_tk3_next = prep.LFSR8(k2_tk3,3); k2p_tk3_next = prep.LFSR8(k2p_tk3,3); 

		for _ in range(numTrails):
			s01 = prep.select(YDDT8[0x21][0x20]); s01p = s01 ^ 0x20
			s20 = prep.select(YDDT8[0x90][0x03]); s20p = s20 ^ 0x03
			
			s15,s15p = prep.Sbox8Pair(s01 ^ k2_tk1 ^ k2_tk2 ^ k2_tk3, s01p ^ k2p_tk1 ^ k2p_tk2 ^ k2p_tk3)
			s34,s34p = prep.Sbox8Pair(s20 ^ 7 ^ k2_tk1 ^ k2_tk2_next ^ k2_tk3_next, s20p ^ 7 ^ k2p_tk1 ^ k2p_tk2_next ^ k2p_tk3_next)

			if s15 ^ s15p != 0x80: continue
			if s34 ^ s34p != 0x20: continue
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
	# 4 (0,0,1,20) (1,4,20,80) 1 (0,0 constant = 1)
	left = [(0,0,1,0x20)]
	right = (1,4,0x20,0x80)
	possibleKeys,distribution = prep.computeLinearProb8(left,right,constant=[1])
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

def key5e(RUN=False):
	index = 0x5
	if RUN == False:
		probs_per_cell[index] = {2.6740719516471674: 0.25390625, 2.7705922966104883: 0.087890625, 2.8158828991035625: 0.01953125, 2.872345713378304: 0.1318359375, 3.1958031965721574: 0.2783203125, 3.4166711753195433: 0.228515625}
		keys_reduction[index] += 0
		return
	# 5 e (2,6,80,3) (1,4,20,80) (1,b,20,80) (2,9,0,0) (3,b,3,20) 0 (1,4 constant = 0)

	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	
	for _ in range(numKeys):
		count = 0
		k5 = randbits(8); k5p = k5
		ke = randbits(8); kep = ke
		for _ in range(numTrails):
			s26 = prep.select(YDDT8[0x80][0x03]); s26p = s26 ^ 0x03
			s14 = prep.select(YDDT8[0x20][0x80]); s14p = s14 ^ 0x80
			s1b = prep.select(YDDT8[0x20][0x80]); s1bp = s1b ^ 0x80

			s29,s29p = prep.Sbox8Pair(s14 ^ 0 ^ s1b ^ ke, s14p ^ 0 ^ s1bp ^ kep)
			s3b,s3bp = prep.Sbox8Pair(s26 ^ s29 ^ k5, s26p ^ s29p ^ k5p)

			if s3b ^ s3bp != 0x20: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_5e.png')
	plt.clf()
	np.savetxt('k5e_counts.txt',counts)





def key8():
	index = 0x8
	# 8 (1,1,20,80) (2,5,80,43) 1
	left = [(1,1,0x20,0x80)]
	right = (2,5,0x80,0x43)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

def key9():
	index = 0x9
	# 9 (3,1,c0,7) (3,b,3,20) (4,d,20,90) 1
	left = [(3,1,0xc0,7),(3,0xb,3,0x20)]
	right = (4,0xd,0x20,0x90)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution


def keybc(RUN=False):
	index = 0xb
	left = [(1,2,0x20,0x80)]; right = (2,6,0x80,0x03)
	pkc, _ = prep.computeLinearProb8(left,right)
	if len(pkc) < 256: keys_reduction[index] += -np.log2(len(pkc)/256)
	if RUN == False:
		probs_per_cell[index] = {5.418494221970053: 1.0} 
		keys_reduction[index] += 1.028456446049228
		return
	# c (1,2,20,80) (2,6,80,3) 1
	# b c (1,5,20,80) (1,2,20,80) (2,a,80,3) (2,e,80,c0)
	# b (1,5,20,80) (2,6,80,3) (2,a,80,3) (2,e,80,c0) # by combining the top 2

	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	
	for _ in range(numKeys):
		count = 0
		kc = prep.select(pkc); kcp = kc
		kb = randbits(8); kbp = kb
		for _ in range(numTrails):
			s12 = prep.select(YDDT8[0x20][0x80]); s12p = s12 ^ 0x80
			s15 = prep.select(YDDT8[0x20][0x80]); s15p = s15 ^ 0x80
			s18 = randbits(8); s18p = s18

			s26,s26p = prep.Sbox8Pair(s12 ^ kc, s12p ^ kcp)
			s2a,s2ap = prep.Sbox8Pair(s15 ^ s18 ^ kb, s15p ^ s18p ^ kbp)
			s2e,s2ep = prep.Sbox8Pair(s12 ^ s18 ^ kc, s12p ^ s18p ^ kcp)

			if s26 ^ s26p != 0x03: continue
			if s2a ^ s2ap != 0x03: continue
			if s2e ^ s2ep != 0xc0: continue
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
	# f (1,3,20,80) (2,7,80,3) 1
	left = [(1,3,0x20,0x80)]
	right = (2,7,0x80,3)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

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