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

def key09(RUN=False):
	index = 0
	if RUN == False:
		probs_per_cell[index] = {2.997516811620603: 0.3381770145310436, 3.41267751295568: 0.28269484808454426, 4.998384201604011: 0.37912813738441214} 
		keys_reduction[index] += 0.43585051001426833
		return
	# 0 9 (2,6,20,93) (1,4,20,80) (1,b,b8,80) (2,9,0,0) (3,b,93,b0) 0 (1,4 constant = 0)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k0 = randbits(8); k0p = k0
		k9 = randbits(8); k9p = k9
		for _ in range(numTrails):
			s26 = prep.select(YDDT8[0x20][0x93]); s26p = s26 ^ 0x93
			s14 = prep.select(YDDT8[0x20][0x80]); s14p = s14 ^ 0x80
			s1b = prep.select(YDDT8[0xb8][0x80]); s1bp = s1b ^ 0x80

			s29,s29p = prep.Sbox8Pair(s14 ^ 0 ^ s1b ^ k9, s14p ^ 0 ^ s1bp ^ k9p)
			s3b,s3bp = prep.Sbox8Pair(s26 ^ s29 ^ k0, s26p ^ s29p ^ k0p)

			if s3b ^ s3bp != 0xb0: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_09.png')
	plt.clf()
	np.savetxt('k09_counts.txt',counts)


def key3(RUN=False):
	index = 3
	left = [(0,1,0x21,0x20)]; right = (1,5,0x20,0x90)
	pk3, _ = prep.computeLinearProb8(left,right)
	if len(pk3) < 256: keys_reduction[index] = -np.log2(len(pk3)/256)
	if RUN == False:
		probs_per_cell[index] = {3.8286049575679852: 0.05813953488372093, 4.085276582313488: 0.1686046511627907, 4.414998930152978: 0.2189922480620155, 4.829498541637042: 0.3430232558139535, 5.422259353908085: 0.16472868217054262, 5.959937339543525: 0.046511627906976744}
		keys_reduction[index] = 0.9887727445767459
		return
	# 3 (0,1,21,20) (1,5,20,90) 1
	# 3n (2,0,90,3) (3,4,3,b0) 1 (2,0 constant = 7)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	
	for _ in range(numKeys):
		count = 0
		k3_tk1 = randbits(8); k3p_tk1 = k3_tk1
		k3_tk2 = prep.select(pk3) ^ k3_tk1; k3p_tk2 = k3_tk2
		k3_tk2_next = prep.LFSR8(k3_tk2,2); k3p_tk2_next = prep.LFSR8(k3p_tk2,2); 
		for _ in range(numTrails):
			s01 = prep.select(YDDT8[0x21][0x20]); s01p = s01 ^ 0x20
			s20 = prep.select(YDDT8[0x90][0x03]); s20p = s20 ^ 0x03
			
			s15,s15p = prep.Sbox8Pair(s01 ^ k3_tk1 ^ k3_tk2, s01p ^ k3p_tk1 ^ k3p_tk2)
			s34,s34p = prep.Sbox8Pair(s20 ^ 7 ^ k3_tk1 ^ k3_tk2_next, s20p ^ 7 ^ k3p_tk1 ^ k3p_tk2_next)

			if s15 ^ s15p != 0x90: continue
			if s34 ^ s34p != 0xb0: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_3.png')
	plt.clf()
	np.savetxt('k3_counts.txt',counts)


def key6(RUN=False):
	index = 6
	left0 = [(0,3,0x22,0x99),(0,9,0x01,0xb8),(0,0xc,0x04,0x01)]; right0 = (1,3,0x20,0x90)
	left1 = [(0,3,0x22,0x99)]; right1 = (1,7,0x99,0x13)
	left2 = [(0,3,0x22,0x99),(0,9,0x01,0xb8)]; right2 = (1,0xf,0x21,0x20)
	possibleKeys0, _ = prep.computeLinearProb8(left0,right0)
	possibleKeys1, _ = prep.computeLinearProb8(left1,right1)
	possibleKeys2, _ = prep.computeLinearProb8(left2,right2)

	pk6 = prep.overlapKeys(possibleKeys0,possibleKeys1,possibleKeys2)
	if len(pk6) < 256: keys_reduction[index] = -np.log2(len(pk6)/256)
	if RUN == False:
		probs_per_cell[index] = {3.6789364275004615: 1.0} 
		keys_reduction[index] += 0
		return
	# 6 (0,3,22,99) (0,9,1,b8) (0,c,4,1) (1,3,20,90) 1
	# 6 (0,3,22,99) (1,7,99,13) 1
	# 6 (0,3,22,99) (0,9,1,b8) (1,f,21,20) 1
	# 6 (0,c,4,1) (0,3,22,99) (1,3,20,90) (1,7,99,13) (1,f,21,20) # this is the sum of the first 3

	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	
	for _ in range(numKeys):
		count = 0
		k6 = prep.select(pk6); k6p = k6
		for _ in range(numTrails):
			s03 = prep.select(YDDT8[0x22][0x99]); s03p = s03 ^ 0x99
			s09 = prep.select(YDDT8[0x01][0xb8]); s09p = s09 ^ 0xb8
			s0c = prep.select(YDDT8[0x04][0x01]); s0cp = s0c ^ 0x01
			
			s13,s13p = prep.Sbox8Pair(s03 ^ s09 ^ s0c ^ k6, s03p ^ s09p ^ s0cp ^ k6p)
			s17,s17p = prep.Sbox8Pair(s03 ^ k6, s03p ^ k6p)
			s1f,s1fp = prep.Sbox8Pair(s03 ^ s09 ^ k6, s03p ^ s09p ^ k6p)
			
			if s13 ^ s13p != 0x90: continue
			if s17 ^ s17p != 0x13: continue
			if s1f ^ s1fp != 0x20: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_6.png')
	plt.clf()
	np.savetxt('k6_counts.txt',counts)




def key7():
	index = 7
	# 7 (0,0,21,20) (1,4,20,80) 1 (0,0 constant = 1)
	left = [(0,0,0x21,0x20)]
	right = (1,4,0x20,0x80)
	possibleKeys,distribution = prep.computeLinearProb8(left,right,constant=[1])
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)

def keyaf(RUN=False):
	index = 0xa
	left = [(1,2,0x01,0x20)]; right = (2,6,0x20,0x93)
	pkf, _ = prep.computeLinearProb8(left,right)
	if len(pkf) < 256: keys_reduction[index] = -np.log2(len(pkf)/256)
	if RUN == False:
		probs_per_cell[index] = {5.990559695999591: 0.4970703125, 7.590451577490543: 0.5029296875}
		keys_reduction[index] += 0
		return
	# f (1,2,1,20) (2,6,20,93) 1
	# a f (1,5,20,90) (1,2,1,20) (2,a,90,3) (2,e,20,80) 
	# a (1,5,20,90) (2,6,20,93) (2,a,90,3) (2,e,20,80) # covered by the first 2

	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	
	for _ in range(numKeys):
		count = 0
		kf = prep.select(pkf); kfp = kf
		ka = randbits(8); kap = ka
		for _ in range(numTrails):
			s12 = prep.select(YDDT8[0x01][0x20]); s12p = s12 ^ 0x20
			s15 = prep.select(YDDT8[0x20][0x90]); s15p = s15 ^ 0x90
			s18 = randbits(8); s18p = s18

			s26,s26p = prep.Sbox8Pair(s12 ^ kf, s12p ^ kfp)
			s2a,s2ap = prep.Sbox8Pair(s15 ^ s18 ^ ka, s15p ^ s18p ^ kap)
			s2e,s2ep = prep.Sbox8Pair(s12 ^ s18 ^ kf, s12p ^ s18p ^ kfp)

			if s26 ^ s26p != 0x93: continue
			if s2a ^ s2ap != 0x03: continue
			if s2e ^ s2ep != 0x80: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_af.png')
	plt.clf()
	np.savetxt('kaf_counts.txt',counts)


def keyb():
	index = 0xb
	# b (b,3,1,20) (c,7,20,80) 1
	left = [(0xb,3,1,0x20)]
	right = (0xc,7,0x20,0x80)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)

def keyc():
	index = 0xc
	# c (1,3,20,90) (2,7,90,3) 1
	left = [(1,3,0x20,0x90)]
	right = (2,7,0x90,3)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)

def keyd():
	index = 0xd
	# d (1,1,20,80) (2,5,80,3) 1
	left = [(1,1,0x20,0x80)]
	right = (2,5,0x80,3)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)

def keye():
	index = 0xe
	# e (3,1,80,3) (3,b,93,b0) (4,d,b0,80) 1
	left = [(3,1,0x80,3),(3,0xb,0x93,0xb0)]
	right = (4,0xd,0xb0,0x80)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution
	print(keys_reduction[index])
	print(distribution)

if __name__ == "__main__":
	key09()
	key3()
	key6()
	key7()
	keyaf()
	keyb()
	keyc()
	keyd()
	keye()
	print(keys_reduction)
	print(probs_per_cell)
	print('sum keyreduction:', sum(keys_reduction))
	print(prep.combineDictionary(probs_per_cell))





