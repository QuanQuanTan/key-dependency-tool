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

def key045cef(): # CPP
	index = 0
	# 4 (0,0,10,40) (1,4,40,6) 1 (0,0 constant=1)
	# c (1,2,4,6) (2,6,6,21) 1
	# c (1,2,4,6) (1,8,5,26) (2,e,20,90) 1 (1,8 constant=2)
	# f (1,3,40,6) (2,7,6,21) 1
	
	# f 5 (1,3,40,6) (0,4,10,40) (0,b,10,40) (1,9,0,0) (2,f,6,21) 0 (0,4 constant)
	# 5n e (2,6,6,21) (1,4,40,6) (1,b,4,6) (2,9,0,0) (3,b,21,20) 0 (1,4 constant)
	# c 4n (1,2,4,6) (1,8,5,26) (1,f,21,20) (2,2,0,0) (2,8,26,21) (3,e,21,20) 0 (1,8 constant, 2,8 constant)
	# 0 (0,7,4,5) (1,4,40,6) (1,8,5,26) (1,c,40,6)
	left0 = [(0,0,0x10,0x40)]; right0 = (1,4,0x40,6) 
	left1 = [(1,2,4,6)]; right1 = (2,6,6,0x21)
	left2 = [(1,2,4,6),(1,8,5,0x26)]; right2 = (2,0xe,0x20,0x90)
	left3 = [(1,3,0x40,6)]; right3 = (2,7,6,0x21)

	pk4, _ = prep.computeLinearProb8(left0,right0,constant=[1])
	tmp1, _ = prep.computeLinearProb8(left1,right1)
	tmp2, _ = prep.computeLinearProb8(left2,right2,constant=[0,2])
	pkf, _ = prep.computeLinearProb8(left3,right3)
	pkc = prep.overlapKeys(tmp1,tmp2)

	if len(pk4) < 256: keys_reduction[index] += -np.log2(len(pk4)/256)
	if len(pkf) < 256: keys_reduction[index] += -np.log2(len(pkf)/256)
	if len(pkc) < 256: keys_reduction[index] += -np.log2(len(pkc)/256)

	probs_per_cell[index] = {17.022385: 0.01050420168067227, 17.563344642857142: 0.029411764705882353, 18.07627215189874: 0.04149159663865546, 18.636602222222226: 0.11817226890756302, 19.148678421052633: 0.09978991596638656, 19.725944782608696: 0.1207983193277311, 20.113318749999998: 0.08403361344537816, 20.39948818181818: 0.05777310924369748, 20.739246111111107: 0.09453781512605042, 21.064258: 0.052521008403361345, 21.5020288: 0.13130252100840337, 22.809225986842105: 0.15966386554621848}
	keys_reduction[index] += 0.6901447374132128
	return


def key1(RUN=False):
	index = 1
	left0 = [(0x0,0x3,0x25,0x25),(0x0,0x9,0x50,0x4),(0x0,0xc,0x15,0x61)]; right0 = (0x1,0x3,0x40,0x06)
	left1 = [(0x0,0x3,0x25,0x25)]; right1 = (0x1,0x7,0x25,0xa5)
	left2 = [(0x0,0x3,0x25,0x25),(0x0,0x9,0x50,0x04)]; right2 = (0x1,0xf,0x21,0x20)
	tmp0, _ = prep.computeLinearProb8(left0,right0)
	tmp1, _ = prep.computeLinearProb8(left1,right1)
	tmp2, _ = prep.computeLinearProb8(left2,right2)
	pk1 = prep.overlapKeys(tmp0,tmp1,tmp2)

	if len(pk1) < 256: keys_reduction[index] = -np.log2(len(pk1)/256)
	if RUN == False:
		probs_per_cell[index] = {6.004652184102134: 1.0}
		keys_reduction[index] += 0
		return
	# 1 (0,3,25,25) (0,9,50,4) (0,c,15,61) (1,3,40,6) 1 
	# 1 (0,3,25,25) (1,7,25,a5) 1
	# 1 (0,3,25,25) (0,9,50,4) (1,f,21,20) 1
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	
	for _ in range(numKeys):
		count = 0
		k1 = prep.select(pk1); k1p = k1
		for _ in range(numTrails):
			s03 = prep.select(YDDT8[0x25][0x25]); s03p = s03 ^ 0x25
			s09 = prep.select(YDDT8[0x50][0x04]); s09p = s09 ^ 0x04
			s0c = prep.select(YDDT8[0x15][0x61]); s0cp = s0c ^ 0x61
			
			s13,s13p = prep.Sbox8Pair(s03 ^ s09 ^ s0c ^ k1, s03p ^ s09p ^ s0cp ^ k1p)
			s17,s17p = prep.Sbox8Pair(s03 ^ k1, s03p ^ k1p)
			s1f,s1fp = prep.Sbox8Pair(s03 ^ s09 ^ k1, s03p ^ s09p ^ k1p)
			
			if s13 ^ s13p != 0x06: continue
			if s17 ^ s17p != 0xa5: continue
			if s1f ^ s1fp != 0x20: continue
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

def key2():
	index = 2
	# 2 (2,0,6,21) (3,4,21,20) 1 (2,0 constant=7)
	left = [(2,0,0x06,0x21)]
	right = (3,4,0x21,0x20)
	possibleKeys,distribution = prep.computeLinearProb8(left,right,constant=[7])
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

def key8():
	index = 8
	# 8 (1,1,40,6) (2,5,6,21) 1 
	left = [(1,1,0x40,6)]
	right = (2,5,6,0x21)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

def key9():
	index = 9
	# 9 (3,1,90,3) (3,b,21,20) (4,d,20,90) 1
	left = [(3,1,0x90,3),(3,0xb,0x21,0x20)]
	right = (4,0xd,0x20,0x90)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

if __name__ == "__main__":
	key045cef() #cpp
	key1()
	key2()
	key8()
	key9()
	print(keys_reduction)
	print(probs_per_cell)
	print('sum keyreduction:', sum(keys_reduction))
	print(prep.combineDictionary(probs_per_cell))