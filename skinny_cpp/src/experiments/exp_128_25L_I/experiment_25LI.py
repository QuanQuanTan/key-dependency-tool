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


def key1(RUN = False):
	index = 1
	# 1 (0,0,21,20) (1,4,20,90) 1 (0,0 constant = 1)
	left = [(0,0,0x21,0x20)]
	right = (1,4,0x20,0x90)
	possibleKeys,distribution = prep.computeLinearProb8(left,right,constant=[1])
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

def key28(RUN = False):
	index = 2
	if RUN == False:
		keys_reduction[index] += 0
		probs_per_cell[index] = {2.678925103893094: 0.244140625, 2.7769369369251846: 0.087890625, 2.862203212120417: 0.1513671875, 3.183370171983341: 0.2392578125, 3.412016061669547: 0.27734375}
		return
	# 2 8 (2,6,80,3) (1,4,20,90) (1,b,20,90) (2,9,0,0) (3,b,3,20) 0 (1,4 constant = 0)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	
	for _ in range(numKeys):
		count = 0
		k2 = randbits(8); k2p = k2
		k8 = randbits(8); k8p = k8
		for _ in range(numTrails):
			s26 = prep.select(YDDT8[0x80][0x03]); s26p = s26 ^ 0x03
			s14 = prep.select(YDDT8[0x20][0x90]); s14p = s14 ^ 0x90
			s1b = prep.select(YDDT8[0x20][0x90]); s1bp = s1b ^ 0x90

			s29,s29p = prep.Sbox8Pair(s14 ^ 0 ^ s1b ^ k8, s14p ^ 0 ^ s1bp ^ k8p)
			s3b,s3bp = prep.Sbox8Pair(s26 ^ s29 ^ k2, s26p ^ s29p ^ k2p)

			if s3b ^ s3bp != 0x20: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_28.png')
	plt.clf()
	np.savetxt('k28_counts.txt',counts)

def key45(RUN = False):
	index = 4
	left0 = [(0,6,0x80,0x02),(0,9,0xb1,0x22)]; right0 = (1,0xb,0x20,0x90)
	left1 = [(0,3,0x90,0x02),(0,9,0xb1,0x22)]; right1 = (1,0xf,0x20,0x80)
	left2 = [(0,3,0x90,0x02)]; right2 = (1,7,0x02,0x09)

	pk4, _ = prep.computeLinearProb8(left0,right0)
	tmp1, _ = prep.computeLinearProb8(left1,right1)
	tmp2, _ = prep.computeLinearProb8(left2,right2)
	pk5 = prep.overlapKeys(tmp1,tmp2)
	if len(pk4) < 256: keys_reduction[index] += -np.log2(len(pk4)/256)
	if len(pk5) < 256: keys_reduction[index] += -np.log2(len(pk5)/256)
	if RUN == False:
		keys_reduction[index] += 0
		probs_per_cell[index] = {5.201048564582327: 0.4921875, 8.007936050820067: 0.5078125}
		return
	# 4 (0,6,80,2) (0,9,b1,22) (1,b,20,90) 1 
	# 5 (0,3,90,2) (0,9,b1,22) (1,f,20,80) 1
	# 5 (0,3,90,2) (1,7,2,9) 1
	# 4 5 (0,9,b1,22) (0,6,80,2) (0,3,90,2) (1,7,2,9) (1,b,20,90) # combined from the constraints above
	# 4 (0,6,80,2) (1,7,2,9) (1,b,20,90) (1,f,20,80) # combined from the constraints above

	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k4 = prep.select(pk4); k4p = k4
		k5 = prep.select(pk5); k5p = k5

		for _ in range(numTrails):
			s06 = prep.select(YDDT8[0x80][0x02]); s06p = s06 ^ 0x02
			s09 = prep.select(YDDT8[0xb1][0x22]); s09p = s09 ^ 0x22
			s03 = prep.select(YDDT8[0x90][0x02]); s03p = s03 ^ 0x02
			
			s1b,s1bp = prep.Sbox8Pair(s06 ^ s09 ^ k4, s06p ^ s09p ^ k4p)
			s1f,s1fp = prep.Sbox8Pair(s03 ^ s09 ^ k5, s03p ^ s09p ^ k5p)
			s17,s17p = prep.Sbox8Pair(s03 ^ k5, s03p ^ k5p)

			if s1b ^ s1bp != 0x90: continue
			if s1f ^ s1fp != 0x80: continue
			if s17 ^ s17p != 0x09: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_4.png')
	plt.clf()
	np.savetxt('k4_counts.txt',counts)

def key7(RUN = False):
	index = 7
	left0 = [(0,1,0x21,0x20)]; right0 = (1,5,0x20,0x80)
	left1 = [(2,0,0x90,0x03)]; right1 = (3,4,0x03,0x20)
	tmp0, _ = prep.computeLinearProb8(left0,right0)
	tmp1, _ = prep.computeLinearProb8(left1,right1,constant=[7])
	pk7 = prep.overlapKeys(tmp0,tmp1)
	if len(pk7) < 256: keys_reduction[index] += -np.log2(len(pk7)/256)

	if RUN == False:
		keys_reduction[index] += 1.0255854101944728
		probs_per_cell[index] = {2.828039969983442: 0.04771371769383698, 3.0959469731862335: 0.22465208747514911, 3.4164681976705444: 0.11133200795228629, 3.764104376513577: 0.3359840954274354, 3.9721601324296905: 0.055666003976143144, 4.413106696689881: 0.1988071570576541, 4.999888345328399: 0.02584493041749503}
		return
	# 7 (0,1,21,20) (1,5,20,80) 1
	# 7n (2,0,90,3) (3,4,3,20) 1 (2,0 constant = 7)
	
	if len(pk7) < 256: keys_reduction[index] = -np.log2(len(pk7)/256)
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k7_tk1 = randbits(8); k7p_tk1 = k7_tk1
		k7_tk2 = randbits(8); k7p_tk2 = k7_tk2
		k7_tk3 = prep.select(pk7) ^ k7_tk1 ^ k7_tk2; k7p_tk3 = k7_tk3

		k7_tk2_next = prep.LFSR8(k7_tk2,2); k7p_tk2_next = prep.LFSR8(k7p_tk2,2)
		k7_tk3_next = prep.LFSR8(k7_tk3,3); k7p_tk3_next = prep.LFSR8(k7p_tk3,3)
		for _ in range(numTrails):
			s01 = prep.select(YDDT8[0x21][0x20]); s01p = s01 ^ 0x20
			s20 = prep.select(YDDT8[0x90][0x03]); s20p = s20 ^ 0x03

			s15,s15p = prep.Sbox8Pair(s01 ^ k7_tk1 ^ k7_tk2 ^ k7_tk3, s01p ^ k7p_tk1 ^ k7p_tk2 ^ k7p_tk3)			
			s34,s34p = prep.Sbox8Pair(s20 ^ 7 ^ k7_tk1 ^ k7_tk2_next ^ k7_tk3_next, s20p ^ 7 ^ k7p_tk1 ^ k7p_tk2_next ^ k7p_tk3_next)
			
			if s15 ^ s15p != 0x80: continue
			if s34 ^ s34p != 0x20: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_7.png')
	plt.clf()
	np.savetxt('k7_counts.txt',counts)


def key9c(RUN = False):
	index = 9
	# 9 (1,2,20,80) (2,6,80,3) 1
	# 9nnnnn (d,1,8,10) (e,5,10,40) 1
	# c 9 (1,5,20,80) (1,2,20,80) (2,a,80,3) (2,e,80,d0)
	# c (1,5,20,80) (2,6,80,3) (2,a,80,3) (2,e,80,d0) # combined from the constraints above 
	left0 = [(0x1,0x2,0x20,0x80)]; right0 = (0x2,0x06,0x80,0x03)
	pk9, _ = prep.computeLinearProb8(left0,right0)
	if len(pk9) < 256: keys_reduction[index] = -np.log2(len(pk9)/256)
	if RUN == False:
		keys_reduction[index] += 2.005646563141142
		probs_per_cell[index] = {5.692995949675065: 0.03529411764705882, 5.849204073082104: 0.20784313725490197, 6.397330488898613: 0.5137254901960784, 7.421745381393787: 0.24313725490196078}
		return
	
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	for _ in range(numKeys):
		count = 0
		k9_tk1 = randbits(8); k9p_tk1 = k9_tk1
		k9_tk2 = randbits(8); k9p_tk2 = k9_tk2
		k9_tk3 = prep.select(pk9) ^ k9_tk1 ^ k9_tk2; k9p_tk3 = k9_tk3

		k9_tk2_nnnnn = prep.LFSR8(prep.LFSR8(prep.LFSR8(prep.LFSR8(prep.LFSR8(k9_tk2,2),2),2),2),2)
		k9_tk3_nnnnn = prep.LFSR8(prep.LFSR8(prep.LFSR8(prep.LFSR8(prep.LFSR8(k9_tk3,3),3),3),3),3)

		k9p_tk2_nnnnn = prep.LFSR8(prep.LFSR8(prep.LFSR8(prep.LFSR8(prep.LFSR8(k9p_tk2,2),2),2),2),2)
		k9p_tk3_nnnnn = prep.LFSR8(prep.LFSR8(prep.LFSR8(prep.LFSR8(prep.LFSR8(k9p_tk3,3),3),3),3),3)

		kc_tk1 = randbits(8); kcp_tk1 = kc_tk1
		kc_tk2 = randbits(8); kcp_tk2 = kc_tk2
		kc_tk3 = randbits(8); kcp_tk3 = kc_tk3

		for _ in range(numTrails):
			s12 = prep.select(YDDT8[0x20][0x80]); s12p = s12 ^ 0x80
			sd1 = prep.select(YDDT8[0x08][0x10]); sd1p = sd1 ^ 0x10
			s15 = prep.select(YDDT8[0x20][0x80]); s15p = s15 ^ 0x80
			s18 = randbits(8); s18p = s18

			s26,s26p = prep.Sbox8Pair(s12 ^ k9_tk1 ^ k9_tk2 ^ k9_tk3, s12p ^ k9p_tk1 ^ k9p_tk2 ^ k9p_tk3)
			s2a,s2ap = prep.Sbox8Pair(s15 ^ s18 ^ kc_tk1 ^ kc_tk2 ^ kc_tk3, s15p ^ s18p ^ kcp_tk1 ^ kcp_tk2 ^ kcp_tk3)
			s2e,s2ep = prep.Sbox8Pair(s12 ^ s18 ^ k9_tk1 ^ k9_tk2 ^ k9_tk3, s12p ^ s18p ^ k9p_tk1 ^ k9p_tk2 ^ k9p_tk3)

			se5,se5p = prep.Sbox8Pair(sd1 ^ k9_tk1 ^ k9_tk2_nnnnn ^ k9_tk3_nnnnn, sd1p ^ k9p_tk1 ^ k9p_tk2_nnnnn ^ k9p_tk3_nnnnn)
			
			if s26 ^ s26p != 0x03: continue
			if se5 ^ se5p != 0x40: continue
			if s2a ^ s2ap != 0x03: continue
			if s2e ^ s2ep != 0xd0: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_9c.png')
	plt.clf()
	np.savetxt('k9c_counts.txt',counts)



def keyb(RUN = False):
	index = 0xb
	# b (1,1,20,90) (2,5,90,3) 1
	left = [(1,1,0x20,0x90)]
	right = (2,5,0x90,3)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

def keyd(RUN = False):
	index = 0xd
	# d (3,1,d0,7) (3,b,3,20) (4,d,20,80) 1
	left = [(3,1,0xd0,7),(3,0xb,3,0x20)]
	right = (2,7,0x80,3)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

def keye(RUN = False):
	index = 0xe
	# e (1,3,20,80) (2,7,80,3) 1
	left = [(1,3,0x20,0x80)]
	right = (2,7,0x80,3)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution



if __name__ == "__main__":
	key1()
	key28()
	key45()
	key7()
	key9c()
	keyb()
	keyd()
	keye()
	print(keys_reduction)
	print(probs_per_cell)
	print('sum keyreduction:', sum(keys_reduction))
	print(prep.combineDictionary(probs_per_cell))





