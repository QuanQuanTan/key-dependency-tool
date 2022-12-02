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

def key0d(RUN=False):
	index = 0
	left = [(1,1,4,5)]; right = (2,5,2,8)
	pk0, _ = prep.computeLinearProb8(left,right)
	if len(pk0) < 256: keys_reduction[index] += -np.log2(len(pk0)/256)
	# 0 (1,1,4,5) (2,5,2,8) 1
	# 0 d (1,1,4,5) (0,6,20,80) (0,9,20,80) (1,b,0,0) (2,d,2,c) 0
	if RUN == False:
		probs_per_cell[index] = {3.4187339682726385: 0.15612903225806452, 4.402806485274654: 0.33161290322580644, 4.996032680838266: 0.1935483870967742, 5.999243099150504: 0.31870967741935485} 
		keys_reduction[index] += 0.4019474998384001
		return
	numTrails = 1 << 14
	numKeys = 1 << 10
	counts = []
	
	for _ in range(numKeys):
		count = 0
		k0 = prep.select(pk0); k0p = k0 ^ 0x7
		kd = randbits(8); kdp = kd
		for _ in range(numTrails):
			s11 = prep.select(YDDT8[0x04][0x05]); s11p = s11 ^ 0x05
			s06 = prep.select(YDDT8[0x20][0x80]); s06p = s06 ^ 0x80
			s09 = prep.select(YDDT8[0x20][0x80]); s09p = s09 ^ 0x80

			s25,s25p = prep.Sbox8Pair(s11 ^ k0, s11p ^ k0p)
			s1b,s1bp = prep.Sbox8Pair(s06 ^ s09 ^ kd, s06p ^ s09p ^ kdp)
			s2d,s2dp = prep.Sbox8Pair(s11 ^ s1b ^ k0, s11p ^ s1bp ^ k0p)

			if s25 ^ s25p != 0x08: continue
			if s2d ^ s2dp != 0x0c: continue
			count += 1
		if count > 0:
			counts.append(-np.log2(count/numTrails))
	if numKeys == len(counts): pass
	else: keys_reduction[index] += -np.log2(len(counts)/numKeys)
	counts.sort()
	plt.plot(counts)
	plt.savefig('experiment_0d.png')
	plt.clf()
	np.savetxt('k0d_counts.txt',counts)


def key28(RUN=False): # CPP
	index = 2
	# 8 (0,1,20,80) (0,b,20,80) (0,e,40,4) (1,1,4,5) 1
	# 8 (0,1,20,80) (1,5,80,2) 1
	# 2 (1,0,80,2) (2,4,2,8) 1 (1,0 constant)
	# 2 8 8n (1,0,80,2) (1,a,80,2) (0,1,20,80) (0,b,20,80) (1,d,0,0) (2,0,0,0) (2,a,2,8) (2,d,2,c) (3,0,4,6) 0 (1,0 constant = 3, 2,0 constant = 7)
	# 2 8 8n (1,0,80,2) (1,a,80,2) (0,1,20,80) (0,b,20,80) (1,d,0,0) (2,0,0,0) (2,a,2,8) (3,c,8,10) 0 (1,0 constant = 3, 2,0 constant = 7)
	
	left0 = [(0,1,0x20,0x80),(0,0xb,0x20,0x80),(0,0xe,0x40,4)]; right0 = (1,1,4,5)
	left1 = [(0,1,0x20,0x80)]; right1 = (1,5,0x80,2)
	possibleKeys0, _ = prep.computeLinearProb8(left0,right0)
	possibleKeys1, _ = prep.computeLinearProb8(left1,right1)
	possibleKeys8 = prep.overlapKeys(possibleKeys0,possibleKeys1)
	left2 = [(1,0,0x80,2)]; right2 =(2,4,2,8)
	possibleKeys2, _ = prep.computeLinearProb8(left2,right2,constant=[3])
	if len(possibleKeys8) < 256: keys_reduction[index] += -np.log2(len(possibleKeys8)/256)
	if len(possibleKeys2) < 256: keys_reduction[index] += -np.log2(len(possibleKeys2)/256)
	probs_per_cell[index] = {9.28750004484305: 0.0272216796875, 9.38039228125: 0.078125, 9.491038515625: 0.015625, 10.170349345703123: 0.25, 10.303799218749997: 0.015625, 11.592329687499998: 0.0078125, 12.215751383928572: 0.2734375, 13.001462304687498: 0.1875, 13.203228124999999: 0.0078125, 14.029278125: 0.0390625, 14.345009375: 0.00390625, 15.004022023809522: 0.08203125, 16.0056587628866: 0.0118408203125}
	keys_reduction[index] += 0
	return

def key5bc(RUN=False): # CPP
	index = 5
	# c (0,2,20,90) (1,6,90,2) 1
	# cn 5 (2,4,2,8) (1,6,90,2) (1,9,80,2) (2,b,0,0) (3,9,8,10) 0 (2,4 constant = 0)
	# b c (0,5,20,80) (0,2,20,90) (1,a,80,2) (1,e,90,2) 
	# b (0,5,20,80) (1,6,90,2) (1,a,80,2) (1,e,90,2) # a combination of 2 constraints above
	if RUN == False:
		left = [(0,2,0x20,0x90)]; right = (1,6,0x90,2)
		possibleKeysc, _ = prep.computeLinearProb8(left,right)
		if len(possibleKeysc) < 256: keys_reduction[index] = -np.log2(len(possibleKeysc)/256)
		probs_per_cell[index] = {5.537881666666667: 0.27692307692307694, 5.7429341379310355: 0.2230769230769231, 6.284784324324325: 0.2846153846153846, 6.66742392857143: 0.2153846153846154} 
		keys_reduction[index] += 1.0313332068047916
		return

def key7():
	index = 7
	# 7 (1,3,80,2) (2,7,2,8) 1 
	left = [(1,3,0x80,2)]
	right = (2,7,2,8)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

def keya():
	index = 0xa
	# a (2,2,2,8) (3,6,8,10) 1
	left = [(2,2,2,8)]
	right = (3,6,8,0x10)
	possibleKeys,distribution = prep.computeLinearProb8(left,right)
	possibleKeys = len(possibleKeys)
	if possibleKeys == 256: keys_reduction[index] = 0
	else: keys_reduction[index] = -np.log2(possibleKeys/256)
	probs_per_cell[index] = distribution

if __name__ == "__main__":
	key0d()
	key28()
	key5bc()
	key7()
	keya()
	print(keys_reduction)
	print(probs_per_cell)
	print('sum keyreduction:', sum(keys_reduction))
	print(prep.combineDictionary(probs_per_cell))




