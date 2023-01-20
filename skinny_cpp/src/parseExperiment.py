import numpy as np
import sys
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift




def main():
	file = sys.argv[1]
	outputfile = file[:file.find('.')]+'_dist.txt'
	counts = np.loadtxt(file,dtype=float)
	counts.sort()
	# N = len(counts)
	# T = 1/N
	# yplot = fftshift(counts)
	# xf = [i for i in range(len(counts))]
	# plt.plot(xf,1/N * np.abs(yplot))
	# plt.plot(counts)
	# plt.savefig('tmp.png')
	# plt.clf()

	# find the difference
	difference = [ counts[i+1] - counts[i] for i in range(len(counts)-1)]
	# we find difference that transit from positive to near 0 to positive again

	peaks, _ = find_peaks(difference, distance=64, prominence=0.015)
	plt.plot(difference)
	plt.plot(counts)
	plt.plot(peaks, counts[peaks], "x")
	plt.savefig(file[:file.find('.')]+'.png')
	plt.clf()

	peaks = peaks.tolist()
	peaks.insert(0,0)
	peaks.append(len(counts))
	with open(outputfile,'w') as f:
		for i in range(len(peaks)-1):
			print(np.average(counts[peaks[i]:peaks[i+1]]),(peaks[i+1]-peaks[i])/len(counts)*100, file=f)


if __name__ == "__main__":
	main()