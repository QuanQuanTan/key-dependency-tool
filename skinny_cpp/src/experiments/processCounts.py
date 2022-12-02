import numpy as np
import os
from collections import Counter
import matplotlib.pyplot as plt

directories = [
'exp_64_17L_I','exp_64_18L_I','exp_64_22L_BDM1','exp_64_22L_BDM2',
'exp_128_18L_I','exp_128_18L_II','exp_128_21L_I','exp_128_22L_I',
'exp_128_25L_I','exp_128_DDH9','exp_128_DDH10']

index = 0
print('directory: ', directories[index])
directory = './' + directories[index]


trails = {}

for filename in os.listdir(directory):
    if '.py' in filename:
        with open(directory + '/' + filename,'r') as f:
            lines = f.readlines()
        capture = False
        for i in range(len(lines)):
            if 'def key' in lines[i]:
                key = lines[i][4:-4]
                capture = True
            if capture == True:
                if 'numKeys' in lines[i]: 
                    trails[key] = int(lines[i].split(' ')[-1])
                    capture = False
                    
print(trails)
distributions = []
keySpaceReduction = []
Redo = []
for filename in os.listdir(directory):
    if 'counts.txt' not in filename: continue
    print(filename)
    counts = np.loadtxt(directory+'/'+filename)
    counts.sort()
    difference = [counts[i+1]-counts[i] for i in range(len(counts)-1)]
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(counts)
    plt.show()
    ok_flag = 'n'
    if len(counts) < 1024: keySpaceReduction.append(-np.log2(1024-len(counts))/1024)
    else: keySpaceReduction.append(0)
    while ok_flag == 'n':
        maxes = []
        peaks = int(input("enter number of peaks:"))
        if peaks == -1: 
            ok_flag = 'y'
            Redo.append(filename)
            del keySpaceReduction[-1]
            continue
        sorted_difference = sorted(difference)[-peaks:]
        for i in range(peaks):
            maxes.append(difference.index(sorted_difference[i]))
        maxes.sort()
        maxes.insert(0,0)
        maxes.append(len(counts))
        distribution = {}
        for i in range(len(maxes)-1):
            size = (maxes[i+1] - maxes[i])/len(counts)
            prob = np.average(counts[maxes[i]:maxes[i+1]])
            distribution[prob] = size
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(counts)
        ax.vlines(maxes, min(counts), max(counts), linestyles='dashed', colors='red')
        plt.show()
        ok_flag = input('is that ok?: ')
        if ok_flag == 'y': distributions.append(distribution)
print(directory)
print(trails)
for i in range(len(distributions)):
    print(distributions[i],keySpaceReduction[i])
print(Redo)