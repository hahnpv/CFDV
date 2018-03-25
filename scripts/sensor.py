#!/usr/bin/python
#shamelessly stolen from stack overflow

import numpy as np
import matplotlib.pyplot as plt
import glob

path = "*.sensor"
labels = list()

# load first line of first sensor file to divine headers and number of columns
headers=open(glob.glob(path)[0]).readline().rstrip().split(" ")
numsubplots = len(headers)-1;
fig, ax = plt.subplots(1, numsubplots)

# load all sensor files and plot data
for filename in glob.glob(path):
    labels.append(filename.split(".")[0])
    data = np.genfromtxt(filename, skip_header=1)
    for i in range(data.shape[1]-1):
        ax[i].plot(data[:,0], data[:,i+1])

for i in range(numsubplots):
    ax[i].set_xlabel(headers[0])
    ax[i].set_ylabel(headers[i+1])
    ax[i].grid()

plt.legend(labels)
plt.show()
