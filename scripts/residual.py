#!/usr/bin/python
#shamelessly stolen from stack overflow

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('RMSError.txt')

fig, ax = plt.subplots(1)

#u'\u03C1 is lower case rho
if data.shape[1]-1 == 3:
    labels = [u'\u03C1',u'\u03C1u',u'\u03C1E']
if data.shape[1]-1 == 4:
    labels = [u'\u03C1',u'\u03C1u',u'\u03C1v',u'\u03C1E']
if data.shape[1]-1 == 5:
    labels = [u'\u03C1',u'\u03C1u',u'\u03C1v',u'\u03C1w',u'\u03C1E']

for i in range(data.shape[1]-1):
    ax.semilogy(data[:,0],abs(data[:,i+1]))

plt.legend(labels)
ax.set_xlabel('Iteration')
ax.set_ylabel('Residual')
plt.grid()
plt.show()
