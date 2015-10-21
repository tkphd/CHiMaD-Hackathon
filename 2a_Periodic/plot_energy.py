#!/usr/bin/python
import matplotlib.pylab as plt
import numpy as np

e,c = np.loadtxt('energy.log',usecols=(0, 1), unpack=True)
x = 0.0025*np.arange(0,len(e))

plt.plot(x, e)
plt.xlabel(r'Time (using $\Delta t=0.0025$)')
plt.ylabel(r'$\mathcal{F}_0=\sum\sum f_0\Delta x\Delta y$ (arb. units)')
plt.title("Periodic Ostwald")
plt.savefig('energy.png', bbox_inches='tight')
