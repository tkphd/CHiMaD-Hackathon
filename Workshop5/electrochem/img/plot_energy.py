#!/usr/bin/python
import matplotlib.pylab as plt
import numpy as np

dt = 0.05

t, e = np.loadtxt('../energy.log',usecols=(0, 1), unpack=True, skiprows=1, delimiter='\t')

#plt.plot(t, e)
plt.semilogy(t, e-np.min(e))
plt.xlabel(r'Time')
plt.ylabel(r'$\mathcal{F}-\mathcal{F}_{min}=\sum\sum f_0\Delta x\Delta y$ (arb. units)')
plt.title("Cahn-Hilliard-Poisson $(\Delta t=0.05)$")
plt.savefig('energy.png', bbox_inches='tight', dpi=300)
plt.close()
