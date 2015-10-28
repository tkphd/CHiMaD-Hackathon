#!/usr/bin/python
import matplotlib.pylab as plt
import numpy as np

dt = 0.00375

e,c = np.loadtxt('energy.log',usecols=(0, 1), unpack=True, skiprows=1, delimiter='\t')
t = dt*np.arange(0,len(e))

#plt.plot(t, e)
plt.loglog(t, e-np.min(e))
plt.xlabel(r'Time (using $\Delta t=0.00375$)')
plt.ylabel(r'$\mathcal{F}-\mathcal{F}_{min}=\sum_{i,j} f(c,\vec{\phi})\Delta x\Delta y$ (arb. units)')
plt.xlim([1, 1e5])
plt.title("T-square Ostwald")
plt.savefig('energy.png', bbox_inches='tight', dpi=300)
plt.close()

#plt.plot(t, c)
plt.semilogy(t, c)
plt.xlabel(r'Time (using $\Delta t=0.00375$)')
plt.ylabel(r'$m=\sum\sum c\Delta x\Delta y$ (arb. units)')
plt.ylim([10, 5e3])
plt.title("T-square Ostwald")
plt.savefig('mass.png', bbox_inches='tight', dpi=300)
plt.close()
