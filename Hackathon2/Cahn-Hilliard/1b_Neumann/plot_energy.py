#!/usr/bin/python
import matplotlib.pylab as plt
import numpy as np

dt = 0.005

e,c = np.loadtxt('energy.log',usecols=(0, 1), unpack=True, skiprows=1, delimiter='\t')
t = dt*np.arange(0,len(e))

#plt.plot(t, e)
plt.semilogy(t, e-np.min(e))
plt.xlabel(r'Time (using $\Delta t=0.005$)')
plt.ylabel(r'$\mathcal{F}-\mathcal{F}_{min}=\sum\sum f_0\Delta x\Delta y$ (arb. units)')
plt.title("Neumann Cahn-Hilliard")
plt.savefig('energy.png', bbox_inches='tight', dpi=300)
plt.close()

#plt.plot(t, c)
plt.semilogy(t, c)
plt.xlabel(r'Time (using $\Delta t=0.005$)')
plt.ylabel(r'$m=\sum\sum c\Delta x\Delta y$ (arb. units)')
plt.ylim([10, 5e4])
plt.title("Neumann Cahn-Hilliard")
plt.savefig('mass.png', bbox_inches='tight', dpi=300)
plt.close()
