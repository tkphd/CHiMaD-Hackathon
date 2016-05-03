#!/usr/bin/python
import matplotlib.pylab as plt
import numpy as np

dt = 0.000439453

e,c = np.loadtxt('energy.log',usecols=(0, 1), unpack=True, skiprows=1, delimiter='\t')
t = dt*np.arange(0,len(e))

#plt.plot(t, e)
plt.semilogy(t, e-np.min(e))
plt.xlabel(r'Time (using $\Delta t=0.0004, Co=0.25$)')
plt.ylabel(r'$\mathcal{F}-\mathcal{F}_{min}=\sum\sum f(c)\Delta x\Delta y$ (arb. units)')
plt.title(r'Periodic Cahn-Hilliard / $\sqrt{2}$')
plt.savefig('energy.png', bbox_inches='tight', dpi=300)
plt.close()

#plt.plot(t, c)
plt.semilogy(t, c)
plt.xlabel(r'Time (using $\Delta t=0.0004$, Co=0.25)')
plt.ylabel(r'$m=\sum\sum c\Delta x\Delta y$ (arb. units)')
plt.ylim([10, 1e5])
plt.title(r'Periodic Cahn-Hilliard / $\sqrt{2}$')
plt.savefig('mass.png', bbox_inches='tight', dpi=300)
plt.close()
