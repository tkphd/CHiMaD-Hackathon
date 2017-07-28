#!/usr/bin/python
import matplotlib.pylab as plt
import numpy as np
from os import path
import glob, re

# Plot energy data

plt.figure(1)

if path.isfile('../energy_dx10.tsv'):
    t, e = np.loadtxt('../energy_dx10.tsv',usecols=(0, 1), unpack=True, skiprows=2, delimiter='\t')
    plt.plot(t, e, label=r"$\Delta x=1.0$")

if path.isfile('../energy_dx05.tsv'):
    t, e = np.loadtxt('../energy_dx05.tsv',usecols=(0, 1), unpack=True, skiprows=2, delimiter='\t')
    plt.plot(t, e, label=r"$\Delta x=0.5$")

if path.isfile('../energy_dx01.tsv'):
    t, e = np.loadtxt('../energy_dx01.tsv',usecols=(0, 1), unpack=True, skiprows=2, delimiter='\t')
    plt.plot(t, e, label=r"$\Delta x=0.1$")

plt.xlabel(r'Time')
plt.ylabel(r'$\mathcal{F}=\sum\sum\left[f_{\mathrm{chem}} + f_{\mathrm{elec}} + \frac{\kappa}{2}|\nabla c|^2\right]\Delta x\Delta y$')
plt.title("Cahn-Hilliard-Poisson Free Energy")
plt.legend(loc='best')
plt.savefig('energy.png', bbox_inches='tight', dpi=400)
plt.close()

# Plot linescan data

plt.figure(2)
fnames = sorted(glob.glob("../square/electrochem_dx10.*.csv"))
if len(fnames) > 0:
    n = 0
    for file in fnames:
        num = int(re.search('[0-9]{3,6}', file).group(0)) / 20
        x, c, u, p = np.loadtxt(file, delimiter=',', skiprows=1, unpack=True)
        plt.plot(x, c+0.5*n, label="t={0}".format(num))
        n += 1
    plt.xlabel("Position $(x)$")
    plt.ylabel(r"Stacked Composition $(c\in[0,1])$")
    plt.title("Cahn-Hilliard-Poisson composition $(\Delta x=1.0, y=50)$")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('comp_dx10.png', bbox_inches='tight', dpi=400)
    plt.close()

plt.figure(3)
fnames = sorted(glob.glob("../square/electrochem_dx05.*.csv"))
if len(fnames) > 0:
    n = 0
    for file in fnames:
        num = int(re.search('[0-9]{3,6}', file).group(0)) / 20
        x, c, u, p = np.loadtxt(file, delimiter=',', skiprows=1, unpack=True)
        plt.plot(x, c+0.5*n, label="t={0}".format(num))
        n += 1
    plt.xlabel("Position $(x)$")
    plt.ylabel(r"Stacked Composition $(c\in[0,1])$")
    plt.title("Cahn-Hilliard-Poisson composition $(\Delta x=0.5, y=50)$")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('comp_dx05.png', bbox_inches='tight', dpi=400)
    plt.close()

plt.figure(4)
fnames = sorted(glob.glob("../square/electrochem_dx01.*.csv"))
if len(fnames) > 0:
    n = 0
    for file in fnames:
        num = int(re.search('[0-9]{3,6}', file).group(0)) / 20
        x, c, u, p = np.loadtxt(file, delimiter=',', skiprows=1, unpack=True)
        plt.plot(x, c+0.5*n, label="t={0}".format(num))
        n += 1
    plt.xlabel("Position $(x)$")
    plt.ylabel(r"Stacked Composition $(c\in[0,1])$")
    plt.title("Cahn-Hilliard-Poisson composition $(\Delta x=0.1, y=50)$")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('comp_dx01.png', bbox_inches='tight', dpi=400)
    plt.close()
