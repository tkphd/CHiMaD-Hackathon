#!/usr/bin/python
from matplotlib import pylab as plt
import numpy as np
from sympy import diff, expand, factor, simplify, symbols

c, ca, cb, k, p, r = symbols("c ca cb k p r")

F = r * (c-ca)**2 * (cb-c)**2 #+ 0.5*k*c*p

f = diff(F, c)

print(simplify(expand(F)))

print(simplify(expand(f)))

def Fcon(c):
    ca = 0.3
    cb = 0.7
    r = 5
    return r * (c**4 + (ca**2 + 4 * ca * cb + cb**2) * c**2 + ca**2 * cb**2)

def Fexp(c):
    ca = 0.3
    cb = 0.7
    r = 5
    k = 0.09
    return r * (-2 * (ca + cb) * c**3 - 2 * ca * cb * (ca + cb) * c) - 0.5 * k * c

x = np.linspace(0, 1, 100)
ycon = Fcon(x)
yexp = Fexp(x)

plt.figure()
plt.plot(x, ycon, label="con")
plt.plot(x, yexp, label="exp")
plt.plot(x, ycon+yexp, label="tot")
plt.xlim([0, 1])
plt.ylim([-0.5, 0.5])
plt.legend(loc='best')
plt.savefig('f.png', bbox_inches='tight', dpi=400)
