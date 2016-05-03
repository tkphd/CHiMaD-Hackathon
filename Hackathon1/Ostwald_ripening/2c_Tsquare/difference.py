#!/users/tnk10/.conda/envs/mmsp/bin/python
# -*- coding: utf-8 -*-
import matplotlib.pylab as plt
from matplotlib import cm
import pandas as pd
import numpy as np

colors = cm.Paired(np.linspace(0,1,12))

d50  = pd.read_csv('lin/ostwald.05000000.csv',encoding='utf-8')
d125 = pd.read_csv('lin/ostwald.12500000.csv',encoding='utf-8')

fig = plt.figure()
ax = plt.subplot(111)

ax.plot(d50['y'], d125['f'] - d50['f'], label=r'$f$', color=colors[0])
ax.plot(d50['y'], d125['p0'] - d50['p0'], label=r'$c$', color=colors[1])
ax.plot(d50['y'], d125['p1'] - d50['p1'], label=r'$\phi_1$', color=colors[2])
ax.plot(d50['y'], d125['p2'] - d50['p2'], label=r'$\phi_2$', color=colors[3])
ax.plot(d50['y'], d125['p3'] - d50['p3'], label=r'$\phi_3$', color=colors[4])
ax.plot(d50['y'], d125['p4'] - d50['p4'], label=r'$\phi_4$', color=colors[5])
ax.plot(d50['y'], d125['p5'] - d50['p5'], label=r'$\phi_5$', color=colors[6])
ax.plot(d50['y'], d125['p6'] - d50['p6'], label=r'$\phi_6$', color=colors[7])
ax.plot(d50['y'], d125['p7'] - d50['p7'], label=r'$\phi_7$', color=colors[8])
ax.plot(d50['y'], d125['p8'] - d50['p8'], label=r'$\phi_8$', color=colors[9])
ax.plot(d50['y'], d125['p9'] - d50['p9'], label=r'$\phi_9$', color=colors[10])
ax.plot(d50['y'], d125['p10'] - d50['p10'], label=r'$\phi_{10}$', color=colors[11])

plt.title("Differences: 12.5M - 5M")
ax.legend(bbox_to_anchor=(1.25,1.0))
plt.savefig('differences.png',bbox_inches='tight',dpi=300)
plt.close()

e50 = d50['f'].sum()
e125 = d125['f'].sum()

print "Total energy: %e (5M), %e (10M)."%(e50, e125)
if (e50<e125):
    print "Energy is unstable: df = %e > 0."%(e125-e50)
