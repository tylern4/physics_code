from __future__ import print_function
import uproot
import math
import sys


import matplotlib.pyplot as plt

tree = uproot.open(sys.argv[1])["h10"]

pid, p, cx, cy, cz = tree.arrays(
    ["id", "p", "cx", "cy", "cz"], outputtype=tuple)

momentum = []
for pid, pi, cxi, cyi, czi in zip(pid, p, cx, cy, cz):
    if(len(pid) > 0 and pid[0] == 11):
        for pj, cxj, cyj, czj in zip(pi, cxi, cyi, czi):
            momentum.append(math.sqrt((pj * cxj)**2 +
                                      (pj * cyj)**2 + (pj * czj)**2))


plt.hist(momentum, bins=500, range=(0, 5))
plt.show()
