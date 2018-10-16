"""
LFV
===

The LFV diagram.
"""
import matplotlib.pyplot as plt
from feynman import Diagram

fig = plt.figure(figsize=(8., 8.))
ax = fig.add_axes([0, 0, 1, 1], frameon=False)

diagram = Diagram(ax)
in1 = diagram.vertex(xy=(.5, .5))
in2 = diagram.vertex(xy=(.65, .5))
#v1 = diagram.vertex(xy=(.1, .5))
#v2 = diagram.vertex(xy=(.5, .35))
#out1 = diagram.vertex(xy=(.9, .65), marker='')
#out2 = diagram.vertex(xy=(.9, .35), marker='')

NStar = diagram.line(in1, in2, arrow=False, style='double')
#nu1 = diagram.line(v1, in2)
#nu2 = diagram.line(in2, v2)
#w = diagram.line(v1, v2, style='wiggly')
#lep = diagram.line(out1, v1)
#tau = diagram.line(v2, out2)

#nu1.text(r"$\nu_\ell$", fontsize=40)
#nu2.text(r"$\nu_\tau$", fontsize=40)
#lep.text(r"$\ell^+$", fontsize=40)
#tau.text(r"$\tau^-$", fontsize=40)
NStar.text(r"$N^*$", fontsize=40)

diagram.plot()
plt.show()
