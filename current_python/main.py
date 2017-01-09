import ROOT
import numpy as np
from ROOT import TLorentzVector,TVector3
import matplotlib.pyplot as plt
import matplotlib.cm as cm

Square = lambda x: x**2
append = lambda _arr,_val: np.append(_arr,_val)
array = lambda _:np.array([])


def all_mom(p,cx,cy,cz):
    calc = lambda _p, _cos: _p * _cos
    return calc(p,cx),calc(p,cy),calc(p,cz)

def fvec(_px,_py,_pz,_mass):
    _vec = TVector3(_px,_py,_pz)
    _4vec = TLorentzVector()
    _4vec.SetVectM(_vec,_mass)
    return _4vec

def fourvec(_p,_cx,_cy,_cz,_mass):
    _px,_py,_pz = all_mom(_p,_cx,_cy,_cz)
    return fvec(_px,_py,_pz,_mass)

my_cmap = cm.get_cmap('viridis')
my_cmap.set_over('w')
color_map = my_cmap
plt.rc('text', usetex=True)

num_bins = 300
fig_size = (16, 9)