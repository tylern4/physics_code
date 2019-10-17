from build.h10 import h10_data
from build.physics_vectors import LorentzVector
import numpy as np


class reaction:
    """A simple reaction class"""
    target = LorentzVector(0, 0, 0, pid=2212)
    beam = LorentzVector(0, 0, 4.81726, pid=11)
    W = np.nan
    Q2 = np.nan

    def __init__(self, *args, **kwargs):
        self._data = args[0]
        if self._data.gpart < 1:
            return
        self.elec = LorentzVector(
            self._data.px[0], self._data.py[0], self._data.pz[0], pid=11)
        self._calc_W()
        self._calc_Q2()

    def _calc_W(self):
        x = ((self.beam - self.elec) + self.target)
        self.W = x.mass

    def _calc_Q2(self):
        x = (self.beam - self.elec)
        self.Q2 = -x.M2
