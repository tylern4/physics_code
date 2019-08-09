# distutils: language = c++
import cython
from libc.stdlib cimport free

cdef dict get_id = {'PROTON': 2212, 'NEUTRON': 2112, 'PIP': 211, 'PIM': -211, 'PI0': 111, 'KP': 321, 'KM': -321, 'PHOTON': 22, 'ELECTRON': 11}
cdef dict part_mass = {11: 0.000511, 211: 0.13957, -211: 0.13957, 2212: 0.93827, 2112: 0.939565, 321: 0.493667, -321: 0.493667, 22: 0}

cdef extern from "constants.hpp":
  cdef cppclass LorentzVector:
    LorentzVector() except +
    LorentzVector(double px, double py, double pz, double m) except +
 
cdef class FourVector:
  cdef LorentzVector*c_LorentzVector
  def __cinit__(LorentzVector self, double px, double py, double pz, **kwargs):
    if "mass" in kwargs:
      self.c_LorentzVector = new LorentzVector(px, py, pz, kwargs["mass"])
    elif "pid" in kwargs:
      self.c_LorentzVector = new LorentzVector(px, py, pz, part_mass[kwargs["pid"]])
