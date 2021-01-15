# distutils: language = c++
import cython
from libc.stdlib cimport free
from libcpp.memory cimport unique_ptr, shared_ptr
from cython.operator cimport dereference as deref
from libcpp cimport bool
import numpy as np
cimport numpy as np
from libc.math cimport sin, cos, sqrt
from scipy import stats
cimport cython

cdef float MP = 0.93827208816
cdef float E0 = 4.81726
cdef float ME = 0.00051099895

cdef float p_targ_px = 0.0
cdef float p_targ_py = 0.0
cdef float p_targ_pz = 0.0
cdef float p_targ_E = MP

cdef float e_beam_px = 0.0
cdef float e_beam_py = 0.0
cdef float e_beam_pz = sqrt(E0**2-ME**2)
cdef float e_beam_E = E0

cdef dict get_id = {'PROTON': 2212, 'NEUTRON': 2112, 'PIP': 211, 'PIM': -211, 'PI0': 111, 'KP': 321, 'KM': -321, 'PHOTON': 22, 'ELECTRON': 11}
cdef dict part_mass = {11: 0.000511, 211: 0.13957, -211: 0.13957, 2212: 0.93827, 2112: 0.939565, 321: 0.493667, -321: 0.493667, 22: 0}

cdef extern from "TLorentzVector.h":
  cdef cppclass TLorentzVector:
    TLorentzVector() except +
    TLorentzVector(double x, double y, double z, double t) except +
    TLorentzVector operator*(double)
    double operator*(TLorentzVector)
    TLorentzVector operator+(TLorentzVector)
    TLorentzVector operator-(TLorentzVector)
    bool operator!=(TLorentzVector)
    bool operator==(TLorentzVector)
    void Boost(double, double, double)
    void SetXYZM(double x, double y, double z, double m)
    void SetXYZT (double x, double y, double z, double t)
    double Px()
    double Py()
    double Pz()
    double P()
    double E()
    double Energy()
    double Theta()
    double CosTheta()
    double Phi()
    double Rho()
    double Perp2()
    double Pt()
    double Perp()
    double Et2()
    double Et()
    double Mag2()
    double M2()
    double Mag()
    double M()
    double Mt2()
    double Mt()
    double Beta()
    double Gamma()
    double Plus()
    double Minus()
    double Rapidity()
    double Eta()
    double PseudoRapidity()

cdef class LorentzVector:
  cdef unique_ptr[TLorentzVector] c_TLorentzVector
  def __cinit__(LorentzVector self, double px, double py, double pz, **kwargs):
    if "energy" in kwargs:
      self.c_TLorentzVector.reset(new TLorentzVector(px, py, pz, kwargs["energy"]))
    elif "mass" in kwargs:
      self.c_TLorentzVector.reset(new TLorentzVector())
      deref(self.c_TLorentzVector).SetXYZM(px, py, pz, kwargs["mass"])
    elif "pid" in kwargs:
      self.c_TLorentzVector.reset(new TLorentzVector())
      deref(self.c_TLorentzVector).SetXYZM(px, py, pz, part_mass[kwargs["pid"]])
    else:
      self.c_TLorentzVector.reset(new TLorentzVector(px, py, pz, 0))
  def __add__(LorentzVector self, LorentzVector other):
    cdef double X = deref(self.c_TLorentzVector).Px() + deref(other.c_TLorentzVector).Px()
    cdef double Y = deref(self.c_TLorentzVector).Py() + deref(other.c_TLorentzVector).Py()
    cdef double Z = deref(self.c_TLorentzVector).Pz() + deref(other.c_TLorentzVector).Pz()
    cdef double E = deref(self.c_TLorentzVector).E() + deref(other.c_TLorentzVector).E()
    return LorentzVector(X, Y, Z, energy=E)
  def __iadd__(LorentzVector self, LorentzVector other):
    cdef double X = deref(self.c_TLorentzVector).Px() + deref(other.c_TLorentzVector).Px()
    cdef double Y = deref(self.c_TLorentzVector).Py() + deref(other.c_TLorentzVector).Py()
    cdef double Z = deref(self.c_TLorentzVector).Pz() + deref(other.c_TLorentzVector).Pz()
    cdef double E = deref(self.c_TLorentzVector).E() + deref(other.c_TLorentzVector).E()
    return LorentzVector(X, Y, Z, energy=E)
  def __sub__(LorentzVector self, LorentzVector other):
    cdef double X = deref(self.c_TLorentzVector).Px() - deref(other.c_TLorentzVector).Px()
    cdef double Y = deref(self.c_TLorentzVector).Py() - deref(other.c_TLorentzVector).Py()
    cdef double Z = deref(self.c_TLorentzVector).Pz() - deref(other.c_TLorentzVector).Pz()
    cdef double E = deref(self.c_TLorentzVector).E() - deref(other.c_TLorentzVector).E()
    return LorentzVector(X, Y, Z, energy=E)
  def __isub__(LorentzVector self, LorentzVector other):
    cdef double X = deref(self.c_TLorentzVector).Px() - deref(other.c_TLorentzVector).Px()
    cdef double Y = deref(self.c_TLorentzVector).Py() - deref(other.c_TLorentzVector).Py()
    cdef double Z = deref(self.c_TLorentzVector).Pz() - deref(other.c_TLorentzVector).Pz()
    cdef double E = deref(self.c_TLorentzVector).E() - deref(other.c_TLorentzVector).E()
    return LorentzVector(X, Y, Z, energy=E)
  def __str__(self):
    return f"Px {self.Px: 0.4f} | Py {self.Py: 0.4f} | Pz {self.Pz: 0.4f} | E {self.E: 0.4f} | M {self.mass: 0.4f}"
  def __repr__(self):
    return self.__str__()
  def SetPxPyPzM(LorentzVector self, double px, double py, double pz, double mass):
    deref(self.c_TLorentzVector).SetXYZM(px, py, pz, mass)
  @property
  def px(LorentzVector self):
    return deref(self.c_TLorentzVector).Px()
  @property
  def py(LorentzVector self):
    return deref(self.c_TLorentzVector).Py()
  @property
  def pz(LorentzVector self):
    return deref(self.c_TLorentzVector).Pz()

  cdef float Px(LorentzVector self):
    return deref(self.c_TLorentzVector).Px()

  cdef float Py(LorentzVector self):
    return deref(self.c_TLorentzVector).Py()

  cdef float Pz(LorentzVector self):
    return deref(self.c_TLorentzVector).Pz()
  @property
  def P(LorentzVector self):
    return deref(self.c_TLorentzVector).P()
  @property
  def E(LorentzVector self):
    return deref(self.c_TLorentzVector).E()
  cdef float Energy(LorentzVector self):
    return deref(self.c_TLorentzVector).E()
  @property
  def Theta(LorentzVector self):
    return deref(self.c_TLorentzVector).Theta()
  @property
  def CosTheta(LorentzVector self):
    return deref(self.c_TLorentzVector).CosTheta()
  @property
  def Phi(LorentzVector self):
    return deref(self.c_TLorentzVector).Phi()
  @property
  def Rho(LorentzVector self):
    return deref(self.c_TLorentzVector).Rho()
  @property
  def Perp2(LorentzVector self):
    return deref(self.c_TLorentzVector).Perp2()
  @property
  def Pt(LorentzVector self):
    return deref(self.c_TLorentzVector).Pt()
  @property
  def Perp(LorentzVector self):
    return deref(self.c_TLorentzVector).Perp()
  @property
  def Et2(LorentzVector self):
    return deref(self.c_TLorentzVector).Et2()
  @property
  def Et(LorentzVector self):
    return deref(self.c_TLorentzVector).Et()
  @property
  def Mag2(LorentzVector self):
    return deref(self.c_TLorentzVector).Mag2()
  @property
  def M2(LorentzVector self):
    return deref(self.c_TLorentzVector).M2()
  @property
  def Mag(LorentzVector self):
    return deref(self.c_TLorentzVector).Mag()
  @property
  def mass(LorentzVector self):
    return deref(self.c_TLorentzVector).M()
  @property
  def M(LorentzVector self):
    return deref(self.c_TLorentzVector).M()
  @property
  def Mt2(LorentzVector self):
    return deref(self.c_TLorentzVector).Mt2()
  @property
  def Mt(LorentzVector self):
    return deref(self.c_TLorentzVector).Mt()
  @property
  def Beta(LorentzVector self):
    return deref(self.c_TLorentzVector).Beta()
  @property
  def Gamma(LorentzVector self):
    return deref(self.c_TLorentzVector).Gamma()
  @property
  def Plus(LorentzVector self):
    return deref(self.c_TLorentzVector).Plus()
  @property
  def Minus(LorentzVector self):
    return deref(self.c_TLorentzVector).Minus()
  @property
  def Rapidity(LorentzVector self):
    return deref(self.c_TLorentzVector).Rapidity()
  @property
  def Eta(LorentzVector self):
    return deref(self.c_TLorentzVector).Eta()
  @property
  def PseudoRapidity(LorentzVector self):
    return deref(self.c_TLorentzVector).PseudoRapidity()


def calc_W(LorentzVector x):
    cdef float e_prime_px = x.Px()
    cdef float e_prime_py = x.Py()
    cdef float e_prime_pz = x.Pz()
    cdef float e_prime_E = x.Energy()
    
    cdef float temp_px = e_beam_px - e_prime_px + p_targ_px
    cdef float temp_py = e_beam_py - e_prime_py + p_targ_py
    cdef float temp_pz = e_beam_pz - e_prime_pz + p_targ_pz
    cdef float temp_E = e_beam_E - e_prime_E + p_targ_E
    
    
    cdef float temp2 = temp_px**2+temp_py**2+temp_pz**2-temp_E**2
    cdef float temp3 = sqrt(-temp2)
    
    return temp3


def calc_q2(LorentzVector x):
    cdef float e_prime_px = x.Px()
    cdef float e_prime_py = x.Py()
    cdef float e_prime_pz = x.Pz()
    cdef float e_prime_E = x.Energy()
    
    cdef float temp_px = e_beam_px - e_prime_px
    cdef float temp_py = e_beam_py - e_prime_py
    cdef float temp_pz = e_beam_pz - e_prime_pz
    cdef float temp_E = e_beam_E - e_prime_E

    cdef float temp2 = temp_px**2+temp_py**2+temp_pz**2-temp_E**2

    return temp2