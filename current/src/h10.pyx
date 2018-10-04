# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cpython cimport array
import numpy as np
cimport numpy as np
import array

cdef dict get_id = {'PROTON': 2212, 'NEUTRON': 2112, 'PIP': 211, 'PIM': -211, 'PI0': 111, 'KP': 321, 'KM': -321, 'PHOTON': 22, 'ELECTRON': 11}

cdef dict part_mass = {11: 0.000511, 211: 0.13957, -211: 0.13957, 2212: 0.93827, 2112: 0.939565, 321: 0.493667, -321: 0.493667, 22: 0}

cdef extern from "TChain.h":
  cdef cppclass TChain:
    TChain(char*) except +
    int Add(char*)
    long GetEntries()
    int GetEntry(long)

cdef extern from "TLorentzVector.h":
  cdef cppclass TLorentzVector:
    TLorentzVector() except +
    TLorentzVector(double x, double y, double z, double t) except +
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

cdef extern from "TVector3.h":
  cdef cppclass TVector3:
    TVector3 () except +
    TVector3 (double x, double y, double z) except +
    void SetXYZ(double x, double y, double z)
    double x()
    double y()
    double z()
    double X()
    double Y()
    double Z()
    double Px()
    double Py()
    double Pz()
    double Phi()
    double Theta()
    double CosTheta()
    double Mag2()
    double Mag()
    double Perp2()
    double Pt()
    double Perp()
    double Eta ()
    void RotateX (double x)
    void RotateY (double x)
    void RotateZ (double x)
    void Print()

cdef class LorentzVector:
  cdef TLorentzVector*c_TLorentzVector
  def __cinit__(LorentzVector self, double px, double py, double pz, **kwargs):
    if "energy" in kwargs:
      self.c_TLorentzVector = new TLorentzVector(px, py, pz, kwargs["energy"])
    elif "mass" in kwargs:
      self.c_TLorentzVector = new TLorentzVector()
      self.c_TLorentzVector.SetXYZM(px, py, pz, kwargs["mass"])
    elif "pid" in kwargs:
      self.c_TLorentzVector = new TLorentzVector()
      self.c_TLorentzVector.SetXYZM(px, py, pz, part_mass[kwargs["pid"]])
    else:
      self.c_TLorentzVector = new TLorentzVector(px, py, pz, 0)
  def __dealloc__(self):
    del self.c_TLorentzVector
  def __add__(LorentzVector self, LorentzVector other):
    cdef double X = self.c_TLorentzVector.Px() + other.c_TLorentzVector.Px()
    cdef double Y = self.c_TLorentzVector.Py() + other.c_TLorentzVector.Py()
    cdef double Z = self.c_TLorentzVector.Pz() + other.c_TLorentzVector.Pz()
    cdef double E = self.c_TLorentzVector.E() + other.c_TLorentzVector.E()
    return LorentzVector(X, Y, Z, energy=E)
  def __iadd__(LorentzVector self, LorentzVector other):
    cdef double X = self.c_TLorentzVector.Px() + other.c_TLorentzVector.Px()
    cdef double Y = self.c_TLorentzVector.Py() + other.c_TLorentzVector.Py()
    cdef double Z = self.c_TLorentzVector.Pz() + other.c_TLorentzVector.Pz()
    cdef double E = self.c_TLorentzVector.E() + other.c_TLorentzVector.E()
    return LorentzVector(X, Y, Z, energy=E)
  def __sub__(LorentzVector self, LorentzVector other):
    cdef double X = self.c_TLorentzVector.Px() - other.c_TLorentzVector.Px()
    cdef double Y = self.c_TLorentzVector.Py() - other.c_TLorentzVector.Py()
    cdef double Z = self.c_TLorentzVector.Pz() - other.c_TLorentzVector.Pz()
    cdef double E = self.c_TLorentzVector.E() - other.c_TLorentzVector.E()
    return LorentzVector(X, Y, Z, energy=E)
  def __isub__(LorentzVector self, LorentzVector other):
    cdef double X = self.c_TLorentzVector.Px() - other.c_TLorentzVector.Px()
    cdef double Y = self.c_TLorentzVector.Py() - other.c_TLorentzVector.Py()
    cdef double Z = self.c_TLorentzVector.Pz() - other.c_TLorentzVector.Pz()
    cdef double E = self.c_TLorentzVector.E() - other.c_TLorentzVector.E()
    return LorentzVector(X, Y, Z, energy=E)
  def __str__(self):
    return "Px {0: 0.2f} | Py {1: 0.2f} | Pz {2: 0.2f} | E {3: 0.2f}".format(self.px,self.py ,self.pz, self.energy)
  def __repr__(self):
    return self.__str__()
  def MomentumVec(LorentzVector self):
    return ThreeVector(self.c_TLorentzVector.Px(), self.c_TLorentzVector.Py(), self.c_TLorentzVector.Pz())
  def SetPxPyPzM(LorentzVector self, double px, double py, double pz, double mass):
    self.c_TLorentzVector.SetXYZM(px, py, pz, mass)
  @property
  def Px(LorentzVector self):
    return self.c_TLorentzVector.Px()
  @property
  def Py(LorentzVector self):
    return self.c_TLorentzVector.Py()
  @property
  def Pz(LorentzVector self):
    return self.c_TLorentzVector.Pz()
  @property
  def P(LorentzVector self):
    return self.c_TLorentzVector.P()
  @property
  def E(LorentzVector self):
    return self.c_TLorentzVector.E()
  @property
  def Energy(LorentzVector self):
    return self.c_TLorentzVector.E()
  @property
  def Theta(LorentzVector self):
    return self.c_TLorentzVector.Theta()
  @property
  def CosTheta(LorentzVector self):
    return self.c_TLorentzVector.CosTheta()
  @property
  def Phi(LorentzVector self):
    return self.c_TLorentzVector.Phi()
  @property
  def Rho(LorentzVector self):
    return self.c_TLorentzVector.Rho()
  @property
  def Perp2(LorentzVector self):
    return self.c_TLorentzVector.Perp2()
  @property
  def Pt(LorentzVector self):
    return self.c_TLorentzVector.Pt()
  @property
  def Perp(LorentzVector self):
    return self.c_TLorentzVector.Perp()
  @property
  def Et2(LorentzVector self):
    return self.c_TLorentzVector.Et2()
  @property
  def Et(LorentzVector self):
    return self.c_TLorentzVector.Et()
  @property
  def Mag2(LorentzVector self):
    return self.c_TLorentzVector.Mag2()
  @property
  def M2(LorentzVector self):
    return self.c_TLorentzVector.M2()
  @property
  def Mag(LorentzVector self):
    return self.c_TLorentzVector.Mag()
  @property
  def M(LorentzVector self):
    return self.c_TLorentzVector.M()
  @property
  def Mt2(LorentzVector self):
    return self.c_TLorentzVector.Mt2()
  @property
  def Mt(LorentzVector self):
    return self.c_TLorentzVector.Mt()
  @property
  def Beta(LorentzVector self):
    return self.c_TLorentzVector.Beta()
  @property
  def Gamma(LorentzVector self):
    return self.c_TLorentzVector.Gamma()
  @property
  def Plus(LorentzVector self):
    return self.c_TLorentzVector.Plus()
  @property
  def Minus(LorentzVector self):
    return self.c_TLorentzVector.Minus()
  @property
  def Rapidity(LorentzVector self):
    return self.c_TLorentzVector.Rapidity()
  @property
  def Eta(LorentzVector self):
    return self.c_TLorentzVector.Eta()
  @property
  def PseudoRapidity(LorentzVector self):
    return self.c_TLorentzVector.PseudoRapidity()

cdef class ThreeVector:
  cdef TVector3*c_TVector3
  def __cinit__(ThreeVector self, double vx, double vy, double vz):
    c_TVector3 = new TVector3(vx, vy, vz)
  def __add__(ThreeVector self, ThreeVector other):
    cdef double X = self.c_TVector3.x() + other.c_TVector3.x()
    cdef double Y = self.c_TVector3.y() + other.c_TVector3.y()
    cdef double Z = self.c_TVector3.z() + other.c_TVector3.z()
    return ThreeVector(X, Y, Z)
  def __str__(ThreeVector self):
    return "Vx {0: 0.2f} | Vy {1: 0.2f} | Vz {2: 0.2f}".format(self.vx,self.vy ,self.vz)
  def __repr__(self):
    return self.__str__()
  @property
  def x(ThreeVector self):
    return self.c_TVector3.x()
  @property
  def y(ThreeVector self):
    return self.c_TVector3.y()
  @property
  def z(ThreeVector self):
    return self.c_TVector3.z()
  @property
  def X(ThreeVector self):
    return self.c_TVector3.X()
  @property
  def Y(ThreeVector self):
    return self.c_TVector3.Y()
  @property
  def Z(ThreeVector self):
    return self.c_TVector3.Z()
  @property
  def Px(ThreeVector self):
    return self.c_TVector3.Px()
  @property
  def Py(ThreeVector self):
    return self.c_TVector3.Py()
  @property
  def Pz(ThreeVector self):
    return self.c_TVector3.Pz()
  @property
  def Phi(ThreeVector self):
    return self.c_TVector3.Phi()
  @property
  def Theta(ThreeVector self):
    return self.c_TVector3.Theta()
  @property
  def CosTheta(ThreeVector self):
    return self.c_TVector3.CosTheta()
  @property
  def Mag2(ThreeVector self):
    return self.c_TVector3.Mag2()
  @property
  def Mag(ThreeVector self):
    return self.c_TVector3.Mag()
  @property
  def Perp2(ThreeVector self):
    return self.c_TVector3.Perp2()
  @property
  def Pt(ThreeVector self):
    return self.c_TVector3.Pt()
  @property
  def Perp(ThreeVector self):
    return self.c_TVector3.Perp()
  @property
  def Eta (ThreeVector self):
    return self.c_TVector3.Eta()

cdef extern from "physics.cpp":
  pass

cdef extern from "physics.hpp" namespace "physics":
  double theta_calc(double cosz)
  double phi_calc(double cosx, double cosy)
  double center_phi_calc(double cosx, double cosy)
  int get_sector(double phi)
  double Get_Mass(int ID)
  double fiducial_phi(double theta_e, double e_p)


cdef extern from "branches.cpp":
  pass

cdef extern from "branches.hpp":
    cdef cppclass Branches:
      Branches(TChain*) except +
      int npart()
      int evstat()
      int intt()
      int evntid()
      int evtype()
      int evntclas()
      int evthel()
      int evntclas2()
      float q_l()
      float t_l()
      float tr_time()
      float rf_time1()
      float rf_time2()
      int gpart()
      vector[int] id()
      vector[int] dc()
      vector[int] cc()
      vector[int] sc()
      vector[int] ec()
      vector[float] p()
      vector[int] q()
      vector[float] b()
      vector[float] cx()
      vector[float] cy()
      vector[float] cz()
      vector[float] vx()
      vector[float] vy()
      vector[float] vz()
      vector[int] dc_sect()
      vector[int] dc_trk()
      vector[int] dc_stat()
      vector[float] dc_vx()
      vector[float] dc_vy()
      vector[float] dc_vz()
      vector[float] dc_vr()
      vector[float] dc_xsc()
      vector[float] dc_ysc()
      vector[float] dc_zsc()
      vector[float] dc_cxsc()
      vector[float] dc_cysc()
      vector[float] dc_czsc()
      vector[float] dc_c2()
      vector[int] ec_stat()
      vector[int] ec_sect()
      vector[int] ec_whol()
      vector[int] ec_inst()
      vector[int] ec_oust()
      vector[float] etot()
      vector[float] ec_ei()
      vector[float] ec_eo()
      vector[float] ec_t()
      vector[float] ec_r()
      vector[float] ech_x()
      vector[float] ech_y()
      vector[float] ech_z()
      vector[float] ec_m2()
      vector[float] ec_m3()
      vector[float] ec_m4()
      vector[float] ec_c2()
      vector[int] sc_sect()
      vector[int] sc_hit()
      vector[int] sc_pd()
      vector[int] sc_stat()
      vector[float] edep()
      vector[float] sc_t()
      vector[float] sc_r()
      vector[float] sc_c2()
      vector[int] cc_sect()
      vector[int] cc_hit()
      vector[int] cc_segm()
      vector[int] nphe()
      vector[float] cc_t()
      vector[float] cc_r()
      vector[float] cc_c2()

cdef char* str_to_char(str name):
  """Convert python string to char*"""
  cdef bytes name_bytes = name.encode()
  cdef char* c_name = name_bytes
  return c_name


cdef class h10:
  cdef int entry
  cdef TChain*c_chain
  cdef Branches*c_branches
  def __cinit__(h10 self, str branch_name, str file_name):
    self.entry = 0
    self.c_chain = new TChain(str_to_char(branch_name))
    self.c_chain.Add(str_to_char(file_name))
    self.c_branches = new Branches(self.c_chain)
  def add(self, file_name):
    self.c_chain.Add(str_to_char(file_name))
  @property
  def num_entries(self):
    return self.c_chain.GetEntries()
  def get_entry(self, num):
    self.c_chain.GetEntry(num)
  def __iter__(self):
      return self
  def next(self):
    if(self.entry <= self.c_chain.GetEntries()):
      self.c_chain.GetEntry(self.entry)
      self.entry += 1
      return self
    else:
      raise StopIteration
  def __next__(self):
    if(self.entry <= self.c_chain.GetEntries()):
      self.c_chain.GetEntry(self.entry)
      self.entry += 1
      return self
    else:
      raise StopIteration
  @property
  def gpart(self):
    return np.array(self.c_branches.gpart())
  @property
  def p(self):
    return np.array(self.c_branches.p())
  @property
  def b(self):
    return np.array(self.c_branches.b())
  @property
  def q(self):
    return np.array(self.c_branches.q())
  @property
  def cx(self):
    return np.array(self.c_branches.cx())
  @property
  def cy(self):
    return np.array(self.c_branches.cy())
  @property
  def cz(self):
    return np.array(self.c_branches.cz())
  @property
  def px(self):
    return np.array(self.c_branches.cx()) * np.array(self.c_branches.p())
  @property
  def py(self):
    return np.array(self.c_branches.cy()) * np.array(self.c_branches.p())
  @property
  def pz(self):
    return np.array(self.c_branches.cz()) * np.array(self.c_branches.p())
  @property
  def vx(self):
    return np.array(self.c_branches.vx())
  @property
  def vy(self):
    return np.array(self.c_branches.vy())
  @property
  def vz(self):
    return np.array(self.c_branches.vz())
  @property
  def id(self):
    return np.array(self.c_branches.id())
  @property
  def dc(self):
    return np.array(self.c_branches.dc())
  @property
  def cc(self):
    return np.array(self.c_branches.cc())
  @property
  def sc(self):
    return np.array(self.c_branches.sc())
  @property
  def ec(self):
    return np.array(self.c_branches.ec())
  @property
  def dc_sect(self):
    return np.array(self.c_branches.dc_sect())
  @property
  def dc_trk(self):
    return np.array(self.c_branches.dc_trk())
  @property
  def dc_stat(self):
    return np.array(self.c_branches.dc_stat())
  @property
  def dc_vx(self):
    return np.array(self.c_branches.dc_vx())
  @property
  def dc_vy(self):
    return np.array(self.c_branches.dc_vy())
  @property
  def dc_vz(self):
    return np.array(self.c_branches.dc_vz())
  @property
  def dc_vr(self):
    return np.array(self.c_branches.dc_vr())
  @property
  def dc_xsc(self):
    return np.array(self.c_branches.dc_xsc())
  @property
  def dc_ysc(self):
    return np.array(self.c_branches.dc_ysc())
  @property
  def dc_zsc(self):
    return np.array(self.c_branches.dc_zsc())
  @property
  def dc_cxsc(self):
    return np.array(self.c_branches.dc_cxsc())
  @property
  def dc_cysc(self):
    return np.array(self.c_branches.dc_cysc())
  @property
  def dc_czsc(self):
    return np.array(self.c_branches.dc_czsc())
  @property
  def dc_c2(self):
    return np.array(self.c_branches.dc_c2())
  @property
  def ec_stat(self):
    return np.array(self.c_branches.ec_stat())
  @property
  def ec_sect(self):
    return np.array(self.c_branches.ec_sect())
  @property
  def ec_whol(self):
    return np.array(self.c_branches.ec_whol())
  @property
  def ec_inst(self):
    return np.array(self.c_branches.ec_inst())
  @property
  def ec_oust(self):
    return np.array(self.c_branches.ec_oust())
  @property
  def etot(self):
    return np.array(self.c_branches.etot())
  @property
  def ec_ei(self):
    return np.array(self.c_branches.ec_ei())
  @property
  def ec_eo(self):
    return np.array(self.c_branches.ec_eo())
  @property
  def ec_t(self):
    return np.array(self.c_branches.ec_t())
  @property
  def ec_r(self):
    return np.array(self.c_branches.ec_r())
  @property
  def ech_x(self):
    return np.array(self.c_branches.ech_x())
  @property
  def ech_y(self):
    return np.array(self.c_branches.ech_y())
  @property
  def ech_z(self):
    return np.array(self.c_branches.ech_z())
  @property
  def ec_m2(self):
    return np.array(self.c_branches.ec_m2())
  @property
  def ec_m3(self):
    return np.array(self.c_branches.ec_m3())
  @property
  def ec_m4(self):
    return np.array(self.c_branches.ec_m4())
  @property
  def ec_c2(self):
    return np.array(self.c_branches.ec_c2())
  @property
  def sc_sect(self):
    return np.array(self.c_branches.sc_sect())
  @property
  def sc_hit(self):
    return np.array(self.c_branches.sc_hit())
  @property
  def sc_pd(self):
    return np.array(self.c_branches.sc_pd())
  @property
  def sc_stat(self):
    return np.array(self.c_branches.sc_stat())
  @property
  def edep(self):
    return np.array(self.c_branches.edep())
  @property
  def sc_t(self):
    return np.array(self.c_branches.sc_t())
  @property
  def sc_r(self):
    return np.array(self.c_branches.sc_r())
  @property
  def sc_c2(self):
    return np.array(self.c_branches.sc_c2())
  @property
  def cc_sect(self):
    return np.array(self.c_branches.cc_sect())
  @property
  def cc_hit(self):
    return np.array(self.c_branches.cc_hit())
  @property
  def cc_segm(self):
    return np.array(self.c_branches.cc_segm())
  @property
  def nphe(self):
    return np.array(self.c_branches.nphe())
  @property
  def cc_t(self):
    return np.array(self.c_branches.cc_t())
  @property
  def cc_r(self):
    return np.array(self.c_branches.cc_r())
  @property
  def cc_c2(self):
    return np.array(self.c_branches.cc_c2())
