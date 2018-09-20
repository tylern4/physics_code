# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cpython cimport array
import numpy as np
cimport numpy as np
import array

cdef extern from "TChain.h":
  cdef cppclass TChain:
    TChain(char*) except +
    int Add(char*)
    long GetEntries()
    int GetEntry(long)

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
    return self.c_branches.gpart()
  @property
  def p(self):
    return self.c_branches.p()
  @property
  def b(self):
    return self.c_branches.b()
  @property
  def q(self):
    return self.c_branches.q()
  @property
  def cx(self):
    return self.c_branches.cx()
  @property
  def cy(self):
    return self.c_branches.cy()
  @property
  def cz(self):
    return self.c_branches.cz()
  @property
  def vx(self):
    return self.c_branches.vx()
  @property
  def vy(self):
    return self.c_branches.vy()
  @property
  def vz(self):
    return self.c_branches.vz()
  @property
  def cx(self):
    return self.c_branches.cx()
  @property
  def id(self):
    return self.c_branches.id()
  @property
  def dc(self):
    return self.c_branches.dc()
  @property
  def cc(self):
    return self.c_branches.cc()
  @property
  def sc(self):
    return self.c_branches.sc()
  @property
  def ec(self):
    return self.c_branches.ec()


cdef extern from "delta_t.cpp":
  pass

cdef extern from "delta_t.hpp":
    cdef cppclass Delta_T:
      Delta_T(double sc_time, double sc_pathlength) except +
      double deltat(double momentum, double sc_t, double sc_r, double mass)

cdef class Dt:
  cdef Delta_T*c_Dt
  def __cinit__(Dt self, double sc_time, double sc_pathlength):
    self.c_Dt = new Delta_T(sc_time, sc_pathlength)
  def deltat(Dt self, double momentum, double sc_t, double sc_r, double mass):
    return self.c_Dt.deltat(momentum, sc_t, sc_r, mass)
