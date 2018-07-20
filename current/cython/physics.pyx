# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cpython cimport array
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
    cdef void getBranches(TChain*)

cdef class h10:
  cdef entry
  cdef TChain*c_chain
  def __cinit__(self, name):
    self.entry = 0
    cdef bytes name_bytes = name.encode()
    cdef char* c_name = name_bytes
    self.c_chain = new TChain(c_name)
  def add(self, file):
    cdef bytes file_bytes = file.encode()
    cdef char* c_file = file_bytes
    self.c_chain.Add(c_file)
  def num_entris(self):
    return self.c_chain.GetEntries()
  def get_entry(self, num):
    self.c_chain.GetEntry(num)
  def get_next(self):
    if(self.entry <= self.c_chain.GetEntries()):
      self.c_chain.GetEntry(self.entry)
      self.entry += 1
      return True
    else:
      return False
  def get_branches(self):
    getBranches(self.c_chain)
  def npart(self):
    return self.npart
  def get_p(self):
    return 
