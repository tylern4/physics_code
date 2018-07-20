from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "../src/skim.hpp" namespace "Skim":
    cdef cppclass Skim:
      Skim(string) except +
      Process()


cdef class py_skim:
  cdef Skim*c_skim
  def __cinit__(self, filename):
    self.c_skim = new Skim(filename)
  #def Process():
  #  self.c_skim.Process()
