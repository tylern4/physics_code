# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector

import glob


cdef extern from "yeilds.hpp":
    cdef cppclass Yeilds:
      Yeilds() except +
      Yeilds(string) except +
      void WriteHeader()
      void Run[T](vector[string])
      void OpenFile(string)
      #int RunNtuple(string, bool)

cdef class yeild:
  cdef Yeilds*c_yeilds
  cdef vector[string] _files
  cdef string output
  def __cinit__(self, list_files, output):
    self._files = glob.glob(list_files)
    self.output = output
    self.c_yeilds = new Yeilds()
    self.c_yeilds.OpenFile(self.output)
  def header(self):
    self.c_yeilds.WriteHeader()
  #def run(self):
  #  self.c_yeilds.Run[CutType](self._files)
  #def ntuple(self):
  #  for f in self._files:
  #   self.c_yeilds.RunNtuple(f)
