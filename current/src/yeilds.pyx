from libcpp.string cimport string
from libcpp.vector cimport vector

from timeit import default_timer as timer
import glob


cdef extern from "yeilds.hpp":
    cdef cppclass Yeilds:
      Yeilds() except +
      Yeilds(string) except +
      void WriteHeader()
      void Run(vector[string])
      void OpenFile(string)

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
  def run(self):
    self.c_yeilds.Run(self._files)
