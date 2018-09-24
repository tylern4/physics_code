# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cpython cimport array
import array
from cython.parallel import prange
from multiprocessing import Pool

ctypedef vector[DataHandeler*] dh_vec
ctypedef vector[Histogram] hist_vec

cdef extern from "histogram.cpp":
  pass

cdef extern from "histogram.hpp":
  cdef cppclass Histogram:
    Histogram() except +
    Histogram(string) except +
    void Write(string)

cdef extern from "datahandeler.cpp":
  pass

cdef extern from "datahandeler.hpp":
  cdef cppclass DataHandeler:
    DataHandeler() except +
    void Run(string, Histogram*)

cdef class handeler:
  cdef dh_vec c_handeler
  cdef Histogram*hists
  cdef list file_names
  cdef string output_file
  def __cinit__(self, list file_names, string output_file):
    self.file_names = file_names
    self.output_file = output_file
    self.c_handeler = dh_vec(len(self.file_names))
    self.hists = new Histogram()
    for i in range(0, len(self.file_names)):
      self.c_handeler[i] = new DataHandeler()
  def run(self):
    cdef i = 0
    cdef l = len(self.file_names)
    for i in range(l):
      self.c_handeler[i].Run(self.file_names[i], self.hists)
    self.hists.Write(self.output_file)
