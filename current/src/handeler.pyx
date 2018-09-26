# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cpython cimport array
import array
from cython.parallel import prange
from multiprocessing import Pool


class colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

cdef extern from "histogram.cpp":
  pass

cdef extern from "histogram.hpp":
  cdef cppclass Histogram:
    Histogram() except +
    void Write(string)

cdef extern from "datahandeler.cpp":
  pass

cdef extern from "datahandeler.hpp":
  cdef cppclass DataHandeler:
    DataHandeler() except +
    void Run(string, Histogram*)

cdef class handeler:
  cdef DataHandeler*c_handeler
  cdef Histogram*hists
  cdef list file_names
  cdef string output_file, input_file
  def __cinit__(self):
    self.c_handeler = new DataHandeler()
    self.hists = new Histogram()
  def run(self, string input_file):
    self.input_file = input_file
    self.c_handeler.Run(self.input_file, self.hists)
  def write(self, string output_file):
    self.output_file = output_file
    self.hists.Write(self.output_file)
