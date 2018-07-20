from libcpp.string cimport string
from libcpp.vector cimport vector
from multiprocessing import Pool
import multiprocessing
from timeit import default_timer as timer
import glob


cdef extern from "skim.hpp":
    cdef cppclass Skim:
      Skim(string) except +
      void Process()

cdef class py_skim:
  cdef Skim*c_skim
  cdef string f
  def __cinit__(self, filename):
    self.f = filename
    self.c_skim = new Skim(filename)
    self.c_skim.Process()
  def __reduce__(self):
    return (self.__class__, (self.f, ))

class skim_files:
  def __init__(self, input):
    self.files = glob.glob(input)
  def run(self):
    start = timer()
    num_cores = multiprocessing.cpu_count()
    pool = Pool(processes=num_cores)
    pool.map(py_skim, (self.files))
    end = timer()
    print("Time: "+ str(end - start))
