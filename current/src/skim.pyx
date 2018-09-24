from libcpp.string cimport string
from libcpp.vector cimport vector
from multiprocessing import Pool
import multiprocessing
from timeit import default_timer as timer
import glob


cdef extern from "skim.hpp":
    cdef cppclass Skim:
      Skim(vector[string], string) except +
      void Basic()
      void Strict()

cdef class py_skim:
  cdef Skim*c_skim
  cdef vector[string] f
  cdef string output
  def __cinit__(self, list_files, output):
    self.f = list_files
    self.output = output
    self.c_skim = new Skim(self.f, self.output)
    self.c_skim.Basic()

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
