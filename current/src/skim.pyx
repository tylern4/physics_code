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
      void Final()

cdef class py_skim:
  cdef Skim*c_skim
  cdef vector[string] f
  cdef string output
  def __cinit__(self, list_files, output):
    self.f = list_files
    self.output = output
  def __cinit__(self, file):
    self.f = [file]
    self.output = file[:-5]+b'_skim.root'
    self.basic()
  def __reduce__(self):
    return (self.output.decode("utf-8"))
  def basic(self):
    self.c_skim = new Skim(self.f, self.output)
    self.c_skim.Basic()
  def strict(self):
    self.c_skim = new Skim(self.f, self.output)
    self.c_skim.Strict()
  def final(self):
    self.c_skim = new Skim(self.f, self.output)
    self.c_skim.Final()

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
