from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from cython.operator cimport dereference as deref
from timeit import default_timer as timer
import glob
from physics_vectors import LorentzVector as LVector

cdef extern from "constants.hpp":
    cdef cppclass LorentzVector:
        float px()
        float py()
        float pz()
        float mass()

cdef extern from "mom_corr.cpp":
  pass

cdef extern from "mom_corr.hpp":
    cdef cppclass MomCorr:
        MomCorr()
        shared_ptr[LorentzVector] CorrectedVector(float, float, float, int)

cdef class MomentumCorrections:
    cdef:
        shared_ptr[MomCorr] c_corrections
        shared_ptr[LorentzVector] c_vec
    def __cinit__(self):
        self.c_corrections.reset(new MomCorr())
    def CorrectVector(self, float px, float py, float pz, int part_type):
        self.c_vec = deref(self.c_corrections).CorrectedVector(px, py, pz, part_type)
        cdef float _px = deref(self.c_vec).px()
        cdef float _py = deref(self.c_vec).py()
        cdef float _pz = deref(self.c_vec).pz()
        cdef float _mass = deref(self.c_vec).mass()

        return LVector(_px, _py, _pz, mass=_mass)