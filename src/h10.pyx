# distutils: language = c++
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: infertypes=True
#cython: initializedcheck=False
#cython: cdivision=True
#cython: embedsignature=True
#cython: profile=True

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr, shared_ptr
from cython.operator cimport dereference as deref
import cython
from libcpp cimport bool
from cpython cimport array
from libc.stdlib cimport free
import numpy as np
cimport numpy as np

class colors:
  RESET = b"\033[0m"
  BLACK = b"\033[30m"
  RED = b"\033[31m"
  GREEN = b"\033[32m"
  YELLOW = b"\033[33m"
  BLUE = b"\033[34m"
  MAGENTA = b"\033[35m"
  CYAN = b"\033[36m"
  WHITE = b"\033[37m"
  BOLDBLACK = b"\033[1m\033[30m"
  BOLDRED = b"\033[1m\033[31m"
  BOLDGREEN = b"\033[1m\033[32m"
  BOLDYELLOW = b"\033[1m\033[33m"
  BOLDBLUE = b"\033[1m\033[34m"
  BOLDMAGENTA = b"\033[1m\033[35m"
  BOLDCYAN = b"\033[1m\033[36m"
  BOLDWHITE = b"\033[1m\033[37m"

cdef extern from "TChain.h":
  cdef cppclass TChain:
    TChain(char*) except +
    int Add(char*)
    long GetEntries()
    int GetEntry(long)

cdef extern from "branches.cpp":
  pass

cdef extern from "branches.hpp":
    cdef cppclass Branches:
      Branches(shared_ptr[TChain]) except +
      Branches(shared_ptr[TChain], bool) except +
      int npart()
      int evstat()
      int intt()
      int evntid()
      int evtype()
      int evntclas()
      int evthel()
      int evntclas2()
      float q_l()
      float t_l()
      float tr_time()
      float rf_time1()
      float rf_time2()
      int gpart()
      vector[int] id()
      vector[int] dc()
      vector[int] cc()
      vector[int] sc()
      vector[int] ec()
      vector[float] p()
      vector[int] q()
      vector[float] b()
      vector[float] cx()
      vector[float] cy()
      vector[float] cz()
      vector[float] vx()
      vector[float] vy()
      vector[float] vz()
      vector[int] dc_sect()
      vector[int] dc_trk()
      vector[int] dc_stat()
      vector[float] dc_vx()
      vector[float] dc_vy()
      vector[float] dc_vz()
      vector[float] dc_vr()
      vector[float] dc_xsc()
      vector[float] dc_ysc()
      vector[float] dc_zsc()
      vector[float] dc_cxsc()
      vector[float] dc_cysc()
      vector[float] dc_czsc()
      vector[float] dc_c2()
      vector[int] ec_stat()
      vector[int] ec_sect()
      vector[int] ec_whol()
      vector[int] ec_inst()
      vector[int] ec_oust()
      vector[float] etot()
      vector[float] ec_ei()
      vector[float] ec_eo()
      vector[float] ec_t()
      vector[float] ec_r()
      vector[float] ech_x()
      vector[float] ech_y()
      vector[float] ech_z()
      vector[float] ec_m2()
      vector[float] ec_m3()
      vector[float] ec_m4()
      vector[float] ec_c2()
      vector[int] sc_sect()
      vector[int] sc_hit()
      vector[int] sc_pd()
      vector[int] sc_stat()
      vector[float] edep()
      vector[float] sc_t()
      vector[float] sc_r()
      vector[float] sc_c2()
      vector[int] cc_sect()
      vector[int] cc_hit()
      vector[int] cc_segm()
      vector[int] nphe()
      vector[float] cc_t()
      vector[float] cc_r()
      vector[float] cc_c2()
      ####################
      vector[int] pidpart()
      vector[float] xpart()
      vector[float] ypart()
      vector[float] zpart()
      vector[float] epart()
      vector[float] pxpart()
      vector[float] pypart()
      vector[float] pzpart()
      vector[float] qpart()


cdef class h10_data:
  cdef:
    int entry
    shared_ptr[TChain] c_chain
    shared_ptr[Branches] c_branches
    int total_num_entries
  def __cinit__(h10_data self):
    self.entry = 0
    self.c_chain.reset(new TChain("h10"))
    self.c_branches.reset(new Branches(self.c_chain))
  def add(self, file_name):
    deref(self.c_chain).Add(file_name.encode())
  @property
  def num_entries(self):
    if not self.total_num_entries:
      self.total_num_entries = deref(self.c_chain).GetEntries()
    return self.total_num_entries
  @property
  def entry(self):
    return self.entry - 1
  def get_entry(self, num):
    deref(self.c_chain).GetEntry(num)
  def reset(self):
      self.entry = 0
  def __iter__(self):
      return self
  def next(self):
    if(self.entry <= self.total_num_entries):
      deref(self.c_chain).GetEntry(self.entry)
      self.entry += 1
      return self
    else:
      raise StopIteration
  def __next__(self):
    if(self.entry <= self.total_num_entries):
      deref(self.c_chain).GetEntry(self.entry)
      self.entry += 1
      return self
    else:
      raise StopIteration
  def __len__(self):
    return deref(self.c_branches).gpart()
  @property
  def gpart(self):
    return deref(self.c_branches).gpart()
  @property
  def p(self):
    return np.array(deref(self.c_branches).p())
  @property
  def b(self):
    return np.array(deref(self.c_branches).b())
  @property
  def q(self):
    return np.array(deref(self.c_branches).q())
  @property
  def cx(self):
    return np.array(deref(self.c_branches).cx())
  @property
  def cy(self):
    return np.array(deref(self.c_branches).cy())
  @property
  def cz(self):
    return np.array(deref(self.c_branches).cz())
  @property
  def px(self):
    return np.array(deref(self.c_branches).cx()) * np.array(deref(self.c_branches).p())
  @property
  def py(self):
    return np.array(deref(self.c_branches).cy()) * np.array(deref(self.c_branches).p())
  @property
  def pz(self):
    return np.array(deref(self.c_branches).cz()) * np.array(deref(self.c_branches).p())
  @property
  def vx(self):
    return np.array(deref(self.c_branches).vx())
  @property
  def vy(self):
    return np.array(deref(self.c_branches).vy())
  @property
  def vz(self):
    return np.array(deref(self.c_branches).vz())
  @property
  def id(self):
    return np.array(deref(self.c_branches).id())
  @property
  def dc(self):
    return np.array(deref(self.c_branches).dc())
  @property
  def cc(self):
    return np.array(deref(self.c_branches).cc())
  @property
  def sc(self):
    return np.array(deref(self.c_branches).sc())
  @property
  def ec(self):
    return np.array(deref(self.c_branches).ec())
  @property
  def dc_sect(self):
    return np.array(deref(self.c_branches).dc_sect())
  @property
  def dc_trk(self):
    return np.array(deref(self.c_branches).dc_trk())
  @property
  def dc_stat(self):
    return np.array(deref(self.c_branches).dc_stat())
  @property
  def dc_vx(self):
    return np.array(deref(self.c_branches).dc_vx())
  @property
  def dc_vy(self):
    return np.array(deref(self.c_branches).dc_vy())
  @property
  def dc_vz(self):
    return np.array(deref(self.c_branches).dc_vz())
  @property
  def dc_vr(self):
    return np.array(deref(self.c_branches).dc_vr())
  @property
  def dc_xsc(self):
    return np.array(deref(self.c_branches).dc_xsc())
  @property
  def dc_ysc(self):
    return np.array(deref(self.c_branches).dc_ysc())
  @property
  def dc_zsc(self):
    return np.array(deref(self.c_branches).dc_zsc())
  @property
  def dc_cxsc(self):
    return np.array(deref(self.c_branches).dc_cxsc())
  @property
  def dc_cysc(self):
    return np.array(deref(self.c_branches).dc_cysc())
  @property
  def dc_czsc(self):
    return np.array(deref(self.c_branches).dc_czsc())
  @property
  def dc_c2(self):
    return np.array(deref(self.c_branches).dc_c2())
  @property
  def ec_stat(self):
    return np.array(deref(self.c_branches).ec_stat())
  @property
  def ec_sect(self):
    return np.array(deref(self.c_branches).ec_sect())
  @property
  def ec_whol(self):
    return np.array(deref(self.c_branches).ec_whol())
  @property
  def ec_inst(self):
    return np.array(deref(self.c_branches).ec_inst())
  @property
  def ec_oust(self):
    return np.array(deref(self.c_branches).ec_oust())
  @property
  def etot(self):
    return np.array(deref(self.c_branches).etot())
  @property
  def ec_ei(self):
    return np.array(deref(self.c_branches).ec_ei())
  @property
  def ec_eo(self):
    return np.array(deref(self.c_branches).ec_eo())
  @property
  def ec_t(self):
    return np.array(deref(self.c_branches).ec_t())
  @property
  def ec_r(self):
    return np.array(deref(self.c_branches).ec_r())
  @property
  def ech_x(self):
    return np.array(deref(self.c_branches).ech_x())
  @property
  def ech_y(self):
    return np.array(deref(self.c_branches).ech_y())
  @property
  def ech_z(self):
    return np.array(deref(self.c_branches).ech_z())
  @property
  def ec_m2(self):
    return np.array(deref(self.c_branches).ec_m2())
  @property
  def ec_m3(self):
    return np.array(deref(self.c_branches).ec_m3())
  @property
  def ec_m4(self):
    return np.array(deref(self.c_branches).ec_m4())
  @property
  def ec_c2(self):
    return np.array(deref(self.c_branches).ec_c2())
  @property
  def sc_sect(self):
    return np.array(deref(self.c_branches).sc_sect())
  @property
  def sc_hit(self):
    return np.array(deref(self.c_branches).sc_hit())
  @property
  def sc_pd(self):
    return np.array(deref(self.c_branches).sc_pd())
  @property
  def sc_stat(self):
    return np.array(deref(self.c_branches).sc_stat())
  @property
  def edep(self):
    return np.array(deref(self.c_branches).edep())
  @property
  def sc_t(self):
    return np.array(deref(self.c_branches).sc_t())
  @property
  def sc_r(self):
    return np.array(deref(self.c_branches).sc_r())
  @property
  def sc_c2(self):
    return np.array(deref(self.c_branches).sc_c2())
  @property
  def cc_sect(self):
    return np.array(deref(self.c_branches).cc_sect())
  @property
  def cc_hit(self):
    return np.array(deref(self.c_branches).cc_hit())
  @property
  def cc_segm(self):
    return np.array(deref(self.c_branches).cc_segm())
  @property
  def nphe(self):
    return np.array(deref(self.c_branches).nphe())
  @property
  def cc_t(self):
    return np.array(deref(self.c_branches).cc_t())
  @property
  def cc_r(self):
    return np.array(deref(self.c_branches).cc_r())
  @property
  def cc_c2(self):
    return np.array(deref(self.c_branches).cc_c2())  

cdef class h10_mc(h10_data):
  def __cinit__(self):
    self.entry = 0
    self.c_chain.reset(new TChain("h10"))
    self.c_branches.reset(new Branches(self.c_chain, True))
  @property
  def pidpart(self):
   return np.array(deref(self.c_branches).pidpart())
  @property
  def xpart(self):
   return np.array(deref(self.c_branches).xpart())
  @property
  def ypart(self):
   return np.array(deref(self.c_branches).ypart())
  @property
  def zpart(self):
   return np.array(deref(self.c_branches).zpart())
  @property
  def epart(self):
   return np.array(deref(self.c_branches).epart())
  @property
  def pxpart(self):
   return np.array(deref(self.c_branches).pxpart())
  @property
  def pypart(self):
   return np.array(deref(self.c_branches).pypart())
  @property
  def pzpart(self):
   return np.array(deref(self.c_branches).pzpart())
  @property
  def qpart(self):
   return np.array(deref(self.c_branches).qpart())

cdef extern from "reaction.cpp":
  pass

cdef extern from "reaction.hpp":
    cdef cppclass Reaction:
      Reaction(shared_ptr[Branches]) except +
      Reaction(shared_ptr[Branches],float) except +
      void SetProton(int)
      void SetPip(int)
      void SetPim(int)
      void SetOther(int)
      void SetNeutron(int)

      bool Reset()
      bool SingleP()

      float W()
      float Q2()

      float MM()
      float MM2()
      float pi0_mass()
      float pi0_mass2()

cdef class reaction:
  cdef:
    shared_ptr[Branches] c_branches
    shared_ptr[Reaction] c_reaction
  def __cinit__(self, h10_data branch, float beam_E):
    self.c_branches = branch.c_branches
    self.c_reaction.reset(new Reaction(self.c_branches, beam_E))
  def __cinit__(self, h10_data branch):
    self.c_branches = branch.c_branches
    self.c_reaction.reset(new Reaction(self.c_branches))
  def reset(self):
    deref(self.c_reaction).Reset()
  def proton(self, int i):
   deref(self.c_reaction).SetProton(i)
  def neutron(self, int i):
   deref(self.c_reaction).SetNeutron(i)
  def pip(self, int i):
   deref(self.c_reaction).SetPip(i)
  def pim(self, int i):
   deref(self.c_reaction).SetPim(i)
  def other(self, int i):
   deref(self.c_reaction).SetOther(i)
  def pi0_M(self):
    return deref(self.c_reaction).pi0_mass()
  def pi0_M2(self):
    return deref(self.c_reaction).pi0_mass2()
  @property
  def w(self):
    return deref(self.c_reaction).W()
  @property
  def q2(self):
    return deref(self.c_reaction).Q2()
  def single_P(self):
    return deref(self.c_reaction).SingleP()

