/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
#include <TFile.h>
#include <TLorentzVector.h>
#include <cstring>
#include <fstream>
#include <future>
#include <string>
#include <thread>
#include <vector>
#include "TChain.h"
#include "color.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "delta_t.hpp"
#include "histogram.hpp"
#include "missing_mass.hpp"
#include "photon_flux.hpp"
#include "physics.hpp"
#include "reaction.hpp"

#if defined(ENABLE_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
#endif

class DataHandeler {
 protected:
  bool _loadbar = false;
  TChain* _chain[4];
  std::vector<std::string> _input_files;
  std::shared_ptr<Histogram> _hists;
  std::shared_ptr<Branches> _data[4];
  bool CUTS = true;

 public:
  DataHandeler();
  DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists);
  ~DataHandeler();
  void setLoadBar(bool load);
  // int Run(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists);
  int Run();
  int Run(int current_event, int thread);
  void loadbar(long x, long n);
};

class mcHandeler : public DataHandeler {
 protected:
  std::shared_ptr<mcHistogram> _mc_hists;

 public:
  mcHandeler() : DataHandeler() {}
  mcHandeler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists);
  int Run();
  int Run(int current_event, int thread);
};

#endif
