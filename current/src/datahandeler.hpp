/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
#include <TFile.h>
#include <cstring>
#include <fstream>
#include <future>
#include <iostream>
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
#include "mom_corr.hpp"
#include "photon_flux.hpp"
#include "physics.hpp"
#include "reaction.hpp"

class DataHandeler {
 protected:
  bool _loadbar = false;
  std::shared_ptr<TChain> _chain;
  std::vector<std::string> _input_files;
  std::shared_ptr<Histogram> _hists;
  std::shared_ptr<Branches> _data;
  std::shared_ptr<BranchesMC> _data_mc;
  std::shared_ptr<MomCorr> _mom_corr = nullptr;
  bool CUTS;

 public:
  DataHandeler();
  DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists,
               const std::shared_ptr<MomCorr>& mom_corr);

  DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists);
  ~DataHandeler();
  void setLoadBar(bool load);

  int Run();
  const void RunEvent(size_t current_event);

  void loadbar(long x, long n);
};

class mcHandeler : public DataHandeler {
 protected:
  std::shared_ptr<mcHistogram> _mc_hists;

 public:
  mcHandeler() : DataHandeler() {}
  mcHandeler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists);
  int Run();
  const void RunEvent(int current_event);
};

#endif
