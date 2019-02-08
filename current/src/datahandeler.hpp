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

class DataHandeler {
 private:
  int PID;
  std::vector<std::string> input_files;
  bool CUTS = true;

 public:
  DataHandeler();
  ~DataHandeler();
  void Run(std::string fin, Histogram *hists);
  void Run(std::vector<std::string> fin, Histogram *hists);
  void loadbar(long x, long n);
};

#endif
