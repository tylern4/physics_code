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
#include <string>
#include <thread>
#include <vector>
#include "TChain.h"
#include "branches.hpp"
#include "color.hpp"
#include "constants.hpp"
#include "histogram_mc.hpp"
#include "missing_mass.hpp"
#include "physics.hpp"

class mcHandeler {
 private:
  int PID;
  std::vector<std::string> input_files;
  double BEAM_ENERGY = E1D_E0;
  bool CUTS = true;
  TLorentzVector *e_mu;

 public:
  mcHandeler();
  ~mcHandeler();
  void Run(std::string fin, mcHistogram *hists);
  void Run(std::vector<std::string> fin, mcHistogram *hists);
  void loadbar(long x, long n);
};

#endif
