/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef Datahandeler_multi_H_GUARD
#define Datahandeler_multi_H_GUARD
#include <TFile.h>
#include <TLorentzVector.h>
#include <cstring>
#include <fstream>
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

// using vec4 = ROOT::Math::PxPyPzMVector;

class Datahandeler_multi {
 private:
  int PID;
  Histogram *hists;
  TFile *RootOutputFile;
  TLorentzVector *e_mu;
  TH1D *mm = new TH1D("mm", "mm", 500, 0.0, 10.0);
  std::vector<std::string> input_files;
  std::vector<std::vector<TLorentzVector>> All_events;
  MissingMass *MM_neutron;
  TCanvas *c1;

 public:
  Datahandeler_multi(std::vector<std::string> fin, std::string RootFile_output);
  ~Datahandeler_multi();
  void load_files(std::string fin);
  void loadbar(long x, long n);
  void run();
};

#endif
