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
#include "color.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "delta_t.hpp"
#include "histogram.hpp"
#include "missing_mass.hpp"
#include "physics.hpp"

class DataHandeler {
 private:
  int PID;
  Histogram *hists;
  TFile *RootOutputFile;
  TLorentzVector *e_mu;
  std::vector<std::string> input_files;
  MissingMass *MM_neutron;
  std::vector<bool> *pip_vec;
  std::vector<bool> *pim_vec;
  std::vector<bool> *proton_vec;
  std::vector<bool> *elec_vec;
  TCanvas *c1;

 public:
  DataHandeler(std::vector<std::string> fin, std::string RootFile_output);
  ~DataHandeler();
  void file_handeler(std::string fin);
  void loadbar(long x, long n);
  void run();
};

#endif
