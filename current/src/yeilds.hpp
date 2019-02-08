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

class Yeilds {
 private:
  int PID;
  std::ofstream csv_output;
  std::vector<std::string> input_files;
  bool CUTS = true;

 public:
  Yeilds();
  Yeilds(std::string output_file_name);
  ~Yeilds();
  void OpenFile(std::string output_file_name);
  void WriteHeader();
  int Run(std::string fin);
  void Run(std::vector<std::string> fin);
};

#endif
