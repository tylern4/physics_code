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
#include "Event.hpp"
#include "Particle.hpp"
#include "TChain.h"
#include "color.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "delta_t.hpp"
#include "histogram.hpp"
#include "missing_mass.hpp"
#include "photon_flux.hpp"
#include "physics.hpp"

class DataHandeler {
 private:
  int PID;
  std::ofstream csv_output;
  std::vector<std::string> input_files;
  std::vector<std::vector<TLorentzVector>> All_events;
  double BEAM_ENERGY = E1D_E0;

 public:
  DataHandeler();
  ~DataHandeler();
  void Run(std::string fin, Histogram *hists);
  void Run(std::vector<std::string> fin, Histogram *hists);
  void loadbar(long x, long n);
};

#endif
