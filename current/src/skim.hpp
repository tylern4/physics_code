/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#ifndef SKIM_H_GUARD
#define SKIM_H_GUARD
#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "branches.hpp"
#include "color.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "delta_t.hpp"
#include "func.hpp"
#include "glob_files.hpp"
#include "missing_mass.hpp"
#include "reaction.hpp"

class Skim {
 private:
  std::shared_ptr<TChain> chain;
  std::string fout;
  std::vector<std::string> fin;
  std::shared_ptr<TFile> RootOutputFile;
  std::shared_ptr<TLorentzVector> e_mu;
  std::shared_ptr<MissingMass> MM_neutron;
  double BEAM_ENERGY = E1D_E0;

 public:
  Skim(const std::vector<std::string> &input, const std::string &output);
  ~Skim();

  float Basic();
  void Strict();
  void Final();

  double sf_top_fit(double P);
  double sf_bot_fit(double P);
  bool sf_cut(double sf, double P);

  double dt_P_bot_fit(double P);
  double dt_P_top_fit(double P);
  bool dt_P_cut(double dt, double P);

  double dt_Pip_bot_fit(double P);
  double dt_Pip_top_fit(double P);
  bool dt_Pip_cut(double dt, double P);
};
#endif
