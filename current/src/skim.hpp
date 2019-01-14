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
  TChain *chain;
  std::shared_ptr<Branches> data;
  std::string fout;
  std::vector<std::string> fin;
  TFile *RootOutputFile;
  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;
  TVector3 particle_3;
  TLorentzVector particle;
  TLorentzVector e_mu;
  TVector3 Particle3;
  TLorentzVector Particle4;
  MissingMass *MM_neutron;
  double BEAM_ENERGY = E1D_E0;

 public:
  Skim(std::vector<std::string> input, std::string output);
  void Basic();
  void Strict();
  void Final();
  ~Skim();
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
