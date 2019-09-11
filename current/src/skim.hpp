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

 public:
  Skim(const std::vector<std::string> &input, const std::string &output);
  ~Skim();

  float Basic();
  void Strict();
  void Final();
};
#endif
