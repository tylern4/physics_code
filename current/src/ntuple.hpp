/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef NTUPLE_H_GUARD
#define NTUPLE_H_GUARD
#include <future>
#include <iostream>
#include <thread>
#include <vector>
#include "constants.hpp"
#include "cuts.hpp"
#include "delta_t.hpp"
#include "physics.hpp"
#include "reaction.hpp"

class Ntuple {
 private:
  Int_t _type;
  Int_t _helicity;
  Float_t _W;
  Float_t _Q2;
  Float_t _MM;
  Float_t _MM2;
  Float_t _Theta_E;
  Float_t _Theta_star;
  Float_t _Phi_star;
  Float_t _theta;
  Float_t _phi;
  Int_t _sector;

  std::shared_ptr<TTree> ntuple = nullptr;
  std::shared_ptr<TFile> rootout = nullptr;

 public:
  Ntuple(const std::string &output_file_name);
  Ntuple(const std::vector<std::string> &infiles, const std::string &output_file_name);
  ~Ntuple();
  size_t Run(const std::shared_ptr<TChain> &chain);
  size_t Run();
  size_t run_file(std::vector<std::string> in, int thread_id);
};

#endif
