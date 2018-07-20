/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#ifndef SKIM_H_GUARD
#define SKIM_H_GUARD
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "color.hpp"
#include "constants.hpp"
#include "glob_files.hpp"
#include <iostream>

class Skim {
private:
  TChain *chain;
  std::string fout;
  std::string fin;
  TFile *RootOutputFile;
  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;
  TLorentzVector e_mu;
  TVector3 Particle3;
  TLorentzVector Particle4;
  int gpart;
  int id[MAX_PARTS];      //[gpart]
  int stat[MAX_PARTS];    //[gpart]
  int dc[MAX_PARTS];      //[gpart]
  int cc[MAX_PARTS];      //[gpart]
  int sc[MAX_PARTS];      //[gpart]
  int ec[MAX_PARTS];      //[gpart]
  float p[MAX_PARTS];     //[gpart]
  float m[MAX_PARTS];     //[gpart]
  int q[MAX_PARTS];       //[gpart]
  float b[MAX_PARTS];     //[gpart]
  float cx[MAX_PARTS];    //[gpart]
  float cy[MAX_PARTS];    //[gpart]
  float cz[MAX_PARTS];    //[gpart]
  float vx[MAX_PARTS];    //[gpart]
  float vy[MAX_PARTS];    //[gpart]
  float vz[MAX_PARTS];    //[gpart]
  int dc_stat[MAX_PARTS]; //[dc_part]
  float etot[MAX_PARTS];  //[ec_part]

  void getSkimBranches();

public:
  Skim(std::string input);
  void Process();
  ~Skim();
};
#endif
