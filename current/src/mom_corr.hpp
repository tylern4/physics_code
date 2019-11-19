#ifndef MOMCORR_H
#define MOMCORR_H
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "branches.hpp"
#include "constants.hpp"

#define Npar 4

class MomCorr {
 private:
  char *_datadir;
  /* Theta Binning for Theta correction */
#define ThetaC_n 144
  Double_t ThetaC_min = 0;
  Double_t ThetaC_max = 144;
  Double_t ThetaC_wid = 1.;

/* Theta Binning for Momentum correction */
#define MomC_T_n 48
  Double_t MomC_T_min = 1;
  Double_t MomC_T_max = 145;
  Double_t MomC_T_wid = 3.;

  Double_t c0_theta[ThetaC_n][NUM_SECTORS];
  Double_t c1_theta[ThetaC_n][NUM_SECTORS];
  Double_t c2_theta[ThetaC_n][NUM_SECTORS];

  Double_t c0_mom[MomC_T_n][NUM_SECTORS][Npar];
  Double_t c1_mom[MomC_T_n][NUM_SECTORS][Npar];

  Double_t d0_mom[MomC_T_n][NUM_SECTORS][Npar];
  Double_t d1_mom[MomC_T_n][NUM_SECTORS][Npar];

  /* REad angle correction parameters */
  void read_theta_par();
  /* Read momentum correction parameters for electrons */
  void read_mom_par();
  /* Read momentum correction parameters for pi+ */
  void read_mom_pip_par();

  /* Angle correction */
  Double_t theta_corr(Double_t, Double_t, Int_t);

  /* momentum correction for electrons */
  Double_t mom_corr(Double_t, Double_t, Double_t, Int_t);

  /* momentum correction for hadrons */
  Double_t mom_corr_pip(Double_t, Double_t, Double_t, Int_t);

  Int_t GetSector(Double_t phi);

 public:
  MomCorr();
  std::shared_ptr<LorentzVector> CorrectedVector(float px, float py, float pz, int particle_type);
};

#endif