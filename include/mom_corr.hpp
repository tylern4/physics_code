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
  std::mutex _readerMutex;
  char *_datadir;
  /* Theta Binning for Theta correction */
#define ThetaC_n 144
  float ThetaC_min = 0;
  float ThetaC_max = 144;
  float ThetaC_wid = 1.;

/* Theta Binning for Momentum correction */
#define MomC_T_n 48
  float MomC_T_min = 1;
  float MomC_T_max = 145;
  float MomC_T_wid = 3.;

  float c0_theta[ThetaC_n][NUM_SECTORS];
  float c1_theta[ThetaC_n][NUM_SECTORS];
  float c2_theta[ThetaC_n][NUM_SECTORS];

  float c0_mom[MomC_T_n][NUM_SECTORS][Npar];
  float c1_mom[MomC_T_n][NUM_SECTORS][Npar];

  float d0_mom[MomC_T_n][NUM_SECTORS][Npar];
  float d1_mom[MomC_T_n][NUM_SECTORS][Npar];

  /* REad angle correction parameters */
  void read_theta_par();
  /* Read momentum correction parameters for electrons */
  void read_mom_par();
  /* Read momentum correction parameters for pi+ */
  void read_mom_pip_par();

  /* Angle correction */
  float theta_corr(float, float, Int_t);

  /* momentum correction for electrons */
  float mom_corr(float, float, float, Int_t);

  /* momentum correction for hadrons */
  float mom_corr_pip(float, float, float, Int_t);

  Int_t GetSector(float phi);

 public:
  MomCorr();
  std::shared_ptr<LorentzVector> CorrectedVector(float px, float py, float pz, int particle_type);
};

#endif