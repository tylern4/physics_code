/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef FITS_HPP
#define FITS_HPP
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "physics.hpp"
#include "Math/MinimizerOptions.h"

class Fits {
 private:
  double par_max, par_mean, par_FWHM;
  double mean, sigma, FWHM;
  double min_value, max_value;
  double a, b, c, d, e = 0.0;
  double min_edge_x = std::nan("0"), max_edge_x = std::nan("0");
  std::string ploy1d = "[0]+[1]*x";
  std::string ploy2d = ploy1d + "+[2]*x*x";
  std::string ploy3d = ploy2d + "+[3]*x*x*x";
  std::string ploy4d = ploy3d + "+[4]*x*x*x*x";
  std::string gaus = "[0]*TMath::Gaus(x,[1],[2],1)";
  double fiducial_phi_lo(double theta_e, double theta_e_min, double k, double m, int c);
  double fiducial_phi_hi(double theta_e, double theta_e_min, double k, double m, int c);

 public:
  Fits();
  ~Fits();
  void FitGaus(TH1D *hist);
  void FitLandau(TH1D *hist);
  void FitPoly_1D(TH1D *hist);
  void FitPoly_2D(TH1D *hist);
  void FitPoly_3D(TH1D *hist);
  void FitPoly_4D(TH1D *hist);
  void FitPoly_fid(TH2D *hist);
  void FitFiducial_lo(TH2D *hist2d);
  void FitFiducial_hi(TH2D *hist2d);
  void FitFiducial(TH2D *hist2d);
  void FitGenNormal(TH1D *hist);
  void Set_min(double val);
  void Set_max(double val);
  double Get_min_edge();
  double Get_max_edge();
  double Get_sigma();
  double Get_mean();
  double Get_FWHM();
};

#endif
