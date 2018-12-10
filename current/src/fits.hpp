/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef FITS_HPP
#define FITS_HPP
#include "Math/MinimizerOptions.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "func.hpp"
#include "physics.hpp"

class Fits {
 private:
  int color = 2;
  double par_max, par_mean, par_FWHM, par_RMS;
  double mean, sigma, FWHM;
  double min_value = -1000, max_value = 1000;
  double a, b, c, d, e = 0.0;
  double left_edge_x = std::nan("0"), right_edge_x = std::nan("0");
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
  TF1 *FitGaus(TH1D *hist);
  TF1 *FitLandauGaus(TH1D *hist);
  TF1 *Fit2Gaus(TH1D *hist);
  TF1 *FitLandau(TH1D *hist);
  TF1 *FitPoly_1D(TH1D *hist);
  TF1 *FitPoly_2D(TH1D *hist);
  TF1 *FitPoly_3D(TH1D *hist);
  TF1 *FitPoly_4D(TH1D *hist);
  TF1 *FitPoly_fid(TGraph *hist);
  TF1 *FitFiducial(TGraph *profile, int sec);
  TF1 *FitFiducial_hi(TH2D *hist2d);
  TF1 *FitFiducial(TH2D *hist2d);
  TF1 *FitGenNormal(TH1D *hist);
  TF1 *FitBreitWigner(TH1D *hist);
  TF1 *FitMissMass(TH1D *hist);
  void Set_min(double val);
  void Set_max(double val);
  void Set_lineColor(int val);
  double Get_left_edge();
  double Get_right_edge();
  double Get_sigma();
  double Get_mean();
  double Get_FWHM();
};

#endif
