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
  double min_edge_x = std::nan("0"), max_edge_x = std::nan("0");
  std::string ploy1d = "[0]+[1]*x";
  std::string ploy2d = ploy1d + "+[2]*x*x";
  std::string ploy3d = ploy2d + "+[3]*x*x*x";
  std::string ploy4d = ploy3d + "+[4]*x*x*x*x";
  std::string gaus = "[0]*TMath::Gaus(x,[1],[2],1)";
  double fiducial_phi_lo(double theta_e, double theta_e_min, double k, double m,
                         int c);
  double fiducial_phi_hi(double theta_e, double theta_e_min, double k, double m,
                         int c);

public:
  double mean, sigma, FWHM;
  double a, b, c, d, e = 0.0;
  void FitGaus(TH1D *hist, double min_value, double max_value);
  void FitLandau(TH1D *hist, double min_value, double max_value);
  void FitPoly_1D(TH1D *hist, double min_value, double max_value);
  void FitPoly_2D(TH1D *hist, double min_value, double max_value);
  void FitPoly_3D(TH1D *hist, double min_value, double max_value);
  void FitPoly_4D(TH1D *hist, double min_value, double max_value);
  void FitPoly_fid(TH2D *hist, double min_value, double max_value);
  void FitFiducial_lo(TH2D *hist2d, double min_value, double max_value);
  void FitFiducial_hi(TH2D *hist2d, double min_value, double max_value);
  void FitFiducial(TH2D *hist2d, double min_value, double max_value);
  void FitGenNormal(TH1D *hist, double min_value, double max_value);
  double Get_min_edge();
  double Get_max_edge();
};

#endif
