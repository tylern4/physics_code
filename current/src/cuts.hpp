/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef CUTS_HPP
#define CUTS_HPP
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

class Cuts {
private:
  double par_max, par_mean, par_FWHM;
  std::string ploy1d = "[0]+[1]*x";
  std::string ploy2d = ploy1d + "+[2]*x*x";
  std::string ploy3d = ploy2d + "+[3]*x*x*x";
  std::string ploy4d = ploy3d + "+[4]*x*x*x*x";
  std::string gaus = "[0]*TMath::Gaus(x,[1],[2],1)";

public:
  double mean, sigma, FWHM;
  double a, b, c, d, e;
  void FitGaus(TH1D *hist, double min_value, double max_value);
  void FitPoly_1D(TH1D *hist, double min_value, double max_value);
  void FitPoly_2D(TH1D *hist, double min_value, double max_value);
  void FitPoly_3D(TH1D *hist, double min_value, double max_value);
  void FitPoly_4D(TH1D *hist, double min_value, double max_value);
};

#endif