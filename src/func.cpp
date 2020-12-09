/**************************************/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "func.hpp"
#include <iostream>

ROOT::Double_v func::genNormal(const ROOT::Double_v *x, const double *par) {
  double frac = (par[1] / (2 * par[0] * TMath::Gamma(1 / par[1])));
  double expo = TMath::Power(TMath::Abs(x[0] - par[2]) / par[0], par[1]);

  double func = par[3] * frac * TMath::Exp(-expo);

  return func;
}

ROOT::Double_v func::ec_fit_func(const ROOT::Double_v *x, const double *par) {
  double func = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] * x[0] * x[0] * x[0] * x[0];
  return func;
}

/*
ROOT::Double_v func::fiducial_phi(double theta, double e_p) {
  /////////// NOTE: Definitly just magic numbers here......... :(
  double theta_min = 9.5 + 17.0 / (e_p + 0.17);
  double k = 0.705 + 1.1 * e_p;
  double m = -63.5 + (-30.0 * e_p);

  return 37.14 * TMath::Power(TMath::Sin((theta - theta_min) * 0.01745), (k + (m / theta) + (1500. / (theta * theta))));
}

ROOT::Double_v func::fiducial_phi(const ROOT::Double_v *x, const double *par) {
  double theta_min = par[3] + par[4] / (x[1] + par[5]);
  double k = par[6] + par[7] * x[1];
  double m = par[8] + (par[9] * x[1]);

  return par[0] * TMath::Power(TMath::Sin((x[0] - theta_min) * par[1]), (k + (m / x[0]) + (par[2] / (x[0] * x[0]))));
}
*/

ROOT::Double_v func::fiducial(const ROOT::Double_v *x, const double *par) {
  double func = 0.0;
  func += par[0] * x[0] * x[0] * x[0] * x[0];
  func += par[1] * x[0] * x[0] * x[0];
  func += par[2] * x[0] * x[0];
  func += par[3] * x[0];
  func += par[4];
  func += TMath::Exp(x[0] * par[5]);

  return func;
}

ROOT::Double_v func::breit_wigner(const ROOT::Double_v *x, const double *par) {
  return par[2] * TMath::BreitWigner(x[0], par[0], par[1]);
}

ROOT::Double_v func::gausian(const ROOT::Double_v *x, const double *par) {
  double g1 = par[0] * TMath::Gaus(x[0], par[1], par[2], false);
  return g1;
}

ROOT::Double_v func::gausian2(const ROOT::Double_v *x, const double *par) {
  double g1 = par[2] * TMath::Gaus(x[0], par[0], par[1], true);
  double g2 = par[5] * TMath::Gaus(x[0], par[3], par[4], true);
  return g1 * g2;
}

ROOT::Double_v func::peak(const ROOT::Double_v *x, const double *par) {
  double g1 = par[0] * TMath::Gaus(x[0], par[1], par[2], true);
  return g1;
}

ROOT::Double_v func::landau(const ROOT::Double_v *x, const double *par) {
  return par[2] * TMath::Landau(x[0], par[0], par[1], true);
}

ROOT::Double_v func::horizontal(const ROOT::Double_v *x, const double *par) { return par[0]; }
ROOT::Double_v func::line(const ROOT::Double_v *x, const double *par) { return x[0] * par[1] + par[0]; }
ROOT::Double_v func::quad(const ROOT::Double_v *x, const double *par) {
  return x[0] * x[0] * par[2] + x[0] * par[1] + par[0];
}

ROOT::Double_v func::pol0(const ROOT::Double_v *x, const double *par) { return par[0]; }
ROOT::Double_v func::pol1(const ROOT::Double_v *x, const double *par) {
  double result = pol0(x, par);
  result += (x[0] * par[1]);
  return result;
}
ROOT::Double_v func::pol2(const ROOT::Double_v *x, const double *par) {
  double result = pol1(x, par);
  result += (x[0] * x[0] * par[2]);
  return result;
}
ROOT::Double_v func::pol3(const ROOT::Double_v *x, const double *par) {
  double result = pol2(x, par);
  result += (x[0] * x[0] * x[0] * par[3]);
  return result;
}
ROOT::Double_v func::pol4(const ROOT::Double_v *x, const double *par) {
  double result = pol3(x, par);
  result += (x[0] * x[0] * x[0] * x[0] * par[4]);
  return result;
}

ROOT::Double_v func::pol5(const ROOT::Double_v *x, const double *par) {
  double result = pol4(x, par);
  result += (x[0] * x[0] * x[0] * x[0] * x[0] * par[5]);
  return result;
}

ROOT::Double_v func::fid(const ROOT::Double_v *x, const double *par) {
  double result = 0;
  for (int i = 0; i <= 2; i++) {
    result += par[i] * TMath::Power(x[0], i);
  }
  return result;
}

ROOT::Double_v func::theta_cc_fit(const ROOT::Double_v *x, const double *par) {
  // intercept, slope, const, exp_const, pol_const
  return par[0] + par[1] * x[0] + par[2] * TMath::Exp(par[3] * x[0]);
  // + TMath::Exp(x[0]) * (par[0] + par[1] * x[0] + par[2] * x[0] * x[0]);
}
ROOT::Double_v func::dt_fit(const ROOT::Double_v *x, const double *par) {
  return par[0] + par[1] * x[0];
}  //+ par[2] * TMath::Exp(x[0]); }

// Lorentzian Peak function
ROOT::Double_v func::missMasspeak(const ROOT::Double_v *x, const double *par) {
  double arg1 = 2 / TMath::Pi();                    // 2 over pi
  double arg2 = par[1] * par[1] * par[2] * par[2];  // Gamma=par[1]  M=par[2]
  double arg3 = ((x[0] * x[0]) - (par[2] * par[2])) * ((x[0] * x[0]) - (par[2] * par[2]));
  double arg4 = x[0] * x[0] * x[0] * x[0] * ((par[1] * par[1]) / (par[2] * par[2]));
  return par[0] * arg1 * arg2 / (arg3 + arg4);
}

ROOT::Double_v func::missMassbackground(const ROOT::Double_v *x, const double *par) {
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] +
         par[4] * x[0] * x[0] * x[0] * x[0];
}

// Sum of background and peak function
ROOT::Double_v func::missMassfitFunction(const ROOT::Double_v *x, const double *par) {
  return missMasspeak(x, par) + missMasspeak(x, &par[3]) + missMasspeak(x, &par[6]) + missMassbackground(x, &par[9]);
}

ROOT::Double_v func::degauss(const ROOT::Double_v *x, const double *par) {
  // The two functions that were convoluted together are:
  //
  // 1. Gaussian with A, mu, and sigma for parameters.
  // Asymmetrical double-exponential function that looks like
  //* exp(lambda1*x) for x<0
  //* exp(-lambda2*x) for x>=0.

  // input
  double xx = x[0];
  // Gaussian smearing parameters
  double A = par[0];
  double mu = par[1];
  double sigma = par[2];
  // Asymmetrical double-exponential parameters
  double lambda1 = par[3];
  double lambda2 = par[4];
  // Useful quantities
  double mu1 = sigma * sigma * lambda1 + xx - mu;
  double mu2 = -sigma * sigma * lambda2 + xx - mu;
  double ret =
      A * 0.5 / (1.0 / lambda1 + 1.0 / lambda2) *
      (exp(0.5 * TMath::Power(sigma * lambda1, 2) + lambda1 * (xx - mu)) * TMath::Erfc(mu1 / (sigma * sqrt(2.0))) +
       exp(0.5 * TMath::Power(sigma * lambda2, 2) - lambda2 * (xx - mu)) * TMath::Erfc(-mu2 / (sigma * sqrt(2.0))));

  return ret;
}

float func::dt_poly4(std::vector<double> params, float p) {
  float y = 0;
  y += params[0] * p * p * p * p;
  y += params[1] * p * p * p;
  y += params[2] * p * p;
  y += params[3] * p;
  y += params[4];

  return y;
}

bool func::fid_chern(float x, float y) {
  float P0 = 55;
  float P1 = -1.75;
  return (x > P0 + P1 * y && x > P0 - P1 * y);
}

float func::log_sqrt_pol1(std::vector<double> params, float p) {
  float y = params[0] * logf(params[1] * p) + params[2] * p * sqrtf(params[3] * p) + params[4] * p + params[5];
  return y;
}