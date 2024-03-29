/**************************************/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "func.hpp"
#include <iostream>

double func::genNormal(const double *x, const double *par) {
  double frac = (par[1] / (2 * par[0] * TMath::Gamma(1 / par[1])));
  double expo = TMath::Power(TMath::Abs(x[0] - par[2]) / par[0], par[1]);

  double func = par[3] * frac * TMath::Exp(-expo);

  return func;
}

double func::ec_fit_func(const double *x, const double *par) {
  double func = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] * x[0] * x[0] * x[0] * x[0];
  return func;
}

/*
double func::fiducial_phi(double theta, double e_p) {
  /////////// NOTE: Definitly just magic numbers here......... :(
  double theta_min = 9.5 + 17.0 / (e_p + 0.17);
  double k = 0.705 + 1.1 * e_p;
  double m = -63.5 + (-30.0 * e_p);

  return 37.14 * TMath::Power(TMath::Sin((theta - theta_min) * 0.01745), (k + (m / theta) + (1500. / (theta * theta))));
}

double func::fiducial_phi(const double *x, const double *par) {
  double theta_min = par[3] + par[4] / (x[1] + par[5]);
  double k = par[6] + par[7] * x[1];
  double m = par[8] + (par[9] * x[1]);

  return par[0] * TMath::Power(TMath::Sin((x[0] - theta_min) * par[1]), (k + (m / x[0]) + (par[2] / (x[0] * x[0]))));
}
*/

double func::fiducial(const double *x, const double *par) {
  double func = 0.0;
  func += par[0] * x[0] * x[0] * x[0] * x[0];
  func += par[1] * x[0] * x[0] * x[0];
  func += par[2] * x[0] * x[0];
  func += par[3] * x[0];
  func += par[4];
  func += TMath::Exp(x[0] * par[5]);

  return func;
}

double func::breit_wigner(const double *x, const double *par) {
  return par[2] * TMath::BreitWigner(x[0], par[0], par[1]);
}

double func::gausian(const double *x, const double *par) {
  double g1 = par[0] * TMath::Gaus(x[0], par[1], par[2], false);
  return g1;
}

double func::gausian2(const double *x, const double *par) {
  double g1 = par[2] * TMath::Gaus(x[0], par[0], par[1], true);
  double g2 = par[5] * TMath::Gaus(x[0], par[3], par[4], true);
  return g1 * g2;
}

double func::peak(const double *x, const double *par) {
  double g1 = par[0] * TMath::Gaus(x[0], par[1], par[2], true);
  return g1;
}

double func::landau(const double *x, const double *par) {
  return par[2] * TMath::Landau(x[0], par[0], par[1], true);
}

double func::horizontal(const double *x, const double *par) { return par[0]; }
double func::line(const double *x, const double *par) { return x[0] * par[1] + par[0]; }
double func::quad(const double *x, const double *par) {
  return x[0] * x[0] * par[2] + x[0] * par[1] + par[0];
}

double func::pol0(const double *x, const double *par) { return par[0]; }
double func::pol1(const double *x, const double *par) {
  double result = pol0(x, par);
  result += (x[0] * par[1]);
  return result;
}
double func::pol2(const double *x, const double *par) {
  double result = pol1(x, par);
  result += (x[0] * x[0] * par[2]);
  return result;
}
double func::pol3(const double *x, const double *par) {
  double result = pol2(x, par);
  result += (x[0] * x[0] * x[0] * par[3]);
  return result;
}
double func::pol4(const double *x, const double *par) {
  double result = pol3(x, par);
  result += (x[0] * x[0] * x[0] * x[0] * par[4]);
  return result;
}

double func::pol5(const double *x, const double *par) {
  double result = pol4(x, par);
  result += (x[0] * x[0] * x[0] * x[0] * x[0] * par[5]);
  return result;
}

double func::fid(const double *x, const double *par) {
  double result = 0;
  for (int i = 0; i <= 2; i++) {
    result += par[i] * TMath::Power(x[0], i);
  }
  return result;
}

double func::theta_cc_fit(const double *x, const double *par) {
  // intercept, slope, const, exp_const, pol_const
  return par[0] + par[1] * x[0] + par[2] * TMath::Exp(par[3] * x[0]);
  // + TMath::Exp(x[0]) * (par[0] + par[1] * x[0] + par[2] * x[0] * x[0]);
}
double func::dt_fit(const double *x, const double *par) {
  return par[0] + par[1] * x[0];
}  //+ par[2] * TMath::Exp(x[0]); }

// Lorentzian Peak function
double func::missMasspeak(const double *x, const double *par) {
  double arg1 = 2 / TMath::Pi();                    // 2 over pi
  double arg2 = par[1] * par[1] * par[2] * par[2];  // Gamma=par[1]  M=par[2]
  double arg3 = ((x[0] * x[0]) - (par[2] * par[2])) * ((x[0] * x[0]) - (par[2] * par[2]));
  double arg4 = x[0] * x[0] * x[0] * x[0] * ((par[1] * par[1]) / (par[2] * par[2]));
  return par[0] * arg1 * arg2 / (arg3 + arg4);
}

double func::missMassbackground(const double *x, const double *par) {
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] +
         par[4] * x[0] * x[0] * x[0] * x[0];
}

// Sum of background and peak function
double func::missMassfitFunction(const double *x, const double *par) {
  return missMasspeak(x, par) + missMasspeak(x, &par[3]) + missMasspeak(x, &par[6]) + missMassbackground(x, &par[9]);
}

double func::degauss(const double *x, const double *par) {
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
  // float P0 = 55;
  // float P1 = -1.75;

  float P0 = 48;
  float P1 = 1.75;

  return ((x > P0 - P1 * y) && (x > P0 + P1 * y));
}

bool func::fid_chern_mc(float x, float y) {
  // float P0 = 55;
  // float P1 = -1.75;

  float P0 = 48;
  float P1 = 1.9;

  return ((x > P0 - P1 * y) && (x > P0 + P1 * y));
}

float func::log_sqrt_pol1(std::vector<double> params, float p) {
  float y = params[0] * logf(params[1] * p) + params[2] * p * sqrtf(params[3] * p) + params[4] * p + params[5];
  return y;
}

float func::log_pol2(std::vector<double> params, float p) {
  float y = params[0] * log(params[1] * p) + params[2] * p * p + params[3] * p + params[4];
  return y;
}

bool func::pip_sec1_cut1(float p, float theta) {
  float s = (300 * p * p * p - 300 * p * p + 497.462 * p + 15) * exp(-2.15 * p) + 14.3;
  if (theta < s) return true;
  return false;
}

bool func::pip_sec2_cut1(float p, float theta) {
  float x = p + 0.1;
  float s = (300 * x * x * x - 300 * x * x + 500 * x - 35) * exp(-2.15 * x) + 14.3;
  if (theta < s) return true;
  return false;
}

bool func::pip_sec3_cut1(float p, float theta) {
  float x = p - 0.01;
  float s = (300 * x * x * x - 300 * x * x + 500 * x - 3) * exp(-2.15 * x) + 14.3;
  if (theta < s) return true;
  return false;
}

bool func::pip_sec4_cut1(float p, float theta) {
  float x = p - 0.01;
  float s = (300 * x * x * x - 300 * x * x + 500 * x - 3) * exp(-2.15 * x) + 14.3;
  if (theta < s) return true;
  return false;
}

bool func::pip_sec5_cut1(float p, float theta) {
  float x = p - 0.01;
  float s = (300 * x * x * x - 300 * x * x + 500 * x - 34) * exp(-2.15 * x) + 14.3;
  if (theta < s) return true;
  return false;
}

bool func::pip_sec6_cut1(float p, float theta) {
  float x = p + 0.03;
  float s = (300 * x * x * x - 300 * x * x + 500 * x - 20) * exp(-2.15 * x) + 14.3;
  if (theta < s) return true;
  return false;
}

bool func::thetaMin(float p, float theta) {
  float par[4] = {2.5, 4, 0.001, 4.8};
  float ans = par[0] + (par[1] / (p + par[2])) + par[3] * p;
  if (ans > theta * DEG2RAD)
    return true;
  else
    return false;
}
