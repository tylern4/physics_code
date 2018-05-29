/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "func.hpp"

double func::genNormal(double *x, double *par) {
  double frac = (par[1] / (2 * par[0] * TMath::Gamma(1 / par[1])));
  double expo = TMath::Power(TMath::Abs(x[0] - par[2]) / par[0], par[1]);

  double func = par[3] * frac * TMath::Exp(-expo);

  return func;
}

double func::fiducial_phi(double theta, double e_p) {
  /////////// NOTE: Definitly just magic numbers here......... :(
  double theta_min = 9.5 + 17.0 / (e_p + 0.17);
  double k = 0.705 + 1.1 * e_p;
  double m = -63.5 + (-30.0 * e_p);

  return 37.14 * TMath::Power(TMath::Sin((theta - theta_min) * 0.01745), (k + (m / theta) + (1500. / (theta * theta))));
}

double func::fiducial_phi(double *x, double *par) {
  /////////// NOTE: Definitly just magic numbers here......... :(
  double theta_min = par[3] + par[4] / (x[1] + par[5]);
  double k = par[6] + par[7] * x[1];
  double m = par[8] + (par[9] * x[1]);

  return par[0] * TMath::Power(TMath::Sin((x[0] - theta_min) * par[1]), (k + (m / x[0]) + (par[2] / (x[0] * x[0]))));
}

double func::breit_wigner(double *x, double *par) { return par[2] * TMath::BreitWigner(x[0], par[0], par[1]); }

double func::gausian(double *x, double *par) {
  double g1 = par[2] * TMath::Gaus(x[0], par[0], par[1], true);
  return g1;
}

double func::gausian2(double *x, double *par) {
  double g1 = par[2] * TMath::Gaus(x[0], par[0], par[1], true);
  double g2 = par[5] * TMath::Gaus(x[0], par[3], par[4], true);
  return g1 + g2;
}

double func::peak(double *x, double *par) {
  // double g1 = par[2] * TMath::Gaus(x[0], par[0], par[1], false);
  // double g2 = par[5] * TMath::Gaus(x[0], par[3], par[4], true);
  return par[0] * TMath::Gaus(x[0], par[1], par[2], true);
}

double func::landau(double *x, double *par) { return par[2] * TMath::Landau(x[0], par[0], par[1], true); }

double func::horizontal(double *x, double *par) { return par[0]; }
double func::line(double *x, double *par) { return x[0] * par[1] + par[0]; }
double func::quad(double *x, double *par) { return x[0] * x[0] * par[1] + x[0] * par[1] + par[0]; }

double func::pol0(double *x, double *par) { return par[0]; }
double func::pol1(double *x, double *par) {
  double result = pol0(x, par);
  result += (x[0] * par[1]);
  return result;
}
double func::pol2(double *x, double *par) {
  double result = pol1(x, par);
  result += (x[0] * x[0] * par[2]);
  return result;
}
double func::pol3(double *x, double *par) {
  double result = pol2(x, par);
  result += (x[0] * x[0] * x[0] * par[3]);
  return result;
}
double func::pol4(double *x, double *par) {
  double result = pol3(x, par);
  result += (x[0] * x[0] * x[0] * x[0] * par[4]);
  return result;
}
double func::pol5(double *x, double *par) {
  double result = pol4(x, par);
  result += (x[0] * x[0] * x[0] * x[0] * x[0] * par[5]);
  return result;
}

double func::theta_cc_fit(double *x, double *par) {
  return par[0] + par[1] * x[0] + par[2] * -1 * TMath::Exp(x[0]);
  // + TMath::Exp(x[0]) * (par[0] + par[1] * x[0] + par[2] * x[0] * x[0]);
}
double func::dt_fit(double *x, double *par) { return par[0] + par[1] * x[0]; }  //+ par[2] * TMath::Exp(x[0]); }
