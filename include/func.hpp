/**************************************/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef FUNC_HPP
#define FUNC_HPP
#include <cmath>
#include "TMath.h"
#include "TSpline.h"
#include "VectorizedTMath.h"
#include "constants.hpp"

namespace func {
ROOT::Double_v ec_fit_func(const ROOT::Double_v *x, const double *par);
ROOT::Double_v ec_fit_func_arjun(const ROOT::Double_v *x, const double *par);
ROOT::Double_v ec_fit_func_invert(const ROOT::Double_v *x, const double *par);
ROOT::Double_v genNormal(const ROOT::Double_v *x, const double *par);
ROOT::Double_v breit_wigner(const ROOT::Double_v *x, const double *par);
ROOT::Double_v fiducial_phi(double theta, double e_p);
ROOT::Double_v fiducial_phi(const ROOT::Double_v *x, const double *par);
ROOT::Double_v fiducial(const ROOT::Double_v *x, const double *par);

ROOT::Double_v gausian(const ROOT::Double_v *x, const double *par);
ROOT::Double_v gausian2(const ROOT::Double_v *x, const double *par);
ROOT::Double_v landau(const ROOT::Double_v *x, const double *par);
ROOT::Double_v peak(const ROOT::Double_v *x, const double *par);

ROOT::Double_v horizontal(const ROOT::Double_v *x, const double *par);
ROOT::Double_v line(const ROOT::Double_v *x, const double *par);
ROOT::Double_v quad(const ROOT::Double_v *x, const double *par);

ROOT::Double_v pol0(const ROOT::Double_v *x, const double *par);
ROOT::Double_v pol1(const ROOT::Double_v *x, const double *par);
ROOT::Double_v pol2(const ROOT::Double_v *x, const double *par);
ROOT::Double_v pol3(const ROOT::Double_v *x, const double *par);
ROOT::Double_v pol4(const ROOT::Double_v *x, const double *par);
ROOT::Double_v pol5(const ROOT::Double_v *x, const double *par);
ROOT::Double_v fid(const ROOT::Double_v *x, const double *par);
ROOT::Double_v dt_fit(const ROOT::Double_v *x, const double *par);
ROOT::Double_v theta_cc_fit(const ROOT::Double_v *x, const double *par);

ROOT::Double_v missMasspeak(const ROOT::Double_v *x, const double *par);
ROOT::Double_v missMassbackground(const ROOT::Double_v *x, const double *par);
ROOT::Double_v missMassfitFunction(const ROOT::Double_v *x, const double *par);
ROOT::Double_v degauss(const ROOT::Double_v *x, const double *par);

float dt_poly4(std::vector<double> params, float p);

bool fid_chern(float x, float y);
bool fid_chern_mc(float x, float y);
float log_sqrt_pol1(std::vector<double> params, float p);
float log_pol2(std::vector<double> params, float p);

bool pip_sec1_cut1(float p, float theta);
bool pip_sec2_cut1(float p, float theta);
bool pip_sec3_cut1(float p, float theta);
bool pip_sec4_cut1(float p, float theta);
bool pip_sec5_cut1(float p, float theta);
bool pip_sec6_cut1(float p, float theta);

bool thetaMin(float p, float theta);

}  // namespace func

#endif
