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
double ec_fit_func(const double *x, const double *par);
double ec_fit_func_arjun(const double *x, const double *par);
double ec_fit_func_invert(const double *x, const double *par);
double genNormal(const double *x, const double *par);
double breit_wigner(const double *x, const double *par);
//double fiducial_phi(double theta, double e_p);
double fiducial_phi(const double *x, const double *par);
double fiducial(const double *x, const double *par);

double gausian(const double *x, const double *par);
double gausian2(const double *x, const double *par);
double landau(const double *x, const double *par);
double peak(const double *x, const double *par);

double horizontal(const double *x, const double *par);
double line(const double *x, const double *par);
double quad(const double *x, const double *par);

double pol0(const double *x, const double *par);
double pol1(const double *x, const double *par);
double pol2(const double *x, const double *par);
double pol3(const double *x, const double *par);
double pol4(const double *x, const double *par);
double pol5(const double *x, const double *par);
double fid(const double *x, const double *par);
double dt_fit(const double *x, const double *par);
double theta_cc_fit(const double *x, const double *par);

double missMasspeak(const double *x, const double *par);
double missMassbackground(const double *x, const double *par);
double missMassfitFunction(const double *x, const double *par);
double degauss(const double *x, const double *par);

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
