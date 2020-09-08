/**************************************/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef FUNC_HPP
#define FUNC_HPP
#include <cmath>
#include "TMath.h"
#include "TSpline.h"

namespace func {
double ec_fit_func(double *x, double *par);
double ec_fit_func_arjun(double *x, double *par);
double ec_fit_func_invert(double *x, double *par);
double genNormal(double *x, double *par);
double breit_wigner(double *x, double *par);
double fiducial_phi(double theta, double e_p);
double fiducial_phi(double *x, double *par);
double fiducial(double *x, double *par);

double gausian(double *x, double *par);
double gausian2(double *x, double *par);
double landau(double *x, double *par);
double peak(double *x, double *par);

double horizontal(double *x, double *par);
double line(double *x, double *par);
double quad(double *x, double *par);

double pol0(double *x, double *par);
double pol1(double *x, double *par);
double pol2(double *x, double *par);
double pol3(double *x, double *par);
double pol4(double *x, double *par);
double pol5(double *x, double *par);
double fid(double *x, double *par);
double dt_fit(double *x, double *par);
double theta_cc_fit(double *x, double *par);

double missMasspeak(double *x, double *par);
double missMassbackground(double *x, double *par);
double missMassfitFunction(double *x, double *par);
double degauss(double *x, double *par);

float dt_poly4(std::vector<double> params, float p);

bool fid_chern(float x, float y);
float log_sqrt_pol1(std::vector<double> params, float p);

}  // namespace func

#endif
