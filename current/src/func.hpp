/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef FUNC_HPP
#define FUNC_HPP
#include "TMath.h"

namespace func {
double ec_fit_func(double *x, double *par);
double genNormal(double *x, double *par);
double breit_wigner(double *x, double *par);
double fiducial_phi(double theta, double e_p);
double fiducial_phi(double *x, double *par);

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
double dt_fit(double *x, double *par);
double theta_cc_fit(double *x, double *par);
}  // namespace func

#endif
