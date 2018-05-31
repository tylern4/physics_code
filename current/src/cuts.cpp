/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "cuts.hpp"

double Cuts::sf_top_fit(double P) { return 0.363901 + -0.00992778 * P + 5.84749e-06 * P * P * P * P * P * P; }
double Cuts::sf_bot_fit(double P) { return 0.103964 + 0.0524214 * P + -3.64355e-05 * P * P * P * P * P * P; }
bool Cuts::sf_cut(double sf, double P) { return ((sf > sf_bot_fit(P)) && (sf < sf_top_fit(P))); }

double Cuts::dt_P_bot_fit(double P) { return -1.509 + 0.4172 * P; }
double Cuts::dt_P_top_fit(double P) { return 1.307 - 0.3473 * P; }
bool Cuts::dt_P_cut(double dt, double P) { return ((dt > dt_P_bot_fit(P)) && (dt < dt_P_top_fit(P))); }

double Cuts::dt_Pip_bot_fit(double P) { return -0.9285 - 0.04094 * P; }
double Cuts::dt_Pip_top_fit(double P) { return 0.9845 - 0.05473 * P; }
bool Cuts::dt_Pip_cut(double dt, double P) { return ((dt > dt_Pip_bot_fit(P)) && (dt < dt_Pip_top_fit(P))); }
