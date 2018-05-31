/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "constants.hpp"

namespace Cuts {
double sf_top_fit(double P);
double sf_bot_fit(double P);
bool sf_cut(double sf, double P);

double dt_P_bot_fit(double P);
double dt_P_top_fit(double P);
bool dt_P_cut(double dt, double P);

double dt_Pip_bot_fit(double P);
double dt_Pip_top_fit(double P);
bool dt_Pip_cut(double dt, double P);
}  // namespace Cuts
#endif
