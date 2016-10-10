/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef FIT_FUNCTIONS_H_GUARD
#define FIT_FUNCTIONS_H_GUARD
#include "main.h"

double Proton_Pos_fit(double x){
	return 8.18213 * exp(-4.47734 * x) + 0.452738 * x +0.233145;
}

double Proton_Neg_fit(double x){
	return -95.1733 * exp(-11.1013 * x) + -0.0244971 * x + -0.840824;
}

double Pip_Pos_fit(double x){
	return 3.53392 * exp(-3.86201 * x) + 0.651278 * x + 0.267719;
}

double Pip_Neg_fit(double x){
	return -3.92061 * exp(-2.58462 * x) + -1.74244 * x + 1.08991;
}



#endif