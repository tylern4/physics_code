\\Auto Generated fit code from e1d
#ifndef FIT_FUNCTIONS_H_GUARD
#define FIT_FUNCTIONS_H_GUARD
#include "main.h"

double Proton_Pos_fit(double x){
	return 21.2822*exp(-5.98669*x)+0.538534*x+0.591178;
}
double Proton_Neg_fit(double x){
	return -112.119*exp(-10.528*x)+-0.122261*x+-1.16527;
}
double Pip_Pos_fit(double x){
	return 5.33217*exp(-3.71301*x)+1.19119*x+0.162472;
}
double Pip_Neg_fit(double x){
	return -5.65144*exp(-2.87401*x)+-2.24482*x+1.12397;
}
#endif

