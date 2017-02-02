//Auto Generated fit code from e1d
#ifndef FIT_FUNCTIONS_H_GUARD
#define FIT_FUNCTIONS_H_GUARD
#include "main.h"

double Proton_Pos_fit(double x){
	return 107.513*exp(-13.0565*x)+-0.532738*x+1.8416;
}
double Proton_Neg_fit(double x){
	return -1181.36*exp(-19.5255*x)+0.608919*x+-1.99309;
}
double Pip_Pos_fit(double x){
	return 2.68699*exp(-3.69242*x)+0.434782*x+0.891199;
}
double Pip_Neg_fit(double x){
	return -3.08468*exp(-4.46416*x)+-0.762092*x+-0.583121;
}
#endif

