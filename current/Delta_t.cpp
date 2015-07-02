#include <cmath>

// speed of light in appropriate units for delta t
const double c_special_units = 29.9792458;

inline double vertex_time(double sc_time, double sc_pathlength, double cut_beta) {
	return sc_time - sc_pathlength/(cut_beta * c_special_units);
}

inline double electron_vertex_time(int* sc, double* sc_t, double* sc_r) {
	return vertex_time(sc_t[sc[0]-1], sc_r[sc[0]-1], 1.0);
}

inline double delta_t(int index, double mass, double* p, int* sc, double* sc_t, double* sc_r){
	const double momentum = p[index];
	const double cut_beta = 1.0/sqrt(1.0 + pow(mass/momentum,2.0));
	return electron_vertex_time(sc,sc_t,sc_r) - vertex_time(sc_t[sc[index]-1],	sc_r[sc[index]-1],	cut_beta);
}
