#include "delta_t.hpp"

inline double delta_t::vertex_time(double sc_time, double sc_pathlength, double cut_beta) {
	return sc_time - sc_pathlength/(cut_beta * c_special_units);
}

inline double delta_t::delta_t_calc(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r){
	double cut_beta = 1.0/sqrt(1.0 + Square(mass/momentum));
	return electron_vertex_time - vertex_time(sc_t,sc_r,cut_beta);
}