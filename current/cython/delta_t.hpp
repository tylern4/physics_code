/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef DELTA_T_H
#define DELTA_T_H
#include <mutex>
#include <thread>
#include "branches.hpp"
#include "constants.hpp"

class Delta_T {
 private:
  const double c_special_units = 29.9792458;
  std::vector<double> masses = {MASS_E, MASS_P, MASS_PIP, MASS_KP};
  double vertex = 0.0;
  double dt_E = 0.0;
  double dt_P = 0.0;
  double dt_Pi = 0.0;
  double dt_K = 0.0;

  double vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta);

 public:
  Delta_T();
  Delta_T(double sc_time, double sc_pathlength);
  ~Delta_T();

  void deltat(double momentum, double sc_t, double sc_r);
  double deltat(double momentum, double sc_t, double sc_r, double mass);
  double Get_dt_E();
  double Get_dt_P();
  double Get_dt_Pi();
  double Get_dt_K();
  double Get_vertex();

  double delta_t(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r);
  std::vector<double> delta_t_array(double mass, int num_parts);
};
#endif
