/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef DELTA_T_H
#define DELTA_T_H
#include "histogram.hpp"
#include "constants.hpp"
#include "branches.hpp"

class Delta_T {
 private:
  const double c_special_units = 29.9792458;
  std::vector<double> masses = {MASS_E, MASS_P, MASS_PIP};
  double vertex = 0.0;
  double dt_E = 0.0;
  double dt_P = 0.0;
  double dt_Pi = 0.0;
  bool first_run = true;

  double vertex_time(double sc_time, double sc_pathlength,
                     double relatavistic_beta);

 public:
  Delta_T(double sc_time, double sc_pathlength);
  ~Delta_T();

  void deltat(double momentum, double sc_t, double sc_r);
  double Get_dt_E();
  double Get_dt_P();
  double Get_dt_Pi();

  double delta_t(double electron_vertex_time, double mass, double momentum,
                 double sc_t, double sc_r);
  void delta_t_hists(Histogram *hists);
  double *delta_t_array(double *dt_array, double mass);
  std::vector<double> delta_t_array(double mass, int num_parts);
};
#endif
