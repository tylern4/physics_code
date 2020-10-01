/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef DELTA_T_H
#define DELTA_T_H
#include "branches.hpp"
#include "constants.hpp"
#include "histogram.hpp"

class Delta_T {
 private:
  std::shared_ptr<Branches> _data = nullptr;
  const float c_special_units = 29.9792458;
  float _vertex = 0.0;
  std::vector<float> _elec_array;
  std::vector<float> _proton_array;
  std::vector<float> _pion_array;
  std::vector<float> _kaon_array;

  friend std::ostream &operator<<(std::ostream &os, Delta_T const &m);
  float vertex_time(float sc_time, float sc_pathlength, float relatavistic_beta);
  std::vector<float> _delta_t_vec(float mass);
  float _delta_t(float mass, float momentum, float sc_t, float sc_r);

 public:
  Delta_T(const std::shared_ptr<Branches> &data);
  ~Delta_T();

  void deltat(float momentum, float sc_t, float sc_r);
  float Get_dt_E(int part);
  float Get_dt_P(int part);
  float Get_dt_Pi(int part);
  float Get_dt_K(int part);
  float Get_vertex(int part);

  void delta_t_hists(const std::shared_ptr<Histogram> &hists);
};
#endif
