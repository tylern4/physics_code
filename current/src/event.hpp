/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef EVENT_T_H
#define EVENT_T_H
#include "physics.hpp"

class Event {
 private:
  TLorentzVector vec;
  double p_e, cx_e, cy_e, cz_e, px_e, py_e, pz_e, mass_e, theta_e, phi_e;
  int sector_e, pid_e;

 public:
  Event();
  Event(double p, double cx, double cy, double cz, int pid, double mass);
  Event(double p, double cx, double cy, double cz, int pid);
  ~Event();
};

#endif
