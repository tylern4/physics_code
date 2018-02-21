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
  double p, cx, cy, cz, px, py, pz, mass, theta, phi;
  int sector, pid;

 public:
  Event();
  Event(double p, double cx, double cy, double cz, int pid, double mass);
  Event(double p, double cx, double cy, double cz, int pid);
  ~Event();
};

#endif
