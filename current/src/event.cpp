/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "event.hpp"

Event::Event() {}
Event::Event(double p, double cx, double cy, double cz, int pid, double mass) {
  p_e = p;
  cx_e = cx;
  cy_e = cy;
  cz_e = cz;
  px_e = p * cx;
  py_e = p * cy;
  pz_e = p * cz;
  pid_e = pid;

  mass_e = mass;

  theta_e = physics::theta_calc(cz_e);
  phi_e = physics::phi_calc(cx_e, cy_e);
  sector_e = physics::get_sector(phi_e);

  vec = physics::fourVec(px_e, py_e, pz_e, physics::Get_Mass(pid_e));
}
Event::Event(double p, double cx, double cy, double cz, int pid) {
  p_e = p;
  cx_e = cx;
  cy_e = cy;
  cz_e = cz;
  px_e = p * cx;
  py_e = p * cy;
  pz_e = p * cz;
  pid_e = pid;

  mass_e = physics::Get_Mass(pid_e);

  theta_e = physics::theta_calc(cz_e);
  phi_e = physics::phi_calc(cx_e, cy_e);
  sector_e = physics::get_sector(phi_e);

  vec = physics::fourVec(px_e, py_e, pz_e, physics::Get_Mass(pid_e));
}
Event::~Event() {}
