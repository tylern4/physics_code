/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "event.hpp"

Event::Event() {}
Event::Event(double p, double cx, double cy, double cz, int pid, double mass) {
  this->p = p;
  this->cx = cx;
  this->cy = cy;
  this->cz = cz;
  this->px = p * cx;
  this->py = p * cy;
  this->pz = p * cz;
  this->pid = pid;

  this->mass = mass;

  this->theta = physics::theta_calc(this->cz);
  this->phi = physics::phi_calc(this->cx, this->cy);
  this->sector = physics::get_sector(this->phi);

  this->vec = physics::fourVec(this->px, this->py, this->pz, physics::Get_Mass(this->pid));
}
Event::Event(double p, double cx, double cy, double cz, int pid) {
  this->p = p;
  this->cx = cx;
  this->cy = cy;
  this->cz = cz;
  this->px = p * cx;
  this->py = p * cy;
  this->pz = p * cz;
  this->pid = pid;

  this->mass = physics::Get_Mass(this->pid);

  this->theta = physics::theta_calc(this->cz);
  this->phi = physics::phi_calc(this->cx, this->cy);
  this->sector = physics::get_sector(this->phi);

  this->vec = physics::fourVec(this->px, this->py, this->pz, physics::Get_Mass(this->pid));
}
Event::~Event() {}
