/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef EVENT_T_H
#define EVENT_T_H

#include <iostream>
#include "Particle.hpp"
#include "physics.hpp"

class Event {
 private:
  Particle _beam;
  Particle _gamma;
  Particle _electron;
  std::vector<Particle> _events;
  std::vector<int> _PID;
  std::vector<int> _event_sig;

 public:
  Event(Particle Electron);
  ~Event();

  void Add_Part(Particle p);

  std::vector<int> Signiture();
  void PrintSigniture();
};

#endif
