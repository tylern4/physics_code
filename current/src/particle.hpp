/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef PARTICLE_T_H
#define PARTICLE_T_H
#include <iostream>
#include "physics.hpp"

class Particle {
 private:
  TLorentzVector _vec4;
  TVector3 _vertex;
  TVector3 _cos_angle;
  float _p = float(NULL);
  float _mass = float(NULL);
  float _beta = float(NULL);
  int _bank_PID = int(NULL);
  int _PID = int(NULL);
  int _charge = int(NULL);

 public:
  Particle();
  Particle(double p, double cx, double cy, double cz, int pid);
  ~Particle();
  void Set_momentum(double p, double cx, double cy, double cz);
  void Set_BankID(int pid);
  void Set_PID(int pid);
  void Set_charge(int charge);
  void Set_beta(float beta);
  void Set_mass(float mass);
  void Set_vertex(double vx, double vy, double vz);
  void setup_4vector();

  TLorentzVector Particle4Vec();
  TVector3 vertex();
  int BankID();
  int PID();
  int charge();
  float beta();
  float mass();
};

#endif
