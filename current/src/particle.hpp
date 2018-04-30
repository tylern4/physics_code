/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef PARTICLE_T_H
#define PARTICLE_T_H
#include "physics.hpp"

class Particle {
 private:
  TLorentzVector *vec;
  double p, cx, cy, cz, px, py, pz, mass, theta, phi;
  int sector, pid;

 public:
  Particle();
  Particle(double p, double cx, double cy, double cz, int id, double mass);
  Particle(double p, double cx, double cy, double cz, int id);
  ~Particle();

  inline void ReCalc();
  double GetPx();
  double GetPy();
  double GetPz();
  double GetP();
  double GetE();
  double GetMass();
  double GetTheta();
  double GetPhi();
  int GetSector();
  int GetPid();
  TLorentzVector Get4Vec();

  double SetPx(double _px);
  double SetPy(double _py);
  double SetPz(double _pz);
  double SetP(double _p);
  double SetMass(double _mass);
  double SetTheta(double _theta);
  double SetPhi(double _phi);
  int SetSector(int _sec);
  int SetPid(int _pid);
};

class Electron : public Particle {
 private:
  double W, Q2;
  TLorentzVector *beam = new TLorentzVector(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)), E1D_E0);
  TLorentzVector *q_mu = new TLorentzVector();

 public:
  Electron(double p, double cx, double cy, double cz);
  ~Electron();
  double GetW();
  double GetQ2();
  TLorentzVector Get_q_mu();
};

#endif
