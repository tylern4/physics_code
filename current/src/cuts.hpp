/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "branches.hpp"
#include "constants.hpp"
#include "func.hpp"
#include "physics.hpp"

class Cuts {
 private:
  std::shared_ptr<Branches> _data;

  int num_phe = 0;
  int gpart = 0;
  int charge = 0;
  double electron_p = 0;
  double samp_frac = 0;
  double _vx = 0;
  double _vy = 0;
  double _vz = 0.3;
  double _theta = 0.0;
  double _phi, _phi_cent = 0.0;
  int _sec = 0;

  bool ec_cut = false;
  bool cc_cut = false;
  bool stat_cut = false;
  bool sc_cut = false;
  bool dc_cut = false;
  bool dc_stat_cut = false;

  bool samp_frac_cut = false;

 public:
  Cuts(std::shared_ptr<Branches> data);
  ~Cuts();
  void Set_elec_fid();
  bool isElecctron();
  bool isStrictElecctron();
  bool Fid_cut();
  bool Beam_cut();

  double sf_top_fit(double P);
  double sf_bot_fit(double P);
  bool sf_cut(double sf, double P);

  double dt_P_bot_fit(double P);
  double dt_P_top_fit(double P);
  bool dt_P_cut(double dt, double P);

  double dt_Pip_bot_fit(double P);
  double dt_Pip_top_fit(double P);
  bool dt_Pip_cut(double dt, double P);

  bool elec_fid_cut();
};

#endif
