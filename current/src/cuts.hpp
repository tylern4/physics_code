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

class Cuts {
 private:
  bool electron_cut = false;
  Branches* _data;

  int num_phe = 0;
  int gpart = 0;
  int charge = 0;
  int electron_id = 0;
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
  Cuts();
  Cuts(Branches* data);
  ~Cuts();
  void Set_num_phe(int set);
  void Set_charge(int set);
  void Set_electron_id(int set);
  void Set_gpart(int set);
  void Set_p(double set);
  void Set_Sf(double set);
  void Set_BeamPosition(double x, double y, double z);
  void Set_elec_fid(double theta, double phi, int sec);

  void Set_ec_cut(bool set);
  void Set_cc_cut(bool set);
  void Set_stat_cut(bool set);
  void Set_sc_cut(bool set);
  void Set_dc_cut(bool set);
  void Set_dc_stat_cut(bool set);

  bool isElecctron();
  bool isStrictElecctron();
  bool CheckElectron();
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
