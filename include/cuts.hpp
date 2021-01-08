/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include <iostream>
#include "branches.hpp"
#include "constants.hpp"
#include "delta_t.hpp"
#include "func.hpp"
#include "physics.hpp"
#include "reaction.hpp"

class Cuts {
 protected:
  std::shared_ptr<Branches> _data = nullptr;
  std::shared_ptr<Delta_T> _dt = nullptr;

  int num_phe = 0;
  int gpart = 0;
  int charge = 0;
  float electron_p = 0;
  float samp_frac = 0;
  float _vx = 0;
  float _vy = 0;
  float _vz = 0.3;
  float _theta = 0.0;
  float _phi, _phi_cent = 0.0;
  int _sec = 0;

  bool ec_cut = false;
  bool cc_cut = false;
  bool stat_cut = false;
  bool sc_cut = false;
  bool dc_cut = false;
  bool dc_stat_cut = false;

  bool samp_frac_cut = false;

 public:
  Cuts(const std::shared_ptr<Branches>& data);
  Cuts(const std::shared_ptr<Branches>& data, const std::shared_ptr<Delta_T>& dt);
  ~Cuts();
  bool check_banks();
  void Set_elec_fid();
  bool isElectron();
  bool isStrictElectron();
  bool Fid_cut();
  bool Beam_cut();

  const std::shared_ptr<Delta_T>& share_dt() { return _dt; }

  bool Electron_fid_arjun();
  bool Hardon_fid_arjun(int part);

  float hardon_fid_phi(int part);
  float hadron_fid_phi_min(float theta, int sec);
  float hadron_fid_phi_max(float theta, int sec);

  bool Pip(int part);
  bool Pim(int part);
  bool Prot(int part);

  bool Pipish(int part);
  bool Pimish(int part);
  bool Protish(int part);

  bool dt_P_cut(int i);
  bool dt_K_cut(int i);
  bool dt_Pip_cut(int i);
  bool elec_fid_cut();
  bool fid_chern_cut();
};

class e1d_Cuts : public Cuts {
 public:
  e1d_Cuts(const std::shared_ptr<Branches>& data) : Cuts(data){};
  e1d_Cuts(const std::shared_ptr<Branches>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){};
  bool isElectron();
  bool Beam_cut();
  float sf_top_fit(float P);
  float sf_bot_fit(float P);
  bool sf_cut(float sf, float P);
};

class e1f_Cuts : public Cuts {
 public:
  e1f_Cuts(const std::shared_ptr<Branches>& data) : Cuts(data){};
  e1f_Cuts(const std::shared_ptr<Branches>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){};
  bool isElectron();
  bool Beam_cut();
  float sf_top_fit(float P);
  float sf_bot_fit(float P);
  bool sf_cut(float sf, float P);
};

class e16_Cuts : public Cuts {
 public:
  e16_Cuts(const std::shared_ptr<Branches>& data) : Cuts(data){};
  e16_Cuts(const std::shared_ptr<Branches>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){};
  bool isElectron();
  bool Beam_cut();
};

#endif