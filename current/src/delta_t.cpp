/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#include "delta_t.hpp"

Delta_T::Delta_T() {}
Delta_T::Delta_T(double sc_time, double sc_pathlength) { vertex = vertex_time(sc_time, sc_pathlength, 1.0); }

Delta_T::~Delta_T() {}

double Delta_T::vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

void Delta_T::deltat(double momentum, double sc_t, double sc_r) {
  double beta = 0.0;
  double mp = (masses.at(0) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  dt_E = vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(1) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  dt_P = vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(2) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  dt_Pi = vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(3) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  dt_K = vertex - vertex_time(sc_t, sc_r, beta);
}

double Delta_T::Get_dt_E() { return dt_E; }
double Delta_T::Get_dt_P() { return dt_P; }
double Delta_T::Get_dt_Pi() { return dt_Pi; }
double Delta_T::Get_dt_K() { return dt_K; }
double Delta_T::Get_vertex() { return vertex; }

void Delta_T::delta_t_hists(Histogram *hists) {
  Cuts *dt_cut = new Cuts();
  double sct, scr, mom;
  int ID, charge, sc_paddle, sc_sector;
  for (int event_number = 0; event_number < gpart; event_number++) {
    sct = (double)sc_t[sc[event_number] - 1];
    scr = (double)sc_r[sc[event_number] - 1];
    mom = (double)p[event_number];
    ID = (int)id[event_number];
    charge = (int)q[event_number];
    sc_paddle = (int)sc_pd[sc[event_number] - 1];
    sc_sector = (int)sc_sect[sc[event_number] - 1];

    deltat(mom, sct, scr);

    if (charge == 1) {
      hists->Fill_deltat_P(mom, dt_P);
      hists->Fill_deltat_PIP(mom, dt_Pi);
      hists->Fill_deltat_kp(mom, dt_K);
      if (dt_cut->dt_P_cut(dt_P, mom)) hists->Fill_deltat_P_PID(mom, dt_P);
      if (dt_cut->dt_Pip_cut(dt_Pi, mom)) hists->Fill_deltat_PIP_PID(mom, dt_Pi);
      if (dt_cut->dt_P_cut(dt_K, mom)) hists->Fill_deltat_kp_PID(mom, dt_K);
    } else if (charge == -1) {
      hists->Fill_deltat_electron(mom, dt_E);
      if (ID == ELECTRON) hists->Fill_deltat_electron_PID(mom, dt_E);
      hists->Fill_deltat_PIM(mom, dt_Pi);
      if (dt_cut->dt_Pip_cut(dt_Pi, mom)) hists->Fill_deltat_PIM(mom, dt_Pi);
    }

    hists->delta_t_Fill(mom, charge, dt_P, dt_Pi, dt_E);
    hists->delta_t_sec_pad(mom, charge, dt_P, dt_Pi, dt_E, sc_sector, sc_paddle);
  }
}

double Delta_T::delta_t(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r) {
  double cut_beta = 1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
  return electron_vertex_time - vertex_time(sc_t, sc_r, cut_beta);
}

double *Delta_T::delta_t_array(double *dt_array, double mass) {
  Delta_T *dt = new Delta_T(sc_t[sc[0] - 1], sc_r[sc[0] - 1]);
  double electron_vertex = dt->Get_vertex();
  double sct, scr, mom;
  int ID, charge, sc_paddle, sc_sector;

  for (int event_number = 0; event_number < gpart; event_number++) {
    sct = (double)sc_t[sc[event_number] - 1];
    scr = (double)sc_r[sc[event_number] - 1];
    mom = (double)p[event_number];
    ID = (int)id[event_number];
    charge = (int)q[event_number];
    sc_paddle = (int)sc_pd[sc[event_number] - 1];
    sc_sector = (int)sc_sect[sc[event_number] - 1];

    dt_array[event_number] = delta_t(electron_vertex, mass, mom, sct, scr);
  }
  delete dt;
  return dt_array;
}

std::vector<double> Delta_T::delta_t_array(double mass, int num_parts) {
  Delta_T *dt = new Delta_T(sc_t[sc[0] - 1], sc_r[sc[0] - 1]);
  std::vector<double> dt_array(num_parts);
  double electron_vertex = dt->Get_vertex();
  double sct, scr, mom;
  int ID, charge, sc_paddle, sc_sector;

  for (int event_number = 0; event_number < gpart; event_number++) {
    sct = (double)sc_t[sc[event_number] - 1];
    scr = (double)sc_r[sc[event_number] - 1];
    mom = (double)p[event_number];
    ID = (int)id[event_number];
    charge = (int)q[event_number];
    sc_paddle = (int)sc_pd[sc[event_number] - 1];
    sc_sector = (int)sc_sect[sc[event_number] - 1];

    dt_array[event_number] = delta_t(electron_vertex, mass, mom, sct, scr);
  }
  delete dt;
  return dt_array;
}
