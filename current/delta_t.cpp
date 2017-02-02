/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "delta_t.hpp"

Delta_T::Delta_T() {}
Delta_T::~Delta_T() {}

double Delta_T::vertex_time(double sc_time, double sc_pathlength,
                            double cut_beta) {
  return sc_time - sc_pathlength / (cut_beta * c_special_units);
}

double Delta_T::delta_t(double electron_vertex_time, double mass,
                        double momentum, double sc_t, double sc_r) {
  double cut_beta = 1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
  return electron_vertex_time - vertex_time(sc_t, sc_r, cut_beta);
}

void Delta_T::delta_t_cut(Histogram *hists, bool first_run) {
  double delta_t_P, delta_t_PIP, delta_t_ELECTRON;
  double electron_vertex = vertex_time(sc_t[sc[0] - 1], sc_r[sc[0] - 1], 1.0);
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

    if (first_run) {
      delta_t_P = delta_t(electron_vertex, MASS_P, mom, sct, scr);
      delta_t_PIP = delta_t(electron_vertex, MASS_PIP, mom, sct, scr);
      delta_t_ELECTRON = delta_t(electron_vertex, MASS_E, mom, sct, scr);
    } else {
      delta_t_P = dt_proton->at(event_number);
      delta_t_PIP = dt_pip->at(event_number);
      delta_t_ELECTRON = delta_t(electron_vertex, MASS_E, mom, sct, scr);
    }

    if (charge == 1) {
      hists->Fill_deltat_P(mom, delta_t_P);
      hists->Fill_deltat_PIP(mom, delta_t_PIP);
      hists->Fill_deltat_positron(mom, delta_t_ELECTRON);
      if (is_proton->at(event_number) && ID == PROTON)
        hists->Fill_deltat_P_PID(mom, delta_t_P);
      if (is_pip->at(event_number) && ID == PIP)
        hists->Fill_deltat_PIP_PID(mom, delta_t_PIP);
      if (ID == -11)
        hists->Fill_deltat_positron_PID(mom, delta_t_ELECTRON);
    } else if (charge == -1) {
      hists->Fill_deltat_electron(mom, delta_t_ELECTRON);
      if (is_electron->at(event_number) && ID == ELECTRON)
        hists->Fill_deltat_electron_PID(mom, delta_t_ELECTRON);
    }

    hists->delta_t_Fill(mom, charge, delta_t_P, delta_t_PIP, delta_t_ELECTRON);
    hists->delta_t_sec_pad(mom, charge, delta_t_P, delta_t_PIP,
                           delta_t_ELECTRON, sc_sector, sc_paddle);
  }
}
