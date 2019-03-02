/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#include "delta_t.hpp"

Delta_T::Delta_T() {}
Delta_T::Delta_T(double sc_time, double sc_pathlength) { this->vertex = vertex_time(sc_time, sc_pathlength, 1.0); }
Delta_T::~Delta_T() {}

double Delta_T::vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

void Delta_T::deltat(double momentum, double sc_t, double sc_r) {
  double beta = 0.0;
  double mp = (masses.at(0) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  this->dt_E = this->vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(1) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  this->dt_P = this->vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(2) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  this->dt_Pi = this->vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(3) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  this->dt_K = this->vertex - vertex_time(sc_t, sc_r, beta);
}

double Delta_T::Get_dt_E() { return dt_E; }
double Delta_T::Get_dt_P() { return dt_P; }
double Delta_T::Get_dt_Pi() { return dt_Pi; }
double Delta_T::Get_dt_K() { return dt_K; }
double Delta_T::Get_vertex() { return vertex; }

void Delta_T::delta_t_hists(std::shared_ptr<Histogram> hists, std::shared_ptr<Branches> data) {
  auto dt_cut = std::make_unique<Cuts>(data);
  double sct, scr, mom;
  int ID, charge, sc_paddle, sc_sector;

  for (int event_number = 1; event_number < data->gpart(); event_number++) {
    sct = data->sc_t(event_number);
    scr = data->sc_r(event_number);
    mom = data->p(event_number);
    ID = data->id(event_number);
    charge = data->q(event_number);
    sc_paddle = data->sc_pd(event_number);
    sc_sector = data->sc_sect(event_number);

    this->deltat(mom, sct, scr);

    if (charge == POSITIVE) {
      hists->Fill_deltat_P(mom, dt_P);
      hists->Fill_deltat_PIP(mom, dt_Pi);
      hists->Fill_deltat_kp(mom, dt_K);
      if (dt_cut->dt_P_cut(dt_P, mom)) hists->Fill_deltat_P_PID(mom, dt_P);
      if (dt_cut->dt_Pip_cut(dt_Pi, mom)) hists->Fill_deltat_PIP_PID(mom, dt_Pi);
      if (dt_cut->dt_P_cut(dt_K, mom)) hists->Fill_deltat_kp_PID(mom, dt_K);
    } else {
      hists->Fill_deltat_electron(mom, dt_E);
      if (abs(dt_E) < 0.04)
        hists->Fill_deltat_electron_PID(mom, dt_E);
      else {
        hists->Fill_deltat_PIM(mom, dt_Pi);
        if (dt_cut->dt_Pip_cut(dt_Pi, mom)) hists->Fill_deltat_PIM_PID(mom, dt_Pi);
      }
    }

    hists->delta_t_Fill(mom, charge, dt_P, dt_Pi, dt_E);
    // if (sc_sector <= 6 && sc_sector != 0 && sc_paddle <= 48 && sc_paddle != 0) {
    //  hists->delta_t_sec_pad(mom, charge, dt_P, dt_Pi, dt_E, sc_sector, sc_paddle);
    //}
  }
}

double Delta_T::delta_t(double mass, double momentum, double sc_t, double sc_r) {
  double cut_beta = 1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
  return this->vertex - vertex_time(sc_t, sc_r, cut_beta);
}

double *Delta_T::delta_t_array(double *dt_array, double mass, std::shared_ptr<Branches> data) {
  double sct, scr, mom;

  for (int event_number = 0; event_number < data->gpart(); event_number++) {
    sct = data->sc_t(event_number);
    scr = data->sc_r(event_number);
    mom = data->p(event_number);

    dt_array[event_number] = delta_t(mass, mom, sct, scr);
  }
  return dt_array;
}

std::vector<double> Delta_T::delta_t_array(double mass, std::shared_ptr<Branches> data) {
  std::vector<double> dt_array(data->gpart());
  double sct, scr, mom;

  for (int event_number = 0; event_number < data->gpart(); event_number++) {
    sct = data->sc_t(event_number);
    scr = data->sc_r(event_number);
    mom = data->p(event_number);

    dt_array[event_number] = delta_t(mass, mom, sct, scr);
  }
  return dt_array;
}

void Delta_T::_delta_t_array(double mass, std::shared_ptr<Branches> data, std::vector<double> *dt_array) {
  dt_array->reserve(data->gpart());
  double sct, scr, mom;
  for (int event_number = 0; event_number < data->gpart(); event_number++) {
    sct = data->sc_t(event_number);
    scr = data->sc_r(event_number);
    mom = data->p(event_number);

    dt_array->at(event_number) = delta_t(mass, mom, sct, scr);
  }
}
