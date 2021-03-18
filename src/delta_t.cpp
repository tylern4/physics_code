/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#include "delta_t.hpp"
#include "cuts.hpp"

Delta_T::Delta_T(const std::shared_ptr<Branches> &data) : _data(data) {
  _vertex = vertex_time(_data->sc_t(0), _data->sc_r(0), 1.0);

  _elec_array = _delta_t_vec(mass_map[ELECTRON]);
  _proton_array = _delta_t_vec(mass_map[PROTON]);
  _pion_array = _delta_t_vec(mass_map[PIP]);
  _kaon_array = _delta_t_vec(mass_map[KP]);
}

Delta_T::~Delta_T() {}

std::ostream &operator<<(std::ostream &os, Delta_T const &dt) {
  auto dt_cut = std::make_unique<Cuts>(dt._data);
  for (int event_number = 1; event_number < dt._data->gpart(); event_number++) {
    if (!std::isnan(dt._vertex))
      os << dt._data->dc_sect(event_number) << "," << physics::theta_rad(dt._data->cz(event_number)) << ","
         << physics::phi_rad(dt._data->cx(event_number), dt._data->cy(event_number)) << "," << dt._data->q(event_number)
         << "," << dt._vertex << "," << dt._data->p(event_number) << "," << dt._data->sc_t(event_number) << ","
         << dt._data->sc_r(event_number) << "\n";
    else {
      os << "";
    }
  }
  return os;
}

float Delta_T::vertex_time(float sc_time, float sc_pathlength, float relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

float Delta_T::Get_dt_E(int part) { return _elec_array.at(part); }
float Delta_T::Get_dt_P(int part) { return _proton_array.at(part); }
float Delta_T::Get_dt_Pi(int part) { return _pion_array.at(part); }
float Delta_T::Get_dt_K(int part) { return _kaon_array.at(part); }
float Delta_T::Get_vertex(int part) { return _vertex; }

void Delta_T::delta_t_hists(const std::shared_ptr<Histogram> &hists) {
  auto dt_cut = std::make_unique<e1d_Cuts>(_data);
  float sct, scr, mom;
  int ID, charge, sc_paddle, sc_sector;

  for (int event_number = 1; event_number < _data->gpart(); event_number++) {
    sct = _data->sc_t(event_number);
    scr = _data->sc_r(event_number);
    mom = _data->p(event_number);
    ID = _data->id(event_number);
    charge = _data->q(event_number);
    sc_paddle = _data->sc_pd(event_number);
    sc_sector = _data->sc_sect(event_number);

    if (charge == POSITIVE) {
      hists->Fill_deltat_P(mom, _proton_array.at(event_number));
      hists->Fill_deltat_PIP(mom, _pion_array.at(event_number));
      hists->Fill_deltat_kp(mom, _kaon_array.at(event_number));
      if (dt_cut->dt_Pip_cut(event_number))
        hists->Fill_deltat_PIP_PID(mom, _pion_array.at(event_number));
      else if (dt_cut->dt_P_cut(event_number))
        hists->Fill_deltat_P_PID(mom, _proton_array.at(event_number));
      else if (dt_cut->dt_K_cut(event_number))
        hists->Fill_deltat_kp_PID(mom, _kaon_array.at(event_number));
    } else {
      hists->Fill_deltat_electron(mom, _elec_array.at(event_number));
      hists->Fill_deltat_PIM(mom, _pion_array.at(event_number));
      if (dt_cut->dt_Pip_cut(event_number)) hists->Fill_deltat_PIM_PID(mom, _pion_array.at(event_number));
    }

    hists->delta_t_Fill(mom, charge, _proton_array.at(event_number), _pion_array.at(event_number),
                        _elec_array.at(event_number));
    if (sc_sector <= 6 && sc_sector != 0 && sc_paddle <= 48 && sc_paddle != 0) {
      hists->delta_t_sec_pad(mom, charge, _proton_array.at(event_number), _pion_array.at(event_number),
                             _elec_array.at(event_number), sc_sector, sc_paddle);
    }
  }
}

float Delta_T::_delta_t(float mass, float momentum, float sc_t, float sc_r) {
  float cut_beta = 1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
  return _vertex - vertex_time(sc_t, sc_r, cut_beta);
}

std::vector<float> Delta_T::_delta_t_vec(float mass) {
  std::vector<float> dt_array(_data->gpart());
  float sct, scr, mom;

  for (int event_number = 0; event_number < _data->gpart(); event_number++) {
    sct = _data->sc_t(event_number);
    scr = _data->sc_r(event_number);
    mom = _data->p(event_number);

    dt_array[event_number] = _delta_t(mass, mom, sct, scr);
  }
  return dt_array;
}
