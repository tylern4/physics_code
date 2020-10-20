/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
#include <TFile.h>
#include <cstring>
#include <fstream>
#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include "TChain.h"
#include "color.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "delta_t.hpp"
#include "histogram.hpp"
#include "missing_mass.hpp"
#include "mom_corr.hpp"
#include "photon_flux.hpp"
#include "physics.hpp"
#include "reaction.hpp"

class DataHandeler {
 protected:
  bool _loadbar = false;
  std::shared_ptr<TChain> _chain;
  std::vector<std::string> _input_files;
  std::shared_ptr<Histogram> _hists;
  std::shared_ptr<Branches> _data;
  std::shared_ptr<MomCorr> _mom_corr = nullptr;
  float _beam_energy = NAN;

 public:
  DataHandeler();
  DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists,
               const std::shared_ptr<MomCorr>& mom_corr);
  DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists);
  ~DataHandeler();
  void setLoadBar(bool load);

  template <class CutType>
  void RunEvent(size_t current_event) {
    _chain->GetEntry(current_event);
    auto check = std::make_unique<CutType>(_data);
    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;
    _hists->Fill_Beam_Position(_data);

    if (!check->isElecctron()) return;

    _hists->EC_inout(_data->ec_ei(0), _data->ec_eo(0));
    _hists->EC_fill(_data->etot(0), _data->p(0));
    _hists->TM_Fill(_data->p(0), physics::theta_calc(_data->cz(0)));
    float theta_cc = TMath::ACos(TMath::Abs(_data->p(0) * _data->cz(0)) / TMath::Abs(_data->p(0))) / D2R;
    _hists->CC_fill(_data->cc_sect(0), (_data->cc_segm(0) % 1000) / 10, _data->cc_segm(0) / 1000 - 1, _data->nphe(0),
                    theta_cc);

    auto event = std::make_shared<Reaction>(_data, _beam_energy, _mom_corr);

    _hists->Fill_E_Prime_fid(event->e_mu_prime());
    _hists->Fill_E_Prime(event->e_mu_prime());

    if (!check->isElecctron()) return;

    auto photon_flux = std::make_unique<PhotonFlux>(event->e_mu(), event->e_mu_prime());
    _hists->Photon_flux_Fill(photon_flux->GetVirtualPhotonFlux());
    _hists->WvsQ2_Fill(event->W(), event->Q2(), _data->dc_sect(0));

    auto dt = check->share_dt();  // std::make_unique<Delta_T>(_data);
    dt->delta_t_hists(_hists);

    float theta = physics::theta_calc(_data->cz(0));
    float phi = physics::phi_calc(_data->cx(0), _data->cy(0));
    int elec_sector = _data->dc_sect(0);
    int sector = 0;

    _hists->Fill_electron_fid(_data, event);

    for (int part_num = 1; part_num < _data->gpart(); part_num++) {
      theta = physics::theta_calc(_data->cz(part_num));
      phi = physics::phi_calc(_data->cx(part_num), _data->cy(part_num));
      sector = _data->dc_sect(part_num);

      _hists->delta_t_sec_pad(_data->p(part_num), _data->q(part_num), dt->Get_dt_P(part_num), dt->Get_dt_Pi(part_num),
                              dt->Get_dt_E(part_num), _data->sc_sect(part_num), _data->sc_pd(part_num));

      _hists->Fill_Target_Vertex(_data->cx(part_num), _data->cy(part_num), _data->cz(part_num), _data->vx(part_num),
                                 _data->vy(part_num), _data->vz(part_num));
      _hists->MomVsBeta_Fill(_data->p(part_num), _data->b(part_num));

      if (_data->q(part_num) == POSITIVE) {
        _hists->MomVsBeta_Fill_pos(_data->p(part_num), _data->b(part_num));
      } else if (_data->q(part_num) == NEGATIVE) {
        _hists->MomVsBeta_Fill_neg(_data->p(part_num), _data->b(part_num));
      } else if (_data->q(part_num) == 0 && _data->id(part_num) == NEUTRON) {
        _hists->MomVsBeta_Fill_neutral(_data->p(part_num), _data->b(part_num));
        _hists->Fill_neutron_fid(_data->cc_c2(part_num), _data->cc_r(part_num), _data->cc_sect(part_num));
      }

      if (check->Pip(part_num) && check->Prot(part_num))
        _hists->Fill_proton_Pi_ID_P(_data->p(part_num), _data->b(part_num));

      if (check->Pip(part_num)) {
        event->SetPip(part_num);
        _hists->Fill_hadron_fid(theta, phi, sector, PIP);
        _hists->Fill_pion_WQ2(event->W(), event->Q2());
        _hists->Fill_Pi_ID_P(_data->p(part_num), _data->b(part_num));
      } else if (check->Prot(part_num)) {
        event->SetProton(part_num);
        _hists->Fill_hadron_fid(theta, phi, sector, PROTON);
        _hists->Fill_proton_WQ2(event->W(), event->Q2());
        _hists->Fill_proton_ID_P(_data->p(part_num), _data->b(part_num));
      } else if (check->Pim(part_num)) {
        _hists->Fill_hadron_fid(theta, phi, sector, PIM);
        event->SetPim(part_num);
      } else
        event->SetOther(part_num);
    }

    if (event->channel()) _hists->EC_cut_fill(_data->etot(0), _data->p(0));
    _hists->FillEvent(event);

    return;
  }

  template <class CutType>
  int Run() {
    if (_hists == nullptr) return 0;
    size_t num_of_events = (size_t)_chain->GetEntries();
    for (size_t current_event = 0; current_event < num_of_events; current_event++) {
      if (_loadbar && current_event % 10000 == 0) DataHandeler::loadbar(current_event, num_of_events);
      DataHandeler::RunEvent<CutType>(current_event);
    }

    _chain->Reset();
    return num_of_events;
  }

  void loadbar(long x, long n);
};

class mcHandeler : public DataHandeler {
 protected:
  std::shared_ptr<BranchesMC> _data_mc;
  std::shared_ptr<mcHistogram> _mc_hists;

 public:
  mcHandeler() : DataHandeler() {}
  mcHandeler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists);
  mcHandeler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists,
             const std::shared_ptr<MomCorr>& mom_corr);

  template <class CutType>
  void RunEvent(int current_event) {
    // _chain->GetEntry(current_event);
    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;
    auto mc_event = std::make_shared<MCReaction>(_data, _beam_energy);
    _mc_hists->Fill_P(_data);
    _mc_hists->Fill_WQ2_MC(mc_event);
    _hists->Fill_ND(mc_event);
    _mc_hists->Fill(mc_event);

    for (int part_num = 1; part_num < _data->gpart(); part_num++) {
      if (_data->pidpart(part_num) == PIP || _data->pidpart(part_num) == 8) mc_event->SetPip(part_num);
    }
    if (mc_event->SinglePip()) _mc_hists->Fill_Missing_Mass(mc_event);

    return;
  }

  template <class CutType>
  int Run() {
    size_t num_of_events = (size_t)_chain->GetEntries();

    for (size_t current_event = 0; current_event < num_of_events; current_event++) {
      if (_loadbar && current_event % 10000 == 0) DataHandeler::loadbar(current_event, num_of_events);
      DataHandeler::RunEvent<CutType>(current_event);
      mcHandeler::RunEvent<CutType>(current_event);
    }

    _chain->Reset();
    return num_of_events;
  }
};

#endif
