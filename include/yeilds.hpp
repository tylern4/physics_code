/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef YEILDS_H_GUARD
#define YEILDS_H_GUARD
#include <TFile.h>
#include <TLorentzVector.h>
#include <cstring>
#include <fstream>
#include <future>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TNtuple.h"
#include "color.hpp"
#include "constants.hpp"
#include "csv_data.hpp"
#include "cuts.hpp"
#include "delta_t.hpp"
#include "histogram.hpp"
#include "missing_mass.hpp"
#include "photon_flux.hpp"
#include "physics.hpp"
#include "reaction.hpp"
#include "syncfile.hpp"

class Yeilds {
 protected:
  int PID;
  std::shared_ptr<TNtuple> ntuple = nullptr;
  std::shared_ptr<TFile> Rootout = nullptr;
  std::shared_ptr<SyncFile> _multi_threaded_csv = nullptr;
  std::shared_ptr<MomCorr> _mom_corr = nullptr;
  std::vector<std::string> input_files;
  float _beam_energy = E1D_E0;

 public:
  Yeilds();
  Yeilds(std::shared_ptr<SyncFile> multi_threaded_csv);
  Yeilds(std::shared_ptr<SyncFile> multi_threaded_csv, std::shared_ptr<MomCorr> mom_corr);
  Yeilds(std::string output_file_name, bool isRoot);
  ~Yeilds();
  void WriteData(const csv_data& toWrite);
  std::string Header();

  template <class CutType>
  int Run(std::string root_file, std::string type) {
    auto chain = std::make_shared<TChain>("h10");
    int num_of_events = 0;
    chain->Add(root_file.c_str());
    num_of_events = (int)chain->GetEntries();
    auto data = std::make_shared<Branches>(chain, false);
    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;
    int total = 0;
    for (int current_event = 0; current_event < num_of_events; current_event++) {
      chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);
      auto event = std::make_shared<Reaction>(data, _beam_energy, _mom_corr);
      int pip_num = 0;
      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (check->Pip(part_num)) {
          pip_num = part_num;
          event->SetPip(part_num);
        } else if (check->Prot(part_num)) {
          event->SetProton(part_num);
        } else if (check->Pim(part_num)) {
          event->SetPim(part_num);
        } else
          event->SetOther(part_num);
      }

      event->boost();

      if ((event->SinglePip() || event->NeutronPip()) &&
          (event->MM2() >= 0.3 && event->MM2() <= 1.5 && event->W() <= 2.0 && event->Q2() <= 4.0)) {
        total++;
        csv_data csv_buffer;
        csv_buffer.electron_sector = data->dc_sect(0);
        csv_buffer.w = event->W();
        csv_buffer.q2 = event->Q2();
        csv_buffer.theta = event->Theta_star();
        csv_buffer.phi = event->Phi_star();
        csv_buffer.mm2 = event->MM2();
        csv_buffer.helicty = data->helicity();
        csv_buffer.photon_flux = event->photon_flux();
        csv_buffer.type = type;
        csv_buffer.hash = 0;  // std::hash<std::string>{}(root_file + std::to_string(current_event));
        WriteData(csv_buffer);
      }
    }

    chain->Reset();

    return total;
  }

  template <class CutType>
  int RunNtuple(const std::shared_ptr<TChain>& chain) {
    size_t num_of_events = (size_t)chain->GetEntries();
    auto data = std::make_shared<Branches>(chain);
    size_t total = 0;

    for (size_t current_event = 0; current_event < num_of_events; current_event++) {
      total++;
      chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);
      auto event = std::make_unique<Reaction>(data);

      if (!check->isElecctron()) continue;
      // total++;
      auto dt = std::make_unique<Delta_T>(data);

      float theta = physics::theta_calc(data->cz(0));
      float phi = physics::phi_calc(data->cx(0), data->cy(0));
      int sector = data->dc_sect(0);

      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (data->q(part_num) == POSITIVE) {
          if (check->dt_Pip_cut(part_num))
            event->SetPip(part_num);
          else if (check->dt_P_cut(part_num))
            event->SetProton(part_num);
          else
            event->SetOther(part_num);
        } else if (data->q(part_num) == NEGATIVE) {
          if (check->dt_Pip_cut(part_num))
            event->SetPim(part_num);
          else
            event->SetOther(part_num);
        } else if (data->q(part_num) == 0) {
          event->SetOther(part_num);
        }
      }

      if (event->W() < 0) continue;
      if (event->Q2() > 6) continue;
      if (event->SinglePip() && event->MM() < 0) continue;

      if (event->SinglePip() || event->NeutronPip()) {
        event->boost();
        ntuple->Fill(event->Type(), event->W(), event->Q2(), event->MM(), event->MM2(), event->Theta_E(),
                     event->Theta_star(), event->Phi_star(), theta, phi, sector);
      }
    }
    chain->Reset();  // delete Tree object

    return total;
  }
};

class mcYeilds : public Yeilds {
 public:
  mcYeilds() : Yeilds(){};
  mcYeilds(std::shared_ptr<SyncFile> multi_threaded_csv) : Yeilds(multi_threaded_csv){};
  mcYeilds(std::string output_file_name, bool isRoot) : Yeilds(output_file_name, isRoot){};
  mcYeilds(std::shared_ptr<SyncFile> multi_threaded_csv, std::shared_ptr<MomCorr> mom_corr)
      : Yeilds(multi_threaded_csv, mom_corr){};

  template <class CutType>
  int RunMC(std::string root_file) {
    int total = 0;
    auto chain = std::make_shared<TChain>("h10");
    int num_of_events = 0;
    chain->Add(root_file.c_str());
    num_of_events = (int)chain->GetEntries();
    auto data = std::make_shared<Branches>(chain, true);
    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;

    for (int current_event = 0; current_event < num_of_events; current_event++) {
      chain->GetEntry(current_event);
      auto mc_event = std::make_shared<MCReaction>(data, _beam_energy);
      int pip_num = 1;
      mc_event->SetPip(pip_num);

      if (!std::isnan(mc_event->W_thrown()) && !std::isnan(mc_event->Q2_thrown())) {
        total++;
        csv_data csv_buffer;
        csv_buffer.electron_sector = data->dc_sect(0);
        csv_buffer.w = mc_event->W_thrown();
        csv_buffer.q2 = mc_event->Q2_thrown();
        csv_buffer.theta = mc_event->Theta_star();
        csv_buffer.phi = mc_event->Phi_star();
        csv_buffer.mm2 = mc_event->MM2();
        csv_buffer.helicty = data->helicity();
        csv_buffer.photon_flux = mc_event->photon_flux();
        csv_buffer.type = "thrown";
        csv_buffer.hash = 0;  // std::hash<std::string>{}(root_file + std::to_string(current_event));
        WriteData(csv_buffer);
      }
    }

    chain->Reset();
    return total;
  }
};
#endif
