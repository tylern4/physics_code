/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
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
#include "cuts.hpp"
#include "delta_t.hpp"
#include "histogram.hpp"
#include "missing_mass.hpp"
#include "photon_flux.hpp"
#include "physics.hpp"
#include "reaction.hpp"

class Yeilds {
 protected:
  int PID;
  std::shared_ptr<TNtuple> ntuple = nullptr;
  std::shared_ptr<TFile> Rootout = nullptr;
  std::ofstream csv_output;
  std::vector<std::string> input_files;
  double _beam_energy = E1D_E0;

 public:
  Yeilds();
  Yeilds(std::string output_file_name);
  Yeilds(std::string output_file_name, bool isRoot);
  ~Yeilds();
  void OpenFile(std::string output_file_name);
  void WriteHeader();

  template <class CutType>
  int Run(std::string root_file, std::string type) {
    auto chain = std::make_shared<TChain>("h10");
    int num_of_events = 0;
    chain->Add(root_file.c_str());
    num_of_events = (int)chain->GetEntries();
    auto data = std::make_shared<Branches>(chain);
    int total = 0;
    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;

    for (int current_event = 0; current_event < num_of_events; current_event++) {
      chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);
      auto event = std::make_shared<Reaction>(data, _beam_energy);
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

      if (event->SinglePip() || event->NeutronPip()) {
        total++;
        std::lock_guard<std::mutex> lk(std::mutex);
        csv_output << std::setprecision(15) << data->dc_sect(0) << "," << event->W() << "," << event->Q2() << ","
                   << event->Theta_star() << "," << event->Phi_star() << "," << event->MM2() << "," << data->p(0) << ","
                   << data->cx(0) << "," << data->cy(0) << "," << data->cz(0) << "," << data->p(pip_num) << ","
                   << data->cx(pip_num) << "," << data->cy(pip_num) << "," << data->cz(pip_num) << ","
                   << data->helicity() << "," << type << std::endl;
      }
    }

    chain->Reset();

    return total;
  }

  template <class CutType>
  int RunNtuple(const std::shared_ptr<TChain> &chain) {
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

  // int Run(std::vector<std::string> fin);
};

class mcYeilds : public Yeilds {
 public:
  mcYeilds() : Yeilds(){};
  mcYeilds(std::string output_file_name) : Yeilds(output_file_name){};
  mcYeilds(std::string output_file_name, bool isRoot) : Yeilds(output_file_name, isRoot){};

  template <class CutType>
  int Run(std::string root_file, std::string type) {
    Yeilds::Run<CutType>(root_file, "mc_rec");

    auto chain = std::make_shared<TChain>("h10");
    int num_of_events = 0;
    chain->Add(root_file.c_str());
    num_of_events = (int)chain->GetEntries();
    auto data = std::make_shared<BranchesMC>(chain);
    int total = 0;
    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;

    for (int current_event = 0; current_event < num_of_events; current_event++) {
      chain->GetEntry(current_event);
      auto mc_event = std::make_shared<MCReaction>(data, _beam_energy);
      int pip_num = 0;
      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (data->pidpart(part_num) == PIP) {
          pip_num = part_num;
          mc_event->SetPip(part_num);
        }
      }

      mc_event->boost();

      if (mc_event->SinglePip()) {
        total++;
        std::lock_guard<std::mutex> lk(std::mutex);
        csv_output << std::setprecision(15) << data->dc_sect(0) << "," << mc_event->W() << "," << mc_event->Q2() << ","
                   << mc_event->Theta_star() << "," << mc_event->Phi_star() << "," << mc_event->MM2() << ","
                   << data->p(0) << "," << data->cx(0) << "," << data->cy(0) << "," << data->cz(0) << ","
                   << data->p(pip_num) << "," << data->cx(pip_num) << "," << data->cy(pip_num) << ","
                   << data->cz(pip_num) << "," << data->helicity() << "," << type << std::endl;
      }
    }

    chain->Reset();

    return total;
  }
};
#endif
