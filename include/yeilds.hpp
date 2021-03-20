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

// 10,000 Ten Thousand
#define WRITEBUFF 10000

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
  void Save();
  std::string Header();

  template <class CutType>
  int Run(const std::shared_ptr<TChain>& chain, const std::string& type, const size_t& thread_id) {
    auto data = std::make_shared<Branches>(chain, false);
    size_t num_of_events = (size_t)chain->GetEntries();

    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;
    size_t total = 0;
    size_t written = 0;
    for (size_t current_event = 0; current_event < num_of_events; current_event++) {
      if (thread_id == 0 && current_event % 100000 == 0)
        std::cerr << "\t" << 100 * (current_event / (float)num_of_events) + 1 << "\r\r" << std::flush;

      chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);

      total++;

      auto event = std::make_shared<Reaction>(data, _beam_energy, _mom_corr);
      if (event->Q2() < 1.1 || event->Q2() > 3.5) continue;

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
      if (cosf(event->Theta_star()) == -1.0) continue;

      bool cut_fid = (check->isElectron() && check->fid_chern_cut());
      short electron_sector = data->dc_sect(0);
      if (electron_sector < 0 || electron_sector > 6) electron_sector = 0;

      // bool cut_angles =
      //     (cos(event->Theta_star()) == -1.0 && event->Phi_star() >= 1.57079 && event->Phi_star() <= 1.57082) ||
      //     (cos(event->Theta_star()) == -1.0 && event->Phi_star() >= 4.71238 && event->Phi_star() <= 4.71240);
      // bool cut_angles = ;

      if ((event->SinglePip() || event->NeutronPip()) && cut_fid &&
          (event->MM2() >= 0.3 && event->MM2() <= 1.5 && event->W() <= 2.0 && event->Q2() <= 4.0)) {
        written++;
        csv_data csv_buffer;
        csv_buffer.electron_sector = electron_sector;
        csv_buffer.w = event->W();
        csv_buffer.q2 = event->Q2();
        csv_buffer.theta = event->Theta_star();
        csv_buffer.phi = event->Phi_star();
        csv_buffer.mm2 = event->MM2();
        csv_buffer.cut_fid = cut_fid;
        csv_buffer.helicty = data->helicity();
        csv_buffer.type = type;
        // csv_buffer.hash = 0;  // std::hash<std::string>{}(root_file + std::to_string(current_event));
        WriteData(csv_buffer);
      }
      if (thread_id == 0 && written % WRITEBUFF == 0) Save();
    }

    // chain->Reset();

    return total;
  }

  template <class CutType>
  int RunPPi0(const std::shared_ptr<TChain>& chain, const std::string& type, const size_t& thread_id) {
    auto data = std::make_shared<Branches>(chain, false);
    size_t num_of_events = (size_t)chain->GetEntries();

    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;
    int total = 0;
    for (size_t current_event = 0; current_event < num_of_events; current_event++) {
      if (thread_id == 0 && current_event % 100000 == 0)
        std::cerr << "\t" << 100 * (current_event / (float)num_of_events) + 1 << "\r\r" << std::flush;

      chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);

      total++;

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
      short electron_sector = data->dc_sect(0);
      if (electron_sector < 0 || electron_sector > 6) electron_sector = 0;
      if (event->Q2() < 1.1 || event->Q2() > 3.5) continue;

      bool cut_fid = (check->isElectron() && check->fid_chern_cut());

      if (event->PPi0()) {
        csv_data csv_buffer;
        csv_buffer.electron_sector = electron_sector;
        csv_buffer.w = event->W();
        csv_buffer.q2 = event->Q2();
        csv_buffer.theta = event->Theta_star();
        csv_buffer.phi = event->Phi_star();
        csv_buffer.mm2 = event->pi0_mass2();
        csv_buffer.cut_fid = cut_fid;
        csv_buffer.helicty = data->helicity();
        csv_buffer.type = type;
        WriteData(csv_buffer);
      }
    }

    // chain->Reset();

    return total;
  }

  template <class CutType>
  int RunElastic(const std::shared_ptr<TChain>& chain, const std::string& type, const size_t& thread_id) {
    auto data = std::make_shared<Branches>(chain, false);
    size_t num_of_events = (size_t)chain->GetEntries();

    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;
    int total = 0;
    for (size_t current_event = 0; current_event < num_of_events; current_event++) {
      if (thread_id == 0 && current_event % 100000 == 0)
        std::cerr << "\t" << 100 * (current_event / (float)num_of_events) + 1 << "\r\r" << std::flush;

      chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);

      total++;

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
      short electron_sector = data->dc_sect(0);
      if (electron_sector < 0 || electron_sector > 6) electron_sector = 0;
      // if (event->Q2() < 1.1 || event->Q2() > 3.5) continue;

      bool cut_fid = (check->isElectron() && check->fid_chern_cut());

      if (event->elastic()) {
        csv_data csv_buffer;
        csv_buffer.electron_sector = electron_sector;
        csv_buffer.w = event->W();
        csv_buffer.q2 = event->Q2();
        csv_buffer.theta = event->Theta_star();
        csv_buffer.phi = event->Phi_star();
        csv_buffer.mm2 = event->MM2();
        csv_buffer.cut_fid = cut_fid;
        csv_buffer.helicty = data->helicity();
        csv_buffer.type = type;
        WriteData(csv_buffer);
      }
    }

    // chain->Reset();

    return total;
  }

  template <class CutType>
  int RunNtuple(const std::shared_ptr<TChain>& chain) {
    size_t num_of_events = (size_t)chain->GetEntries();
    auto data = std::make_shared<Branches>(chain);
    size_t total = 0;

    for (size_t current_event = 0; current_event < num_of_events; current_event++) {
      chain->GetEntry(current_event);
      auto check = std::make_unique<CutType>(data);
      auto event = std::make_unique<Reaction>(data);

      if (!check->isElectron()) continue;
      total++;
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
      if (event->Q2() < 1.1 || event->Q2() > 3.5) continue;
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

  template <class CutType>
  size_t Run_fidFits(const std::shared_ptr<TChain>& chain, const size_t& thread_id) {
    auto start = std::chrono::high_resolution_clock::now();
    auto data = std::make_shared<Branches>(chain, false);
    size_t num_of_events = (size_t)chain->GetEntries();

    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;
    size_t total = 0;
    size_t current_event = 0;
    for (current_event = 0; current_event < num_of_events; current_event++) {
      if (thread_id == 0 && current_event % 100000 == 0)
        std::cerr << "\t" << 100 * (current_event / (float)num_of_events) + 1 << "\r\r" << std::flush;

      chain->GetEntry(current_event);
      if (std::isnan(data->p(0))) continue;
      auto event = std::make_shared<Reaction>(data, _beam_energy, nullptr);
      if (event->cc_x() == 0 || event->cc_y() == 0) continue;
      auto check = std::make_unique<CutType>(data);

      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (check->Pip(part_num)) {
          event->SetPip(part_num);
        } else if (check->Prot(part_num)) {
          event->SetProton(part_num);
        } else if (check->Pim(part_num)) {
          event->SetPim(part_num);
        } else
          event->SetOther(part_num);
      }

      std::string type = "other";
      if (event->channel()) type = "channel";
      if (event->PPi0()) type = "PPi0";

      total++;
      fid_data csv_buffer;
      csv_buffer.sector = data->dc_sect(0);
      csv_buffer.e_p = data->p(0);
      csv_buffer.e_theta = physics::theta_calc_rad(data->cz(0));
      csv_buffer.e_phi = physics::phi_calc_rad(data->cx(0), data->cy(0));
      csv_buffer.theta = physics::theta_calc_rad(data->cz(0));
      csv_buffer.phi = physics::center_phi_calc_rad(data->cx(0), data->cy(0));
      csv_buffer.x = event->cc_x();
      csv_buffer.y = event->cc_y();
      csv_buffer.type = type;

      _multi_threaded_csv->write(csv_buffer.print());
      // std::cout << csv_buffer << std::endl;
    }

    chain->Reset();

    return total;
  }

  template <class CutType>
  size_t Run_fidFitsPip(const std::shared_ptr<TChain>& chain, const size_t& thread_id) {
    auto start = std::chrono::high_resolution_clock::now();
    auto data = std::make_shared<Branches>(chain, false);
    size_t num_of_events = (size_t)chain->GetEntries();

    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;
    size_t total = 0;
    size_t current_event = 0;
    for (current_event = 0; current_event < num_of_events; current_event++) {
      if (thread_id == 0 && current_event % 100000 == 0)
        std::cerr << "\t" << 100 * (current_event / (float)num_of_events) + 1 << "\r\r" << std::flush;

      chain->GetEntry(current_event);

      short pip_num = -1;
      auto check = std::make_unique<CutType>(data);
      if (!check->isElectron()) continue;
      if (!check->fid_chern_cut()) continue;

      auto event = std::make_shared<Reaction>(data, _beam_energy, nullptr);
      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (check->Pip(part_num)) {
          pip_num = part_num;
          event->SetPip(part_num);
          continue;
        }
      }

      short e_sector = data->dc_sect(0);
      if (e_sector == 0 || e_sector > NUM_SECTORS || e_sector == int(NULL)) continue;

      // if (event->W() < 1.0 || event->W() > 1.8) continue;

      event->boost();

      total++;
      fid_data_pip csv_buffer;
      csv_buffer.e_sector = data->dc_sect(0);
      csv_buffer.e_p = data->p(0);
      csv_buffer.e_theta = physics::theta_calc_rad(data->cz(0));
      csv_buffer.e_phi = physics::phi_calc_rad(data->cx(0), data->cy(0));
      if (pip_num != -1) {
        short pip_sector = data->dc_sect(pip_num);
        if (pip_sector == 0 || pip_sector > NUM_SECTORS || pip_sector == int(NULL)) continue;
        if (data->p(pip_num) == 0) continue;
        csv_buffer.pip_sector = pip_sector;
        csv_buffer.pip_p = data->p(pip_num);
        csv_buffer.pip_theta = physics::theta_calc_rad(data->cz(pip_num));
        csv_buffer.pip_phi = physics::phi_calc_rad(data->cx(pip_num), data->cy(pip_num));
        csv_buffer.pip_theta_star = event->Theta_star();
        csv_buffer.pip_phi_star = event->Phi_star();
      } else {
        csv_buffer.pip_sector = -1;
        csv_buffer.pip_p = NAN;
        csv_buffer.pip_theta = NAN;
        csv_buffer.pip_phi = NAN;
        csv_buffer.pip_theta_star = NAN;
        csv_buffer.pip_phi_star = NAN;
      }

      _multi_threaded_csv->write(csv_buffer.print());
    }

    chain->Reset();

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
  int RunMC(const std::shared_ptr<TChain>& chain, const size_t& thread_id) {
    auto start = std::chrono::high_resolution_clock::now();
    auto data = std::make_shared<Branches>(chain, true);
    size_t num_of_events = (size_t)chain->GetEntries();
    // PRINT_TIMEING(start, "Got number of events from chain " << num_of_events << " : ");

    // auto data = std::make_shared<Branches>(chain, true);
    _beam_energy = std::is_same<CutType, e1f_Cuts>::value ? E1F_E0 : E1D_E0;

    int total = 0;
    for (size_t current_event = 0; current_event < num_of_events; current_event++) {
      if (thread_id == 0 && current_event % 100000 == 0)
        std::cerr << "\t" << 100 * (current_event / (float)num_of_events) + 1 << "\r\r" << std::flush;
      chain->GetEntry(current_event);
      auto mc_event = std::make_shared<MCReaction>(data, _beam_energy);
      if (mc_event->Q2_thrown() < 1.1 || mc_event->Q2_thrown() > 3.5) continue;

      int pip_num = 1;
      mc_event->SetPip(pip_num);

      short electron_sector = data->dc_sect(0);
      if (electron_sector < 0 || electron_sector > 6) electron_sector = 0;

      if (!std::isnan(mc_event->W_thrown()) && !std::isnan(mc_event->Q2_thrown())) {
        total++;
        csv_data csv_buffer;
        csv_buffer.electron_sector = electron_sector;
        csv_buffer.w = mc_event->W_thrown();
        csv_buffer.q2 = mc_event->Q2_thrown();
        csv_buffer.theta = mc_event->Theta_star();
        csv_buffer.phi = mc_event->Phi_star();
        csv_buffer.mm2 = mc_event->MM2();
        csv_buffer.cut_fid = true;
        csv_buffer.helicty = data->helicity();
        csv_buffer.type = "thrown";
        // csv_buffer.hash = 0;  // std::hash<std::string>{}(root_file + std::to_string(current_event));
        WriteData(csv_buffer);
      }
      if (thread_id == 0 && total % WRITEBUFF == 0) Save();
    }

    // chain->Reset();
    return total;
  }
};
#endif
