/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "yeilds.hpp"

Yeilds::Yeilds() {}
Yeilds::Yeilds(std::string output_file_name) { csv_output.open(output_file_name); }
Yeilds::Yeilds(std::string output_file_name, bool isRoot = true) {
  Rootout = new TFile(output_file_name.c_str(), "RECREATE");
  ntuple = new TNtuple("ntuple", "", "type:W:Q2:MM:MM2:theta_e:theta_star:phi_star:theta_lab:phi_lab:sector");
}
Yeilds::~Yeilds() {
  if (ntuple) ntuple->Write();
}

void Yeilds::OpenFile(std::string output_file_name) { csv_output.open(output_file_name); }

void Yeilds::WriteHeader() {
  csv_output << "W,Q2,MM,MM2,theta_e,theta_star,phi_star,theta_lab,phi_lab,sector" << std::endl;
}

int Yeilds::Run(std::vector<std::string> fin) {
  std::cout.imbue(std::locale(""));
  int total = 0;
  int f = 0;
  for (f = 0; f < fin.size(); f++) {
    total += Run(fin.at(f));
  }
  return total;
}

int Yeilds::Run(std::string root_file) {
  auto chain = std::make_unique<TChain>("h10");
  int num_of_events = 0;
  chain->Add(root_file.c_str());
  num_of_events = (int)chain->GetEntries();
  auto data = std::make_shared<Branches>(chain.get());
  int total = 0;

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
    auto event = std::make_unique<Reaction>(data);

    if (check->isElecctron()) {
      total++;
      auto dt = std::make_unique<Delta_T>(data->sc_t(0), data->sc_r(0));
      std::vector<double> dt_proton = dt->delta_t_array(MASS_P, data);
      std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, data);

      float theta = physics::theta_calc(data->cz(0));
      float phi = physics::phi_calc(data->cx(0), data->cy(0));
      int sector = data->dc_sect(0);

      // auto photon_flux = std::make_unique<PhotonFlux>(event->e_mu(), event->e_mu_prime());
      // photon_flux->GetVirtualPhotonFlux();

      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (data->q(part_num) == POSITIVE) {
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num)))
            event->SetPip(part_num);
          else if (check->dt_P_cut(dt_proton.at(part_num), data->p(part_num)))
            event->SetProton(part_num);
          else
            event->SetOther(part_num);

        } else if (data->q(part_num) == NEGATIVE) {
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num)))
            event->SetPim(part_num);
          else
            event->SetOther(part_num);
        } else if (data->q(part_num) == 0)
          event->SetOther(part_num);
      }

      if ((event->SinglePip() || event->NeutronPip()))
        csv_output << std::setprecision(15) << event->W() << "," << event->Q2() << "," << event->MM() << ","
                   << event->MM2() << "," << event->Theta_E() << "," << event->Theta_star() << "," << event->Phi_star()
                   << "," << theta << "," << phi << "," << sector << std::endl;
    }
  }
  chain->Reset();  // delete Tree object

  return total;
}

int Yeilds::RunNtuple(std::string root_file, bool python) {
  auto chain = std::make_unique<TChain>("h10");
  int num_of_events = 0;
  chain->Add(root_file.c_str());
  num_of_events = (int)chain->GetEntries();
  auto data = std::make_shared<Branches>(chain.get());
  int total = 0;

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
    auto event = std::make_unique<Reaction>(data);

    if (check->isElecctron()) {
      total++;
      auto dt = std::make_unique<Delta_T>(data->sc_t(0), data->sc_r(0));
      std::vector<double> dt_proton = dt->delta_t_array(MASS_P, data);
      std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, data);

      float theta = physics::theta_calc(data->cz(0));
      float phi = physics::phi_calc(data->cx(0), data->cy(0));
      int sector = data->dc_sect(0);

      // auto photon_flux = std::make_unique<PhotonFlux>(event->e_mu(), event->e_mu_prime());
      // photon_flux->GetVirtualPhotonFlux();

      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        if (data->q(part_num) == POSITIVE) {
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num)))
            event->SetPip(part_num);
          else if (check->dt_P_cut(dt_proton.at(part_num), data->p(part_num)))
            event->SetProton(part_num);
          else
            event->SetOther(part_num);

        } else if (data->q(part_num) == NEGATIVE) {
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num)))
            event->SetPim(part_num);
          else
            event->SetOther(part_num);
        } else if (data->q(part_num) == 0)
          event->SetOther(part_num);
      }

      if (event->W() < 0) continue;
      if (event->Q2() > 6) continue;
      if (event->SinglePip() || event->NeutronPip() || event->ProtonPim() || event->SingleP() || event->TwoPion()) {
        ntuple->Fill(event->Type(), event->W(), event->Q2(), event->MM(), event->MM2(), event->Theta_E(),
                     event->Theta_star(), event->Phi_star(), theta, phi, sector);
      }
    }
  }
  chain->Reset();  // delete Tree object

  return total;
}
