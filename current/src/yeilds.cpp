/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "yeilds.hpp"

Yeilds::Yeilds() {}
Yeilds::Yeilds(std::string output_file_name) { csv_output.open(output_file_name); }
Yeilds::Yeilds(std::string output_file_name, bool isRoot = true) {
  Rootout = std::make_shared<TFile>(output_file_name.c_str(), "RECREATE");
  ntuple =
      std::make_shared<TNtuple>("ntuple", "", "type:W:Q2:MM:MM2:theta_e:theta_star:phi_star:theta_lab:phi_lab:sector");
}
Yeilds::~Yeilds() {
  if (ntuple) ntuple->Write();
}

void Yeilds::OpenFile(std::string output_file_name) { csv_output.open(output_file_name); }

void Yeilds::WriteHeader() {
  csv_output << "electron_sector,e_p,e_cx,e_cy,e_cz,pip_p,pip_cx,pip_cy,pip_cz" << std::endl;
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
  auto chain = std::make_shared<TChain>("h10");
  int num_of_events = 0;
  chain->Add(root_file.c_str());
  num_of_events = (int)chain->GetEntries();
  auto data = std::make_shared<Branches>(chain);
  int total = 0;

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
    auto event = std::make_unique<Reaction>(data);
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

    if ((event->SinglePip() || event->NeutronPip()) && pip_num != 0) {
      total++;
      csv_output << std::setprecision(15) << data->dc_sect(0) << "," << data->p(0) << "," << data->cx(0) << ","
                 << data->cy(0) << "," << data->cz(0) << "," << data->p(pip_num) << "," << data->cx(pip_num) << ","
                 << data->cy(pip_num) << "," << data->cz(pip_num) << "," << std::endl;
    }
  }

  chain->Reset();

  return total;
}

int Yeilds::RunNtuple(const std::shared_ptr<TChain> &chain) {
  size_t num_of_events = (size_t)chain->GetEntries();
  auto data = std::make_shared<Branches>(chain);
  size_t total = 0;

  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    total++;
    chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
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
    if (event->SinglePip() || event->NeutronPip() || event->ProtonPim() || event->SingleP() || event->TwoPion()) {
      event->boost();
      ntuple->Fill(event->Type(), event->W(), event->Q2(), event->MM(), event->MM2(), event->Theta_E(),
                   event->Theta_star(), event->Phi_star(), theta, phi, sector);
    }
  }
  chain->Reset();  // delete Tree object

  return total;
}
