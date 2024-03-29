/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "ntuple.hpp"

Ntuple::Ntuple(const std::string &output_file_name) {
  rootout = std::make_shared<TFile>(output_file_name.c_str(), "RECREATE");
  rootout->SetCompressionAlgorithm(404);
  ntuple = std::make_shared<TTree>("events", "events");
  ntuple->Branch("type", &_type);
  ntuple->Branch("helicity", &_helicity);
  ntuple->Branch("W", &_W);
  ntuple->Branch("Q2", &_Q2);
  ntuple->Branch("MM", &_MM);
  ntuple->Branch("MM2", &_MM2);
  ntuple->Branch("Theta_E", &_Theta_E);
  ntuple->Branch("Theta_star", &_Theta_star);
  ntuple->Branch("Phi_star", &_Phi_star);
  ntuple->Branch("theta", &_theta);
  ntuple->Branch("phi", &_phi);
  ntuple->Branch("sector", &_sector);
}

Ntuple::Ntuple(const std::vector<std::string> &infiles, const std::string &output_file_name) {
  ROOT::EnableThreadSafety();
  rootout = std::make_shared<TFile>(output_file_name.c_str(), "RECREATE");
  rootout->SetCompressionAlgorithm(404);
  ntuple = std::make_shared<TTree>("events", "events");
  ntuple->Branch("type", &_type);
  ntuple->Branch("helicity", &_helicity);
  ntuple->Branch("W", &_W);
  ntuple->Branch("Q2", &_Q2);
  ntuple->Branch("MM", &_MM);
  ntuple->Branch("MM2", &_MM2);
  ntuple->Branch("Theta_E", &_Theta_E);
  ntuple->Branch("Theta_star", &_Theta_star);
  ntuple->Branch("Phi_star", &_Phi_star);
  ntuple->Branch("theta", &_theta);
  ntuple->Branch("phi", &_phi);
  ntuple->Branch("sector", &_sector);
}

Ntuple::~Ntuple() {
  if (ntuple) ntuple->Write();
}

size_t Ntuple::Run(const std::shared_ptr<TChain> &chain) {
  size_t num_of_events = (size_t)chain->GetEntries();
  auto data = std::make_shared<Branches>(chain);
  size_t total = 0;

  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    if (current_event % 100000 == 0) std::cerr << "\t" << 100.0 * current_event / num_of_events << "\r\r";
    chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
    auto event = std::make_unique<Reaction>(data);

    if (!check->isElectron()) continue;

    auto dt = std::make_unique<Delta_T>(data);

    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      if (check->Pip(part_num))
        event->SetPip(part_num);
      else if (check->Prot(part_num))
        event->SetProton(part_num);
      else if (check->Prot(part_num))
        event->SetPim(part_num);
      else
        event->SetOther(part_num);
    }

    if (event->W() < 0) continue;
    if (event->W() > 3) continue;
    if (event->Q2() < 0) continue;
    if (event->Q2() > 5) continue;
    if (event->SinglePip() && event->MM() <= 0) continue;

    event->boost();

    _type = event->Type();
    _helicity = data->helicity();
    _W = event->W();
    _Q2 = event->Q2();
    _MM = event->MM();
    _MM2 = event->MM2();
    _Theta_E = event->Theta_E();
    _Theta_star = event->Theta_star();
    _Phi_star = event->Phi_star();
    _theta = physics::theta_calc(data->cz(0));
    _phi = physics::phi_calc(data->cx(0), data->cy(0));
    _sector = event->sector();

    ntuple->Fill();

    total++;
  }
  chain->Reset();  // delete Tree object

  return total;
}
