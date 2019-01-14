/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "datahandeler_mc.hpp"

mcHandeler::mcHandeler() {
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  } else {
    BEAM_ENERGY = E1D_E0;
  }
  e_mu = new TLorentzVector(0.0, 0.0, sqrt(BEAM_ENERGY * BEAM_ENERGY - MASS_E * MASS_E), BEAM_ENERGY);
}
mcHandeler::~mcHandeler() {}

void mcHandeler::Run(std::vector<std::string> fin, mcHistogram *hists) {
  int i = 1;
  for (auto f : fin) {
    this->loadbar(i++, fin.size());
    Run(f, hists);
  }
}

void mcHandeler::Run(std::string fin, mcHistogram *hists) {
  auto *chain = new TChain("h10");
  chain->Add(fin.c_str());
  auto *data = new Branches(chain, true);
  int num_of_events = (int)chain->GetEntries();

  int current_event = 0;
  for (current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto check = std::make_shared<Cuts>(data);
    if (!check->isElecctron()) continue;

    auto event = std::make_shared<Reaction>();
    auto mc_event = std::make_shared<Reaction>();
    event->SetElec(data->p(0), data->cx(0), data->cy(0), data->cz(0));
    mc_event->SetElec(data->p(0), data->cx(0), data->cy(0), data->cz(0));
    hists->Fill_P(data);
    hists->Fill_WQ2(event->W(), event->Q2());
    hists->Fill_WQ2_MC(mc_event->W(), mc_event->Q2());

    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      if (data->pidpart(part_num) == PIP)
        event->SetPip(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num));
    }
    hists->Fill_Missing_Mass(event->MM(), event->MM2());
  }
  chain->Reset();  // delete Tree object
}

void mcHandeler::loadbar(long x, long n) {
  int w = 50;
  if ((x > n) || n == 1) return;
  if ((x != n) && (x % (n / 100 + 1) != 0)) return;

  double ratio = x / (double)n;
  int c = (int)ratio * w;

  std::cout << BLUE << " [";
  for (int z = 0; z < c; z++) std::cout << GREEN << "=" << DEF;
  std::cout << GREEN << ">" << DEF;
  for (int y = c; y < w; y++) std::cout << " ";
  std::cout << BLUE << (int)(ratio * 100) << "%]\r" << DEF << std::flush;
}
