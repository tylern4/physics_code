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
  MissingMass *MM_neutron = new MissingMass(MASS_P, 0.0);
  double W, Q2;
  bool electron_cuts;

  TChain *chain = new TChain("h10");
  chain->Add(fin.c_str());
  Branches *data = new Branches(chain, true);
  int num_of_events = (int)chain->GetEntries();

  int current_event = 0;
  for (current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    electron_cuts = true;
    electron_cuts &= (data->q(0) == -1);
    electron_cuts &= (data->ec(0) > 0);
    electron_cuts &= (data->cc(0) > 0);
    electron_cuts &= (data->sc(0) > 0);
    electron_cuts &= (data->dc(0) > 0);
    if (!electron_cuts) continue;
    // Setup scattered electron 4 vector
    TLorentzVector e_mu_prime = physics::fourVec(data->p(0), data->cx(0), data->cy(0), data->cz(0), MASS_E);
    TLorentzVector gamma_mu = (*e_mu - e_mu_prime);

    TLorentzVector e_mu_prime_mc = physics::fourVec(data->pxpart(0), data->pypart(0), data->pzpart(0), MASS_E);
    TLorentzVector gamma_mu_mc = (*e_mu - e_mu_prime_mc);

    W = physics::W_calc(*e_mu, e_mu_prime);
    Q2 = physics::Q2_calc(*e_mu, e_mu_prime);
    hists->Fill_WQ2(W, Q2);

    W = physics::W_calc(*e_mu, e_mu_prime_mc);
    Q2 = physics::Q2_calc(*e_mu, e_mu_prime_mc);
    hists->Fill_WQ2_MC(W, Q2);

    TLorentzVector particle;
    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      particle =
          physics::fourVec(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num), MASS_PIP);
      MM_neutron->Set_4Vec(particle);
      MM_neutron->missing_mass(gamma_mu_mc);
      hists->Fill_Missing_Mass(MM_neutron);
    }
  }
  chain->Reset();  // delete Tree object
}

void mcHandeler::loadbar(long x, long n) {
  int w = 50;
  if ((x > n) || n == 1) return;
  if ((x != n) && (x % (n / 100 + 1) != 0)) return;

  double ratio = x / (double)n;
  int c = ratio * w;

  std::cout << BLUE << " [";
  for (int x = 0; x < c; x++) std::cout << GREEN << "=" << DEF;
  std::cout << GREEN << ">" << DEF;
  for (int x = c; x < w; x++) std::cout << " ";
  std::cout << BLUE << (int)(ratio * 100) << "%]\r" << DEF << std::flush;
}
