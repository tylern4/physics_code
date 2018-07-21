#include "skim.hpp"

void Skim::getSkimBranches() {
  chain->SetBranchAddress("gpart", &gpart);
  chain->SetBranchAddress("id", id);
  chain->SetBranchAddress("stat", stat);
  chain->SetBranchAddress("dc", dc);
  chain->SetBranchAddress("cc", cc);
  chain->SetBranchAddress("sc", sc);
  chain->SetBranchAddress("ec", ec);
  chain->SetBranchAddress("p", p);
  chain->SetBranchAddress("m", m);
  chain->SetBranchAddress("q", q);
  chain->SetBranchAddress("b", b);
  chain->SetBranchAddress("cx", cx);
  chain->SetBranchAddress("cy", cy);
  chain->SetBranchAddress("cz", cz);
  chain->SetBranchAddress("vx", vx);
  chain->SetBranchAddress("vy", vy);
  chain->SetBranchAddress("vz", vz);
  chain->SetBranchAddress("dc_stat", dc_stat);
  chain->SetBranchAddress("etot", etot);
  chain->SetBranchStatus("*", 1);
}

Skim::Skim(std::string input) {
  fin = input;
  chain = new TChain("h10");
  chain->AddFile(fin.c_str());

  fout = fin.substr(0, fin.size() - 5) + "_skim.root";
  RootOutputFile = new TFile(fout.c_str(), "RECREATE");

  e_mu.SetPxPyPzE(0.0, 0.0, sqrt((E1D_E0 * E1D_E0) - (MASS_E * MASS_E)), E1D_E0);
}
Skim::~Skim() {}

void Skim::Process() {
  int num_of_events;
  bool electron_cuts;
  std::cout << BLUE << "Analyzing file " << GREEN << fin << DEF << std::endl;
  getSkimBranches();
  num_of_events = (int)chain->GetEntries();
  TTree *skim = chain->CloneTree(0);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);

    electron_cuts = true;
    // electron cuts
    electron_cuts &= (id[0] == ELECTRON);       // First particle is electron
    electron_cuts &= (gpart > 0);               // Number of good particles is greater than 0
    electron_cuts &= (stat[0] > 0);             // First Particle hit stat
    electron_cuts &= ((int)q[0] == -1);         // First particle is negative Q
    electron_cuts &= (sc[0] > 0);               // First Particle hit sc
    electron_cuts &= (dc[0] > 0);               // ``` ``` ``` d
    electron_cuts &= (ec[0] > 0);               // ``` ``` ``` ec
    electron_cuts &= (dc_stat[dc[0] - 1] > 0);  //??
    if (electron_cuts) {
      electron_cuts &= (etot[ec[0] - 1] / p[0]) < 0.4;
      electron_cuts &= (etot[ec[0] - 1] / p[0]) > 0.2;
    }

    // e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
    // e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

    if (electron_cuts) {
      skim->Fill();  // Fill the banks after the skim
    }
  }
  chain->Reset();  // delete Tree object
  delete chain;

  RootOutputFile->cd();
  RootOutputFile->Write();
  RootOutputFile->Close();
  delete RootOutputFile;
}
