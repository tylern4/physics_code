/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef CSV_MAKER_H_GUARD
#define CSV_MAKER_H_GUARD
#include "main.h"
// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//
void make_electron_csv(char *fin) {
  char *csv_name_output = "511_lab.csv";
  char *root_name_output = "511_lab.root";
  const char *progress = "-\\|/";
  int num_of_events;
  bool electron_cuts;
  double _p, _cx, _cy, _cz;
  TFile *OutputFile = new TFile(root_name_output, "RECREATE");
  TTree *lab = new TTree("lab", "lab");

  lab->Branch("p", &_p);
  lab->Branch("cx", &_cx);
  lab->Branch("cy", &_cy);
  lab->Branch("cz", &_cz);
  // in main.h now
  // ofstream cut_outputs;
  csv_output.open(csv_name_output);
  csv_output << "p,cx,cy,cz" << endl;

  // Load chain from branch h10
  TChain chain("h10");
  chain.Add(fin);
  getBranches(&chain);

  num_of_events = (int)chain.GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);

    cout << "\t[ " << progress[((current_event / 100) % 4)] << " ]\t\t[" << (current_event / num_of_events) << "]\r\r"
         << flush;

    int n_prot = 0;
    int n_other = 0;
    for (int x = 0; x < gpart; x++)
      if (id[x] == PROTON)
        n_prot++;
      else
        n_other++;

    if (n_prot != 1 || n_other > 1) continue;

    electron_cuts = true;
    // electron cuts
    electron_cuts &= ((int)ec[0] > 0);                             // ``` ``` ``` ec
    electron_cuts &= ((int)id[0] == ELECTRON || (int)id[0] == 0);  // First particle is electron`
    electron_cuts &= ((int)gpart > 0);                             // Number of good particles is gt 0
    electron_cuts &= ((int)stat[0] > 0);                           // First Particle hit stat
    electron_cuts &= ((int)q[0] == -1);                            // First particle is negative Q
    electron_cuts &= ((int)sc[0] > 0);                             // First Particle hit sc
    electron_cuts &= ((int)dc[0] > 0);                             // ``` ``` ``` dc
    electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);

    // if (electron_cuts) electron_cuts &= (p[0] > MIN_P_CUT);  // Minimum Momentum cut
    if (!electron_cuts) continue;

    // Setup scattered electron 4 vector
    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)), E1D_E0);

    e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

    W = physics::W_calc(e_mu, e_mu_prime);
    Q2 = physics::Q2_calc(e_mu, e_mu_prime);

    if (W < 4.0 && Q2 < 4.0) {
      _p = (double)p[0];
      _cx = (double)cx[0];
      _cy = (double)cy[0];
      _cz = (double)cz[0];
      lab->Fill();
      csv_output << (double)p[0] << ",";
      csv_output << (double)cx[0] << ",";
      csv_output << (double)cy[0] << ",";
      csv_output << (double)cz[0];
      csv_output << endl;
    }
  }
  //
  // end stuff
  chain.Reset();
  OutputFile->cd();
  lab->Write();
  OutputFile->Close();

  csv_output.close();
}

void analyze(char *fin, const char *fout) {
  double _p, _cx, _cy, _cz;
  TH2D *hist = new TH2D("WvsQ2", "WvsQ2", 500, 0, 2, 500, 0, 4);
  TH1D *hist_W = new TH1D("W", "W", 500, 0, 2);
  TH1D *hist_Q2 = new TH1D("Q2", "Q2", 500, 0, 4);
  // Load chain from branch h10
  TChain lab("lab");
  lab.Add(fin);
  TFile *OutputFile = new TFile(fout, "RECREATE");
  lab.SetBranchAddress("p", &_p);
  lab.SetBranchAddress("cx", &_cx);
  lab.SetBranchAddress("cy", &_cy);
  lab.SetBranchAddress("cz", &_cz);

  int num_of_events = (int)lab.GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    lab.GetEntry(current_event);
    // Setup scattered electron 4 vector
    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)), E1D_E0);

    e_mu_prime_3.SetXYZ(_p * _cx, _p * _cy, _p * _cz);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

    W = physics::W_calc(e_mu, e_mu_prime);
    Q2 = physics::Q2_calc(e_mu, e_mu_prime);
    hist_W->Fill(W);
    hist_Q2->Fill(Q2);
    hist->Fill(W, Q2);
  }
  //
  // end stuff
  lab.Reset();
  OutputFile->cd();
  hist->Write();
  hist_W->Write();
  hist_Q2->Write();
  OutputFile->Close();
}
#endif
