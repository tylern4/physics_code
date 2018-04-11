/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef CSV_MAKER_H_GUARD
#define CSV_MAKER_H_GUARD
#include "main.h"
#include "TLorentzVector.h"
// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//

void make_electron_csv(char *fin) {
  char *csv_name_output = "511_lab_E_data.csv";
  char *root_name_output = "511_lab_E_data.root";
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
    if (current_event % 100 == 0)
      cout << "\t[ " << progress[((current_event / 100) % 4)] << " ]\t\t["
           << 100 * ((float)current_event / (float)num_of_events) << "]\r\r" << flush;

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

void make_mm_csv(char *fin) {
  char *csv_name_output = "511_lab_E_PIP_data.csv";
  char *root_name_output = "511_lab_E_PIP_data.root";
  const char *progress = "-\\|/";
  int num_of_events;
  bool electron_cuts;
  double e_p, e_cx, e_cy, e_cz;
  double pip_p, pip_cx, pip_cy, pip_cz;
  TFile *OutputFile = new TFile(root_name_output, "RECREATE");
  TTree *lab = new TTree("lab", "lab");

  lab->Branch("e_p", &e_p);
  lab->Branch("e_cx", &e_cx);
  lab->Branch("e_cy", &e_cy);
  lab->Branch("e_cz", &e_cz);
  lab->Branch("pip_p", &pip_p);
  lab->Branch("pip_cx", &pip_cx);
  lab->Branch("pip_cy", &pip_cy);
  lab->Branch("pip_cz", &pip_cz);

  // in main.h now
  // ofstream cut_outputs;
  csv_output.open(csv_name_output);
  csv_output << "e_p,e_cx,e_cy,e_cz,pip_p,pip_cx,pip_cy,pip_cz" << endl;

  // Load chain from branch h10
  TChain chain("h10");
  chain.Add(fin);
  getBranches(&chain);

  num_of_events = (int)chain.GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    if (current_event % 100 == 0)
      cout << "\t[ " << progress[((current_event / 100) % 4)] << " ]\t\t["
           << 100 * ((float)current_event / (float)num_of_events) << "]\r\r" << flush;

    int n_pion = 0;
    int n_other = 0;
    int pion_is = 0;
    for (int x = 0; x < gpart; x++)
      if (id[x] == PIP) {
        n_pion++;
        pion_is = x;
      } else
        n_other++;

    if (n_pion != 1 || n_other > 1) continue;

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
    if (electron_cuts) {
      e_p = (double)p[0];
      e_cx = (double)cx[0];
      e_cy = (double)cy[0];
      e_cz = (double)cz[0];
      pip_p = (double)p[pion_is];
      pip_cx = (double)cx[pion_is];
      pip_cy = (double)cy[pion_is];
      pip_cz = (double)cz[pion_is];
      lab->Fill();
      csv_output << (double)p[0] << ",";
      csv_output << (double)cx[0] << ",";
      csv_output << (double)cy[0] << ",";
      csv_output << (double)cz[0] << ",";
      csv_output << (double)p[pion_is] << ",";
      csv_output << (double)cx[pion_is] << ",";
      csv_output << (double)cy[pion_is] << ",";
      csv_output << (double)cz[pion_is];
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

void analyze_wq2(char *fin, const char *fout) {
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

double missing_mass(TLorentzVector gamma_mu, TLorentzVector pip_mu) {
  TVector3 target_3;
  TLorentzVector target;
  // Set target vector
  target_3.SetXYZ(0.0, 0.0, 0.0);
  target.SetVectM(target_3, MASS_P);
  TLorentzVector reaction;
  reaction = (gamma_mu + target - pip_mu);

  return reaction.M();
}

double Breit(double *x, double *par) { return par[2] * TMath::BreitWigner(x[0], par[0], par[1]); }

void analyze_MM(char *fin, const char *fout) {
  double e_p, e_cx, e_cy, e_cz;
  double pip_p, pip_cx, pip_cy, pip_cz;
  TH2D *hist_2 = new TH2D("WvsQ2_mm", "WvsQ2_mm", 500, 0, 3, 500, 0, 4);
  TH1D *MM = new TH1D("mm", "mm", 500, 0, 3);
  TH1D *hist_W_2 = new TH1D("W_mm", "W_mm", 500, 0, 3);
  TH1D *hist_Q2_2 = new TH1D("Q2_mm", "Q2_mm", 500, 0, 4);
  TF1 *bw = new TF1("bw", Breit, 0, 2, 3);

  // Load chain from branch h10
  TChain lab("lab");
  lab.Add(fin);
  TFile *OutputFile = new TFile(fout, "RECREATE");
  lab.SetBranchAddress("e_p", &e_p);
  lab.SetBranchAddress("e_cx", &e_cx);
  lab.SetBranchAddress("e_cy", &e_cy);
  lab.SetBranchAddress("e_cz", &e_cz);
  lab.SetBranchAddress("pip_p", &pip_p);
  lab.SetBranchAddress("pip_cx", &pip_cx);
  lab.SetBranchAddress("pip_cy", &pip_cy);
  lab.SetBranchAddress("pip_cz", &pip_cz);

  int num_of_events = (int)lab.GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    lab.GetEntry(current_event);
    // Setup scattered electron 4 vector
    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TVector3 pip_mu_prime_3;
    TLorentzVector pip_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)), E1D_E0);

    e_mu_prime_3.SetXYZ(e_p * e_cx, e_p * e_cy, e_p * e_cz);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

    W = physics::W_calc(e_mu, e_mu_prime);
    Q2 = physics::Q2_calc(e_mu, e_mu_prime);
    hist_W_2->Fill(W);
    hist_Q2_2->Fill(Q2);
    hist_2->Fill(W, Q2);

    pip_mu_prime_3.SetXYZ(pip_p * pip_cx, pip_p * pip_cy, pip_p * pip_cz);
    pip_mu_prime.SetVectM(pip_mu_prime_3, MASS_PIP);

    TLorentzVector gamma_mu;
    gamma_mu = e_mu - e_mu_prime;
    double mm = missing_mass(gamma_mu, pip_mu_prime);
    MM->Fill(mm);
  }
  //
  // end stuff
  lab.Reset();
  OutputFile->cd();
  hist_2->Write();
  hist_W_2->Write();
  hist_Q2_2->Write();
  bw->SetParName(0, "Mean");
  bw->SetParName(1, "Width");
  bw->SetParName(2, "Const");
  bw->SetParameter(0, 1);
  bw->SetParameter(1, 1);
  bw->SetParameter(2, 1);
  MM->Fit("bw", "M", "", 0.8, 1.0);
  MM->Write();
  OutputFile->Close();
}
#endif
