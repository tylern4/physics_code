/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef CSV_MAKER_H_GUARD
#define CSV_MAKER_H_GUARD
#include "TFile.h"
#include "branches.hpp"
#include "main.h"
#include "physics.hpp"

double BEAM_ENERGY = 0.0;

void make_electron_csv(std::string fin) {
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  } else {
    BEAM_ENERGY = E1D_E0;
  }
  std::string csv_name_output = "511_lab_E_data.csv";
  std::string root_name_output = "511_lab_E_data.root";
  const char *progress = "-\\|/";
  int num_of_events;
  bool electron_cuts;
  double _p, _cx, _cy, _cz;
  TFile *OutputFile = new TFile(root_name_output.c_str(), "RECREATE");
  OutputFile->SetCompressionSettings(404);
  auto chain = std::make_shared<TChain>("h10");
  chain->Add(fin.c_str());
  auto data = std::make_shared<Branches>(chain);

  TTree *lab = new TTree("lab", "lab");

  lab->Branch("p", &_p);
  lab->Branch("cx", &_cx);
  lab->Branch("cy", &_cy);
  lab->Branch("cz", &_cz);
  // in main.h now
  // ofstream cut_outputs;
  csv_output.open(csv_name_output);
  csv_output << "p,cx,cy,cz" << endl;

  num_of_events = (int)chain->GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (current_event % 1000000 == 0)
      cout << "\t[ " << progress[((current_event / 1000000) % 4)] << " ]\t\t"
           << 100 * ((float)current_event / (float)num_of_events) << "\r\r" << flush;

    electron_cuts = true;
    // electron cuts
    electron_cuts &= (data->gpart() > 0);  // Number of good particles is gt 0
    electron_cuts &= (data->stat(0) > 0);  // First Particle hit stat
    electron_cuts &= (data->q(0) == -1);   // First particle is negative Q
    electron_cuts &= (data->sc(0) > 0);    // First Particle hit sc
    electron_cuts &= (data->dc(0) > 0);    // ``` ``` ``` dc
    electron_cuts &= (data->ec(0) > 0);    // ``` ``` ``` ec
    electron_cuts &= (data->dc_stat(0) > 0);

    if (!electron_cuts) continue;

    int n_prot = 0;
    int n_other = 0;
    for (int x = 1; x < data->gpart(); x++)
      if (data->id(x) == PROTON)
        n_prot++;
      else
        n_other++;

    if (n_prot == 1 && n_other == 0) {
      // Setup scattered electron 4 vector
      LorentzVector e_mu_prime(data->px(0), data->py(0), data->pz(0), MASS_E);
      LorentzVector e_mu(0.0, 0.0, sqrt(Square(BEAM_ENERGY) - Square(MASS_E)), BEAM_ENERGY);
      double W = physics::W_calc(e_mu, e_mu_prime);
      double Q2 = physics::Q2_calc(e_mu, e_mu_prime);

      if (W < 4.0 && Q2 < 4.0) {
        _p = data->p(0);
        _cx = data->cx(0);
        _cy = data->cy(0);
        _cz = data->cz(0);
        lab->Fill();
        csv_output << data->p(0) << ",";
        csv_output << data->cx(0) << ",";
        csv_output << data->cy(0) << ",";
        csv_output << data->cz(0);
        csv_output << endl;
      }
    }
  }
  //
  // end stuff
  chain->Reset();
  OutputFile->cd();
  lab->Write();
  OutputFile->Close();

  csv_output.close();
}

void make_mm_csv(std::string fin) {
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  } else {
    BEAM_ENERGY = E1D_E0;
  }
  std::string csv_name_output = "511_lab_E_PIP_data.csv";
  std::string root_name_output = "511_lab_E_PIP_data.root";
  const char *progress = "-\\|/";
  int num_of_events;
  bool electron_cuts;
  double e_p, e_cx, e_cy, e_cz;
  double pip_p, pip_cx, pip_cy, pip_cz;
  TFile *OutputFile = new TFile(root_name_output.c_str(), "RECREATE");
  OutputFile->SetCompressionSettings(404);
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
  auto chain = std::make_shared<TChain>("h10");
  chain->Add(fin.c_str());
  auto data = std::make_shared<Branches>(chain);

  num_of_events = (int)chain->GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (current_event % 10000 == 0)
      cout << "\t[ " << progress[((current_event / 10000) % 4)] << " ]\t\t"
           << 100 * ((float)current_event / (float)num_of_events) << "\r\r" << flush;

    electron_cuts = true;
    // electron cuts
    electron_cuts &= (data->gpart() > 0);  // Number of good particles is gt 0
    electron_cuts &= (data->stat(0) > 0);  // First Particle hit stat
    electron_cuts &= (data->q(0) == -1);   // First particle is negative Q
    electron_cuts &= (data->sc(0) > 0);    // First Particle hit sc
    electron_cuts &= (data->dc(0) > 0);    // ``` ``` ``` dc
    electron_cuts &= (data->ec(0) > 0);    // ``` ``` ``` ec
    electron_cuts &= (data->dc_stat(0) > 0);
    if (!electron_cuts) continue;

    int n_pion = 0;
    int n_other = 0;
    int pion_is = 0;
    for (int x = 1; x < data->gpart(); x++)
      if (data->id(x) == PIP) {
        n_pion++;
        pion_is = x;
      } else
        n_other++;

    if (electron_cuts && n_pion == 1 && n_other == 0) {
      e_p = data->p(0);
      e_cx = data->cx(0);
      e_cy = data->cy(0);
      e_cz = data->cz(0);
      pip_p = data->p(pion_is);
      pip_cx = data->cx(pion_is);
      pip_cy = data->cy(pion_is);
      pip_cz = data->cz(pion_is);
      lab->Fill();
      csv_output << data->p(0) << ",";
      csv_output << data->cx(0) << ",";
      csv_output << data->cy(0) << ",";
      csv_output << data->cz(0) << ",";
      csv_output << data->p(pion_is) << ",";
      csv_output << data->cx(pion_is) << ",";
      csv_output << data->cy(pion_is) << ",";
      csv_output << data->cz(pion_is);
      csv_output << endl;
    }
  }
  //
  // end stuff
  chain->Reset();
  OutputFile->cd();
  lab->Write();
  OutputFile->Close();

  csv_output.close();
}

void analyze_wq2(std::string fin, std::string fout) {
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  } else {
    BEAM_ENERGY = E1D_E0;
  }
  double _p, _cx, _cy, _cz;
  TH2D *hist = new TH2D("WvsQ2", "WvsQ2", 500, 0, 2, 500, 0, 4);
  TH1D *hist_W = new TH1D("W", "W", 500, 0, 2);
  TH1D *hist_Q2 = new TH1D("Q2", "Q2", 500, 0, 4);
  // Load chain from branch h10
  TChain lab("lab");
  lab.Add(fin.c_str());
  TFile *OutputFile = new TFile(fout.c_str(), "RECREATE");
  lab.SetBranchAddress("p", &_p);
  lab.SetBranchAddress("cx", &_cx);
  lab.SetBranchAddress("cy", &_cy);
  lab.SetBranchAddress("cz", &_cz);

  int num_of_events = (int)lab.GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    lab.GetEntry(current_event);
    // Setup scattered electron 4 vector
    LorentzVector e_mu_prime(_p * _cx, _p * _cy, _p * _cz, MASS_E);
    LorentzVector e_mu(0.0, 0.0, sqrt(Square(BEAM_ENERGY) - Square(MASS_E)), BEAM_ENERGY);

    double W = physics::W_calc(e_mu, e_mu_prime);
    double Q2 = physics::Q2_calc(e_mu, e_mu_prime);
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

double missing_mass(LorentzVector gamma_mu, LorentzVector pip_mu) {
  LorentzVector target(0.0, 0.0, 0.0, MASS_P);
  LorentzVector reaction;
  reaction = (gamma_mu + target - pip_mu);
  return reaction.mag();
}

double Breit(double *x, double *par) { return par[2] * TMath::BreitWigner(x[0], par[0], par[1]); }

void analyze_MM(std::string fin, std::string fout) {
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  } else {
    BEAM_ENERGY = E1D_E0;
  }
  double e_p, e_cx, e_cy, e_cz;
  double pip_p, pip_cx, pip_cy, pip_cz;
  TH2D *hist_2 = new TH2D("WvsQ2_mm", "WvsQ2_mm", 500, 0, 3, 500, 0, 4);
  TH2D *w_vs_mm = new TH2D("w_vs_mm", "w_vs_mm", 500, 0, 4, 500, -2, 2);
  TH1D *MM = new TH1D("mm", "mm", 500, 0, 3);
  TH1D *hist_W_2 = new TH1D("W_mm", "W_mm", 500, 0, 3);
  TH1D *hist_Q2_2 = new TH1D("Q2_mm", "Q2_mm", 500, 0, 4);

  TH2D *hist_2_after = new TH2D("WvsQ2_mm_after", "WvsQ2_mm_after", 500, 0, 3, 500, 0, 4);
  TH1D *MM_after = new TH1D("mm_after", "mm_after", 500, 0, 3);
  TH1D *hist_W_2_after = new TH1D("W_mm_after", "W_mm_after", 500, 0, 3);
  TH1D *hist_Q2_2_after = new TH1D("Q2_mm_after", "Q2_mm_after", 500, 0, 4);

  TF1 *bw = new TF1("bw", Breit, 0, 2, 3);

  // Load chain from branch h10
  TChain lab("lab");
  lab.Add(fin.c_str());
  TFile *OutputFile = new TFile(fout.c_str(), "RECREATE");
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

    LorentzVector e_mu_prime(e_p * e_cx, e_p * e_cy, e_p * e_cz, MASS_E);

    LorentzVector pip_mu_prime(pip_p * pip_cx, pip_p * pip_cy, pip_p * pip_cz, MASS_PIP);
    LorentzVector e_mu(0.0, 0.0, sqrt(Square(BEAM_ENERGY) - Square(MASS_E)), BEAM_ENERGY);

    LorentzVector gamma_mu;
    gamma_mu = e_mu - e_mu_prime;
    double mm = missing_mass(gamma_mu, pip_mu_prime);
    MM->Fill(mm);
    double W = physics::W_calc(e_mu, e_mu_prime);
    double Q2 = physics::Q2_calc(e_mu, e_mu_prime);
    hist_W_2->Fill(W);
    hist_Q2_2->Fill(Q2);
    hist_2->Fill(W, Q2);

    w_vs_mm->Fill(W, mm);

    if (mm < 1.1) {
      MM_after->Fill(mm);
      hist_W_2_after->Fill(W);
      hist_Q2_2_after->Fill(Q2);
      hist_2_after->Fill(W, Q2);
    }
  }
  //
  // end stuff
  lab.Reset();
  OutputFile->cd();
  hist_2->Write();
  hist_W_2->Write();
  hist_Q2_2->Write();
  w_vs_mm->Write();
  bw->SetParName(0, "Mean");
  bw->SetParName(1, "Width");
  bw->SetParName(2, "Const");
  bw->SetParameter(0, 1);
  bw->SetParameter(1, 1);
  bw->SetParameter(2, 1000);
  bw->SetParLimits(2, 100.0, 10000.0);
  MM->Fit("bw", "M", "", 0.5, 1.1);
  MM->Write();
  MM_after->Write();
  hist_2_after->Write();
  hist_W_2_after->Write();
  hist_Q2_2_after->Write();
  OutputFile->Close();
}
#endif
