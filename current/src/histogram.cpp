/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "histogram.hpp"

// using namespace std;

Histogram::Histogram() {
  makeHists_delta_t();
  makeHists_CC();
  makeHists_fid();
}

Histogram::~Histogram() {
  delete WvsQ2_hist;
  delete W_hist;
  delete Q2_hist;
  delete E_prime_hist;
  delete Q2_vs_xb;
  delete WvsQ2_proton;
  delete W_proton;
  delete Q2_proton;
  delete WvsQ2_pion;
  delete W_pion;
  delete Q2_pion;
  delete WvsQ2_single_pi;
  delete W_single_pi;
  delete Q2_single_pi;
  delete WvsQ2_single_proton;
  delete W_single_proton;
  delete Q2_single_proton;
  delete MomVsBeta_hist;
  delete MomVsBeta_hist_pos;
  delete MomVsBeta_hist_neg;
  delete Mom;
  delete Energy_hist;
  delete MomVsBeta_proton_ID;
  delete MomVsBeta_Pi_ID;
  delete MomVsBeta_proton_Pi_ID;
  delete delta_t_mass_P;
  delete delta_t_mass_P_PID;
  delete delta_t_mass_PIP;
  delete delta_t_mass_PIP_PID;
  delete delta_t_mass_electron;
  delete delta_t_mass_electron_PID;
  delete delta_t_mass_positron;
  delete delta_t_mass_positron_PID;
  delete fid_hist;
  delete EC_sampling_fraction;
  delete Missing_Mass;
  delete Missing_Mass_square;
  delete Theta_CC;
}

// W and Q^2
void Histogram::Fill_proton_WQ2(double W, double Q2) {
  WvsQ2_proton->Fill(W, Q2);
  W_proton->Fill(W);
  Q2_proton->Fill(Q2);
}

void Histogram::Fill_single_pi_WQ2(double W, double Q2) {
  WvsQ2_single_pi->Fill(W, Q2);
  W_single_pi->Fill(W);
  Q2_single_pi->Fill(Q2);
}

void Histogram::Fill_single_proton_WQ2(double W, double Q2) {
  WvsQ2_single_proton->Fill(W, Q2);
  W_single_proton->Fill(W);
  Q2_single_proton->Fill(Q2);
}

void Histogram::WvsQ2_Fill(double E_prime, double W, double Q2, double xb) {
  E_prime_hist->Fill(E_prime);
  WvsQ2_hist->Fill(W, Q2);
  W_hist->Fill(W);
  Q2_hist->Fill(Q2);
  Q2_vs_xb->Fill(xb, Q2);
}

void Histogram::Fill_pion_WQ2(double W, double Q2) {
  WvsQ2_pion->Fill(W, Q2);
  W_pion->Fill(W);
  Q2_pion->Fill(Q2);
}

void Histogram::WvsQ2_Write() {
  WvsQ2_hist->SetXTitle("W (GeV)");
  WvsQ2_hist->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_hist->Write();

  W_hist->SetXTitle("W (GeV)");
  W_hist->Write();

  Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
  Q2_hist->Write();

  E_prime_hist->SetXTitle("Energy (GeV)");
  E_prime_hist->Write();

  Q2_vs_xb->SetXTitle("x_{b}");
  Q2_vs_xb->SetYTitle("Q^{2}");
  Q2_vs_xb->Write();

  WvsQ2_proton->SetXTitle("W (GeV)");
  WvsQ2_proton->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_proton->Write();

  W_proton->SetXTitle("W (GeV)");
  W_proton->Write();

  Q2_proton->SetXTitle("Q^{2} (GeV^{2})");
  Q2_proton->Write();

  WvsQ2_pion->SetXTitle("W (GeV)");
  WvsQ2_pion->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_pion->Write();

  W_pion->SetXTitle("W (GeV)");
  W_pion->Write();

  Q2_pion->SetXTitle("Q^{2} (GeV^{2})");
  Q2_pion->Write();

  WvsQ2_single_pi->SetXTitle("W (GeV)");
  WvsQ2_single_pi->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_single_pi->Write();

  W_single_pi->SetXTitle("W (GeV)");
  W_single_pi->Write();

  Q2_single_pi->SetXTitle("Q^{2} (GeV^{2})");
  Q2_single_pi->Write();

  WvsQ2_single_proton->SetXTitle("W (GeV)");
  WvsQ2_single_proton->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_single_proton->Write();

  W_single_proton->SetXTitle("W (GeV)");
  W_single_proton->Write();

  Q2_single_proton->SetXTitle("Q^{2} (GeV^{2})");
  Q2_single_proton->Write();
}
// W and Q^2

// P and E
void Histogram::MomVsBeta_Fill_pos(double P, double Beta) {
  MomVsBeta_hist_pos->Fill(P, Beta);
}

void Histogram::MomVsBeta_Fill_neg(double P, double Beta) {
  MomVsBeta_hist_neg->Fill(P, Beta);
}

void Histogram::Fill_proton_ID_P(double p, double beta) {
  MomVsBeta_proton_ID->Fill(p, beta);
}

void Histogram::Fill_Pi_ID_P(double p, double beta) {
  MomVsBeta_Pi_ID->Fill(p, beta);
}

void Histogram::Fill_proton_Pi_ID_P(double p, double beta) {
  MomVsBeta_proton_Pi_ID->Fill(p, beta);
}

void Histogram::MomVsBeta_Fill(double Energy, double P, double Beta) {
  Energy_hist->Fill(Energy);
  MomVsBeta_hist->Fill(P, Beta);
  Mom->Fill(P);
}

void Histogram::MomVsBeta_Write() {
  MomVsBeta_hist->SetXTitle("Momentum (GeV)");
  MomVsBeta_hist->SetYTitle("#beta");
  MomVsBeta_hist_pos->SetXTitle("Momentum (GeV)");
  MomVsBeta_hist_pos->SetYTitle("#beta");
  MomVsBeta_hist_neg->SetXTitle("Momentum (GeV)");
  MomVsBeta_hist_neg->SetYTitle("#beta");
  Mom->SetXTitle("Momentum (GeV)");

  MomVsBeta_proton_ID->SetXTitle("Momentum (GeV)");
  MomVsBeta_proton_ID->SetYTitle("#beta");
  MomVsBeta_proton_ID->Write();

  MomVsBeta_Pi_ID->SetXTitle("Momentum (GeV)");
  MomVsBeta_Pi_ID->SetYTitle("#beta");
  MomVsBeta_Pi_ID->Write();

  MomVsBeta_proton_Pi_ID->SetXTitle("Momentum (GeV)");
  MomVsBeta_proton_Pi_ID->SetYTitle("#beta");
  MomVsBeta_proton_Pi_ID->Write();

  Energy_hist->Write();
  MomVsBeta_hist->Write();
  MomVsBeta_hist_pos->Write();
  MomVsBeta_hist_neg->Write();
  Mom->Write();
}

// Missing Mass
void Histogram::Fill_Missing_Mass(double miss_mass) {
  Missing_Mass->Fill(miss_mass);
}

void Histogram::Fill_Mass(double mass) { Mass->Fill(mass); }

void Histogram::Fill_Missing_Mass_square(double miss_mass_2) {
  Missing_Mass_square->Fill(miss_mass_2);
}

void Histogram::Write_Missing_Mass() {
  Missing_Mass->SetXTitle("Mass (GeV)");
  Missing_Mass->Write();

  Mass->SetXTitle("Mass (GeV)");
  Mass->Write();

  Missing_Mass_square->SetXTitle("Mass (GeV)");
  Missing_Mass_square->Write();
}

void Histogram::makeHists_delta_t() {
  for (int jj = 0; jj < num_points; jj++) {
    sprintf(hname, "delta_t_p_%d", jj);
    sprintf(htitle, "#Deltat P %d", jj);
    delta_t_hist[0][jj] = new TH1D(hname, htitle, bins, Dt_min, Dt_max);

    sprintf(hname, "delta_t_pip_%d", jj);
    sprintf(htitle, "#Deltat #pi^{+} %d", jj);
    delta_t_hist[1][jj] = new TH1D(hname, htitle, bins, Dt_min, Dt_max);

    sprintf(hname, "delta_t_electron_%d", jj);
    sprintf(htitle, "#Deltat electron %d", jj);
    delta_t_hist[2][jj] = new TH1D(hname, htitle, bins, Dt_min, Dt_max);
  }

  for (int jj = 0; jj < sc_sector_num; jj++) {
    for (int jjj = 0; jjj < sc_paddle_num; jjj++) {
      sprintf(hname, "delta_t_p_sec%d_pad%d", jj + 1, jjj + 1);
      sprintf(htitle, "#Deltat P Sector %d Paddle %d", jj + 1, jjj + 1);
      delta_t_sec_pad_hist[0][jj][jjj] = new TH2D(
          hname, htitle, bins / 2, p_min, p_max, bins / 2, Dt_min, Dt_max);

      sprintf(hname, "delta_t_pip_sec%d_pad%d", jj + 1, jjj + 1);
      sprintf(htitle, "#Deltat #pi^{+} Sector %d Paddle %d", jj + 1, jjj + 1);
      delta_t_sec_pad_hist[1][jj][jjj] = new TH2D(
          hname, htitle, bins / 2, p_min, p_max, bins / 2, Dt_min, Dt_max);

      sprintf(hname, "delta_t_electron_sec%d_pad%d", jj + 1, jjj + 1);
      sprintf(htitle, "#Deltat electron Sector %d Paddle %d", jj + 1, jjj + 1);
      delta_t_sec_pad_hist[2][jj][jjj] = new TH2D(
          hname, htitle, bins / 2, p_min, p_max, bins / 2, Dt_min, Dt_max);
    }
  }
}

void Histogram::Fill_deltat_P(double momentum, double delta_t) {
  delta_t_mass_P->Fill(momentum, delta_t);
}

void Histogram::Fill_deltat_P_PID(double momentum, double delta_t) {
  delta_t_mass_P_PID->Fill(momentum, delta_t);
}

void Histogram::Fill_deltat_PIP(double momentum, double delta_t) {
  delta_t_mass_PIP->Fill(momentum, delta_t);
}

void Histogram::Fill_deltat_PIP_PID(double momentum, double delta_t) {
  delta_t_mass_PIP_PID->Fill(momentum, delta_t);
}

void Histogram::Fill_deltat_electron(double momentum, double delta_t) {
  delta_t_mass_electron->Fill(momentum, delta_t);
}

void Histogram::Fill_deltat_electron_PID(double momentum, double delta_t) {
  delta_t_mass_electron_PID->Fill(momentum, delta_t);
}

void Histogram::Fill_deltat_positron(double momentum, double delta_t) {
  delta_t_mass_positron->Fill(momentum, delta_t);
}

void Histogram::Fill_deltat_positron_PID(double momentum, double delta_t) {
  delta_t_mass_positron_PID->Fill(momentum, delta_t);
}

void Histogram::delta_t_slice_fit() {
  Header *fit_functions = new Header("../src/fit_functions.hpp", "FF");

  TF1 *peak = new TF1("peak", "gaus", -1, 1);
  //[0]*exp(-[1]*x) +
  char *func = "[0]*exp(-[1]*x) + [2]*x*x + [3]*x + [4]";
  delta_t_mass_P->FitSlicesY(peak, 0, -1, 10, "QRG5");
  TH1D *delta_t_mass_P_0 = (TH1D *)gDirectory->Get("delta_t_mass_P_0");
  TH1D *delta_t_mass_P_1 = (TH1D *)gDirectory->Get("delta_t_mass_P_1");
  TH1D *delta_t_mass_P_2 = (TH1D *)gDirectory->Get("delta_t_mass_P_2");
  double x[500];
  double y_plus[500];
  double y_minus[500];
  int num = 0;
  for (int i = 0; i < 500; i++) {
    if (delta_t_mass_P_1->GetBinContent(i) != 0) {
      // Get momentum from bin center
      x[num] = (double)delta_t_mass_P_1->GetBinCenter(i);
      // mean + 3sigma
      y_plus[num] = (double)delta_t_mass_P_1->GetBinContent(i) +
                    N_SIGMA * (double)delta_t_mass_P_2->GetBinContent(i);
      // mean - 3simga
      y_minus[num] = (double)delta_t_mass_P_1->GetBinContent(i) -
                     N_SIGMA * (double)delta_t_mass_P_2->GetBinContent(i);
      num++;
    }
  }

  TGraph *P = new TGraph(num, x, y_plus);
  TGraph *M = new TGraph(num, x, y_minus);
  TF1 *Proton_Pos_fit = new TF1("Proton_Pos_fit", func);
  TF1 *Proton_Neg_fit = new TF1("Proton_Neg_fit", func);
  P->Fit(Proton_Pos_fit, "QRG5", "", 0.2, 2);
  P->Write();
  M->Fit(Proton_Neg_fit, "QRG5", "", 0.2, 2);
  M->Write();
  Proton_Pos_fit->Write();
  Proton_Neg_fit->Write();
  P->Draw("Same");
  M->Draw("Same");
  Proton_Pos_fit->Draw("Same");
  Proton_Neg_fit->Draw("Same");

  fit_functions->NewFunction();
  fit_functions->Set_RetrunType("double");
  fit_functions->Set_FuncName("Proton_Pos_fit");
  fit_functions->Set_FuncInputs("double x");
  fit_functions->Set_Function(Proton_Pos_fit->GetExpFormula("P"));
  fit_functions->WriteFunction();

  fit_functions->NewFunction();
  fit_functions->Set_RetrunType("double");
  fit_functions->Set_FuncName("Proton_Neg_fit");
  fit_functions->Set_FuncInputs("double x");
  fit_functions->Set_Function(Proton_Neg_fit->GetExpFormula("P"));
  fit_functions->WriteFunction();

  delta_t_mass_PIP->FitSlicesY(peak, 0, -1, 10, "QRG5");
  TH1D *delta_t_mass_PIP_0 = (TH1D *)gDirectory->Get("delta_t_mass_PIP_0");
  TH1D *delta_t_mass_PIP_1 = (TH1D *)gDirectory->Get("delta_t_mass_PIP_1");
  TH1D *delta_t_mass_PIP_2 = (TH1D *)gDirectory->Get("delta_t_mass_PIP_2");
  double x_pip[500];
  double y_plus_pip[500];
  double y_minus_pip[500];
  num = 0;
  for (int i = 0; i < 500; i++) {
    if (delta_t_mass_PIP_1->GetBinContent(i) != 0) {
      // Get momentum from bin center
      x_pip[num] = (double)delta_t_mass_PIP_1->GetBinCenter(i);
      // mean + 3sigma
      y_plus_pip[num] = (double)delta_t_mass_PIP_1->GetBinContent(i) +
                        N_SIGMA * (double)delta_t_mass_PIP_2->GetBinContent(i);
      // mean - 3simga
      y_minus_pip[num] = (double)delta_t_mass_PIP_1->GetBinContent(i) -
                         N_SIGMA * (double)delta_t_mass_PIP_2->GetBinContent(i);
      num++;
    }
  }

  TGraph *P_pip = new TGraph(num, x_pip, y_plus_pip);
  TGraph *M_pip = new TGraph(num, x_pip, y_minus_pip);
  TF1 *Pip_Pos_fit = new TF1("Pip_Pos_fit", func);
  TF1 *Pip_Neg_fit = new TF1("Pip_Neg_fit", func);
  P_pip->Fit(Pip_Pos_fit, "QRG5", "", 0.1, 1.75);
  P_pip->Write();
  M_pip->Fit(Pip_Neg_fit, "QRG5", "", 0.1, 1.75);
  M_pip->Write();
  Pip_Pos_fit->Write();
  Pip_Neg_fit->Write();
  P_pip->Draw("Same");
  M_pip->Draw("Same");
  Pip_Pos_fit->Draw("Same");
  Pip_Neg_fit->Draw("Same");

  fit_functions->NewFunction();
  fit_functions->Set_RetrunType("double");
  fit_functions->Set_FuncName("Pip_Pos_fit");
  fit_functions->Set_FuncInputs("double x");
  fit_functions->Set_Function(Pip_Pos_fit->GetExpFormula("P"));
  fit_functions->WriteFunction();

  fit_functions->NewFunction();
  fit_functions->Set_RetrunType("double");
  fit_functions->Set_FuncName("Pip_Neg_fit");
  fit_functions->Set_FuncInputs("double x");
  fit_functions->Set_Function(Pip_Neg_fit->GetExpFormula("P"));
  fit_functions->WriteFunction();

  fit_functions->NewFunction();
  fit_functions->Set_RetrunType("bool");
  fit_functions->Set_FuncName("Between_Pip_fit");
  fit_functions->Set_FuncInputs("double dt, double p");
  fit_functions->AddLine("bool between = true");
  fit_functions->AddLine("between &= (dt >= Pip_Neg_fit(p))");
  fit_functions->AddLine("between &= (dt <= Pip_Pos_fit(p))");
  fit_functions->Set_Function("between");
  fit_functions->WriteFunction();

  fit_functions->NewFunction();
  fit_functions->Set_RetrunType("bool");
  fit_functions->Set_FuncName("Between_Proton_fit");
  fit_functions->Set_FuncInputs("double dt, double p");
  fit_functions->AddLine("bool between = true");
  fit_functions->AddLine("between &= (dt >= Proton_Neg_fit(p))");
  fit_functions->AddLine("between &= (dt <= Proton_Pos_fit(p))");
  fit_functions->Set_Function("between");
  fit_functions->WriteFunction();

  delete fit_functions;
}

void Histogram::delta_t_Write() {
  delta_t_mass_P->SetXTitle("Momentum (GeV)");
  delta_t_mass_P->SetYTitle("#Deltat");
  delta_t_mass_PIP->SetXTitle("Momentum (GeV)");
  delta_t_mass_PIP->SetYTitle("#Deltat");
  delta_t_mass_P_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_P_PID->SetYTitle("#Deltat");
  delta_t_mass_PIP_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_PIP_PID->SetYTitle("#Deltat");
  delta_t_mass_electron->SetXTitle("Momentum (GeV)");
  delta_t_mass_electron->SetYTitle("#Deltat");
  delta_t_mass_electron_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_electron_PID->SetYTitle("#Deltat");
  delta_t_mass_positron->SetXTitle("Momentum (GeV)");
  delta_t_mass_positron->SetYTitle("#Deltat");
  delta_t_mass_positron_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_positron_PID->SetYTitle("#Deltat");

  delta_t_slice_fit();

  delta_t_mass_P->Write();
  delta_t_mass_P_PID->Write();

  delta_t_mass_PIP->Write();
  delta_t_mass_PIP_PID->Write();
  delta_t_mass_electron->Write();
  delta_t_mass_electron_PID->Write();
  delta_t_mass_positron->Write();
  delta_t_mass_positron_PID->Write();
}

void Histogram::delta_t_Fill(double momentum, int charge, double delta_t_proton,
                             double delta_t_pip, double delta_t_electron) {
  for (int jj = 0; jj < num_points; jj++) {
    if (momentum > jj * bin_width && momentum <= (jj + 1) * bin_width) {
      if (charge == 1 && !std::isnan(delta_t_proton) &&
          !std::isnan(delta_t_pip)) {
        delta_t_hist[0][jj]->Fill(delta_t_proton);
        delta_t_hist[1][jj]->Fill(delta_t_pip);
      }
      if (charge == -1)
        delta_t_hist[2][jj]->Fill(delta_t_electron);
    }
  }
}

void Histogram::delta_t_slices_Write() {
  Fits delta_t_cut[3][num_points];
  double fit_dt_min = -1.0;
  double fit_dt_max = 1.0;
  for (int j = 0; j < 3; j++) {
    for (int jj = 0; jj < num_points; jj++) {
      if (j != 2) {
        delta_t_cut[j][num_points].FitGaus(delta_t_hist[j][jj], fit_dt_min,
                                           fit_dt_max);
        // cout << j << ',' << jj << ',' << delta_t_cut[j][num_points].mean <<
        // ',' << delta_t_cut[j][num_points].sigma << std::endl;
      }

      delta_t_hist[j][jj]->SetYTitle("#Deltat");
      delta_t_hist[j][jj]->Write();
    }
  }
}

void Histogram::delta_t_sec_pad(double momentum, int charge,
                                double delta_t_proton, double delta_t_pip,
                                double delta_t_electron, int sc_sector,
                                int sc_paddle) {
  if (charge == 1) {
    delta_t_sec_pad_hist[0][sc_sector - 1][sc_paddle - 1]->Fill(momentum,
                                                                delta_t_proton);
    delta_t_sec_pad_hist[1][sc_sector - 1][sc_paddle - 1]->Fill(momentum,
                                                                delta_t_pip);
  }
  if (charge == -1)
    delta_t_sec_pad_hist[2][sc_sector - 1][sc_paddle - 1]->Fill(
        momentum, delta_t_electron);
}

void Histogram::delta_t_sec_pad_Write() {
  for (int j = 0; j < 3; j++) {
    for (int jj = 0; jj < sc_sector_num; jj++) {
      for (int jjj = 0; jjj < sc_paddle_num; jjj++) {
        delta_t_sec_pad_hist[j][jj][jjj]->SetYTitle("#Deltat");
        delta_t_sec_pad_hist[j][jj][jjj]->Write();
      }
    }
  }
}

void Histogram::delta_T_canvas() {
  TCanvas *can_dt[sc_sector_num][3];
  char can_name[50];
  char *P_PIP_E;
  for (int particle_i = 0; particle_i < 3; particle_i++) {
    for (int sec_i = 0; sec_i < sc_sector_num; sec_i++) {
      if (particle_i == 0)
        P_PIP_E = "Proton";
      else if (particle_i == 1)
        P_PIP_E = "Pip";
      else if (particle_i == 2)
        P_PIP_E = "Electron";

      sprintf(can_name, "Sector %d %s", sec_i + 1, P_PIP_E);
      can_dt[sec_i][particle_i] = new TCanvas(can_name, can_name, 1200, 800);
      can_dt[sec_i][particle_i]->Divide(6, 8);
      for (int pad_i = 0; pad_i < sc_paddle_num; pad_i++) {
        can_dt[sec_i][particle_i]->cd((int)pad_i + 1);
        delta_t_sec_pad_hist[particle_i][sec_i][pad_i]->Draw("same"
                                                             "colz");
      }
      can_dt[sec_i][particle_i]->Write();
    }
  }
}

void Histogram::CC_fill(int cc_sector, int cc_segment, int cc_pmt, int cc_nphe,
                        double theta_cc) {
  x_cc_sparse[0] = cc_sector;
  x_cc_sparse[1] = cc_segment;
  x_cc_sparse[2] = cc_pmt;
  x_cc_sparse[3] = cc_nphe;
  cc_sparse->Fill(x_cc_sparse);

  if (cc_pmt == -1)
    cc_pmt = 2;
  cc_hist[cc_sector - 1][cc_segment - 1][cc_pmt]->Fill(cc_nphe);
  cc_hist_allSeg[cc_sector - 1][cc_pmt]->Fill(cc_nphe);

  Theta_CC->Fill(cc_segment, theta_cc);
}

void Histogram::makeHists_CC() {
  cc_sparse->GetAxis(0)->SetTitle(" cc_sector ");
  cc_sparse->GetAxis(1)->SetTitle(" cc_segment ");
  cc_sparse->GetAxis(2)->SetTitle(" cc_pmt ");
  cc_sparse->GetAxis(3)->SetTitle(" cc_nphe ");

  for (int sec_i = 0; sec_i < sector; sec_i++) {
    for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
      if (pmt_i == 0)
        L_R_C = "both";
      if (pmt_i == 1)
        L_R_C = "right";
      if (pmt_i == 2)
        L_R_C = "left";
      sprintf(hname, "CC_sec%d_%s", sec_i + 1, L_R_C);
      sprintf(htitle, "CC sector %d %s", sec_i + 1, L_R_C);
      cc_hist_allSeg[sec_i][pmt_i] =
          new TH1D(hname, htitle, bins_CC, CC_min, CC_max);
      for (int seg_i = 0; seg_i < segment; seg_i++) {
        sprintf(hname, "CC_sec%d_seg%d_%s", sec_i + 1, seg_i + 1, L_R_C);
        sprintf(htitle, "CC sector %d segment %d %s", sec_i + 1, seg_i + 1,
                L_R_C);
        cc_hist[sec_i][seg_i][pmt_i] =
            new TH1D(hname, htitle, bins_CC, CC_min, CC_max);
      }
    }
  }
}

void Histogram::CC_Write() {
  Theta_CC->SetXTitle("CC segment");
  Theta_CC->SetYTitle("#theta_CC");
  Theta_CC->Write();

  cc_sparse->Write();

  for (int sec_i = 0; sec_i < sector; sec_i++) {
    for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
      cc_hist_allSeg[sec_i][pmt_i]->SetYTitle("number photoelectrons");
      cc_hist_allSeg[sec_i][pmt_i]->Write();
      for (int seg_i = 0; seg_i < segment; seg_i++) {
        cc_hist[sec_i][seg_i][pmt_i]->SetYTitle("number photoelectrons");
        cc_hist[sec_i][seg_i][pmt_i]->Write();
      }
    }
  }
}

void Histogram::CC_canvas() {
  TCanvas *can[sector][PMT];
  char can_name[50];
  for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
    for (int sec_i = 0; sec_i < sector; sec_i++) {
      if (pmt_i == 0)
        L_R_C = "both";
      if (pmt_i == 1)
        L_R_C = "right";
      if (pmt_i == 2)
        L_R_C = "left";

      sprintf(can_name, "Sector %d %s", sec_i + 1, L_R_C);
      can[sec_i][pmt_i] = new TCanvas(can_name, can_name, 1200, 800);
      can[sec_i][pmt_i]->Divide(6, 3);
      for (int seg_i = 0; seg_i < segment; seg_i++) {
        can[sec_i][pmt_i]->cd((int)seg_i + 1);
        cc_hist[sec_i][seg_i][pmt_i]->Draw("same");
      }
      can[sec_i][pmt_i]->Write();
    }
  }
}

void Histogram::makeHists_fid() {
  fid_sec_hist.reserve(sector_num);
  for (int sec = 0; sec < sector_num; sec++) {
    sprintf(hname, "fid_sec%d", sec + 1);
    sprintf(htitle, "fid_sec%d", sec + 1);
    fid_sec_hist[sec] = new TH2D(hname, htitle, bins, min_phi[sec],
                                 max_phi[sec], bins, theta_min, theta_max);
  }
}

void Histogram::Fill_fid(double theta, double phi, int sector) {
  fid_hist->Fill(phi, theta);
  fid_sec_hist[sector]->Fill(phi, theta);
}

void Histogram::Fid_Write() {
  fid_hist->SetYTitle("#theta");
  fid_hist->SetXTitle("#phi");

  fid_hist->Write();
  Fits fid_sec_lo[sector_num];
  Fits fid_sec_hi[sector_num];
  Fits fid_sec[sector_num];

  for (int sec = 0; sec < sector_num; sec++) {
    // fid_sec_lo[sec].FitFiducial_lo(fid_sec_hist[sec], min_phi[sec],
    //                               max_phi[sec]);
    // fid_sec[sec].FitPoly_fid(fid_sec_hist[sec], min_phi[sec], max_phi[sec]);
    // fid_sec_hi[sec].FitFiducial_hi(fid_sec_hist[sec], min_phi[sec],
    //                               max_phi[sec]);
    fid_sec_hist[sec]->SetYTitle("#theta");
    fid_sec_hist[sec]->SetXTitle("#phi");
    fid_sec_hist[sec]->Write();
  }
}

void Histogram::EC_fill(double etot, double momentum) {
  double sampling_frac = etot / momentum;
  EC_sampling_fraction->Fill(momentum, sampling_frac);
}

void Histogram::EC_Write() {
  EC_sampling_fraction->SetXTitle("Momentum (GeV)");
  EC_sampling_fraction->SetYTitle("Sampling Fraction");

  EC_sampling_fraction->Write();
}

void Histogram::Fill_Beam_Position(double vertex_x, double vertex_y) {
  Beam_Position->Fill(vertex_x, vertex_y);
  Beam_Position_X->Fill(vertex_x);
  Beam_Position_Y->Fill(vertex_y);
}

void Histogram::Beam_Position_Write() {
  Beam_Position->SetXTitle("X");
  Beam_Position->SetYTitle("Y");
  Beam_Position->Write();

  Beam_Position_X->SetXTitle("X");
  Beam_Position_X->Write();

  Beam_Position_Y->SetXTitle("Y");
  Beam_Position_Y->Write();
}
