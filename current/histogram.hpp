/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HISTOGRAM_H
#define HISTOGRAM_H
#include "TH2.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TGraph.h"
#include "TF1.h"

class Histogram {
private:
  void makeHists_fid();
  void makeHists_delta_t();
  void makeHists_CC();
  const int bins = 500;
  const double p_min = 0.0;
  const double p_max = 5.0;
  char hname[50];
  char htitle[500];

  // W and Q^2
  double w_min = 0;
  double w_max = 3.25;
  double q2_min = 0;
  double q2_max = 10;

  TH2D *WvsQ2_hist = new TH2D("WvsQ2_hist", "W vs Q^{2}", bins, w_min, w_max,
                              bins, q2_min, q2_max);
  TH1D *W_hist = new TH1D("W", "W", bins, w_min, w_max);
  TH1D *Q2_hist = new TH1D("Q2", "Q^{2}", bins, q2_min, q2_max);
  TH1D *E_prime_hist =
      new TH1D("E_prime", "Scattered Electron Energy", bins, 0.0, 5.0);
  TH2D *Q2_vs_xb =
      new TH2D("Q2_vs_xb", "Q^{2} vs x_{b}", bins, 0.1, 0.6, bins, 1.0, 3.5);
  TH2D *WvsQ2_proton = new TH2D("WvsQ2_proton", "W vs Q^{2} P", bins, w_min,
                                w_max, bins, q2_min, q2_max);
  TH1D *W_proton = new TH1D("W_proton", "W P", bins, w_min, w_max);
  TH1D *Q2_proton = new TH1D("Q2_proton", "Q^{2} P", bins, q2_min, q2_max);
  TH2D *WvsQ2_pion = new TH2D("WvsQ2_pion", "W vs Q^{2} #pi^{+} only", bins,
                              w_min, w_max, bins, q2_min, q2_max);
  TH1D *W_pion = new TH1D("W_pion", "W #pi^{+} only", bins, w_min, w_max);
  TH1D *Q2_pion =
      new TH1D("Q2_pion", "Q^{2} #pi^{+} only", bins, q2_min, q2_max);
  TH2D *WvsQ2_single_pi = new TH2D("WvsQ2_single_pi", "W vs Q^{2} #pi^{+}",
                                   bins, w_min, w_max, bins, q2_min, q2_max);
  TH1D *W_single_pi = new TH1D("W_single_pi", "W #pi^{+}", bins, w_min, w_max);
  TH1D *Q2_single_pi =
      new TH1D("Q2_single_pi", "Q^{2} #pi^{+}", bins, q2_min, q2_max);
  TH2D *WvsQ2_single_proton =
      new TH2D("WvsQ2_single_proton", "W vs Q^{2} P", bins, w_min, w_max, bins,
               q2_min, q2_max);
  TH1D *W_single_proton =
      new TH1D("W_single_proton", "W P", bins, w_min, w_max);
  TH1D *Q2_single_proton =
      new TH1D("Q2_single_proton", "Q^{2} P", bins, q2_min, q2_max);
  // W and Q^2

  // P and E
  double b_min = 0.1;
  double b_max = 1.2;
  TH2D *MomVsBeta_hist = new TH2D("MomVsBeta", "Momentum versus #beta", bins,
                                  p_min, p_max, bins, b_min, b_max);
  TH2D *MomVsBeta_hist_pos =
      new TH2D("MomVsBeta_pos", "Momentum versus #beta Positive", bins, p_min,
               p_max, bins, b_min, b_max);
  TH2D *MomVsBeta_hist_neg =
      new TH2D("MomVsBeta_neg", "Momentum versus #beta Negative", bins, p_min,
               p_max, bins, b_min, b_max);
  TH1D *Mom = new TH1D("Momentum", "Momentum", bins, 0, 2.5);
  TH1D *Energy_hist = new TH1D("Energy_hist", "Energy_hist", bins, 0.0, 2.5);
  TH2D *MomVsBeta_proton_ID =
      new TH2D("MomVsBeta_proton_ID", "Momentum versus #beta p", bins, p_min,
               p_max, bins, b_min, b_max);
  TH2D *MomVsBeta_Pi_ID =
      new TH2D("MomVsBeta_Pi_ID", "Momentum versus #beta #pi^{+}", bins, p_min,
               p_max, bins, b_min, b_max);
  TH2D *MomVsBeta_proton_Pi_ID =
      new TH2D("MomVsBeta_proton_Pi_ID", "Momentum versus #beta P #pi^{+}",
               bins, p_min, p_max, bins, b_min, b_max);
  // P and E

  // Missing Mass
  int bins_MM = 200;
  double MM_min = 0.0;
  double MM_max = 3.0;
  TH1D *Mass = new TH1D("Mass", "Mass", 600, 0, 6);
  // Missing Mass

  // Delta T
  // Histogram declarations, fills, and write
  // j -> type: 0=>Proton,1=>Pip,2=>Electron
  // jj -> Fit point
  static const int N_SIGMA = 3;
  double Dt_min = -10;
  double Dt_max = 10;
  static const int num_points = 20;
  TH1D *delta_t_hist[3][num_points];
  const double bin_width = (p_max - p_min) / num_points;

  static const int sc_sector_num = 6;
  static const int sc_paddle_num = 48;
  TH2D *delta_t_sec_pad_hist[3][sc_sector_num][sc_paddle_num];
  TH2D *delta_t_mass_P =
      new TH2D("delta_t_mass_P", "#Deltat assuming mass of proton", bins, p_min,
               p_max, bins, Dt_min, Dt_max);
  TH2D *delta_t_mass_P_PID = new TH2D(
      "delta_t_mass_P_PID", "#Deltat assuming mass of proton with PID proton",
      bins, p_min, p_max, bins, Dt_min, Dt_max);

  TH2D *delta_t_mass_PIP =
      new TH2D("delta_t_mass_PIP", "#Deltat assuming mass of #pi^{+}", bins,
               p_min, p_max, bins, Dt_min, Dt_max);
  TH2D *delta_t_mass_PIP_PID =
      new TH2D("delta_t_mass_PIP_PID",
               "#Deltat assuming mass of #pi^{+} with PID #pi^{+}", bins, p_min,
               p_max, bins, Dt_min, Dt_max);

  TH2D *delta_t_mass_electron =
      new TH2D("delta_t_mass_electron", "#Deltat assuming mass of e^{-}", bins,
               p_min, p_max, bins, Dt_min, Dt_max);
  TH2D *delta_t_mass_electron_PID =
      new TH2D("delta_t_mass_electron_PID",
               "#Deltat assuming mass of e^{-} with PID e^{-}", bins, p_min,
               p_max, bins, Dt_min, Dt_max);

  TH2D *delta_t_mass_positron =
      new TH2D("delta_t_mass_postitron", "#Deltat assuming mass of e^{+}", bins,
               p_min, p_max, bins, Dt_min, Dt_max);
  TH2D *delta_t_mass_positron_PID =
      new TH2D("delta_t_mass_postitron_PID",
               "#Deltat assuming mass of e^{+} with PID e^{+}", bins, p_min,
               p_max, bins, Dt_min, Dt_max);
  // Delta T

  // cc hist
  int bins_CC = 50;
  double CC_min = 0;
  double CC_max = 250;
  char *L_R_C;
  int sector = 6;
  int segment = 18;
  int PMT = 3;
  TH1D *cc_hist[6][18][3];
  TH1D *cc_hist_allSeg[6][3];
  static const Int_t ndims_cc_sparse = 4;
  Int_t bins_cc_sparse[ndims_cc_sparse] = {sector, segment, PMT, bins_CC};
  Double_t xmin_cc_sparse[ndims_cc_sparse] = {0.0, 0.0, -2.0, CC_min};
  Double_t xmax_cc_sparse[ndims_cc_sparse] = {sector + 1.0, segment + 1.0, 1.0,
                                              CC_max};
  Double_t x_cc_sparse[ndims_cc_sparse];
  THnSparse *cc_sparse =
      new THnSparseD("cc_sparse", "Histogram", ndims_cc_sparse, bins_cc_sparse,
                     xmin_cc_sparse, xmax_cc_sparse);
  // cc hist

  // fiducial
  double theta_min = 0;
  double theta_max = 90;
  double phi_min = -360 / 2.0;
  double phi_max = 360 / 2.0;

  static const int sector_num = 6;
  TH2D *fid_sec_hist[sector_num];
  TH2D *fid_hist = new TH2D("fid", "fid", bins, phi_min, phi_max, bins,
                            theta_min, theta_max);
  // fiducial

  // EC hists
  double EC_min = 0;
  double EC_max = 1;
  TH2D *EC_sampling_fraction =
      new TH2D("EC_sampling_fraction", "EC_sampling_fraction", bins, p_min,
               p_max, bins, EC_min, EC_max);
  // EC hists

public:
  Histogram();
  ~Histogram();
  TH1D *Missing_Mass =
      new TH1D("Missing_Mass", "Missing Mass", bins_MM, MM_min, MM_max);
  TH1D *Missing_Mass_square =
      new TH1D("Missing_Mass_square", "Missing Mass square", bins_MM, MM_min,
               MM_max * MM_max);
  // W and Q^2
  void Fill_proton_WQ2(double W, double Q2);
  void Fill_single_pi_WQ2(double W, double Q2);
  void Fill_single_proton_WQ2(double W, double Q2);
  void WvsQ2_Fill(double E_prime, double W, double Q2, double xb);
  void Fill_pion_WQ2(double W, double Q2);
  void WvsQ2_Write();

  // P and E
  void MomVsBeta_Fill_pos(double P, double Beta);
  void MomVsBeta_Fill_neg(double P, double Beta);
  void Fill_proton_ID_P(double p, double beta);
  void Fill_Pi_ID_P(double p, double beta);
  void Fill_proton_Pi_ID_P(double p, double beta);
  void MomVsBeta_Fill(double Energy, double P, double Beta);
  void MomVsBeta_Write();

  // Missing Mass
  void Fill_Missing_Mass(double miss_mass);
  void Fill_Mass(double mass);
  void Fill_Missing_Mass_square(double miss_mass_2);
  void Write_Missing_Mass();

  // Delta T
  void Fill_deltat_P(double momentum, double delta_t);
  void Fill_deltat_P_PID(double momentum, double delta_t);
  void Fill_deltat_PIP(double momentum, double delta_t);
  void Fill_deltat_PIP_PID(double momentum, double delta_t);
  void Fill_deltat_electron(double momentum, double delta_t);
  void Fill_deltat_electron_PID(double momentum, double delta_t);
  void Fill_deltat_positron(double momentum, double delta_t);
  void Fill_deltat_positron_PID(double momentum, double delta_t);
  void delta_t_slice_fit();
  void delta_t_Write();
  void delta_t_Fill(double momentum, int charge, double delta_t_proton,
                    double delta_t_pip, double delta_t_electron);
  void delta_t_slices_Write();
  void delta_t_sec_pad(double momentum, int charge, double delta_t_proton,
                       double delta_t_pip, double delta_t_electron,
                       int sc_sector, int sc_paddle);
  void delta_t_sec_pad_Write();
  void delta_T_canvas();

  // cc hist
  void CC_fill(int cc_sector, int cc_segment, int cc_pmt, int cc_nphe);
  void CC_Write();
  void CC_canvas();

  // fiducial hist
  void Fill_fid(double theta, double phi, int sector);
  void Fid_Write();

  // EC hists
  void EC_fill(double etot, double momentum);
  void EC_Write();
};

#endif
