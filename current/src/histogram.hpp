/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef HISTOGRAM_H
#define HISTOGRAM_H
#include <cmath>
#include <fstream>
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "color.hpp"
#include "fits.hpp"
#include "missing_mass.hpp"

#define W_bins 20
#define Q2_bins 10
#define theta_bins 100
#define phi_bins 100
#define w_binned_min 1.0
#define w_binned_max 3.0
#define q2_binned_min 1.0
#define q2_binned_max 3.0

class Histogram {
 private:
 public:
  Histogram();
  Histogram(std::string output_file);
  ~Histogram();
  void Write(std::string output_file);

  TFile *RootOutputFile;
  TCanvas *def;
  void makeHists_fid();
  void makeHists_deltat();
  void makeHists_CC();
  void makeHists_WvsQ2();
  int bins = 500;
  double p_min = 0.0;
  double p_max = 5.0;
  char hname[50];
  char htitle[500];

  static const int sector = 6;

  // Missing Mass
  int bins_MM = 300;
  double MM_min = 0.0;
  double MM_max = 3.0;
  // Missing Mass

  // W and Q^2
  double w_min = 0;
  double w_max = 3.25;
  double q2_min = 0;
  double q2_max = 5;

  double W_width = (w_binned_max - w_binned_min) / (double)W_bins;
  double Q2_width = (q2_binned_max - q2_binned_min) / (double)Q2_bins;

  std::unique_ptr<TH2D> WvsQ2_hist =
      std::make_unique<TH2D>("WvsQ2_hist", "W vs Q^{2}", bins, w_min, w_max, bins, q2_min, q2_max);
  std::unique_ptr<TH1D> W_hist = std::make_unique<TH1D>("W", "W", bins, w_min, w_max);
  std::unique_ptr<TH1D> Q2_hist = std::make_unique<TH1D>("Q2", "Q^{2}", bins, q2_min, q2_max);
  std::unique_ptr<TH1D> E_prime_hist = std::make_unique<TH1D>("E_prime", "Scattered Electron Energy", bins, 0.0, 5.0);
  std::unique_ptr<TH1D> photon_flux_hist = std::make_unique<TH1D>("photon_flux", "Photon Flux", bins, -0.1, 0.1);
  std::unique_ptr<TH2D> Q2_vs_xb = std::make_unique<TH2D>("Q2_vs_xb", "Q^{2} vs x_{b}", bins, 0.1, 0.6, bins, 1.0, 3.5);
  std::unique_ptr<TH2D> WvsQ2_proton =
      std::make_unique<TH2D>("WvsQ2_proton", "W vs Q^{2} P", bins, w_min, w_max, bins, q2_min, q2_max);
  std::unique_ptr<TH1D> W_proton = std::make_unique<TH1D>("W_proton", "W P", bins, w_min, w_max);
  std::unique_ptr<TH1D> Q2_proton = std::make_unique<TH1D>("Q2_proton", "Q^{2} P", bins, q2_min, q2_max);
  std::unique_ptr<TH2D> WvsQ2_pion =
      std::make_unique<TH2D>("WvsQ2_pion", "W vs Q^{2} #pi^{+} only", bins, w_min, w_max, bins, q2_min, q2_max);
  std::unique_ptr<TH1D> W_pion = std::make_unique<TH1D>("W_pion", "W #pi^{+} only", bins, w_min, w_max);
  std::unique_ptr<TH1D> Q2_pion = std::make_unique<TH1D>("Q2_pion", "Q^{2} #pi^{+} only", bins, q2_min, q2_max);
  std::unique_ptr<TH2D> WvsQ2_single_pi =
      std::make_unique<TH2D>("WvsQ2_single_pi", "W vs Q^{2} #pi^{+}", bins, w_min, w_max, bins, q2_min, q2_max);
  std::unique_ptr<TH1D> W_single_pi = std::make_unique<TH1D>("W_single_pi", "W #pi^{+}", bins, w_min, w_max);
  std::unique_ptr<TH1D> Q2_single_pi = std::make_unique<TH1D>("Q2_single_pi", "Q^{2} #pi^{+}", bins, q2_min, q2_max);
  std::unique_ptr<TH2D> WvsQ2_MM =
      std::make_unique<TH2D>("WvsQ2_MM", "W vs Q^{2} mm N cut", bins, w_min, w_max, bins, q2_min, q2_max);
  std::unique_ptr<TH1D> W_MM = std::make_unique<TH1D>("W_MM", "W mm N cut", bins, w_min, w_max);
  std::unique_ptr<TH1D> Q2_MM = std::make_unique<TH1D>("Q2_MM", "Q^{2} mm N cut", bins, q2_min, q2_max);
  std::unique_ptr<TH2D> WvsQ2_channel =
      std::make_unique<TH2D>("WvsQ2_channel", "W vs Q^{2} #pi^{+} N", bins, w_min, w_max, bins, q2_min, q2_max);
  std::unique_ptr<TH1D> W_channel = std::make_unique<TH1D>("W_channel", "W #pi^{+} N", bins, w_min, w_max);
  std::unique_ptr<TH1D> Q2_channel = std::make_unique<TH1D>("Q2_channel", "Q^{2} #pi^{+} N", bins, q2_min, q2_max);
  std::unique_ptr<TH2D> WvsQ2_Ppi0 =
      std::make_unique<TH2D>("WvsQ2_Ppi0", "W vs Q^{2} P #pi^{0}", bins, w_min, w_max, bins, q2_min, q2_max);
  std::unique_ptr<TH1D> W_Ppi0 = std::make_unique<TH1D>("W_Ppi0", "W P #pi^{0}", bins, w_min, w_max);
  std::unique_ptr<TH1D> Q2_Ppi0 = std::make_unique<TH1D>("Q2_Ppi0", "Q^{2} P #pi^{0}", bins, q2_min, q2_max);
  std::unique_ptr<TH2D> WvsQ2_single_proton =
      std::make_unique<TH2D>("WvsQ2_single_proton", "W vs Q^{2} P", bins, w_min, w_max, bins, q2_min, q2_max);
  std::unique_ptr<TH1D> W_single_proton = std::make_unique<TH1D>("W_single_proton", "W P", bins, w_min, w_max);
  std::unique_ptr<TH1D> Q2_single_proton = std::make_unique<TH1D>("Q2_single_proton", "Q^{2} P", bins, q2_min, q2_max);
  std::unique_ptr<TH2D> WvsQ2_binned =
      std::make_unique<TH2D>("WvsQ2_hist_binned", "W vs Q^{2} binned", W_bins, w_binned_min, w_binned_max, Q2_bins,
                             q2_binned_min, q2_binned_max);

  std::unique_ptr<TH1D> W_binned[Q2_bins];
  std::unique_ptr<TH1D> Q2_binned[W_bins];

  std::unique_ptr<TH1D> Missing_Mass_WBinned[W_bins];
  Fits *Fit_Missing_Mass_WBinned[W_bins];
  std::unique_ptr<TH1D> Missing_Mass_WBinned_square[W_bins];

  static const int ndims_pip_N = 5;
  int bins_pip_N[ndims_pip_N] = {W_bins, Q2_bins, sector, bins, bins};                  // theta_bins, phi_bins, };
  double xmin_pip_N[ndims_pip_N] = {w_binned_min, q2_binned_min, 0.0, MM_min, MM_min};  // 0.0, -TMath::Pi(), };
  double xmax_pip_N[ndims_pip_N] = {w_binned_max, q2_binned_max, sector, MM_max,
                                    MM_max *MM_max};  // TMath::Pi() / 2.0,TMath::Pi(),  };
  double x_pip_N[ndims_pip_N];
  THnD *pip_N = new THnD("pip_N", "WvsQ2_NPIP", ndims_pip_N, bins_pip_N, xmin_pip_N, xmax_pip_N);

  // W and Q^2

  // P and E
  double b_min = 0.1;
  double b_max = 1.2;
  std::unique_ptr<TH2D> MomVsBeta_hist =
      std::make_unique<TH2D>("MomVsBeta", "Momentum versus #beta", bins, p_min, p_max, bins, b_min, b_max);
  std::unique_ptr<TH2D> MomVsBeta_hist_pos =
      std::make_unique<TH2D>("MomVsBeta_pos", "Momentum versus #beta Positive", bins, p_min, p_max, bins, b_min, b_max);
  std::unique_ptr<TH2D> MomVsBeta_hist_neg =
      std::make_unique<TH2D>("MomVsBeta_neg", "Momentum versus #beta Negative", bins, p_min, p_max, bins, b_min, b_max);
  std::unique_ptr<TH1D> Mom = std::make_unique<TH1D>("Momentum", "Momentum", bins, 0, 2.5);
  std::unique_ptr<TH1D> Energy_hist = std::make_unique<TH1D>("Energy_hist", "Energy_hist", bins, 0.0, 2.5);
  std::unique_ptr<TH2D> MomVsBeta_proton_ID =
      std::make_unique<TH2D>("MomVsBeta_proton_ID", "Momentum versus #beta p", bins, p_min, p_max, bins, b_min, b_max);
  std::unique_ptr<TH2D> MomVsBeta_Pi_ID = std::make_unique<TH2D>("MomVsBeta_Pi_ID", "Momentum versus #beta #pi^{+}",
                                                                 bins, p_min, p_max, bins, b_min, b_max);
  std::unique_ptr<TH2D> MomVsBeta_proton_Pi_ID = std::make_unique<TH2D>(
      "MomVsBeta_proton_Pi_ID", "Momentum versus #beta P #pi^{+}", bins, p_min, p_max, bins, b_min, b_max);
  // P and E

  // Delta T
  // Histogram declarations, fills, and write
  // j -> type: 0=>Proton,1=>Pip,2=>Electron
  // jj -> Fit point
  static const int N_SIGMA = 3;
  double Dt_min = -10;
  double Dt_max = 10;
  static const int num_points = 20;
  std::unique_ptr<TH1D> delta_t_hist[3][num_points];
  double bin_width = (p_max - p_min) / num_points;

  static constexpr int sc_paddle_num = 48;
  std::unique_ptr<TH2D> delta_t_sec_pad_hist[3][sector][sc_paddle_num];
  std::unique_ptr<TH2D> delta_t_mass_P = std::make_unique<TH2D>("delta_t_mass_P", "#Deltat assuming mass of proton",
                                                                bins, p_min, p_max, bins, Dt_min, Dt_max);
  std::unique_ptr<TH2D> delta_t_mass_P_PID =
      std::make_unique<TH2D>("delta_t_mass_P_PID", "#Deltat assuming mass of proton with PID proton", bins, p_min,
                             p_max, bins, Dt_min, Dt_max);

  std::unique_ptr<TH2D> delta_t_mass_PIP = std::make_unique<TH2D>(
      "delta_t_mass_PIP", "#Deltat assuming mass of #pi^{+}", bins, p_min, p_max, bins, Dt_min, Dt_max);
  std::unique_ptr<TH2D> delta_t_mass_PIP_PID =
      std::make_unique<TH2D>("delta_t_mass_PIP_PID", "#Deltat assuming mass of #pi^{+} with PID #pi^{+}", bins, p_min,
                             p_max, bins, Dt_min, Dt_max);

  std::unique_ptr<TH2D> delta_t_mass_PIM = std::make_unique<TH2D>(
      "delta_t_mass_PIM", "#Deltat assuming mass of #pi^{-}", bins, p_min, p_max, bins, Dt_min, Dt_max);
  std::unique_ptr<TH2D> delta_t_mass_PIM_PID =
      std::make_unique<TH2D>("delta_t_mass_PIM_PID", "#Deltat assuming mass of #pi^{-} with PID #pi^{-}", bins, p_min,
                             p_max, bins, Dt_min, Dt_max);

  std::unique_ptr<TH2D> delta_t_mass_electron = std::make_unique<TH2D>(
      "delta_t_mass_electron", "#Deltat assuming mass of e^{-}", bins, p_min, p_max, bins, Dt_min, Dt_max);
  std::unique_ptr<TH2D> delta_t_mass_electron_PID =
      std::make_unique<TH2D>("delta_t_mass_electron_PID", "#Deltat assuming mass of e^{-} with PID e^{-}", bins, p_min,
                             p_max, bins, Dt_min, Dt_max);

  std::unique_ptr<TH2D> delta_t_mass_kp = std::make_unique<TH2D>("delta_t_mass_kp", "#Deltat assuming mass of k^{+}",
                                                                 bins, p_min, p_max, bins, Dt_min, Dt_max);
  std::unique_ptr<TH2D> delta_t_mass_kp_PID = std::make_unique<TH2D>(
      "delta_t_mass_kp_PID", "#Deltat assuming mass of k^{+} with PID k^{+}", bins, p_min, p_max, bins, Dt_min, Dt_max);
  // Delta T

  // cc hist
  int bins_CC = 50;
  double CC_min = 0;
  double CC_max = 250;
  std::string L_R_C;
  static const int segment = 18;
  static const int PMT = 3;

  std::unique_ptr<TH1D> cc_hist[sector][segment][PMT];
  std::unique_ptr<TH1D> cc_hist_allSeg[sector][PMT];
  static const int ndims_cc_sparse = 4;
  int bins_cc_sparse[ndims_cc_sparse] = {sector, segment, PMT, bins_CC};
  double xmin_cc_sparse[ndims_cc_sparse] = {0.0, 0.0, -2.0, CC_min};
  double xmax_cc_sparse[ndims_cc_sparse] = {sector + 1.0, segment + 1.0, 1.0, CC_max};
  double x_cc_sparse[ndims_cc_sparse];
  THnD *cc_sparse = new THnD("cc_sparse", "Histogram", ndims_cc_sparse, bins_cc_sparse, xmin_cc_sparse, xmax_cc_sparse);

  std::unique_ptr<TH2D> Theta_CC = std::make_unique<TH2D>("Theta_CC", "Theta_CC", 20, 0.0, 20.0, 60, 0.0, 60.0);
  std::unique_ptr<TH2D> Theta_CC_Sec[sector];
  std::unique_ptr<TH2D> Theta_CC_Sec_cut[sector];
  // cc hist

  // fiducial
  double theta_min = 0;
  double theta_max = 60;
  double phi_min = -360 / 2.0;
  double phi_max = 360 / 2.0;

  static const int fid_slices = 10;
  static const int start_slice = 0;
  std::vector<std::unique_ptr<TH2D> > electron_fid_sec_hist;
  TH1D *electron_fid_sec_slice[sector][fid_slices];
  std::unique_ptr<TH2D> electron_fid_hist =
      std::make_unique<TH2D>("electron_fid", "electron_fid", bins, phi_min, phi_max, bins, theta_min, theta_max);

  std::unique_ptr<TH2D> hadron_fid_sec_hist[3][sector];
  std::unique_ptr<TH1D> hadron_fid_sec_slice[sector][fid_slices];
  std::unique_ptr<TH2D> hadron_fid_hist[3];
  // fiducial

  // EC hists
  double EC_min = 0;
  double EC_max = 1;
  std::unique_ptr<TH2D> EC_sampling_fraction =
      std::make_unique<TH2D>("EC_sampling_fraction", "EC_sampling_fraction", bins, p_min, p_max, bins, EC_min, EC_max);
  std::unique_ptr<TH2D> EC_sampling_fraction_cut = std::make_unique<TH2D>(
      "EC_sampling_fraction_cut", "EC_sampling_fraction_cut", bins, p_min, p_max, bins, EC_min, EC_max);
  std::unique_ptr<TH1D> EC_hist[num_points];
  std::unique_ptr<TH1D> EC_hist_cut[num_points];
  std::unique_ptr<TH2D> Theta_vs_mom =
      std::make_unique<TH2D>("Theta_vs_mom", "Theta_vs_mom", bins, p_min, p_max, bins, 0, 100);
  std::unique_ptr<TH2D> ECin_ECout = std::make_unique<TH2D>("ECin_ECout", "ECin_ECout", bins, 0.0, 0.5, bins, 0.0, 0.5);
  // EC hists

  // Beam Position
  double x_y_min_max = 0.5;
  std::unique_ptr<TH1D> Beam_Position_X =
      std::make_unique<TH1D>("Beam_Position_X", "Beam_Position_X", bins, -x_y_min_max, x_y_min_max);
  std::unique_ptr<TH1D> Beam_Position_Y =
      std::make_unique<TH1D>("Beam_Position_Y", "Beam_Position_Y", bins, -x_y_min_max, x_y_min_max);
  std::unique_ptr<TH1D> Beam_Position_Z = std::make_unique<TH1D>("Beam_Position_Z", "Beam_Position_Z", bins, -5.0, 5.0);

  std::unique_ptr<TH2D> Beam_Position = std::make_unique<TH2D>("Beam_Position", "Beam_Position", bins, -x_y_min_max,
                                                               x_y_min_max, bins, -x_y_min_max, x_y_min_max);

  // Beam Position

  // Vertex
  std::unique_ptr<TH1D> target_vertex_X = std::make_unique<TH1D>("Target_vertex_X", "Target_vertex_X", bins, -6.0, 6.0);
  std::unique_ptr<TH1D> target_vertex_Y = std::make_unique<TH1D>("Target_vertex_Y", "Target_vertex_Y", bins, -6.0, 6.0);
  std::unique_ptr<TH1D> target_vertex_Z = std::make_unique<TH1D>("Target_vertex_Z", "Target_vertex_Z", bins, -6.0, 6.0);

  std::unique_ptr<TH2D> target_vertex_xy =
      std::make_unique<TH2D>("Target_vertex_xy", "Target_vertex_xy", bins, -6.0, 6.0, bins, -6.0, 6.0);
  std::unique_ptr<TH2D> target_vertex_zy =
      std::make_unique<TH2D>("Target_vertex_zy", "Target_vertex_zy", bins, -6.0, 6.0, bins, -6.0, 6.0);
  std::unique_ptr<TH2D> target_vertex_zx =
      std::make_unique<TH2D>("Target_vertex_zx", "Target_vertex_zx", bins, -6.0, 6.0, bins, -6.0, 6.0);

  std::unique_ptr<TH1D> Missing_Mass = std::make_unique<TH1D>("Missing_Mass", "Missing Mass", bins_MM, MM_min, MM_max);
  std::unique_ptr<TH1D> Missing_Mass_square =
      std::make_unique<TH1D>("Missing_Mass_square", "Missing Mass square", bins_MM, MM_min, MM_max *MM_max);

  std::unique_ptr<TH1D> Missing_Mass_strict =
      std::make_unique<TH1D>("Missing_Mass_strict", "Missing Mass", bins_MM, MM_min, MM_max);
  std::unique_ptr<TH1D> Missing_Mass_square_strict =
      std::make_unique<TH1D>("Missing_Mass_square_strict", "Missing Mass square", bins_MM, MM_min, MM_max *MM_max);

  std::unique_ptr<TH1D> Missing_Mass_pi0 =
      std::make_unique<TH1D>("Missing_Mass_pi0", "Missing Mass #pi^{0}", bins_MM, -3, 3);
  std::unique_ptr<TH1D> Missing_Mass_square_pi0 =
      std::make_unique<TH1D>("Missing_Mass_pi0_2", "MM^{2} #pi^{0}", bins_MM, -3, 3);

  std::unique_ptr<TH1D> Missing_Mass_2pi =
      std::make_unique<TH1D>("Missing_Mass_2pi", "Missing Mass 2 #pi", bins_MM, MM_min, MM_max);
  std::unique_ptr<TH1D> Missing_Mass_square_2pi =
      std::make_unique<TH1D>("Missing_Mass_square_2pi", "Missing Mass 2 #pi", bins_MM, MM_min, MM_max *MM_max);

  std::unique_ptr<TH1D> Missing_Mass_nutron_no2pi =
      std::make_unique<TH1D>("Missing_Mass_nutron_no2pi", "Missing Mass removing 2 #pi", bins_MM, MM_min, MM_max);

  std::unique_ptr<TH1D> energy_no_cuts =
      std::make_unique<TH1D>("Energy_no_cuts", "Scattered electron energy", 500, 0.0, 5.0);
  std::unique_ptr<TH1D> energy_fid_cuts =
      std::make_unique<TH1D>("Energy_fid_cuts", "Scattered electron energy after fiducial cuts", 500, 0.0, 5.0);
  std::unique_ptr<TH1D> energy_channel_cuts =
      std::make_unique<TH1D>("Energy_channel_cuts", "Scattered electron energy for N #pi^{+} events", 500, 0.0, 5.0);

  // W and Q^2
  void Fill_proton_WQ2(double W, double Q2);
  void Fill_P_PI0(double W, double Q2);
  void Fill_single_pi_WQ2(double W, double Q2);
  void Fill_MM_WQ2(double W, double Q2);
  void Fill_channel_WQ2(double W, double Q2, TLorentzVector e_prime, double mm, double mm2, int sec);
  void Fill_single_proton_WQ2(double W, double Q2);
  void WvsQ2_Fill(double W, double Q2);
  void Fill_pion_WQ2(double W, double Q2);
  void WvsQ2_Write();
  void WvsQ2_binned_Write();

  // P and E
  void MomVsBeta_Fill_pos(double P, double Beta);
  void MomVsBeta_Fill_neg(double P, double Beta);
  void Fill_proton_ID_P(double p, double beta);
  void Fill_Pi_ID_P(double p, double beta);
  void Fill_proton_Pi_ID_P(double p, double beta);
  void MomVsBeta_Fill(double P, double Beta);
  void Photon_flux_Fill(double photon_flux);
  void MomVsBeta_Write();

  // Missing Mass
  void Fill_Missing_Mass(double miss_mass);
  void Fill_Missing_Mass(float mm, float mm2);
  void Fill_Missing_Mass_strict(float mm, float mm2);
  // void Fill_Mass(double mass);
  void Fill_Missing_Mass_square(double miss_mass_2);
  void Write_Missing_Mass();

  void Fill_Missing_Mass_pi0(float mm, float mm2);
  void Fill_Missing_Mass_twoPi(float mm, float mm2);
  void Fill_W_Missing_Mass(double W, float mm, float mm2);

  // Delta T
  void Fill_deltat_P(double momentum, double delta_t);
  void Fill_deltat_P_PID(double momentum, double delta_t);
  void Fill_deltat_PIP(double momentum, double delta_t);
  void Fill_deltat_PIP_PID(double momentum, double delta_t);
  void Fill_deltat_PIM(double momentum, double delta_t);
  void Fill_deltat_PIM_PID(double momentum, double delta_t);
  void Fill_deltat_electron(double momentum, double delta_t);
  void Fill_deltat_electron_PID(double momentum, double delta_t);
  void Fill_deltat_kp(double momentum, double delta_t);
  void Fill_deltat_kp_PID(double momentum, double delta_t);
  void delta_t_slice_fit();
  void delta_t_Write();
  void delta_t_Fill(double momentum, int charge, double delta_t_proton, double delta_t_pip, double delta_t_electron);
  void delta_t_slices_Write();
  void delta_t_sec_pad(double momentum, int charge, double delta_t_proton, double delta_t_pip, double delta_t_electron,
                       int sc_sector, int sc_paddle);
  void delta_t_sec_pad_Write();
  void delta_T_canvas();

  // cc hist
  void CC_fill(int cc_sector, int cc_segment, int cc_pmt, int cc_nphe, double theta_cc);
  void CC_Write();
  void Theta_CC_Write();
  void theta_cc_slice_fit();
  void CC_canvas();

  // fiducial hist
  void Fill_electron_fid(double theta, double phi, int sector);
  void Fill_hadron_fid(double theta, double phi, int sector, int id);
  void Fid_Write();
  void fid_canvas();

  // EC hists
  void makeHists_EC();
  void EC_slices_Write();
  void EC_fill(double etot, double momentum);
  void EC_inout(double Ein, double Eout);
  void EC_cut_fill(double etot, double momentum);
  void EC_slice_fit();
  void EC_Write();

  void TM_Fill(double momentum, double theta);

  // Beam Position
  void Fill_Beam_Position(double vertex_x, double vertex_y, double vertex_z);
  void Beam_Position_Write();

  void Fill_Target_Vertex(double vertex_x, double vertex_y, double vertex_z);
  void Target_Vertex_Write();

  void Fill_E_Prime(TLorentzVector e_prime);
  void Fill_E_Prime_fid(TLorentzVector e_prime);
  void Fill_E_Prime_channel(TLorentzVector e_prime);
  void E_Prime_Write();
};

#endif
