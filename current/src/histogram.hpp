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

#define W_BINS 20
#define Q2_BINS 10
#define THETA_BINS 100
#define PHI_BINS 100
#define NUM_SECTORS 6
#define NDIMS_PIP_N 7
#define SC_PADDLE_NUM 48

#define BINS 500
#define BINS_MM 300

class Histogram {
 private:
 public:
  Histogram();
  Histogram(std::string output_file);
  ~Histogram();
  void Write(std::string output_file);

  TFile* RootOutputFile;
  TCanvas* def;
  void makeHists_fid();
  void makeHists_deltat();
  void makeHists_CC();
  void makeHists_WvsQ2();

  float p_min = 0.0;
  float p_max = 5.0;
  char hname[50];
  char htitle[500];

  float w_binned_min = 1.0;
  float w_binned_max = 3.0;
  float q2_binned_min = 1.0;
  float q2_binned_max = 3.0;

  // Missing Mass
  float MM_min = 0.0;
  float MM_max = 3.0;
  // Missing Mass

  // W and Q^2
  float w_min = 0;
  float w_max = 3.25;
  float q2_min = 0;
  float q2_max = 5;

  float W_width = (w_binned_max - w_binned_min) / (float)W_BINS;
  float Q2_width = (q2_binned_max - q2_binned_min) / (float)Q2_BINS;

  TH2D* WvsQ2_hist = new TH2D("WvsQ2_hist", "W vs Q^{2}", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1D* W_hist = new TH1D("W", "W", BINS, w_min, w_max);
  TH1D* Q2_hist = new TH1D("Q2", "Q^{2}", BINS, q2_min, q2_max);
  TH1D* E_prime_hist = new TH1D("E_prime", "Scattered Electron Energy", BINS, 0.0, 5.0);
  TH1D* photon_flux_hist = new TH1D("photon_flux", "Photon Flux", BINS, -0.1, 0.1);
  TH2D* Q2_vs_xb = new TH2D("Q2_vs_xb", "Q^{2} vs x_{b}", BINS, 0.1, 0.6, BINS, 1.0, 3.5);
  TH2D* WvsQ2_proton = new TH2D("WvsQ2_proton", "W vs Q^{2} P", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1D* W_proton = new TH1D("W_proton", "W P", BINS, w_min, w_max);
  TH1D* Q2_proton = new TH1D("Q2_proton", "Q^{2} P", BINS, q2_min, q2_max);
  TH2D* WvsQ2_pion = new TH2D("WvsQ2_pion", "W vs Q^{2} #pi^{+} only", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1D* W_pion = new TH1D("W_pion", "W #pi^{+} only", BINS, w_min, w_max);
  TH1D* Q2_pion = new TH1D("Q2_pion", "Q^{2} #pi^{+} only", BINS, q2_min, q2_max);
  TH2D* WvsQ2_NeutronPip =
      new TH2D("WvsQ2_NeutronPip", "W vs Q^{2} N #pi^{+}", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1D* W_NeutronPip = new TH1D("W_NeutronPip", "W N #pi^{+}", BINS, w_min, w_max);
  TH1D* Q2_NeutronPip = new TH1D("Q2_NeutronPip", "Q^{2} N #pi^{+}", BINS, q2_min, q2_max);
  TH2D* WvsQ2_MM = new TH2D("WvsQ2_MM", "W vs Q^{2} mm N cut", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1D* W_MM = new TH1D("W_MM", "W mm N cut", BINS, w_min, w_max);
  TH1D* Q2_MM = new TH1D("Q2_MM", "Q^{2} mm N cut", BINS, q2_min, q2_max);
  TH2D* WvsQ2_channel = new TH2D("WvsQ2_channel", "W vs Q^{2} #pi^{+} N", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1D* W_channel = new TH1D("W_channel", "W #pi^{+} N", BINS, w_min, w_max);
  TH1D* Q2_channel = new TH1D("Q2_channel", "Q^{2} #pi^{+} N", BINS, q2_min, q2_max);
  TH2D* WvsQ2_Ppi0 = new TH2D("WvsQ2_Ppi0", "W vs Q^{2} P #pi^{0}", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1D* W_Ppi0 = new TH1D("W_Ppi0", "W P #pi^{0}", BINS, w_min, w_max);
  TH1D* Q2_Ppi0 = new TH1D("Q2_Ppi0", "Q^{2} P #pi^{0}", BINS, q2_min, q2_max);
  TH2D* WvsQ2_single_proton = new TH2D("WvsQ2_single_proton", "W vs Q^{2} P", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1D* W_single_proton = new TH1D("W_single_proton", "W P", BINS, w_min, w_max);
  TH1D* Q2_single_proton = new TH1D("Q2_single_proton", "Q^{2} P", BINS, q2_min, q2_max);
  TH2D* WvsQ2_binned = new TH2D("WvsQ2_hist_binned", "W vs Q^{2} binned", W_BINS, w_binned_min, w_binned_max, Q2_BINS,
                                q2_binned_min, q2_binned_max);

  TH1D* W_binned[Q2_BINS];
  TH1D* Q2_binned[W_BINS];

  TH1D* Missing_Mass_WBinned[W_BINS];
  Fits* Fit_Missing_Mass_WBinned[W_BINS];
  TH1D* Missing_Mass_WBinned_square[W_BINS];

  /*
    int bins_pip_N[NDIMS_PIP_N] = {W_BINS, Q2_BINS, NUM_SECTORS, BINS_MM, BINS_MM, THETA_BINS, PHI_BINS};
    float xmin_pip_N[NDIMS_PIP_N] = {w_binned_min,     q2_binned_min,   0.0, MM_min, MM_min,
                                     -2 * TMath::Pi(), -2 * TMath::Pi()};
    float xmax_pip_N[NDIMS_PIP_N] = {w_binned_max, q2_binned_max,   NUM_SECTORS,    MM_max,
                                     MM_max,       2 * TMath::Pi(), 2 * TMath::Pi()};
    THnF* pip_N = new THnF("pip_N", "WvsQ2_NPIP", NDIMS_PIP_N, bins_pip_N, xmin_pip_N, xmax_pip_N);
  */

  // W and Q^2

  // P and E
  float b_min = 0.1;
  float b_max = 1.2;
  TH2D* MomVsBeta_hist = new TH2D("MomVsBeta", "Momentum versus #beta", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2D* MomVsBeta_hist_pos =
      new TH2D("MomVsBeta_pos", "Momentum versus #beta Positive", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2D* MomVsBeta_hist_neg =
      new TH2D("MomVsBeta_neg", "Momentum versus #beta Negative", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2D* MomVsBeta_hist_neutral =
      new TH2D("MomVsBeta_Fill_neutral", "Momentum versus #beta neutral", BINS, p_min, p_max, BINS, b_min, b_max);
  TH1D* Mom = new TH1D("Momentum", "Momentum", BINS, 0, 2.5);
  TH1D* Energy_hist = new TH1D("Energy_hist", "Energy_hist", BINS, 0.0, 2.5);
  TH2D* MomVsBeta_proton_ID =
      new TH2D("MomVsBeta_proton_ID", "Momentum versus #beta p", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2D* MomVsBeta_Pi_ID =
      new TH2D("MomVsBeta_Pi_ID", "Momentum versus #beta #pi^{+}", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2D* MomVsBeta_proton_Pi_ID =
      new TH2D("MomVsBeta_proton_Pi_ID", "Momentum versus #beta P #pi^{+}", BINS, p_min, p_max, BINS, b_min, b_max);
  // P and E

  // Delta T
  // Histogram declarations, fills, and write
  // j -> type: 0=>Proton,1=>Pip,2=>Electron
  // jj -> Fit point
  static const int N_SIGMA = 3;
  float Dt_min = -10;
  float Dt_max = 10;
  static const int num_points = 20;
  TH1D* delta_t_hist[3][num_points];
  float bin_width = (p_max - p_min) / num_points;

  TH2D* delta_t_sec_pad_hist[3][NUM_SECTORS][SC_PADDLE_NUM];
  TH2D* delta_t_mass_P =
      new TH2D("delta_t_mass_P", "#Deltat assuming mass of proton", BINS, p_min, p_max, BINS, Dt_min, Dt_max);
  TH2D* delta_t_mass_P_PID = new TH2D("delta_t_mass_P_PID", "#Deltat assuming mass of proton with PID proton", BINS,
                                      p_min, p_max, BINS, Dt_min, Dt_max);

  TH2D* delta_t_mass_PIP =
      new TH2D("delta_t_mass_PIP", "#Deltat assuming mass of #pi^{+}", BINS, p_min, p_max, BINS, Dt_min, Dt_max);
  TH2D* delta_t_mass_PIP_PID = new TH2D("delta_t_mass_PIP_PID", "#Deltat assuming mass of #pi^{+} with PID #pi^{+}",
                                        BINS, p_min, p_max, BINS, Dt_min, Dt_max);

  TH2D* delta_t_mass_PIM =
      new TH2D("delta_t_mass_PIM", "#Deltat assuming mass of #pi^{-}", BINS, p_min, p_max, BINS, Dt_min, Dt_max);
  TH2D* delta_t_mass_PIM_PID = new TH2D("delta_t_mass_PIM_PID", "#Deltat assuming mass of #pi^{-} with PID #pi^{-}",
                                        BINS, p_min, p_max, BINS, Dt_min, Dt_max);

  TH2D* delta_t_mass_electron =
      new TH2D("delta_t_mass_electron", "#Deltat assuming mass of e^{-}", BINS, p_min, p_max, BINS, Dt_min, Dt_max);
  TH2D* delta_t_mass_electron_PID =
      new TH2D("delta_t_mass_electron_PID", "#Deltat assuming mass of e^{-} with PID e^{-}", BINS, p_min, p_max, BINS,
               Dt_min, Dt_max);

  TH2D* delta_t_mass_kp =
      new TH2D("delta_t_mass_kp", "#Deltat assuming mass of k^{+}", BINS, p_min, p_max, BINS, Dt_min, Dt_max);
  TH2D* delta_t_mass_kp_PID = new TH2D("delta_t_mass_kp_PID", "#Deltat assuming mass of k^{+} with PID k^{+}", BINS,
                                       p_min, p_max, BINS, Dt_min, Dt_max);
  // Delta T

  // cc hist
  int bins_CC = 50;
  float CC_min = 0;
  float CC_max = 250;
  std::string L_R_C;
  static const int segment = 18;
  static const int PMT = 3;

  TH1D* cc_hist[NUM_SECTORS][segment][PMT];
  TH1D* cc_hist_allSeg[NUM_SECTORS][PMT];
  /*
  static const int ndims_cc_sparse = 4;
  int bins_cc_sparse[ndims_cc_sparse] = {NUM_SECTORS, segment, PMT, bins_CC};
  float xmin_cc_sparse[ndims_cc_sparse] = {0.0, 0.0, -2.0, CC_min};
  float xmax_cc_sparse[ndims_cc_sparse] = {NUM_SECTORS + 1.0, segment + 1.0, 1.0, CC_max};
  float x_cc_sparse[ndims_cc_sparse];
  THnF* cc_sparse = new THnF("cc_sparse", "Histogram", ndims_cc_sparse, bins_cc_sparse, xmin_cc_sparse, xmax_cc_sparse);
  */

  TH2D* Theta_CC = new TH2D("Theta_CC", "Theta_CC", 20, 0.0, 20.0, 60, 0.0, 60.0);
  TH2D* Theta_CC_Sec[NUM_SECTORS];
  TH2D* Theta_CC_Sec_cut[NUM_SECTORS];
  // cc hist

  // fiducial
  float theta_min = 0;
  float theta_max = 60;
  float phi_min = -360 / 2.0;
  float phi_max = 360 / 2.0;

  static const int fid_slices = 10;
  static const int start_slice = 0;
  std::vector<TH2D*> electron_fid_sec_hist;
  TH1D* electron_fid_sec_slice[NUM_SECTORS][fid_slices];
  TH2D* electron_fid_hist =
      new TH2D("electron_fid", "electron_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);

  TH2D* hadron_fid_sec_hist[3][NUM_SECTORS];
  TH1D* hadron_fid_sec_slice[NUM_SECTORS][fid_slices];
  TH2D* hadron_fid_hist[3];
  // fiducial

  // EC hists
  float EC_min = 0;
  float EC_max = 1;
  TH2D* EC_sampling_fraction =
      new TH2D("EC_sampling_fraction", "EC_sampling_fraction", BINS, p_min, p_max, BINS, EC_min, EC_max);
  TH2D* EC_sampling_fraction_cut =
      new TH2D("EC_sampling_fraction_cut", "EC_sampling_fraction_cut", BINS, p_min, p_max, BINS, EC_min, EC_max);
  TH1D* EC_hist[num_points];
  TH1D* EC_hist_cut[num_points];
  TH2D* Theta_vs_mom = new TH2D("Theta_vs_mom", "Theta_vs_mom", BINS, p_min, p_max, BINS, 0, 100);
  TH2D* ECin_ECout = new TH2D("ECin_ECout", "ECin_ECout", BINS, 0.0, 0.5, BINS, 0.0, 0.5);
  // EC hists

  // Beam Position
  float x_y_min_max = 0.5;
  TH1D* Beam_Position_X = new TH1D("Beam_Position_X", "Beam_Position_X", BINS, -x_y_min_max, x_y_min_max);
  TH1D* Beam_Position_Y = new TH1D("Beam_Position_Y", "Beam_Position_Y", BINS, -x_y_min_max, x_y_min_max);
  TH1D* Beam_Position_Z = new TH1D("Beam_Position_Z", "Beam_Position_Z", BINS, -5.0, 5.0);

  TH2D* Beam_Position =
      new TH2D("Beam_Position", "Beam_Position", BINS, -x_y_min_max, x_y_min_max, BINS, -x_y_min_max, x_y_min_max);

  // Beam Position

  // Vertex
  TH1D* target_vertex_X = new TH1D("Target_vertex_X", "Target_vertex_X", BINS, -6.0, 6.0);
  TH1D* target_vertex_Y = new TH1D("Target_vertex_Y", "Target_vertex_Y", BINS, -6.0, 6.0);
  TH1D* target_vertex_Z = new TH1D("Target_vertex_Z", "Target_vertex_Z", BINS, -6.0, 6.0);

  TH2D* target_vertex_xy = new TH2D("Target_vertex_xy", "Target_vertex_xy", BINS, -6.0, 6.0, BINS, -6.0, 6.0);
  TH2D* target_vertex_zy = new TH2D("Target_vertex_zy", "Target_vertex_zy", BINS, -6.0, 6.0, BINS, -6.0, 6.0);
  TH2D* target_vertex_zx = new TH2D("Target_vertex_zx", "Target_vertex_zx", BINS, -6.0, 6.0, BINS, -6.0, 6.0);

  TH1D* Missing_Mass = new TH1D("Missing_Mass", "Missing Mass", BINS_MM, MM_min, MM_max);
  TH1D* Missing_Mass_square = new TH1D("Missing_Mass_square", "Missing Mass square", BINS_MM, MM_min, MM_max* MM_max);

  TH1D* Missing_Mass_strict = new TH1D("Missing_Mass_strict", "Missing Mass", BINS_MM, MM_min, MM_max);
  TH1D* Missing_Mass_square_strict =
      new TH1D("Missing_Mass_square_strict", "Missing Mass square", BINS_MM, MM_min, MM_max* MM_max);

  TH1D* Missing_Mass_pi0 = new TH1D("Missing_Mass_pi0", "Missing Mass #pi^{0}", BINS_MM, -3, 3);
  TH1D* Missing_Mass_square_pi0 = new TH1D("Missing_Mass_pi0_2", "MM^{2} #pi^{0}", BINS_MM, -3, 3);

  TH1D* Missing_Mass_2pi = new TH1D("Missing_Mass_2pi", "Missing Mass 2 #pi", BINS_MM, MM_min, MM_max);
  TH1D* Missing_Mass_square_2pi =
      new TH1D("Missing_Mass_square_2pi", "Missing Mass 2 #pi", BINS_MM, MM_min, MM_max* MM_max);

  TH1D* Missing_Mass_nutron_no2pi =
      new TH1D("Missing_Mass_nutron_no2pi", "Missing Mass removing 2 #pi", BINS_MM, MM_min, MM_max);

  TH1D* energy_no_cuts = new TH1D("Energy_no_cuts", "Scattered electron energy", 500, 0.0, 5.0);
  TH1D* energy_fid_cuts = new TH1D("Energy_fid_cuts", "Scattered electron energy after fiducial cuts", 500, 0.0, 5.0);
  TH1D* energy_channel_cuts =
      new TH1D("Energy_channel_cuts", "Scattered electron energy for N #pi^{+} events", 500, 0.0, 5.0);

  // W and Q^2
  void Fill_proton_WQ2(float W, float Q2);
  void Fill_P_PI0(float W, float Q2);
  void Fill_NeutronPip_WQ2(float W, float Q2);
  void Fill_MM_WQ2(float W, float Q2);
  void Fill_channel_WQ2(float W, float Q2, TLorentzVector e_prime, float mm, float mm2, int sec);
  void Fill_single_proton_WQ2(float W, float Q2);
  void WvsQ2_Fill(float W, float Q2);
  void Fill_pion_WQ2(float W, float Q2);
  void WvsQ2_Write();
  void WvsQ2_binned_Write();

  // P and E
  void MomVsBeta_Fill_pos(float P, float Beta);
  void MomVsBeta_Fill_neg(float P, float Beta);
  void MomVsBeta_Fill_neutral(float P, float Beta);
  void Fill_proton_ID_P(float p, float beta);
  void Fill_Pi_ID_P(float p, float beta);
  void Fill_proton_Pi_ID_P(float p, float beta);
  void MomVsBeta_Fill(float P, float Beta);
  void Photon_flux_Fill(float photon_flux);
  void MomVsBeta_Write();

  // Missing Mass
  void Fill_Missing_Mass(float miss_mass);
  void Fill_Missing_Mass(float mm, float mm2);
  void Fill_Missing_Mass_strict(float mm, float mm2);
  // void Fill_Mass(float mass);
  void Fill_Missing_Mass_square(float miss_mass_2);
  void Write_Missing_Mass();

  void Fill_Missing_Mass_pi0(float mm, float mm2);
  void Fill_Missing_Mass_twoPi(float mm, float mm2);
  void Fill_W_Missing_Mass(float W, float mm, float mm2);

  // Delta T
  void Fill_deltat_P(float momentum, float delta_t);
  void Fill_deltat_P_PID(float momentum, float delta_t);
  void Fill_deltat_PIP(float momentum, float delta_t);
  void Fill_deltat_PIP_PID(float momentum, float delta_t);
  void Fill_deltat_PIM(float momentum, float delta_t);
  void Fill_deltat_PIM_PID(float momentum, float delta_t);
  void Fill_deltat_electron(float momentum, float delta_t);
  void Fill_deltat_electron_PID(float momentum, float delta_t);
  void Fill_deltat_kp(float momentum, float delta_t);
  void Fill_deltat_kp_PID(float momentum, float delta_t);
  void delta_t_slice_fit();
  void delta_t_Write();
  void delta_t_Fill(float momentum, int charge, float delta_t_proton, float delta_t_pip, float delta_t_electron);
  void delta_t_slices_Write();
  void delta_t_sec_pad(float momentum, int charge, float delta_t_proton, float delta_t_pip, float delta_t_electron,
                       int sc_sector, int sc_paddle);
  void delta_t_sec_pad_Write();
  void delta_T_canvas();

  // cc hist
  void CC_fill(int cc_sector, int cc_segment, int cc_pmt, int cc_nphe, float theta_cc);
  void CC_Write();
  void Theta_CC_Write();
  void theta_cc_slice_fit();
  void CC_canvas();

  // fiducial hist
  void Fill_electron_fid(float theta, float phi, int sector);
  void Fill_hadron_fid(float theta, float phi, int sector, int id);
  void Fid_Write();
  void fid_canvas();

  // EC hists
  void makeHists_EC();
  void EC_slices_Write();
  void EC_fill(float etot, float momentum);
  void EC_inout(float Ein, float Eout);
  void EC_cut_fill(float etot, float momentum);
  void EC_slice_fit();
  void EC_Write();

  void TM_Fill(float momentum, float theta);

  // Beam Position
  void Fill_Beam_Position(float vertex_x, float vertex_y, float vertex_z);
  void Beam_Position_Write();

  void Fill_Target_Vertex(float vertex_x, float vertex_y, float vertex_z);
  void Target_Vertex_Write();

  void Fill_E_Prime(TLorentzVector e_prime);
  void Fill_E_Prime_fid(TLorentzVector e_prime);
  void Fill_E_Prime_channel(TLorentzVector e_prime);
  void E_Prime_Write();
};

#endif
