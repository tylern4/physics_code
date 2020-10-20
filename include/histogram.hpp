/***************************************/
/*									   */
/*  Created by Nick Tyler              */
/*	University Of South Carolina       */
/***************************************/

#ifndef HISTOGRAM_H
#define HISTOGRAM_H
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <unordered_map>
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "THnSparse.h"
#include "TThread.h"
#include "branches.hpp"
#include "color.hpp"
#include "constants.hpp"
#include "fits.hpp"
#include "makeHeader.hpp"
#include "missing_mass.hpp"
#include "reaction.hpp"

#define DIMENSIONS 4
#define W_BINS 20
#define Q2_BINS 10
#define THETA_BINS 100
#define PHI_BINS 100
#define NUM_SECTORS 6
#define NDIMS_PIP_N 7
#define SC_PADDLE_NUM 48

#define BINS 500
#define BINS_MM 250

#define FID_SLICES 20

#define N_SIGMA 3
#define NUM_POINTS 20

using TH2F_ptr = std::shared_ptr<TH2F>;
using TH1F_ptr = std::shared_ptr<TH1F>;

class Histogram {
 protected:
  std::mutex mutex;
  std::shared_ptr<TFile> RootOutputFile;
  std::shared_ptr<TCanvas> def;

  /*
Garys binning
  W: 1.1 to 1.9 GeV in 32 bins
  Q^2: 0.4 to 1.0 GeV in 3 bins
  cos(theta) range: -1 to 1 in 10 bins
  phi range: 0 to 360 (degrees) in 6, 8, or 9 bins
*/
  //////////////// W , Q2, Theta_star_pip, Phi_star_pip
  //////////////// int nbins[DIMENSIONS] = {28, 5, 10, 12};
  //////////////// float xmin[DIMENSIONS] = {1.1, 1.0, -1.0, 0.0};
  //////////////// float xmax[DIMENSIONS] = {1.8, 2.0, 1.0, (M_PI * 2.0)};

  int nbins[DIMENSIONS] = {24, 4, 10, 20};
  double xmin[DIMENSIONS] = {1.2, 1.5, -1.0, 0};
  double xmax[DIMENSIONS] = {1.8, 3.5, 1.0, (M_PI * 2.0)};

  std::unique_ptr<THn> ndhist_nPip;
  std::unique_ptr<THn> ndhist_protPi0;

  TH1F* final_hists[24][4][10][20];

  float p_min = 0.0;
  float p_max = 5.0;

  float w_binned_min = 1.3;
  float w_binned_max = 1.8;
  float q2_binned_min = 1.0;
  float q2_binned_max = 4.5;

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

  TH2F_ptr WvsQ2_hist = std::make_shared<TH2F>("WvsQ2_hist", "W vs Q^{2}", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH2F_ptr ThetaVsPhi_hist = std::make_shared<TH2F>("ThetaVsPhi_hist", "Theta vs Phi", BINS, 0, PI, BINS, 0, 2 * PI);
  TH2F_ptr CosThetaVsPhi_hist =
      std::make_shared<TH2F>("CosThetaVsPhi_hist", "Cos Theta vs Phi", BINS, -1.0, 1.0, BINS, 0, 2 * PI);
  TH1F_ptr W_hist = std::make_shared<TH1F>("W", "W", BINS, w_min, w_max);

  std::vector<TH2F_ptr> WvsQ2_sec;
  std::vector<TH1F_ptr> W_sec;

  std::vector<TH2F_ptr> WvsQ2_channel_sec;
  std::vector<TH1F_ptr> W_channel_sec;

  TH1F_ptr Q2_hist = std::make_shared<TH1F>("Q2", "Q^{2}", BINS, q2_min, q2_max);
  TH1F_ptr E_prime_hist = std::make_shared<TH1F>("E_prime", "Scattered Electron Energy", BINS, 0.0, 5.0);
  TH1F_ptr photon_flux_hist = std::make_shared<TH1F>("photon_flux", "Photon Flux", 100, -0.0001, 0.001);
  TH2F_ptr Q2_vs_xb = std::make_shared<TH2F>("Q2_vs_xb", "Q^{2} vs x_{b}", BINS, 0.1, 0.6, BINS, 1.0, 3.5);
  TH2F_ptr WvsQ2_proton =
      std::make_shared<TH2F>("WvsQ2_proton", "W vs Q^{2} P", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1F_ptr W_proton = std::make_shared<TH1F>("W_proton", "W P", BINS, w_min, w_max);
  TH1F_ptr Q2_proton = std::make_shared<TH1F>("Q2_proton", "Q^{2} P", BINS, q2_min, q2_max);
  TH2F_ptr WvsQ2_pion =
      std::make_shared<TH2F>("WvsQ2_pion", "W vs Q^{2} #pi^{+} only", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1F_ptr W_pion = std::make_shared<TH1F>("W_pion", "W #pi^{+} only", BINS, w_min, w_max);
  TH1F_ptr Q2_pion = std::make_shared<TH1F>("Q2_pion", "Q^{2} #pi^{+} only", BINS, q2_min, q2_max);
  TH2F_ptr WvsQ2_NeutronPip =
      std::make_shared<TH2F>("WvsQ2_NeutronPip", "W vs Q^{2} N #pi^{+}", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH2F_ptr WvsMM_NeutronPip =
      std::make_shared<TH2F>("WvsMM_NeutronPip", "W vs MM N #pi^{+}", BINS, w_min, w_max, BINS, -q2_max, q2_max);
  TH2F_ptr WvsMM2_NeutronPip =
      std::make_shared<TH2F>("WvsMM2_NeutronPip", "W vs MM^{2} N #pi^{+}", BINS, w_min, w_max, BINS, -q2_max, q2_max);
  TH1F_ptr W_NeutronPip = std::make_shared<TH1F>("W_NeutronPip", "W N #pi^{+}", BINS, w_min, w_max);
  TH1F_ptr Q2_NeutronPip = std::make_shared<TH1F>("Q2_NeutronPip", "Q^{2} N #pi^{+}", BINS, q2_min, q2_max);
  TH2F_ptr WvsQ2_MM =
      std::make_shared<TH2F>("WvsQ2_MM", "W vs Q^{2} mm N cut", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1F_ptr W_MM = std::make_shared<TH1F>("W_MM", "W mm N cut", BINS, w_min, w_max);
  TH1F_ptr Q2_MM = std::make_shared<TH1F>("Q2_MM", "Q^{2} mm N cut", BINS, q2_min, q2_max);
  TH2F_ptr WvsQ2_channel =
      std::make_shared<TH2F>("WvsQ2_channel", "W vs Q^{2} #pi^{+} N", BINS, w_min, w_max, BINS, q2_min, q2_max);

  TH2F_ptr ThetaVsPhi_channel =
      std::make_shared<TH2F>("ThetaVsPhi_channel", "Theta vs Phi #pi^{+} N", 100, 0, PI, 100, 0, 2 * PI);
  TH2F_ptr CosThetaVsPhi_channel =
      std::make_shared<TH2F>("CosThetaVsPhi_channel", "Cos Theta vs Phi #pi^{+} N", 100, -1.0, 1.0, 100, 0, 2 * PI);

  TH1F_ptr W_channel = std::make_shared<TH1F>("W_channel", "W #pi^{+} N", BINS, w_min, w_max);
  TH1F_ptr Q2_channel = std::make_shared<TH1F>("Q2_channel", "Q^{2} #pi^{+} N", BINS, q2_min, q2_max);
  TH2F_ptr WvsQ2_Ppi0 =
      std::make_shared<TH2F>("WvsQ2_Ppi0", "W vs Q^{2} P #pi^{0}", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1F_ptr W_Ppi0 = std::make_shared<TH1F>("W_Ppi0", "W P #pi^{0}", BINS, w_min, w_max);
  TH1F_ptr Q2_Ppi0 = std::make_shared<TH1F>("Q2_Ppi0", "Q^{2} P #pi^{0}", BINS, q2_min, q2_max);
  TH2F_ptr WvsQ2_elastic =
      std::make_shared<TH2F>("WvsQ2_elastic", "W vs Q^{2} P", BINS, w_min, 2.0, BINS, q2_min, q2_max);
  TH1F_ptr W_elastic = std::make_shared<TH1F>("W_elastic", "W P", BINS, w_min, 2.0);
  TH1F_ptr Q2_elastic = std::make_shared<TH1F>("Q2_elastic", "Q^{2} P", BINS, q2_min, q2_max);

  TH1F_ptr elastic_MM = std::make_shared<TH1F>("Elastic_MM", "Elastic_MM", BINS, -0.2, 0.2);
  TH1F_ptr elastic_phi = std::make_shared<TH1F>("Elastic_phi", "Elastic_phi", BINS, 2, 4);
  TH2F_ptr elastic_thetaVsP = std::make_shared<TH2F>("elastic_thetaVsP", "elastic_thetaVsP", 500, 0.0, 3.5, 500, 0, 90);

  TH2F_ptr WvsQ2_binned = std::make_shared<TH2F>("WvsQ2_hist_binned", "W vs Q^{2} binned", W_BINS, w_binned_min,
                                                 w_binned_max, Q2_BINS, q2_binned_min, q2_binned_max);

  std::vector<TH1F_ptr> W_binned;
  std::vector<TH1F_ptr> Q2_binned;
  std::vector<TH1F_ptr> Missing_Mass_WBinned;
  std::vector<TH1F_ptr> Missing_Mass_WBinned_square;
  std::vector<Fits*> Fit_Missing_Mass_WBinned;
  // W and Q^2

  // P and E
  float b_min = 0.1;
  float b_max = 1.2;
  TH2F_ptr MomVsBeta_hist =
      std::make_shared<TH2F>("MomVsBeta", "Momentum versus #beta", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2F_ptr MomVsBeta_hist_pos =
      std::make_shared<TH2F>("MomVsBeta_pos", "Momentum versus #beta Positive", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2F_ptr MomVsBeta_hist_neg =
      std::make_shared<TH2F>("MomVsBeta_neg", "Momentum versus #beta Negative", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2F_ptr MomVsBeta_hist_neutral = std::make_shared<TH2F>("MomVsBeta_Fill_neutral", "Momentum versus #beta neutral",
                                                           BINS, p_min, p_max, BINS, b_min, b_max);
  TH1F_ptr Mom = std::make_shared<TH1F>("Momentum", "Momentum", BINS, 0, 2.5);
  TH1F_ptr Energy_hist = std::make_shared<TH1F>("Energy_hist", "Energy_hist", BINS, 0.0, 2.5);
  TH2F_ptr MomVsBeta_proton_ID =
      std::make_shared<TH2F>("MomVsBeta_proton_ID", "Momentum versus #beta p", BINS, p_min, p_max, BINS, b_min, b_max);
  TH2F_ptr MomVsBeta_Pi_ID = std::make_shared<TH2F>("MomVsBeta_Pi_ID", "Momentum versus #beta #pi^{+}", BINS, p_min,
                                                    p_max, BINS, b_min, b_max);
  TH2F_ptr MomVsBeta_proton_Pi_ID = std::make_shared<TH2F>("MomVsBeta_proton_Pi_ID", "Momentum versus #beta P #pi^{+}",
                                                           BINS, p_min, p_max, BINS, b_min, b_max);
  // P and E

  // Delta T
  // Histogram declarations, fills, and write
  // j -> type: 0=>Proton,1=>Pip,2=>Electron
  // jj -> Fit point
  float Dt_min = -10;
  float Dt_max = 10;
  TH1F* delta_t_hist[3][NUM_POINTS];
  float bin_width = (p_max - p_min) / NUM_POINTS;

  TH2F* delta_t_sec_pad_hist[3][NUM_SECTORS][SC_PADDLE_NUM];
  TH2F_ptr delta_t_mass_P = std::make_shared<TH2F>("delta_t_mass_P", "#Deltat assuming mass of proton", BINS, p_min,
                                                   p_max, BINS, Dt_min, Dt_max);
  TH2F_ptr delta_t_mass_P_PID =
      std::make_shared<TH2F>("delta_t_mass_P_PID", "#Deltat assuming mass of proton with PID proton", BINS, p_min,
                             p_max, BINS, Dt_min, Dt_max);

  TH2F_ptr delta_t_mass_PIP = std::make_shared<TH2F>("delta_t_mass_PIP", "#Deltat assuming mass of #pi^{+}", BINS,
                                                     p_min, p_max, BINS, Dt_min, Dt_max);
  TH2F_ptr delta_t_mass_PIP_PID =
      std::make_shared<TH2F>("delta_t_mass_PIP_PID", "#Deltat assuming mass of #pi^{+} with PID #pi^{+}", BINS, p_min,
                             p_max, BINS, Dt_min, Dt_max);

  TH2F_ptr delta_t_mass_PIM = std::make_shared<TH2F>("delta_t_mass_PIM", "#Deltat assuming mass of #pi^{-}", BINS,
                                                     p_min, p_max, BINS, Dt_min, Dt_max);
  TH2F_ptr delta_t_mass_PIM_PID =
      std::make_shared<TH2F>("delta_t_mass_PIM_PID", "#Deltat assuming mass of #pi^{-} with PID #pi^{-}", BINS, p_min,
                             p_max, BINS, Dt_min, Dt_max);

  TH2F_ptr delta_t_mass_electron = std::make_shared<TH2F>("delta_t_mass_electron", "#Deltat assuming mass of e^{-}",
                                                          BINS, p_min, p_max, BINS, Dt_min, Dt_max);
  TH2F_ptr delta_t_mass_electron_PID =
      std::make_shared<TH2F>("delta_t_mass_electron_PID", "#Deltat assuming mass of e^{-} with PID e^{-}", BINS, p_min,
                             p_max, BINS, Dt_min, Dt_max);

  TH2F_ptr delta_t_mass_kp = std::make_shared<TH2F>("delta_t_mass_kp", "#Deltat assuming mass of k^{+}", BINS, p_min,
                                                    p_max, BINS, Dt_min, Dt_max);
  TH2F_ptr delta_t_mass_kp_PID = std::make_shared<TH2F>(
      "delta_t_mass_kp_PID", "#Deltat assuming mass of k^{+} with PID k^{+}", BINS, p_min, p_max, BINS, Dt_min, Dt_max);
  // Delta T

  // cc hist
  int bins_CC = 50;
  float CC_min = 0;
  float CC_max = 250;
  std::unordered_map<int, std::string> L_R_C = {{0, "both"}, {1, "right"}, {2, "left"}};
  static const int segment = 18;
  static const int PMT = 3;

  TH1F_ptr cc_hist[NUM_SECTORS][segment][PMT];
  TH1F_ptr cc_hist_allSeg[NUM_SECTORS][PMT];

  TH2F_ptr fid_xy_hist;
  std::vector<TH2F_ptr> fid_xy;  //[NUM_SECTORS];

  TH2F_ptr Theta_CC = std::make_shared<TH2F>("Theta_CC", "Theta_CC", 20, 0.0, 20.0, 60, 0.0, 60.0);
  std::vector<TH2F_ptr> Theta_CC_Sec;      //[NUM_SECTORS];
  std::vector<TH2F_ptr> Theta_CC_Sec_cut;  //[NUM_SECTORS];
  // cc hist

  // fiducial
  float theta_min = 0;
  float theta_max = 60;
  float phi_min = -360 / 2.0;
  float phi_max = 360 / 2.0;

  static const int start_slice = 0;
  std::vector<TH2F_ptr> electron_fid_sec_hist;  //[NUM_SECTORS];
  TH1F* electron_fid_sec_slice[NUM_SECTORS][FID_SLICES];
  TH2F_ptr electron_fid_hist =
      std::make_shared<TH2F>("electron_fid", "electron_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
  TH2F_ptr neutron_fid_hist = std::make_shared<TH2F>("neutron_fid", "neutron_fid", BINS, -360, 360, BINS, -360, 360);

  TH2F* hadron_fid_sec_hist[3][NUM_SECTORS];
  TH1F* hadron_fid_sec_slice[NUM_SECTORS][FID_SLICES];
  TH2F* hadron_fid_hist[3];

  std::vector<TH2F_ptr> cerenkov_fid;

  // fiducial

  // EC hists
  float EC_min = 0;
  float EC_max = 1;
  TH2F_ptr EC_sampling_fraction =
      std::make_shared<TH2F>("EC_sampling_fraction", "EC_sampling_fraction", BINS, p_min, p_max, BINS, EC_min, EC_max);
  TH2F_ptr EC_sampling_fraction_cut = std::make_shared<TH2F>("EC_sampling_fraction_cut", "EC_sampling_fraction_cut",
                                                             BINS, p_min, p_max, BINS, EC_min, EC_max);
  TH1F* EC_hist[NUM_POINTS];
  TH1F* EC_hist_cut[NUM_POINTS];
  TH2F_ptr Theta_vs_mom = std::make_shared<TH2F>("Theta_vs_mom", "Theta_vs_mom", BINS, p_min, p_max, BINS, 0, 100);
  TH2F_ptr ECin_ECout = std::make_shared<TH2F>("ECin_ECout", "ECin_ECout", BINS, 0.0, 0.5, BINS, 0.0, 0.5);

  TH1F_ptr EC_tot_energy = std::make_shared<TH1F>("EC_tot_energy", "EC_tot_energy", BINS, 0, 1.5);
  TH2F_ptr EC_etot_vs_P = std::make_shared<TH2F>("EC_etot_vs_P", "EC_etot_vs_P", BINS, 0, 5.0, BINS, 0, 1.5);
  // EC hists

  // Beam Position
  float x_y_min_max = 0.5;
  TH1F_ptr Beam_Position_X =
      std::make_shared<TH1F>("Beam_Position_X", "Beam_Position_X", BINS, -x_y_min_max, x_y_min_max);
  TH1F_ptr Beam_Position_Y =
      std::make_shared<TH1F>("Beam_Position_Y", "Beam_Position_Y", BINS, -x_y_min_max, x_y_min_max);
  TH1F_ptr Beam_Position_Z = std::make_shared<TH1F>("Beam_Position_Z", "Beam_Position_Z", 10 * BINS, -10, 15.0);

  TH2F_ptr Beam_Position = std::make_shared<TH2F>("Beam_Position", "Beam_Position", BINS, -x_y_min_max, x_y_min_max,
                                                  BINS, -x_y_min_max, x_y_min_max);

  // Beam Position

  // Vertex
  TH1F_ptr target_vertex_X = std::make_shared<TH1F>("Target_vertex_X", "Target_vertex_X", BINS, -6.0, 6.0);
  TH1F_ptr target_vertex_Y = std::make_shared<TH1F>("Target_vertex_Y", "Target_vertex_Y", BINS, -6.0, 6.0);
  TH1F_ptr target_vertex_Z = std::make_shared<TH1F>("Target_vertex_Z", "Target_vertex_Z", BINS, -6.0, 6.0);

  TH2F_ptr target_vertex_xy =
      std::make_shared<TH2F>("Target_vertex_xy", "Target_vertex_xy", BINS, -6.0, 6.0, BINS, -6.0, 6.0);
  TH2F_ptr target_vertex_zy =
      std::make_shared<TH2F>("Target_vertex_zy", "Target_vertex_zy", BINS, -6.0, 6.0, BINS, -6.0, 6.0);
  TH2F_ptr target_vertex_zx =
      std::make_shared<TH2F>("Target_vertex_zx", "Target_vertex_zx", BINS, -6.0, 6.0, BINS, -6.0, 6.0);

  TH2F_ptr target_vertex_vz_phi =
      std::make_shared<TH2F>("target_vertex_vz_phi", "target_vertex_vz_phi", BINS, -10, 15.0, BINS, -200.0, 200.0);

  TH1F_ptr Missing_Mass = std::make_shared<TH1F>("Missing_Mass", "Missing Mass", BINS_MM, MM_min, MM_max);
  TH1F_ptr Missing_Mass_small = std::make_shared<TH1F>("Missing_Mass_small", "e(p,#pi^{+} X)e'", BINS_MM, 0.8, 1.3);
  std::vector<TH1F_ptr> Missing_Mass_small_sec;
  std::vector<TH1F_ptr> Missing_Mass_Sq_small_sec;
  TH1F_ptr Missing_Mass_square =
      std::make_shared<TH1F>("Missing_Mass_square", "Missing Mass square", BINS_MM, 0.7, 1.5);

  TH1F_ptr Missing_Mass_strict = std::make_shared<TH1F>("Missing_Mass_strict", "Missing Mass", BINS_MM, MM_min, MM_max);
  TH1F_ptr Missing_Mass_square_strict =
      std::make_shared<TH1F>("Missing_Mass_square_strict", "Missing Mass square", BINS_MM, 0.7, 1.5);

  TH1F_ptr Missing_Mass_pi0 = std::make_shared<TH1F>("Missing_Mass_pi0", "Missing Mass #pi^{0}", BINS_MM, -1, 1);
  TH1F_ptr Missing_Mass_square_pi0 = std::make_shared<TH1F>("Missing_Mass_pi0_2", "MM^{2} #pi^{0}", BINS_MM, -1, 1);

  TH1F_ptr Mass_pi0 = std::make_shared<TH1F>("Mass_pi0", "Mass #pi^{0}", BINS_MM, 0, 0.5);
  TH1F_ptr Mass_square_pi0 = std::make_shared<TH1F>("Mass_pi0_2", "#pi^{0} mass^{2}", BINS_MM, 0, 0.05);

  TH1F_ptr Mass_eta = std::make_shared<TH1F>("Mass_eta", "Mass #eta", BINS_MM, 0.2, 1.0);
  TH1F_ptr Mass_square_eta = std::make_shared<TH1F>("Mass_eta_2", "#eta mass^{2}", BINS_MM, 0.2, 1.0);

  TH1F_ptr Missing_Mass_pi0_otherCut =
      std::make_shared<TH1F>("Missing_Mass_pi0_otherCut", "Missing Mass #pi^{0}", BINS_MM, -1, 1);
  TH1F_ptr Missing_Mass_square_pi0_otherCut =
      std::make_shared<TH1F>("Missing_Mass_pi0_2_otherCut", "MM^{2} #pi^{0}", BINS_MM, -1, 1);

  TH1F_ptr Mass_pi0_otherCut = std::make_shared<TH1F>("Mass_pi0_otherCut", "Mass #pi^{0}", BINS_MM, 0, 0.5);
  TH1F_ptr Mass_square_pi0_otherCut =
      std::make_shared<TH1F>("Mass_pi0_2_otherCut", "#pi^{0} mass^{2}", BINS_MM, 0, 0.05);

  TH1F_ptr Missing_Mass_2pi = std::make_shared<TH1F>("Missing_Mass_2pi", "Missing Mass 2 #pi", BINS_MM, MM_min, MM_max);
  TH1F_ptr Missing_Mass_square_2pi =
      std::make_shared<TH1F>("Missing_Mass_square_2pi", "Missing Mass 2 #pi", BINS_MM, MM_min, MM_max* MM_max);

  TH1F_ptr energy_no_cuts = std::make_shared<TH1F>("Energy_no_cuts", "Scattered electron energy", 500, 0.0, 5.0);
  TH1F_ptr energy_fid_cuts =
      std::make_shared<TH1F>("Energy_fid_cuts", "Scattered electron energy after fiducial cuts", 500, 0.0, 5.0);
  TH1F_ptr energy_channel_cuts =
      std::make_shared<TH1F>("Energy_channel_cuts", "Scattered electron energy for N #pi^{+} events", 500, 0.0, 5.0);

  void makeHists_fid();
  void makeHists_deltat();
  void makeHists_CC();
  void makeHists_WvsQ2();

 public:
  Histogram();
  Histogram(const std::string& output_file);
  ~Histogram();

  void FillEvent(const std::shared_ptr<Reaction>& event);

  void Write(const std::string& output_file);
  void Write();

  // W and Q^2
  void Fill_proton_WQ2(float W, float Q2);
  void Fill_P_PI0(float W, float Q2);
  void Fill_NeutronPip_WQ2(float W, float Q2, float MM, float MM2);
  void Fill_MM_WQ2(float W, float Q2);
  void Fill_channel_WQ2(const std::shared_ptr<Reaction>& event);
  void Fill_elastic(const std::shared_ptr<Reaction>& event);
  void WvsQ2_Fill(float W, float Q2, int sector);
  void Fill_pion_WQ2(float W, float Q2);
  void WvsQ2_Write();
  void WvsQ2_sec_Write();
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
  void Fill_Missing_Mass(const std::shared_ptr<Reaction>& event);
  void Fill_Missing_Mass_strict(float mm, float mm2);
  // void Fill_Mass(float mass);
  void Fill_Missing_Mass_square(float miss_mass_2);

  void Fill_Mass_photons(std::shared_ptr<Reaction> _e);
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
  void Fill_electron_fid(const std::shared_ptr<Branches>& _data, const std::shared_ptr<Reaction>& _r);
  void Fill_neutron_fid(float theta, float phi, int sector);
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
  void Fill_Beam_Position(const std::shared_ptr<Branches>& _data);
  void Beam_Position_Write();
  void Fill_Target_Vertex(float px, float py, float pz, float vertex_x, float vertex_y, float vertex_z);
  void Target_Vertex_Write();
  void Fill_E_Prime(const LorentzVector& e_prime);
  void Fill_E_Prime_fid(const LorentzVector& e_prime);
  void Fill_E_Prime_channel(const LorentzVector& e_prime);
  void E_Prime_Write();

  void Fill_ND(const std::shared_ptr<Reaction>& event);
};

class mcHistogram : public Histogram {
 private:
  std::unique_ptr<THnSparse> ndhist_mc;
  TH2F_ptr WvsQ2_MC =
      std::make_shared<TH2F>("WvsQ2_MC", "W vs Q^{2} #pi^{+} N", BINS, w_min, w_max, BINS, q2_min, q2_max);
  TH1F_ptr W_MC = std::make_shared<TH1F>("W_MC", "W #pi^{+} N", BINS, w_min, w_max);

  TH2F_ptr WvsQ2_binned_MC = std::make_shared<TH2F>("WvsQ2_hist_binned_MC", "W vs Q^{2} binned", W_BINS, w_binned_min,
                                                    w_binned_max, Q2_BINS, q2_binned_min, q2_binned_max);

  TH1F_ptr W_binned_MC[Q2_BINS];
  TH1F_ptr delta_p[4];
  TH1F_ptr delta_p_electron[4];
  TH2F_ptr delta_px_py_electron =
      std::make_shared<TH2F>("delta_px_py_electron", "#DeltaP_x vs #DeltaP_y", 500, -0.1, 0.1, 500, -0.1, 0.1);

 public:
  mcHistogram() : Histogram() {}
  mcHistogram(const std::string& output_file) : Histogram(output_file) { makeMCHists(); }
  ~mcHistogram();

  // W and Q^2
  void makeMCHists();
  void Fill_WQ2_MC(const std::shared_ptr<MCReaction>& _e);
  void Fill(const std::shared_ptr<MCReaction>& _e);
  void Fill_P(const std::shared_ptr<Branches>& d);
  void Write();
  void Write_DeltaP();
  void WvsQ2_MC_Write();
};

#endif
