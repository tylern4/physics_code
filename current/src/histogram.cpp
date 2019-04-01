/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram() {
  def = new TCanvas("def");
  makeHists_WvsQ2();
  makeHists_deltat();
  makeHists_EC();
  makeHists_CC();
  makeHists_fid();
  hadron_fid_hist[0] = new TH2D("hadron_fid", "hadron_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
  hadron_fid_hist[1] = new TH2D("proton_fid", "proton_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
  hadron_fid_hist[2] = new TH2D("pip_fid", "pip_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
}

Histogram::Histogram(std::string output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  def = new TCanvas("def");
  makeHists_WvsQ2();
  makeHists_deltat();
  makeHists_EC();
  makeHists_CC();
  makeHists_fid();
  hadron_fid_hist[0] = new TH2D("hadron_fid", "hadron_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
  hadron_fid_hist[1] = new TH2D("proton_fid", "proton_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
  hadron_fid_hist[2] = new TH2D("pip_fid", "pip_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
}

Histogram::~Histogram() {}
void Histogram::Write(const std::string &output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  Write();
}
void Histogram::Write() {
  std::cout << GREEN << "\nFitting" << DEF << std::endl;
  // Start of cuts
  auto MM_neutron_cut = std::make_unique<Fits>();
  MM_neutron_cut->FitMissMass(Missing_Mass.get());

  auto MissingMassSquare_cut = std::make_unique<Fits>();
  MissingMassSquare_cut->Set_max(1.1);
  MissingMassSquare_cut->Set_min(0.7);
  MissingMassSquare_cut->FitBreitWigner(Missing_Mass_square.get());
  MissingMassSquare_cut->Get_sigma();

  RootOutputFile->cd();
  std::cerr << BOLDBLUE << "EC_Write()" << DEF << std::endl;
  TDirectory *EC_folder = RootOutputFile->mkdir("EC_hists");
  EC_folder->cd();
  EC_Write();
  std::cerr << BOLDBLUE << "EC_slices()" << DEF << std::endl;
  TDirectory *EC_slices = RootOutputFile->mkdir("EC_slices");
  EC_slices->cd();
  EC_slices_Write();
  std::cerr << BOLDBLUE << "Beam_Position()" << DEF << std::endl;
  TDirectory *Beam_Folder = RootOutputFile->mkdir("Beam Position");
  Beam_Folder->cd();
  Beam_Position_Write();
  Target_Vertex_Write();
  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  WvsQ2_folder->cd();
  WvsQ2_Write();
  TDirectory *WvsQ2_sec_folder = RootOutputFile->mkdir("W_vs_Q2_sec");
  WvsQ2_sec_folder->cd();
  WvsQ2_sec_Write();
  TDirectory *W_Q2_binned = RootOutputFile->mkdir("W_Q2_binned");
  W_Q2_binned->cd();
  WvsQ2_binned_Write();
  std::cerr << BOLDBLUE << "MomVsBeta_Fill()" << DEF << std::endl;
  TDirectory *MomVsBeta_folder = RootOutputFile->mkdir("Momentum vs beta");
  MomVsBeta_folder->cd();
  MomVsBeta_Write();
  std::cerr << BOLDBLUE << "Write_Missing_Mass()" << DEF << std::endl;
  // Missing Mass Write
  TDirectory *MissMass = RootOutputFile->mkdir("Missing_Mass");
  MissMass->cd();
  Write_Missing_Mass();
  std::cerr << BOLDBLUE << "delta_t_Write()" << DEF << std::endl;
  // Delta T Write
  TDirectory *DeltaT = RootOutputFile->mkdir("Delta_T");
  DeltaT->cd();
  delta_t_Write();
  std::cerr << BOLDBLUE << "delta_t_slices_Write()" << DEF << std::endl;
  TDirectory *DeltaT_slices = RootOutputFile->mkdir("Delta_T_slices");
  DeltaT_slices->cd();
  delta_t_slices_Write();
  std::cerr << BOLDBLUE << "delta_t_sec_pad_Write()" << DEF << std::endl;
  TDirectory *DeltaT_sec_pad = RootOutputFile->mkdir("Delta_T_sec_pad");
  DeltaT_sec_pad->cd();
  delta_t_sec_pad_Write();
  std::cerr << BOLDBLUE << "delta_T_canvas()" << DEF << std::endl;
  TDirectory *Delta_T_canvases = RootOutputFile->mkdir("Delta_T_canvases");
  Delta_T_canvases->cd();
  delta_T_canvas();
  std::cerr << BOLDBLUE << "Theta_CC_Write()" << DEF << std::endl;
  TDirectory *Theta_CC_hists = RootOutputFile->mkdir("Theta_CC_hists");
  Theta_CC_hists->cd();
  Theta_CC_Write();
  std::cerr << BOLDBLUE << "CC_Write()" << DEF << std::endl;
  TDirectory *CC_hists = RootOutputFile->mkdir("CC_hists");
  CC_hists->cd();
  CC_Write();
  std::cerr << BOLDBLUE << "CC_canvas()" << DEF << std::endl;
  TDirectory *CC_canvases = RootOutputFile->mkdir("CC_canvases");
  CC_canvases->cd();
  CC_canvas();
  std::cerr << BOLDBLUE << "Fid_Write()" << DEF << std::endl;
  TDirectory *Fid_cuts = RootOutputFile->mkdir("Fid_cuts");
  Fid_cuts->cd();
  Fid_Write();
  std::cerr << BOLDBLUE << "fid_canvas()" << DEF << std::endl;
  TDirectory *Fid_canvas = RootOutputFile->mkdir("Fid_canvas");
  Fid_canvas->cd();
  fid_canvas();
  std::cerr << BOLDBLUE << "E_Prime_Write()" << DEF << std::endl;
  TDirectory *E_Prime = RootOutputFile->mkdir("E_Prime");
  E_Prime->cd();
  E_Prime_Write();
  std::cerr << "MM:mean,+3,-3" << std::endl;
  std::cerr << MM_neutron_cut->Get_mean() << ",";
  std::cerr << MM_neutron_cut->Get_mean() + 3 * MM_neutron_cut->Get_sigma() << ",";
  std::cerr << MM_neutron_cut->Get_mean() - 3 * MM_neutron_cut->Get_sigma() << std::endl;
  std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::Write(const std::string &output_file, bool multi) {
  _multi = multi;
  RootOutputFile = std::make_unique<TFile>(output_file.c_str(), "RECREATE");
  RootOutputFile->cd();

  TDirectory *EC_folder = RootOutputFile->mkdir("EC_hists");
  EC_folder->cd();
  EC_Write();

  TDirectory *Beam_Folder = RootOutputFile->mkdir("Beam Position");
  Beam_Folder->cd();
  Beam_Position_Write();
  Target_Vertex_Write();

  TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  WvsQ2_folder->cd();
  WvsQ2_Write();

  TDirectory *WvsQ2_sec_folder = RootOutputFile->mkdir("W_vs_Q2_sec");
  WvsQ2_sec_folder->cd();
  WvsQ2_sec_Write();
  TDirectory *W_Q2_binned = RootOutputFile->mkdir("W_Q2_binned");
  W_Q2_binned->cd();
  WvsQ2_binned_Write();

  TDirectory *MomVsBeta_folder = RootOutputFile->mkdir("Momentum vs beta");
  MomVsBeta_folder->cd();
  MomVsBeta_Write();

  TDirectory *MissMass = RootOutputFile->mkdir("Missing_Mass");
  MissMass->cd();
  Write_Missing_Mass();

  TDirectory *DeltaT = RootOutputFile->mkdir("Delta_T");
  DeltaT->cd();
  delta_t_Write();

  TDirectory *DeltaT_slices = RootOutputFile->mkdir("Delta_T_slices");
  DeltaT_slices->cd();
  delta_t_slices_Write();

  TDirectory *DeltaT_sec_pad = RootOutputFile->mkdir("Delta_T_sec_pad");
  DeltaT_sec_pad->cd();
  delta_t_sec_pad_Write();

  TDirectory *Theta_CC_hists = RootOutputFile->mkdir("Theta_CC_hists");
  Theta_CC_hists->cd();
  Theta_CC_Write();

  TDirectory *Fid_cuts = RootOutputFile->mkdir("Fid_cuts");
  Fid_cuts->cd();
  Fid_Write();

  RootOutputFile->Close();
}

// W and Q^2
void Histogram::makeHists_WvsQ2() {
  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    sprintf(hname, "W_vs_Q2_sec_%d", sec + 1);
    sprintf(htitle, "W vs Q^{2} Sector: %d", sec + 1);
    WvsQ2_sec[sec] = new TH2D(hname, htitle, BINS, w_min, w_max, BINS, q2_min, q2_max);

    sprintf(hname, "W_sec_%d", sec + 1);
    sprintf(htitle, "W Sector: %d", sec + 1);
    W_sec[sec] = new TH1D(hname, htitle, BINS, w_min, w_max);
  }

  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    sprintf(hname, "W_vs_Q2_channel_sec_%d", sec + 1);
    sprintf(htitle, "W vs Q^{2} N #pi^{+} Sector: %d", sec + 1);
    WvsQ2_channel_sec[sec] = new TH2D(hname, htitle, BINS, w_min, w_max, BINS, q2_min, q2_max);

    sprintf(hname, "W_channel_sec_%d", sec + 1);
    sprintf(htitle, "W N #pi^{+} Sector: %d", sec + 1);
    W_channel_sec[sec] = new TH1D(hname, htitle, BINS, w_min, w_max);
  }

  for (short y = 0; y < Q2_BINS; y++) {
    sprintf(hname, "W_%0.3f_%0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    sprintf(htitle, "W hist\nQ^{2} %0.3f %0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    W_binned[y] = new TH1D(hname, htitle, BINS, w_binned_min, w_binned_max);
  }
  for (short x = 0; x < W_BINS; x++) {
    sprintf(hname, "Q2_%0.3f_%0.3f", w_binned_min + (W_width * x), w_binned_min + (W_width * (x + 1)));
    sprintf(htitle, "Q^{2} hist\nW %0.3f %0.3f", w_binned_min + (W_width * x), w_binned_min + (W_width * (x + 1)));
    Q2_binned[x] = new TH1D(hname, htitle, BINS, q2_binned_min, q2_binned_max);

    sprintf(hname, "MM_W_%0.3f_%0.3f", w_binned_min + (W_width * x), w_binned_min + (W_width * (x + 1)));
    sprintf(htitle, "Missing Mass\nW %0.3f %0.3f", w_binned_min + (W_width * x), w_binned_min + (W_width * (x + 1)));
    Missing_Mass_WBinned[x] = new TH1D(hname, htitle, BINS, MM_min, MM_max);

    sprintf(hname, "MM2_W_%0.3f_%0.3f", w_binned_min + (W_width * x), w_binned_min + (W_width * (x + 1)));
    sprintf(htitle, "Missing Mass^{2}\nW %0.3f %0.3f", w_binned_min + (W_width * x),
            w_binned_min + (W_width * (x + 1)));
    Missing_Mass_WBinned_square[x] = new TH1D(hname, htitle, BINS, MM_min * MM_min, MM_max * MM_max);
  }
}

void Histogram::Fill_proton_WQ2(float W, float Q2) {
  WvsQ2_proton->Fill(W, Q2);
  W_proton->Fill(W);
  Q2_proton->Fill(Q2);
}
void Histogram::Fill_P_PI0(float W, float Q2) {
  WvsQ2_Ppi0->Fill(W, Q2);
  W_Ppi0->Fill(W);
  Q2_Ppi0->Fill(Q2);
}
void Histogram::Fill_NeutronPip_WQ2(float W, float Q2, float MM, float MM2) {
  WvsQ2_NeutronPip->Fill(W, Q2);
  W_NeutronPip->Fill(W);
  Q2_NeutronPip->Fill(Q2);

  WvsMM_NeutronPip->Fill(W, MM);
  WvsMM2_NeutronPip->Fill(W, MM2);
}

void Histogram::Fill_MM_WQ2(float W, float Q2) {
  WvsQ2_MM->Fill(W, Q2);
  W_MM->Fill(W);
  Q2_MM->Fill(Q2);
}

void Histogram::Fill_channel_WQ2(float W, float Q2, int sector, TLorentzVector e_prime, float mm, float mm2) {
  /*
  float x_pip_N[NDIMS_PIP_N];
  x_pip_N[0] = W;
  x_pip_N[1] = Q2;
  x_pip_N[2] = (float)sec;
  x_pip_N[3] = mm;
  x_pip_N[4] = mm2;
  x_pip_N[5] = e_prime.Theta();
  x_pip_N[6] = e_prime.Phi();


  // pip_N->Fill(x_pip_N);
  // std::cout << "Nope" << '\n';
  */
  E_prime_hist->Fill(e_prime.E());
  Q2_vs_xb->Fill(physics::xb_calc(Q2, e_prime.E()), Q2);

  WvsQ2_channel->Fill(W, Q2);
  W_channel->Fill(W);
  Q2_channel->Fill(Q2);

  WvsQ2_binned->Fill(W, Q2);

  for (int y = 0; y < Q2_BINS; y++) {
    if (q2_binned_min + (Q2_width * y) <= Q2 && q2_binned_min + (Q2_width * (y + 1)) >= Q2) {
      W_binned[y]->Fill(W);
      continue;
    }
  }

  for (int x = 0; x < W_BINS; x++) {
    if (w_binned_min + (W_width * x) <= W && w_binned_min + (W_width * (x + 1)) >= W) {
      Q2_binned[x]->Fill(Q2);
      continue;
    }
  }
  if (sector == 0 || sector > NUM_SECTORS) return;
  WvsQ2_channel_sec[sector - 1]->Fill(W, Q2);
  W_channel_sec[sector - 1]->Fill(W);
}

void Histogram::Fill_single_proton_WQ2(float W, float Q2) {
  WvsQ2_single_proton->Fill(W, Q2);
  W_single_proton->Fill(W);
  Q2_single_proton->Fill(Q2);
}

void Histogram::WvsQ2_Fill(float W, float Q2, int sector) {
  WvsQ2_hist->Fill(W, Q2);
  W_hist->Fill(W);
  Q2_hist->Fill(Q2);
  if (sector == 0 || sector > NUM_SECTORS) return;
  WvsQ2_sec[sector - 1]->Fill(W, Q2);
  W_sec[sector - 1]->Fill(W);
}

void Histogram::Fill_pion_WQ2(float W, float Q2) {
  WvsQ2_pion->Fill(W, Q2);
  W_pion->Fill(W);
  Q2_pion->Fill(Q2);
}

void Histogram::WvsQ2_Write() {
  /*
    pip_N->GetAxis(0)->SetTitle("W");
    pip_N->GetAxis(1)->SetTitle("Q^{2}");
    pip_N->GetAxis(2)->SetTitle("sector");
    pip_N->GetAxis(3)->SetTitle("MissMass");
    pip_N->GetAxis(4)->SetTitle("MissMass2");
    pip_N->GetAxis(5)->SetTitle("#theta");
    pip_N->GetAxis(6)->SetTitle("#phi");
    pip_N->Write();
  */
  WvsQ2_channel->SetXTitle("W (GeV)");
  WvsQ2_channel->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_channel->SetOption("COLZ");
  WvsQ2_channel->Write();

  W_channel->SetXTitle("W (GeV)");
  W_channel->Write();

  Q2_channel->SetXTitle("Q^{2} (GeV^{2})");
  Q2_channel->Write();

  E_prime_hist->SetXTitle("Energy (GeV)");
  E_prime_hist->Write();

  Q2_vs_xb->SetXTitle("x_{b}");
  Q2_vs_xb->SetYTitle("Q^{2}");
  Q2_vs_xb->SetOption("COLZ");
  Q2_vs_xb->Write();

  WvsQ2_hist->SetXTitle("W (GeV)");
  WvsQ2_hist->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_hist->SetOption("COLZ");
  WvsQ2_hist->Write();

  W_hist->SetXTitle("W (GeV)");
  W_hist->Write();

  Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
  Q2_hist->Write();

  WvsQ2_proton->SetXTitle("W (GeV)");
  WvsQ2_proton->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_proton->SetOption("COLZ");
  WvsQ2_proton->Write();

  W_proton->SetXTitle("W (GeV)");
  W_proton->Write();

  Q2_proton->SetXTitle("Q^{2} (GeV^{2})");
  Q2_proton->Write();

  WvsQ2_pion->SetXTitle("W (GeV)");
  WvsQ2_pion->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_pion->SetOption("COLZ");
  WvsQ2_pion->Write();

  W_pion->SetXTitle("W (GeV)");
  W_pion->Write();

  Q2_pion->SetXTitle("Q^{2} (GeV^{2})");
  Q2_pion->Write();

  WvsQ2_NeutronPip->SetXTitle("W (GeV)");
  WvsQ2_NeutronPip->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_NeutronPip->SetOption("COLZ");
  WvsQ2_NeutronPip->Write();

  WvsMM_NeutronPip->SetXTitle("W (GeV)");
  WvsMM_NeutronPip->SetYTitle("MM (GeV)");
  WvsMM_NeutronPip->SetOption("COLZ");
  WvsMM_NeutronPip->Write();

  WvsMM2_NeutronPip->SetXTitle("W (GeV)");
  WvsMM2_NeutronPip->SetYTitle("MM^{2} (GeV^{2})");
  WvsMM2_NeutronPip->SetOption("COLZ");
  WvsMM2_NeutronPip->Write();

  W_NeutronPip->SetXTitle("W (GeV)");
  W_NeutronPip->Write();

  Q2_NeutronPip->SetXTitle("Q^{2} (GeV^{2})");
  Q2_NeutronPip->Write();

  WvsQ2_single_proton->SetXTitle("W (GeV)");
  WvsQ2_single_proton->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_single_proton->SetOption("COLZ");
  WvsQ2_single_proton->Write();

  W_single_proton->SetXTitle("W (GeV)");
  W_single_proton->Write();

  Q2_single_proton->SetXTitle("Q^{2} (GeV^{2})");
  Q2_single_proton->Write();

  WvsQ2_Ppi0->SetXTitle("W (GeV)");
  WvsQ2_Ppi0->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_Ppi0->SetOption("COLZ");
  WvsQ2_Ppi0->Write();

  W_Ppi0->SetXTitle("W (GeV)");
  W_Ppi0->Write();

  Q2_Ppi0->SetXTitle("Q^{2} (GeV^{2})");
  Q2_Ppi0->Write();

  WvsQ2_MM->SetXTitle("W (GeV)");
  WvsQ2_MM->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_MM->SetOption("COLZ");
  WvsQ2_MM->Write();

  W_MM->SetXTitle("W (GeV)");
  W_MM->Write();

  Q2_MM->SetXTitle("Q^{2} (GeV^{2})");
  Q2_MM->Write();
}

void Histogram::WvsQ2_sec_Write() {
  for (auto &sec : WvsQ2_channel_sec) {
    sec->SetXTitle("W (GeV)");
    sec->SetYTitle("Q^{2} (GeV^{2})");
    sec->SetOption("COLZ");
    sec->Write();
  }
  for (auto &sec : W_channel_sec) {
    sec->SetXTitle("W (GeV)");
    sec->Write();
  }

  if (!_multi) {
    auto canWQ2_channel = std::make_unique<TCanvas>("w_q2_sec_channel_can", "W vs Q2 N #pi^{+} Sec", 1600, 900);
    canWQ2_channel->Divide(2, 3);
    for (short sec = 0; sec < NUM_SECTORS; sec++) {
      canWQ2_channel->cd(sec + 1);
      WvsQ2_channel_sec[sec]->SetXTitle("W (GeV)");
      WvsQ2_channel_sec[sec]->SetYTitle("Q^{2} (GeV^{2})");
      WvsQ2_channel_sec[sec]->SetOption("COLZ");
      WvsQ2_channel_sec[sec]->Draw();
    }
    canWQ2_channel->Write();

    auto canW_channel = std::make_unique<TCanvas>("w_channel_sec_can", "W Sec", 1600, 900);
    canW_channel->Divide(2, 3);
    for (short sec = 0; sec < NUM_SECTORS; sec++) {
      canW_channel->cd(sec + 1);
      W_channel_sec[sec]->SetXTitle("W (GeV)");
      W_channel_sec[sec]->Draw();
    }
    canW_channel->Write();
  }

  for (auto &sec : WvsQ2_sec) {
    sec->SetXTitle("W (GeV)");
    sec->SetYTitle("Q^{2} (GeV^{2})");
    sec->SetOption("COLZ");
    sec->Write();
  }

  for (auto &sec : W_sec) {
    sec->SetXTitle("W (GeV)");
    sec->Write();
  }

  if (!_multi) {
    auto canWQ2 = std::make_unique<TCanvas>("w_q2_sec_can", "W vs Q2 Sec", 1600, 900);
    canWQ2->Divide(2, 3);
    for (short sec = 0; sec < NUM_SECTORS; sec++) {
      canWQ2->cd(sec + 1);
      WvsQ2_sec[sec]->SetXTitle("W (GeV)");
      WvsQ2_sec[sec]->SetYTitle("Q^{2} (GeV^{2})");
      WvsQ2_sec[sec]->SetOption("COLZ");
      WvsQ2_sec[sec]->Draw();
    }
    canWQ2->Write();

    auto canW = std::make_unique<TCanvas>("w_sec_can", "W Sec", 1600, 900);
    canW->Divide(2, 3);
    for (short sec = 0; sec < NUM_SECTORS; sec++) {
      canW->cd(sec + 1);
      W_sec[sec]->SetXTitle("W (GeV)");
      W_sec[sec]->Draw();
    }
    canW->Write();
  }
}

void Histogram::WvsQ2_binned_Write() {
  WvsQ2_binned->SetXTitle("W (GeV)");
  WvsQ2_binned->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_binned->SetOption("COLZ");
  WvsQ2_binned->Write();

  for (auto &sec : W_binned) {
    sec->SetXTitle("W (GeV)");
    sec->Write();
  }
  for (auto &sec : Q2_binned) {
    sec->SetXTitle("Q^{2} (GeV^{2})");
    sec->Write();
  }
}
// W and Q^2

// P and E
void Histogram::MomVsBeta_Fill_pos(float P, float Beta) { MomVsBeta_hist_pos->Fill(P, Beta); }

void Histogram::MomVsBeta_Fill_neg(float P, float Beta) { MomVsBeta_hist_neg->Fill(P, Beta); }

void Histogram::MomVsBeta_Fill_neutral(float P, float Beta) { MomVsBeta_hist_neutral->Fill(P, Beta); }

void Histogram::Fill_proton_ID_P(float p, float beta) { MomVsBeta_proton_ID->Fill(p, beta); }

void Histogram::Fill_Pi_ID_P(float p, float beta) { MomVsBeta_Pi_ID->Fill(p, beta); }

void Histogram::Fill_proton_Pi_ID_P(float p, float beta) { MomVsBeta_proton_Pi_ID->Fill(p, beta); }

void Histogram::MomVsBeta_Fill(float P, float Beta) {
  MomVsBeta_hist->Fill(P, Beta);
  Mom->Fill(P);
}

void Histogram::Photon_flux_Fill(float photon_flux) { photon_flux_hist->Fill(photon_flux); }

void Histogram::MomVsBeta_Write() {
  MomVsBeta_hist->SetXTitle("Momentum (GeV)");
  MomVsBeta_hist->SetYTitle("#beta");
  MomVsBeta_hist->SetOption("COLZ");
  MomVsBeta_hist_pos->SetXTitle("Momentum (GeV)");
  MomVsBeta_hist_pos->SetYTitle("#beta");
  MomVsBeta_hist_pos->SetOption("COLZ");
  MomVsBeta_hist_neg->SetXTitle("Momentum (GeV)");
  MomVsBeta_hist_neg->SetYTitle("#beta");
  MomVsBeta_hist_neg->SetOption("COLZ");

  MomVsBeta_hist_neutral->SetXTitle("Momentum (GeV)");
  MomVsBeta_hist_neutral->SetYTitle("#beta");
  MomVsBeta_hist_neutral->SetOption("COLZ");
  Mom->SetXTitle("Momentum (GeV)");

  MomVsBeta_proton_ID->SetXTitle("Momentum (GeV)");
  MomVsBeta_proton_ID->SetYTitle("#beta");
  MomVsBeta_proton_ID->SetOption("COLZ");
  MomVsBeta_proton_ID->Write();

  MomVsBeta_Pi_ID->SetXTitle("Momentum (GeV)");
  MomVsBeta_Pi_ID->SetYTitle("#beta");
  MomVsBeta_Pi_ID->SetOption("COLZ");
  MomVsBeta_Pi_ID->Write();

  MomVsBeta_proton_Pi_ID->SetXTitle("Momentum (GeV)");
  MomVsBeta_proton_Pi_ID->SetYTitle("#beta");
  MomVsBeta_proton_Pi_ID->SetOption("COLZ");
  MomVsBeta_proton_Pi_ID->Write();

  Energy_hist->Write();
  MomVsBeta_hist->Write();
  MomVsBeta_hist_pos->Write();
  MomVsBeta_hist_neg->Write();
  MomVsBeta_hist_neutral->Write();
  Mom->Write();

  photon_flux_hist->Write();
}

// Missing Mass
void Histogram::Fill_Missing_Mass(float miss_mass) { Missing_Mass->Fill(miss_mass); }

void Histogram::Fill_Missing_Mass(float mm, float mm2) {
  Missing_Mass->Fill(mm);
  Missing_Mass_square->Fill(mm2);
}
void Histogram::Fill_Missing_Mass_strict(float mm, float mm2) {
  Missing_Mass_strict->Fill(mm);
  Missing_Mass_square_strict->Fill(mm2);
}

void Histogram::Fill_Missing_Mass_pi0(float mm, float mm2) {
  if (mm != 0) {
    Missing_Mass_pi0->Fill(mm);
    Missing_Mass_square_pi0->Fill(mm2);
  }
}

void Histogram::Fill_Missing_Mass_twoPi(float mm, float mm2) {
  Missing_Mass_2pi->Fill(mm);
  Missing_Mass_square_2pi->Fill(mm2);
}

void Histogram::Fill_W_Missing_Mass(float W, float mm, float mm2) {
  for (int x = 0; x < W_BINS; x++) {
    if (w_binned_min + (W_width * x) <= W && w_binned_min + (W_width * (x + 1)) >= W) {
      Missing_Mass_WBinned[x]->Fill(mm);
      Missing_Mass_WBinned_square[x]->Fill(mm2);
      continue;
    }
  }
}

// void Histogram::Fill_Mass(float mass) { Mass->Fill(mass); }

void Histogram::Fill_Missing_Mass_square(float miss_mass_2) { Missing_Mass_square->Fill(miss_mass_2); }

void Histogram::Write_Missing_Mass() {
  Missing_Mass->SetXTitle("Mass (GeV)");
  Missing_Mass->Write();

  Missing_Mass_square->SetXTitle("Mass^{2} (GeV^{2})");
  Missing_Mass_square->Write();

  Missing_Mass_strict->SetXTitle("Mass (GeV)");
  Missing_Mass_strict->Write();

  Missing_Mass_square_strict->SetXTitle("Mass^{2} (GeV^{2})");
  Missing_Mass_square_strict->Write();

  Missing_Mass_pi0->SetXTitle("Mass (GeV)");
  Missing_Mass_pi0->Write();
  Missing_Mass_square_pi0->SetXTitle("Mass^{2} (GeV^{2})");
  Missing_Mass_square_pi0->Write();

  Missing_Mass_2pi->SetXTitle("Mass (GeV)");
  Missing_Mass_2pi->Write();
  Missing_Mass_square_2pi->SetXTitle("Mass^{2} (GeV^{2})");
  Missing_Mass_square_2pi->Write();

  TDirectory *mm_binned = RootOutputFile->mkdir("MM_binned");
  mm_binned->cd();
  for (int y = 0; y < W_BINS; y++) {
    if (!_multi) {
      Fit_Missing_Mass_WBinned[y] = new Fits();
      Fit_Missing_Mass_WBinned[y]->FitMissMass(Missing_Mass_WBinned[y]);
    }
    Missing_Mass_WBinned[y]->SetXTitle("Mass (GeV)");
    Missing_Mass_WBinned[y]->Write();
  }
  TDirectory *mm2_binned = RootOutputFile->mkdir("MM2_binned");
  mm2_binned->cd();
  for (int y = 0; y < W_BINS; y++) {
    Missing_Mass_WBinned_square[y]->SetXTitle("Mass^{2} (GeV^{2})");
    Missing_Mass_WBinned_square[y]->Write();
  }
}

void Histogram::makeHists_deltat() {
  for (int jj = 0; jj < NUM_POINTS; jj++) {
    sprintf(hname, "delta_t_p_%d", jj);
    sprintf(htitle, "#Deltat P %d", jj);
    delta_t_hist[0][jj] = new TH1D(hname, htitle, BINS, Dt_min, Dt_max);

    sprintf(hname, "delta_t_pip_%d", jj);
    sprintf(htitle, "#Deltat #pi^{+} %d", jj);
    delta_t_hist[1][jj] = new TH1D(hname, htitle, BINS, Dt_min, Dt_max);

    sprintf(hname, "delta_t_electron_%d", jj);
    sprintf(htitle, "#Deltat electron %d", jj);
    delta_t_hist[2][jj] = new TH1D(hname, htitle, BINS, Dt_min, Dt_max);
  }

  for (int jj = 0; jj < NUM_SECTORS; jj++) {
    for (int jjj = 0; jjj < SC_PADDLE_NUM; jjj++) {
      sprintf(hname, "delta_t_p_sec%d_pad%d", jj + 1, jjj + 1);
      sprintf(htitle, "#Deltat P Sector %d Paddle %d", jj + 1, jjj + 1);
      delta_t_sec_pad_hist[0][jj][jjj] = new TH2D(hname, htitle, BINS / 2, p_min, p_max, BINS / 2, Dt_min, Dt_max);

      sprintf(hname, "delta_t_pip_sec%d_pad%d", jj + 1, jjj + 1);
      sprintf(htitle, "#Deltat #pi^{+} Sector %d Paddle %d", jj + 1, jjj + 1);
      delta_t_sec_pad_hist[1][jj][jjj] = new TH2D(hname, htitle, BINS / 2, p_min, p_max, BINS / 2, Dt_min, Dt_max);

      sprintf(hname, "delta_t_electron_sec%d_pad%d", jj + 1, jjj + 1);
      sprintf(htitle, "#Deltat electron Sector %d Paddle %d", jj + 1, jjj + 1);
      delta_t_sec_pad_hist[2][jj][jjj] = new TH2D(hname, htitle, BINS / 2, p_min, p_max, BINS / 2, Dt_min, Dt_max);
    }
  }
}

void Histogram::Fill_deltat_P(float momentum, float delta_t) { delta_t_mass_P->Fill(momentum, delta_t); }

void Histogram::Fill_deltat_P_PID(float momentum, float delta_t) { delta_t_mass_P_PID->Fill(momentum, delta_t); }

void Histogram::Fill_deltat_PIP(float momentum, float delta_t) { delta_t_mass_PIP->Fill(momentum, delta_t); }

void Histogram::Fill_deltat_PIP_PID(float momentum, float delta_t) { delta_t_mass_PIP_PID->Fill(momentum, delta_t); }

void Histogram::Fill_deltat_PIM(float momentum, float delta_t) { delta_t_mass_PIM->Fill(momentum, delta_t); }

void Histogram::Fill_deltat_PIM_PID(float momentum, float delta_t) { delta_t_mass_PIM_PID->Fill(momentum, delta_t); }

void Histogram::Fill_deltat_electron(float momentum, float delta_t) { delta_t_mass_electron->Fill(momentum, delta_t); }

void Histogram::Fill_deltat_electron_PID(float momentum, float delta_t) {
  delta_t_mass_electron_PID->Fill(momentum, delta_t);
}

void Histogram::Fill_deltat_kp(float momentum, float delta_t) { delta_t_mass_kp->Fill(momentum, delta_t); }

void Histogram::Fill_deltat_kp_PID(float momentum, float delta_t) { delta_t_mass_kp_PID->Fill(momentum, delta_t); }

void Histogram::delta_t_slice_fit() {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

  TF1 *peak = new TF1("peak", "gaus", -1, 1);
  // TF1 *peak = new TF1("peak", func::peak, -1, 1, 3);
  peak->SetParNames("constant", "mean", "#sigma");
  // Bin 50 = 0.5GeV, Bin 300 = 3 GeV
  delta_t_mass_P->FitSlicesY(peak, 50, 300, 10, "QRG5");
  TH1D *delta_t_mass_P_const = (TH1D *)gDirectory->Get("delta_t_mass_P_0");
  TH1D *delta_t_mass_P_mean = (TH1D *)gDirectory->Get("delta_t_mass_P_1");
  TH1D *delta_t_mass_P_sigma = (TH1D *)gDirectory->Get("delta_t_mass_P_2");
  // delta_t_mass_P_const->Write();
  // delta_t_mass_P_mean->Write();
  // delta_t_mass_P_sigma->Write();

  float x[500];
  float y_plus[500];
  float y_minus[500];
  int num = 0;
  for (int i = 0; i < 500; i++) {
    if (delta_t_mass_P_mean->GetBinContent(i) != 0) {
      // Get momentum from bin center
      x[num] = (float)delta_t_mass_P_mean->GetBinCenter(i);
      // mean + 3sigma
      y_plus[num] =
          (float)delta_t_mass_P_mean->GetBinContent(i) + N_SIGMA * (float)delta_t_mass_P_sigma->GetBinContent(i);
      // mean - 3simga
      y_minus[num] =
          (float)delta_t_mass_P_mean->GetBinContent(i) - N_SIGMA * (float)delta_t_mass_P_sigma->GetBinContent(i);
      num++;
    }
  }

  TGraph *P = new TGraph(num, x, y_plus);
  TGraph *M = new TGraph(num, x, y_minus);
  P->SetName("Proton_Pos_graph");
  M->SetName("Proton_Neg_graph");
  TF1 *Proton_Pos_fit = new TF1("Proton_Pos_fit", func::dt_fit, 0.1, 3.0, 2);
  TF1 *Proton_Neg_fit = new TF1("Proton_Neg_fit", func::dt_fit, 0.1, 3.0, 2);
  P->Fit(Proton_Pos_fit, "QRG5", "", 0.2, 2);
  M->Fit(Proton_Neg_fit, "QRG5", "", 0.2, 2);
  // P->Write();
  // M->Write();
  // Proton_Pos_fit->Write();
  // Proton_Neg_fit->Write();

  TCanvas *dt_proton_canvas = new TCanvas("dt_proton_canvas", "#dt Proton", 1280, 720);
  dt_proton_canvas->cd();
  delta_t_mass_P->Draw();
  Proton_Pos_fit->Draw("same");
  Proton_Neg_fit->Draw("same");
  P->Draw("*same");
  M->Draw("*same");
  dt_proton_canvas->Write();

  delta_t_mass_PIP->FitSlicesY(peak, 50, 300, 10, "QRG5");
  TH1D *delta_t_mass_PIP_const = (TH1D *)gDirectory->Get("delta_t_mass_PIP_0");
  TH1D *delta_t_mass_PIP_mean = (TH1D *)gDirectory->Get("delta_t_mass_PIP_1");
  TH1D *delta_t_mass_PIP_sigma = (TH1D *)gDirectory->Get("delta_t_mass_PIP_2");
  float x_pip[500];
  float y_plus_pip[500];
  float y_minus_pip[500];
  num = 0;
  for (int i = 0; i < 500; i++) {
    if (delta_t_mass_PIP_mean->GetBinContent(i) != 0) {
      // Get momentum from bin center
      x_pip[num] = (float)delta_t_mass_PIP_mean->GetBinCenter(i);
      // mean + 3sigma
      y_plus_pip[num] =
          (float)delta_t_mass_PIP_mean->GetBinContent(i) + N_SIGMA * (float)delta_t_mass_PIP_sigma->GetBinContent(i);
      // mean - 3simga
      y_minus_pip[num] =
          (float)delta_t_mass_PIP_mean->GetBinContent(i) - N_SIGMA * (float)delta_t_mass_PIP_sigma->GetBinContent(i);
      num++;
    }
  }

  TGraph *P_pip = new TGraph(num, x_pip, y_plus_pip);
  TGraph *M_pip = new TGraph(num, x_pip, y_minus_pip);
  P_pip->SetName("Pip_Pos_graph");
  M_pip->SetName("Pip_Neg_graph");
  TF1 *Pip_Pos_fit = new TF1("Pip_Pos_fit", func::dt_fit, 0.1, 3.0, 2);
  TF1 *Pip_Neg_fit = new TF1("Pip_Neg_fit", func::dt_fit, 0.1, 3.0, 2);
  P_pip->Fit(Pip_Pos_fit, "QRG5", "", 0.1, 1.75);
  M_pip->Fit(Pip_Neg_fit, "QRG5", "", 0.1, 1.75);
  // P_pip->Write();
  // M_pip->Write();
  // Pip_Pos_fit->Write();
  // Pip_Neg_fit->Write();

  TCanvas *dt_Pip_canvas = new TCanvas("dt_Pip_canvas", "#dt #pi^{+}", 1280, 720);
  dt_Pip_canvas->cd();
  delta_t_mass_PIP->Draw();
  Pip_Pos_fit->Draw("same");
  Pip_Neg_fit->Draw("same");
  P_pip->Draw("*same");
  M_pip->Draw("*same");
  dt_Pip_canvas->Write();
}

void Histogram::delta_t_Write() {
  delta_t_mass_P->SetXTitle("Momentum (GeV)");
  delta_t_mass_P->SetYTitle("#Deltat");
  delta_t_mass_P->SetOption("COLZ");
  delta_t_mass_PIP->SetXTitle("Momentum (GeV)");
  delta_t_mass_PIP->SetYTitle("#Deltat");
  delta_t_mass_PIP->SetOption("COLZ");
  delta_t_mass_PIM->SetXTitle("Momentum (GeV)");
  delta_t_mass_PIM->SetYTitle("#Deltat");
  delta_t_mass_PIM->SetOption("COLZ");
  delta_t_mass_P_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_P_PID->SetYTitle("#Deltat");
  delta_t_mass_P_PID->SetOption("COLZ");
  delta_t_mass_PIP_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_PIP_PID->SetYTitle("#Deltat");
  delta_t_mass_PIP_PID->SetOption("COLZ");
  delta_t_mass_PIM_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_PIM_PID->SetYTitle("#Deltat");
  delta_t_mass_PIM_PID->SetOption("COLZ");
  delta_t_mass_electron->SetXTitle("Momentum (GeV)");
  delta_t_mass_electron->SetYTitle("#Deltat");
  delta_t_mass_electron->SetOption("COLZ");
  delta_t_mass_electron_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_electron_PID->SetYTitle("#Deltat");
  delta_t_mass_electron_PID->SetOption("COLZ");
  delta_t_mass_kp->SetXTitle("Momentum (GeV)");
  delta_t_mass_kp->SetYTitle("#Deltat");
  delta_t_mass_kp->SetOption("COLZ");
  delta_t_mass_kp_PID->SetXTitle("Momentum (GeV)");
  delta_t_mass_kp_PID->SetYTitle("#Deltat");
  delta_t_mass_kp_PID->SetOption("COLZ");

  if (!_multi) delta_t_slice_fit();

  delta_t_mass_P->Write();
  delta_t_mass_P_PID->Write();

  delta_t_mass_PIP->Write();
  delta_t_mass_PIP_PID->Write();
  delta_t_mass_PIM->Write();
  delta_t_mass_PIM_PID->Write();
  delta_t_mass_electron->Write();
  delta_t_mass_electron_PID->Write();
  delta_t_mass_kp->Write();
  delta_t_mass_kp_PID->Write();
}

void Histogram::delta_t_Fill(float momentum, int charge, float delta_t_proton, float delta_t_pip,
                             float delta_t_electron) {
  for (int jj = 0; jj < NUM_POINTS; jj++) {
    if (momentum > jj * bin_width && momentum <= (jj + 1) * bin_width) {
      if (charge == 1 && !std::isnan(delta_t_proton) && !std::isnan(delta_t_pip)) {
        delta_t_hist[0][jj]->Fill(delta_t_proton);
        delta_t_hist[1][jj]->Fill(delta_t_pip);
      }
      if (charge == -1) delta_t_hist[2][jj]->Fill(delta_t_electron);
    }
  }
}

void Histogram::delta_t_slices_Write() {
  std::unique_ptr<Fits> delta_t_cut[3][NUM_POINTS];
  float fit_dt_min = -1.0;
  float fit_dt_max = 1.0;
  for (int j = 0; j < 3; j++) {
    for (int jj = 0; jj < NUM_POINTS; jj++) {
      if (j != 2 && !_multi) {
        delta_t_cut[j][jj] = std::make_unique<Fits>();
        delta_t_cut[j][jj]->Set_min(fit_dt_min);
        delta_t_cut[j][jj]->Set_max(fit_dt_max);
        delta_t_cut[j][jj]->FitGaus(delta_t_hist[j][jj]);
      }

      delta_t_hist[j][jj]->SetYTitle("#Deltat");
      delta_t_hist[j][jj]->Write();
    }
  }
}

void Histogram::delta_t_sec_pad(float momentum, int charge, float delta_t_proton, float delta_t_pip,
                                float delta_t_electron, int sc_sector, int sc_paddle) {
  if (std::isinf(delta_t_proton) || std::isnan(delta_t_proton)) return;
  if (std::isinf(delta_t_pip) || std::isnan(delta_t_pip)) return;
  if (std::isinf(delta_t_electron) || std::isnan(delta_t_electron)) return;
  if (sc_sector == 0 || sc_sector > NUM_SECTORS || sc_sector < 0) return;
  if (sc_paddle == 0 || sc_paddle > SC_PADDLE_NUM || sc_paddle < 0) return;

  if (charge == 1) {
    delta_t_sec_pad_hist[0][sc_sector - 1][sc_paddle - 1]->Fill(momentum, delta_t_proton);
    delta_t_sec_pad_hist[1][sc_sector - 1][sc_paddle - 1]->Fill(momentum, delta_t_pip);
  } else if (charge == -1)
    delta_t_sec_pad_hist[2][sc_sector - 1][sc_paddle - 1]->Fill(momentum, delta_t_electron);
}

void Histogram::delta_t_sec_pad_Write() {
  for (int j = 0; j < 3; j++) {
    for (int jj = 0; jj < NUM_SECTORS; jj++) {
      for (int jjj = 0; jjj < SC_PADDLE_NUM; jjj++) {
        delta_t_sec_pad_hist[j][jj][jjj]->SetYTitle("#Deltat");
        delta_t_sec_pad_hist[j][jj][jjj]->SetOption("COLZ");
        delta_t_sec_pad_hist[j][jj][jjj]->Write();
      }
    }
  }
}

void Histogram::delta_T_canvas() {
  TCanvas *can_dt[NUM_SECTORS][3];
  char can_name[50];
  std::string P_PIP_E;
  for (int particle_i = 0; particle_i < 3; particle_i++) {
    for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
      if (particle_i == 0)
        P_PIP_E = "Proton";
      else if (particle_i == 1)
        P_PIP_E = "Pip";
      else if (particle_i == 2)
        P_PIP_E = "Electron";

      sprintf(can_name, "Sector %d %s", sec_i + 1, P_PIP_E.c_str());
      can_dt[sec_i][particle_i] = new TCanvas(can_name, can_name, 1200, 800);
      can_dt[sec_i][particle_i]->Divide(6, 8);
      for (int pad_i = 0; pad_i < SC_PADDLE_NUM; pad_i++) {
        can_dt[sec_i][particle_i]->cd((int)pad_i + 1);
        delta_t_sec_pad_hist[particle_i][sec_i][pad_i]->Draw("sameCOLZ");
      }
      can_dt[sec_i][particle_i]->Write();
    }
  }
}

void Histogram::CC_fill(int cc_sector, int cc_segment, int cc_pmt, int cc_nphe, float theta_cc) {
  if (cc_pmt == -1) cc_pmt = 2;
  /*
  x_cc_sparse[0] = cc_sector;
  x_cc_sparse[1] = cc_segment;
  x_cc_sparse[2] = cc_pmt;
  x_cc_sparse[3] = cc_nphe;
  cc_sparse->Fill(x_cc_sparse);
  */

  if (cc_sector <= NUM_SECTORS && cc_segment <= segment && cc_pmt < PMT) {
    cc_hist[cc_sector - 1][cc_segment - 1][cc_pmt]->Fill(cc_nphe);
    cc_hist_allSeg[cc_sector - 1][cc_pmt]->Fill(cc_nphe);

    Theta_CC->Fill(cc_segment, theta_cc);
    Theta_CC_Sec[cc_sector - 1]->Fill(cc_segment, theta_cc);
    Theta_CC_Sec_cut[cc_sector - 1]->Fill(cc_segment, theta_cc);
  }
}

void Histogram::makeHists_CC() {
  /*
  cc_sparse->GetAxis(0)->SetTitle(" cc_sector ");
  cc_sparse->GetAxis(1)->SetTitle(" cc_segment ");
  cc_sparse->GetAxis(2)->SetTitle(" cc_pmt ");
  cc_sparse->GetAxis(3)->SetTitle(" cc_nphe ");
  */

  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    sprintf(hname, "Theta_CC_sec%d", sec_i + 1);
    sprintf(htitle, "Theta CC sector %d", sec_i + 1);
    Theta_CC_Sec[sec_i] = new TH2D(hname, htitle, 20, 0.0, 20.0, 60, 0.0, 60.0);
    sprintf(hname, "Theta_CC_sec_cut%d", sec_i + 1);
    sprintf(htitle, "Theta CC sector cut %d", sec_i + 1);
    Theta_CC_Sec_cut[sec_i] = new TH2D(hname, htitle, 20, 0.0, 20.0, 60, 0.0, 60.0);
    for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
      if (pmt_i == 0) L_R_C = "both";
      if (pmt_i == 1) L_R_C = "right";
      if (pmt_i == 2) L_R_C = "left";
      sprintf(hname, "CC_sec%d_%s", sec_i + 1, L_R_C.c_str());
      sprintf(htitle, "CC sector %d %s", sec_i + 1, L_R_C.c_str());
      cc_hist_allSeg[sec_i][pmt_i] = new TH1D(hname, htitle, bins_CC, CC_min, CC_max);
      for (int seg_i = 0; seg_i < segment; seg_i++) {
        sprintf(hname, "CC_sec%d_seg%d_%s", sec_i + 1, seg_i + 1, L_R_C.c_str());
        sprintf(htitle, "CC sector %d segment %d %s", sec_i + 1, seg_i + 1, L_R_C.c_str());
        cc_hist[sec_i][seg_i][pmt_i] = new TH1D(hname, htitle, bins_CC, CC_min, CC_max);
      }
    }
  }
}

void Histogram::Theta_CC_Write() {
  Theta_CC->SetXTitle("CC segment");
  Theta_CC->SetYTitle("#theta_CC");
  Theta_CC->SetOption("COLZ");
  Theta_CC->Write();
  if (!_multi) {
    for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
      Theta_CC_Sec[sec_i]->SetXTitle("CC segment");
      Theta_CC_Sec[sec_i]->SetYTitle("#theta_CC");
      Theta_CC_Sec[sec_i]->SetOption("COLZ");
      Theta_CC_Sec[sec_i]->Write();

      Theta_CC_Sec_cut[sec_i]->SetXTitle("CC segment");
      Theta_CC_Sec_cut[sec_i]->SetYTitle("#theta_CC");
      Theta_CC_Sec_cut[sec_i]->SetOption("COLZ");
      Theta_CC_Sec_cut[sec_i]->Write();
    }

    theta_cc_slice_fit();
  }
}

void Histogram::CC_Write() {
  // cc_sparse->Write();
  std::unique_ptr<Fits> cc_fits[NUM_SECTORS][segment][PMT];
  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
      cc_hist_allSeg[sec_i][pmt_i]->SetYTitle("number photoelectrons");
      cc_hist_allSeg[sec_i][pmt_i]->Write();
      for (int seg_i = 0; seg_i < segment; seg_i++) {
        cc_fits[sec_i][seg_i][pmt_i] = std::make_unique<Fits>();
        /*
        cc_fits[sec_i][seg_i][pmt_i]->FitLandauGaus(cc_hist[sec_i][seg_i][pmt_i]);

        cc_fits[sec_i][seg_i][pmt_i]->Set_lineColor(9);
        cc_fits[sec_i][seg_i][pmt_i]->Set_min(0.0);
        cc_fits[sec_i][seg_i][pmt_i]->Set_max(30.0);
        cc_fits[sec_i][seg_i][pmt_i]->FitLandau(cc_hist[sec_i][seg_i][pmt_i]);
        */
        cc_fits[sec_i][seg_i][pmt_i]->Set_lineColor(8);
        // cc_fits[sec_i][seg_i][pmt_i]->Set_min(30.0);
        cc_fits[sec_i][seg_i][pmt_i]->Set_min(0.0);
        cc_fits[sec_i][seg_i][pmt_i]->Set_max(250.0);
        cc_fits[sec_i][seg_i][pmt_i]->FitGaus(cc_hist[sec_i][seg_i][pmt_i]);

        cc_hist[sec_i][seg_i][pmt_i]->SetYTitle("number photoelectrons");
        cc_hist[sec_i][seg_i][pmt_i]->Write();
      }
    }
  }
}

void Histogram::theta_cc_slice_fit() {
  TCanvas *Theta_CC_canvas[NUM_SECTORS];

  TF1 *peak[NUM_SECTORS];
  TH1D *Theta_CC_0[NUM_SECTORS];
  TH1D *Theta_CC_1[NUM_SECTORS];
  TH1D *Theta_CC_2[NUM_SECTORS];

  TGraph *CC_P[NUM_SECTORS];
  TGraph *CC_M[NUM_SECTORS];
  TF1 *Theta_CC_Pos_fit[NUM_SECTORS];
  TF1 *Theta_CC_Neg_fit[NUM_SECTORS];
  char get_name[100];
  char can_name[100];
  char can_title[100];

  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    // TF1 *peak = new TF1("peak", func::landau, 0, 60, 3);
    peak[sec_i] = new TF1("peak", "landau", 10, 60);

    Theta_CC_Sec[sec_i]->FitSlicesY(peak[sec_i], 0, -1, 20, "QRM+");
    sprintf(get_name, "Theta_CC_sec%d_%d", sec_i + 1, 0);
    Theta_CC_0[sec_i] = (TH1D *)gDirectory->Get(get_name);
    sprintf(get_name, "Theta_CC_sec%d_%d", sec_i + 1, 1);
    Theta_CC_1[sec_i] = (TH1D *)gDirectory->Get(get_name);
    sprintf(get_name, "Theta_CC_sec%d_%d", sec_i + 1, 2);
    Theta_CC_2[sec_i] = (TH1D *)gDirectory->Get(get_name);

    float x[20];
    float y_plus[20];
    float y_minus[20];
    int num = 0;
    for (int i = 0; i < 20; i++) {
      if (Theta_CC_1[sec_i]->GetBinContent(i) != 0) {
        // Get momentum from bin center
        x[num] = (float)Theta_CC_1[sec_i]->GetBinCenter(i);
        // mean + 3sigma
        y_plus[num] =
            (float)Theta_CC_1[sec_i]->GetBinContent(i) + (3 * (float)Theta_CC_2[sec_i]->GetBinContent(i));  //(N_SIGMA)
        // mean - 3simga
        y_minus[num] =
            (float)Theta_CC_1[sec_i]->GetBinContent(i) - (3 * (float)Theta_CC_2[sec_i]->GetBinContent(i));  //(N_SIGMA)
        num++;
      }
    }

    CC_P[sec_i] = new TGraph(num, x, y_plus);
    CC_M[sec_i] = new TGraph(num, x, y_minus);
    sprintf(get_name, "Theta_CC_sec%d_pos_fit", sec_i + 1);
    // Theta_CC_Pos_fit[sec_i] = new TF1(get_name, func::pol1, 0, 14, 2);
    Theta_CC_Pos_fit[sec_i] = new TF1(get_name, func::theta_cc_fit, 1, 14, 4);
    Theta_CC_Pos_fit[sec_i]->SetParNames("intercept", "slope", "a", "exp^c*x");
    Theta_CC_Pos_fit[sec_i]->SetParLimits(3, 0, 10);

    sprintf(get_name, "Theta_CC_sec%d_neg_fit", sec_i + 1);
    // Theta_CC_Neg_fit[sec_i] = new TF1(get_name, func::pol1, 0, 14, 2);
    Theta_CC_Neg_fit[sec_i] = new TF1(get_name, func::theta_cc_fit, 1, 14, 4);
    Theta_CC_Neg_fit[sec_i]->SetParNames("intercept", "slope", "a", "exp^c*x");
    Theta_CC_Neg_fit[sec_i]->SetParLimits(3, 0, 10);
    CC_P[sec_i]->Fit(Theta_CC_Pos_fit[sec_i], "QRM+", "", 1, 14);
    CC_M[sec_i]->Fit(Theta_CC_Neg_fit[sec_i], "QRM+", "", 1, 14);

    sprintf(can_name, "Theta_CC_sec%d", sec_i + 1);
    sprintf(can_title, "Theta CC sec_i %d", sec_i + 1);
    Theta_CC_canvas[sec_i] = new TCanvas(can_name, can_title, 1280, 720);
    Theta_CC_canvas[sec_i]->cd();
    Theta_CC_Sec[sec_i]->Draw();
    Theta_CC_Pos_fit[sec_i]->Draw("same");
    Theta_CC_Neg_fit[sec_i]->Draw("same");
    CC_P[sec_i]->Draw("*same");
    CC_M[sec_i]->Draw("*same");
    Theta_CC_canvas[sec_i]->Write();

    delete peak[sec_i];
    delete Theta_CC_0[sec_i];
    delete Theta_CC_1[sec_i];
    delete Theta_CC_2[sec_i];
  }
}

void Histogram::CC_canvas() {
  TCanvas *can[NUM_SECTORS][PMT];
  char can_name[50];
  for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
    for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
      if (pmt_i == 0) L_R_C = "both";
      if (pmt_i == 1) L_R_C = "right";
      if (pmt_i == 2) L_R_C = "left";

      sprintf(can_name, "Sector %d %s", sec_i + 1, L_R_C.c_str());
      can[sec_i][pmt_i] = new TCanvas(can_name, can_name, 1200, 800);
      can[sec_i][pmt_i]->Divide(6, 3);
      for (int seg_i = 0; seg_i < segment; seg_i++) {
        can[sec_i][pmt_i]->cd((int)seg_i + 1);
        // cc_hist[sec_i][seg_i][pmt_i]->Fit("gaus");
        cc_hist[sec_i][seg_i][pmt_i]->Draw("same");
      }
      can[sec_i][pmt_i]->Write();
    }
  }
}

void Histogram::makeHists_fid() {
  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    sprintf(hname, "electron_fid_sec%d", sec_i + 1);
    sprintf(htitle, "electron_fid_sec%d", sec_i + 1);
    electron_fid_sec_hist[sec_i] =
        new TH2D(hname, htitle, BINS, min_phi[sec_i], max_phi[sec_i], BINS, theta_min, theta_max);

    for (int t = 0; t < 3; t++) {
      sprintf(hname, "hadron_fid_sec%d_%d", sec_i + 1, t);
      sprintf(htitle, "hadron_fid_sec%d_%d", sec_i + 1, t);
      hadron_fid_sec_hist[t][sec_i] =
          new TH2D(hname, htitle, BINS, min_phi[sec_i], max_phi[sec_i], BINS, theta_min, theta_max);
    }
  }
}

void Histogram::Fill_electron_fid(float theta, float phi, int sector) {
  electron_fid_hist->Fill(phi, theta);
  if (sector == 0 || sector > NUM_SECTORS) return;
  electron_fid_sec_hist[sector - 1]->Fill(phi, theta);
}

void Histogram::Fill_hadron_fid(float theta, float phi, int sector, int id) {
  if (sector == 0 || sector > NUM_SECTORS) return;
  if (id == PROTON) {
    hadron_fid_hist[0]->Fill(phi, theta);
    hadron_fid_sec_hist[0][sector - 1]->Fill(phi, theta);
    hadron_fid_hist[1]->Fill(phi, theta);
    hadron_fid_sec_hist[1][sector - 1]->Fill(phi, theta);
  } else if (id == PIP) {
    hadron_fid_hist[0]->Fill(phi, theta);
    hadron_fid_sec_hist[0][sector - 1]->Fill(phi, theta);
    hadron_fid_hist[2]->Fill(phi, theta);
    hadron_fid_sec_hist[2][sector - 1]->Fill(phi, theta);
  }
}

void Histogram::Fid_Write() {
  electron_fid_hist->SetYTitle("#theta");
  electron_fid_hist->SetXTitle("#phi");
  electron_fid_hist->SetOption("COLZ");
  electron_fid_hist->Write();

  for (int t = 0; t < 3; t++) {
    hadron_fid_hist[t]->SetYTitle("#theta");
    hadron_fid_hist[t]->SetXTitle("#phi");
    hadron_fid_hist[t]->SetOption("COLZ");
    hadron_fid_hist[t]->Write();
  }
  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    electron_fid_sec_hist[sec_i]->SetXTitle("#phi");
    electron_fid_sec_hist[sec_i]->SetYTitle("#theta");
    electron_fid_sec_hist[sec_i]->SetOption("COLZ");
    electron_fid_sec_hist[sec_i]->Write();
    for (int t = 0; t < 3; t++) {
      hadron_fid_sec_hist[t][sec_i]->SetYTitle("#theta");
      hadron_fid_sec_hist[t][sec_i]->SetXTitle("#phi");
      hadron_fid_sec_hist[t][sec_i]->SetOption("COLZ");
      hadron_fid_sec_hist[t][sec_i]->Write();
    }
  }

  if (!_multi) {
    float slice_width = ((float)100 / (float)FID_SLICES);
    float y_width = (60.0 / (float)100);
    std::unique_ptr<Fits> SliceFit[NUM_SECTORS][FID_SLICES];
    TGraph *fid[NUM_SECTORS];
    std::unique_ptr<Fits> FidGraph[NUM_SECTORS];

    TCanvas *electron_fid_can[NUM_SECTORS];
    float x[FID_SLICES * 2];
    float y[FID_SLICES * 2];

    // std::cout << "sec,y,x" << '\n';
    for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
      sprintf(hname, "electron_fid_sector_%d", sec_i + 1);
      sprintf(htitle, "electron_fid_sector_%d", sec_i + 1);
      electron_fid_can[sec_i] = new TCanvas(hname, htitle, 1280, 720);
      for (int slice = start_slice; slice < FID_SLICES; slice++) {
        sprintf(hname, "electron_fid_sec_%d_%d", sec_i + 1, slice + 1);
        electron_fid_sec_slice[sec_i][slice] = (TH1D *)electron_fid_sec_hist[sec_i]->ProjectionY(
            hname, slice_width * slice, slice_width * slice + (slice_width - 1));
        electron_fid_sec_slice[sec_i][slice]->Rebin(4);
        SliceFit[sec_i][slice] = std::make_unique<Fits>();
        SliceFit[sec_i][slice]->Set_min(min_phi[sec_i]);
        SliceFit[sec_i][slice]->Set_max(max_phi[sec_i]);
        SliceFit[sec_i][slice]->FitGenNormal(electron_fid_sec_slice[sec_i][slice]);

        if (SliceFit[sec_i][slice]->Get_left_edge() == SliceFit[sec_i][slice]->Get_left_edge() &&
            SliceFit[sec_i][slice]->Get_right_edge() == SliceFit[sec_i][slice]->Get_right_edge()) {
          x[slice + FID_SLICES] = y_width * slice_width * slice;
          y[slice + FID_SLICES] = SliceFit[sec_i][slice]->Get_left_edge();
          x[slice] = y_width * slice_width * slice;
          y[slice] = SliceFit[sec_i][slice]->Get_right_edge();

          // std::cout << sec_i + 1 << "," << y[slice] << "," << x[slice + FID_SLICES] << '\n';
          // std::cout << sec_i + 1 << "," << y[slice] << "," << x[slice] << '\n';
        }
      }

      fid[sec_i] = new TGraph(FID_SLICES * 2, x, y);
      FidGraph[sec_i] = std::make_unique<Fits>();

      FidGraph[sec_i]->Set_min(min_phi[sec_i]);
      FidGraph[sec_i]->Set_max(max_phi[sec_i]);
      // FidGraph[sec_i]->FitFiducial(fid[sec_i], sec_i);
      // FidGraph[sec_i]->FitPoly_fid(fid[sec_i]);

      electron_fid_can[sec_i]->cd();

      electron_fid_sec_hist[sec_i]->Draw();
      fid[sec_i]->Draw("*same");
      electron_fid_can[sec_i]->Write();
    }
  }
}

void Histogram::fid_canvas() {
  TCanvas *can[NUM_SECTORS];
  char can_name[50];

  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    sprintf(can_name, "Electron Fid Sector %d Slices", sec_i + 1);
    can[sec_i] = new TCanvas(can_name, can_name, 1600, 900);
    can[sec_i]->Divide(10, (FID_SLICES - 10) / 10);
    int x = 1;
    for (int slice = start_slice; slice < FID_SLICES; slice++) {
      can[sec_i]->cd(x++);
      electron_fid_sec_slice[sec_i][slice]->Draw("same");
    }
    can[sec_i]->Write();
  }
}

void Histogram::makeHists_EC() {
  for (int n = 0; n < NUM_POINTS; n++) {
    sprintf(hname, "ec_%d", n);
    sprintf(htitle, "Sampling Fraction %d", n);
    EC_hist[n] = new TH1D(hname, htitle, BINS, EC_min, EC_max);

    sprintf(hname, "ec_cut_%d", n);
    sprintf(htitle, "Sampling Fraction cut %d", n);
    EC_hist_cut[n] = new TH1D(hname, htitle, BINS, EC_min, EC_max);
  }
}

void Histogram::EC_fill(float etot, float momentum) {
  float sampling_frac = etot / momentum;
  EC_sampling_fraction->Fill(momentum, sampling_frac);

  for (int n = 0; n < NUM_POINTS; n++) {
    if (momentum > n * bin_width && momentum <= (n + 1) * bin_width) {
      EC_hist[n]->Fill(sampling_frac);
    }
  }
}

void Histogram::EC_inout(float Ein, float Eout) { ECin_ECout->Fill(Ein, Eout); }

void Histogram::TM_Fill(float momentum, float theta) { Theta_vs_mom->Fill(momentum, theta); }

void Histogram::EC_cut_fill(float etot, float momentum) {
  float sampling_frac = etot / momentum;
  EC_sampling_fraction_cut->Fill(momentum, sampling_frac);

  for (int n = 0; n < NUM_POINTS; n++) {
    if (momentum > n * bin_width && momentum <= (n + 1) * bin_width) {
      EC_hist_cut[n]->Fill(sampling_frac);
    }
  }
}

void Histogram::EC_slice_fit() {
  if (EC_sampling_fraction->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  // Header *fit_functions = new Header("EC_fit_functions.hpp", "FF");

  TF1 *peak = new TF1("peak", "gaus", 0.2, 0.4);
  EC_sampling_fraction->FitSlicesY(peak, 0, -1, 0, "QRG5");
  TH1D *EC_sampling_fraction_0 = (TH1D *)gDirectory->Get("EC_sampling_fraction_0");
  TH1D *EC_sampling_fraction_1 = (TH1D *)gDirectory->Get("EC_sampling_fraction_1");
  TH1D *EC_sampling_fraction_2 = (TH1D *)gDirectory->Get("EC_sampling_fraction_2");
  float x[BINS];
  float y_plus[BINS];
  float y_minus[BINS];
  int num = 0;
  for (int i = 0; i < BINS; i++) {
    if (EC_sampling_fraction_1->GetBinContent(i) != 0) {
      // Get momentum from bin center
      x[num] = (float)EC_sampling_fraction_1->GetBinCenter(i);
      // mean + 3sigma
      y_plus[num] =
          (float)EC_sampling_fraction_1->GetBinContent(i) + 2.0 * (float)EC_sampling_fraction_2->GetBinContent(i);
      // mean - 3simga
      y_minus[num] =
          (float)EC_sampling_fraction_1->GetBinContent(i) - 2.0 * (float)EC_sampling_fraction_2->GetBinContent(i);
      num++;
    }
  }

  TGraph *EC_P = new TGraph(num, x, y_plus);
  TGraph *EC_M = new TGraph(num, x, y_minus);
  EC_P->SetName("Positive_EC_graph");
  EC_M->SetName("Negative_EC_graph");

  TF1 *EC_P_fit = new TF1("EC_P_fit", func::ec_fit_func, 0.25, 4.0, 3);
  TF1 *EC_M_fit = new TF1("EC_M_fit", func::ec_fit_func, 0.25, 4.0, 3);

  /*
    EC_P_fit->SetParameters(0.3296, 0.002571, 4.8e-7);
    EC_M_fit->SetParameters(0.1715, 0.02044, -1.581e-5);
    EC_P_fit->SetParLimits(2, 1.0e-7, 5.0e-7);
    EC_M_fit->SetParLimits(2, -2.0e-5, -1.0e-5);
  */

  EC_P->Fit(EC_P_fit, "QMRG+", "", 0.75, 3.75);
  EC_M->Fit(EC_M_fit, "QMRG+", "", 0.75, 3.75);

  EC_P_fit->Write();
  EC_M_fit->Write();
  EC_P->Write();
  EC_M->Write();

  TCanvas *EC_canvas = new TCanvas("EC_canvas", "EC canvas", 1280, 720);
  EC_canvas->cd();
  EC_sampling_fraction->Draw();
  EC_P_fit->Draw("same");
  EC_M_fit->Draw("same");
  EC_P->Draw("*same");
  EC_M->Draw("*same");
  EC_canvas->Write();
  /*
    fit_functions->NewFunction();
    fit_functions->Set_RetrunType("float");
    fit_functions->Set_FuncName("EC_P_fit");
    fit_functions->Set_FuncInputs("float x");
    fit_functions->Set_Function(EC_P_fit->GetExpFormula("P"));
    fit_functions->WriteFunction();

    fit_functions->NewFunction();
    fit_functions->Set_RetrunType("float");
    fit_functions->Set_FuncName("EC_M_fit");
    fit_functions->Set_FuncInputs("float x");
    fit_functions->Set_Function(EC_M_fit->GetExpFormula("P"));
    fit_functions->WriteFunction();

    delete fit_functions;
  */
}

void Histogram::EC_slices_Write() {
  std::unique_ptr<Fits> EC_fit[NUM_POINTS];
  float fit_ec_min = 0.2;
  float fit_ec_max = 0.4;
  for (int n = 0; n < NUM_POINTS; n++) {
    EC_fit[n] = std::make_unique<Fits>();
    EC_fit[n]->Set_min(fit_ec_min);
    EC_fit[n]->Set_max(fit_ec_max);
    EC_fit[n]->FitGaus(EC_hist[n]);
    EC_hist[n]->SetYTitle("Sampling Fraction");
    EC_hist[n]->Write();
    EC_hist_cut[n]->SetYTitle("Sampling Fraction");
    EC_hist_cut[n]->Write();
  }
}

void Histogram::EC_Write() {
  EC_sampling_fraction->SetXTitle("Momentum (GeV)");
  EC_sampling_fraction->SetYTitle("Sampling Fraction");
  EC_sampling_fraction->SetOption("COLZ");
  EC_sampling_fraction->Write();

  EC_sampling_fraction_cut->SetXTitle("Momentum (GeV)");
  EC_sampling_fraction_cut->SetYTitle("Sampling Fraction");
  EC_sampling_fraction_cut->SetOption("COLZ");
  EC_sampling_fraction_cut->Write();
  if (!_multi) EC_slice_fit();

  Theta_vs_mom->SetXTitle("Momentum (GeV)");
  Theta_vs_mom->SetYTitle("Theta #theta");
  Theta_vs_mom->SetOption("COLZ");
  Theta_vs_mom->Write();

  ECin_ECout->SetXTitle("EC_{inner}");
  ECin_ECout->SetYTitle("EC_{outer}");
  ECin_ECout->SetOption("COLZ");
  ECin_ECout->Write();
}

void Histogram::Fill_Beam_Position(float vertex_x, float vertex_y, float vertex_z) {
  Beam_Position->Fill(vertex_x, vertex_y);
  Beam_Position_X->Fill(vertex_x);
  Beam_Position_Y->Fill(vertex_y);
  Beam_Position_Z->Fill(vertex_z);

  // Phi vs vertex
}

void Histogram::Beam_Position_Write() {
  Beam_Position->SetXTitle("X");
  Beam_Position->SetYTitle("Y");
  Beam_Position->SetOption("COLZ");
  Beam_Position->Write();

  Beam_Position_X->SetXTitle("X");
  Beam_Position_X->Write();

  Beam_Position_Y->SetXTitle("Y");
  Beam_Position_Y->Write();

  Beam_Position_Z->SetXTitle("Z");
  Beam_Position_Z->Write();
}

void Histogram::Fill_Target_Vertex(float vertex_x, float vertex_y, float vertex_z) {
  if (0 == vertex_x) return;
  if (0 == vertex_y && 0 == vertex_z) return;
  target_vertex_X->Fill(vertex_x);
  target_vertex_Y->Fill(vertex_y);
  target_vertex_Z->Fill(vertex_z);
  target_vertex_xy->Fill(vertex_x, vertex_y);
  target_vertex_zy->Fill(vertex_z, vertex_y);
  target_vertex_zx->Fill(vertex_z, vertex_x);
}

void Histogram::Target_Vertex_Write() {
  target_vertex_X->SetXTitle("X");
  target_vertex_X->Write();

  target_vertex_Y->SetXTitle("Y");
  target_vertex_Y->Write();

  target_vertex_Z->SetXTitle("Z");
  target_vertex_Z->Write();

  target_vertex_xy->SetXTitle("X");
  target_vertex_xy->SetYTitle("Y");
  target_vertex_xy->SetOption("COLZ");
  target_vertex_xy->Write();

  target_vertex_zy->SetXTitle("Z");
  target_vertex_zy->SetYTitle("Y");
  target_vertex_zy->SetOption("COLZ");
  target_vertex_zy->Write();

  target_vertex_zx->SetXTitle("Z");
  target_vertex_zx->SetYTitle("X");
  target_vertex_zx->SetOption("COLZ");
  target_vertex_zx->Write();
}

void Histogram::Fill_E_Prime(TLorentzVector e_prime) {
  if (e_prime.E() > 0.1) energy_no_cuts->Fill(e_prime.E());
}
void Histogram::Fill_E_Prime_fid(TLorentzVector e_prime) {
  if (e_prime.E() > 0.1) energy_fid_cuts->Fill(e_prime.E());
}
void Histogram::Fill_E_Prime_channel(TLorentzVector e_prime) {
  if (e_prime.E() > 0.1) energy_channel_cuts->Fill(e_prime.E());
}

void Histogram::E_Prime_Write() {
  energy_no_cuts->SetXTitle("Energy (GeV)");
  energy_no_cuts->Write();

  energy_fid_cuts->SetXTitle("Energy (GeV)");
  energy_fid_cuts->Write();

  energy_channel_cuts->SetXTitle("Energy (GeV)");
  energy_channel_cuts->Write();
}

mcHistogram::~mcHistogram() {}

void mcHistogram::Write() {
  std::cerr << GREEN << "\nWriting" << DEF << std::endl;
  Histogram::Write();
  RootOutputFile->cd();
  // Start of cuts
  // Fits *MM_neutron_cut = new Fits();
  // MM_neutron_cut->Set_min(0.8);
  // MM_neutron_cut->Set_max(1.2);
  // MM_neutron_cut->FitBreitWigner(Missing_Mass.get());
  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  // TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2 MC");
  // WvsQ2_folder->cd();
  WvsQ2_MC_Write();

  TDirectory *delta_mom = RootOutputFile->mkdir("delta_mom");
  delta_mom->cd();
  Write_DeltaP();
  std::cerr << BOLDBLUE << "Done!!!" << DEF << std::endl;
}

void mcHistogram::makeMCHists() {
  std::string xyz[4] = {"X", "Y", "Z", "all"};
  for (int i = 0; i < 4; i++) {
    sprintf(hname, "dPvsP_%s", xyz[i].c_str());
    sprintf(htitle, "#DeltaP/P_{rec} vs P_{%s}", xyz[i].c_str());
    delta_p[i] = std::make_unique<TH1D>(hname, htitle, 500, -0.5, 0.5);
  }

  for (int y = 0; y < Q2_BINS; y++) {
    sprintf(hname, "W_MC_%0.3f_%0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    sprintf(htitle, "W hist from true MC\nQ^{2} %0.3f %0.3f", q2_binned_min + (Q2_width * y),
            q2_binned_min + (Q2_width * (y + 1)));
    W_binned_MC[y] = std::make_unique<TH1D>(hname, htitle, BINS, w_binned_min, w_binned_max);
  }
}

void mcHistogram::Fill_WQ2_MC(double W, double Q2) {
  WvsQ2_MC->Fill(W, Q2);
  W_MC->Fill(W);
  WvsQ2_binned_MC->Fill(W, Q2);
  for (int y = 0; y < Q2_BINS; y++) {
    if (q2_binned_min + (Q2_width * y) <= Q2 && q2_binned_min + (Q2_width * (y + 1)) >= Q2) {
      W_binned_MC[y]->Fill(W);
      continue;
    }
  }
}

void mcHistogram::Fill_P(const std::shared_ptr<Branches> &d) {
  double P = 0;
  for (int part_num = 0; part_num < d->gpart(); part_num++) {
    double px = d->p(part_num) * d->cx(part_num);
    delta_p[0]->Fill((px - d->pxpart(part_num)) / px);
    double py = d->p(part_num) * d->cy(part_num);
    delta_p[1]->Fill((py - d->pypart(part_num)) / py);
    double pz = d->p(part_num) * d->cz(part_num);
    delta_p[2]->Fill((pz - d->pzpart(part_num)) / pz);
    P = TMath::Sqrt(d->pxpart(part_num) * d->pxpart(part_num) + d->pypart(part_num) * d->pypart(part_num) +
                    d->pzpart(part_num) * d->pzpart(part_num));
    delta_p[3]->Fill((d->p(part_num) - P) / d->p(part_num));
  }
}

void mcHistogram::WvsQ2_MC_Write() {
  WvsQ2_MC->SetXTitle("W (GeV)");
  WvsQ2_MC->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_MC->SetOption("COLZ");
  WvsQ2_MC->Write();

  W_MC->SetXTitle("W (GeV)");
  W_MC->Write();
}

void mcHistogram::Write_DeltaP() {
  TCanvas *dp_canvas = new TCanvas("dp_canvas", "#Delta P", 1280, 720);
  dp_canvas->Divide(2, 2);
  for (int i = 0; i < 4; i++) {
    dp_canvas->cd(i + 1);
    delta_p[i]->SetXTitle("#Delta P (GeV)");
    delta_p[i]->Fit("gaus", "QM+", "", -0.1, 0.1);
    delta_p[i]->Draw("same");
  }
  dp_canvas->Write();
  for (auto &dp : delta_p) {
    dp->SetXTitle("#Delta P (GeV)");
    dp->Write();
  }
}
