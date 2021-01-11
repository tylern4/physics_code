/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram() {
  def = std::make_shared<TCanvas>();
  makeHists_WvsQ2();
  makeHists_deltat();
  makeHists_EC();
  makeHists_CC();
  makeHists_fid();

  ndhist_nPip = std::make_unique<THnD>("ndhist", "ndhist", DIMENSIONS, nbins, xmin, xmax);
  ndhist_nPip->GetAxis(0)->SetName("W");
  ndhist_nPip->GetAxis(1)->SetName("Q2");
  ndhist_nPip->GetAxis(2)->SetName("cos_Theta_star");
  ndhist_nPip->GetAxis(3)->SetName("Phi_star");

  ndhist_protPi0 = std::make_unique<THnD>("ndhist_protPi0", "ndhist_protPi0", DIMENSIONS, nbins, xmin, xmax);
  ndhist_protPi0->GetAxis(0)->SetName("W");
  ndhist_protPi0->GetAxis(1)->SetName("Q2");
  ndhist_protPi0->GetAxis(2)->SetName("cos_Theta_star");
  ndhist_protPi0->GetAxis(3)->SetName("Phi_star");
}

Histogram::Histogram(const std::string &output_file) : Histogram() {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
}

Histogram::~Histogram() {}

void Histogram::FillEvent(const std::shared_ptr<Reaction> &event) {
  this->Fill_ND(event);
  this->Fill_Mass_photons(event);
  ThetaVsPhi_hist->Fill(event->Theta_star(), event->Phi_star());
  CosThetaVsPhi_hist->Fill(cosf(event->Theta_star()), event->Phi_star());

  if (event->SinglePip()) {
    this->Fill_Missing_Mass(event);
    this->Fill_W_Missing_Mass(event->W(), event->MM(), event->MM2());
  }

  if (event->MM_cut()) this->Fill_MM_WQ2(event->W(), event->Q2());

  if (event->channel()) {
    this->Fill_channel_WQ2(event);
    this->Fill_Missing_Mass_strict(event->MM(), event->MM2());
    this->Fill_E_Prime_channel(event->e_mu_prime());
    ThetaVsPhi_channel->Fill(event->Theta_star(), event->Phi_star());
    CosThetaVsPhi_channel->Fill(cosf(event->Theta_star()), event->Phi_star());
  }

  if ((event->SinglePip() || event->NeutronPip()))
    this->Fill_NeutronPip_WQ2(event->W(), event->Q2(), event->MM(), event->MM2());

  if (event->SingleP()) {
    this->Fill_Missing_Mass_pi0(event->MM(), event->MM2());
    this->Fill_elastic(event);
    if (event->MM() >= 0.05 && event->MM() <= 0.3) {
      Mass_pi0_otherCut->Fill(event->pi0_mass());
      Mass_square_pi0_otherCut->Fill(event->pi0_mass2());
    }
    if (event->pi0_mass() >= 0.05 && event->pi0_mass() <= 0.2) {
      Missing_Mass_pi0_otherCut->Fill(event->MM());
      Missing_Mass_square_pi0_otherCut->Fill(event->MM2());
    }
  }

  if (event->PPi0()) this->Fill_P_PI0(event->W(), event->Q2());
  if (event->TwoPion()) this->Fill_Missing_Mass_twoPi(event->MM(), event->MM2());
}

void Histogram::Fill_ND(const std::shared_ptr<Reaction> &event) {
  std::lock_guard<std::mutex> lk(mutex);
  bool _good = true;
  _good &= !std::isnan(event->W());
  if (!_good) return;
  _good &= !std::isnan(event->Q2());
  if (!_good) return;
  _good &= !std::isnan(event->Theta_star());
  if (!_good) return;
  _good &= !std::isnan(event->Phi_star());
  if (!_good) return;

  std::array<double, DIMENSIONS> to_fill = {event->W(), event->Q2(), cos(event->Theta_star()), event->Phi_star()};

  if (_good && event->channel()) {
    ndhist_nPip->Fill(to_fill.data());
  } else if (_good && event->PPi0()) {
    ndhist_protPi0->Fill(to_fill.data());
  }
}

void Histogram::Write(const std::string &output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  Write();
}

void Histogram::Write() {
  ndhist_nPip->Write();
  ndhist_protPi0->Write();
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
  /*
    std::cerr << "MM:mean,+3,-3" << std::endl;
    std::cerr << MM_neutron_cut->Get_mean() << ",";
    std::cerr << MM_neutron_cut->Get_mean() + 3 * MM_neutron_cut->Get_sigma() << ",";
    std::cerr << MM_neutron_cut->Get_mean() - 3 * MM_neutron_cut->Get_sigma() << std::endl;
  */
  std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

// W and Q^2
void Histogram::makeHists_WvsQ2() {
  WvsQ2_sec.resize(NUM_SECTORS);
  W_sec.resize(NUM_SECTORS);
  WvsQ2_channel_sec.resize(NUM_SECTORS);
  W_channel_sec.resize(NUM_SECTORS);
  Missing_Mass_small_sec.resize(NUM_SECTORS);
  Missing_Mass_Sq_small_sec.resize(NUM_SECTORS);

  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    WvsQ2_sec[sec] = std::make_shared<TH2D>(Form("W_vs_Q2_sec_%d", sec + 1), Form("W vs Q^{2} Sector: %d", sec + 1),
                                            BINS, w_min, w_max, BINS, q2_min, q2_max);

    W_sec[sec] = std::make_shared<TH1D>(Form("W_sec_%d", sec + 1), Form("W Sector: %d", sec + 1), BINS, w_min, w_max);

    WvsQ2_channel_sec[sec] = std::make_shared<TH2D>(Form("W_vs_Q2_channel_sec_%d", sec + 1),
                                                    Form("W vs Q^{2} N #pi^{+} Sector: %d", sec + 1), BINS, w_min,
                                                    w_max, BINS / 2, q2_min, q2_max);

    W_channel_sec[sec] = std::make_shared<TH1D>(Form("W_channel_sec_%d", sec + 1),
                                                Form("W N #pi^{+} Sector: %d", sec + 1), BINS / 2, w_min, w_max);
    Missing_Mass_small_sec[sec] = std::make_shared<TH1D>(Form("Missing_Mass_small_%d", sec),
                                                         Form("e(p,#pi^{+} X)e' Sector %d", sec), BINS_MM, 0.8, 1.3);
    Missing_Mass_Sq_small_sec[sec] = std::make_shared<TH1D>(Form("Missing_Mass_Sq_small_%d", sec),
                                                            Form("e(p,#pi^{+} X)e' Sector %d", sec), BINS_MM, 0.7, 1.5);
  }
  W_binned.resize(Q2_BINS);
  for (short y = 0; y < Q2_BINS; y++) {
    float _min = q2_binned_min + (Q2_width * y);
    float _max = q2_binned_min + (Q2_width * (y + 1));
    W_binned[y] =
        std::make_shared<TH1D>(Form("W_%0.3f_%0.3f", _min, _max), Form("W hist\nQ^{2} %0.3f %0.3f", _min, _max), BINS,
                               w_binned_min, w_binned_max);
  }

  Q2_binned.resize(W_BINS);
  Missing_Mass_WBinned.resize(W_BINS);
  Missing_Mass_WBinned_square.resize(W_BINS);
  for (short x = 0; x < W_BINS; x++) {
    float _min = w_binned_min + (W_width * x);
    float _max = w_binned_min + (W_width * (x + 1));

    Q2_binned[x] =
        std::make_shared<TH1D>(Form("Q2_%0.3f_%0.3f", _min, _max), Form("Q^{2} hist\nW %0.3f %0.3f", _min, _max), BINS,
                               q2_binned_min, q2_binned_max);

    Missing_Mass_WBinned[x] = std::make_shared<TH1D>(Form("MM_W_%0.3f_%0.3f", _min, _max),
                                                     Form("Missing Mass\nW %0.3f %0.3f", _min, _max), BINS, 0.8, 1.5);

    Missing_Mass_WBinned_square[x] =
        std::make_shared<TH1D>(Form("MM2_W_%0.3f_%0.3f", _min, _max),
                               Form("Missing Mass^{2}\nW %0.3f %0.3f", _min, _max), BINS, 0.8 * 0.8, 1.5 * 1.5);
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

void Histogram::Fill_channel_WQ2(const std::shared_ptr<Reaction> &event) {
  float W = event->W();
  float Q2 = event->Q2();
  short sector = event->sector();
  E_prime_hist->Fill(event->E_prime());
  Q2_vs_xb->Fill(event->xb(), Q2);

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

void Histogram::Fill_elastic(const std::shared_ptr<Reaction> &event) {
  if (event->W() > 2.0) return;
  elastic_phi->Fill(event->phi_diff());
  if (event->phi_diff() > (3.125) && event->phi_diff() < (3.165)) {
    elastic_MM->Fill(event->MM2());
    if (abs(event->MM2()) < 0.005) elastic_thetaVsP->Fill(event->P_Mom(), event->P_Theta());
  }

  if (event->elastic()) {
    WvsQ2_elastic->Fill(event->W(), event->Q2());
    W_elastic->Fill(event->W());
    Q2_elastic->Fill(event->Q2());
  }
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
  WvsQ2_channel->SetXTitle("W (GeV/c^{2})");
  WvsQ2_channel->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_channel->SetOption("COLZ");
  WvsQ2_channel->Write();

  ThetaVsPhi_channel->SetXTitle("#Theta");
  ThetaVsPhi_channel->SetYTitle("#phi");
  ThetaVsPhi_channel->SetOption("COLZ");
  ThetaVsPhi_channel->Write();

  CosThetaVsPhi_channel->SetXTitle("Cos(#Theta)");
  CosThetaVsPhi_channel->SetYTitle("#phi");
  CosThetaVsPhi_channel->SetOption("COLZ");
  CosThetaVsPhi_channel->Write();

  W_channel->SetXTitle("W (GeV/c^{2})");
  W_channel->Write();

  Q2_channel->SetXTitle("Q^{2} (GeV^{2})");
  Q2_channel->Write();

  E_prime_hist->SetXTitle("Energy (GeV)");
  E_prime_hist->Write();

  Q2_vs_xb->SetXTitle("x_{b}");
  Q2_vs_xb->SetYTitle("Q^{2}");
  Q2_vs_xb->SetOption("COLZ");
  Q2_vs_xb->Write();

  WvsQ2_hist->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
  WvsQ2_hist->SetYTitle("Q^{2} [Momentum Transfer] (GeV^{2})");
  WvsQ2_hist->SetOption("COLZ");
  WvsQ2_hist->Write();

  ThetaVsPhi_hist->SetXTitle("#Theta");
  ThetaVsPhi_hist->SetYTitle("#phi");
  ThetaVsPhi_hist->SetOption("COLZ");
  ThetaVsPhi_hist->Write();

  CosThetaVsPhi_hist->SetXTitle("Cos(#Theta)");
  CosThetaVsPhi_hist->SetYTitle("#phi");
  CosThetaVsPhi_hist->SetOption("COLZ");
  CosThetaVsPhi_hist->Write();

  W_hist->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
  W_hist->Write();

  Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
  Q2_hist->Write();

  WvsQ2_proton->SetXTitle("W (GeV/c^{2})");
  WvsQ2_proton->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_proton->SetOption("COLZ");
  WvsQ2_proton->Write();

  W_proton->SetXTitle("W (GeV/c^{2})");
  W_proton->Write();

  Q2_proton->SetXTitle("Q^{2} (GeV^{2})");
  Q2_proton->Write();

  WvsQ2_pion->SetXTitle("W (GeV/c^{2})");
  WvsQ2_pion->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_pion->SetOption("COLZ");
  WvsQ2_pion->Write();

  W_pion->SetXTitle("W (GeV/c^{2})");
  W_pion->Write();

  Q2_pion->SetXTitle("Q^{2} (GeV^{2})");
  Q2_pion->Write();

  WvsQ2_NeutronPip->SetXTitle("W (GeV/c^{2})");
  WvsQ2_NeutronPip->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_NeutronPip->SetOption("COLZ");
  WvsQ2_NeutronPip->Write();

  WvsMM_NeutronPip->SetXTitle("W (GeV/c^{2})");
  WvsMM_NeutronPip->SetYTitle("MM (GeV)");
  WvsMM_NeutronPip->SetOption("COLZ");
  WvsMM_NeutronPip->Write();

  WvsMM2_NeutronPip->SetXTitle("W (GeV/c^{2})");
  WvsMM2_NeutronPip->SetYTitle("MM^{2} (GeV^{2})");
  WvsMM2_NeutronPip->SetOption("COLZ");
  WvsMM2_NeutronPip->Write();

  W_NeutronPip->SetXTitle("W (GeV/c^{2})");
  W_NeutronPip->Write();

  Q2_NeutronPip->SetXTitle("Q^{2} (GeV^{2})");
  Q2_NeutronPip->Write();

  WvsQ2_elastic->SetXTitle("W (GeV/c^{2})");
  WvsQ2_elastic->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_elastic->SetOption("COLZ");
  WvsQ2_elastic->Write();

  W_elastic->SetXTitle("W (GeV/c^{2})");
  W_elastic->Write();

  Q2_elastic->SetXTitle("Q^{2} (GeV^{2})");
  Q2_elastic->Write();

  elastic_MM->Write();
  elastic_phi->Write();
  elastic_thetaVsP->Write();

  WvsQ2_Ppi0->SetXTitle("W (GeV/c^{2})");
  WvsQ2_Ppi0->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_Ppi0->SetOption("COLZ");
  WvsQ2_Ppi0->Write();

  W_Ppi0->SetXTitle("W (GeV/c^{2})");
  W_Ppi0->Write();

  Q2_Ppi0->SetXTitle("Q^{2} (GeV^{2})");
  Q2_Ppi0->Write();

  WvsQ2_MM->SetXTitle("W (GeV/c^{2})");
  WvsQ2_MM->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_MM->SetOption("COLZ");
  WvsQ2_MM->Write();

  W_MM->SetXTitle("W (GeV/c^{2})");
  W_MM->Write();

  Q2_MM->SetXTitle("Q^{2} (GeV^{2})");
  Q2_MM->Write();
}

void Histogram::WvsQ2_sec_Write() {
  for (auto &sec : WvsQ2_channel_sec) {
    sec->SetXTitle("W (GeV/c^{2})");
    sec->SetYTitle("Q^{2} (GeV^{2})");
    sec->SetOption("COLZ");
    sec->Write();
  }
  for (auto &sec : W_channel_sec) {
    sec->SetXTitle("W (GeV/c^{2})");
    sec->Write();
  }

  auto canWQ2_channel = std::make_unique<TCanvas>("w_q2_sec_channel_can", "W vs Q2 N #pi^{+} Sec", 1600, 900);
  canWQ2_channel->Divide(2, 3);
  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    canWQ2_channel->cd(sec + 1);
    WvsQ2_channel_sec[sec]->SetXTitle("W (GeV/c^{2})");
    WvsQ2_channel_sec[sec]->SetYTitle("Q^{2} (GeV^{2})");
    WvsQ2_channel_sec[sec]->SetOption("COLZ");
    WvsQ2_channel_sec[sec]->Draw();
  }
  canWQ2_channel->Write();

  auto canW_channel = std::make_unique<TCanvas>("w_channel_sec_can", "W Sec", 1600, 900);
  canW_channel->Divide(2, 3);
  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    canW_channel->cd(sec + 1);
    W_channel_sec[sec]->SetXTitle("W (GeV/c^{2})");
    W_channel_sec[sec]->Draw();
  }
  canW_channel->Write();

  for (auto &sec : WvsQ2_sec) {
    sec->SetXTitle("W (GeV/c^{2})");
    sec->SetYTitle("Q^{2} (GeV^{2})");
    sec->SetOption("COLZ");
    sec->Write();
  }

  for (auto &sec : W_sec) {
    sec->SetXTitle("W (GeV/c^{2})");
    sec->Write();
  }

  auto canWQ2 = std::make_unique<TCanvas>("w_q2_sec_can", "W vs Q2 Sec", 1600, 900);
  canWQ2->Divide(2, 3);
  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    canWQ2->cd(sec + 1);
    WvsQ2_sec[sec]->SetXTitle("W (GeV/c^{2})");
    WvsQ2_sec[sec]->SetYTitle("Q^{2} (GeV^{2})");
    WvsQ2_sec[sec]->SetOption("COLZ");
    WvsQ2_sec[sec]->Draw();
  }
  canWQ2->Write();

  auto canW = std::make_unique<TCanvas>("w_sec_can", "W Sec", 1600, 900);
  canW->Divide(2, 3);
  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    canW->cd(sec + 1);
    W_sec[sec]->SetXTitle("W (GeV/c^{2})");
    W_sec[sec]->Draw();
  }
  canW->Write();
}

void Histogram::WvsQ2_binned_Write() {
  WvsQ2_binned->SetXTitle("W (GeV/c^{2})");
  WvsQ2_binned->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_binned->SetOption("COLZ");
  WvsQ2_binned->Write();

  for (auto &sec : W_binned) {
    sec->SetXTitle("W (GeV/c^{2})");
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
void Histogram::Fill_Missing_Mass(const std::shared_ptr<Reaction> &event) {
  Missing_Mass->Fill(event->MM());
  Missing_Mass_square->Fill(event->MM2());
  Missing_Mass_small->Fill(event->MM());
  if (event->sector() == 0 || event->sector() > NUM_SECTORS) {
    // std::cerr << "Error filling electron fid = " << event->sector() << std::endl;
    return;
  }
  Missing_Mass_small_sec[event->sector() - 1]->Fill(event->MM());
  Missing_Mass_Sq_small_sec[event->sector() - 1]->Fill(event->MM2());
}
void Histogram::Fill_Missing_Mass_strict(float mm, float mm2) {
  Missing_Mass_strict->Fill(mm);
  Missing_Mass_square_strict->Fill(mm2);
}

void Histogram::Fill_Missing_Mass_pi0(float mm, float mm2) {
  Missing_Mass_pi0->Fill(mm);
  Missing_Mass_square_pi0->Fill(mm2);
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

void Histogram::Fill_Mass_photons(const std::shared_ptr<Reaction> &_e) {
  Mass_pi0->Fill(_e->pi0_mass());
  Mass_square_pi0->Fill(_e->pi0_mass2());

  for (auto &&_m : _e->pair_mass()) {
    Mass_eta->Fill(_m);
    Mass_square_eta->Fill(_m * _m);
  }
}

void Histogram::Fill_Missing_Mass_square(float miss_mass_2) { Missing_Mass_square->Fill(miss_mass_2); }

void Histogram::Write_Missing_Mass() {
  Missing_Mass->SetXTitle("Mass (GeV)");
  Missing_Mass->Write();

  Missing_Mass_small->SetXTitle("Mass (GeV)");
  Missing_Mass_small->Write();

  for (auto &&sec : Missing_Mass_small_sec) {
    auto fit = std::make_shared<Fits>();
    fit->Set_min(0.8);
    fit->Set_max(1.1);
    fit->FitDeGauss(sec);
    sec->SetXTitle("Mass (GeV)");
    sec->Write();
  }

  for (auto &&sec : Missing_Mass_Sq_small_sec) {
    auto fit = std::make_shared<Fits>();
    fit->Set_min(0.7);
    fit->Set_max(1.0);
    fit->FitDeGauss(sec);
    sec->SetXTitle("Mass^2 (GeV)");
    sec->Write();
  }

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

  Mass_pi0->SetXTitle("Mass (GeV)");
  Mass_square_pi0->SetXTitle("Mass (GeV)");
  auto fit = std::make_shared<Fits>();
  fit->Set_min(0.1);
  fit->Set_max(0.2);
  fit->FitDeGauss(Mass_pi0);
  Mass_pi0->Write();

  auto fit_sq = std::make_shared<Fits>();
  fit_sq->Set_min(0.01);
  fit_sq->Set_max(0.03);
  fit_sq->FitDeGauss(Mass_square_pi0);
  Mass_square_pi0->Write();

  Mass_eta->SetXTitle("Mass (GeV)");
  Mass_square_eta->SetXTitle("Mass (GeV)");
  Mass_eta->Write();
  Mass_square_eta->Write();

  Missing_Mass_pi0_otherCut->Write();
  Missing_Mass_square_pi0_otherCut->Write();
  Mass_pi0_otherCut->Write();
  Mass_square_pi0_otherCut->Write();

  Missing_Mass_2pi->SetXTitle("Mass (GeV)");
  Missing_Mass_2pi->Write();
  Missing_Mass_square_2pi->SetXTitle("Mass^{2} (GeV^{2})");
  Missing_Mass_square_2pi->Write();
  Fit_Missing_Mass_WBinned.resize(W_BINS);
  TDirectory *mm_binned = RootOutputFile->mkdir("MM_binned");
  mm_binned->cd();
  for (int y = 0; y < W_BINS; y++) {
    // std::cout << Missing_Mass_WBinned[y]->GetEntries() << std::endl;
    Fit_Missing_Mass_WBinned[y] = new Fits();
    Fit_Missing_Mass_WBinned[y]->FitBreitWigner(Missing_Mass_WBinned[y]);
    // Fit_Missing_Mass_WBinned[y]->FitMissMass(Missing_Mass_WBinned[y]);

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
    delta_t_hist[0][jj] = new TH1D(Form("delta_t_p_%d", jj), Form("#Deltat P %d", jj), BINS, Dt_min, Dt_max);
    delta_t_hist[1][jj] = new TH1D(Form("delta_t_pip_%d", jj), Form("#Deltat #pi^{+} %d", jj), BINS, Dt_min, Dt_max);
    delta_t_hist[2][jj] =
        new TH1D(Form("delta_t_electron_%d", jj), Form("#Deltat electron %d", jj), BINS, Dt_min, Dt_max);
  }

  for (int jj = 0; jj < NUM_SECTORS; jj++) {
    for (int jjj = 0; jjj < SC_PADDLE_NUM; jjj++) {
      delta_t_sec_pad_hist[0][jj][jjj] = new TH2D(Form("delta_t_p_sec%d_pad%d", jj + 1, jjj + 1),
                                                  Form("#Deltat P Sector %d Paddle %d", jj + 1, jjj + 1), BINS / 2,
                                                  p_min, p_max, BINS / 2, Dt_min, Dt_max);
      delta_t_sec_pad_hist[1][jj][jjj] = new TH2D(Form("delta_t_pip_sec%d_pad%d", jj + 1, jjj + 1),
                                                  Form("#Deltat #pi^{+} Sector %d Paddle %d", jj + 1, jjj + 1),
                                                  BINS / 2, p_min, p_max, BINS / 2, Dt_min, Dt_max);

      delta_t_sec_pad_hist[2][jj][jjj] = new TH2D(Form("delta_t_electron_sec%d_pad%d", jj + 1, jjj + 1),
                                                  Form("#Deltat electron Sector %d Paddle %d", jj + 1, jjj + 1),
                                                  BINS / 2, p_min, p_max, BINS / 2, Dt_min, Dt_max);
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
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

  TF1 *peak = new TF1("peak", "gaus", -1, 1);
  // TF1 *peak = new TF1("peak", func::peak, -1, 1, 3);
  peak->SetParNames("constant", "mean", "#sigma");
  // Bin 50 = 0.5GeV, Bin 300 = 3 GeV
  delta_t_mass_P->FitSlicesY(peak, 0, 300, 10, "QRG5");
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
  P->Write();
  M->Write();
  Proton_Pos_fit->Write();
  Proton_Neg_fit->Write();
  delete P;
  delete M;

  TCanvas *dt_proton_canvas = new TCanvas("dt_proton_canvas", "#dt Proton", 1280, 720);
  dt_proton_canvas->cd();
  delta_t_mass_P->Draw();
  Proton_Pos_fit->Draw("same");
  Proton_Neg_fit->Draw("same");
  P->Draw("*same");
  M->Draw("*same");
  dt_proton_canvas->Write();
  delete dt_proton_canvas;

  delta_t_mass_PIP->FitSlicesY(peak, 0, 300, 10, "QRG5");
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
  P_pip->Write();
  M_pip->Write();
  Pip_Pos_fit->Write();
  Pip_Neg_fit->Write();
  delete P_pip;
  delete M_pip;

  TCanvas *dt_Pip_canvas = new TCanvas("dt_Pip_canvas", "#dt #pi^{+}", 1280, 720);
  dt_Pip_canvas->cd();
  delta_t_mass_PIP->Draw();
  Pip_Pos_fit->Draw("same");
  Pip_Neg_fit->Draw("same");
  P_pip->Draw("*same");
  M_pip->Draw("*same");
  dt_Pip_canvas->Write();
  delete dt_Pip_canvas;
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

  // delta_t_slice_fit();

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
      delta_t_cut[j][jj] = std::make_unique<Fits>();
      delta_t_cut[j][jj]->Set_min(fit_dt_min);
      delta_t_cut[j][jj]->Set_max(fit_dt_max);
      delta_t_cut[j][jj]->FitGaus(delta_t_hist[j][jj]);

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

      can_dt[sec_i][particle_i] = new TCanvas(Form("Sector %d %s", sec_i + 1, P_PIP_E.c_str()),
                                              Form("Sector %d %s", sec_i + 1, P_PIP_E.c_str()), 1200, 800);
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

  if (cc_sector <= NUM_SECTORS && cc_segment <= segment && cc_pmt < PMT) {
    cc_hist[cc_sector - 1][cc_segment - 1][cc_pmt]->Fill(cc_nphe);
    cc_hist_allSeg[cc_sector - 1][cc_pmt]->Fill(cc_nphe);

    Theta_CC->Fill(cc_segment, theta_cc);
    Theta_CC_Sec[cc_sector - 1]->Fill(cc_segment, theta_cc);
    Theta_CC_Sec_cut[cc_sector - 1]->Fill(cc_segment, theta_cc);
  }
}

void Histogram::makeHists_CC() {
  Theta_CC_Sec.reserve(6);
  Theta_CC_Sec_cut.reserve(6);
  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    Theta_CC_Sec[sec_i] = std::make_shared<TH2D>(Form("Theta_CC_sec%d", sec_i + 1),
                                                 Form("Theta CC sector %d", sec_i + 1), 20, 0.0, 20.0, 60, 0.0, 60.0);
    Theta_CC_Sec_cut[sec_i] = std::make_shared<TH2D>(
        Form("Theta_CC_sec_cut%d", sec_i + 1), Form("Theta CC sector cut %d", sec_i + 1), 20, 0.0, 20.0, 60, 0.0, 60.0);
    for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
      cc_hist_allSeg[sec_i][pmt_i] =
          std::make_shared<TH1D>(Form("CC_sec%d_%s", sec_i + 1, L_R_C[pmt_i].c_str()),
                                 Form("CC sector %d %s", sec_i + 1, L_R_C[pmt_i].c_str()), bins_CC, CC_min, CC_max);
      for (int seg_i = 0; seg_i < segment; seg_i++) {
        cc_hist[sec_i][seg_i][pmt_i] = std::make_shared<TH1D>(
            Form("CC_sec%d_seg%d_%s", sec_i + 1, seg_i + 1, L_R_C[pmt_i].c_str()),
            Form("CC sector %d segment %d %s", sec_i + 1, seg_i + 1, L_R_C[pmt_i].c_str()), bins_CC, CC_min, CC_max);
      }
    }
  }
}

void Histogram::Theta_CC_Write() {
  Theta_CC->SetXTitle("CC segment");
  Theta_CC->SetYTitle("#theta_CC");
  Theta_CC->SetOption("COLZ");
  Theta_CC->Write();

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

void Histogram::CC_Write() {
  // cc_sparse->Write();
  std::unique_ptr<Fits> cc_fits[NUM_SECTORS][segment][PMT];
  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    for (int pmt_i = 1; pmt_i < PMT; pmt_i++) {
      cc_hist_allSeg[sec_i][pmt_i]->SetYTitle("number photoelectrons");
      cc_hist_allSeg[sec_i][pmt_i]->Write();
      for (int seg_i = 0; seg_i < segment; seg_i++) {
        cc_fits[sec_i][seg_i][pmt_i] = std::make_unique<Fits>();

        // cc_fits[sec_i][seg_i][pmt_i]->Set_lineColor(9);
        // cc_fits[sec_i][seg_i][pmt_i]->Set_min(0.0);
        // cc_fits[sec_i][seg_i][pmt_i]->Set_max(30.0);
        // cc_fits[sec_i][seg_i][pmt_i]->FitLandau(cc_hist[sec_i][seg_i][pmt_i].get());

        // cc_fits[sec_i][seg_i][pmt_i]->Set_lineColor(8);
        // cc_fits[sec_i][seg_i][pmt_i]->Set_min(30.0);
        // cc_fits[sec_i][seg_i][pmt_i]->Set_min(0.0);
        // cc_fits[sec_i][seg_i][pmt_i]->Set_max(250.0);
        // cc_fits[sec_i][seg_i][pmt_i]->FitGaus(cc_hist[sec_i][seg_i][pmt_i].get());

        cc_fits[sec_i][seg_i][pmt_i]->Set_lineColor(10);
        cc_fits[sec_i][seg_i][pmt_i]->FitLandauGaus(cc_hist[sec_i][seg_i][pmt_i].get());

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
    Theta_CC_0[sec_i] = (TH1D *)gDirectory->Get(Form("Theta_CC_sec%d_%d", sec_i + 1, 0));
    Theta_CC_1[sec_i] = (TH1D *)gDirectory->Get(Form("Theta_CC_sec%d_%d", sec_i + 1, 1));
    Theta_CC_2[sec_i] = (TH1D *)gDirectory->Get(Form("Theta_CC_sec%d_%d", sec_i + 1, 2));

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

    // Theta_CC_Pos_fit[sec_i] = new TF1(get_name, func::pol1, 0, 14, 2);
    Theta_CC_Pos_fit[sec_i] = new TF1(Form("Theta_CC_sec%d_pos_fit", sec_i + 1), func::theta_cc_fit, 1, 14, 4);
    Theta_CC_Pos_fit[sec_i]->SetParNames("intercept", "slope", "a", "exp^c*x");
    Theta_CC_Pos_fit[sec_i]->SetParLimits(3, 0, 10);

    // Theta_CC_Neg_fit[sec_i] = new TF1(get_name, func::pol1, 0, 14, 2);
    Theta_CC_Neg_fit[sec_i] = new TF1(Form("Theta_CC_sec%d_neg_fit", sec_i + 1), func::theta_cc_fit, 1, 14, 4);
    Theta_CC_Neg_fit[sec_i]->SetParNames("intercept", "slope", "a", "exp^c*x");
    Theta_CC_Neg_fit[sec_i]->SetParLimits(3, 0, 10);
    CC_P[sec_i]->Fit(Theta_CC_Pos_fit[sec_i], "QRM+", "", 1, 14);
    CC_M[sec_i]->Fit(Theta_CC_Neg_fit[sec_i], "QRM+", "", 1, 14);

    Theta_CC_canvas[sec_i] =
        new TCanvas(Form("Theta_CC_sec%d", sec_i + 1), Form("Theta CC sec_i %d", sec_i + 1), 1280, 720);
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
  for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
    for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
      can[sec_i][pmt_i] = new TCanvas(Form("Sector %d %s", sec_i + 1, L_R_C[pmt_i].c_str()),
                                      Form("Sector %d %s", sec_i + 1, L_R_C[pmt_i].c_str()), 1200, 800);
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
  hadron_fid_hist[fid_part::hadron] =
      std::make_unique<TH2D>("hadron_fid", "hadron_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
  hadron_fid_hist[fid_part::proton] =
      std::make_unique<TH2D>("proton_fid", "proton_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);
  hadron_fid_hist[fid_part::pip] =
      std::make_unique<TH2D>("pip_fid", "pip_fid", BINS, phi_min, phi_max, BINS, theta_min, theta_max);

  hadron_fid_xy_hist[fid_part::hadron] =
      std::make_unique<TH2D>("hadron_fid_xy", "hadron_fid_xy", BINS, -200, 200, BINS, 0, 500);
  hadron_fid_xy_hist[fid_part::proton] =
      std::make_unique<TH2D>("proton_fid_xy", "proton_fid_xy", BINS, -200, 200, BINS, 0, 500);
  hadron_fid_xy_hist[fid_part::pip] = std::make_unique<TH2D>("pip_fid_xy", "pip_fid_xy", BINS, -200, 200, BINS, 0, 500);

  cerenkov_fid.reserve(NUM_SECTORS);
  cerenkov_fid_cut.reserve(NUM_SECTORS);
  cerenkov_fid_anti.reserve(NUM_SECTORS);

  fid_xy.reserve(NUM_SECTORS);
  fid_xy_cut.reserve(NUM_SECTORS);
  fid_xy_anti.reserve(NUM_SECTORS);

  electron_fid_sec_hist.reserve(NUM_SECTORS);
  electron_fid_sec_cut_hist.reserve(NUM_SECTORS);
  electron_fid_sec_anti_hist.reserve(NUM_SECTORS);

  fid_xy_hist = std::make_shared<TH2D>("fid_dc_xy", "fid_dc_xy", BINS, -150, 150, BINS, 0, 300);
  fid_xy_cut_hist = std::make_shared<TH2D>("fid_dc_xy_cut", "fid_dc_xy_cut", BINS, -150, 150, BINS, 0, 300);
  fid_xy_anti_hist = std::make_shared<TH2D>("fid_dc_xy_anti", "fid_dc_xy_anti", BINS, -150, 150, BINS, 0, 300);
  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    cerenkov_fid[sec_i] = std::make_shared<TH2D>(Form("fid_cher_xy_%d", sec_i + 1), Form("fid_cher_xy_%d", sec_i + 1),
                                                 BINS, -150, 150, BINS, 0, 300);
    cerenkov_fid_cut[sec_i] = std::make_shared<TH2D>(
        Form("fid_cher_xy_cut_%d", sec_i + 1), Form("fid_cher_xy_cut_%d", sec_i + 1), BINS, -150, 150, BINS, 0, 300);
    cerenkov_fid_anti[sec_i] = std::make_shared<TH2D>(
        Form("fid_cher_xy_anti_%d", sec_i + 1), Form("fid_cher_xy_anti_%d", sec_i + 1), BINS, -150, 150, BINS, 0, 300);
    fid_xy[sec_i] = std::make_shared<TH2D>(Form("fid_dc_xy_%d", sec_i + 1), Form("fid_dc_xy_%d", sec_i + 1), BINS, -150,
                                           150, BINS, 0, 300);
    fid_xy_cut[sec_i] = std::make_shared<TH2D>(Form("fid_dc_xy_cut_%d", sec_i + 1), Form("fid_dc_xy_cut_%d", sec_i + 1),
                                               BINS, -150, 150, BINS, 0, 300);
    fid_xy_anti[sec_i] = std::make_shared<TH2D>(Form("fid_dc_xy_anti_%d", sec_i + 1),
                                                Form("fid_dc_xy_anti_%d", sec_i + 1), BINS, -150, 150, BINS, 0, 300);

    electron_fid_sec_hist[sec_i] =
        std::make_shared<TH2D>(Form("electron_fid_sec%d", sec_i + 1), Form("electron_fid_sec%d", sec_i + 1), BINS,
                               min_phi[sec_i], max_phi[sec_i], BINS, theta_min, theta_max);

    electron_fid_sec_cut_hist[sec_i] =
        std::make_shared<TH2D>(Form("electron_fid_sec_cut_%d", sec_i + 1), Form("electron_fid_sec_cut_%d", sec_i + 1),
                               BINS, min_phi[sec_i], max_phi[sec_i], BINS, theta_min, theta_max);

    electron_fid_sec_anti_hist[sec_i] =
        std::make_shared<TH2D>(Form("electron_fid_sec_anti_%d", sec_i + 1), Form("electron_fid_sec_anti_%d", sec_i + 1),
                               BINS, min_phi[sec_i], max_phi[sec_i], BINS, theta_min, theta_max);

    for (int t = 0; t < 3; t++) {
      hadron_fid_sec_hist[t][sec_i] =
          std::make_shared<TH2D>(Form("hadron_fid_sec%d_%d", sec_i + 1, t), Form("hadron_fid_sec%d_%d", sec_i + 1, t),
                                 BINS, min_phi[sec_i], max_phi[sec_i], BINS, theta_min, theta_max);
      hadron_fid_xy_sec_hist[t][sec_i] =
          std::make_shared<TH2D>(Form("hadron_fid_xy_sec%d_%d", sec_i + 1, t),
                                 Form("hadron_fid_xy_sec%d_%d", sec_i + 1, t), BINS, -200, 200, BINS, 0, 500);
    }
  }
}

void Histogram::Fill_electron_fid(const std::shared_ptr<Branches> &_data, const std::shared_ptr<Reaction> &_r) {
  if ((_r->cc_y() == 0 && _r->cc_x() == 0) || (_data->dc_ysc(0) == 0 && _data->dc_xsc(0) == 0)) return;

  float theta = physics::theta_calc(_data->cz(0));
  float phi = physics::phi_calc(_data->cx(0), _data->cy(0));
  int sector = _data->dc_sect(0);

  electron_fid_hist->Fill(phi, theta);
  fid_xy_hist->Fill(_data->dc_ysc(0), _data->dc_xsc(0));
  if (sector == 0 || sector > NUM_SECTORS) {
    // std::cerr << "Error filling electron fid = " << sector << std::endl;
    return;
  }

  cerenkov_fid[sector - 1]->Fill(_r->cc_y(), _r->cc_x());
  fid_xy[sector - 1]->Fill(_data->dc_ysc(0), _data->dc_xsc(0));
  electron_fid_sec_hist[sector - 1]->Fill(phi, theta);
}

void Histogram::Fill_electron_fid_cut(const std::shared_ptr<Branches> &_data, const std::shared_ptr<Reaction> &_r) {
  if ((_r->cc_y() == 0 && _r->cc_x() == 0) || (_data->dc_ysc(0) == 0 && _data->dc_xsc(0) == 0)) return;
  float theta = physics::theta_calc(_data->cz(0));
  float phi = physics::phi_calc(_data->cx(0), _data->cy(0));
  int sector = _data->dc_sect(0);

  electron_fid_cut_hist->Fill(phi, theta);
  fid_xy_cut_hist->Fill(_data->dc_ysc(0), _data->dc_xsc(0));
  if (sector == 0 || sector > NUM_SECTORS) {
    // std::cerr << "Error filling electron fid = " << sector << std::endl;
    return;
  }

  cerenkov_fid_cut[sector - 1]->Fill(_r->cc_y(), _r->cc_x());
  fid_xy_cut[sector - 1]->Fill(_data->dc_ysc(0), _data->dc_xsc(0));
  electron_fid_sec_cut_hist[sector - 1]->Fill(phi, theta);
}

void Histogram::Fill_electron_fid_anti(const std::shared_ptr<Branches> &_data, const std::shared_ptr<Reaction> &_r) {
  if ((_r->cc_y() == 0 && _r->cc_x() == 0) || (_data->dc_ysc(0) == 0 && _data->dc_xsc(0) == 0)) return;
  float theta = physics::theta_calc(_data->cz(0));
  float phi = physics::phi_calc(_data->cx(0), _data->cy(0));
  int sector = _data->dc_sect(0);

  electron_fid_anti_hist->Fill(phi, theta);
  fid_xy_anti_hist->Fill(_data->dc_ysc(0), _data->dc_xsc(0));
  if (sector == 0 || sector > NUM_SECTORS) {
    // std::cerr << "Error filling electron fid = " << sector << std::endl;
    return;
  }

  cerenkov_fid_anti[sector - 1]->Fill(_r->cc_y(), _r->cc_x());
  fid_xy_anti[sector - 1]->Fill(_data->dc_ysc(0), _data->dc_xsc(0));
  electron_fid_sec_anti_hist[sector - 1]->Fill(phi, theta);
}

void Histogram::Fill_neutron_fid(float theta, float phi, int sector) {
  neutron_fid_hist->Fill(phi, theta);
  if (sector == 0 || sector > NUM_SECTORS) {
    // std::cerr << "Error filling electron fid = " << sector << std::endl;
    return;
  }
  electron_fid_sec_hist[sector - 1]->Fill(phi, theta);
}

void Histogram::Fill_hadron_fid(const std::shared_ptr<Branches> &_data, int part_num, int id) {
  short sector = _data->dc_sect(part_num);
  if (sector == 0 || sector > NUM_SECTORS) return;

  float phi = physics::phi_calc(_data->cx(part_num), _data->cy(part_num));
  float theta = physics::theta_calc(_data->cz(part_num));

  hadron_fid_hist[fid_part::hadron]->Fill(phi, theta);
  hadron_fid_sec_hist[fid_part::hadron][sector - 1]->Fill(phi, theta);

  hadron_fid_xy_hist[fid_part::hadron]->Fill(_data->dc_ysc(part_num), _data->dc_xsc(part_num));
  hadron_fid_xy_sec_hist[fid_part::hadron][sector - 1]->Fill(_data->dc_ysc(part_num), _data->dc_xsc(part_num));

  if (id == PROTON) {
    hadron_fid_hist[fid_part::proton]->Fill(phi, theta);
    hadron_fid_sec_hist[fid_part::proton][sector - 1]->Fill(phi, theta);
    hadron_fid_xy_hist[fid_part::proton]->Fill(_data->dc_ysc(part_num), _data->dc_xsc(part_num));
    hadron_fid_xy_sec_hist[fid_part::proton][sector - 1]->Fill(_data->dc_ysc(part_num), _data->dc_xsc(part_num));

  } else if (id == PIP) {
    hadron_fid_hist[fid_part::pip]->Fill(phi, theta);
    hadron_fid_sec_hist[fid_part::pip][sector - 1]->Fill(phi, theta);
    hadron_fid_xy_hist[fid_part::pip]->Fill(_data->dc_ysc(part_num), _data->dc_xsc(part_num));
    hadron_fid_xy_sec_hist[fid_part::pip][sector - 1]->Fill(_data->dc_ysc(part_num), _data->dc_xsc(part_num));
  }
}

void Histogram::Fid_Write() {
  electron_fid_hist->SetYTitle("#theta");
  electron_fid_hist->SetXTitle("#phi");
  electron_fid_hist->SetOption("COLZ");
  electron_fid_hist->Write();

  electron_fid_cut_hist->SetYTitle("#theta");
  electron_fid_cut_hist->SetXTitle("#phi");
  electron_fid_cut_hist->SetOption("COLZ");
  electron_fid_cut_hist->Write();

  electron_fid_anti_hist->SetYTitle("#theta");
  electron_fid_anti_hist->SetXTitle("#phi");
  electron_fid_anti_hist->SetOption("COLZ");
  electron_fid_anti_hist->Write();

  for (int t = 0; t < 3; t++) {
    hadron_fid_hist[t]->SetYTitle("#theta");
    hadron_fid_hist[t]->SetXTitle("#phi");
    hadron_fid_hist[t]->SetOption("COLZ");
    hadron_fid_hist[t]->Write();

    hadron_fid_xy_hist[t]->SetYTitle("#theta");
    hadron_fid_xy_hist[t]->SetXTitle("#phi");
    hadron_fid_xy_hist[t]->SetOption("COLZ");
    hadron_fid_xy_hist[t]->Write();
  }

  fid_xy_hist->SetXTitle("X");
  fid_xy_hist->SetYTitle("Y");
  fid_xy_hist->SetOption("COLZ");
  fid_xy_hist->Write();

  fid_xy_cut_hist->SetXTitle("X");
  fid_xy_cut_hist->SetYTitle("Y");
  fid_xy_cut_hist->SetOption("COLZ");
  fid_xy_cut_hist->Write();

  fid_xy_anti_hist->SetXTitle("X");
  fid_xy_anti_hist->SetYTitle("Y");
  fid_xy_anti_hist->SetOption("COLZ");
  fid_xy_anti_hist->Write();

  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    cerenkov_fid[sec_i]->SetXTitle("X");
    cerenkov_fid[sec_i]->SetYTitle("Y");
    cerenkov_fid[sec_i]->SetOption("COLZ");
    cerenkov_fid[sec_i]->Write();

    cerenkov_fid_cut[sec_i]->SetXTitle("X");
    cerenkov_fid_cut[sec_i]->SetYTitle("Y");
    cerenkov_fid_cut[sec_i]->SetOption("COLZ");
    cerenkov_fid_cut[sec_i]->Write();

    cerenkov_fid_anti[sec_i]->SetXTitle("X");
    cerenkov_fid_anti[sec_i]->SetYTitle("Y");
    cerenkov_fid_anti[sec_i]->SetOption("COLZ");
    cerenkov_fid_anti[sec_i]->Write();

    fid_xy[sec_i]->SetXTitle("X");
    fid_xy[sec_i]->SetYTitle("Y");
    fid_xy[sec_i]->SetOption("COLZ");
    fid_xy[sec_i]->Write();

    fid_xy_cut[sec_i]->SetXTitle("X");
    fid_xy_cut[sec_i]->SetYTitle("Y");
    fid_xy_cut[sec_i]->SetOption("COLZ");
    fid_xy_cut[sec_i]->Write();

    fid_xy_anti[sec_i]->SetXTitle("X");
    fid_xy_anti[sec_i]->SetYTitle("Y");
    fid_xy_anti[sec_i]->SetOption("COLZ");
    fid_xy_anti[sec_i]->Write();
  }

  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    electron_fid_sec_hist[sec_i]->SetXTitle("#phi");
    electron_fid_sec_hist[sec_i]->SetYTitle("#theta");
    electron_fid_sec_hist[sec_i]->SetOption("COLZ");
    electron_fid_sec_hist[sec_i]->Write();

    electron_fid_sec_cut_hist[sec_i]->SetXTitle("#phi");
    electron_fid_sec_cut_hist[sec_i]->SetYTitle("#theta");
    electron_fid_sec_cut_hist[sec_i]->SetOption("COLZ");
    electron_fid_sec_cut_hist[sec_i]->Write();
    for (int t = 0; t < 3; t++) {
      hadron_fid_sec_hist[t][sec_i]->SetYTitle("#theta");
      hadron_fid_sec_hist[t][sec_i]->SetXTitle("#phi");
      hadron_fid_sec_hist[t][sec_i]->SetOption("COLZ");
      hadron_fid_sec_hist[t][sec_i]->Write();

      hadron_fid_xy_sec_hist[t][sec_i]->SetYTitle("#theta");
      hadron_fid_xy_sec_hist[t][sec_i]->SetXTitle("#phi");
      hadron_fid_xy_sec_hist[t][sec_i]->SetOption("COLZ");
      hadron_fid_xy_sec_hist[t][sec_i]->Write();
    }
  }

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
    electron_fid_can[sec_i] =
        new TCanvas(Form("electron_fid_sector_%d", sec_i + 1), Form("electron_fid_sector_%d", sec_i + 1), 1280, 720);

    for (int slice = start_slice; slice < FID_SLICES; slice++) {
      // std::cout << "Fid Write " << sec_i << "\t" << slice << std::endl;
      electron_fid_sec_slice[sec_i][slice] = (TH1D *)electron_fid_sec_hist[sec_i]->ProjectionX(
          Form("electron_fid_sec_%d_%d", sec_i + 1, slice + 1), slice_width * slice,
          slice_width * slice + (slice_width - 1));
      electron_fid_sec_slice[sec_i][slice]->Rebin(4);
      /*
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
      */
    }

    fid[sec_i] = new TGraph(FID_SLICES * 2, x, y);
    // FidGraph[sec_i] = std::make_unique<Fits>();

    // FidGraph[sec_i]->Set_min(min_phi[sec_i]);
    // FidGraph[sec_i]->Set_max(max_phi[sec_i]);
    // FidGraph[sec_i]->FitFiducial(fid[sec_i], sec_i);
    // FidGraph[sec_i]->FitPoly_fid(fid[sec_i]);

    electron_fid_can[sec_i]->cd();

    electron_fid_sec_hist[sec_i]->Draw();
    fid[sec_i]->Draw("*same");
    electron_fid_can[sec_i]->Write();
  }
}

void Histogram::fid_canvas() {
  TCanvas *can[NUM_SECTORS];
  char can_name[500];

  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    can[sec_i] = new TCanvas(Form("Electron Fid Sector %d Slices", sec_i + 1),
                             Form("Electron Fid Sector %d Slices", sec_i + 1), 1600, 900);
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
    EC_hist[n] = new TH1D(Form("ec_%d", n), Form("Sampling Fraction %d", n), BINS, EC_min, EC_max);
    EC_hist_cut[n] = new TH1D(Form("ec_cut_%d", n), Form("Sampling Fraction cut %d", n), BINS, EC_min, EC_max);
  }

  Theta_vs_mom_sec.reserve(NUM_SECTORS);
  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    Theta_vs_mom_sec[sec] = std::make_shared<TH2D>(Form("Theta_vs_mom_%d", sec + 1), Form("Theta_vs_mom_%d", sec + 1),
                                                   BINS, p_min, p_max, BINS, 0, 70);
  }
}

void Histogram::EC_fill(float etot, float momentum) {
  EC_etot_vs_P->Fill(momentum, etot);
  float sampling_frac = etot / momentum;
  EC_sampling_fraction->Fill(momentum, sampling_frac);

  EC_tot_energy->Fill(etot);

  for (int n = 0; n < NUM_POINTS; n++) {
    if (momentum > n * bin_width && momentum <= (n + 1) * bin_width) {
      EC_hist[n]->Fill(sampling_frac);
    }
  }
}

void Histogram::EC_inout(float Ein, float Eout) {
  if (Eout > 0 && Ein > 0) ECin_ECout->Fill(Ein, Eout);
}

void Histogram::TM_Fill(float momentum, float theta) { Theta_vs_mom->Fill(momentum, theta); }
void Histogram::Theta_vs_p_Fill(const std::shared_ptr<Branches> &_data) {
  short sector = _data->dc_sect(0);
  if (sector == 0 || sector > NUM_SECTORS || sector == int(NULL)) return;
  Theta_vs_mom_sec[sector - 1]->Fill(_data->p(0), physics::theta_calc(_data->cz(0)));
}

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
  if (EC_sampling_fraction_cut->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  // Header *fit_functions = new Header("EC_fit_functions.hpp", "FF");

  TF1 *peak = new TF1("peak", "gaus", 0.1, 0.5);
  EC_sampling_fraction_cut->FitSlicesY(peak, 0, -1, 0, "QRG5");
  TH1D *EC_sampling_fraction_cut_0 = (TH1D *)gDirectory->Get("EC_sampling_fraction_cut_0");
  TH1D *EC_sampling_fraction_cut_1 = (TH1D *)gDirectory->Get("EC_sampling_fraction_cut_1");
  TH1D *EC_sampling_fraction_cut_2 = (TH1D *)gDirectory->Get("EC_sampling_fraction_cut_2");
  float x[BINS];
  float y_plus[BINS];
  float y_minus[BINS];
  int num = 0;
  for (int i = 0; i < BINS; i++) {
    if (EC_sampling_fraction_cut_1->GetBinContent(i) != 0) {
      // Get momentum from bin center
      x[num] = (float)EC_sampling_fraction_cut_1->GetBinCenter(i);
      // mean + 3sigma
      y_plus[num] = (float)EC_sampling_fraction_cut_1->GetBinContent(i) +
                    4.0 * (float)EC_sampling_fraction_cut_2->GetBinContent(i);
      // mean - 3simga
      y_minus[num] = (float)EC_sampling_fraction_cut_1->GetBinContent(i) -
                     4.0 * (float)EC_sampling_fraction_cut_2->GetBinContent(i);
      num++;
    }
  }

  TGraph *EC_P = new TGraph(num, x, y_plus);
  TGraph *EC_M = new TGraph(num, x, y_minus);
  EC_P->SetName("Positive_EC_graph");
  EC_M->SetName("Negative_EC_graph");

  TF1 *EC_P_fit = new TF1("EC_P_fit", func::ec_fit_func, 0.25, 4.0, 3);
  TF1 *EC_M_fit = new TF1("EC_M_fit", func::ec_fit_func, 0.25, 4.0, 3);

  EC_P_fit->SetParameters(0.3296, 0.002571, 4.8e-7);
  EC_M_fit->SetParameters(0.1715, 0.02044, -1.581e-5);
  EC_P_fit->SetParLimits(2, 1.0e-7, 5.0e-7);
  EC_M_fit->SetParLimits(2, -2.0e-5, -1.0e-5);

  EC_P->Fit(EC_P_fit, "QMRG+", "", 0.75, 3.75);
  EC_M->Fit(EC_M_fit, "QMRG+", "", 0.75, 3.75);

  EC_P_fit->Write();
  EC_M_fit->Write();
  EC_P->Write();
  EC_M->Write();

  TCanvas *EC_canvas = new TCanvas("EC_canvas", "EC canvas", 1280, 720);
  EC_canvas->cd();
  EC_sampling_fraction_cut->Draw();
  // EC_sampling_fraction->Draw();
  EC_P_fit->Draw("same");
  // std::cout << "EC_P_fit " << EC_P_fit->GetParameter(0) << "," << EC_P_fit->GetParameter(1) << ","
  //          << EC_P_fit->GetParameter(2) << std::endl;
  EC_M_fit->Draw("same");
  // std::cout << "EC_M_fit " << EC_M_fit->GetParameter(0) << "," << EC_M_fit->GetParameter(1) << ","
  //          << EC_M_fit->GetParameter(2) << std::endl;
  EC_P->Draw("*same");
  EC_M->Draw("*same");
  EC_canvas->Write();
  delete EC_canvas;
  delete EC_P;
  delete EC_M;
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
    if (EC_hist[n]->GetEntries() < 1000) continue;
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
  EC_tot_energy->SetXTitle("Energy (GeV)");
  EC_tot_energy->Write();

  EC_etot_vs_P->SetYTitle("Energy (GeV)");
  EC_etot_vs_P->SetXTitle("Momentum (GeV)");
  EC_etot_vs_P->SetOption("COLZ");
  EC_etot_vs_P->Write();

  EC_sampling_fraction->SetXTitle("Momentum (GeV)");
  EC_sampling_fraction->SetYTitle("Sampling Fraction");
  EC_sampling_fraction->SetOption("COLZ");
  EC_sampling_fraction->Write();

  EC_sampling_fraction_cut->SetXTitle("Momentum (GeV)");
  EC_sampling_fraction_cut->SetYTitle("Sampling Fraction");
  EC_sampling_fraction_cut->SetOption("COLZ");
  EC_sampling_fraction_cut->Write();
  EC_slice_fit();

  Theta_vs_mom->SetXTitle("Momentum (GeV)");
  Theta_vs_mom->SetYTitle("Theta #theta");
  Theta_vs_mom->SetOption("COLZ");
  Theta_vs_mom->Write();

  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    Theta_vs_mom_sec[sec_i]->SetXTitle("Momentum (GeV)");
    Theta_vs_mom_sec[sec_i]->SetYTitle("Theta #theta");
    Theta_vs_mom_sec[sec_i]->SetOption("COLZ LOGZ");
    Theta_vs_mom_sec[sec_i]->Write();
  }

  ECin_ECout->SetXTitle("EC_{inner}");
  ECin_ECout->SetYTitle("EC_{outer}");
  ECin_ECout->SetOption("COLZ");
  ECin_ECout->Write();
}

void Histogram::Fill_Beam_Position_cut(const std::shared_ptr<Branches> &_data) {
  Beam_Position_cut->Fill(_data->dc_vx(0), _data->dc_vy(0));
  Beam_Position_X_cut->Fill(_data->dc_vx(0));
  Beam_Position_Y_cut->Fill(_data->dc_vy(0));
  Beam_Position_Z_cut->Fill(_data->dc_vz(0));
}

void Histogram::Fill_Beam_Position(const std::shared_ptr<Branches> &_data) {
  Beam_Position->Fill(_data->dc_vx(0), _data->dc_vy(0));
  Beam_Position_X->Fill(_data->dc_vx(0));
  Beam_Position_Y->Fill(_data->dc_vy(0));
  Beam_Position_Z->Fill(_data->dc_vz(0));

  // Phi vs vertex
  target_vertex_vz_phi->Fill(_data->dc_vz(0), physics::phi_calc(_data->cx(0), _data->cy(0)));
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

  Beam_Position_cut->SetXTitle("X");
  Beam_Position_cut->SetYTitle("Y");
  Beam_Position_cut->SetOption("COLZ");
  Beam_Position_cut->Write();

  Beam_Position_X_cut->SetXTitle("X");
  Beam_Position_X_cut->Write();

  Beam_Position_Y_cut->SetXTitle("Y");
  Beam_Position_Y_cut->Write();

  Beam_Position_Z_cut->SetXTitle("Z");
  Beam_Position_Z_cut->Write();
}

void Histogram::Fill_Target_Vertex(const std::shared_ptr<Branches> &_data) {
  for (int part_num = 1; part_num < _data->gpart(); part_num++) {
    float vertex_x = _data->dc_vx(part_num);
    float vertex_y = _data->dc_vy(part_num);
    float vertex_z = _data->dc_vz(part_num);

    if (0 == vertex_x) return;
    if (0 == vertex_y && 0 == vertex_z) return;

    target_vertex_X->Fill(vertex_x);
    target_vertex_Y->Fill(vertex_y);
    target_vertex_Z->Fill(vertex_z);
    target_vertex_xy->Fill(vertex_x, vertex_y);
    target_vertex_zy->Fill(vertex_z, vertex_y);
    target_vertex_zx->Fill(vertex_z, vertex_x);
  }
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

  target_vertex_vz_phi->SetXTitle("vz");
  target_vertex_vz_phi->SetYTitle("#phi");
  target_vertex_vz_phi->SetOption("COLZ");
  target_vertex_vz_phi->Write();
}

void Histogram::Fill_E_Prime(const LorentzVector &e_prime) {
  if (e_prime.E() > 0.1) energy_no_cuts->Fill(e_prime.E());
}
void Histogram::Fill_E_Prime_fid(const LorentzVector &e_prime) {
  if (e_prime.E() > 0.1) energy_fid_cuts->Fill(e_prime.E());
}
void Histogram::Fill_E_Prime_channel(const LorentzVector &e_prime) {
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
  RootOutputFile->cd();
  ndhist_mc->Write();
  // THn *acceptance = (THn *)ndhist_mc->Clone("Acceptance");
  // acceptance->Scale(1 / ndhist_mc->GetEntries());
  // ndhist->Scale(1 / ndhist->GetEntries());
  // acceptance->Divide(ndhist.get());
  // acceptance->Write();

  Histogram::Write();
  // Start of cuts
  auto MM_neutron_cut = std::make_unique<Fits>();
  MM_neutron_cut->Set_min(0.8);
  MM_neutron_cut->Set_max(1.2);
  MM_neutron_cut->FitBreitWigner(Missing_Mass.get());
  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;

  TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2 MC");
  WvsQ2_folder->cd();
  WvsQ2_MC_Write();

  TDirectory *delta_mom = RootOutputFile->mkdir("delta_mom");
  delta_mom->cd();
  Write_DeltaP();
  std::cerr << BOLDBLUE << "Done!!!" << DEF << std::endl;
}

void mcHistogram::makeMCHists() {
  ndhist_mc = std::make_unique<THnSparseD>("ndhist_mc", "ndhist_mc", DIMENSIONS, nbins, xmin, xmax);
  ndhist_mc->GetAxis(0)->SetName("W");
  ndhist_mc->GetAxis(1)->SetName("Q2");
  ndhist_mc->GetAxis(2)->SetName("cos_Theta_star");
  ndhist_mc->GetAxis(3)->SetName("Phi_star");

  std::string xyz[4] = {"X", "Y", "Z", "all"};
  for (int i = 0; i < 4; i++) {
    delta_p[i] = std::make_unique<TH1D>(Form("dPvsP_%s", xyz[i].c_str()),
                                        Form("#DeltaP/P_{rec} vs P_{%s}", xyz[i].c_str()), 500, -0.5, 0.5);
    delta_p_electron[i] =
        std::make_unique<TH1D>(Form("dPvsP_electron_%s", xyz[i].c_str()),
                               Form("Electron #DeltaP/P_{rec} vs P_{%s}", xyz[i].c_str()), 500, -0.5, 0.5);
  }

  for (int y = 0; y < Q2_BINS; y++) {
    W_binned_MC[y] = std::make_unique<TH1D>(
        Form("W_MC_%0.3f_%0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1))),
        Form("W hist from true MC\nQ^{2} %0.3f %0.3f", q2_binned_min + (Q2_width * y),
             q2_binned_min + (Q2_width * (y + 1))),
        BINS, w_binned_min, w_binned_max);
  }
}

void mcHistogram::Fill_WQ2_MC(const std::shared_ptr<MCReaction> &_e) {
  WvsQ2_MC->Fill(_e->W_thrown(), _e->Q2_thrown());
  W_MC->Fill(_e->W_thrown());
  WvsQ2_binned_MC->Fill(_e->W_thrown(), _e->Q2_thrown());
  for (int y = 0; y < Q2_BINS; y++) {
    if (q2_binned_min + (Q2_width * y) <= _e->Q2_thrown() && q2_binned_min + (Q2_width * (y + 1)) >= _e->Q2_thrown()) {
      W_binned_MC[y]->Fill(_e->W_thrown());
      continue;
    }
  }
}

void mcHistogram::Fill(const std::shared_ptr<MCReaction> &event) {
  // std::lock_guard<std::mutex> lk(mutex);
  // std::array<double, DIMENSIONS> to_fill = {event->W_thrown(), event->Q2_thrown(), cos(event->Theta_star()),
  //                                           event->Phi_star() * RAD2DEG};
  // ndhist_mc->Fill(to_fill.data());
}

void mcHistogram::Fill_P(const std::shared_ptr<Branches> &d) {
  if (d->gpart() > 0) {
    delta_p_electron[0]->Fill((d->px(0) - d->pxpart(0)) / d->px(0));
    delta_p_electron[1]->Fill((d->py(0) - d->pypart(0)) / d->py(0));
    delta_p_electron[2]->Fill((d->pz(0) - d->pzpart(0)) / d->pz(0));
    delta_p_electron[3]->Fill((d->p(0) - TMath::Sqrt(d->pxpart(0) * d->pxpart(0) + d->pypart(0) * d->pypart(0) +
                                                     d->pzpart(0) * d->pzpart(0))) /
                              d->p(0));
    delta_px_py_electron->Fill(d->px(0) - d->pxpart(0), d->py(0) - d->pypart(0));
  }

  for (int part_num = 1; part_num < d->gpart(); part_num++) {
    delta_p[0]->Fill((d->px(part_num) - d->pxpart(part_num)) / d->px(part_num));
    delta_p[1]->Fill((d->py(part_num) - d->pypart(part_num)) / d->py(part_num));
    delta_p[2]->Fill((d->pz(part_num) - d->pzpart(part_num)) / d->pz(part_num));
    delta_p[3]->Fill((d->p(part_num) - TMath::Sqrt(d->pxpart(part_num) * d->pxpart(part_num) +
                                                   d->pypart(part_num) * d->pypart(part_num) +
                                                   d->pzpart(part_num) * d->pzpart(part_num))) /
                     d->p(part_num));
  }
}

void mcHistogram::WvsQ2_MC_Write() {
  WvsQ2_MC->SetXTitle("W (GeV/c^{2})");
  WvsQ2_MC->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_MC->SetOption("COLZ");
  WvsQ2_MC->Write();

  W_MC->SetXTitle("W (GeV/c^{2})");
  W_MC->Write();
}

void mcHistogram::Write_DeltaP() {
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  for (auto &dp : delta_p_electron) {
    dp->Fit("gaus", "QM+", "", -0.1, 0.1);
    dp->SetXTitle("#Delta P (GeV)");
    dp->Write();
  }
  auto f2 = std::make_unique<TF2>("f2", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", -0.1, 0.1, -0.1, 0.1);
  f2->SetParameters(10, 0.0, 0.1, 0.0, 0.1);
  delta_px_py_electron->Fit(f2.get(), "MR+", "");
  delta_px_py_electron->SetOption("SURF3");
  delta_px_py_electron->SetXTitle("Px");
  delta_px_py_electron->SetYTitle("Py");
  delta_px_py_electron->Write();

  auto dp_canvas = std::make_unique<TCanvas>("dp_canvas", "#Delta P", 1280, 720);
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
