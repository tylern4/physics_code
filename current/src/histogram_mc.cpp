/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram_mc.hpp"

mcHistogram::mcHistogram(std::string output_file) {
  RootOutputFile = new TFile(output_file.c_str(), "RECREATE");
  def = new TCanvas("def");
  makeHists();
}

mcHistogram::~mcHistogram() {
  std::cerr << GREEN << "\nWriting" << DEF << std::endl;
  RootOutputFile->cd();
  // Start of cuts
  Fits *MM_neutron_cut = new Fits();
  MM_neutron_cut->Set_min(0.8);
  MM_neutron_cut->Set_max(1.2);
  MM_neutron_cut->FitBreitWigner(Missing_Mass);

  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  WvsQ2_folder->cd();
  WvsQ2_Write();
  TDirectory *W_Q2_binned = RootOutputFile->mkdir("W_Q2_binned");
  W_Q2_binned->cd();
  WvsQ2_binned_Write();
  TDirectory *MM_folder = RootOutputFile->mkdir("MM_folder");
  MM_folder->cd();
  Write_Missing_Mass();
  TDirectory *delta_mom = RootOutputFile->mkdir("delta_mom");
  delta_mom->cd();
  Write_DeltaP();
  RootOutputFile->Close();
  std::cerr << BOLDBLUE << "Done!!!" << DEF << std::endl;
}

// W and Q^2
void mcHistogram::makeHists() {
  std::string xyz[4] = {"X", "Y", "Z", "all"};
  for (int i = 0; i < 4; i++) {
    sprintf(hname, "dPvsP_%s", xyz[i].c_str());
    sprintf(htitle, "#DeltaP/P_{rec} vs P_{%s}", xyz[i].c_str());
    delta_p[i] = std::make_unique<TH1D>(hname, htitle, 500, -0.5, 0.5);
  }

  for (int y = 0; y < Q2_bins; y++) {
    sprintf(hname, "W_%0.3f_%0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    sprintf(htitle, "W hist\nQ^{2} %0.3f %0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    W_binned[y] = std::make_unique<TH1D>(hname, htitle, bins, w_binned_min, w_binned_max);

    sprintf(hname, "W_MC_%0.3f_%0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    sprintf(htitle, "W hist from true MC\nQ^{2} %0.3f %0.3f", q2_binned_min + (Q2_width * y),
            q2_binned_min + (Q2_width * (y + 1)));
    W_binned_MC[y] = std::make_unique<TH1D>(hname, htitle, bins, w_binned_min, w_binned_max);
  }
}

void mcHistogram::Fill_WQ2(double W, double Q2) {
  WvsQ2_hist->Fill(W, Q2);
  W_hist->Fill(W);
  WvsQ2_binned->Fill(W, Q2);
  for (int y = 0; y < Q2_bins; y++) {
    if (q2_binned_min + (Q2_width * y) <= Q2 && q2_binned_min + (Q2_width * (y + 1)) >= Q2) {
      W_binned[y]->Fill(W);
      continue;
    }
  }
}

void mcHistogram::Fill_WQ2_MC(double W, double Q2) {
  WvsQ2_MC->Fill(W, Q2);
  W_MC->Fill(W);
  WvsQ2_binned_MC->Fill(W, Q2);
  for (int y = 0; y < Q2_bins; y++) {
    if (q2_binned_min + (Q2_width * y) <= Q2 && q2_binned_min + (Q2_width * (y + 1)) >= Q2) {
      W_binned_MC[y]->Fill(W);
      continue;
    }
  }
}

void mcHistogram::Fill_P(std::shared_ptr<Branches> d) {
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

void mcHistogram::WvsQ2_Write() {
  WvsQ2_hist->SetXTitle("W (GeV)");
  WvsQ2_hist->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_hist->SetOption("COLZ");
  WvsQ2_hist->Write();

  W_hist->SetXTitle("W (GeV)");
  W_hist->Write();

  WvsQ2_MC->SetXTitle("W (GeV)");
  WvsQ2_MC->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_MC->SetOption("COLZ");
  WvsQ2_MC->Write();

  W_MC->SetXTitle("W (GeV)");
  W_MC->Write();
}

void mcHistogram::WvsQ2_binned_Write() {
  WvsQ2_binned->SetXTitle("W (GeV)");
  WvsQ2_binned->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_binned->SetOption("COLZ");
  WvsQ2_binned->Write();

  WvsQ2_binned_MC->SetXTitle("W (GeV)");
  WvsQ2_binned_MC->SetYTitle("Q^{2} (GeV^{2})");
  WvsQ2_binned_MC->SetOption("COLZ");
  WvsQ2_binned_MC->Write();

  for (int x = 0; x < Q2_bins; x++) {
    W_binned[x]->SetXTitle("W (GeV)");
    W_binned[x]->Write();
  }
  for (int x = 0; x < Q2_bins; x++) {
    W_binned_MC[x]->SetXTitle("W (GeV)");
    W_binned_MC[x]->Write();
  }
}
// W and Q^2

// Missing Mass
void mcHistogram::Fill_Missing_Mass(double mm, double mm2) {
  Missing_Mass->Fill(mm);
  Missing_Mass_square->Fill(mm2);
}

void mcHistogram::Write_Missing_Mass() {
  Missing_Mass->SetXTitle("Mass (GeV)");
  Missing_Mass->Write();

  Missing_Mass_square->SetXTitle("Mass^{2} (GeV^{2})");
  Missing_Mass_square->Write();
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

  for (int i = 0; i < 4; i++) delta_p[i]->Write();
}
