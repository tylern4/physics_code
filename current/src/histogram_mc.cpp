/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram_mc.hpp"

mcHistogram::mcHistogram(std::string output_file) {
  RootOutputFile = new TFile(output_file.c_str(), "RECREATE");
  def = new TCanvas("def");
  makeHists_W();
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

  RootOutputFile->Close();
  std::cerr << BOLDBLUE << "Done!!!" << DEF << std::endl;
}

// W and Q^2
void mcHistogram::makeHists_W() {
  for (int y = 0; y < Q2_bins; y++) {
    sprintf(hname, "W_%0.3f_%0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    sprintf(htitle, "W hist\nQ^{2} %0.3f %0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    W_binned[y] = new TH1D(hname, htitle, bins, w_binned_min, w_binned_max);

    sprintf(hname, "W_MC_%0.3f_%0.3f", q2_binned_min + (Q2_width * y), q2_binned_min + (Q2_width * (y + 1)));
    sprintf(htitle, "W hist from true MC\nQ^{2} %0.3f %0.3f", q2_binned_min + (Q2_width * y),
            q2_binned_min + (Q2_width * (y + 1)));
    W_binned_MC[y] = new TH1D(hname, htitle, bins, w_binned_min, w_binned_max);
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
void mcHistogram::Fill_Missing_Mass(MissingMass *miss_mass) {
  Missing_Mass->Fill(miss_mass->Get_MM());
  Missing_Mass_square->Fill(miss_mass->Get_MM2());
}

void mcHistogram::Write_Missing_Mass() {
  Missing_Mass->SetXTitle("Mass (GeV)");
  Missing_Mass->Write();

  Missing_Mass_square->SetXTitle("Mass^{2} (GeV^{2})");
  Missing_Mass_square->Write();
}
