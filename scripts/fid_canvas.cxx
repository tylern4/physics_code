#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"

TCanvas *fid_canvas(std::string infile = "v2_all.root") {
  TFile *f = new TFile(infile.c_str());

  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);

  c1->Divide(2, 3);

  TH2D *h1 = (TH2D *)f->Get("Fid_cuts/electron_fid_sec1");

  TH1D *electron_fid_sec_slice;
  double slice_width = (500 / 200);
  char *hname;
  /*
    for (int slice = 0; slice < 200; slice++) {
      sprintf(hname, "electron_fid_slice_%d", slice + 1);
      electron_fid_sec_slice[slice] =
          h1->ProjectionX(hname, slice_width * slice, slice_width * slice + (slice_width - 1));
      electron_fid_sec_slice[slice]->Rebin(10);
      electron_fid_sec_slice[slice]->Draw();
    }
  */
  c1->cd(2);
  electron_fid_sec_slice = h1->ProjectionX("name", 200, 300);
  electron_fid_sec_slice->Rebin(10);
  electron_fid_sec_slice->Draw();

  c1->cd(1);
  h1->Draw();

  TH2D *h2 = (TH2D *)f->Get("Fid_cuts/electron_fid_sec2");
  c1->cd(2);
  h2->Draw();

  TH2D *h3 = (TH2D *)f->Get("Fid_cuts/electron_fid_sec3");
  c1->cd(3);
  h3->Draw();

  TH2D *h4 = (TH2D *)f->Get("Fid_cuts/electron_fid_sec4");
  c1->cd(4);
  h4->Draw();

  TH2D *h5 = (TH2D *)f->Get("Fid_cuts/electron_fid_sec5");
  c1->cd(5);
  h5->Draw();

  TH2D *h6 = (TH2D *)f->Get("Fid_cuts/electron_fid_sec6");
  c1->cd(6);
  h6->Draw();

  return c1;
}

#if !defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tfid_canvas data.root" << std::endl;
    exit(1);
  }

  auto can = fid_canvas(argv[1]);
  can->SaveAs("fid_canvas.pdf");
  return 0;
}
#endif