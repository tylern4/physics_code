// Show the slice of a TH2 following the mouse position

static constexpr int W_bins = 20;
double w_binned_min = 1.0;
double w_binned_max = 3.0;
double W_width = (w_binned_max - w_binned_min) / (double)W_bins;
TH1D *Missing_Mass_WBinned[W_bins];

void WvsMM(string name = "v2_all.root") {
  // Create a new canvas.
  TCanvas *c1 = new TCanvas("c1", "Dynamic Slice Example", 10, 10, 700, 500);
  TFile *f = new TFile(name.c_str());

  char hname[80];
  for (int x = 0; x < W_bins; x++) {
    sprintf(hname, "MM_binned/MM_W_%0.3f_%0.3f", w_binned_min + (W_width * x),
            w_binned_min + (W_width * (x + 1)));
    Missing_Mass_WBinned[x] = (TH1D *)f->Get(hname);
  }

  // c1->SetFillColor(42);
  // create a 2-d histogram, fill and draw it
  TH2D *wq2 = (TH2D *)f->Get("W vs Q2/WvsQ2_channel");
  wq2->Draw("col");

  // Add a TExec object to the canvas
  c1->AddExec("dynamic", "DynamicExec()");
}

void DynamicExec() {
  TObject *select = gPad->GetSelected();
  if (!select)
    return;
  if (!select->InheritsFrom(TH2::Class())) {
    gPad->SetUniqueID(0);
    return;
  }
  TH2 *h = (TH2 *)select;
  gPad->GetCanvas()->FeedbackMode(kTRUE);

  // erase old position and draw a line at current position
  int pyold = gPad->GetUniqueID();
  int px = gPad->GetEventX();
  int py = gPad->GetEventY();
  float uxmin = gPad->GetUxmin();
  float uxmax = gPad->GetUxmax();
  int pxmin = gPad->XtoAbsPixel(uxmin);
  int pxmax = gPad->XtoAbsPixel(uxmax);
  if (pyold)
    gVirtualX->DrawLine(pxmin, pyold, pxmax, pyold);
  gVirtualX->DrawLine(pxmin, py, pxmax, py);
  gPad->SetUniqueID(py);

  Float_t upy = gPad->AbsPixeltoY(py);
  Float_t y = gPad->PadtoY(upy);

  Float_t upx = gPad->AbsPixeltoX(px);
  Float_t W = gPad->PadtoX(upx);

  // create or set the new canvas c2
  TVirtualPad *padsav = gPad;
  TCanvas *c2 = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("c2");
  if (c2)
    delete c2->GetPrimitive("Projection");
  else
    c2 = new TCanvas("c2", "Projection Canvas", 710, 10, 700, 500);
  c2->SetGrid();
  c2->cd();

  // draw slice corresponding to mouse position
  char title[80];
  for (int i = 0; i < W_bins; i++) {
    if (w_binned_min + (W_width * i) <= W &&
        w_binned_min + (W_width * (i + 1)) >= W) {
      sprintf(title, "MM for W = %0.3f -> %0.3f", w_binned_min + (W_width * i),
              w_binned_min + (W_width * (i + 1)));
      Missing_Mass_WBinned[i]->SetFillColor(38);
      Missing_Mass_WBinned[i]->SetTitle(title);
      Missing_Mass_WBinned[i]->Draw();
      continue;
    }
  }

  c2->Update();
  padsav->cd();
}
