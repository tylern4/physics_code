{
  TFile *f = new TFile("v2_all.root");

  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
  c1->SetLogz();
  TH2D *h1 = (TH2D *)f->Get("Delta_T/delta_t_mass_PIP");
  gStyle->SetPalette(kCividis);
  // TColor::InvertPalette();
  h1->GetXaxis()->SetRange(100, 200);
  h1->Draw();
}
