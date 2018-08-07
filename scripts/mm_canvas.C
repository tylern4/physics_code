{
  TFile *f = new TFile("v2_all.root");

  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);

  TH2D *h1 = (TH2D *)f->Get("Missing_Mass/Missing_Mass");

  double l = 0.8;
  double h = 1.1;

  TLine *lower = new TLine(l, 0, l, 180000);
  TLine *upper = new TLine(h, 0, h, 180000);
  lower->SetLineColor(46);
  upper->SetLineColor(46);

  lower->SetLineWidth(2);
  upper->SetLineWidth(2);

  h1->Draw();
  lower->Draw("same");
  upper->Draw("same");
}
