{
  TFile *f_wFid = new TFile("v2_all.root");
  TFile *f_noFid = new TFile("v2_all_noFid.root");

  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);

  c1->Divide(0, 3);

  TH2D *h_noFid = (TH2D *)f_noFid->Get("Fid_cuts/electron_fid");
  TH2D *h_wFid = (TH2D *)f_wFid->Get("Fid_cuts/electron_fid");
  TH2D h_sub = *h_noFid;
  h_sub.Add(h_wFid, -1);

  c1->cd(1);
  h_noFid->Draw("same");
  c1->cd(2);
  h_wFid->Draw();
  c1->cd(3);
  h_sub.Draw();
}
