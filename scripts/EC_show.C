{
  TCanvas c1;

  TFile in_file("v2_all.root");
  TH2D* h;
  in_file.GetObject("EC_hists/EC_sampling_fraction", h);

  TF1* fit1;
  in_file.GetObject("EC_hists/EC_P_fit", fit1);

  TF1* fit2;
  in_file.GetObject("EC_hists/EC_M_fit", fit2);

  TGraph* P;
  in_file.GetObject("EC_hists/Positive_EC_graph", P);

  TGraph* M;
  in_file.GetObject("EC_hists/Negative_EC_graph", M);

  h->Draw();

  fit1->Draw("same");
  fit2->Draw("same");
  P->Draw("*same");
  M->Draw("*same");
}
