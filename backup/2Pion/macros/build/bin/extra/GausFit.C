{

   gStyle->SetPalette(1);
   gStyle->SetOptFit(1);
   using namespace std;
   Int_t i;

   TCanvas *c1 = new TCanvas("c1", "c1",0,0,1280,960);
   TCanvas *c2 = new TCanvas("c2", "c2",0,0,1280,960);
   TCanvas *c3 = new TCanvas("c3", "c3",0,0,1280,960);
   TCanvas *c4 = new TCanvas("c4", "c4",0,0,1280,960);
   TCanvas *c5 = new TCanvas("c5", "c5",0,0,1280,960);
   TCanvas *c6 = new TCanvas("c6", "c6",0,0,1280,960);


   char Hist[20];
   TH1D *h1D;
   TFile *infile = new TFile("AnaTwoPion.root", "READ");
   TF1 *f1 = new TF1("f1","gaus",0.489,0.505);

      sprintf(Hist,"MassCutAroundTgt_C_0");
      c1->cd();
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Fit("f1","R");
      h1D->Draw();
      c1->SaveAs("gaus/MassCutAroundCGaus.gif");
      c1->Close();

      sprintf(Hist,"MassCutAroundTgt_FeTi_0");
      c2->cd();
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Fit("f1","R");
      h1D->Draw();
      c2->SaveAs("gaus/MassCutAroundFeTiGaus.gif");
      c2->Close();

      sprintf(Hist,"MassCutAroundTgt_Pb_0");
      c3->cd();
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Fit("f1","R");
      h1D->Draw();
      c3->SaveAs("gaus/MassCutAroundPbGaus.gif");
      c3->Close();

      sprintf(Hist,"MassCutAfterTgt_C_0");
      c4->cd();
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Fit("f1","R");
      h1D->Draw();
      c4->SaveAs("gaus/MassCutAfterCGaus.gif");
      c4->Close();

      sprintf(Hist,"MassCutAfterTgt_FeTi_0");
      c5->cd();
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Fit("f1","R");
      h1D->Draw();
      c5->SaveAs("gaus/MassCutAfterFeTiGaus.gif");
      c5->Close();

      sprintf(Hist,"MassCutAfterTgt_Pb_0");
      c6->cd();
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Fit("f1","R");
      h1D->Draw();
      c6->SaveAs("gaus/MassCutAfterPbGaus.gif");
      c6->Close();



}

