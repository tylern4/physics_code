{
   gStyle->SetPalette(1);
   
   using namespace std;
   
   Int_t i;   
   char Hist[20];
   TH1D *h1D;  
    
   TFile *infile = new TFile("AnaTwoPion.root", "READ");
   
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,1280,960);


   for(i=1; i<=6; i++){

      c1->cd();      
	 sprintf(Hist,"MassCutAfter_%i",i);

      h1D = (TH1D*)infile->Get(Hist);
      //h1D->
      h1D->Draw();
      }
      
      
   c1->SaveAs("results/Targets.gif");
   c1->Close();









}
