{
//
//Macro to place the mass cut histograms from AnaTwoPion.root to a single canvas and create gifs of these canvases.
//

   using namespace std;
   
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,1280,960);
   c1->Divide(3,2);
   
   TCanvas *c2 = new TCanvas("c2", "c2",0,0,1280,960);
   c2->Divide(3,2);
   
   TCanvas *c3 = new TCanvas("c3", "c3",0,0,1280,960);
   c3->Divide(3,2);
   
   TCanvas *c4 = new TCanvas("c4", "c4",0,0,1280,960);
   c4->Divide(3,2);
   
   TCanvas *c5 = new TCanvas("c5", "c5",0,0,1280,960);
   c5->Divide(3,2);
 
   TCanvas *c6 = new TCanvas("c6", "c6",0,0,1280,960);
   c6->Divide(3,2);
  
   TCanvas *c7 = new TCanvas("c7", "c7",0,0,1280,960);
   c7->Divide(3,2);
   
   TCanvas *c8 = new TCanvas("c8", "c8",0,0,1280,960);
   c8->Divide(3,2);
   
   TCanvas *c9 = new TCanvas("c9", "c9",0,0,1280,960);
   c9->Divide(3,2);
   
   TCanvas *c10 = new TCanvas("c10", "c10",0,0,1280,960);
   c10->Divide(3,2);
   
   TCanvas *c11 = new TCanvas("c11", "c11",0,0,1280,960);
   c11->Divide(3,2);
   
   TCanvas *c12 = new TCanvas("c12", "c12",0,0,1280,960);
   c12->Divide(3,2);
   
   TCanvas *c13 = new TCanvas("c13", "c13",0,0,1280,960);
   c13->Divide(3,2);
   
   TCanvas *c14 = new TCanvas("c14", "c14",0,0,1280,960);
   c14->Divide(3,2);
   
   TCanvas *c15 = new TCanvas("c15", "c15",0,0,1280,960);
   c15->Divide(3,2);
   
   TCanvas *c16 = new TCanvas("c16", "c16",0,0,1280,960);
   c16->Divide(3,2);
   
   TCanvas *c17 = new TCanvas("c17","c17",0,0,1280,960);
   c17->Divide(2,2);
   
   TCanvas *c18 = new TCanvas("c18","c18",0,0,1280,960);
   c18->Divide(2,2);
   
   TCanvas *c19 = new TCanvas("c19", "c19",0,0,1280,960);
   c19->Divide(3,2);
   
   TCanvas *c20 = new TCanvas("c20", "c20",0,0,1280,960);
   c20->Divide(3,2);
   
   TFile *infile = new TFile("AnaTwoPion.root", "READ");

   Int_t i;
   

   char Hist[20];
   TH1D *h1D;
 
   for(i=1; i<=6; i++){
      c1->cd(i);
	 sprintf(Hist,"MassCutAfter%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c1->SaveAs("results/MassCutAfter.gif");
   c1->Close();
     
   for(i=1; i<=6; i++){
      c2->cd(i);
	 sprintf(Hist,"MassCutAround%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c2->SaveAs("results/MassCutAround.gif");
   c2->Close();
   
   for(i=1; i<=6; i++){
      c3->cd(i);
	 sprintf(Hist,"MassCutAfterTgt_2H_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c3->SaveAs("results/MassCutAfter2H.gif");
   c3->Close();
   
   for(i=1; i<=6; i++){
      c4->cd(i);
	 sprintf(Hist,"MassCutAfterTgt_C_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c4->SaveAs("results/MassCutAfterC.gif");
   c4->Close();
   
   for(i=1; i<=6; i++){
      c5->cd(i);
	 sprintf(Hist,"MassCutAfterTgt_FeTi_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c5->SaveAs("results/MassCutAfterFeTi.gif");
   c5->Close();
  
   for(i=1; i<=6; i++){
      c6->cd(i);
	 sprintf(Hist,"MassCutAfterTgt_Pb_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c6->SaveAs("results/MassCutAfterPb.gif");
   c6->Close();
  
   for(i=1; i<=6; i++){
      c7->cd(i);
	 sprintf(Hist,"MassCutAroundTgt_2H_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c7->SaveAs("results/MassCutAround2H.gif");
   c7->Close();
   
   for(i=1; i<=6; i++){
      c8->cd(i);
	 sprintf(Hist,"MassCutAroundTgt_C_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c8->SaveAs("results/MassCutAroundC.gif");
   c8->Close();
   
   
   for(i=1; i<=6; i++){
      c9->cd(i);
	 sprintf(Hist,"MassCutAroundTgt_FeTi_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c9->SaveAs("results/MassCutAroundFeTi.gif");
   c9->Close();
   
   for(i=1; i<=6; i++){
      c10->cd(i);
	 sprintf(Hist,"MassCutAroundTgt_Pb_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c10->SaveAs("results/MassCutAroundPb.gif");
   c10->Close();
   
   for(i=1; i<=6; i++){
      c11->cd(i);
	 sprintf(Hist,"ErrorMassCutAroundSideband_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c11->SaveAs("results/ErrorMassCutAroundSideband.gif");
   c11->Close();

   for(i=1; i<=6; i++){
      c12->cd(i);
	 sprintf(Hist,"ErrorMassCutAfterSideband_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c12->SaveAs("results/ErrorMassCutAfterSideband.gif");
   c12->Close();


   for(i=1; i<=6; i++){
      c13->cd(i);
	 sprintf(Hist,"SidebandMassCutAfter_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c13->SaveAs("results/SidebandMassCutAfter.gif");
   c13->Close();
   
   for(i=1; i<=6; i++){
      c14->cd(i);
	 sprintf(Hist,"SidebandMassCutAround_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c14->SaveAs("results/SidebandMassCutAround.gif");
   c14->Close();
      
   for(i=1; i<=6; i++){
      c15->cd(i);
	 sprintf(Hist,"YeildMassCutAroundSideband_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c15->SaveAs("results/YeildMassCutAroundSideband.gif");
   c15->Close();
   
   for(i=1; i<=6; i++){
      c16->cd(i);
	 sprintf(Hist,"YeildMassCutAfterSideband_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c16->SaveAs("results/YeildMassCutAfterSideband.gif");
   c16->Close();
   
   c17->cd(1);
   h1D = (TH1D*) infile->Get("MassCutAfterTgt_2H_0");
   h1D->Draw();
   c17->cd(2);
   h1D = (TH1D*) infile->Get("MassCutAfterTgt_C_0");
   h1D->Draw();
   c17->cd(3);
   h1D = (TH1D*) infile->Get("MassCutAfterTgt_FeTi_0");
   h1D->Draw();
   c17->cd(4);
   h1D = (TH1D*) infile->Get("MassCutAfterTgt_Pb_0");
   h1D->Draw();
   
   c17->SaveAs("results/MassCutAfterSideband.gif");
   c17->Close();
   
   c18->cd(1);
   h1D = (TH1D*) infile->Get("MassCutAroundTgt_2H_0");
   h1D->Draw();
   c18->cd(2);
   h1D = (TH1D*) infile->Get("MassCutAroundTgt_C_0");
   h1D->Draw();
   c18->cd(3);
   h1D = (TH1D*) infile->Get("MassCutAroundTgt_FeTi_0");
   h1D->Draw();
   c18->cd(4);
   h1D = (TH1D*) infile->Get("MassCutAroundTgt_Pb_0");
   h1D->Draw();
   
   c18->SaveAs("results/MassCutAroundSideband.gif");
   c18->Close();
   
   for(i=1; i<=6; i++){
      c19->cd(i);
	 sprintf(Hist,"RatioSidebandMassCutAround_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c19->SaveAs("results/RatioMassCutAroundSideband.gif");
   c19->Close();
   
   for(i=1; i<=6; i++){
      c20->cd(i);
	 sprintf(Hist,"RatioSidebandMassCutAfter_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->Draw();
   }
   c20->SaveAs("results/RatioMassCutAfterSideband.gif");
   c20->Close();
   

   for(i=1; i<=6; i++){
      c19->cd(i);
	 sprintf(Hist,"RatioSidebandMassCutAround_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->SetMarkerStyle(20);
      h1D->SetMarkerColor(i);
	 h1D->Draw("P");
	 h1D->Draw("SAME");
   }
   c19->SaveAs("results/RatioMassCutAroundSideband.gif");
   c19->Close();
   
   for(i=1; i<=6; i++){
      c20->cd(i);
      sprintf(Hist,"RatioSidebandMassCutAfter_%i",i);
      h1D = (TH1D*)infile->Get(Hist);
      h1D->SetMarkerStyle(20);
      h1D->SetMarkerColor(i);
      h1D->Draw("P");
      h1D->Draw("SAME");
   }
   c20->SaveAs("results/RatioMassCutAfterSideband.gif");
   c20->Close();
    
   
   
}
