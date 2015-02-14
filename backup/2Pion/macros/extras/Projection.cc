//
//   Projection 
//
//
//   Makes either an X or Y projection histogram from 2D histograms.
//


#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH2.h"
#include <TLorentzVector.h>
#include "TFile.h"
#include <iostream>
#include "TObject.h"
#include <stdio.h>



using namespace std;


int main(int argc, char **argv){

	float time1 = clock();
	
	
  	TH1D *h1D;
  	TH2F *h2D;  	
  	  	
	//char OutCan = "title";
	Int_t binLo = 0;
	Int_t binHi = -1;
	  
  	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
  	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  	c1->SetBorderSize(5); 
  	gStyle->SetOptStat(0);
  	gStyle->SetPalette(1);
  	c1->SetFillStyle(4000);

  	char strname[50];

     char infilename[128];
     char outfilename[128];
    	char XorY[128];


	if(argc < 4){
             cout<<"Error : Too Few Arguments" <<endl;
             cout<<"To run: ./Test [Infilename] [Outfilename] [Axis]"<<endl;
             exit(0);}


	sprintf(outfilename,"%s",argv[2]);
	sprintf(infilename, "%s",argv[1]);
	sprintf(XorY, "%s", argv[3]);

	TFile *myFile = new TFile(outfilename,"RECREATE");
	TFile *infile = new TFile(infilename,"READ"); 
	
	
	
  	c1->cd();
  	gPad->SetLeftMargin(1.5);
  	gPad->SetRightMargin(1.5);
  	gPad->SetFillColor(0);


	cout<<"Analyzing file "<<infilename<<endl;

	if (!strcmp(XorY, "X")or!strcmp(XorY,"x") or !strcmp(XorY, "b")){

		h2D = (TH2F*)infile->Get("hIM_P0_0");
		sprintf(strname,"Pi+MassAround");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
 		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi^{+} Mass Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P0_1");
		sprintf(strname,"Pi-MassAround");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi^{-} Mass Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P0_2");
		sprintf(strname,"PionPairMassAround");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi Pair Mass Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P1_0");
		sprintf(strname,"Pi+MassAfter");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi^{+} Mass After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P1_1");
		sprintf(strname,"Pi-MassAfter");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi^{-} Mass After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P1_2");
		sprintf(strname,"PionPairMassAfter");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi pair Mass After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();}


	if (!strcmp(XorY, "Y")or!strcmp(XorY,"y")or !strcmp(XorY, "b")){
		
		h2D = (TH2F*)infile->Get("hIM_P0_0");
		sprintf(strname,"Pi+MomentumAround");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi^{+} Momentum Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();

			
		h2D = (TH2F*)infile->Get("hIM_P0_1");
		sprintf(strname,"Pi-MomentumAround");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi^{-} Momentum Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P0_2");
		sprintf(strname,"PionPairMomentumAround");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi Pair Momentum Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P1_0");
		sprintf(strname,"Pi+MomentumAfter");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi^{+} Momentum After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P1_1");
		sprintf(strname,"Pi-MomentumAfter");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi^{-} Momentum After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();			
			
		h2D = (TH2F*)infile->Get("hIM_P1_2");
		sprintf(strname,"PionMomentumAfter");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi Mometum After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();}


	/*else if{
		cout<<"No axis "<<XorY<<endl;
		cout<<"Try x or y"<<endl;
		exit(0);}*/

	
	myFile->Write(outfilename);
	myFile->Close(outfilename);


	//time to complete function
  	float minutes = 0;
	float seconds = 0;
  	float time2 = clock();
  	minutes = (time2 - time1)/1000000;
  	minutes = (minutes)/60;
  	seconds = fmod(minutes,1);
	minutes = minutes-seconds;
	seconds = seconds*60;


  	if (minutes==0) cout<<endl<<"Completed in "<<seconds<<" seconds."<<endl<<endl;
  	else cout<<endl<<"Completed in "<<minutes<<" minutes and "<<seconds<<" seconds."<<endl<<endl;


	return 0;


}
