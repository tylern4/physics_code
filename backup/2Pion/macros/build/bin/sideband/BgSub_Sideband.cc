// 
// BgSub_Sideband - background subtraction routine. The sidebands around the peak form the
//                  background shape and are subtracted the yield under the peak.
//                  The output is a list of yields.
//                  
//                  fname = ROOT file containing the histogram
//                  hname = histogram name
//					yieldFile = name of output file without the suffix
//
// M. H. Wood, Canisius College
//--------------------------------------------------------------------------------------------
#include "TTree.h"
#include "TROOT.h"
#include "TH2.h"
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include "TStyle.h"
#include "TCanvas.h"
#include "TObject.h"
#include <stdio.h>



void BgSub_Sideband(char *rootFile, char *hname, char* yieldFile)
{
     gROOT->Reset();

	Float_t Lmar = 0.125;
	Float_t Rmar = 0.125;
	Float_t yoff = 1.5;
	
	Float_t xLo = 0.35;  // lower value of x-axis for drawing histogram
	Float_t xHi = 0.6;   // upper value of x-axis for drawing histogram
	Float_t PeakLo = 0.48; // lower limit on the peak range
	Float_t PeakHi = 0.51; // upper limit on the peak range
	Float_t Ks_width = PeakHi - PeakLo; // width of the peak range
	Float_t SidebandLo = PeakLo - 0.5*Ks_width; // lower limit of left sideband
	Float_t SidebandHi = PeakHi + 0.5*Ks_width; // upper limit of right sideband
	
	Int_t i, j, k;
	Float_t x;
	Float_t yield, yield_peak, yield_sb;
	
	TH1F *hist; // original histogram
	TH1F *hPeak; // histogram of peak range
	TH1F *hSideband; // histogram of sidebands
		
	// create canvas
	char xtitle[100];
	char ytitle[100];
	char title[100];
	sprintf(title,"Sideband Analysis"); // canvas title
	TCanvas *can1 = new TCanvas("can1",title,0,0,600,600); // create the canvas
	
	gStyle->SetOptStat(1111);
	can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
	can1->SetBorderSize(5);  
	can1->SetFillStyle(4000); 
	
	sprintf(xtitle,"Invariant Mass (GeV)"); // set the x-axis title
	sprintf(ytitle,"Counts"); // set the y-axis title
	
	// data files contain the trees
	printf("Analyzing file %s\n",rootFile);  
	TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
	hist = (TH1F*)fd->Get(hname); // get the histogram from the ROOT file
					
	// clone the original histogram for the sideband histogram
	sprintf(hname,"hSideband");
	hSideband = (TH1F*)hist->Clone(hname);    
		
	// clone the original histogram for the peak histogram
	sprintf(hname,"hPeak");
	hPeak = (TH1F*)hist->Clone(hname);    
		
	// loop to fill the peak and sideband histograms
	for(k=1;k<=hist->GetNbinsX();k++){
		x = hist->GetBinCenter(k); // x-axis value for bin k

		// fill peak histogram with counts from original histogram for bins in range
		if(x>=PeakLo && x<PeakHi){
			hPeak->SetBinContent(k,hist->GetBinContent(k));
		}else{
			hPeak->SetBinContent(k,0.0);
		}
		
		// fill sideband histogram when bin is in the proper ranges
		if(x>=SidebandLo && x<PeakLo){
			hSideband->SetBinContent(k,hist->GetBinContent(k));
		}else if(x>=PeakHi && x<SidebandHi){
			hSideband->SetBinContent(k,hist->GetBinContent(k));
		}else{
			hSideband->SetBinContent(k,0.0);
		}
	}
		
	// set up the Pad parameters
	gPad->SetLeftMargin(Lmar);
	gPad->SetRightMargin(Rmar);
	gPad->SetFillColor(0);
	
	// draw the original histogram
	hist->SetTitle(0);
	hist->GetXaxis()->SetTitle(xtitle);
	hist->GetXaxis()->CenterTitle();
	hist->GetYaxis()->SetTitle(ytitle);
	hist->GetYaxis()->CenterTitle();
	hist->GetYaxis()->SetTitleOffset(yoff);
	hist->GetXaxis()->SetRangeUser(xLo,xHi);
	hist->SetLineWidth(2);
	hist->Draw(); 
	
	// draw the peak histogram
	hPeak->SetLineWidth(2);
	hPeak->SetFillColor(2);
	hPeak->Draw("same");
	
	// draw the sideband histogram
	hSideband->SetLineWidth(2);
	hSideband->SetFillColor(4);
	hSideband->Draw("same");
	

	// create the image files
	char OutCan[100];
	sprintf(OutCan,"%s.gif",yieldFile);
	can1->Print(OutCan);
	sprintf(OutCan,"%s.eps",yieldFile);
	can1->Print(OutCan);

	// open text file for the yields
	char OutFile[100];
	sprintf(OutFile,"%s.yld",yieldFile);
	ofstream fout(OutFile); 

	// sum counts in the peak range
	yield_peak = hPeak->Integral(hPeak->FindBin(PeakLo),hPeak->FindBin(PeakHi));

    // sum counts in sidebands
	yield_sb = hSideband->Integral(hSideband->FindBin(SidebandLo),hSideband->FindBin(SidebandHi));

	yield = yield_peak - yield_sb; // subtract background from peak sum

	fout<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;
	fout.close(); // close the text file

}
