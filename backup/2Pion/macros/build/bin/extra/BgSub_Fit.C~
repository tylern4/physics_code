// 
// BgSub_Fit - background subtraction routine where the peak is fit with a well-known function
//			   like a gaussian.  The backgorund is fit with another well-known function like a
//             polynomial.
//             The output is a list of yields.
//                  
//                  fname = ROOT file containing the histogram
//                  hname = histogram name
//		    yieldFile = name of output file without the suffix
//
// M. H. Wood, Canisius College
//--------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();

Int_t PolNum = 10;
Double_t par[10];

void BgSub_Fit(char *rootFile, char *hname, char* yieldFile)
{


	Float_t Lmar = 0.125; // set the left margin
	Float_t Rmar = 0.125; // set the right margin
	Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values
	
	Double_t fPeak, fWidth, SumLo, SumHi;
	Float_t xLo = 0.41;  // lower value of x-axis for drawing histogram
	Float_t xHi = 0.6;   // upper value of x-axis for drawing histogram
	Float_t PeakLo = 0.49; // lower limit on the peak range
	Float_t PeakHi = 0.51; // upper limit on the peak range
	Float_t Nsigma = 3.0;  // number of sigma of the peak
	
	Int_t ibg = 3; // parameter index in the par[] array where the background parameters start

	Int_t i, k;
	Float_t x, pval;
	Float_t yield;
	
	TH1F *hist; // original histogram
	TH1F *hBgFit; // histogram of background
		
	// create canvas
	char title[100];
	char xtitle[100];
	char ytitle[100];
	sprintf(title,"Fitting Analysis"); // canvas title
	TCanvas *can1 = new TCanvas("can1",title,0,0,1280,720); // create the canvas
	
//	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(1111);
	can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
	can1->SetBorderSize(5);  
	can1->SetFillStyle(4000); 
	
	sprintf(xtitle,"Invariant Mass (GeV)"); // set the x-axis title
	sprintf(ytitle,"Counts"); // set the y-axis title
	
	// data files contain the trees
	cout<<"Analyzing file "<<rootFile<<endl;  
	TFile *fd = new TFile(rootFile,"READ"); // open up the ROOT file
	hist = (TH1F*)fd->Get(hname); // get the histogram from the ROOT file
					
	// histogram of background
	hBgFit = (TH1F*)hist->Clone("hBgFit"); // clone original hist. into background temp. hist.
	hBgFit->SetName("hBgFit");
	hBgFit->SetTitle("Background");

	// fit the histogram

	for(Int_t i=0; i<=PolNum; i++){par[i] = 1.0;}

    TF1 *g1 = new TF1("g1",gaussFit,PeakLo,PeakHi,3); // declare fit fcn
    TF1 *pol = new TF1("pol",polFit,xLo,xHi,PolNum-3);
    TF1 *t1 = new TF1("t1",totFit,xLo,xHi,PolNum);
		
	g1->SetParameters(&par[0]);    // set parameters for initial peak fit
	hist->Fit("g1","R");           // fit the peak
	g1->GetParameters(&par[0]);    // get parameters from initial peak fit

	hist->Fit("pol","R+");         // fit the background
	pol->GetParameters(&par[ibg]); // get parameters fromt background fit

	t1->SetParameters(par);       // set the parameters from initial fits for total fit
	t1->SetLineWidth(5);          // make fit line thicker
	t1->SetLineColor(1);          // set the fit line color
	hist->Fit("t1","R");          // fit spectrum with total fit function
	t1->GetParameters(&par[0]);   //get final parameters

	SumLo = 0.486;
	SumHi = 0.507;


			
	pol->SetParameters(&par[ibg]); // set the pfinal parameters for the background function
	hBgFit->Add(pol,-1.0); // subtract background function
	for(k=1; k<hBgFit->GetNbinsX(); k++){
		x = hBgFit->GetBinCenter(k); // read bin center value on the x-axis
		if(x>=SumLo && x<=SumHi){ // check that x is in the peak summation region
				pval = hBgFit->GetBinContent(k); // get the number of counts in the bin
			if(pval<0.0) pval = 0.0; // if neg. counts, set to zero
		}else{
			pval =0.0;
		}
		hBgFit->SetBinContent(k,pval); // refill the histogram
	}
	
	yield = hBgFit->Integral(hBgFit->FindBin(SumLo),hBgFit->FindBin(SumHi));//sum total counts in peak

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
	hist->GetXaxis()->SetRangeUser(xLo,xHi);  // set the x-axis range for the plot
	hist->SetLineWidth(2);
	hist->SetMinimum(0); // start the y-axis at zero.
	hist->Draw();
	
	// draw the background-subtracted histogram
	hBgFit->SetLineWidth(2);
	hBgFit->SetFillColor(4);
	hBgFit->Draw("same");
	
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

	fout<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;
	fout.close(); // close the text file
	cout<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;

}

// peak is Gaussian
Double_t gaussFit(Double_t *x, Double_t *par){
	return TMath::Max(1.e-10,par[0]*TMath::Gaus(x[0],par[1],par[2]));
}
// peak is Breit Wigner
Double_t breitwigner(Double_t *x, Double_t *par){
	return TMath::Max(1.e-10,par[0]*TMath::BreitWigner(x[0],par[1],par[2]));
}

// background function is polynomial
Double_t polFit(Double_t *x, Double_t *par){
	Int_t nmax = PolNum-4;
	Double_t bck = 0.0;
	for (Int_t i = 0; i<=nmax; i++){
		bck += par[i]*pow(x[0],i);
	}
	return bck;
}

// Sum of background and peak function
// peak uses par[0-2]
// background uses par[3-5]
Double_t totFit(Double_t *x, Double_t *par) {
	return gaussFit(x,par) + polFit(x,&par[3]);
}

