// 
// PolSubtraction - background subtraction routine where the peak is fit with a well-known function
//			   like a gaussian.  The backgorund is fit with another well-known function like a
//             polynomial.
//             The output is a list of yields.
//                  
//                  fname = ROOT file containing the histogram
//                  hname = histogram name
//		    yieldFile = name of output file without the suffix
//
//	to run
//	.x PolSubtraction ("exmaple.root", "exmaple_histogram_name", "outfile_name")
//
//
// M. H. Wood, Canisius College
// Nick Tyler Canisius College
//--------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();
gStyle->SetPalette(1);  


Int_t PolNum = 10;
Double_t par[10];

void PolSubtractionAngle2(char *rootFile, char* yieldFile)
{
	char *tar;
	char hname[50];
	Double_t ratio[37];
	Double_t ratios;


	Float_t Lmar = 0.125; // set the left margin
	Float_t Rmar = 0.125; // set the right margin
	Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values
	
	Double_t fPeak, fWidth, SumLo, SumHi;
	Float_t xLo = 0.41;  // lower value of x-axis for drawing histogram
	Float_t xHi = 0.6;   // upper value of x-axis for drawing histogram
	Float_t PeakLo = 0.49; // lower limit on the peak range
	Float_t PeakHi = 0.51; // upper limit on the peak range
	
	Int_t ibg = 1; // parameter index in the par[] array where the background parameters start

	Int_t i, j,k;
	Float_t x, pval;
	Float_t yield;
	
	TH1F *hist; // original histogram
	TH1F *hBgFit; // histogram of background
	TH1F *Peak;  //histogram of peak
	TH1F *RatioHist[4];
		
	// create canvas
	char title[100];
	char xtitle[100];
	char ytitle[100];
	sprintf(title,"Fitting Analysis"); // canvas title
	TCanvas *can1 = new TCanvas("can1",title,0,0,1280,720); // create the canvas
	
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(1111);
	can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
	can1->SetBorderSize(5);  
	can1->SetFillStyle(4000); 
	
	sprintf(xtitle,"Invariant Mass (GeV)"); 	// set the x-axis title
	sprintf(ytitle,"Counts"); 			// set the y-axis title
	
	// data files contain the trees
	cout<<"Analyzing file "<<rootFile<<endl;  
	TFile *fd = new TFile(rootFile,"UPDATE"); 	// open up the ROOT file

	for(j=0;j<=3;j++){
	   if(j==0) tar = "2H";
	   else if(j==1) tar = "C";
	   else if(j==2) tar = "FeTi";
	   else if(j==3) tar = "Pb";
	   cout<<j<<" "<<tar<<endl;


	for(i=0;i<=36;i++){
//for(i=0;i<=3;i++){
	   sprintf(hname, "AngleCutAroundTgt_range_%d_tgt_%s",i, tar);
	   hist = (TH1F*)fd->Get(hname); 			// get the histogram from the ROOT file
					
	   // histogram of background
	   hBgFit = (TH1F*)hist->Clone("hBgFit"); // clone original hist. into background temp. hist.
	   hBgFit->SetName("hBgFit");
	   hBgFit->SetTitle("Background");

	   Peak = (TH1F*)hist->Clone("Peak");
	   Peak->SetName("Peak");
	   Peak->SetTitle("Peak");


/* 
   use SumHi and SumLo to remove the bin contents from this region from the Bghist
   find the polinomial relating to the Bghist without the peak
   use this new function to subtract from the original histogram with the peak leaving only the peak value
*/

	   SumLo = 0.486;
	   SumHi = 0.507;			
	   Int_t Xbins = hBgFit->GetNbinsX();
	   Double_t Sub;

	   for(k=1; k<Xbins; k++){
	      if (hBgFit->GetBinCenter(k)<=SumHi && hBgFit->GetBinCenter(k)>=SumLo){ 
	         hBgFit->SetBinContent(k,0.0);//Set all bins inside peak range to be 0
	      }
	   }

	   //TF1 *pol = new TF1("pol",polFit,xLo,xHi,PolNum);
TF1 *pol = new TF1("pol",totFit,xLo,xHi,PolNum);
	   hBgFit->Fit("pol","R+");         // fit the background
	   pol->GetParameters(&par[ibg]);
	

	   pol->SetParameters(&par[ibg]); 	// set the pfinal parameters for the background function
	   Peak->Add(pol,-1.0); 		// subtract background function
	   for(k=1; k<Peak->GetNbinsX(); k++){
		x = Peak->GetBinCenter(k); 		// read bin center value on the x-axis
		if(x>=SumLo && x<=SumHi){ 		// check that x is in the peak summation region
			pval = Peak->GetBinContent(k);	// get the number of counts in the bin
			if(pval<0.0) pval = 0.0; 	// if neg. counts, set to zero
		}
		else{
			pval =0.0;
		}
		Peak->SetBinContent(k,pval); // refill the histogram
	   }
	   yield = Peak->Integral(Peak->FindBin(SumLo),Peak->FindBin(SumHi)); //sum total counts in peak


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
	   hist->SetLineWidth(1);
	   hist->SetMinimum(0); // start the y-axis at zero.
	   hist->Draw();

	   // draw the polinomial fit
	   pol->Draw("same");


	   //draw peak
	   Peak->SetFillColor(4);
	   Peak->Draw("same");
	
	
	   // create the image files
	   char OutCan[100];
	   sprintf(OutCan,"pol/%s%d.gif",yieldFile,i);
	   can1->Print(OutCan);
	   sprintf(OutCan,"pol/%s%d.eps",yieldFile,i);
	   can1->Print(OutCan);

	   // open text file for the yields
	   char OutFile[100];

	   sprintf(OutFile,"pol/%s%d.yld",yieldFile,i);
	   ofstream fout(OutFile); 

           if(i == 0) ratios = yield;
	   ratio[i] = yield/ratios;

	   fout<<"\t"<<yield<<"\t"<<sqrt(yield)<<"\t"<<ratio[i]<<endl;
	   cout<<"\t"<<yield<<"\t"<<sqrt(yield)<<"\t"<<ratio[i]<<endl;

  }

	/*sprintf(title,"Ratio of %s for different angles", tar); // canvas title
	TCanvas *can2 = new TCanvas("can2",title,0,0,1280,720); // create the canvas*/


        sprintf(hname,"RatioHist_%s",tar);
        sprintf(title,"ratio around %s",tar);
        RatioHist[j] = new TH1F(hname,title, 360, 0, 180 );

	for(i=0;i<=36;i++) RatioHist[j]->SetBinContent(i*10, ratio[i]);
//for(i=0;i<=3;i++) RatioHist[j]->SetBinContent(i*10, ratio[i]);
}

	   gPad->SetLogy();
	   RatioHist[0]->GetXaxis()->SetTitle("angle in degrees");
	   RatioHist[0]->GetXaxis()->CenterTitle();
	   RatioHist[0]->SetTitle("Ratio of Yeilds to 2H");  	   
  	   RatioHist[3]->SetMarkerSize(1.25);
  	   RatioHist[3]->SetMarkerStyle(22);
	   RatioHist[3]->SetMarkerColor(5);
	   RatioHist[3]->Draw("p");
  	   RatioHist[2]->SetMarkerSize(1.25);
  	   RatioHist[2]->SetMarkerStyle(22);
	   RatioHist[2]->SetMarkerColor(4);  	   
	   RatioHist[2]->Draw("p""same");
  	   RatioHist[1]->SetMarkerSize(1.25);
  	   RatioHist[1]->SetMarkerStyle(22);
	   RatioHist[1]->SetMarkerColor(3);  	   
	   RatioHist[1]->Draw("p""same");
  	   RatioHist[0]->SetMarkerSize(1.25);
  	   RatioHist[0]->SetMarkerStyle(22);
  	   RatioHist[0]->SetMarkerColor(2);
	   RatioHist[0]->Draw("p""same");

	   sprintf(OutCan,"pol/Ratio%s.gif",yieldFile);
	   can1->Print(OutCan);
	   sprintf(OutCan,"pol/Ratio%s.eps",yieldFile);
	   can1->Print(OutCan);


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
	Int_t nmax = PolNum-1;
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

