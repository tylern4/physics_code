/* 
Quick and dirty picture maker 
Change lines with *** to change file and histograms to make
use myFile->ls(); to find the hists in the file
*/

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include "TStyle.h"

using namespace std;
int main(int argc, char **argv){

	TFile *myFile;
	char  outfilename[128];
	TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);

	myFile = new TFile("../outputFiles/release_5-22_1-42.root","READ"); /***/

	//myFile->ls();  /***/
	TH2D *h1 = (TH2D*)myFile->Get("WvsQ2_hist"); /***/

	c1->cd();
	h1->Draw("color");
	for (int i = 0; i < 200; ++i)
	{
		gStyle->SetPalette(i);
		sprintf(outfilename,"%s_%d.png","WvsQ2",i);
		//cout << outfilename << endl;
		c1->Print(outfilename); /***/	
	}


	myFile->Close();

	return 0;
}
