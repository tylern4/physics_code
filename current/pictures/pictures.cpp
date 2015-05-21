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

using namespace std;
int main(int argc, char **argv){

	TFile *myFile;
	TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);

	myFile = new TFile("now.root","READ"); /***/

	//myFile->ls();  /***/
	TH2D *h1 = (TH2D*)myFile->Get("WvsQ2_hist"); /***/

	c1->cd();
	h1->Draw("color");
	c1->Print("Q2.pdf"); /***/

	myFile->Close();

	return 0;
}
