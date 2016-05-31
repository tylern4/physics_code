/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef CC_HIST_H_GUARD
#define CC_HIST_H_GUARD
#include "main.h"
#include "THnSparse.h"
//
//Histogram declarations, fills, and write
//
//

int bins_CC = 100;
float CC_min = 0;
float CC_max = 200;
char* L_R_C;
int sector = 6;
int segment = 18;
int PMT = 3;

TH1D *cc_hist[6][18][3];
TH1D *cc_hist_allSeg[6][3];


const Int_t ndims_cc_sparse = 4;
Int_t bins_cc_sparse[ndims_cc_sparse] = {sector, segment, PMT, bins_CC};
Double_t xmin_cc_sparse[ndims_cc_sparse] = {0.0,		0.0,			-2.0,	CC_min};
Double_t xmax_cc_sparse[ndims_cc_sparse] = {sector+1.0,	segment+1.0,	1.0,	CC_max};
Double_t x_cc_sparse[ndims_cc_sparse];

//THnSparse* cc_sparse = new THnSparseD("cc_sparse", "Histogram", ndims_cc_sparse, bins_cc_sparse, xmin_cc_sparse, xmax_cc_sparse); //
THnSparse* cc_sparse = new THnSparseD("cc_sparse", "Histogram", ndims_cc_sparse, bins_cc_sparse, xmin_cc_sparse, xmax_cc_sparse); //

void CC_fill(int cc_sector,int cc_segment,int cc_pmt,int cc_nphe){
	x_cc_sparse[0] = cc_sector;
	x_cc_sparse[1] = cc_segment;
	x_cc_sparse[2] = cc_pmt;
	x_cc_sparse[3] = cc_nphe;
	cc_sparse->Fill(x_cc_sparse);

	if(cc_pmt == -1 ) cc_pmt = 2;
	cc_hist[cc_sector-1][cc_segment-1][cc_pmt]->Fill(cc_nphe);
	cc_hist_allSeg[cc_sector-1][cc_pmt]->Fill(cc_nphe);
}

void makeHists_CC(){
	cc_sparse->GetAxis(0)->SetTitle(" cc_sector ");
	cc_sparse->GetAxis(1)->SetTitle(" cc_segment ");
	cc_sparse->GetAxis(2)->SetTitle(" cc_pmt ");
	cc_sparse->GetAxis(3)->SetTitle(" cc_nphe ");

	for (int j = 0; j < sector; j++) {
		for (int jjj = 0; jjj < PMT; jjj++) {
			if(jjj == 0) L_R_C = "both";
			if(jjj == 1) L_R_C = "right";
			if(jjj == 2) L_R_C = "left";
			sprintf(hname,"CC_sec%d_%s",j+1,L_R_C);
			sprintf(htitle,"CC sector %d %s",j+1,L_R_C);
			cc_hist_allSeg[j][jjj] = new TH1D(hname,htitle, bins_CC, CC_min, CC_max);
			for (int jj = 0; jj < segment; jj++) {
				sprintf(hname,"CC_sec%d_seg%d_%s",j+1,jj+1,L_R_C);
				sprintf(htitle,"CC sector %d segment %d %s",j+1,jj+1,L_R_C);
				cc_hist[j][jj][jjj] = new TH1D(hname,htitle, bins_CC, CC_min, CC_max);
			}
		}
	}
}


void CC_Write(){
	cc_sparse->Write();

	for (int j = 0; j < sector; j++) {
		for (int jjj = 0; jjj < PMT; jjj++) {
			cc_hist_allSeg[j][jjj]->SetYTitle("number photoelectrons");
			cc_hist_allSeg[j][jjj]->Write();
			for (int jj = 0; jj < segment; jj++){
				cc_hist[j][jj][jjj]->SetYTitle("number photoelectrons");
				cc_hist[j][jj][jjj]->Write();
			}
		}
	}

}

void CC_canvas(){
	TCanvas* c1 = new TCanvas("c1");
	for (int j = 0; j < sector; j++) {
		for (int jj = 0; jj < 1; jj++) {
			for (int jjj = 0; jjj < 1; jjj++) {
				cc_hist[j][jj][jjj]->Draw();
			}
		}
	}
}
#endif
