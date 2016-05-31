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
Int_t bins_cc_sparse[ndims_cc_sparse] = 	{sector, 		segment, 		PMT, 	bins_CC};
Double_t xmin_cc_sparse[ndims_cc_sparse] = 	{0.0,			0.0,			-2.0,	CC_min};
Double_t xmax_cc_sparse[ndims_cc_sparse] = 	{sector+1.0,	segment+1.0,	1.0,	CC_max};
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

	for (int sec_i = 0; sec_i < sector; sec_i++) {
		for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
			if(pmt_i == 0) L_R_C = "both";
			if(pmt_i == 1) L_R_C = "right";
			if(pmt_i == 2) L_R_C = "left";
			sprintf(hname,"CC_sec%d_%s",sec_i+1,L_R_C);
			sprintf(htitle,"CC sector %d %s",sec_i+1,L_R_C);
			cc_hist_allSeg[sec_i][pmt_i] = new TH1D(hname,htitle, bins_CC, CC_min, CC_max);
			for (int seg_i = 0; seg_i < segment; seg_i++) {
				sprintf(hname,"CC_sec%d_seg%d_%s",sec_i+1,seg_i+1,L_R_C);
				sprintf(htitle,"CC sector %d segment %d %s",sec_i+1,seg_i+1,L_R_C);
				cc_hist[sec_i][seg_i][pmt_i] = new TH1D(hname,htitle, bins_CC, CC_min, CC_max);
			}
		}
	}
}


void CC_Write(){
	cc_sparse->Write();
	for (int sec_i = 0; sec_i < sector; sec_i++) {
		for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
			cc_hist_allSeg[sec_i][pmt_i]->SetYTitle("number photoelectrons");
			cc_hist_allSeg[sec_i][pmt_i]->Write();
			for (int seg_i = 0; seg_i < segment; seg_i++){
				cc_hist[sec_i][seg_i][pmt_i]->SetYTitle("number photoelectrons");
				cc_hist[sec_i][seg_i][pmt_i]->Write();
			}
		}
	}

}

void CC_canvas(){
	TCanvas* can[sector][PMT];
	char can_name[50];
	for (int pmt_i = 0; pmt_i < PMT; pmt_i++) {
		for (int sec_i = 0; sec_i < sector; sec_i++) {
			if(pmt_i == 0) L_R_C = "both";
			if(pmt_i == 1) L_R_C = "right";
			if(pmt_i == 2) L_R_C = "left";

			sprintf(can_name, "Sector %d %s",sec_i+1,L_R_C);
			can[sec_i][pmt_i] = new TCanvas(can_name,can_name,1200,800);
			can[sec_i][pmt_i]->Divide(6, 3);
			for (int seg_i = 0; seg_i < segment; seg_i++) {
				if(seg_i == 0) can[sec_i][pmt_i]->cd(0);
				else can[sec_i][pmt_i]->cd((int)seg_i);
				cc_hist[sec_i][seg_i][pmt_i]->Draw("same");
			}
			can[sec_i][pmt_i]->Write();
		}
	}
}
#endif
