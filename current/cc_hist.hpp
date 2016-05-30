/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef CC_HIST_H_GUARD
#define CC_HIST_H_GUARD
#include "main.h"
//
//Histogram declarations, fills, and write
//
//

int bins_CC = 500;
float CC_min = 0;
float CC_max = 200;
char* L_R_C;
int sector = 6;
int segment = 18;
int PMT = 3;

TH1D *cc_hist[6][18][3];

void CC_fill(int cc_sector,int cc_segment,int cc_pmt,int cc_nphe){
	if(cc_pmt == -1 ) cc_pmt = 2;
	cc_hist[cc_sector-1][cc_segment-1][cc_pmt]->Fill(cc_nphe);
}

void makeHists_CC(){
	for (int j = 0; j < sector; j++) {
		for (int jj = 0; jj < segment; jj++) {
			for (int jjj = 0; jjj < PMT; jjj++) {
				if(jjj == 0) L_R_C = "center";
				if(jjj == 1) L_R_C = "right";
				if(jjj == 2) L_R_C = "left";
				sprintf(hname,"CC_sec%d_seg%d_%s",j+1,jj+1,L_R_C);
				sprintf(htitle,"CC sector %d segment %d %s",j+1,jj+1,L_R_C);
				cc_hist[j][jj][jjj] = new TH1D(hname,htitle, bins_CC, CC_min, CC_max);
			}
		}
	}
}


void CC_Write(){
	for (int j = 0; j < sector; j++) {
		for (int jj = 0; jj < segment; jj++) {
			for (int jjj = 0; jjj < PMT; jjj++) {
				cc_hist[j][jj][jjj]->SetYTitle("number photoelectrons");
				cc_hist[j][jj][jjj]->Write();
			}
		}
	}
}
#endif
