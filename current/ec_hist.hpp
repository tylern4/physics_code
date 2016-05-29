/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef EC_HIST_H_GUARD
#define EC_HIST_H_GUARD
//
//Histogram declarations, fills, and write
//
//

int bins_EC = 500;
int bins_p_ec = 500;
float P_ec_min = 0;
float P_ec_max = 5;
float EC_min = 0;
float EC_max = 1;

TH2D *EC_sampling_fraction = new TH2D("EC_sampling_fraction","EC_sampling_fraction", 
	bins_p_ec, P_ec_min, P_ec_max, bins_EC, EC_min, EC_max);

void EC_fill(double etot, double momentum){
	double sampling_frac = etot/momentum;
	EC_sampling_fraction->Fill(momentum,sampling_frac);
}

void EC_Write(){
	EC_sampling_fraction->SetXTitle("Momentum (GeV)");
	EC_sampling_fraction->SetYTitle("Sampling Fraction");

	EC_sampling_fraction->Write();

}
#endif
