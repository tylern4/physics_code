/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef CUTS_HPP_GUARD
#define CUTS_HPP_GUARD
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
// The idea for this class would be to:
// 1) Take in the histograms I have made in pervious routines
// 2) Perform fits on the histograms
//		make a few TF1's and choose which one to fit based on type
// 3) Place the mean and sigma values for each cut into a variables
//		Which i should be able to get back later as: 
//		cut_delta_t.mean, cut_delta_t.sigma, etc.

/*
get root file with needed hists (should be the output from other part of program)
fit hists and place values in mean and sigma


*/


class Cuts
{

public:
	//Cuts();
	//~Cuts();
	double mean;
	double sigma;

	inline void CutFit(TH1D *hist, double min_value, double max_value, double *parGuess){
		//TF1 *fitFunc = new TF1("fitFunc","gaus", min_value-5, max_value+5);
		TCanvas c1;
		c1.cd();
		TF1 *fitFunc = new TF1("fitFunc","gaus", min_value, max_value);
		hist->Fit("fitFunc","","", min_value, max_value);

		mean = fitFunc->GetParameter(1);
		sigma = fitFunc->GetParameter(2); 

	} //

	
};

#endif