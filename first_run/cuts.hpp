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
#include "TStyle.h"
#include <math.h>
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
	double mean;
	double sigma;
	double width;

	inline void FitGaus(TH1D *hist, double min_value, double max_value){
		TF1 *fitFunc = new TF1("fitFunc","[0]*TMath::Gaus(x,[1],[2],1)", min_value, max_value);
		fitFunc->SetParameter(0, hist->GetMaximum());
		fitFunc->SetParameter(1, hist->GetMean());
		fitFunc->SetParameter(2, 1);
		fitFunc->SetParNames("height","mass","width");

		hist->Fit("fitFunc","qM0+","", min_value, max_value);

		fitFunc->SetParameter(0, hist->GetMaximum());
		fitFunc->SetParameter(1, fitFunc->GetParameter("mass"));
		fitFunc->SetParameter(2, fitFunc->GetParameter("width"));
		hist->Fit("fitFunc","qM+","", min_value, max_value);

		mean = fitFunc->GetParameter("mass");
		width = fitFunc->GetParameter("width");
		sigma = fitFunc->GetParameter("width") / (2 * sqrt(2 * log(2))); // 2.35482004503; 
		gStyle->SetOptFit(1111);

	} //

	
};

#endif