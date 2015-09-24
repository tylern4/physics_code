/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef CUTS_HPP_GUARD
#define CUTS_HPP_GUARD
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
	double mean;
	double sigma;

public:
	//Cuts();
	//~Cuts();

	
};

#endif