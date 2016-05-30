/************************************************************************/
/*  AnaTwoPion.cc                                                       */
/*                                                                      */
/*  Created by Nick Tyler, Canisius College                             */
/*                                                                      */
/************************************************************************/

#define PI 3.14159265;

#include "TTree.h"
#include "TROOT.h"
#include "TH2.h"
#include <TLorentzVector.h>
#include <TFile.h>
#include "TStyle.h"
#include "TCanvas.h"
#include <fstream>
#include "TF1.h"

//#include <omp.h>
//#include <iostream>
//#include "TObject.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include "TH1.h"


using namespace std;

  const Int_t nTgt = 6;  			//initialize number of targets types (H2,C,Fe,Ti,Pb,no Pb)
  const Int_t nFoils = 8;			//initialize number of foils
  const Int_t tgtTot = nTgt + nFoils;
  char *target[nTgt] = {"All","2H","C","Fe-Ti","Pb","No Pb"}; //label targets
  char *Foil[nFoils] = {"2H","C1","Fe","C2","Pb","C3","Ti","C4"}; //label foils

  Float_t zcut[nFoils+1] ={-22.0,-14.0,-11.0,-8.5,-6.0,-3.5,-1.0,1.5,4.0};

  Float_t zhi = 5.0;

  char hname[50];
  char htitle[500];
  char title[500];
  char ctitle[500];
  char xtitle[100];
  char ytitle[100];

  const Int_t nPart = 2;
  const Int_t nType = 2;
  char *VzAna[nType] = {"Around Targets","After Targets"};
  TH2F *hEbeam;
  TH2F *ht;
  TH2F *hVz;
  TH2F *hIM_cuts;
  TH2F *hIM_tgt[nType];
  TH2F *hPcomp[nType];
  TH2F *hIM_P[nType][3];
 
  TH1F *MassCutAround[7];
  TH1F *MassCutAfter[7];
  TH1F *MomentumCutAfter[7];
  TH1F *MomentumCutAngle[7];
  TH1F *MomentumCutAround[7];
  TH1F *MassCutAroundTgt[5][7];
  TH1F *MassCutAngleTgt[5][7];
  TH1F *MassCutAfterTgt[5][7];

  TH1F *MassCutAfterTarget[5]; 
  TH1F *MassCutAroundTarget[5];

  TH1F *AngleCutAroundTarget[5][37];
  TH1F *AngleCutAroundTargetCarbon[5][37];
  TH1F *MassCutAroundTargetAngle[5][37];
  TH1F *MassCutAfterTargetAngle[5][37];

  TH1F *hYldVsA;
  TH1F *hYldVsR;

  TH1F *Angle;
  TH1F *Mass[2];

  Double_t zFoil_mean[nFoils] = {-18.0,-12.3,-9.89,-7.42,-4.93,-2.46,0.02,2.44}; // foil peak mean
  Double_t zFoil_sigma[nFoils] = {1.17,0.258,0.249,0.252,0.235,0.234,0.245,0.231}; // foil peak sigma


  Double_t zFoil_mean2[nFoils+6] = {-21.0,-20.0,-19.0,-18.0,-17.0,-16.0,-15.0,-12.3,-9.89,-7.42,-4.93,-2.46,0.02,2.44}; // foil peak mean
  Double_t zFoil_sigma2[nFoils+6] = {0.17,0.167,0.167,0.167,0.167,0.167,0.167,0.258,0.249,0.252,0.235,0.234,0.245,0.231}; // foil peak sigma

  Int_t PolNum = 4;
  Double_t par[10];


  //declarations of functions
  Double_t polFit(Double_t *x, Double_t *par);
  void PolSubtractionAngle(char* rootFile, char* yieldFile);
  void MomPolSubtraction(char* rootFile, char* yieldFile);

  Int_t FindFoilIndex(Float_t vz);
  Int_t FindVzBehindFoilIndex(Float_t vz);
  Int_t FindTargetIndex(Float_t vz, Int_t method);
  Int_t FindNoPb(Float_t vz, Int_t method);

  Int_t GetSign(Double_t fValue);

  void WriteHist(char *RootFile);
  void BookHists_Yields();
  void WriteHists_Yields(char *RootFile);
  void BookHist();
  void Projection(char *infilename,char *outfilename);



  class ParticleTarget {
     public:
     Double_t zdiff[nFoils+6];  // z difference between target foil and decay vertex
     Double_t rdiff[nFoils+6];  // radial difference between target foil and decay vertex
     Double_t decayL[nFoils+6]; // decay length between target foil and decay vertex
     Double_t ctau[nFoils+6];   // c*tau for the decay length
     Double_t rtgt[nFoils+6];   // parent particle r for each target foil position
     Double_t rtgt_doca;      // parent particle r closest to z-axis 
     Double_t decayL_doca;    // decay length corresponding to rtgt_doca
     Double_t ctau_doca;      // c*tau corresponding to rtgt_doca
     Double_t zdiff_doca;     // z difference corresponding to rtgt_doca
     Double_t rdiff_doca;     // radial difference corresponding to rtgt_doca
     Int_t tgt_doca;          // target foil index corresponding to rtgt_doca
     Int_t nTgtTrav;          // number of traversed targets, fabs(r)<rlim
     Int_t beamlineCrossing;
     Double_t theta;};  // pos. - does not cross beamline, neg. - crosses beamline


// 
// GetParticleTargetInfo - return info with info about particle trajectory and target foils
//                  
// vz = position where particle particle intersects the z-axis
// vr = target radial position of the parent particle
// V4 = reconstructed 4-vector of the parent particle
//

ParticleTarget GetParticleTargetInfo(Float_t vz, Float_t vr, TLorentzVector V4){

    Int_t i;
    Int_t rSign;
    Int_t firstSign = 0;
    Int_t lastSign = 0;
    ParticleTarget temp; // create the ParticleTarget class

    Double_t theta = V4.Theta();
    Double_t beta = V4.Beta();

    temp.theta = (theta*180)/PI;
 

    memset(&temp.rtgt,-99,sizeof(temp.rtgt));		// initialize rtgt array elements to -99
    memset(&temp.rdiff,-99,sizeof(temp.rdiff)); 	// initialize rdiff array elements to -99
    memset(&temp.zdiff,-99,sizeof(temp.zdiff)); 	// initialize zdiff array elements to -99
    memset(&temp.decayL,-99,sizeof(temp.decayL)); 	// initialize decayL array elements to -99
    memset(&temp.ctau,-99,sizeof(temp.ctau)); 		// initialize ctau array elements to -99

    temp.tgt_doca = 0;
    temp.nTgtTrav = 0; 					// initialize coutner to zero
    temp.beamlineCrossing = 0; 				// initialize beamline crossing flag


    //for(i=0;i<=nFoils;i++){
    for(i=0;i<=nFoils+6;i++){

       //if(vz>zFoil_mean[i]){
          //temp.zdiff[i] = vz - zFoil_mean[i]; 		// z difference between decay vertex and target position
       if(vz>zFoil_mean2[i]){
          temp.zdiff[i] = vz - zFoil_mean2[i];     
          temp.rdiff[i] = temp.zdiff[i]*tan(theta); 	// radial difference between decay vertex and target
          temp.rtgt[i] = vr - temp.rdiff[i]; 		// radial distance of parent particle in target foil
          temp.decayL[i] = sqrt(temp.zdiff[i]*temp.zdiff[i] + temp.rdiff[i]*temp.rdiff[i]); // decay length
          if(beta) temp.ctau[i] = temp.decayL[i]/beta; 	// c*tau for ith target foil


         // find the smallest value of rtgt or radial position of the parent particle
         // that is closest to the z-axis
         if(fabs(temp.rtgt[i])<fabs(temp.rtgt[i-1])){
            temp.tgt_doca = i;
            temp.rtgt_doca = temp.rtgt[i];
            temp.decayL_doca = temp.decayL[i];
            temp.zdiff_doca = temp.zdiff[i];
            temp.rdiff_doca = temp.rdiff[i];
            temp.ctau_doca = temp.ctau[i];}


         // get the sign of rtgt
         rSign = GetSign(temp.rtgt[i]);
         if(firstSign==0) firstSign = rSign; // sign of first target in assembly
         if(rSign!=0) lastSign = rSign; // sign of last target in assembly

         // count number of targets traversed by parent particle before it decays
         if(fabs(temp.rtgt[i])<0.6){
         temp.nTgtTrav++;}}}


    // check if the particle vector crosses the beamline
    // firstSign - positive or negative sign of rtgt for 1st target
    // lastSign - positive or negative sign of rtgt for last target
    // If product is neg., vector crosses beamline.  If pos., vector does not cross beamline
    temp.beamlineCrossing = firstSign*lastSign;
    return temp;}



//
//  AnaTwoPion
//
//  Called in main funtion.
//  Opens a root file for analysis and sets up the branches of the tree.
//  Goes through the data and creates a root file with 2D histograms displaying mass vs momentum.
//
//  Also creates a set of projections from these histograms displaying the mass and momentum seperatly.
//  While also creating histogram cuts for the mass and momentum graphs based on the range of momentum.
//
//

void AnaTwoPion(char *fin="TwoPion.lis", char *RootFile="AnaTwoPion.root", Int_t MaxEvents=0, Int_t dEvents=10000){
  gROOT->Reset();
  Int_t i, j, k, ii, kk;
  Int_t nEvents;	
  Int_t TotEvents = 0;
  
  Float_t Plim = 0.5;
  Float_t dtlim = 1.002;
  Float_t dzlim = 3.0;
  Float_t rlim = 2.0;

  Int_t tgt[nType];
  Int_t NoPbtgt[nType];
  Int_t Ftgt[nType];

  Int_t vz_cut;
  Int_t dt_cut;
  Int_t dz_cut;
  Int_t r_cut;
  Int_t P_cut;
  Int_t isPip;
  Int_t isPim;
  Int_t cuts;

  Double_t r_pair, dz, pairM;

  Float_t Z,R;		     //floats for Get Particle Target Info
  
  TVector3 V3[nPart];        // vertex 3-vectors
  TLorentzVector V4[nPart];  // particle 4-vectors
  TLorentzVector pairV4;     // particle pair 4-vector
  TLorentzVector Beam;       // beam 4-vector
  TLorentzVector A;          // beam-pair 4-vector


  Double_t theta;
  ParticleTarget TgtInfo;
  Int_t tar = 0;
  Float_t PeakLo = 0.49; // lower limit on the peak range
  Float_t PeakHi = 0.51; // upper limit on the peak range


  // declarations for branch addresses
  Float_t eGAM, vTIME, tPHO, vX, vY, vZ, pairX, pairY, pairZ;  
  Int_t numPI, partQ[nPart], PID[nPart], SEC[nPart], isLEP[nPart], isG7[nPart];
  Int_t isEC[nPart], isECT[nPart], isCC[nPart], isMM2[nPart];
  Int_t gold[nPart], ecFID[nPart], scID[nPart], ccHITSTAT[nPart],ccSTAT[nPart];
  Int_t ccPMT[nPart], ccSEG[nPart], ccPHI[nPart], lacSTAT[nPart];
  Float_t partX[nPart], partY[nPart], partZ[nPart], pE[nPart], pX[nPart];
  Float_t pY[nPart], pZ[nPart], tPROP[nPart], ecTIME[nPart], ecTOT[nPart];
  Float_t ecIN[nPart], ecOUT[nPart], ecX[nPart], ecY[nPart], ecZ[nPart];
  Float_t scTIME[nPart], scLEN[nPart], timing[nPart], timing_e[nPart];
  Float_t scVTpi[nPart], scVTpi_pi[nPart], scVTpi_e[nPart], tofM2[nPart];
  Float_t scX[nPart], scY[nPart], scZ[nPart], ccTIME[nPart], ccTHETA[nPart];
  Float_t ccDTHETA[nPart], ccNPE[nPart], ccQF[nPart], ccCS[nPart], ccX[nPart];
  Float_t ccY[nPart], ccZ[nPart], lacTOT[nPart], lacIN[nPart], lacX[nPart];
  Float_t lacY[nPart], lacZ[nPart];

  // data files contain the trees  
  cout << "Analyzing file " << fin << endl;

  BookHist(); // declare histograms

  TFile *myFile;
  TTree *myTree;
  Int_t ncols;
  Int_t nfiles = 0;
  char rootFile[500];
  FILE *in1 = fopen(fin,"r"); 

  while (1){
	
	ncols = fscanf(in1,"%s",rootFile); 

	if (ncols<0) break;

	myFile = new TFile(rootFile); 
	// declare the tree
	myTree = (TTree*)myFile->Get("PipPimTree");
  
     myTree->SetBranchAddress("egam",&eGAM);
     myTree->SetBranchAddress("vTime",&vTIME);
     myTree->SetBranchAddress("tpho",&tPHO);
     myTree->SetBranchAddress("vx",&vX);
     myTree->SetBranchAddress("vy",&vY);
     myTree->SetBranchAddress("vz",&vZ);
     myTree->SetBranchAddress("vx_pair",&pairX);
     myTree->SetBranchAddress("vy_pair",&pairY);
     myTree->SetBranchAddress("vz_pair",&pairZ);

     myTree->SetBranchAddress("Npi",&numPI);
     myTree->SetBranchAddress("x",&partX);
     myTree->SetBranchAddress("y",&partY);
     myTree->SetBranchAddress("z",&partZ);
     myTree->SetBranchAddress("Q",&partQ);
     myTree->SetBranchAddress("Pid",&PID);
     myTree->SetBranchAddress("Sec",&SEC);
     myTree->SetBranchAddress("E",&pE);
     myTree->SetBranchAddress("Px",&pX);
     myTree->SetBranchAddress("Py",&pY);
     myTree->SetBranchAddress("Pz",&pZ);
     myTree->SetBranchAddress("Tprop",&tPROP);
     myTree->SetBranchAddress("IsLep",&isLEP);
     myTree->SetBranchAddress("IsLepG7",&isG7);
     myTree->SetBranchAddress("IsLepG7ec",&isEC);
     myTree->SetBranchAddress("IsLepG7ect",&isECT);
     myTree->SetBranchAddress("IsLepG7cc",&isCC);
     myTree->SetBranchAddress("IsLepG7mm2",&isMM2);
     myTree->SetBranchAddress("Golden",&gold);
     myTree->SetBranchAddress("EC_time",&ecTIME);
     myTree->SetBranchAddress("EC",&ecTOT);
     myTree->SetBranchAddress("ECin",&ecIN);
     myTree->SetBranchAddress("ECout",&ecOUT);
     myTree->SetBranchAddress("ECx",&ecX);
     myTree->SetBranchAddress("ECy",&ecY);
     myTree->SetBranchAddress("ECz",&ecZ);
     myTree->SetBranchAddress("ECfid",&ecFID);
     myTree->SetBranchAddress("SC_time",&scTIME);
     myTree->SetBranchAddress("scLen",&scLEN);
     myTree->SetBranchAddress("scId",&scID);
     myTree->SetBranchAddress("Timing",&timing);
     myTree->SetBranchAddress("Timing_e",&timing_e);
     myTree->SetBranchAddress("scvT_pi",&scVTpi);
     myTree->SetBranchAddress("scvT_pi_pi",&scVTpi_pi);
     myTree->SetBranchAddress("scvT_pi_e",&scVTpi_e);
     myTree->SetBranchAddress("TOF_MassSq",&tofM2);
     myTree->SetBranchAddress("SCx",&scX);
     myTree->SetBranchAddress("SCy",&scY);
     myTree->SetBranchAddress("SCz",&scZ);
     myTree->SetBranchAddress("CCstat",&ccSTAT);
     myTree->SetBranchAddress("CC_time",&ccTIME);
     myTree->SetBranchAddress("CChit_stat",&ccHITSTAT);
     myTree->SetBranchAddress("CCtheta",&ccTHETA);
     myTree->SetBranchAddress("CCdtheta",&ccDTHETA);
     myTree->SetBranchAddress("CCnpe",&ccNPE);
     myTree->SetBranchAddress("CC_QF",&ccQF);
     myTree->SetBranchAddress("CC_CS",&ccCS);
     myTree->SetBranchAddress("CCpmt",&ccPMT);
     myTree->SetBranchAddress("CCseg",&ccSEG);
     myTree->SetBranchAddress("CCphimatch",&ccPHI);
     myTree->SetBranchAddress("CCx",&ccX);
     myTree->SetBranchAddress("CCy",&ccY);
     myTree->SetBranchAddress("CCz",&ccZ);
     myTree->SetBranchAddress("LAC",&lacTOT);
     myTree->SetBranchAddress("LACin",&lacIN);
     myTree->SetBranchAddress("LACstat",&lacSTAT);
     myTree->SetBranchAddress("LACx",&lacX);
     myTree->SetBranchAddress("LACy",&lacY);
     myTree->SetBranchAddress("LACz",&lacZ);

	// loop over events
	nEvents = (Int_t)myTree->GetEntries(); 
	i = 0; // inititalize the event counter for each file

    	while(i<nEvents && (!MaxEvents || TotEvents<MaxEvents)){

      	// Initialize cuts to false
      	   vz_cut = 0;
      	   dt_cut = 0;
      	   dz_cut = 0;
      	   r_cut = 0;
      	   isPip = 0;
      	   isPim = 0;
      	   P_cut = 0;
      	   cuts = 0;

      	   memset(&tgt,0,sizeof(tgt));
      	   memset(&NoPbtgt,0,sizeof(NoPbtgt));
      	   memset(&Ftgt,0,sizeof(Ftgt));

     	   myTree->GetEntry(i); // retrieve the event from the tree

      	   // print statement for certain eevent intervals
      	   if(!(i%dEvents)) cerr<<i<<" / "<<nEvents<<" / "<<nfiles<<"\r"; 

	        Beam.SetPxPyPzE(0,0,eGAM,eGAM); // beam 4-vector

             // pair kinematics
             for(ii=0;ii<nPart;ii++){
                V3[ii].SetXYZ(partX[ii],partY[ii],partZ[ii]);
                V4[ii].SetPxPyPzE(pX[ii],pY[ii],pZ[ii],pE[ii]);
             }

      	   pairV4 = V4[0] + V4[1];                              // particle pair 4-vector
	   pairM = pairV4.M();                                  // pair invariant mass



	   A = Beam - pairV4;                                   // beam minus pair 4-vector
	
      	   dz = V3[0].Z() - V3[1].Z();                          // z vertex difference
      	   r_pair = sqrt(pairX*pairX + pairY*pairY);            // radial vertex (x,y) for the particle pair

      	   P_cut = (V4[0].P()>=Plim && V4[1].P()>=Plim);        // momentum cut
      	   vz_cut = (pairZ>=zcut[0] && pairZ>=zcut[nFoils]);    // vertex z cut
      
      	   dz_cut = (fabs(dz)<=dzlim);                          // vertex z difference cut
      	   dt_cut = (fabs(scTIME[1]-scTIME[0])<=dtlim);         // vertex time difference cut
      	   r_cut = (r_pair<=rlim);                              // target radius cut 
     	   isPim = (PID[0]==9);                                 // pi- GEANT id
      	   isPip = (PID[1]==8);                                 // pi+ GEANT id

      	   cuts=(dz_cut && r_cut && dt_cut && isPim && isPip);  // total cuts

	  	   hEbeam->Fill(Beam.E(),0);                            // beam histogram
	  	   ht->Fill(A.Mag2(),0);                                // histogram of t (mom. transfer)

      	   if(isPim && isPip){                                       // cut histograms
        	      hIM_cuts->Fill(pairM,1);                               // no cuts
        	      if(dz_cut) hIM_cuts->Fill(pairM,2);                    // dz cut only
        	      if(r_cut) hIM_cuts->Fill(pairM,3);                     // r cut only
        	      if(dt_cut) hIM_cuts->Fill(pairM,4);                    // dt cut only
        	      if(dz_cut && r_cut) hIM_cuts->Fill(pairM,5);           // dz and r cuts
        	      if(dz_cut && r_cut && dt_cut) hIM_cuts->Fill(pairM,6); // dz, r, and dt cuts
      	   }

      for(j=0;j<nType;j++){                                     // loop over vz type, [0,in target], [1,after target]
         tgt[j] = FindTargetIndex(pairZ,j);                     // target index (1-4)
         NoPbtgt[j] = FindNoPb(pairZ,j);                        // no Pb index (5)
         if(NoPbtgt[j]) NoPbtgt[j] = nTgt-1; 
         Ftgt[j] = FindFoilIndex(pairZ);                        // foil index (6-13)
         if(Ftgt[j]) Ftgt[j] += nTgt-1;
         	         
	 if (cuts && tgt[j]){
	    if(j==0){
	       for(k=0; k<=36; k++){
	          if(TgtInfo.theta >= k*5){
	      	  MassCutAroundTargetAngle[0][k]->Fill(pairM);
  	     	  MassCutAroundTargetAngle[tgt[j]][k]->Fill(pairM);}
		}

		MassCutAroundTarget[(tgt[j])]->Fill(pairM);

		if(pairV4.P()>0 and pairV4.P()<=0.35){ 
  		   MomentumCutAround[1]->Fill(pairV4.P());
  		   MassCutAround[1]->Fill(pairM);
		   MassCutAroundTgt[(tgt[j])][1]->Fill(pairM);}

  		else if(pairV4.P()>0.35 and pairV4.P()<=0.70){
		   MomentumCutAround[2]->Fill(pairV4.P());
		   MassCutAround[2]->Fill(pairM);
		   MassCutAroundTgt[(tgt[j])][2]->Fill(pairM);}

		else if(pairV4.P()>0.70 and pairV4.P()<=1.05){
		   MomentumCutAround[3]->Fill(pairV4.P());
		   MassCutAround[3]->Fill(pairM);
		   MassCutAroundTgt[(tgt[j])][3]->Fill(pairM);}

		else if(pairV4.P()>1.05 and pairV4.P()<=1.40){
		   MomentumCutAround[4]->Fill(pairV4.P());
		   MassCutAround[4]->Fill(pairM);
		   MassCutAroundTgt[(tgt[j])][4]->Fill(pairM);}

		else if(pairV4.P()>1.40 and pairV4.P()<=1.75){
		   MomentumCutAround[5]->Fill(pairV4.P());
		   MassCutAround[5]->Fill(pairM);
		   MassCutAroundTgt[(tgt[j])][5]->Fill(pairM);}

		else if(pairV4.P()>1.75 and pairV4.P()<=2.10){
		   MomentumCutAround[6]->Fill(pairV4.P());
		   MassCutAround[6]->Fill(pairM);
		   MassCutAroundTgt[(tgt[j])][6]->Fill(pairM);}		      
		   }

	    else if(j==1){
	       for(k=0; k<=36; k++){
	          if(TgtInfo.theta >= k*5){
  	             MassCutAfterTargetAngle[0][k]->Fill(pairM);
                     MassCutAfterTargetAngle[tgt[j]][k]->Fill(pairM);}
		}

		MassCutAfterTarget[(tgt[j])]->Fill(pairM);

		if(pairV4.P()>0 and pairV4.P()<=0.35){ 
  		   MomentumCutAfter[1]->Fill(pairV4.P());
  		   MassCutAfter[1]->Fill(pairM);
		   MassCutAfterTgt[(tgt[j])][1]->Fill(pairM);}

  		else if(pairV4.P()>0.35 and pairV4.P()<=0.70){
		   MomentumCutAfter[2]->Fill(pairV4.P());
		   MassCutAfter[2]->Fill(pairM);
		   MassCutAfterTgt[(tgt[j])][2]->Fill(pairM);}

		else if(pairV4.P()>0.70 and pairV4.P()<=1.05){
		   MomentumCutAfter[3]->Fill(pairV4.P());
		   MassCutAfter[3]->Fill(pairM);
		   MassCutAfterTgt[(tgt[j])][3]->Fill(pairM);}

		else if(pairV4.P()>1.05 and pairV4.P()<=1.40){
		   MomentumCutAfter[4]->Fill(pairV4.P());
		   MassCutAfter[4]->Fill(pairM);
		   MassCutAfterTgt[(tgt[j])][4]->Fill(pairM);}

		else if(pairV4.P()>1.40 and pairV4.P()<=1.75){
		   MomentumCutAfter[5]->Fill(pairV4.P());
		   MassCutAfter[5]->Fill(pairM);
		   MassCutAfterTgt[(tgt[j])][5]->Fill(pairM);}

		else if(pairV4.P()>1.75 and pairV4.P()<=2.10){
		   MomentumCutAfter[6]->Fill(pairV4.P());
		   MassCutAfter[6]->Fill(pairM);
		   MassCutAfterTgt[(tgt[j])][6]->Fill(pairM);}		      
		}
	}



	   TgtInfo = GetParticleTargetInfo(pairZ, sqrt(pairX*pairX+pairY*pairY), pairV4);

	         switch(TgtInfo.tgt_doca){
  		    case 1:
   		    case 2:
  		    case 3:
  		    case 4:
  		    case 5:
  		    case 6:
  		    case 7:
        	       tar = 1; break;
  		    case 8: 
  		    case 10: 
  		    case 12: 
  		    case 14: 
    	   	       tar = 2; break;
  		    case 9: 
  		    case 13: 
    	   	       tar = 3; break;
  		    case 11: 
    		       tar = 4; break;
  		    default: 
    	 	       tar = 0; break;}

           if (cuts && tar){

	      if (pairV4.M()>=PeakLo && pairV4.M()<=PeakHi)Angle->Fill(TgtInfo.theta);// Angle between the kaon and z axis for mass range
              Mass[1]->Fill(pairV4.M());

		if(pairV4.P()>0 and pairV4.P()<=0.35){ 
  		   MomentumCutAngle[1]->Fill(pairV4.P());
		   MassCutAngleTgt[tar][1]->Fill(pairV4.M());}

  		else if(pairV4.P()>0.35 and pairV4.P()<=0.70){
		   MomentumCutAngle[2]->Fill(pairV4.P());
		   MassCutAngleTgt[tar][2]->Fill(pairV4.M());}

		else if(pairV4.P()>0.70 and pairV4.P()<=1.05){
		   MomentumCutAngle[3]->Fill(pairV4.P());
		   MassCutAngleTgt[tar][3]->Fill(pairV4.M());}

		else if(pairV4.P()>1.05 and pairV4.P()<=1.40){
		   MomentumCutAngle[4]->Fill(pairV4.P());
		   MassCutAngleTgt[tar][4]->Fill(pairV4.M());}

		else if(pairV4.P()>1.40 and pairV4.P()<=1.75){
		   MomentumCutAngle[5]->Fill(pairV4.P());
		   MassCutAngleTgt[tar][5]->Fill(pairV4.M());}

		else if(pairV4.P()>1.75 and pairV4.P()<=2.10){
		   MomentumCutAngle[6]->Fill(pairV4.P());
		   MassCutAngleTgt[tar][6]->Fill(pairV4.M());}

	      for(k=0; k<=36; k++){

	         if(TgtInfo.theta >= k*5){
		    AngleCutAroundTarget[0][k]->Fill(pairV4.M());
	            AngleCutAroundTarget[tar][k]->Fill(pairV4.M());}}
	  }


//carbon target stuff
	  if (cuts && j==0){
	      for(k=0; k<=36; k++){

	         switch(TgtInfo.tgt_doca){
  		    case 8: 
        	       tar = 1; break;
  		    case 10: 
        	       tar = 2; break;
  		    case 12: 
        	       tar = 3; break;
  		    case 14: 
        	       tar = 4; break;
  		    default: 
    	 	       tar = 0; break;}

	         if(TgtInfo.theta >= k*5){
		    AngleCutAroundTargetCarbon[0][k]->Fill(pairM);
	            AngleCutAroundTargetCarbon[tar][k]->Fill(pairM);}}
}

  
         if(cuts && tgt[j]){ // histograms with std cuts parsed by target
		   hVz->Fill(pairZ,j);
		   hPcomp[j]->Fill(V4[0].P(),V4[1].P());
		   hIM_tgt[j]->Fill(pairM,0);
		   hIM_tgt[j]->Fill(pairM,tgt[j]);
		   if(NoPbtgt[j]) hIM_tgt[j]->Fill(pairM,NoPbtgt[j]);
		   hIM_tgt[j]->Fill(pairM,Ftgt[j]);
		   hIM_P[j][0]->Fill(pairM,V4[0].P());
		   hIM_P[j][1]->Fill(pairM,V4[1].P());
		   hIM_P[j][2]->Fill(pairM,pairV4.P());
	  	}
	}
	  i++; 		// increment event counter
	  TotEvents++; // increment total event counter 
    }
    	cerr<<endl;
	myTree->Delete(); 						// delete Tree object
	myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
    	nfiles++; 							// increment file counter
  }

  fclose(in1); 							// close file with input file list
  cout<<TotEvents<<" events in "<<nfiles<< " files."<<endl; // print out stats

  WriteHist(RootFile); 						// write histograms to a file
  
}



void BgSub_Sideband(char *rootFile, char* yieldFile)
{

	Float_t Lmar = 0.125;
	Float_t Rmar = 0.125;
	Float_t yoff = 1.5;
	
	Float_t xLo = 0.41;  // lower value of x-axis for drawing histogram
	Float_t xHi = 0.6;   // upper value of x-axis for drawing histogram
	Float_t PeakLo = 0.49; // lower limit on the peak range
	Float_t PeakHi = 0.51; // upper limit on the peak range
	Float_t Ks_width = PeakHi - PeakLo; // width of the peak range
	Float_t SidebandLo = PeakLo - 0.5*Ks_width; // lower limit of left sideband
	Float_t SidebandHi = PeakHi + 0.5*Ks_width; // upper limit of right sideband
	
	Int_t i, j, k, t, m, a, n;
	Float_t x;
	Float_t yield, yield_peak, yield_sb, yield_2H;
	
	TH1F *hist; // original histogram
	TH1F *hPeak; // histogram of peak range
	TH1F *hSideband; // histogram of sidebands
	
	char *AfterOrAround[3];
	AfterOrAround[1] = "After";
	AfterOrAround[2] = "Around";
	char *Target[4];
	Target [1] = "2H";
	Target [2] = "C";
	Target [3] = "FeTi";
	Target [4] = "Pb";
	char inName[200];
	char outName[200];
	char name[200];
	char YieldHistName[200];
	char ErrorHistName[200];
	Double_t error;
	Double_t function;
	
	TH1F *Sideband[3][7];
	TH1F *SidebandYield[3][7];
	TH1F *SidebandError[3][7];
	TH1F *SidebandRatio[3][7];
		
	// create canvas
	char xtitle[100];
	char ytitle[100];
	char title[100];
	sprintf(title,"Sideband Analysis"); // canvas title
	TCanvas *can1 = new TCanvas("can1",title,0,0,600,600); // create the canvas
	
	gStyle->SetOptStat(1111);
	can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
	can1->SetBorderSize(5);  
	can1->SetFillStyle(4000); 
	
	sprintf(xtitle,"Invariant Mass (GeV)"); // set the x-axis title
	sprintf(ytitle,"Counts"); // set the y-axis title
	
	// data files contain the trees
	printf("Analyzing file %s\n",rootFile);  
	TFile *fd = new TFile(rootFile,"UPDATE"); //open up the ROOT file
	

     for (a = 1; a <=2 ; a++){   //after = 1, around = 2
       for (m = 0; m <= 6; m++){  // mometum ranges

       sprintf(outName,"SidebandMassCut%s_%i",AfterOrAround[a],m);
       sprintf(name,"Mass Cut %s Sideband range %i",AfterOrAround[a],m);
       Sideband[a][m] = new TH1F(outName,name,220, 0.5, 220.5);
       
       sprintf(outName,"RatioSidebandMassCut%s_%i",AfterOrAround[a],m);
       sprintf(name,"Mass Cut %s Sideband range %i ratio",AfterOrAround[a],m);
       SidebandRatio[a][m] = new TH1F(outName,name,220, 0.5, 220.5);

       sprintf(YieldHistName,"YeildMassCut%sSideband_%i",AfterOrAround[a],m);       
       SidebandYield[a][m] = new TH1F(YieldHistName,YieldHistName,220,0.5 ,220.5 );
     
       sprintf(ErrorHistName,"ErrorMassCut%sSideband_%i",AfterOrAround[a],m);
       SidebandError[a][m] = new TH1F(ErrorHistName,ErrorHistName,220,0.5 ,220.5 );



       for(t=1;t<=4;t++){    //target ranges
     
         if (t==1) n = 2;
         else if (t == 2) n = 12;
         else if (t == 3) n = 56;
         else if (t == 4) n = 208;
          	    
         sprintf(inName,"MassCut%sTgt_%s_%i", AfterOrAround[a],Target[t],m);

	    hist = (TH1F*)fd->Get(inName); // get the histogram from the ROOT file

	    hSideband = (TH1F*)hist->Clone(inName);  // clone the original histogram for the sideband histogram  

	    hPeak = (TH1F*)hist->Clone(inName);  // clone the original histogram for the peak histogram  
	
	    for(k=1;k<=hist->GetNbinsX();k++){   // loop to fill the peak and sideband histograms
	  	    x = hist->GetBinCenter(k); // x-axis value for bin k

		    // fill peak histogram with counts from original histogram for bins in range
		    if(x>=PeakLo && x<PeakHi){
			    hPeak->SetBinContent(k,hist->GetBinContent(k));
		    }
		    else{
			    hPeak->SetBinContent(k,0.0);
		    }
		
		    // fill sideband histogram when bin is in the proper ranges
		    if(x>=SidebandLo && x<PeakLo){
			    hSideband->SetBinContent(k,hist->GetBinContent(k));
		    }else if(x>=PeakHi && x<SidebandHi){
			    hSideband->SetBinContent(k,hist->GetBinContent(k));
		    }else{
		 	    hSideband->SetBinContent(k,0.0);
		    }
	    }


	   
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
	    hist->GetXaxis()->SetRangeUser(xLo,xHi);
	    hist->SetLineWidth(2);
	    hist->Draw();
	    hist->Write("hist");
	
	    // draw the peak histogram
	    hPeak->SetLineWidth(2);
	    hPeak->SetFillColor(2);
	    hPeak->Draw("same");
	    hPeak->Write("hPeak");
 
	
	    // draw the sideband histogram
	    hSideband->SetLineWidth(2);
	    hSideband->SetFillColor(4);
	    hSideband->Draw("same");
	    hSideband->Write("hSideband");




	    //create the image files
	    char OutCan[200];
	    sprintf(OutCan,"pictures/Sideband%s.gif",inName);
	    can1->Print(OutCan);
	   
	    sprintf(OutCan,"results/%s.eps",name);
	    can1->Print(OutCan);

	    // sum counts in the peak range
	    yield_peak = hPeak->Integral(hPeak->FindBin(PeakLo),hPeak->FindBin(PeakHi));

            // sum counts in sidebands
	    yield_sb = hSideband->Integral(hSideband->FindBin(SidebandLo),hSideband->FindBin(SidebandHi));

	    yield = yield_peak - yield_sb; // subtract background from peak sum

	// open text file for the yields
	char OutFile[100];
	sprintf(OutFile,"results/%s.yld",inName);
	ofstream fout(OutFile); 

	fout<<"\t"<<yield<<"\t"<<sqrt(yield)<<endl;


         if (t == 1) yield_2H = yield;
	
         fd->cd();

         SidebandYield[a][m]->SetBinContent(n,yield);
        
         function = ((2*yield)/(n*yield_2H));

         if (function > 0) Sideband[a][m]->SetBinContent(n,function);
         
         SidebandRatio[a][m]->SetBinContent(n,(yield/yield_2H));
             
         error = (function*sqrt((1/yield_2H)+(1/yield)));
         SidebandError[a][m]->SetBinError(n,error);
         }
      
     SidebandRatio[a][m]->Write();   
     SidebandYield[a][m]->Write();
     SidebandError[a][m]->Write();
	Sideband[a][m]->Write();}}

	fd->Close();

	
}


int main(int argc, char **argv){


     float time1 = clock();


     if(argc < 3){
        cout<<"Error : Too Few Arguments" <<endl;
        cout<<"To run: ./AnaTwoPion [Infilename] [Outfilename]"<<endl;
        exit(0);}


     char  infilename[128];
     char  outfilename[128];


     sprintf(infilename,"%s",argv[1]);
     sprintf(outfilename,"%s",argv[2]);



	//AnaTwoPion(infilename,outfilename);
	//Projection(outfilename,outfilename);
	BgSub_Sideband(outfilename,"BgSub_Sideband");
   PolSubtractionAngle(outfilename, "PolSubtractionAngle");
	MomPolSubtraction (outfilename, "MomPolSubtraction");
    
	     

     //time to complete function
  	float minutes = 0;
	float seconds = 0;
  	float time2 = clock();
  	minutes = (time2 - time1)/1000000;
  	minutes = (minutes)/60;
  	seconds = fmod(minutes,1);
	minutes = minutes-seconds;
	seconds = seconds*60;


  	if (minutes==0) cout<<endl<<"Completed in "<<seconds<<" seconds."<<endl<<endl;
  	else cout<<endl<<"Completed in "<<minutes<<" minutes and "<<seconds<<" seconds."<<endl<<endl;


	return 0;
}




Int_t GetSign(Double_t fValue){
   Double_t sign = 0;
   if(fValue) sign = fValue/fabs(fValue);
   return (Int_t)sign;}




//
// Makes projections from 2D histograms
//
//
//
//
//
void Projection(char *infilename="AnaTwoPion.root",char *outfilename="AnaTwoPion.root"){

	
  	TH1D *h1D;
  	TH2F *h2D;  	
  	  	
	Int_t binLo = 0;
	Int_t binHi = -1;
	  
  	TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
  	c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  	c1->SetBorderSize(5); 
  	gStyle->SetOptStat(0);
  	gStyle->SetPalette(1);
  	c1->SetFillStyle(4000);

  	char strname[50];


	TFile *myFile = new TFile(outfilename,"UPDATE");
	TFile *infile = new TFile(infilename,"READ"); 
	
	
	
  	c1->cd();
  	gPad->SetLeftMargin(1.5);
  	gPad->SetRightMargin(1.5);
  	gPad->SetFillColor(0);


	cout<<endl<<"Analyzing file "<<infilename<<endl;
			
		h2D = (TH2F*)infile->Get("hIM_P0_2");
		sprintf(strname,"ProjectionPionPairMassAround");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
		h1D->SetTitle("#pi Pair Mass Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();
			
		h2D = (TH2F*)infile->Get("hIM_P1_2");
		sprintf(strname,"ProjectionPionPairMassAfter");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
		h1D->SetTitle("#pi Pair Mass After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();

		h2D = (TH2F*)infile->Get("hIM_P0_2");
		sprintf(strname,"ProjectionPionPairMomentumAround");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
		h1D->SetTitle("#pi Pair Momentum Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();
  					
     	h2D = (TH2F*)infile->Get("hIM_P1_2");
		sprintf(strname,"ProjectionPionPairMomentumAfter");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
		h1D->SetTitle("#pi Pair Mometum After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();


	
	myFile->Write(outfilename);
	myFile->Close(outfilename);

}

// 
// FindVzBehindFoilIndex - return index number of target foil where vz is directly behind
//                  
//				vz = z vertex position
//
Int_t FindVzBehindFoilIndex(Float_t vz){
	Int_t i;
	Int_t ret = 0;

  
	Float_t z1, z2;

	for(i=0;i<nFoils;i++){  
           z1 = zFoil_mean[i] - 3.0*zFoil_sigma[i];
	   z2 = (i==7) ? 4.9 : zFoil_mean[i+1] - 3.0*zFoil_sigma[i+1];
	   if(vz>=z1 && vz<z2) ret = i+1;
        }
  return ret;
}

// 
// FindTargetIndex - return index number of target
//                  
//				vz = z vertex position
//				method = [0 = target],[1 = after target]
//
Int_t FindTargetIndex(Float_t vz, Int_t method){
	Int_t ret = 0;
	Int_t index = (method==0) ? FindFoilIndex(vz) : FindVzBehindFoilIndex(vz);

	switch(index){
  		case 1: 
        	   ret = 1; break;
  		case 2: 
  		case 4: 
  		case 6: 
  		case 8: 
    	   	   ret = 2; break;
  		case 3: 
  		case 7: 
    	   	   ret = 3; break;
  		case 5: 
    		   ret = 4; break;
  		default: 
    	 	   ret = 0; break;    
  	}
  	return ret;
}

// 
// FindFoilIndex - return index number of target foils
//                  
//				vz = z vertex position
//
Int_t FindFoilIndex(Float_t vz){
	Int_t i;
  	Int_t ret = 0;

	Float_t z1, z2;
	
  
  	for(i=0;i<nFoils;i++){
        z1 = zFoil_mean[i] - 3.0*zFoil_sigma[i];
	z2 = zFoil_mean[i] + 3.0*zFoil_sigma[i];
        if(vz>=z1 && vz<=z2) ret = i+1;

     }
  return ret;
}

// 
// FindNoPb - determine if target is not Pb
//                  
//				vz = z vertex position
//				method = [0 = target],[1 = after target]
//
Int_t FindNoPb(Float_t vz, Int_t method){
	Int_t ret = 0;
  	Int_t ii = (method==0) ? FindFoilIndex(vz) : FindVzBehindFoilIndex(vz);
  	
	if(ii>0 && ii!=5) ret = 1;
  	return ret;
}

// 
// WriteHist - routine to write histograms to the output file
//             Called in the AnaTwoPion() routine
//
void WriteHist(char *RootFile){

  Int_t i=0, j=0 ,k=0 , l=0, jj=0;
  
  TFile *fa = new TFile(RootFile,"RECREATE");
  fa->cd();
  
  hEbeam->SetXTitle("Photon Energy (GeV)");  
  hEbeam->Write();
  ht->SetXTitle("Momentum Transfer");
  ht->Write();
  hIM_cuts->SetXTitle("Mass (GeV)");
  hIM_cuts->Write();
  hVz->SetXTitle("Z vertex");
  hVz->Write();


  Angle->Write(); 
  Mass[1]->Write();

 
  for(jj=0; jj<=4; jj++){
     for (j=0;j<=36;j++){
        AngleCutAroundTarget[jj][j]->SetXTitle("Mass (GeV)");
        AngleCutAroundTarget[jj][j]->GetXaxis()->CenterTitle();
        AngleCutAroundTarget[jj][j]->Write();

        AngleCutAroundTargetCarbon[jj][j]->SetXTitle("Mass (GeV)");
        AngleCutAroundTargetCarbon[jj][j]->GetXaxis()->CenterTitle();
        AngleCutAroundTargetCarbon[jj][j]->Write();

        MassCutAroundTargetAngle[jj][j]->SetXTitle("Mass (GeV)");
        MassCutAroundTargetAngle[jj][j]->GetXaxis()->CenterTitle();
        MassCutAroundTargetAngle[jj][j]->Write();

        MassCutAfterTargetAngle[jj][j]->SetXTitle("Mass (GeV)");
        MassCutAfterTargetAngle[jj][j]->GetXaxis()->CenterTitle();
        MassCutAfterTargetAngle[jj][j]->Write();}
     }

  for(k=1;k<=4;k++){
     MassCutAfterTarget[k]->SetXTitle("Mass (GeV)");
     MassCutAfterTarget[k]->GetXaxis()->CenterTitle();
     MassCutAfterTarget[k]->Write();
   
     MassCutAroundTarget[k]->SetXTitle("Mass (GeV)");
     MassCutAroundTarget[k]->GetXaxis()->CenterTitle();
     MassCutAroundTarget[k]->Write();}


  for(i=1;i<=6;i++){
    MassCutAround[i]->SetXTitle("Mass (GeV)");
    MassCutAround[i]->GetXaxis()->CenterTitle();
    MassCutAround[i]->Write();

    MassCutAfter[i]->SetXTitle("Mass (GeV)");
    MassCutAfter[i]->GetXaxis()->CenterTitle();
    MassCutAfter[i]->Write();

    MomentumCutAfter[i]->SetXTitle("Momentum (GeV)");
    MomentumCutAfter[i]->GetXaxis()->CenterTitle();
    MomentumCutAfter[i]->Write();  

    MomentumCutAngle[i]->SetXTitle("Momentum (GeV)");
    MomentumCutAngle[i]->GetXaxis()->CenterTitle();
    MomentumCutAngle[i]->Write();
  
    MomentumCutAround[i]->SetXTitle("Momentum (GeV)");
    MomentumCutAround[i]->GetXaxis()->CenterTitle();
    MomentumCutAround[i]->Write();
  
    for(k=1;k<=4;k++){  
       MassCutAroundTgt[k][i]->SetXTitle("Mass (GeV)");
       MassCutAroundTgt[k][i]->GetXaxis()->CenterTitle();
       MassCutAroundTgt[k][i]->Write();

       MassCutAngleTgt[k][i]->SetXTitle("Mass (GeV)");
       MassCutAngleTgt[k][i]->GetXaxis()->CenterTitle();
       MassCutAngleTgt[k][i]->Write();
         
       MassCutAfterTgt[k][i]->SetXTitle("Mass (GeV)");
       MassCutAfterTgt[k][i]->GetXaxis()->CenterTitle();
       MassCutAfterTgt[k][i]->Write();}
   }


  for(i=0;i<nPart;i++){
    hIM_tgt[i]->SetXTitle("Mass (GeV)");
    hIM_tgt[i]->SetYTitle("Target");
    hIM_tgt[i]->Write();
    hPcomp[i]->SetXTitle("Momentum (GeV)");
    hPcomp[i]->SetYTitle("Momentum (GeV)");
    hPcomp[i]->Write();
    for(j=0;j<3;j++){
      hIM_P[i][j]->SetXTitle("Mass (GeV)");
      hIM_P[i][j]->SetYTitle("Momentum (GeV)");
      hIM_P[i][j]->Write();
    }
  }
  fa->Close();
}

// 
// BookHists_Yields - create histograms for plotYields()
//                  
//
void BookHists_Yields(){
  Int_t i;
  
  Float_t Alo = -20.0;
  Float_t Ahi = 220.0;
  Int_t Abins = 10*Int_t(Ahi-Alo);
  Float_t Abwid = (Ahi-Alo)/float(Abins);

  Float_t Rlo = 0.0;
  Float_t Rhi = 10.0;
  Int_t Rbins = 10*Int_t(Rhi-Rlo);
  Float_t Rbwid = (Rhi-Rlo)/float(Abins);
  
  Alo = 0.5*Abwid;
  Ahi = 0.5*Abwid;

  Rlo = 0.5*Rbwid;
  Rhi = 0.5*Rbwid;  

  sprintf(hname,"hYldVsA");  
  sprintf(htitle,"K_{s}, Yield vs A");
  hYldVsA = new TH1F(hname,htitle,Abins,Alo,Ahi);

  sprintf(hname,"hYldVsR");  
  sprintf(htitle,"K_{s} Results vs Radius");
  hYldVsR = new TH1F(hname,htitle,Rbins,Rlo,Rhi);
}

// 
// WriteHists_Yields - write histograms to ROOT file for plotYields()
//                  
//
void WriteHists_Yields(char *RootFile){
  	Int_t i;
  	Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
     	Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};

  	TFile *fa = new TFile(RootFile,"RECREATE");
  	fa->cd();

  	hYldVsA->SetMarkerSize(1.25);
  	hYldVsA->SetMarkerStyle(mkr[0]);
  	hYldVsA->SetMarkerColor(lcol[0]);
  	hYldVsA->Write();

  	hYldVsR->SetMarkerSize(1.25);
  	hYldVsR->SetMarkerStyle(mkr[0]);
  	hYldVsR->SetMarkerColor(lcol[0]);
  	hYldVsR->Write();
}


// 
// BookHist - routine to set up histograms for the AnaTwoPion() routine
//
void BookHist(){

  Int_t i, j, jj, k;
  char *between, *target;
  double range, range2;
  
  Int_t nIM = 500;
  Float_t IMl = 0.0;
  Float_t IMh = 1.2;
  
  Int_t nP = 500;
  Float_t Pl = 0.0;
  Float_t Ph = 3.0;
  
  char *Plabel[3] = {"#pi^{+}","#pi^{-}","#pi^{+}#pi^{-}"};
  
  Int_t nVz = 270;
  Float_t zlo = -22.0;
  
  
    for (i=1;i<=4;i++){
        
      if(i==1) target = "2H";
      else if(i==2) target = "C";
      else if(i==3) target = "FeTi";
      else if(i==4) target = "Pb";
      
      sprintf(hname,"MassCutAfterTgt_%s_0",target);
      sprintf(title,"Mass Cut After Target %s",target,between);
      MassCutAfterTarget[i] = new TH1F(hname,title,10000, 0, 2);
         	            
      sprintf(hname,"MassCutAroundTgt_%s_0",target);
      sprintf(title,"Mass Cut Around Target %s",target,between);
      MassCutAroundTarget[i] = new TH1F(hname,title,10000, 0, 2);}
  
  sprintf(hname,"hEbeam");
  sprintf(title,"Photon Energy (GeV)");
  hEbeam = new TH2F(hname,title,100,0.5,4.0,2,-0.5,1.5);

  sprintf(hname,"ht");
  sprintf(title,"Momentum Transfer");
  ht = new TH2F(hname,title,100,-2.5,0.0,2,-0.5,1.5);

  sprintf(hname,"hIM_cuts");
  sprintf(title,"IM by cuts");
  hIM_cuts = new TH2F(hname,title,nIM,IMl,IMh,6,0.5,6.5);

  Mass[1] = new TH1F("Mass1","Mass1",1000,0,2);
  Angle = new TH1F("Angle","Angle",360,-10,360);


  for(j=0;j<=36;j++){
    for(jj=0;jj<=4;jj++){     
       if(jj==0) target = "all";
       else if(jj==1) target = "2H";
       else if(jj==2) target = "C";
       else if(jj==3) target = "FeTi";
       else if(jj==4) target = "Pb";  

        sprintf(hname,"AngleCutAroundTgt_%s_range_%d",target,j);
        sprintf(title,"Angle Cut Around %s angle range %d",target,j);
        AngleCutAroundTarget[jj][j] = new TH1F(hname,title,10000, 0, 2);

        sprintf(hname,"MassCutAroundTgt_%s_range_%d",target,j);
        sprintf(title,"Mass Cut Around %s  angle range %d",target,j);
        MassCutAroundTargetAngle[jj][j] = new TH1F(hname,title,10000, 0, 2);

        sprintf(hname,"MassCutAfterTgt_%s_range_%d",target,j);
        sprintf(title,"Mass Cut After %s  angle range %d",target,j);
        MassCutAfterTargetAngle[jj][j] = new TH1F(hname,title,10000, 0, 2);

        sprintf(hname,"Carbon%d_AngleCutAround_range_%d",jj,j);
        sprintf(title,"Angle Cut Around %s angle range %d",target,j);
        AngleCutAroundTargetCarbon[jj][j] = new TH1F(hname,title,10000, 0, 2);}
  }


  sprintf(hname,"hVz");
  sprintf(title,"z vertex");
  hVz = new TH2F(hname,title,nVz,zlo,zhi,2,-0.5,1.5);
  

   for(k=1;k<=6;k++){

     if(k==1) {between = "0-0.35"; range = 0.0; range2 = 0.35;}
     else if(k==2) {between = "0.35-0.70"; range = 0.35; range2 = 0.70;}
     else if(k==3) {between = "0.70-1.05"; range = 0.70; range2 = 1.05;} 
     else if(k==4) {between = "1.05-1.40"; range = 1.05; range2 = 1.40;} 
     else if(k==5) {between = "1.40-1.75"; range = 1.40; range2 = 1.75;}
     else if(k==6) {between = "1.75-2.10"; range = 1.75; range2 = 2.10;}
    
     sprintf(hname,"MassCutAround%i",k);
     sprintf(title,"Mass Cut Around Target Momentums between %s",between);
     MassCutAround[k] = new TH1F(hname,title,10000, 0, 2);

     sprintf(hname,"MassCutAfter%i",k);
     sprintf(title,"Mass Cut After Target Momentums between %s",between);
     MassCutAfter[k] = new TH1F(hname,title,10000, 0, 2);
        
     sprintf(hname,"MomentumCutAfter%i",k);
     sprintf(title,"Momentum Cut After Target Momentums between %s",between);
     MomentumCutAfter[k] = new TH1F(hname,title,1000, range, range2);

     sprintf(hname,"MomentumCutAngle%i",k);
     sprintf(title,"Momentum Angle Cut Target Momentums between %s",between);
     MomentumCutAngle[k] = new TH1F(hname,title,1000, range, range2);
        
     sprintf(hname,"MomentumCutAround%i",k);
     sprintf(title,"Momentum Cut Around Target Momentums between %s",between);
     MomentumCutAround[k] = new TH1F(hname,title,1000, range, range2);
        
     for (i=1;i<=4;i++){
        
       if(i==1) target = "2H";
       else if(i==2) target = "C";
       else if(i==3) target = "FeTi";
       else if(i==4) target = "Pb";       	 
        	 
       sprintf(hname,"MassCutAroundTgt_%s_%i",target,k);
       sprintf(title,"Mass Cut Around Target %s Momentums between %s",target,between);
       MassCutAroundTgt[i][k] = new TH1F(hname,title,10000, 0, 2);

       sprintf(hname,"MassCutAngleTgt_%s_%i",target,k);
       sprintf(title,"Mass Cut Around Target Angle method %s Momentums between %s",target,between);
       MassCutAngleTgt[i][k] = new TH1F(hname,title,10000, 0, 2);
        
       sprintf(hname,"MassCutAfterTgt_%s_%i",target,k);
       sprintf(title,"Mass Cut After Target %s Momentums between %s",target,between);
       MassCutAfterTgt[i][k] = new TH1F(hname,title,10000, 0, 2);}
     }

 
  for(i=0;i<nType;i++){  
    sprintf(hname,"hIM_tgt%i",i);
    sprintf(title,"IM by Target, %s",VzAna[i]);
    hIM_tgt[i] = new TH2F(hname,title,nIM,IMl,IMh,tgtTot,-0.5,tgtTot-0.5);   
    sprintf(hname,"hPcomp%i",i);
    sprintf(title,"Momentum Comparison, %s",VzAna[i]);
    hPcomp[i] = new TH2F(hname,title,nP,Pl,Ph,nP,Pl,Ph);
    for(j=0;j<3;j++){
      sprintf(hname,"hIM_P%i_%i",i,j);
      sprintf(title,"IM vs Momentum of %s, %s",Plabel[j],VzAna[i]);
      hIM_P[i][j] = new TH2F(hname,title,nIM,IMl,IMh,nP,Pl,Ph);
    }
  }
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


void PolSubtractionAngle(char *rootFile, char* yieldFile){
	char *tar;
	char hname[50];
	char hname2[50];
	char name[50];
	char AAA[10];
	Double_t ratio[37];
	Double_t Chi[37];
	Double_t yields[4][4][37];
	Int_t mass;
	Int_t Xbins;
	Double_t NDF;

	Double_t R;
	Int_t A[4] = {2,12,52,208};
	Double_t d[4] = {1.015,0.88,0.933,1.094};
	Double_t Aw[4] = {2.0140,21.011,51.856,207.298};
	
	Float_t Lmar = 0.125; // set the left margin
	Float_t Rmar = 0.125; // set the right margin
	Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values
	
	Float_t SumLo, SumHi;
	Float_t xLo = 0.41;  // lower value of x-axis for drawing histogram
	Float_t xHi = 0.6;   // upper value of x-axis for drawing histogram
	Float_t PeakLo = 0.49; // lower limit on the peak range
	Float_t PeakHi = 0.51; // upper limit on the peak range
	
	Int_t ibg = 1; // parameter index in the par[] array where the background parameters start

	Int_t i, meth_c, tar_c, k;
	Int_t count = 36;
	Float_t x, pval;
	
	TH1F *hist; // original histogram
	TH1F *hBgFit; // histogram of background
	TH1F *Peak;  //histogram of peak
	TH1F *YieldHist[8];
	TH1F *ChiHist[8];
	TH1F *RatioHist[8];
	TH1F *Ratio[8];
	TH1F *Stuff[8];

		
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

	for(meth_c=0; meth_c<=3; meth_c++){		//for loop for differnt methods (angle, around, after)

	   for(tar_c=0;tar_c<=3;tar_c++){		//for loop for different targets

	      switch (tar_c){ 
	         case 0: tar = "2H";
	                 break;
	         case 1: tar = "C";
	                 break;
	         case 2: tar = "FeTi";
	                 break;
	         case 3: tar = "Pb";
	                 break;
	      }





	for(i=0;i<=count;i++){		//for loop for difffernt angle ranges

	   switch (meth_c){ //switch for naming histograms with differnt target methods
	      case 0: sprintf(hname, "AngleCutAroundTgt_%s_range_%d", tar,i);
		            sprintf(name, "AngleCutAroundTgt_%s_", tar);
	              sprintf(AAA,"Angle"); 
	              break;
	      case 1: sprintf(hname, "MassCutAroundTgt_%s_range_%d",tar,i);
		            sprintf(name, "MassCutAroundTgt_%s",tar);
	              sprintf(AAA,"Around"); 
	              break;
	      case 2: sprintf(hname, "MassCutAfterTgt_%s_range_%d",tar,i);
		            sprintf(name, "MassCutAfterTgt_%s",tar);
	              sprintf(AAA,"After"); 
	              break;
	      case 3: sprintf(hname, "Carbon%d_AngleCutAround_range_%d",tar_c+1,i);
		            sprintf(name, "Carbon%d_AngleCutAround",tar_c+1);
	              sprintf(AAA,"Carbon");
	              break;}


	   hist = (TH1F*)fd->Get(hname); 			// get the histogram from the ROOT file
					
	   // histogram of background
	   hBgFit = (TH1F*)hist->Clone("hBgFit"); // clone original hist. into background temp. hist.
	   hBgFit->SetName("hBgFit");
	   hBgFit->SetTitle("Background");

	   Peak = (TH1F*)hist->Clone("Peak");
	   Peak->SetName("Peak");
	   Peak->SetTitle("Peak");

	   SumLo = 0.486;
	   SumHi = 0.507;			
	   Xbins = hBgFit->GetNbinsX();


	   for(k=1; k<Xbins; k++){
	      if (hBgFit->GetBinCenter(k)<=SumHi && hBgFit->GetBinCenter(k)>=SumLo){ 
	         hBgFit->SetBinContent(k,0.0);//Set all bins inside peak range to be 0
	      }
	   }

	   TF1 *pol = new TF1("pol",polFit,xLo,xHi,PolNum);

	   hBgFit->Fit("pol","R+");         // fit the background
	   pol->GetParameters(&par[ibg]);
	   pol->SetParameters(&par[ibg]); 	// set the pfinal parameters for the background function
	   NDF = pol->GetNDF();
	   Chi[i] = pol->GetChisquare();

	   if (NDF) Chi[i] = Chi[i]/NDF;
	   else Chi[i] = 0;

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

	      yields[tar_c][meth_c][i] = Peak->Integral(Peak->FindBin(SumLo),Peak->FindBin(SumHi)); //sum total counts in peak


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
	      sprintf(OutCan,"pol/%s_%s_%s_%d.gif",yieldFile,tar,AAA,i);
	      can1->Print(OutCan);
	      sprintf(OutCan,"pol/%s_%s_%s_%d.eps",yieldFile,tar,AAA,i);
	      can1->Print(OutCan);

	      // open text file for the yields
	      char OutFile[100];
	      sprintf(OutFile,"pol/%s_%s_%s_%d.yld",yieldFile,tar,AAA,i);
	      ofstream fout(OutFile);

	      if(yields[tar_c][meth_c][i] >= 1 && yields[0][meth_c][i] >= 1) ratio[i] = yields[tar_c][meth_c][i]/yields[0][meth_c][i];
	      if(yields[tar_c][meth_c][i] <1 && yields[0][meth_c][i]<1) ratio[i] = 0;






	   }//end for i range


           //sprintf(hname,"YieldHist_%s_%s",tar,AAA);
	   sprintf(hname2,"Yield_%s",name);
           sprintf(title,"Yield around %s %s for different #circ 's",tar,AAA);
	   YieldHist[tar_c] = new TH1F(hname2,title, 380, 0, 190 );
	   YieldHist[tar_c]->GetYaxis()->SetTitle("counts");
	   YieldHist[tar_c]->GetYaxis()->CenterTitle();
	   YieldHist[tar_c]->GetXaxis()->SetTitle("Angle");
	   YieldHist[tar_c]->GetXaxis()->CenterTitle();


           //sprintf(hname,"ChiHist_%s_%s",tar,AAA);
	   sprintf(hname2,"Chi_%s",name);
           sprintf(title,"Chi-Square for target %s, %s",tar,AAA);
	   ChiHist[tar_c] = new TH1F(hname2,title, 380, 0, 190 );
	   ChiHist[tar_c]->GetYaxis()->SetTitle("Chi-Square Value");
	   ChiHist[tar_c]->GetYaxis()->CenterTitle();
	   ChiHist[tar_c]->GetXaxis()->SetTitle("Angle");
	   ChiHist[tar_c]->GetXaxis()->CenterTitle();

	   char OutFile[100];
	   sprintf(OutFile,"pol/%s_%s_%s.yld",yieldFile,tar,AAA);
	   ofstream fout(OutFile);



	   for(i=0;i<=count;i++){
	      ChiHist[tar_c]->SetBinContent((i*10)+2, Chi[i]);
	      YieldHist[tar_c]->SetBinContent((i*10)+2, yields[tar_c][meth_c][i]);
	      fout<<i<<"\t\t"<<yields[tar_c][meth_c][i]<<"\t\t"<<yields[0][meth_c][i]<<"\t\t"<<ratio[i]<<"\t\t"<<Chi[i]<<endl;}

           //sprintf(hname,"RatioHist_%s_%s",tar,AAA);
	   sprintf(hname2,"Ratio_%s",name);
           sprintf(title,"Ratio for %s %s target for differnt #circ 's",tar,AAA);
	   RatioHist[tar_c] = new TH1F(hname2,title, 380, 0, 190 ); 
	   RatioHist[tar_c]->GetYaxis()->SetTitle("ratio");
	   RatioHist[tar_c]->GetYaxis()->CenterTitle();
	   RatioHist[tar_c]->GetXaxis()->SetTitle("Angle");
	   RatioHist[tar_c]->GetXaxis()->CenterTitle();


	   for(i=0;i<=count;i++){
	      RatioHist[tar_c]->SetBinContent((i*10)+2, ratio[i]);}

	   RatioHist[tar_c]->Scale((A[0]*d[0]*Aw[tar_c])/(A[tar_c]*d[tar_c]*Aw[0]));	
	   RatioHist[tar_c]->Write();
           YieldHist[tar_c]->Write();
	   ChiHist[tar_c]->Write();

	}//end for tar_c targets

        //sprintf(hname,"%s_Ratio",AAA);
	sprintf(hname2,"Ratio_%s",AAA);
        sprintf(title,"Ratio %s",AAA);
	Ratio[meth_c] = new TH1F(hname2,title, 210, 0, 210 );
	Ratio[meth_c]->GetYaxis()->SetTitle("Ratio");
	Ratio[meth_c]->GetYaxis()->CenterTitle();
	Ratio[meth_c]->GetXaxis()->SetTitle("Mass");
	Ratio[meth_c]->GetXaxis()->CenterTitle();

        //sprintf(hname,"%s_Yield",AAA);
	sprintf(hname2,"Yield_%s",AAA);
        sprintf(title,"Yield %s",AAA);
	Stuff[meth_c] = new TH1F(hname2,title, 210, 0, 210 );
	Stuff[meth_c]->GetYaxis()->SetTitle("Counts");
	Stuff[meth_c]->GetYaxis()->CenterTitle();
	Stuff[meth_c]->GetXaxis()->SetTitle("Mass");
	Stuff[meth_c]->GetXaxis()->CenterTitle();


	for(tar_c=0;tar_c<=3;tar_c++){
	   R = (A[0]*d[0]*Aw[tar_c]*yields[tar_c][meth_c][0])/(A[tar_c]*d[tar_c]*Aw[0]*yields[0][meth_c][0]);
	   if(meth_c == 3) R = (yields[tar_c][meth_c][0]/yields[0][meth_c][0]);
	   Ratio[meth_c]->SetBinContent(A[tar_c], R);
	   Stuff[meth_c]->SetBinContent(A[tar_c], yields[tar_c][meth_c][0]);}

	Ratio[meth_c]->Write();
	Stuff[meth_c]->Write();


	}//end for meth_c method

	   gPad->SetLogy();
	   YieldHist[0]->GetXaxis()->SetTitle("angle in degrees");
	   YieldHist[0]->GetXaxis()->CenterTitle();
	   YieldHist[0]->SetTitle("Ratio of Yeilds to 2H");  	   
  	   YieldHist[3]->SetMarkerSize(1.25);
  	   YieldHist[3]->SetMarkerStyle(22);
	   YieldHist[3]->SetMarkerColor(5);
	   YieldHist[3]->Draw("p");
  	   YieldHist[2]->SetMarkerSize(1.25);
  	   YieldHist[2]->SetMarkerStyle(22);
	   YieldHist[2]->SetMarkerColor(4);  	   
	   YieldHist[2]->Draw("p""same");
  	   YieldHist[1]->SetMarkerSize(1.25);
  	   YieldHist[1]->SetMarkerStyle(22);
	   YieldHist[1]->SetMarkerColor(3);  	   
	   YieldHist[1]->Draw("p""same");
  	   YieldHist[0]->SetMarkerSize(1.25);
  	   YieldHist[0]->SetMarkerStyle(22);
  	   YieldHist[0]->SetMarkerColor(2);
	   YieldHist[0]->Draw("p""same");

	   char OutCan[100];
	   sprintf(OutCan,"pol/Yield%s.gif",yieldFile);
	   can1->Print(OutCan);
	   sprintf(OutCan,"pol/Yield%s.eps",yieldFile);
	   can1->Print(OutCan);

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

	   //char OutCan[100];
	   sprintf(OutCan,"pol/Ratio%s.gif",yieldFile);
	   can1->Print(OutCan);
	   sprintf(OutCan,"pol/Ratio%s.eps",yieldFile);
	   can1->Print(OutCan);
	   fd->Close();

}


void MomPolSubtraction(char *rootFile, char* yieldFile){
	char *tar;
	char hname[50];
	char hname2[50];
	char name[50];
	char AAA[10];
	Double_t ratio[37];
	Double_t Chi[37];
	Double_t yields[4][4][37];
	Int_t mass;
	Int_t Xbins;
	Double_t NDF;

	Double_t R;
	Int_t A[4] = {2,12,52,208};
	Double_t d[4] = {1.015,0.88,0.933,1.094};
	Double_t Aw[4] = {2.0140,21.011,51.856,207.298};
	
	Float_t Lmar = 0.125; // set the left margin
	Float_t Rmar = 0.125; // set the right margin
	Float_t yoff = 1.75;  // set the offset between the y-axis label and the axis values
	
	Float_t SumLo, SumHi;
	Float_t xLo = 0.41;  // lower value of x-axis for drawing histogram
	Float_t xHi = 0.6;   // upper value of x-axis for drawing histogram
	Float_t PeakLo = 0.49; // lower limit on the peak range
	Float_t PeakHi = 0.51; // upper limit on the peak range
	
	Int_t ibg = 1; // parameter index in the par[] array where the background parameters start

	Int_t i, meth_c, tar_c, k;
	Int_t count = 6;
	Float_t x, pval;
	
	TH1F *Momhist; // original histogram
	TH1F *MomhBgFit; // histogram of background
	TH1F *MomPeak;  //histogram of peak
	TH1F *MomYieldHist[8];
	TH1F *MomChiHist[8];
	TH1F *MomRatioHist[8];
	TH1F *MomRatio[8];
	TH1F *MomStuff[8];

		
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

	for(meth_c=0; meth_c<=2; meth_c++){		//for loop for differnt methods (angle, around, after)

	   for(tar_c=0;tar_c<=3;tar_c++){		//for loop for different targets

	      switch (tar_c){ 
	         case 0: tar = "2H";
	                 break;
	         case 1: tar = "C";
	                 break;
	         case 2: tar = "FeTi";
	                 break;
	         case 3: tar = "Pb";
	                 break;
	      }


	for(i=1;i<=count;i++){		//for loop for difffernt momentum ranges

	   switch (meth_c){ //switch for naming histograms with differnt target methods
	      case 0: sprintf(AAA,"Angle"); 
	              break;
	      case 1: sprintf(AAA,"Around"); 
	              break;
	      case 2: sprintf(AAA,"After"); 
	              break;}

	sprintf(hname, "MassCut%sTgt_%s_%d",AAA,tar,i);
	sprintf(name, "MassCut%sTgt_%s",AAA,tar);


	   Momhist = (TH1F*)fd->Get(hname); 			// get the histogram from the ROOT file
					
	   // histogram of background
	   MomhBgFit = (TH1F*)Momhist->Clone("hBgFit"); // clone original hist. into background temp. hist.
	   MomhBgFit->SetName("hBgFit");
	   MomhBgFit->SetTitle("Background");

	   MomPeak = (TH1F*)Momhist->Clone("Peak");
	   MomPeak->SetName("Peak");
	   MomPeak->SetTitle("Peak");

	   SumLo = 0.486;
	   SumHi = 0.507;			
	   Xbins = MomhBgFit->GetNbinsX();


	   for(k=1; k<Xbins; k++){
	      if (MomhBgFit->GetBinCenter(k)<=SumHi && MomhBgFit->GetBinCenter(k)>=SumLo){ 
	         MomhBgFit->SetBinContent(k,0.0);//Set all bins inside peak range to be 0
	      }
	   }

	   TF1 *pol = new TF1("pol",polFit,xLo,xHi,PolNum);

	   MomhBgFit->Fit("pol","R+");         // fit the background
	   pol->GetParameters(&par[ibg]);
	   pol->SetParameters(&par[ibg]); 	// set the pfinal parameters for the background function
	   NDF = pol->GetNDF();
	   Chi[i] = pol->GetChisquare();

	   if (NDF) Chi[i] = Chi[i]/NDF;
	   else Chi[i] = 0;

	   MomPeak->Add(pol,-1.0); 		// subtract background function

	   for(k=1; k<MomPeak->GetNbinsX(); k++){
		x = MomPeak->GetBinCenter(k); 		// read bin center value on the x-axis
		if(x>=SumLo && x<=SumHi){ 		// check that x is in the peak summation region
			pval = MomPeak->GetBinContent(k);	// get the number of counts in the bin
			if(pval<0.0) pval = 0.0; 	// if neg. counts, set to zero
		}
		else{
			pval =0.0;
		}
		MomPeak->SetBinContent(k,pval); // refill the histogram
	      }

	      yields[tar_c][meth_c][i] = MomPeak->Integral(MomPeak->FindBin(SumLo),MomPeak->FindBin(SumHi)); //sum total counts in peak


	      // set up the Pad parameters
	      gPad->SetLeftMargin(Lmar);
	      gPad->SetRightMargin(Rmar);
	      gPad->SetFillColor(0);

	      // draw the original histogram
	      Momhist->SetTitle(0);
	      Momhist->GetXaxis()->SetTitle(xtitle);
	      Momhist->GetXaxis()->CenterTitle();
	      Momhist->GetYaxis()->SetTitle(ytitle);
	      Momhist->GetYaxis()->CenterTitle();
	      Momhist->GetYaxis()->SetTitleOffset(yoff);
	      Momhist->GetXaxis()->SetRangeUser(xLo,xHi);  // set the x-axis range for the plot
	      Momhist->SetLineWidth(1);
	      Momhist->SetMinimum(0); // start the y-axis at zero.
	      Momhist->Draw();

	      // draw the polinomial fit
	      pol->Draw("same");


	      //draw peak
	      MomPeak->SetFillColor(4);
	      MomPeak->Draw("same");
	
	
	      // create the image files
	      char OutCan[100];
	      sprintf(OutCan,"pol/%s_%s_%s_%d.gif",yieldFile,tar,AAA,i);
	      can1->Print(OutCan);
	      sprintf(OutCan,"pol/%s_%s_%s_%d.eps",yieldFile,tar,AAA,i);
	      can1->Print(OutCan);


	      // open text file for the yields
	      char OutFile[100];
	      sprintf(OutFile,"pol/%s_%s_%s_%d.yld",yieldFile,tar,AAA,i);
	      ofstream fout(OutFile);


	      if(yields[tar_c][meth_c][i] >= 1 && yields[0][meth_c][i] >= 1) ratio[i] = yields[tar_c][meth_c][i]/yields[0][meth_c][i];

	      if(yields[tar_c][meth_c][i] <1 && yields[0][meth_c][i]<1)ratio[i] = 0;






	   }//end for i range


           //sprintf(hname,"YieldHist_%s_%s",tar,AAA);
	   sprintf(hname2,"Yield_mom_%s",name,tar_c,AAA);
     sprintf(title,"Yield around %s %s for different momentum ranges",tar,AAA);
	   MomYieldHist[tar_c] = new TH1F(hname2,title, 30, 0, 15 );
	   MomYieldHist[tar_c]->GetYaxis()->SetTitle("counts");
	   MomYieldHist[tar_c]->GetYaxis()->CenterTitle();
	   MomYieldHist[tar_c]->GetXaxis()->SetTitle("momentum range");
	   MomYieldHist[tar_c]->GetXaxis()->CenterTitle();


           //sprintf(hname,"ChiHist_%s_%s",tar,AAA);
	   sprintf(hname2,"Chi_mom_%s",name,tar_c);
     sprintf(title,"Chi-Square for target %s, %s",tar,AAA);
	   MomChiHist[tar_c] = new TH1F(hname2,title, 30, 0, 15 );
	   MomChiHist[tar_c]->GetYaxis()->SetTitle("Chi-Square Value");
	   MomChiHist[tar_c]->GetYaxis()->CenterTitle();
	   MomChiHist[tar_c]->GetXaxis()->SetTitle("momentum range");
	   MomChiHist[tar_c]->GetXaxis()->CenterTitle();

	   char OutFile[100];
	   sprintf(OutFile,"pol/%s_%s_%s.yld",yieldFile,tar,AAA);
	   ofstream fout(OutFile);



	   for(i=0;i<=count;i++){
	      MomChiHist[tar_c]->SetBinContent((i*4)+1, Chi[i]);
	      MomYieldHist[tar_c]->SetBinContent((i*4)+1, yields[tar_c][meth_c][i]);
	      fout<<i<<"\t\t"<<yields[tar_c][meth_c][i]<<"\t\t"<<yields[0][meth_c][i]<<"\t\t"<<ratio[i]<<"\t\t"<<Chi[i]<<endl;}

           //sprintf(hname,"RatioHist_%s_%s",tar,AAA);
	   sprintf(hname2,"Ratio_mom_%s",name);
           sprintf(title,"Ratio for %s %s target for differnt momentum ranges",tar,AAA);
	   MomRatioHist[tar_c] = new TH1F(hname2,title, 30, 0, 15 ); 
	   MomRatioHist[tar_c]->GetYaxis()->SetTitle("ratio");
	   MomRatioHist[tar_c]->GetYaxis()->CenterTitle();
	   MomRatioHist[tar_c]->GetXaxis()->SetTitle("momentum range");
	   MomRatioHist[tar_c]->GetXaxis()->CenterTitle();


	   for(i=0;i<=count;i++) MomRatioHist[tar_c]->SetBinContent((i*4)+1, ratio[i]); 

	   MomRatioHist[tar_c]->Scale((A[0]*d[0]*Aw[tar_c])/(A[tar_c]*d[tar_c]*Aw[0]));	
	   MomRatioHist[tar_c]->Write();
     MomYieldHist[tar_c]->Write();
	   MomChiHist[tar_c]->Write();

	}//end for tar_c targets

	sprintf(hname2,"%s_Ratio",AAA);
  sprintf(title,"Ratio %s",AAA);
	MomRatio[meth_c] = new TH1F(hname2,title, 420, 0, 210 );
	MomRatio[meth_c]->GetYaxis()->SetTitle("Ratio");
	MomRatio[meth_c]->GetYaxis()->CenterTitle();
	MomRatio[meth_c]->GetXaxis()->SetTitle("Mass");
	MomRatio[meth_c]->GetXaxis()->CenterTitle();

	sprintf(hname2,"%s_Yield",AAA);
  sprintf(title,"Yield %s",AAA);
	MomStuff[meth_c] = new TH1F(hname2,title, 420, 0, 210 );
	MomStuff[meth_c]->GetYaxis()->SetTitle("Counts");
	MomStuff[meth_c]->GetYaxis()->CenterTitle();
	MomStuff[meth_c]->GetXaxis()->SetTitle("Mass");
	MomStuff[meth_c]->GetXaxis()->CenterTitle();

/*
	for(tar_c=0;tar_c<=3;tar_c++){
	   R = (A[0]*d[0]*Aw[tar_c]*yields[tar_c][meth_c][0])/(A[tar_c]*d[tar_c]*Aw[0]*yields[0][meth_c][0]);
	   MomRatio[meth_c]->SetBinContent(A[tar_c], R);
	   MomStuff[meth_c]->SetBinContent(A[tar_c], yields[tar_c][meth_c][0]);}

	//MomRatio[meth_c]->Write();
	//MomStuff[meth_c]->Write();*/

	}//end for meth_c method

	   gPad->SetLogy();
	   MomYieldHist[0]->GetXaxis()->SetTitle("angle in degrees");
	   MomYieldHist[0]->GetXaxis()->CenterTitle();
	   MomYieldHist[0]->SetTitle("Ratio of Yeilds to 2H");  	   
  	 MomYieldHist[3]->SetMarkerSize(1.25);
  	 MomYieldHist[3]->SetMarkerStyle(22);
	   MomYieldHist[3]->SetMarkerColor(5);
	   MomYieldHist[3]->Draw("p");
  	 MomYieldHist[2]->SetMarkerSize(1.25);
  	 MomYieldHist[2]->SetMarkerStyle(22);
	   MomYieldHist[2]->SetMarkerColor(4);  	   
	   MomYieldHist[2]->Draw("p""same");
  	 MomYieldHist[1]->SetMarkerSize(1.25);
  	 MomYieldHist[1]->SetMarkerStyle(22);
	   MomYieldHist[1]->SetMarkerColor(3);  	   
	   MomYieldHist[1]->Draw("p""same");
  	 MomYieldHist[0]->SetMarkerSize(1.25);
  	 MomYieldHist[0]->SetMarkerStyle(22);
  	 MomYieldHist[0]->SetMarkerColor(2);
	   MomYieldHist[0]->Draw("p""same");

	   char OutCan[100];
	   sprintf(OutCan,"pol/Yield%s.gif",yieldFile);
	   can1->Print(OutCan);
	   sprintf(OutCan,"pol/Yield%s.eps",yieldFile);
	   can1->Print(OutCan);

	   gPad->SetLogy();
	   MomRatioHist[0]->GetXaxis()->SetTitle("angle in degrees");
	   MomRatioHist[0]->GetXaxis()->CenterTitle();
	   MomRatioHist[0]->SetTitle("Ratio of Yeilds to 2H");  	   
  	 MomRatioHist[3]->SetMarkerSize(1.25);
  	 MomRatioHist[3]->SetMarkerStyle(22);
	   MomRatioHist[3]->SetMarkerColor(5);
	   MomRatioHist[3]->Draw("p");
  	 MomRatioHist[2]->SetMarkerSize(1.25);
  	 MomRatioHist[2]->SetMarkerStyle(22);
	   MomRatioHist[2]->SetMarkerColor(4);  	   
	   MomRatioHist[2]->Draw("p""same");
  	 MomRatioHist[1]->SetMarkerSize(1.25);
  	 MomRatioHist[1]->SetMarkerStyle(22);
	   MomRatioHist[1]->SetMarkerColor(3);  	   
	   MomRatioHist[1]->Draw("p""same");
  	 MomRatioHist[0]->SetMarkerSize(1.25);
  	 MomRatioHist[0]->SetMarkerStyle(22);
  	 MomRatioHist[0]->SetMarkerColor(2);
	   MomRatioHist[0]->Draw("p""same");

	   //char OutCan[100];
	   sprintf(OutCan,"pol/Ratio%s.gif",yieldFile);
	   can1->Print(OutCan);
	   sprintf(OutCan,"pol/Ratio%s.eps",yieldFile);
	   can1->Print(OutCan);
	   fd->Close();

}      
