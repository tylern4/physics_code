/************************************************************************/
/*  AnaTwoPion.cc                                                       */
/*                                                                      */
/*  Created by Nick Tyler, Canisius College                             */
/*                                                                      */
/************************************************************************/



#include "TTree.h"
#include "TROOT.h"
#include "TH2.h"
#include <TLorentzVector.h>
#include <TFile.h>
#include <iostream>
#include "TStyle.h"
#include "TCanvas.h"
#include "TObject.h"
#include <stdio.h>
#include <fstream>
#define PI 3.14159265;

using namespace std;

  const Int_t nTgt = 6;
  const Int_t nFoils = 8;
  const Int_t tgtTot = nTgt + nFoils;
  char *target[nTgt] = {"All","2H","C","Fe-Ti","Pb","No Pb"};
  char *Foil[nFoils] = {"2H","C1","Fe","C2","Pb","C3","Ti","C4"};

  Float_t zcut[nFoils+1] ={-22.0,-14.0,-11.0,-8.5,-6.0,-3.5,-1.0,1.5,4.0};

  Float_t zhi = 5.0;

  char hname[50];
  char htitle[500];
  char title[500];
  char cname[50];
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
  TH1F *MomentumCutAround[7];
  TH1F *MassCutAroundTgt[5][7];
  TH1F *MassCutAfterTgt[5][7];

  TH1F *MassCutAfterTarget[5]; 
  TH1F *MassCutAroundTarget[5];

  TH1F *AngleCutAroundTarget[37];

  TH1F *hYldVsA;
  TH1F *hYldVsR;

TH1F *Angle;

  struct labels {
    char *x;
    char *y;
    char *hPrefix;
    char *fPrefix;};



  Double_t zFoil_mean[nFoils] = {-18.0,-12.3,-9.89,-7.42,-4.93,-2.46,0.02,2.44}; // foil peak mean
  Double_t zFoil_sigma[nFoils] = {1.17,0.258,0.249,0.252,0.235,0.234,0.245,0.231}; // foil peak sigma



  //declarations of functions

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
     Double_t zdiff[nFoils];  // z difference between target foil and decay vertex
     Double_t rdiff[nFoils];  // radial difference between target foil and decay vertex
     Double_t decayL[nFoils]; // decay length between target foil and decay vertex
     Double_t ctau[nFoils];   // c*tau for the decay length
     Double_t rtgt[nFoils];   // parent particle r for each target foil position
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
    ParticleTarget temp; // create the ParticleTarget struct

    Double_t theta = V4.Theta();
    Double_t beta = V4.Beta();


    temp.theta = (theta*180)/PI;
 

    memset(&temp.rtgt,-99,sizeof(temp.rtgt));		// initialize rtgt array elements to -99
    memset(&temp.rdiff,-99,sizeof(temp.rdiff)); 	// initialize rdiff array elements to -99
    memset(&temp.zdiff,-99,sizeof(temp.zdiff)); 	// initialize zdiff array elements to -99
    memset(&temp.decayL,-99,sizeof(temp.decayL)); 	// initialize decayL array elements to -99
    memset(&temp.ctau,-99,sizeof(temp.ctau)); 		// initialize ctau array elements to -99

    temp.nTgtTrav = 0; 					// initialize coutner to zero
    temp.beamlineCrossing = 0; 				// initialize beamline crossing flag


    for(i=1;i<nFoils;i++){
       if(vz>zFoil_mean[i]){
          temp.zdiff[i] = vz - zFoil_mean[i]; // z difference between decay vertex and target position      
          temp.rdiff[i] = temp.zdiff[i]*tan(theta); // radial difference between decay vertex and target
          temp.rtgt[i] = vr - temp.rdiff[i]; // radial distance of parent particle in target foil
          temp.decayL[i] = sqrt(temp.zdiff[i]*temp.zdiff[i] + temp.rdiff[i]*temp.rdiff[i]); // decay length
          if(beta) temp.ctau[i] = temp.decayL[i]/beta; // c*tau for ith target foil


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
  Int_t i, j, k, ii;
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
  
  TVector3 V3[nPart];        // vertex 3-vectors
  TLorentzVector V4[nPart];  // particle 4-vectors
  TLorentzVector pairV4;     // particle pair 4-vector
  TLorentzVector Beam;       // beam 4-vector
  TLorentzVector A;          // beam-pair 4-vector


Double_t theta;
ParticleTarget TgtInfo;
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


	   TgtInfo = GetParticleTargetInfo(pZ[ii], sqrt(pX[ii]*pX[ii]+pY[ii]*pY[ii]), pairV4);

	   if (pairM>=PeakLo && pairM<=PeakHi) Angle->Fill(TgtInfo.theta);// Angle between the kaon and z axis for mass range


	   for(k=0; k<=36; k++){
	      if(TgtInfo.theta >= k*5) AngleCutAroundTarget[k]->Fill(pairM);}




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




/*	    //create the image files
	    char OutCan[200];
	    sprintf(OutCan,"pictures/Sideband%s.gif",inName);
	    can1->Print(OutCan);
	   
	    sprintf(OutCan,"results/%s.eps",name);
	    can1->Print(OutCan);
*/
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
	
	AnaTwoPion(infilename,outfilename);
	Projection(outfilename,outfilename);
	BgSub_Sideband(outfilename,"test");
     
	     

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
	
	Double_t zFoil_mean[8] = {-18.0,-12.3,-9.89,-7.42,-4.93,-2.46,0.02,2.44};      // foil peak mean
     Double_t zFoil_sigma[8] = {1.17,0.258,0.249,0.252,0.235,0.234,0.245,0.231};    // foil peak sigma
  
	Float_t z1, z2;

	for(i=0;i<nFoils;i++){
    	   z1 = (i==0) ? zcut[i] : zFoil_mean[i] + 3.0*zFoil_sigma[i];    
        //    z1 = zFoil_mean[i] + 3.0*zFoil_sigma[i];    
        z2 = (i+1<nFoils) ? zFoil_mean[i+1] - 3.0*zFoil_sigma[i+1] : zhi;
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
  
  	for(i=0;i<nFoils;i++){
        if(vz>=zcut[i] && vz<zcut[i+1]) ret = i+1;
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

  Int_t i=0, j=0 ,k=0 , l=0;
  
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

for (j=0;j<=36;j++){
     AngleCutAroundTarget[j]->SetXTitle("Mass (GeV)");
     AngleCutAroundTarget[j]->GetXaxis()->CenterTitle();
     AngleCutAroundTarget[j]->Write();}
  
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
  
    MomentumCutAround[i]->SetXTitle("Momentum (GeV)");
    MomentumCutAround[i]->GetXaxis()->CenterTitle();
    MomentumCutAround[i]->Write();
  
    for(k=1;k<=4;k++){  
       MassCutAroundTgt[k][i]->SetXTitle("Mass (GeV)");
       MassCutAroundTgt[k][i]->GetXaxis()->CenterTitle();
       MassCutAroundTgt[k][i]->Write();
         
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

  Int_t i, j, k;
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

Angle = new TH1F("Angle","Angle",360,-10,360);

for(j=0;j<=36;j++){
      sprintf(hname,"AngleCutAfterTgt_range_%d", j);
      sprintf(title,"AngleCutAfterTgt_range_%d",j);
      AngleCutAroundTarget[j] = new TH1F(hname,title,10000, 0, 2);}
  
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
     MassCutAround[k] = new TH1F(hname,title,400, 0, 2);

     sprintf(hname,"MassCutAfter%i",k);
     sprintf(title,"Mass Cut After Target Momentums between %s",between);
     MassCutAfter[k] = new TH1F(hname,title,400, 0, 2);
        
     sprintf(hname,"MomentumCutAfter%i",k);
     sprintf(title,"Momentum Cut After Target Momentums between %s",between);
     MomentumCutAfter[k] = new TH1F(hname,title,60, range, range2);
        
     sprintf(hname,"MomentumCutAround%i",k);
     sprintf(title,"Momentum Cut Around Target Momentums between %s",between);
     MomentumCutAround[k] = new TH1F(hname,title,60, range, range2);
        
     for (i=1;i<=4;i++){
        
       if(i==1) target = "2H";
       else if(i==2) target = "C";
       else if(i==3) target = "FeTi";
       else if(i==4) target = "Pb";       	 
        	 
       sprintf(hname,"MassCutAroundTgt_%s_%i",target,k);
       sprintf(title,"Mass Cut Around Target %s Momentums between %s",target,between);
       MassCutAroundTgt[i][k] = new TH1F(hname,title,400, 0, 2);

        
       sprintf(hname,"MassCutAfterTgt_%s_%i",target,k);
       sprintf(title,"Mass Cut After Target %s Momentums between %s",target,between);
       MassCutAfterTgt[i][k] = new TH1F(hname,title,400, 0, 2);}
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
