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


using namespace std;

  Float_t Alo = -20.0;
  Float_t Ahi = 220.0;
  Int_t Abins = 10*Int_t(Ahi-Alo);
  Float_t Abwid = (Ahi-Alo)/float(Abins);

  Float_t Rlo = 0.0;
  Float_t Rhi = 10.0;
  Int_t Rbins = 10*Int_t(Rhi-Rlo);
  Float_t Rbwid = (Rhi-Rlo)/float(Abins);

  Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
  Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};

  Float_t Lmar = 0.125;
  Float_t Rmar = 0.125;
  Float_t yoff = 1.5;

  Float_t Plim = 0.5;
  Float_t dtlim = 1.002;
  Float_t dzlim = 3.0;
  Float_t rlim = 2.0;


  Int_t nIM = 240;
  Float_t IMl = 0.0;
  Float_t IMh = 1.2;

  Int_t nP = 300;
  Float_t Pl = 0.0;
  Float_t Ph = 3.0;

  char *Plabel[3] = {"#pi^{+}","#pi^{-}","#pi^{+}#pi^{-}"};

  const Int_t nTgt = 6;
  const Int_t nFoils = 8;
  const Int_t tgtTot = nTgt + nFoils;
  char *target[nTgt] = {"All","2H","C","Fe-Ti","Pb","No Pb"};
  char *Foil[nFoils] = {"2H","C1","Fe","C2","Pb","C3","Ti","C4"};

  Float_t zcut[nFoils+1] ={-22.0,-14.0,-11.0,-8.5,-6.0,-3.5,-1.0,1.5,4.0};

  Double_t zFoil_mean[8] = {-18.0,-12.3,-9.89,-7.42,-4.93,-2.46,0.02,2.44};      // foil peak mean
  Double_t zFoil_sigma[8] = {1.17,0.258,0.249,0.252,0.235,0.234,0.245,0.231};    // foil peak sigma

  Int_t nVz = 270;
  Float_t zlo = -22.0;
  Float_t zhi = 5.0;


  char hname[50];
  char htitle[500];
  char title[500];
  char cname[50];
  char ctitle[500];
  char xtitle[100];
  char ytitle[100];
  char OutCan[100];
  char OutText[100];

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
/*////*/
  TH1F *MomentumCutAround1;
  TH1F *MomentumCutAround2;
  TH1F *MomentumCutAround3;
  TH1F *MomentumCutAround4;
  TH1F *MomentumCutAround5;
  TH1F *MomentumCutAround6;
////
  TH1F *MomentumCutAfter1;
  TH1F *MomentumCutAfter2;
  TH1F *MomentumCutAfter3;
  TH1F *MomentumCutAfter4;
  TH1F *MomentumCutAfter5;
  TH1F *MomentumCutAfter6;
////
  TH1F *MassCutAround1;
  TH1F *MassCutAround2;
  TH1F *MassCutAround3;
  TH1F *MassCutAround4;
  TH1F *MassCutAround5;
  TH1F *MassCutAround6;
////
  TH1F *MassCutAfter1;
  TH1F *MassCutAfter2;
  TH1F *MassCutAfter3;
  TH1F *MassCutAfter4;
  TH1F *MassCutAfter5;
  TH1F *MassCutAfter6;
/*////*/


  TH1F *hYldVsA;
  TH1F *hYldVsR;

  struct labels {
    char *x;
    char *y;
    char *hPrefix;
    char *fPrefix;};


  //declarations of functions

  Int_t FindFoilIndex(Float_t vz);
  Int_t FindVzBehindFoilIndex(Float_t vz);
  Int_t FindTargetIndex(Float_t vz, Int_t method);
  Int_t FindNoPb(Float_t vz, Int_t method);
  void WriteHist(char *RootFile);
  void BookHists_Yields();
  void WriteHists_Yields(char *RootFile);
  struct labels SelectHistogram(Int_t flag);
  void BookHist();
  void Projection(char *infilename,char *outfilename);




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
  Float_t Plim;

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
      	   ////if(!(i%dEvents)) cerr<<i<<" / "<<nEvents<<" / "<<nfiles<<"\r"; 

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
         
         
/*////*/
		if (tgt[j] != 0){
		   if(j==0){
		      if(pairV4.P()>0 and pairV4.P()<0.5){ 
  			    MomentumCutAround1->Fill(pairV4.P());
  			    MassCutAround1->Fill(pairM);}
  		      if(pairV4.P()>0.5 and pairV4.P()<1){
		         MomentumCutAround2->Fill(pairV4.P());
		         MassCutAround2->Fill(pairM);}
		      if(pairV4.P()>1 and pairV4.P()<1.5){
		         MomentumCutAround3->Fill(pairV4.P());
		         MassCutAround3->Fill(pairM);}
		      if(pairV4.P()>1.5 and pairV4.P()<2){
		         MomentumCutAround4->Fill(pairV4.P());
		         MassCutAround4->Fill(pairM);}
		      if(pairV4.P()>2 and pairV4.P()<2.5){
		         MomentumCutAround5->Fill(pairV4.P());
		         MassCutAround5->Fill(pairM);}
		      if(pairV4.P()>2.5 and pairV4.P()<3){
		         MomentumCutAround6->Fill(pairV4.P());
		         MassCutAround6->Fill(pairM);}
		   }
/*////*/
/*////*/
		   if(j==1){
		      if(pairV4.P()>0 and pairV4.P()<0.5){ 
		         MomentumCutAfter1->Fill(pairV4.P());
		         MassCutAfter1->Fill(pairM);}
  		      if(pairV4.P()>0.5 and pairV4.P()<1){
		         MomentumCutAfter2->Fill(pairV4.P());
		         MassCutAfter2->Fill(pairM);}
		      if(pairV4.P()>1 and pairV4.P()<1.5){
		         MomentumCutAfter3->Fill(pairV4.P());
		         MassCutAfter3->Fill(pairM);}
		      if(pairV4.P()>1.5 and pairV4.P()<2){
		         MomentumCutAfter4->Fill(pairV4.P());
		         MassCutAfter4->Fill(pairM);}
		      if(pairV4.P()>2 and pairV4.P()<2.5){
		         MomentumCutAfter5->Fill(pairV4.P());
		         MassCutAfter5->Fill(pairM);}
		      if(pairV4.P()>2.5 and pairV4.P()<3){
		         MomentumCutAfter6->Fill(pairV4.P());
		         MassCutAfter6->Fill(pairM);}
		      }
		   }
/*////*/
      
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
    	////cerr<<endl;
	myTree->Delete(); 						// delete Tree object
	myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
    	nfiles++; 							// increment file counter
  }
  fclose(in1); 							// close file with input file list
  cout<<TotEvents<<" events in "<<nfiles<< " files."<<endl; // print out stats

  WriteHist(RootFile); 						// write histograms to a file
  
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


	Alo = 0.5*Abwid;
	Ahi = 0.5*Abwid;

	Rlo = 0.5*Rbwid;
	Rhi = 0.5*Rbwid;
	
	AnaTwoPion(infilename,outfilename);
	Projection(outfilename,outfilename);
     
	     

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

/*////*/
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
		sprintf(strname,"PionPairMassAround");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi Pair Mass Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();
			
		h2D = (TH2F*)infile->Get("hIM_P1_2");
		sprintf(strname,"PionPairMassAfter");
		h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Mass (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi pair Mass After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();

		h2D = (TH2F*)infile->Get("hIM_P0_2");
		sprintf(strname,"PionPairMomentumAround");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi Pair Momentum Around target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();
  					
     	h2D = (TH2F*)infile->Get("hIM_P1_2");
		sprintf(strname,"PionPairMomentumAfter");
		h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,"");	
  		h1D->SetXTitle("Momentum (GeV)");
  		h1D->GetXaxis()->CenterTitle();
  		h1D->SetYTitle("Counts");
  		h1D->GetYaxis()->CenterTitle();
  		h1D->GetYaxis()->SetTitleOffset(1.2);
  		h1D->SetLineWidth(2);
		h1D->SetTitle("#pi Pair Mometum After target");
  		h1D->Draw();
  		myFile->cd();
  		h1D->Write();


	
	myFile->Write(outfilename);
	myFile->Close(outfilename);



}

/*////*/





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

  Int_t i, j;
  
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
  

/*////*/
  MomentumCutAround1->SetXTitle("Momentum (GeV)");
  MomentumCutAround1->GetXaxis()->CenterTitle();
  MomentumCutAround1->SetYTitle("Counts");
  MomentumCutAround1->GetYaxis()->CenterTitle();
  MomentumCutAround1->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAround1->Write();
  
  MomentumCutAround2->SetXTitle("Momentum (GeV)");
  MomentumCutAround2->GetXaxis()->CenterTitle();
  MomentumCutAround2->SetYTitle("Counts");
  MomentumCutAround2->GetYaxis()->CenterTitle();
  MomentumCutAround2->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAround2->Write();
  
  MomentumCutAround3->SetXTitle("Momentum (GeV)");
  MomentumCutAround3->GetXaxis()->CenterTitle();
  MomentumCutAround3->SetYTitle("Counts");
  MomentumCutAround3->GetYaxis()->CenterTitle();
  MomentumCutAround3->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAround3->Write();
  
  MomentumCutAround4->SetXTitle("Momentum (GeV)");
  MomentumCutAround4->GetXaxis()->CenterTitle();
  MomentumCutAround4->SetYTitle("Counts");
  MomentumCutAround4->GetYaxis()->CenterTitle();
  MomentumCutAround4->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAround4->Write();
  
  MomentumCutAround5->SetXTitle("Momentum (GeV)");
  MomentumCutAround5->GetXaxis()->CenterTitle();
  MomentumCutAround5->SetYTitle("Counts");
  MomentumCutAround5->GetYaxis()->CenterTitle();
  MomentumCutAround5->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAround5->Write();
  
  MomentumCutAround6->SetXTitle("Momentum (GeV)");
  MomentumCutAround6->GetXaxis()->CenterTitle();
  MomentumCutAround6->SetYTitle("Counts");
  MomentumCutAround6->GetYaxis()->CenterTitle();
  MomentumCutAround6->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAround6->Write();
  
////
  MomentumCutAfter1->SetXTitle("Momentum (GeV)");
  MomentumCutAfter1->GetXaxis()->CenterTitle();
  MomentumCutAfter1->SetYTitle("Counts");
  MomentumCutAfter1->GetYaxis()->CenterTitle();
  MomentumCutAfter1->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAfter1->Write();
  
  MomentumCutAfter2->SetXTitle("Momentum (GeV)");
  MomentumCutAfter2->GetXaxis()->CenterTitle();
  MomentumCutAfter2->SetYTitle("Counts");
  MomentumCutAfter2->GetYaxis()->CenterTitle();
  MomentumCutAfter2->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAfter2->Write();
  
  MomentumCutAfter3->SetXTitle("Momentum (GeV)");
  MomentumCutAfter3->GetXaxis()->CenterTitle();
  MomentumCutAfter3->SetYTitle("Counts");
  MomentumCutAfter3->GetYaxis()->CenterTitle();
  MomentumCutAfter3->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAfter3->Write();
  
  MomentumCutAfter4->SetXTitle("Momentum (GeV)");
  MomentumCutAfter4->GetXaxis()->CenterTitle();
  MomentumCutAfter4->SetYTitle("Counts");
  MomentumCutAfter4->GetYaxis()->CenterTitle();
  MomentumCutAfter4->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAfter4->Write();
  
  MomentumCutAfter5->SetXTitle("Momentum (GeV)");
  MomentumCutAfter5->GetXaxis()->CenterTitle();
  MomentumCutAfter5->SetYTitle("Counts");
  MomentumCutAfter5->GetYaxis()->CenterTitle();
  MomentumCutAfter5->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAfter5->Write();

  MomentumCutAfter6->SetXTitle("Momentum (GeV)");
  MomentumCutAfter6->GetXaxis()->CenterTitle();
  MomentumCutAfter6->SetYTitle("Counts");
  MomentumCutAfter6->GetYaxis()->CenterTitle();
  MomentumCutAfter6->GetYaxis()->SetTitleOffset(1.2);
  MomentumCutAfter6->Write();
////
  MassCutAround1->SetXTitle("Mass (GeV)");
  MassCutAround1->GetXaxis()->CenterTitle();
  MassCutAround1->SetYTitle("Counts");
  MassCutAround1->GetYaxis()->CenterTitle();
  MassCutAround1->GetYaxis()->SetTitleOffset(1.2);
  MassCutAround1->Write();

  MassCutAround2->SetXTitle("Mass (GeV)");
  MassCutAround2->GetXaxis()->CenterTitle();
  MassCutAround2->SetYTitle("Counts");
  MassCutAround2->GetYaxis()->CenterTitle();
  MassCutAround2->GetYaxis()->SetTitleOffset(1.2);
  MassCutAround2->Write();
  
  MassCutAround3->SetXTitle("Mass (GeV)");
  MassCutAround3->GetXaxis()->CenterTitle();
  MassCutAround3->SetYTitle("Counts");
  MassCutAround3->GetYaxis()->CenterTitle();
  MassCutAround3->GetYaxis()->SetTitleOffset(1.2);
  MassCutAround3->Write();
  
  MassCutAround4->SetXTitle("Mass (GeV)");
  MassCutAround4->GetXaxis()->CenterTitle();
  MassCutAround4->SetYTitle("Counts");
  MassCutAround4->GetYaxis()->CenterTitle();
  MassCutAround4->GetYaxis()->SetTitleOffset(1.2);
  MassCutAround4->Write();
  
  MassCutAround5->SetXTitle("Mass (GeV)");
  MassCutAround5->GetXaxis()->CenterTitle();
  MassCutAround5->SetYTitle("Counts");
  MassCutAround5->GetYaxis()->CenterTitle();
  MassCutAround5->GetYaxis()->SetTitleOffset(1.2);
  MassCutAround5->Write();
  
  MassCutAround6->SetXTitle("Mass (GeV)");
  MassCutAround6->GetXaxis()->CenterTitle();
  MassCutAround6->SetYTitle("Counts");
  MassCutAround6->GetYaxis()->CenterTitle();
  MassCutAround6->GetYaxis()->SetTitleOffset(1.2);
  MassCutAround6->Write();


////
  MassCutAfter1->SetXTitle("Mass (GeV)");
  MassCutAfter1->GetXaxis()->CenterTitle();
  MassCutAfter1->SetYTitle("Counts");
  MassCutAfter1->GetYaxis()->CenterTitle();
  MassCutAfter1->GetYaxis()->SetTitleOffset(1.2);
  MassCutAfter1->Write();
  
  MassCutAfter2->SetXTitle("Mass (GeV)");
  MassCutAfter2->GetXaxis()->CenterTitle();
  MassCutAfter2->SetYTitle("Counts");
  MassCutAfter2->GetYaxis()->CenterTitle();
  MassCutAfter2->GetYaxis()->SetTitleOffset(1.2);
  MassCutAfter2->Write();
  
  MassCutAfter3->SetXTitle("Mass (GeV)");
  MassCutAfter3->GetXaxis()->CenterTitle();
  MassCutAfter3->SetYTitle("Counts");
  MassCutAfter3->GetYaxis()->CenterTitle();
  MassCutAfter3->GetYaxis()->SetTitleOffset(1.2);
  MassCutAfter3->Write();
  
  MassCutAfter4->SetXTitle("Mass (GeV)");
  MassCutAfter4->GetXaxis()->CenterTitle();
  MassCutAfter4->SetYTitle("Counts");
  MassCutAfter4->GetYaxis()->CenterTitle();
  MassCutAfter4->GetYaxis()->SetTitleOffset(1.2);
  MassCutAfter4->Write();
  
  MassCutAfter5->SetXTitle("Mass (GeV)");
  MassCutAfter5->GetXaxis()->CenterTitle();
  MassCutAfter5->SetYTitle("Counts");
  MassCutAfter5->GetYaxis()->CenterTitle();
  MassCutAfter5->GetYaxis()->SetTitleOffset(1.2);
  MassCutAfter5->Write();
  
  MassCutAfter6->SetXTitle("Mass (GeV)");
  MassCutAfter6->GetXaxis()->CenterTitle();
  MassCutAfter6->SetYTitle("Counts");
  MassCutAfter6->GetYaxis()->CenterTitle();
  MassCutAfter6->GetYaxis()->SetTitleOffset(1.2);
  MassCutAfter6->Write();
/*////*/

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

struct labels SelectHistogram(Int_t flag){
  struct labels temp;

  switch(flag){
  case 0:
    temp.x = "A";
    temp.y = "Y/N_{t}#epsilon";
    temp.hPrefix = "hXsnVsA";
    temp.fPrefix = "Ana_Xsn_A";
    break;
  case 1:
    temp.x = "Radius (fm)";
    temp.y = "Y/N_{t}#epsilon";
    temp.hPrefix = "hXsnVsR";
    temp.fPrefix = "Ana_Xsn_R";
    break;
  case 2:
    temp.x = "A";
    temp.y = "T";
    temp.hPrefix = "hTVsA";
    temp.fPrefix = "Ana_T_A";
    break;
  case 3:
    temp.x = "Radius (fm)";
    temp.y = "T";
    temp.hPrefix = "hTVsR";
    temp.fPrefix = "Ana_T_R";
    break;
  case 4:
    temp.x = "A";
    temp.y = "T_{A}/T_{^{2}H}";
    temp.hPrefix = "hRatVsA";
    temp.fPrefix = "Ana_Rat_A";
    break;
  case 5:
    temp.x = "Radius (fm)";
    temp.y = "T_{A}/T_{^{2}H}";
    temp.hPrefix = "hRatVsR";
    temp.fPrefix = "Ana_Rat_R";
    break;
  case 6:
    temp.x = "A";
    temp.y = "Yield";
    temp.hPrefix = "hYldVsA";
    temp.fPrefix = "Ana_Yld_A";
    break;
  case 7:
    temp.x = "Radius (fm)";
    temp.y = "Yield";
    temp.hPrefix = "hYldVsR";
    temp.fPrefix = "Ana_Yld_R";
    break;
  case 8:
    temp.x = "A";
    temp.y = "T";
    temp.hPrefix = "hTVsA";
    temp.fPrefix = "Ana_T_A_C";
    break;
  case 9:
    temp.x = "Radius (fm)";
    temp.y = "T";
    temp.hPrefix = "hTVsR";
    temp.fPrefix = "Ana_T_R_C";
    break;
  case 10:
    temp.x = "A";
    temp.y = "T_{A}/T_{^{12}C}";
    temp.hPrefix = "hRatVsA_C";
    temp.fPrefix = "Ana_Rat_A_C";
    break;
  case 11:
    temp.x = "Radius (fm)";
    temp.y = "T_{A}/T_{^{12}C}";
    temp.hPrefix = "hRatVsR_C";
    temp.fPrefix = "Ana_Rat_R_C";
    break;
  case 12:
    temp.x = "A";
    temp.y = "#epsilon";
    temp.hPrefix = "hAccVsA";
    temp.fPrefix = "Ana_Acc_A";
    break;
  case 13:
    temp.x = "A";
    temp.y = "#epsilon";
    temp.hPrefix = "hAccVsA_fin";
    temp.fPrefix = "Ana_Acc_A_fin";
    break;
  case 14:
    temp.x = "A";
    temp.y = "N_{t}";
    temp.hPrefix = "hNtVsA";
    temp.fPrefix = "Ana_Nt_A";
    break;
  case 15:
    temp.x = "A";
    temp.y = "N_{t}";
    temp.hPrefix = "hNtVsA_fin";
    temp.fPrefix = "Ana_Nt_A_fin";
    break;
  case 16:
    temp.x = "A";
    temp.y = "d(g/cm^{2})";
    temp.hPrefix = "hThkVsA";
    temp.fPrefix = "Ana_Thk_A";
    break;
  case 17:
    temp.x = "A";
    temp.y = "d(g/cm^{2})";
    temp.hPrefix = "hThkVsA_fin";
    temp.fPrefix = "Ana_Thk_A_fin";
    break;
  case 18:
    temp.x = "A";
    temp.y = "N_{t}#epsilon";
    temp.hPrefix = "hNormVsA";
    temp.fPrefix = "Ana_Norm_A";
    break;
  default:
    cout<<"Invalid flag "<<flag<<endl;
    exit(1);
    break;
  }
  return temp;
}


// 
// BookHist - routine to set up histograms for the AnaTwoPion() routine
//
void BookHist(){
  Int_t i, j;

  sprintf(hname,"hEbeam");
  sprintf(title,"Photon Energy (GeV)");
  hEbeam = new TH2F(hname,title,100,0.5,4.0,2,-0.5,1.5);

  sprintf(hname,"ht");
  sprintf(title,"Momentum Transfer");
  ht = new TH2F(hname,title,100,-2.5,0.0,2,-0.5,1.5);

  sprintf(hname,"hIM_cuts");
  sprintf(title,"IM by cuts");
  hIM_cuts = new TH2F(hname,title,nIM,IMl,IMh,6,0.5,6.5);
  
  sprintf(hname,"hVz");
  sprintf(title,"z vertex");
  hVz = new TH2F(hname,title,nVz,zlo,zhi,2,-0.5,1.5);
  
/*////*/
  MomentumCutAround1 = new TH1F("MomentumCutAround1","Momentum Cut Around Target 0-0.5",60, 0, 0.5);
  MomentumCutAround2 = new TH1F("MomentumCutAround2","Momentum Cut Around Target 0.5-1",60 ,0.5, 1);
  MomentumCutAround3 = new TH1F("MomentumCutAround3","Momentum Cut Around Target 1-1.5",60 ,1, 1.5);
  MomentumCutAround4 = new TH1F("MomentumCutAround4","Momentum Cut Around Target 1.5-2",60 ,1.5, 2);
  MomentumCutAround5 = new TH1F("MomentumCutAround5","Momentum Cut Around Target 2-2.5",60 ,2, 2.5);
  MomentumCutAround6 = new TH1F("MomentumCutAround6","Momentum Cut Around Target 2.5-3",60 ,2.5, 3);
////
  MomentumCutAfter1 = new TH1F("MomentumCutAfter1","Momentum Cut After Target 0-0.5",60, 0, 0.5);
  MomentumCutAfter2 = new TH1F("MomentumCutAfter2","Momentum Cut After Target 0.5-1",60, 0.5, 1);
  MomentumCutAfter3 = new TH1F("MomentumCutAfter3","Momentum Cut After Target 1-1.5",60, 1, 1.5);
  MomentumCutAfter4 = new TH1F("MomentumCutAfter4","Momentum Cut After Target 1.5-2",60, 1.5, 2);
  MomentumCutAfter5 = new TH1F("MomentumCutAfter5","Momentum Cut After Target 2-2.5",60, 2, 2.5);
  MomentumCutAfter6 = new TH1F("MomentumCutAfter6","Momentum Cut After Target 2.5-3",60, 2.5, 3);
////
  MassCutAfter1 = new TH1F("MassCutAfter1","Mass Cut After Target Momentums between 0-0.5",400, 0, 2);
  MassCutAfter2 = new TH1F("MassCutAfter2","Mass Cut After Target Momentums between 0.5-1",400 ,0, 2);
  MassCutAfter3 = new TH1F("MassCutAfter3","Mass Cut After Target Momentums between 1-1.5",400 ,0, 2);
  MassCutAfter4 = new TH1F("MassCutAfter4","Mass Cut After Target Momentums between 1.5-2",400 ,0, 2);
  MassCutAfter5 = new TH1F("MassCutAfter5","Mass Cut After Target Momentums between 2-2.5",400 ,0, 2);
  MassCutAfter6 = new TH1F("MassCutAfter6","Mass Cut After Target Momentums between 2.5-3",400 ,0, 2);
////
  MassCutAround1 = new TH1F("MassCutAround1","Mass Cut Around Target Momentums between 0-0.5",400, 0, 2);
  MassCutAround2 = new TH1F("MassCutAround2","Mass Cut Around Target Momentums between 0.5-1",400 ,0, 2);
  MassCutAround3 = new TH1F("MassCutAround3","Mass Cut Around Target Momentums between 1-1.5",400 ,0, 2);
  MassCutAround4 = new TH1F("MassCutAround4","Mass Cut Around Target Momentums between 1.5-2",400 ,0, 2);
  MassCutAround5 = new TH1F("MassCutAround5","Mass Cut Around Target Momentums between 2-2.5",400 ,0, 2);
  MassCutAround6 = new TH1F("MassCutAround6","Mass Cut Around Target Momentums between 2.5-3",400 ,0, 2);
/*////*/
  
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
