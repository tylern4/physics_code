//AnaTwoPion.C
//
// 
// M. H. Wood, Canisius College
//
// Notes: macro to study pi+pi-
//
//--------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();
 
Float_t Alo = -20.0;
Float_t Ahi = 220.0;
Int_t Abins = 10*(Ahi-Alo);
Float_t Abwid = (Ahi-Alo)/float(Abins);
Alo -= 0.5*Abwid;
Ahi -= 0.5*Abwid;

Float_t Rlo = 0.0;
Float_t Rhi = 10.0;
Int_t Rbins = 10*(Rhi-Rlo);
Float_t Rbwid = (Rhi-Rlo)/float(Abins);
Rlo -= 0.5*Rbwid;
Rhi -= 0.5*Rbwid;

Int_t lcol[10] = {1,2,4,6,7,8,9,10,14,15};
Int_t mkr[10] = {20,21,22,24,23,25,26,27,28,29};
char *fSame[10] = {"","same","same","same","same","same","same","same","same","same"};

Float_t R0 = 1.2;
Float_t AvoNum = 6.022137E23; // Avigadro's number [#/mol]

Float_t Lmar = 0.125;
Float_t Rmar = 0.125;
Float_t yoff = 1.5;

Float_t Plim = 0.5;
Float_t dtlim = 1.002;
Float_t dzlim = 3.0;
Float_t rlim = 2.0;

const Int_t Psteps = 9;
Double_t Pwidth = 0.1;

Int_t nIM = 240;
Float_t IMl = 0.0;
Float_t IMh = 1.2;

Int_t nP = 200;
Float_t Pl = 0.0;
Float_t Ph = 2.0;

char *Plabel[3] = {"#pi^{+}","#pi^{-}","#pi^{+}#pi^{-}"};

const Int_t nTgt = 6;
const Int_t nFoils = 8;
const Int_t tgtTot = nTgt + nFoils;
char *target[nTgt] = {"All","2H","C","Fe-Ti","Pb","No Pb"};
char *Foil[nFoils] = {"2H","C1","Fe","C2","Pb","C3","Ti","C4"};

Float_t zcut[nFoils+1] ={-22.0,-14.0,-11.0,-8.5,-6.0,-3.5,-1.0,1.5,4.0};

Double_t zFoil_mean[8] = {-18.0,-12.3,-9.89,-7.42,-4.93,-2.46,0.02,2.44}; // foil peak mean
Double_t zFoil_sigma[8] = {1.17,0.258,0.249,0.252,0.235,0.234,0.245,0.231}; // foil peak sigma

Int_t nVz = 270;
Float_t zlo = -22.0;
Float_t zhi = 5.0;

Float_t Mlo_Ks = 0.48;
Float_t Mhi_Ks = 0.51;

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

TH1F *hYldVsA;
TH1F *hYldVsR;

struct labels {
  char *x;
  char *y;
  char *hPrefix;
  char *fPrefix;
};

//
//	AnaTwoPion - main analysis routine. Read in filtered ROOT data, apply cuts, and plot invariant mass.
//
//               fin = text file with a list of filtered ROOT files
//				 RootFile = output ROOT file
//               MaxEvents = max. number of events to analyze
//               dEvents = event counter for print statements
//
void AnaTwoPion(char *fin="TwoPion.lis", char *RootFile="AnaTwoPion.root", Int_t MaxEvents=0, Int_t dEvents=10000)
{
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
  
  TVector3 V3[nPart]; // vertex 3-vectors
  TLorentzVector V4[nPart]; // particle 4-vectors
  TLorentzVector pairV4; // particle pair 4-vector
  TLorentzVector Beam; // beam 4-vector
  TLorentzVector A; // beam-pair 4-vector

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
  printf("Analyzing file %s\n",fin);  

  BookHist(); // declare histograms

  TFile *myFile;
  TTree *myTree;
  Int_t ncols;
  Int_t nfiles = 0;
  char rootFile[500];
  FILE *in1 = fopen(fin,"r");
  while (1) {
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
      pairV4 = V4[0] + V4[1]; // particle pair 4-vector
	  pairM = pairV4.M(); // pair invariant mass
	
	  A = Beam - pairV4; // beam minus pair 4-vector
	
      dz = V3[0].Z() - V3[1].Z(); // z vertex difference
      r_pair = sqrt(pairX*pairX + pairY*pairY); // radial vertex (x,y) for the particle pair

      P_cut = (V4[0].P()>=Plim && V4[1].P()>=Plim); // momentum cut
      vz_cut = (pairZ>=zcut[0] && pairZ>=zcut[nFoils]); // vertex z cut
      dz_cut = (abs(dz)<=dzlim); // vertex z difference cut
      dt_cut = (abs(scTIME[1]-scTIME[0])<=dtlim); // vertex time difference cut
      r_cut = (r_pair<=rlim); // target radius cut 
      isPim = (PID[0]==9); // pi- GEANT id
      isPip = (PID[1]==8); // pi+ GEANT id

      cuts=(dz_cut && r_cut && dt_cut && isPim && isPip);  // total cuts

	  hEbeam->Fill(Beam.E(),0); // beam histogram
	  ht->Fill(A.Mag2(),0); // histogram of t (mom. transfer)

      if(isPim && isPip){ // cut histograms
        hIM_cuts->Fill(pairM,1); // no cuts
        if(dz_cut) hIM_cuts->Fill(pairM,2); // dz cut only
        if(r_cut) hIM_cuts->Fill(pairM,3); // r cut only
        if(dt_cut) hIM_cuts->Fill(pairM,4); // dt cut only
        if(dz_cut && r_cut) hIM_cuts->Fill(pairM,5); // dz and r cuts
        if(dz_cut && r_cut && dt_cut) hIM_cuts->Fill(pairM,6); // dz, r, and dt cuts
      }

      for(j=0;j<nType;j++){ // loop over vz type, [0,in target], [1,after target]
        tgt[j] = FindTargetIndex(pairZ,j); // target index (1-4)
        NoPbtgt[j] = FindNoPb(pairZ,j); // no Pb index (5)
        if(NoPbtgt[j]) NoPbtgt[j] = nTgt-1; 
        Ftgt[j] = FindFoilIndex(pairZ); // foil index (6-13)
        if(Ftgt[j]) Ftgt[j] += nTgt-1;
      
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
	  i++; // increment event counter
	  TotEvents++; // increment total event counter
    }
    cerr<<endl;
	myTree->Delete(); // delete Tree object
	myFile->Close("R"); // close input ROOT file.  The R flag deletes TProcessIDs
    nfiles++; // increment file counter
  }
  fclose(in1); // close file with input file list
  cout<<"Total "<<TotEvents<<" / "<<nfiles<<endl; // print out stats
	
  WriteHist(RootFile); // write histograms to a file
}

// 
// BookHist - routine to set up histograms for the AnaTwoPion() routine
//
void BookHist()
{
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

// 
// WriteHist - routine to write histograms to the output file
//             Called in the AnaTwoPion() routine
//
void WriteHist(char *RootFile)
{
  Int_t i, j;

  TFile *fa = new TFile(RootFile,"RECREATE");
  fa->cd();
  hEbeam->Write();
  ht->Write();
  hIM_cuts->Write();
  hVz->Write();
  for(i=0;i<nPart;i++){
    hIM_tgt[i]->Write();
    hPcomp[i]->Write();
    for(j=0;j<3;j++) hIM_P[i][j]->Write();
  }
  fa->Close();
}

// 
// plotHist2D - plot 2D histogram with labels
//                  
//                  fAna = output from AnaTwoPion()
//                  hname = 2D histogram name
//                  xname = x-axis label
//                  yname = y-axis label
//
void plotHist2D(char *fAna="Ana.root", char *hname="hEbeam", char *xname="E_{#gamma}", char *yname="Cut Index")
{
  TH2F *h2D;

  // Canvas to compare IM by each cut
  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
  c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  c1->SetBorderSize(5); 
  gStyle->SetOptStat(0);
  c1->SetFillStyle(4000);

  // data files contain the trees
  printf("Analyzing file %s\n",fAna);  
  TFile *fm = new TFile(fAna,"READ");

  c1->cd();
  gPad->SetLeftMargin(Lmar);
  gPad->SetRightMargin(Rmar);
  gPad->SetFillColor(0);
    
  h2D = (TH2F*)fm->Get(hname);
  h2D->SetTitle(0);
  h2D->SetXTitle(xname);
  h2D->GetXaxis()->CenterTitle();
  h2D->SetYTitle(yname);
  h2D->GetYaxis()->CenterTitle();
  h2D->GetYaxis()->SetTitleOffset(yoff);
  h2D->SetLineWidth(2);
  h2D->Draw("colz");

  sprintf(OutCan,"AnaTwoPion_%s.gif",hname);
  c1->Print(OutCan);
  sprintf(OutCan,"AnaTwoPion_%s.eps",hname);
  c1->Print(OutCan);
}

// 
// projectX_1Dfrom2D - project 1D histogram along the x-axis from a 2D histogram
//                  
//                  fAna = output from AnaTwoPion()
//                  WhichAxis = axis to project upon (x or y)
//                  binLo = set lower bin of the range
//                  binHi = set higher bin of the range
//                  hname = 2D histogram name
//                  xname = x-axis label
//                  yname = y-axis label
//
void project_1Dfrom2D(char *fAna="Ana.root", char *WhichAxis="x", Int_t binLo=-1, Int_t binHi=-1, char *hname="hEbeam", char *xname="E_{#gamma}", char *yname="Counts")
{
  char strname[50];

  TH1D *h1D;
  TH2F *h2D;

  // Canvas to compare IM by each cut
  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
  c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  c1->SetBorderSize(5); 
  gStyle->SetOptStat(0);
  c1->SetFillStyle(4000);

  // data files contain the trees
  printf("Analyzing file %s\n",fAna);  
  TFile *fm = new TFile(fAna,"READ");

  c1->cd();
  gPad->SetLeftMargin(Lmar);
  gPad->SetRightMargin(Rmar);
  gPad->SetFillColor(0);
    
  h2D = (TH2F*)fm->Get(hname);
  sprintf(strname,"h1D_p");
  switch(WhichAxis){
  case "x":
  case "X":
    h1D = (TH1D*)h2D->ProjectionX(strname,binLo,binHi,""); break;  
  case "y":
  case "Y":
    h1D = (TH1D*)h2D->ProjectionY(strname,binLo,binHi,""); break;
  default:
    cout<<"No axis "<<WhichAxis<<endl;
	cout<<"Try x or y"<<endl;
	exit(0);
  }
  h1D->SetTitle(0);
  h1D->SetXTitle(xname);
  h1D->GetXaxis()->CenterTitle();
  h1D->SetYTitle(yname);
  h1D->GetYaxis()->CenterTitle();
  h1D->GetYaxis()->SetTitleOffset(yoff);
  h1D->SetLineWidth(2);
  h1D->Draw();

  sprintf(OutCan,"AnaTwoPion_%s_%s%i_%i.gif",hname,WhichAxis,binLo,binHi);
  c1->Print(OutCan);
  sprintf(OutCan,"AnaTwoPion_%s_%s%i_%i.eps",hname,WhichAxis,binLo,binHi);
  c1->Print(OutCan);
}

// 
// compIMvsTarget - compare the inv. mass by target.  Each histogram is plotted in its own pad
//                  
//                  fAna = output from AnaTwoPion()
//                  iType = target region [0], after target [1]
//
void compIMvsTarget(char *fAna="Ana.root", Int_t iType=0)
{
  Int_t i;
  char strname[50];

  TH1D *h1D[nTgt];
  TH2F *h2D;

  if(iType<0 || iType>=nType){
    cout<<"Unknown analysis Type "<<iType<<endl;
    exit;
  }

  // Canvas to compare IM by each cut
  TCanvas *c1 = new TCanvas("c1","c1",0,0,1200,700);
  c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  c1->SetBorderSize(5); 
  gStyle->SetOptStat(0);
  c1->SetFillStyle(4000);
  c1->Divide(3,2);

  // data files contain the trees
  printf("Analyzing file %s\n",fAna);  
  TFile *fm = new TFile(fAna,"READ");

  sprintf(xtitle,"Invariant Mass (GeV)");
  sprintf(ytitle,"Counts");

  for(i=0;i<nTgt;i++){
    c1->cd(i+1);
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    sprintf(hname,"hIM_tgt%i",iType);
    h2D = (TH2F*)fm->Get(hname);
    sprintf(strname,"h1D_px%i",i);
    h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
    h1D[i]->SetTitle(target[i]);
    h1D[i]->SetXTitle(xtitle);
    h1D[i]->GetXaxis()->CenterTitle();
    h1D[i]->SetYTitle(ytitle);
    h1D[i]->GetYaxis()->CenterTitle();
    h1D[i]->GetYaxis()->SetTitleOffset(yoff);
    h1D[i]->SetLineWidth(2);
    h1D[i]->Draw();
  }

  sprintf(OutCan,"IMvsTgt_comp_PipPim%i.gif",iType);
  c1->Print(OutCan);
  sprintf(OutCan,"IMvsTgt_comp_PipPim%i.eps",iType);
  c1->Print(OutCan);
}

// 
// overlayIMvsTarget - overlay the inv. mass for each target in a single pad
//                  
//                  fAna = output from AnaTwoPion()
//                  iType = target region [0], after target [1]
//
void overlayIMvsTarget(char *fAna="Ana.root", Int_t iType=0)
{
  Int_t i, j, xbin;
  char strname[50];

  TH1D *h1D[nTgt];
  TH2F *h2D;

  // Canvas to compare IM by each cut
  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
  c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  c1->SetBorderSize(5); 
  gStyle->SetOptStat(0);
  c1->SetFillStyle(4000);
 
  // data files contain the trees
  printf("Analyzing file %s\n",fAna);  
  TFile *fm = new TFile(fAna,"READ");

  sprintf(xtitle,"Invariant Mass (GeV)");
  sprintf(ytitle,"Counts");

  TLegend *leg = new TLegend(0.2,0.625,0.65,0.865);
  for(i=0;i<nTgt;i++){
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    sprintf(hname,"hIM_tgt%i",iType);
    h2D = (TH2F*)fm->Get(hname);
    sprintf(strname,"h1D_px%i",i);
    h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
    if(h1D[i]->Integral()) h1D[i]->Scale(h1D[0]->Integral()/h1D[i]->Integral());
    h1D[i]->SetTitle(0);
    h1D[i]->Sumw2();
    h1D[i]->SetXTitle(xtitle);
    h1D[i]->GetXaxis()->CenterTitle();
    h1D[i]->SetYTitle(ytitle);
    h1D[i]->GetYaxis()->CenterTitle();
    h1D[i]->GetYaxis()->SetTitleOffset(yoff);
    h1D[i]->SetLineColor(lcol[i]);
    h1D[i]->SetLineWidth(2);
    h1D[i]->Draw(fSame[i]);
    leg->AddEntry(h1D[i],target[i],"l");
  }
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader("Targets:");
  leg->Draw();

  sprintf(OutCan,"IMvsTgt_overlay_PipPim%i.gif",iType);
  c1->Print(OutCan);
  sprintf(OutCan,"IMvsTgt_overlay_PipPim%i.eps",iType);
  c1->Print(OutCan);
}

// 
// overlayVz - overlay the Vz for each iType in a single pad
//                  
//                  fAna = output from AnaTwoPion()
//
void overlayVz(char *fAna="Ana.root")
{
  Int_t i, j, xbin;
  char strname[50];

  TH1D *h1D[nType];
  TH2F *h2D;

  // Canvas to compare IM by each cut
  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,600);
  c1->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  c1->SetBorderSize(5); 
  gStyle->SetOptStat(0);
  c1->SetFillStyle(4000);
 
  // data files contain the trees
  printf("Analyzing file %s\n",fAna);  
  TFile *fm = new TFile(fAna,"READ");

  sprintf(xtitle,"z (cm)");
  sprintf(ytitle,"Counts");

  TLegend *leg = new TLegend(0.2,0.635,0.5,0.875);
  for(i=0;i<nType;i++){
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    
    sprintf(hname,"hVz");
    h2D = (TH2F*)fm->Get(hname);
    sprintf(strname,"h1D_px%i",i);
    h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
    h1D[i]->SetTitle(0);
    //    h1D[i]->Sumw2();
    h1D[i]->SetXTitle(xtitle);
    h1D[i]->GetXaxis()->CenterTitle();
    h1D[i]->SetYTitle(ytitle);
    h1D[i]->GetYaxis()->CenterTitle();
    h1D[i]->GetYaxis()->SetTitleOffset(yoff);
    h1D[i]->SetLineColor(lcol[i]);
    h1D[i]->SetFillColor(lcol[i]);
    h1D[i]->SetFillStyle(3001);    
    h1D[i]->SetLineWidth(2);
    h1D[i]->Draw(fSame[i]);
    leg->AddEntry(h1D[i],VzAna[i],"l");
  }
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader("Vertex Cuts:");
  leg->Draw();

  sprintf(OutCan,"Vz_overlay_PipPim.gif");
  c1->Print(OutCan);
  sprintf(OutCan,"Vz_overlay_PipPim.eps");
  c1->Print(OutCan);
}

// 
// BgSubIM - background subtraction routine (under development)
//                  
//                  fAna = output from AnaTwoPion()
//                  suffix = output file suffix
//                  iType = target region [0], after target [1]
//
void BgSubIM(char *fAna="Ana.root", char *suffix="100bins", Int_t iType=0)
{
  Int_t i, j, k;
  Float_t cts, err;

  char strname[50];
  char tgtname[50];

  TH1D *h1D[tgtTot];
  TH2F *h2D;
  TCanvas *can[tgtTot];

  // data files contain the trees
  printf("Analyzing file %s\n",fAna);  
  TFile *fd = new TFile(fAna,"READ");
  sprintf(hname,"hIM_tgt%i",iType);
  h2D = (TH2F*)fd->Get(hname);

  for(i=0;i<1;i++){
    //  for(i=0;i<tgtTot;i++){
    if(i<=nTgt){
      strcpy(tgtname,target[i]);
    }else{
      strcpy(tgtname,Foil[i-nTgt-1]);
    }

    sprintf(cname,"can%i",i);
    sprintf(ctitle,"Canvas %i",i);
    can[i] = new TCanvas(cname,ctitle,0,i*25,900,500);
    can[i]->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
    can[i]->SetBorderSize(5); 
    can[i]->Divide(2,1);
    gStyle->SetOptStat(0);

    can[i]->cd(1);
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);

    sprintf(strname,"h1D_px%i",i);
    h1D[i] = (TH1D*)h2D->ProjectionX(strname,i+1,i+1,"");
    h1D[i]->Sumw2();
    h1D[i]->SetLineColor(2);
    h1D[i]->SetLineWidth(2);
    h1D[i]->Draw();

    can[i]->cd(2);
    sprintf(hname,"hIMsub%i",i);
    hIMsub[i] = (TH1F*)h1D[i]->Clone();
    //    hIMsub[i]->Add(h1D[i],-1);
    hIMsub[i]->Draw();

    sprintf(OutCan,"BgSubIM_PipPim_%i.gif",i);
    can[i]->Print(OutCan);
    sprintf(OutCan,"BgSubIM_PipPim_%i.eps",i);
    can[i]->Print(OutCan);

    sprintf(OutText,"PipPim_%s-%s.dat",tgtname,suffix);
    ofstream fout(OutText);
    for(k=1;k<=hIMsub[i]->GetNbinsX();k++){
      cts = (hIMsub[i]->GetBinContent(k)>0) ? hIMsub[i]->GetBinContent(k) : 0;
      err = (hIMsub[i]->GetBinContent(k)>0) ? hIMsub[i]->GetBinError(k) : 1;
      fout<<hIMsub[i]->GetBinCenter(k)<<"\t"<<cts<<"\t"<<err<<endl;
    }
    fout.close();
  }
}

// 
// KsMeson_justBgd - background subtraction routine. The bins defining the Ks meson are removed
//               and the neighboring bins are fit to polynomial functions for a background
//               shape.  The background shape is subtracted from the original spectrum.
//               The output is a list of yields.
//                  
//                  fAna = output from AnaTwoPion()
//                  iTgt = target index
//                  suffix = output file suffix
//                  iType = target region [0], after target [1]
//
void KsMeson_justBgd(char *fAna="Ana.root", Int_t iTgt=0, char *suffix="std", Int_t iType=0)
{
  Float_t fitLo = 0.35;
  Float_t fitHi = 0.6;
  Float_t PeakLo = 0.48;
  Float_t PeakHi = 0.51;

  Int_t i, j, k;
  Int_t tgtA;
  Float_t x, pval;
  const Int_t nFits = 4;
  
  char tgtname[50];
  char strName[50];
  char gName[50];
  char pName[50];
  char tName[50];
  char pFunc[50];
  char tFunc[50];
  Float_t yield[nFits];

  TH1D *h1D;
  TH2F *h2D;

  if(iTgt<=nTgt){
    strcpy(tgtname,target[iTgt]);
  }else{
    strcpy(tgtname,Foil[iTgt-nTgt-1]);
  }

  // data files contain the trees
  printf("Analyzing file %s\n",fAna);  
  TFile *fd = new TFile(fAna,"READ");
  sprintf(hname,"hIM_tgt%i",iType);
  h2D = (TH2F*)fd->Get(hname);

  sprintf(strName,"h1D_px");
  h1D = (TH1D*)h2D->ProjectionX(strName,iTgt+1,iTgt+1,"");

  TH1F *hBgFit[nFits];
  TH1F *MgrNoPeak[nFits];
  for(j=0;j<nFits;j++){
    sprintf(hname,"hBgFit%i",j);
    hBgFit[j] = (TH1F*)h1D->Clone(hname);    

    sprintf(hname,"MgrNoPeak%i",j);
    MgrNoPeak[j] = (TH1F*)h1D->Clone(hname);    
    for(k=1;k<=MgrNoPeak[j]->GetNbinsX();k++){
      x = MgrNoPeak[j]->GetBinCenter(k+1);
      if(x>=PeakLo && x<PeakHi){
	MgrNoPeak[j]->SetBinContent(k,0.0);
	MgrNoPeak[j]->SetBinError(k,0.0);
      }
    }

    yield[i] = 0.0;
  }

  sprintf(hname,"MgrZoom");
  sprintf(title,"g7 K_{s} meson, %s",tgtname);
  MgrZoom = (TH1F*)h1D->Clone(hname);    
  MgrZoom->SetTitle(title);

  // open canvas
  // create canvas
  char xtitle[100];
  char ytitle[100];
  sprintf(title,"g7 K_{s} meson Bgd Comp, %s",tgtname);
  TCanvas *can1 = new TCanvas("can1",title,0,0,800,800);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
  can1->SetBorderSize(5); 
  can1->SetFillStyle(4000); 
  can1->Divide(2,2);

  sprintf(xtitle,"%s Invariant Mass (GeV)",Plabel[2]);
  sprintf(ytitle,"Counts");

  // fit the final spectrum
  Double_t par[8]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  TF1 *pol[nFits];
  for(i=0;i<nFits;i++){
    can1->cd(i+1);
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);

    sprintf(pName,"pol%i",i+1);
    sprintf(pFunc,"pol%i",i+1);
    pol[i] = new TF1(pName,pFunc,fitLo,fitHi);
    pol[i]->SetLineColor(lcol[i+1]);

    MgrNoPeak[i]->Fit(pName,"R");
    //    pol[i]->GetParameters(par);    
    //    pol[i]->SetParameters(par);
    hBgFit[i]->Add(pol[i],-1.0);
    for(k=1; k<hBgFit[i]->GetNbinsX(); k++){
      x = hBgFit[i]->GetBinCenter(k);
      if(x>=fitLo && x<=fitHi){
	pval = hBgFit[i]->GetBinContent(k);
	if(pval<0.0) pval = 0.0;
      }else{
	pval =0.0;
      }
      hBgFit[i]->SetBinContent(k,pval);
      if(pval==0.0) hBgFit[i]->SetBinError(k,0.0);
    }
    yield[i]=hBgFit[i]->Integral(hBgFit[i]->FindBin(PeakLo),hBgFit[i]->FindBin(PeakHi));
  }

  char OutCan[100];
  sprintf(OutCan,"fitKsMeson_justBgd_full%i_%s-%s.gif",iType,tgtname,suffix);
  can1->Print(OutCan);
  sprintf(OutCan,"fitKsMeson_justBgd_full%i_%s-%s.eps",iType,tgtname,suffix);
  can1->Print(OutCan);

  sprintf(title,"g7 K_{s} meson fit, %s",tgtname);
  TCanvas *can2 = new TCanvas("can2",title,50,50,1000,600);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  can2->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
  can2->SetBorderSize(5); 
  can2->SetFillStyle(4000); 
  can2->Divide(2,1);

  TLegend *leg = new TLegend(0.45,0.15,0.85,0.45);
  can2->cd(1);
  gPad->SetLeftMargin(Lmar);
  gPad->SetRightMargin(Rmar);
  gPad->SetFillColor(0);
  MgrZoom->SetTitle(0);
  MgrZoom->GetXaxis()->SetTitle(xtitle);
  MgrZoom->GetXaxis()->CenterTitle();
  MgrZoom->GetYaxis()->SetTitle(ytitle);
  MgrZoom->GetYaxis()->CenterTitle();
  MgrZoom->GetYaxis()->SetTitleOffset(yoff);
  MgrZoom->SetMarkerSize(0.75);
  MgrZoom->SetMarkerStyle(20);
  MgrZoom->GetXaxis()->SetRangeUser(fitLo,fitHi);
  MgrZoom->Draw(); 

  for(i=0;i<nFits;i++){
    pol[i]->SetLineColor(lcol[i+1]);
    pol[i]->Draw("same");
    sprintf(title,"Order %i polynomial",i+1);
    leg->AddEntry(pol[i],title,"l");
  }
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader("Background:");
  leg->Draw();

  can2->cd(2);
  gPad->SetLeftMargin(Lmar);
  gPad->SetRightMargin(Rmar);
  gPad->SetFillColor(0);
  MgrZoom->SetLineWidth(2);
  MgrZoom->SetMinimum(0);
  MgrZoom->Draw(); 
  for(i=0;i<nFits;i++){
    hBgFit[i]->SetTitle(0);
    hBgFit[i]->GetXaxis()->SetTitle(xtitle);
    hBgFit[i]->GetXaxis()->CenterTitle();
    hBgFit[i]->GetYaxis()->SetTitle(ytitle);
    hBgFit[i]->GetYaxis()->CenterTitle();
    hBgFit[i]->GetYaxis()->SetTitleOffset(yoff);
    hBgFit[i]->SetMarkerSize(0.75);
    hBgFit[i]->SetMarkerStyle(21+i);
    hBgFit[i]->SetMarkerColor(lcol[i+1]);
    hBgFit[i]->SetLineColor(lcol[i+1]);
    hBgFit[i]->Draw(fSame[i+1]); 
  }

  TLine *lLow = new TLine(PeakLo,0.0,PeakLo,hBgFit[0]->GetMaximum());
  lLow->SetLineWidth(2);
  lLow->Draw();
  TLine *lHigh = new TLine(PeakHi,0.0,PeakHi,hBgFit[0]->GetMaximum());
  lHigh->SetLineWidth(2);
  lHigh->Draw();

  sprintf(OutCan,"fitKsMeson_justBgd%i_%s-%s.gif",iType,tgtname,suffix);
  can2->Print(OutCan);
  sprintf(OutCan,"fitKsMeson_justBgd%i_%s-%s.eps",iType,tgtname,suffix);
  can2->Print(OutCan);

  char OutFile[100];
  sprintf(OutFile,"fitKsMeson_justBgd%i_%s-%s.lis",iType,tgtname,suffix);
  ofstream fout(OutFile);
  fout<<TargetIndex2A(iTgt);
  for(i=0;i<nFits;i++) fout<<"\t"<<yield[i];
  fout<<endl;
  fout.close();

}

// 
// KsMeson_sbBgd - background subtraction routine. The sidebands around the Ks meson form the
//               background shape and are subtracted the yield under the peak.
//               The output is a list of yields.
//                  
//                  fAna = output from AnaTwoPion()
//                  iTgt = target index
//                  suffix = output file suffix
//                  iType = target region [0], after target [1]
//
void KsMeson_sbBgd(char *fAna="Ana.root", char *suffix="std", Int_t iType=0)
{
  Float_t fitLo = 0.35;
  Float_t fitHi = 0.6;
  Float_t PeakLo = 0.48;
  Float_t PeakHi = 0.51;
  Float_t Ks_width = PeakHi - PeakLo;
  Float_t SidebandLo = PeakLo - 0.5*Ks_width;
  Float_t SidebandHi = PeakHi + 0.5*Ks_width;
  
  Int_t i, j, k;
  Int_t iTgt, tgtA;
  Float_t x, pval;
  Float_t yield, yield_peak, yield_sb;
  
  char tgtname[50];
  char strName[50];

  TH1D *h1D[nTgt-2];
  TH1F *hPeak[nTgt-2];
  TH1F *hSideband[nTgt-2];
  TH2F *h2D;

  char OutFile[100];
  sprintf(OutFile,"fitKsMeson_sbBgd%i_%s.lis",iType,suffix);
  ofstream fout(OutFile);

  // create canvas
  char xtitle[100];
  char ytitle[100];
  sprintf(title,"g7 K_{s} meson, Sideband Analysis");
  TCanvas *can1 = new TCanvas("can1",title,0,0,800,800);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  can1->SetBorderMode(0);  //Bordermode (-1=down, 0 = no border, 1=up)
  can1->SetBorderSize(5); 
  can1->SetFillStyle(4000); 
  can1->Divide(2,2);
  
  sprintf(xtitle,"%s Invariant Mass (GeV)",Plabel[2]);
  sprintf(ytitle,"Counts");

  // data files contain the trees
  printf("Analyzing file %s\n",fAna);  
  TFile *fd = new TFile(fAna,"READ");
  sprintf(hname,"hIM_tgt%i",iType);
  h2D = (TH2F*)fd->Get(hname);

  for(j=0;j<nTgt-2;j++){
    iTgt = j+1;

    if(iTgt<=nTgt){
      strcpy(tgtname,target[iTgt]);
    }else{
      strcpy(tgtname,Foil[iTgt-nTgt-1]);
    }

    sprintf(strName,"h1D_px%i",j);
    h1D[j] = (TH1D*)h2D->ProjectionX(strName,iTgt+1,iTgt+1,"");

    sprintf(hname,"hSideband%i",j);
    hSideband[j] = (TH1F*)h1D[j]->Clone(hname);    

    sprintf(hname,"hPeak%i",j);
    hPeak[j] = (TH1F*)h1D[j]->Clone(hname);    

    for(k=1;k<=h1D[j]->GetNbinsX();k++){
	  x = h1D[j]->GetBinCenter(k);
	  if(x>=PeakLo && x<PeakHi){
	    hPeak[j]->SetBinContent(k,h1D[j]->GetBinContent(k));
      }else{
	    hPeak[j]->SetBinContent(k,0.0);
	  }
	  if(x>=SidebandLo && x<PeakLo){
	    hSideband[j]->SetBinContent(k,h1D[j]->GetBinContent(k));
	  }else if(x>=PeakHi && x<SidebandHi){
	    hSideband[j]->SetBinContent(k,h1D[j]->GetBinContent(k));
      }else{
	    hSideband[j]->SetBinContent(k,0.0);
	  }
    }

    yield_peak = hPeak[j]->Integral(hPeak[j]->FindBin(PeakLo),hPeak[j]->FindBin(PeakHi));
    yield_sb = hSideband[j]->Integral(hSideband[j]->FindBin(SidebandLo),hSideband[j]->FindBin(SidebandHi));
    yield = yield_peak - yield_sb;
    fout<<TargetIndex2A(iTgt)<<"\t"<<yield<<"\t"<<sqrt(yield)<<"\t"<<sqrt(yield)<<endl;
 
    can1->cd(j+1);
    gPad->SetLeftMargin(Lmar);
    gPad->SetRightMargin(Rmar);
    gPad->SetFillColor(0);
    h1D[j]->SetTitle(tgtname);
    h1D[j]->GetXaxis()->SetTitle(xtitle);
    h1D[j]->GetXaxis()->CenterTitle();
    h1D[j]->GetYaxis()->SetTitle(ytitle);
    h1D[j]->GetYaxis()->CenterTitle();
    h1D[j]->GetYaxis()->SetTitleOffset(yoff);
    h1D[j]->GetXaxis()->SetRangeUser(fitLo,fitHi);
    h1D[j]->SetLineWidth(2);
    h1D[j]->Draw(); 
    hPeak[j]->SetLineWidth(2);
    hPeak[j]->SetFillColor(lcol[1]);
    hPeak[j]->Draw("same");
    hSideband[j]->SetLineWidth(2);
    hSideband[j]->SetFillColor(lcol[2]);
    hSideband[j]->Draw("same");
  }
  fout.close();
  
  char OutCan[100];
  sprintf(OutCan,"fitKsMeson_sbBgd%i_%s.gif",iType,suffix);
  can1->Print(OutCan);
  sprintf(OutCan,"fitKsMeson_sbBgd%i_%s.eps",iType,suffix);
  can1->Print(OutCan);
}

// 
// ListYields - read yields for all targets from text files, calculate errors, and print to new file 
//                  
//                  suffix = output file suffix
//                  iType = target region [0], after target [1]
//
void ListYields(Int_t Nbins=100, char *suffix="std", Int_t iType=0)
{
  const Int_t nPol = 4;

  Int_t i,j;
  Int_t n = 0;

  Int_t tgtA[nTgt-2];
  Float_t Cts[nTgt-2][nPol];
  Float_t eCts[nTgt-2][nPol];
  
  char fname[500];
  char oname[500];

  for(i=0;i<nTgt-2;i++){
    sprintf(fname,"KsMeson_justBgd/%ibins/fitKsMeson_justBgd%i_%s-%s.lis",Nbins,iType,target[i+1],suffix);
    ifstream inPhi(fname);
    while(inPhi>>tgtA[i]>>Cts[i][0]>>Cts[i][1]>>Cts[i][2]>>Cts[i][3]){
      cout<<"Analyzing file "<<fname<<endl;
    }
    inPhi.close();
    
    for(j=0;j<nPol;j++){
      eCts[i][j] = sqrt(Cts[i][j]);
    }
  }

  for(i=0;i<nPol;i++){
    sprintf(oname,"mFit_yields-bc-%s-%i-pol%i.lis",suffix,iType,i+1);
    ofstream fout(oname);
    for(j=0;j<nTgt-2;j++){
      fout<<tgtA[j]<<"\t";
      fout<<Cts[i][j]<<"\t"<<eCts[i][j]<<"\t"<<eCts[i][j]<<"\t";
      fout<<endl;
    }
    fout.close();
  }
}

// 
// plotYields - read yields for all targets from a single text file and plot vs target 
//                  
//                  fname = input text file with yields from ListYields()
//                  RootFile = output ROOT file
//
void plotYields(char *fname="mFit.lis",char *RootFile="mFit_yields_Ks.root")
{
  Float_t tgt;
  Float_t err;
  Float_t y;
  Float_t y_eLo;
  Float_t y_eHi;
  Float_t radius;
  Int_t i, j;
  Int_t ibin;
  Int_t ntgt = 0;

  BookHists_Yields();

  ifstream fin(fname);
  while(fin>>tgt>>y>>y_eLo>>y_eHi){
    radius = R0*pow(tgt,1./3.);
    hYldVsA->Fill(tgt,y);
    ibin = hYldVsA->FindBin(tgt);
    hYldVsA->SetBinError(ibin,y_eLo);

    hYldVsR->Fill(radius,y);
    ibin = hYldVsR->FindBin(radius);
    hYldVsR->SetBinError(ibin,y);
  }
  WriteHists_Yields(RootFile);  
}

// 
// BookHists_Yields - create histograms for plotYields()
//                  
//
void BookHists_Yields()
{
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
void WriteHists_Yields(char *RootFile)
{
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

struct labels SelectHistogram(Int_t flag)
{
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
    exit;
    break;
  }
  return temp;
}

// 
// PrintHistogramOptions - list histograms available for SelectHistogram()
//                  
//
void PrintHistogramOptions()
{
  Int_t i;

  struct labels newLabels;

  for(i=0;i<19;i++){
    newLabels = SelectHistogram(i);
    cout<<"Flag "<<i<<"\t"<<newLabels.fPrefix<<"\t"<<newLabels.hPrefix<<endl;
  }
}

// 
// PrintHists - print contents of histogram created by plotYields()
//                  
//				RootFile = output from plotYields() 
//				hFlag = histogram selector
//
void PrintHists(char *RootFile="Ana.root", Int_t hFlag=0)
{
  Int_t i, j;

  char hist[50];
  
  TH1F *h[nP];
  TFile *f2 = new TFile(RootFile ,"READ");

  struct labels newLabels = SelectHistogram(hFlag);
  cout<<"Data for "<<newLabels.hPrefix<<" (flag = "<<hFlag<<")"<<endl;

  for(i=0;i<nP;i++){
    sprintf(hist,"%s%i",newLabels.hPrefix,i);
    h[i] = (TH1F*)f2->Get(hist);
    cout<<endl<<"Channel "<<i+1<<endl;
    for(j=0;j<h[i]->GetNbinsX();j++){
      if(h[i]->GetBinContent(j)>0.0){
	cout<<h[i]->GetBinCenter(j)<<"\t";
	cout<<h[i]->GetBinContent(j)<<" +/- "<<h[i]->GetBinError(j)<<endl;
      }
    }
  }
}

// 
// PlotHists - plot histogram created by plotYields()
//                  
//				RootFile = output from plotYields() 
//				hFlag = histogram selector
//
void PlotHists(char *RootFile="Ana.root",char *fSuffix="std",Int_t hFlag=0)
{
  Int_t i;

  char hist[50];
  
  TH1F *htmp;
  TFile *f2 = new TFile(RootFile ,"READ");

  struct labels newLabels = SelectHistogram(hFlag);
  cout<<"Data for "<<newLabels.hPrefix<<" (flag = "<<hFlag<<")"<<endl;

  sprintf(cname,"can");
  TCanvas *can = new TCanvas(cname,"Plot Hists",0,0,600,600);
  
  gStyle->SetOptStat(0);
  can->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  can->SetBorderSize(5); 
  can->SetFillStyle(4000); 
  
  sprintf(hist,"%s",newLabels.hPrefix);
  htmp = (TH1F*)f2->Get(hist);
  
  gPad->SetLeftMargin(Lmar);
  gPad->SetRightMargin(Rmar);
  gPad->SetFillColor(0);
  htmp->SetTitle(0);
  htmp->GetXaxis()->SetTitle(newLabels.x);
  htmp->GetXaxis()->CenterTitle();
  htmp->GetYaxis()->SetTitle(newLabels.y);
  htmp->GetYaxis()->CenterTitle();
  htmp->GetYaxis()->SetTitleOffset(yoff);
  htmp->SetMarkerSize(1.75);
  htmp->SetMarkerStyle(mkr[0]);
  htmp->SetMarkerColor(lcol[0]);
  htmp->Draw("PE1");
  
  sprintf(OutCan,"%s_%s.gif",newLabels.fPrefix,fSuffix);
  can->Print(OutCan);
  sprintf(OutCan,"%s_%s.eps",newLabels.fPrefix,fSuffix);
  can->Print(OutCan);
}

// 
// FindFoilIndex - return index number of target foils
//                  
//				vz = z vertex position
//
Int_t FindFoilIndex(Float_t vz)
{
  Int_t i;
  Int_t ret = 0;
  
  for(i=0;i<nFoils;i++){
    if(vz>=zcut[i] && vz<zcut[i+1]) ret = i+1;
  }
  return ret;
}

// 
// FindVzBehindFoilIndex - return index number of target foil where vz is directly behind
//                  
//				vz = z vertex position
//
Int_t FindVzBehindFoilIndex(Float_t vz)
{
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
Int_t FindTargetIndex(Float_t vz, Int_t method)
{
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
// FindNoPb - determine if target is not Pb
//                  
//				vz = z vertex position
//				method = [0 = target],[1 = after target]
//
Int_t FindNoPb(Float_t vz, Int_t method)
{
  Int_t ret = 0;
  Int_t ii = (method==0) ? FindFoilIndex(vz) : FindVzBehindFoilIndex(vz);
  if(ii>0 && ii!=5) ret = 1;
  return ret;
}

// 
// TargetIndex2A - return target mass number from target index
//                  
//				index = target index
//
Int_t TargetIndex2A(Int_t index)
{
  Int_t ret = 0;
  switch(index){
  case 1: 
  case 6: 
    ret = 2; break;
  case 2: 
  case 7: 
  case 9: 
  case 11: 
  case 13: 
    ret = 12; break;
  case 12: 
    ret = 48; break;
  case 3: 
  case 8: 
    ret = 56; break;
  case 4: 
  case 10: 
    ret = 208; break;
  case 0: 
    ret = (2+12+48+56+208)/5; break;
  case 5: 
    ret = (2+12+48+56)/4; break;
  default: 
    ret = 0; break;    
  }
  return ret;
}
