// nt2_plot_TwoPion.C
//
// macro to plot data trees for each g7 pi+pi- events
// 
// M H Wood, University of Massachusetts, Amherst
//
// Notes:
//
//
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();

Float_t Lmar = 0.125;
Float_t Rmar = 0.125;
Float_t yoff = 1.5;

Float_t dtlim = 1.002;
Float_t dzlim = 3.0;
Float_t rlim = 2.0;
Float_t zlo = -22.0;
Float_t zhi = 5.0;
Float_t Plim = 0.5;
Float_t Mlo_Ks = 0.48;
Float_t Mhi_Ks = 0.51;

Double_t zFoil_mean[8] = {-18.0,-12.3,-9.89,-7.42,-4.93,-2.46,0.02,2.44}; // foil peak mean
Double_t zFoil_sigma[8] = {1.17,0.258,0.249,0.252,0.235,0.234,0.245,0.231}; // foil peak sigma

const Int_t MHITS= 2;
char *piName[MHITS] = {"#pi^{+}","#pi^{-}"};
char *lepName[MHITS] = {"e^{+}","e^{-}"};

struct labelsTH1 {
  Int_t bins;
  Float_t lo;
  Float_t hi;
  char x[50];
  char y[50];
  char hname[50];
  char title[100];
  char fPrefix[50];
  char itemTree[50];
};

struct labelsTH2 {
  Int_t binsX;
  Float_t xlo;
  Float_t xhi;
  Int_t binsY;
  Float_t ylo;
  Float_t yhi;
  char x[50];
  char y[50];
  char hname[50];
  char title[100];
  char fPrefix[50];
  char itemTree[50];
};

void nt2_plot_TwoPion(char *fin="TwoPion.lis",char *item="Theta-Phi", Int_t iType=0, Int_t tnum=1, Int_t iDim=2, Int_t dEvents=1000, Int_t MaxEvents=1000)
{
  // set up cuts for projection from tree
  Int_t H_cut;
  Int_t C1_cut;
  Int_t C2_cut;
  Int_t C3_cut;
  Int_t C4_cut;
  Int_t Fe_cut;
  Int_t Ti_cut;
  Int_t Pb_cut;
  Int_t tgt_cut;
  Int_t posP_cut;
  Int_t negP_cut;
  Int_t vz_cut;
  Int_t dt_cut;
  Int_t dz_cut;
  Int_t r_cut;
  Int_t sec_cut;
  Int_t gold_cut;
  Int_t M_Ks_cut;
  Int_t cuts;
  Int_t acuts;

  Int_t i, j, ii, max;
  Int_t fnum = 0;
  Int_t ilogz = 0;

  Double_t z1, z2;
  Double_t value, r_pair, dz;
  
  TVector3 V3[MHITS]; // vertex 3-vectors
  TLorentzVector V4[MHITS]; // particle 4-vectors
  TLorentzVector pairV4; // particle pair 4-vector

  char outgif[100];      // generic output gif filename
  
  if(iDim==1){
    struct labelsTH1 myLabels;
    if(iType==0 || iType==1){
      myLabels = PionVariable(item,iType);
    }else{
      myLabels = PairVariable(item);
    }
  }else if(iDim==2){
      struct labelsTH2 myLabels = xyVariable(item,iType);
      TH2F *hist = new TH2F(myLabels.hname,myLabels.title,myLabels.binsX,myLabels.xlo,myLabels.xhi,myLabels.binsY,myLabels.ylo,myLabels.yhi);
  }else{
    cout<<"Wrong number of histogram dimensions!!!"<<endl;
    cout<<"Try again."<<endl;
    exit(0);
  }

  // data files contain the trees
  printf("Analyzing file %s\n",fin);  

  // declare the chain
  TChain chain("PipPimTree");

  Int_t ncols;
  Int_t nfiles = 0;
  char rootFile[500];
  FILE *in1 = fopen(fin,"r");
  while (1) {
    ncols = fscanf(in1,"%s",rootFile);
    if (ncols<0) break;
    chain.Add(rootFile);
    nfiles++;
    if(!(nfiles % 100)) cerr << nfiles << "\r";
  }
  fclose(in1);
  cerr<<"Total number of files analyzed = "<<nfiles<<endl;

  Int_t Evt, runno, nTRACKS, TRIG, nPART, nPOS, nNEG;
  Int_t nNEUT, nPIM, nPIP, nKM, nKP, nPROTON, nGAMMA;
  chain.SetBranchAddress("Event",&Evt);
  chain.SetBranchAddress("Run",&runno);
  chain.SetBranchAddress("Tracks",&nTRACKS);
  chain.SetBranchAddress("trig",&TRIG);
  chain.SetBranchAddress("nPart",&nPART);
  chain.SetBranchAddress("nQpos",&nPOS);
  chain.SetBranchAddress("nQneg",&nNEG);
  chain.SetBranchAddress("nQneut",&nNEUT);
  chain.SetBranchAddress("nPim",&nPIM);
  chain.SetBranchAddress("nPip",&nPIP);
  chain.SetBranchAddress("nKm",&nKM);
  chain.SetBranchAddress("nKp",&nKP);
  chain.SetBranchAddress("nProton",&nPROTON);
  chain.SetBranchAddress("nGamma",&nGAMMA);

  Float_t eGAM, vTIME, tPHO, vX, vY, vZ, pairX, pairY, pairZ;
  chain.SetBranchAddress("egam",&eGAM);
  chain.SetBranchAddress("vTime",&vTIME);
  chain.SetBranchAddress("tpho",&tPHO);
  chain.SetBranchAddress("vx",&vX);
  chain.SetBranchAddress("vy",&vY);
  chain.SetBranchAddress("vz",&vZ);
  chain.SetBranchAddress("vx_pair",&pairX);
  chain.SetBranchAddress("vy_pair",&pairY);
  chain.SetBranchAddress("vz_pair",&pairZ);

  Int_t numPI, partQ[MHITS], PID[MHITS], SEC[MHITS], isLEP[MHITS], isG7[MHITS];
  Int_t isEC[MHITS], isECT[MHITS], isCC[MHITS], isMM2[MHITS];
  Int_t gold[MHITS], ecFID[MHITS], scID[MHITS], ccHITSTAT[MHITS],ccSTAT[MHITS];
  Int_t ccPMT[MHITS], ccSEG[MHITS], ccPHI[MHITS], lacSTAT[MHITS];
  Float_t partX[MHITS], partY[MHITS], partZ[MHITS], pE[MHITS], pX[MHITS];
  Float_t pY[MHITS], pZ[MHITS], tPROP[MHITS], ecTIME[MHITS], ecTOT[MHITS];
  Float_t ecIN[MHITS], ecOUT[MHITS], ecX[MHITS], ecY[MHITS], ecZ[MHITS];
  Float_t scTIME[MHITS], scLEN[MHITS], timing[MHITS], timing_e[MHITS];
  Float_t scVTpi[MHITS], scVTpi_pi[MHITS], scVTpi_e[MHITS], tofM2[MHITS];
  Float_t scX[MHITS], scY[MHITS], scZ[MHITS], ccTIME[MHITS], ccTHETA[MHITS];
  Float_t ccDTHETA[MHITS], ccNPE[MHITS], ccQF[MHITS], ccCS[MHITS], ccX[MHITS];
  Float_t ccY[MHITS], ccZ[MHITS], lacTOT[MHITS], lacIN[MHITS], lacX[MHITS];
  Float_t lacY[MHITS], lacZ[MHITS];
  chain.SetBranchAddress("Npi",&numPI);
  chain.SetBranchAddress("x",&partX);
  chain.SetBranchAddress("y",&partY);
  chain.SetBranchAddress("z",&partZ);
  chain.SetBranchAddress("Q",&partQ);
  chain.SetBranchAddress("Pid",&PID);
  chain.SetBranchAddress("Sec",&SEC);
  chain.SetBranchAddress("E",&pE);
  chain.SetBranchAddress("Px",&pX);
  chain.SetBranchAddress("Py",&pY);
  chain.SetBranchAddress("Pz",&pZ);
  chain.SetBranchAddress("Tprop",&tPROP);
  chain.SetBranchAddress("IsLep",&isLEP);
  chain.SetBranchAddress("IsLepG7",&isG7);
  chain.SetBranchAddress("IsLepG7ec",&isEC);
  chain.SetBranchAddress("IsLepG7ect",&isECT);
  chain.SetBranchAddress("IsLepG7cc",&isCC);
  chain.SetBranchAddress("IsLepG7mm2",&isMM2);
  chain.SetBranchAddress("Golden",&gold);
  chain.SetBranchAddress("EC_time",&ecTIME);
  chain.SetBranchAddress("EC",&ecTOT);
  chain.SetBranchAddress("ECin",&ecIN);
  chain.SetBranchAddress("ECout",&ecOUT);
  chain.SetBranchAddress("ECx",&ecX);
  chain.SetBranchAddress("ECy",&ecY);
  chain.SetBranchAddress("ECz",&ecZ);
  chain.SetBranchAddress("ECfid",&ecFID);
  chain.SetBranchAddress("SC_time",&scTIME);
  chain.SetBranchAddress("scLen",&scLEN);
  chain.SetBranchAddress("scId",&scID);
  chain.SetBranchAddress("Timing",&timing);
  chain.SetBranchAddress("Timing_e",&timing_e);
  chain.SetBranchAddress("scvT_pi",&scVTpi);
  chain.SetBranchAddress("scvT_pi_pi",&scVTpi_pi);
  chain.SetBranchAddress("scvT_pi_e",&scVTpi_e);
  chain.SetBranchAddress("TOF_MassSq",&tofM2);
  chain.SetBranchAddress("SCx",&scX);
  chain.SetBranchAddress("SCy",&scY);
  chain.SetBranchAddress("SCz",&scZ);
  chain.SetBranchAddress("CCstat",&ccSTAT);
  chain.SetBranchAddress("CC_time",&ccTIME);
  chain.SetBranchAddress("CChit_stat",&ccHITSTAT);
  chain.SetBranchAddress("CCtheta",&ccTHETA);
  chain.SetBranchAddress("CCdtheta",&ccDTHETA);
  chain.SetBranchAddress("CCnpe",&ccNPE);
  chain.SetBranchAddress("CC_QF",&ccQF);
  chain.SetBranchAddress("CC_CS",&ccCS);
  chain.SetBranchAddress("CCpmt",&ccPMT);
  chain.SetBranchAddress("CCseg",&ccSEG);
  chain.SetBranchAddress("CCphimatch",&ccPHI);
  chain.SetBranchAddress("CCx",&ccX);
  chain.SetBranchAddress("CCy",&ccY);
  chain.SetBranchAddress("CCz",&ccZ);
  chain.SetBranchAddress("LAC",&lacTOT);
  chain.SetBranchAddress("LACin",&lacIN);
  chain.SetBranchAddress("LACstat",&lacSTAT);
  chain.SetBranchAddress("LACx",&lacX);
  chain.SetBranchAddress("LACy",&lacY);
  chain.SetBranchAddress("LACz",&lacZ);

  // loop over events
  Int_t nevents = (Int_t)chain.GetEntries();
  printf("events = %i\n",nevents);

  if(MaxEvents){
    max = MaxEvents;
  } else{
    max = nevents;
  }

  // create the histogram object
  if(iDim==1){
    TH1F *hist = new TH1F(myLabels.hname,myLabels.title,myLabels.bins,myLabels.lo,myLabels.hi);
  }else if(iDim==2){
      TH2F *hist = new TH2F(myLabels.hname,myLabels.title,myLabels.binsX,myLabels.xlo,myLabels.xhi,myLabels.binsY,myLabels.ylo,myLabels.yhi);
  }else{
    cout<<"Wrong number of histogram dimensions!!!"<<endl;
    cout<<"Try again."<<endl;
    exit(0);
  }

  for(j=0;j<max;j++){
    chain.GetEntry(j); // retrieve the event from the tree
    if(!(j%dEvents)) cerr<<j<<" of "<<nevents<<" total\r"; // event counter

    // z vertex difference
    dz = partZ[0] - partZ[1]; 

    // radial vertex (x,y) for the particle pair
    r_pair = sqrt(pairX*pairX + pairY*pairY);

    // pair kinematics
    for(ii=0;ii<MHITS;ii++){
      V3[ii].SetXYZ(partX[ii],partY[ii],partZ[ii]);
      V4[ii].SetPxPyPzE(pX[ii],pY[ii],pZ[ii],pE[ii]);
    }
    pairV4 = V4[0] + V4[1];
    
    // set up the cuts
    dz_cut = (abs(dz) <= dzlim);
    r_cut = (r_pair <= rlim);
    dt_cut = (abs(scTIME[0]-scTIME[1])<= dtlim);
    sec_cut = (SEC[0]!=SEC[1]);
    gold_cut = (gold[0]==1 && gold[1]==1);
    negP_cut = (V4[0].P() >= Plim);
    posP_cut = (V4[1].P() >= Plim);
    M_Ks_cut = (pairV4.M()>=Mlo_Ks && pairV4.M()<=Mhi_Ks);

    vz_cut = (vZ>=zlo && vZ<=zhi);
    for(i=0;i<8;i++){
      z1 = zFoil_mean[i] - 3.0*zFoil_sigma[i];
      z2 = zFoil_mean[i] + 3.0*zFoil_sigma[i];
      switch(i){
      case 0: H_cut = (vZ>=z1 && vZ<=z2); break;
      case 1: C1_cut = (vZ>=z1 && vZ<=z2); break;
      case 2: Fe_cut = (vZ>=z1 && vZ<=z2); break;
      case 3: C2_cut = (vZ>=z1 && vZ<=z2); break;
      case 4: Pb_cut = (vZ>=z1 && vZ<=z2); break;
      case 5: C3_cut = (vZ>=z1 && vZ<=z2); break;
      case 6: Ti_cut = (vZ>=z1 && vZ<=z2); break;
      case 7: C4_cut = (vZ>=z1 && vZ<=z2); break;
      }
    }
    tgt_cut=(H_cut||C1_cut||Fe_cut||C2_cut||Pb_cut||C3_cut||Ti_cut||C4_cut);
    
    // make the cuts
    switch (tnum){
    case 0: // golden cuts
      cuts = gold_cut; break;
    case 1: // all standard cuts
      cuts = dz_cut && r_cut && sec_cut && dt_cut && gold_cut; break;
    case 2: // no dz cut
      cuts = r_cut && sec_cut && dt_cut && gold_cut; break;
    case 3: // no r cuts
      cuts = dz_cut && sec_cut && dt_cut && gold_cut; break;
    case 4: // no timing cut
      cuts = r_cut && sec_cut && dz_cut && gold_cut; break;
    case 5: // no sector cut
      cuts = r_cut && dt_cut && dz_cut && gold_cut; break;
    case 6: // no timing and sector cuts
      cuts = r_cut && dz_cut && gold_cut; break;
    case 7: // all standard cuts with mom. cuts
      cuts = dz_cut && r_cut && sec_cut && dt_cut && gold_cut && posP_cut && negP_cut; break;
    case 8: // all standard cuts with anit-mom. cuts
      cuts = dz_cut && r_cut && sec_cut && dt_cut && gold_cut && posP_cut && negP_cut; break;
    case 9: // all standard cuts with Ks mass cut
      cuts = dz_cut && r_cut && sec_cut && dt_cut && gold_cut && M_Ks_cut; break;
    case 10: // all standard cuts with target cuts
      cuts = vz_cut && dz_cut && r_cut && sec_cut && dt_cut && gold_cut && tgt_cut; break;
    case 11: // all standard cuts with anti-target cuts
      cuts = vz_cut && dz_cut && r_cut && sec_cut && dt_cut && gold_cut && !tgt_cut; break;
    default:
      cuts = dz_cut && r_cut && sec_cut && dt_cut && gold_cut; break;
    }

    switch(item){
    case "vz": hist->Fill(vZ); break;
    case "vz_pair": hist->Fill(pairZ); break;
    case "rmPair": hist->Fill(pairV4.M()); break;
    case "rmPair_e": hist->Fill(pairV4.M()); break;
    case "r_pair": hist->Fill(r_pair); break;
    case "dz": hist->Fill(dz); break;
    case "diffSCtime": hist->Fill(scTime[0] - scTime[1]); break;
	case "Theta": hist->Fill(V3[iType].Theta()); break;
	case "z": hist->Fill(V3[iType].Z()); break;
	case "Phi": hist->Fill(V3[iType].Phi()); break;
	case "Sec": hist->Fill(SEC[iType]); break;
	case "E": hist->Fill(V4[iType].E()); break;
	case "P": hist->Fill(V4[iType].P()); break;
	case "CCnpe": hist->Fill(ccNPE[iType]); break;
	case "CC_QF": hist->Fill(ccQF[iType]); break;
	case "CC_CS": hist->Fill(ccCS[iType]); break;
	case "Timing": hist->Fill(timing[iType]); break;
	default:
      cout<<"Invalid variable name: "<<item<<endl;
      cout<<"Try again."<<endl;
      exit(0);
      break;
    }
  }

  // create canvas
  TCanvas *kcan = new TCanvas("kcan","g7 Data",0,0,600,600);
  kcan->SetBorderMode(1);  //Bordermode (-1=down, 0 = no border, 1=up)
  kcan->SetBorderSize(5);  
  gStyle->SetOptStat(0);

  //draw hist with vaious options
  gPad->SetLeftMargin(Lmar);
  gPad->SetRightMargin(Rmar);
  hist->SetFillColor(2);
  if(ilogz) gPad->SetLogz();
  hist->SetXTitle(myLabels.x);
  hist->GetXaxis()->CenterTitle();
  hist->SetYTitle(myLabels.y);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleOffset(yoff);
  if(iDim==2){
    hist->Draw("colz");
  }else{
    cout<<"Plotting for Dim="<<iDim<<endl;
    hist->Draw();
  }

  if(fnum==1){
    TF1 *fcut1 = new TF1("fcut1","(x-70.0)/1.75",70.0,ch1Hi);
    fcut1->SetLineColor(2);
    fcut1->Draw("same");
    TF1 *fcut2 = new TF1("fcut2","(-x+70.0)/1.75",70.0,ch1Hi);
    fcut2->SetLineColor(2);
    fcut2->Draw("same");
  }

  sprintf(outgif,"nt_TwoPion_%s_%d.gif",myLabels.fPrefix,tnum);
  kcan->Print(outgif);
  sprintf(outgif,"nt_TwoPion_%s_%d.eps",myLabels.fPrefix,tnum);
  kcan->Print(outgif);

} 

// peak is Breit Wigner
Double_t breitwigner(Double_t *x, Double_t *par){
  return TMath::Max(1.e-10,par[0]*TMath::BreitWigner(x[0],par[2],par[1]));
}

struct labelsTH2 xyVariable(char *item, Int_t iType)
{
  Float_t ch1Lo,ch1Hi,ch2Lo,ch2Hi;
  Int_t Nch1,Nch2;
 
  struct labelsTH2 temp;

  sprintf(temp.fPrefix,"%s%i",item,iType); // rename item for output file

  switch(item){
  case "Theta-Phi":
  case "Phi-Theta":
    sprintf(temp.itemTree,"Theta[%i]:Phi[%i]",iType,iType);
    Nch1=180; ch1Lo=-180.0; ch1Hi=180.0;
    Nch2=80; ch2Lo=0.0; ch2Hi=80.0;
    sprintf(temp.title,"%s  #theta_{lab.} vs #phi_{lab.}",piName[iType]);
    sprintf(temp.x,"%s #phi_{lab.} (deg)",piName[iType]);
    sprintf(temp.y,"%s #theta_{lab.} (deg)",piName[iType]);
    break;
  case "compSec":
    sprintf(temp.itemTree,"Sec[0]:Sec[1]");
    Nch1=6; ch1Lo=0.5; ch1Hi=6.0;
    Nch2=6; ch2Lo=0.5; ch2Hi=6.0;
    sprintf(temp.title,"Sector comparison");
    sprintf(temp.x,"%s Sector",piName[1]);
    sprintf(temp.y,"%s Sector",piName[0]);
    break;
  case "compP":
    sprintf(temp.itemTree,"P[0]:P[1]");
    Nch1=100; ch1Lo=0.0; ch1Hi=2.0;
    Nch2=100; ch2Lo=0.0; ch2Hi=2.0;
    sprintf(temp.title,"Momentum comparison");
    sprintf(temp.x,"%s Momentum",piName[1]);
    sprintf(temp.y,"%s Momentum",piName[0]);
    break;
  case "compSCtime":
    sprintf(temp.itemTree,"SC_time[0]:SC_time[1]");
    Nch1=500; ch1Lo=20.0; ch1Hi=120.0;
    Nch2=500; ch2Lo=20.0; ch2Hi=120.0;
    sprintf(temp.title,"TOF Timing comparison");
    sprintf(temp.x,"%s SC Time (ns)",piName[1]);
    sprintf(temp.y,"%s SC Time (ns)",piName[0]);
    break;
  case "compscvT_pi_pi":
    sprintf(temp.itemTree,"scvT_pi_pi[0]:scvT_pi_pi[1]");
    Nch1=500; ch1Lo=20.0; ch1Hi=120.0;
    Nch2=500; ch2Lo=20.0; ch2Hi=120.0;
    sprintf(temp.title,"TOF Timing comparison w/ Mass Constrained");
    sprintf(temp.x,"%s SC Vertex Time (ns)"),piName[1];
    sprintf(temp.y,"%s SC Vertex Time (ns)",piName[0]);
    break;
  case "scvT_pi_pi-scvT_pi_e": 
  case "scvT_pi_e-scvT_pi_pi":
    sprintf(temp.itemTree,"scvT_pi_e[%i]:scvT_pi_pi[%i]",iType,iType);
    Nch1=500; ch1Lo=20.0; ch1Hi=120.0;
    Nch2=500; ch2Lo=20.0; ch2Hi=120.0;
    sprintf(temp.title,"TOF Timing comparison w/ Mass Constrained");
    sprintf(temp.x,"%s SC Vertex Time (ns)",piName[iType]);
    sprintf(temp.y,"e^{%s} SC Vertex Time (ns)",piName[iType]);
    break;
  case "compscvT":
    sprintf(temp.itemTree,"scvT_pi[0]:scvT_pi[1]");
    Nch1=500; ch1Lo=20.0; ch1Hi=120.0;
    Nch2=500; ch2Lo=20.0; ch2Hi=120.0;
    sprintf(temp.title,"TOF Timing comparison");
    sprintf(temp.x,"%s SC Vertex Time (ns)",piName[1]);
    sprintf(temp.y,"%s SC Vertex Time (ns)",piName[0]);
    break;
  case "pTime-eTime": 
  case "eTime-pTime":
    sprintf(temp.itemTree,"Timing[%i]:Timing_e[%i]",iType,iType);
    Nch1=500; ch1Lo=-20.0; ch1Hi=20.0;
    Nch2=500; ch2Lo=-20.0; ch2Hi=20.0;
    sprintf(temp.title,"Vertex Timing comparison");
    sprintf(temp.x,"%s vTime (ns)",piName[iType]);
    sprintf(temp.y,"e^{%s} vTime (ns)",piName[iType]);
    break;
  case "EC-P": 
  case "P-EC":
    sprintf(temp.itemTree,"EC[%i]:P[%i]",iType,iType);
    Nch1=250; ch1Lo=0.0; ch1Hi=2.5;
    Nch2=200; ch2Lo=0.0; ch2Hi=1.0;
    sprintf(temp.title,"%s EC vs P",piName[iType]);
    sprintf(temp.x,"Momentum");
    sprintf(temp.y,"EC energy");
    break;
  case "ECP-P": 
  case "P-ECP":
    sprintf(temp.itemTree,"EC[%i]/P[%i]:P[%i]",iType,iType,iType);
    Nch1=250; ch1Lo=0.0; ch1Hi=2.5;
    Nch2=100; ch2Lo=0.0; ch2Hi=0.5;
    sprintf(temp.title,"%s EC vs P",piName[iType]);
    sprintf(temp.x,"Momentum");
    sprintf(temp.y,"EC/P");
    break;
  case "ECi-ECo": 
  case "ECo-ECi":
    sprintf(temp.itemTree,"ECout[%i]:ECin[%i]",iType,iType);
    Nch1=100; ch1Lo=0.01; ch1Hi=0.51;
    Nch2=150; ch2Lo=0.01; ch2Hi=0.31;
    sprintf(temp.title,"%s Results",piName[iType]);
    sprintf(temp.x,"EC_{in}");
    sprintf(temp.y,"EC_{out}");
    break;
  case "ECiP-ECoP":
  case "ECoP-ECiP":
    sprintf(temp.itemTree,"ECout[%i]/P[%i]:ECin[%i]/P[%i]",iType,iType,iType,iType);
    Nch1=100; ch1Lo=0.01; ch1Hi=0.31;
    Nch2=100; ch2Lo=0.01; ch2Hi=0.36;
    sprintf(temp.title,"%s Results",piName[iType]);
    sprintf(temp.x,"EC_{in}/P");
    sprintf(temp.y,"EC_{out}/P");
  case "ECx-ECy": 
  case "ECy-ECx":
    sprintf(temp.itemTree,"ECy[%i]:ECx[%i]",iType,iType);
    Nch1=100; ch1Lo=0.0; ch1Hi=400.0;
    Nch2=75; ch2Lo=-150.0; ch2Hi=150.0;
    sprintf(temp.title,"%s EC x vs y",piName[iType]);
    sprintf(temp.x,"x_{EC}");
    sprintf(temp.y,"y_{EC}");
  case "CCx-CCy": 
  case "CCy-CCx":
    sprintf(temp.itemTree,"CCy[%i]:CCx[%i]",iType,iType);
    Nch1=100; ch1Lo=0.0; ch1Hi=400.0;
    Nch2=75; ch2Lo=-150.0; ch2Hi=150.0;
    sprintf(temp.title,"%s CC x vs y",piName[iType]);
    sprintf(temp.x,"x_{CC}");
    sprintf(temp.y,"y_{CC}");
  case "SCx-SCy": 
  case "SCy-SCx":
    sprintf(temp.itemTree,"SCy[%i]:SCx[%i]",iType,iType);
    Nch1=100; ch1Lo=0.0; ch1Hi=400.0;
    Nch2=75; ch2Lo=-150.0; ch2Hi=150.0;
    sprintf(temp.title,"%s SC x vs y",piName[iType]);
    sprintf(temp.x,"x_{SC}");
    sprintf(temp.y,"y_{SC}");
    break;
  case "LACx-LACy": 
  case "LACy-LACx":
    sprintf(temp.itemTree,"LACy[%i]:LACx[%i]",iType,iType);
    Nch1=100; ch1Lo=0.0; ch1Hi=400.0;
    Nch2=75; ch2Lo=-150.0; ch2Hi=150.0;
    sprintf(temp.title,"%s LAC x vs y",piName[iType]);
    sprintf(temp.x,"x_{LAC}");
    sprintf(temp.y,"y_{LAC}");
    break;
  default:
    cout<<"Invalid item "<<item<<". Try again!!!"<<endl;
    exit(0);
    break;
  }

  temp.binsX = Nch1;
  temp.xlo = ch1Lo;
  temp.xhi = ch1Hi;
  temp.binsY = Nch2;
  temp.ylo = ch2Lo;
  temp.yhi = ch2Hi;

  sprintf(temp.hname,"hist");

  return temp;
}

struct labelsTH1 PionVariable(char *item, Int_t iType)
{
  struct labelsTH1 temp;

  strcpy(temp.y,"Counts");
  sprintf(temp.fPrefix,"%s%i",item,iType); // rename item for output file

  switch(item){
  case "Theta":
    sprintf(temp.itemTree,"Theta[%i]",iType);
    temp.lo=0.0;        
    temp.hi=180.0;
    temp.bins=180;
    sprintf(temp.title,"%s  #theta_{lab.}",piName[iType]);
    sprintf(temp.x,"%s #theta_{lab.} (deg)",piName[iType]);
    break;
  case "z":
    sprintf(temp.itemTree,"z[%i]",iType);
    temp.lo=-40.0;        
    temp.hi=10.0;
    temp.bins=200;
    sprintf(temp.title,"%s  z Vertex",piName[iType]);
    sprintf(temp.x,"z (cm)");
    break;
  case "Phi":
    sprintf(temp.itemTree,"Phi[%i]",iType);
    temp.lo=-180.0;        
    temp.hi=180.0;
    temp.bins=180;
    sprintf(temp.title,"%s  #phi_{lab.}",piName[iType]);
    sprintf(temp.x,"%s #phi_{lab.} (deg)",piName[iType]);
    break;
  case "Sec":
    sprintf(temp.itemTree,"Sec[%i]",iType);
    temp.lo=0.5;        
    temp.hi=6.5;
    temp.bins=6;
    sprintf(temp.title,"%s Sector",piName[iType]);
    sprintf(temp.x,"%s Sector",piName[iType]);
    break;
  case "E":
    sprintf(temp.itemTree,"E[%i]",iType);
    temp.lo=0.0;        
    temp.hi=2.0;
    temp.bins=100;
    sprintf(temp.title,"%s Energy (GeV)",piName[iType]);
    sprintf(temp.x,"%s Energy (GeV)",piName[iType]);
    break;
  case "P":
    sprintf(temp.itemTree,"P[%i]",iType);
    temp.lo=0.0;        
    temp.hi=2.0;
    temp.bins=100;
    sprintf(temp.title,"%s Momentum",piName[iType]);
    sprintf(temp.x,"%s Momentum",piName[iType]);
    break;
  case "CCnpe":
    sprintf(temp.itemTree,"CCnpe[%i]",iType);
    temp.lo=-1.0;        
    temp.hi=50;
    temp.bins=101;
    sprintf(temp.title,"%s CC",piName[iType]);
    sprintf(temp.x,"%s Number of Photoelectrons",piName[iType]);
    break;
  case "CC_QF":
    sprintf(temp.itemTree,"CC_QF[%i]",iType);
    temp.lo=0.0;        
    temp.hi=1.0;
    temp.bins=100;
    sprintf(temp.title,"%s CC",piName[iType]);
    sprintf(temp.x,"%s CC Quality Factor (rad.)",piName[iType]);
    break;
  case "CC_CS":
    sprintf(temp.itemTree,"CC_CS[%i]",iType);
    temp.lo=-1.0;        
    temp.hi=1.0;
    temp.bins=100;
    sprintf(temp.title,"%s CC",piName[iType]);
    sprintf(temp.x,"%s CS factor",piName[iType]);
    break;
  case "Timing":
    sprintf(temp.itemTree,"Timing[%i]",iType);
    temp.lo=-20.0;        
    temp.hi=20.0;
    temp.bins=240;
    sprintf(temp.title,"%s Timing",piName[iType]);
    sprintf(temp.x,"%s Vertex Timing (ns)",piName[iType]);
    break;
  default:
    cout<<"Invalid variable name: "<<item<<endl;
    cout<<"Try again."<<endl;
    exit(0);
    break;
  }

  sprintf(temp.hname,"hist");

  return temp;
}

struct labelsTH1 PairVariable(char *item)
{
  struct labelsTH1 temp;

  strcpy(temp.y,"Counts");
  sprintf(temp.fPrefix,"%s",item); // rename item for output file

  switch(item){
  case "vz":
    sprintf(temp.itemTree,"%s",item);
    temp.lo=-40.0;        
    temp.hi=10.0;
    temp.bins=200;
    sprintf(temp.title,"z Vertex");
    sprintf(temp.x,"z (cm)");
    break;
  case "vz_pair":
    sprintf(temp.itemTree,"%s",item);
    temp.lo=-40.0;        
    temp.hi=10.0;
    temp.bins=200;
    sprintf(temp.title,"%s%s z Vertex",piName[0],piName[1]);
    sprintf(temp.x,"z (cm)");
    break;
  case "rmPair":
    sprintf(temp.itemTree,"%s",item);
    temp.lo=0.0;        
    temp.hi=1.0;
    temp.bins=200;
    sprintf(temp.title,"%s%s IM",piName[0],piName[1]);
    sprintf(temp.x,"%s%s Invariant Mass (GeV)",piName[0],piName[1]);
    break;
  case "rmPair_e":
    sprintf(temp.itemTree,"%s",item);
    temp.lo=0.0;        
    temp.hi=1.0;
    temp.bins=200;
    sprintf(temp.title,"%s%s IM (with electron masses)",piName[0],piName[1]);
    sprintf(temp.x,"Invariant Mass (GeV)");
    break;
  case "r_pair":
    sprintf(temp.itemTree,"%s",item);
    temp.lo=0.0;        
    temp.hi=6.0;
    temp.bins=100;
    sprintf(temp.title,"%s%s radial vertex",piName[0],piName[1]);
    sprintf(temp.x,"%s%s Radial Vertex (cm)",piName[0],piName[1]);
    break;
  case "dz":
    sprintf(temp.itemTree,"%s",item);
    temp.lo=-6.0;        
    temp.hi=6.0;
    temp.bins=200;
    sprintf(temp.title,"%s%s #Deltaz",piName[0],piName[1]);
    sprintf(temp.x,"%s%s #Deltaz (cm)",piName[0],piName[1]);
    break;
  case "diffSCtime":
    sprintf(temp.itemTree,"SC_time[0]-SC_time[1]");
    temp.lo=-6.0;        
    temp.hi=6.0;
    temp.bins=200;
    sprintf(temp.title,"%s%s #Deltat",piName[0],piName[1]);
    sprintf(temp.x,"%s%s #Deltat (ns)",piName[0],piName[1]);
    break;
  default:
    cout<<"Invalid variable name: "<<item<<endl;
    cout<<"Try again."<<endl;
    exit(0);
    break;
  }

  sprintf(temp.hname,"hist");

  return temp;
}
