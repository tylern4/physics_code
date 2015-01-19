#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TMath.h"
#include "TCanvas.h"
//#include "vector"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

void count_after_cut() {
  ifstream infile("dir_filename"); //get input files
  string s,str_runnum,str_filenum; //make variables s="string name from input file list,temperary name for each time through loop" str_runnum=?? str_filenum=??
  TFile *rootfile;
  TTree *tree;

//if there is no input throw an error
  if(!infile) {  
    cout<<"Can't open infile "<<endl;
    exit(1);
    }

//make output file 
  ofstream outfile("count_event",ios::out);

//while loop to read through input file. Last line says end on it to break the loop
  while(infile>>s){
    
    if(s == "end") { //Break loop
      cout<<"end of list"<<endl;
      break;
    }

    cout<<s<<endl;
    rootfile=new TFile(s.c_str());
    tree=(TTree*)rootfile->Get("h10");
    //Set variables for branches
    //Trip_flag: identify if events are good or bad
    Float_t         q_l;
    Float_t         t_l;
    Float_t         tr_time;
    Float_t         rf_time1;
    Int_t           latch1;
    Int_t           hlsc;
    Int_t           intt;
    Int_t           gpart;
    Short_t         id[20];   //[gpart]
    Char_t          stat[20];   //[gpart]
    UChar_t         dc[20];   //[gpart]
    UChar_t         cc[20];   //[gpart]
    UChar_t         sc[20];   //[gpart]
    UChar_t         ec[20];   //[gpart]
    UChar_t         lec[20];   //[gpart]
    Float_t         p[20];   //[gpart]
    Char_t          q[20];   //[gpart]
    Float_t         b[20];   //[gpart]
    Float_t         cx[20];   //[gpart]
    Float_t         cy[20];   //[gpart]
    Float_t         cz[20];   //[gpart]
    Float_t         vx[20];   //[gpart]
    Float_t         vy[20];   //[gpart]
    Float_t         vz[20];   //[gpart]

    tree->SetBranchAddress("q_l", &q_l);
    tree->SetBranchAddress("t_l", &t_l);
    tree->SetBranchAddress("tr_time", &tr_time);
    tree->SetBranchAddress("rf_time1", &rf_time1);
    tree->SetBranchAddress("latch1", &latch1);
    tree->SetBranchAddress("hlsc", &hlsc);
    tree->SetBranchAddress("intt", &intt);
    tree->SetBranchAddress("gpart", &gpart);
    tree->SetBranchAddress("id", id);
    tree->SetBranchAddress("stat", stat);
    tree->SetBranchAddress("dc", dc);
    tree->SetBranchAddress("cc", cc);
    tree->SetBranchAddress("sc", sc);
    tree->SetBranchAddress("ec", ec);
    tree->SetBranchAddress("lec", lec);
    tree->SetBranchAddress("p", p);
    tree->SetBranchAddress("q", q);
    tree->SetBranchAddress("b", b);
    tree->SetBranchAddress("cx", cx);
    tree->SetBranchAddress("cy", cy);
    tree->SetBranchAddress("cz", cz);
    tree->SetBranchAddress("vx", vx);
    tree->SetBranchAddress("vy", vy);
    tree->SetBranchAddress("vz", vz);

    //Define variables
    float totalQ=0;
    float qcurr ;
    float qprev ;
    float deltaq ;
    float q_temp;
    TLorentzVector ni_4vec(0, 0, 0, 0);
    TVector3 ni_3vec(0, 0, 0);
    TLorentzVector n_4vec(0,0,0,0);
    TVector3 n_3vec(0,0,0);
    TVector3 q_3vec(0,0,0);
    ni_4vec.SetVectM(ni_3vec, 0.9396);
    TLorentzVector ei_4vec(0, 0, 0, 0);
    TVector3 ei_3vec(0, 0, 2.039);
    ei_4vec.SetVectM(ei_3vec, 0.000511);
    TVector3 ef_3vec(0, 0, 0);
    TLorentzVector ef_4vec(0, 0, 0, 0);
    TVector3 pionf_3vec(0, 0, 0);
    TLorentzVector pionf_4vec(0, 0, 0, 0);
    TVector3 pf_3vec(0, 0, 0);
    TLorentzVector pf_4vec(0, 0, 0, 0);
    TVector3 protonf_3vec(0, 0, 0);
    TLorentzVector protonf_4vec(0, 0, 0, 0);
    TLorentzVector w_4vec(0,0,0,0);
    TLorentzVector w_n_q_4vec(0,0,0,0);
    TLorentzVector ZERO_4vec(0,0,0,0);
    TLorentzVector pionf_CM_4vec(0,0,0,0);
    //TLorentzVector cmq_4vec(0,0,0,0);
    // TLorentzVector cmPip_4vec;
    TVector3 b_3vec(0,0,0);
    TLorentzVector m_4vec(0,0,0,0);

    Long64_t nentries=tree->GetEntries();

    int n_event=0;

    for(Long64_t i=0;i<nentries;i++){

      //Get the ith entry
      tree->GetEntry(i);
      //Get total farda cup charge
        q_temp=q_l;
        qcurr=q_temp;
        //cout<<"q_l="<<q_l<<endl;
           
        if ((q_temp>0.)){
          // cout<<"q_l"<<q_temp<<"qcurr"<<qcurr<<endl;        
          if(( qcurr>qprev)){
            deltaq = qcurr-qprev;
            totalQ += deltaq;
            cout<<"qcurr="<<qcurr<<"qprev="<<qprev<<"deltaq="<<deltaq<<endl;
            // cout<<"q_l"<<q_temp<<"qcurr"<<qcurr<<"qprev"<<qprev<<"deltaq"<<deltaq<<endl;
          }
          qprev = qcurr;                    
        } 
        //Cut through trip_flag
        for(int j=0;j<gpart;j++){
	      //Cut delta time of good photons for SC && ST
          if(id[0]== 11){
            ef_3vec.SetXYZ(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0]);
            ef_4vec.SetVectM(ef_3vec, 0.000511);
            //cout<<"electron"<<endl;
          }
          else if(id[j]== 2212){
            pf_3vec.SetXYZ(p[j]*cx[j],p[j]*cy[j],p[j]*cz[j]);
            pf_4vec.SetVectM(pf_3vec,0.9383);
            // cout<<"proton"<<endl;
          }
          else if(id[j] == -211 ){
            //if(dc[i]>0 && cc[i]>0 && sc[i]>0){
            pionf_3vec.SetXYZ(p[j]*cx[j],p[j]*cy[j],p[j]*cz[j]);
            pionf_4vec.SetVectM(pionf_3vec, 0.13957);
            // cout<<"pion"<<endl;
            //Pmom_pif_hist->Fill(p[i]);
          }
          if( ef_4vec != ZERO_4vec && pf_4vec != ZERO_4vec && pionf_4vec != ZERO_4vec){
            /* TLorentzVector q_4vec = ei_4vec - ef_4vec;
            q_3vec = ei_3vec - ef_3vec;
            if(pionf_4vec != ZERO_4vec){
            TLorentzVector n_4vec = q_4vec-pionf_4vec -pf_4vec ;
            n_3vec = q_3vec-pionf_3vec -pf_3vec ;
            TLorentzVector w_4vec = pf_4vec + pionf_4vec;  
            TLorentzVector w_n_q_4vec = ni_4vec + q_4vec;	*/   
            // cout<<"ef_E"<<ef_4vec.E()<<endl;  
            // if( (im_cor<1.116-3*0.001761) || (im_cor>1.116+3*0.001761) ) continue;
            //if( (mm<0.9404-3*0.01231) || (mm>0.9404+3*0.01231) ) continue;

            n_event++;
          } 
        }
      }	
      tree->Delete();
      rootfile->Close();
      str_runnum=s.substr(s.length()-24,5);
      str_filenum=s.substr(s.length()-11,2);
      cout<<str_runnum<<"    "<<str_filenum<<"    "<<n_event<<"totalQ="<<totalQ<<endl;
      if(n_event != 0 && totalQ != 0){ 
        outfile<<setw(20)<<setiosflags(ios::left)<<str_runnum<<setw(20)<<setiosflags(ios::left)<<str_filenum
        <<setw(20)<<setiosflags(ios::left)<<n_event<<setw(20)<<setiosflags(ios::left)<<totalQ<<setw(20)<<setiosflags(ios::left)<<n_event/totalQ<<endl;
      }    
    }
  outfile.close();
}
