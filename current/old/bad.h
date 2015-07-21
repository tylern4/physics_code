#include "TTree.h"
#include "TROOT.h"
#include "TH2.h"
#include <TLorentzVector.h>
#include <TFile.h>
#include "TStyle.h"
#include "TCanvas.h"
#include <fstream>
#include "TF1.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "TGraph.h"
#include "TLine.h"
#include <iomanip>
using namespace std;

void golden_run(char *fin="all.lis", char *RootFile="outFile.root", Int_t MaxEvents=0, Int_t dEvents=10000)
{
	gStyle->SetPalette(1);
	gStyle->SetOptStat(10);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetStatY(0.90);
	gStyle->SetStatX(0.90);
	gStyle->SetStatW(0.35);
	gStyle->SetStatH(0.5);
	gStyle->SetLabelSize(0.04,"X");
	gStyle->SetLabelFont(62,"X");
	gStyle->SetTitleSize(0.05,"X");
	gStyle->SetTitleFont(62,"X");
	gStyle->SetTitleOffset(0.85,"X")
	gStyle->SetLabelSize(0.03,"Y");
	gStyle->SetLabelFont(62,"Y");
	gStyle->SetTitleSize(0.05,"Y");
	gStyle->SetTitleFont(62,"Y");
	gStyle->SetTitleOffset(0.75,"Y");
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleY(0.9);
	//gStyle->SetTitleTextColor(kRed);
	//gStyle->SetFillColor(kYellow);
	gStyle->SetHistFillColor(kYellow);
	
	ifstream infile("good_file_table");
	int runnum,filenum,file;
	double count_event,count_phot,ratio;

	double runfilenum_array[2500], ratio_array[2500];
	TH2F *hist=new TH2F("Ratio Vs Filenum_hist","Ratio vs Filenum after cut;File num;Event/Charge",50,3650000,3662000,50,10000,20000);
	TH1F *hist1 = new TH1F("Ratio","Ratio after cut;Event/Charge",50,10000,20000);

	for(int i=0;i<2500;i++) {
		infile>>runnum>>filenum>>count_event>>count_phot>>ratio;
		// cout<<"filenum="<<filenum<<endl;
		file = runnum*100+filenum;
		//cout<<"file"<<file<<endl;
		hist->Fill(file,ratio);
		hist1->Fill(ratio);
		runfilenum_array[i]=runnum*100+filenum;
		ratio_array[i]=ratio;
		//cout<<"run"<<runfilenum_array[i]<<endl;
	}

	TGraph *gr=new TGraph(2500,runfilenum_array,ratio_array);
	gr->SetNameTitle("Ratio Vs Filenum","Ratio Vs Filenum after cut");
	gr->SetMinimum(10000);
	gr->SetMaximum(20000);
	gr->GetXaxis()->SetLabelSize(0.04);
	gr->GetXaxis()->SetLabelFont(22);
	gr->GetXaxis()->SetTitle("File number");
	gr->GetXaxis()->SetTitleSize(0.06);
	gr->GetXaxis()->SetTitleOffset(0.75);
	gr->GetXaxis()->CenterTitle();
	gr->GetXaxis()->SetTitleFont(22);
	gr->GetYaxis()->SetLabelSize(0.05);
	gr->GetYaxis()->SetLabelFont(22);
	gr->GetYaxis()->SetTitle("Event/Charge");
	gr->GetYaxis()->SetTitleSize(0.06);
	gr->GetYaxis()->SetTitleOffset(0.75);
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitleFont(22);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetMarkerColor(kRed);

	TF1 *gaus1=new TF1("gaus1","gaus",10000,20000);
	double mean1,sigma1;
	int n_sigma=3;

	TCanvas *canv=new TCanvas("Ratio Vs Filenum","Ratio Vs Filenum after cut",1000,600);
	gPad->SetLogz();
	hist->Draw("colz");
	canv->SaveAs("ratio vs filenum-----after cut.ps");

	TCanvas *canv_hist1=new TCanvas("Ratio dustribution","Ratio distribution",1000,600);
	gPad->SetLogz();
	hist1->Fit(gaus1,"R");
	hist1->Draw();
	mean1=(double)gaus1->GetParameter(1);
	sigma1=(double)gaus1->GetParameter(2);
	TLine *line1_1=new TLine(mean1-n_sigma*sigma1,0,mean1-n_sigma*sigma1,100);
	line1_1->SetLineColor(kBlue);
	line1_1->SetLineWidth(2);
	line1_1->Draw();
	TLine *line1_2=new TLine(mean1+n_sigma*sigma1,0,mean1+n_sigma*sigma1,100);
	line1_2->SetLineColor(kBlue);
	line1_2->SetLineWidth(2);
	line1_2->Draw();
	canv_hist1->SaveAs("Ratio distribution after cut.png");

	TCanvas *canv1=new TCanvas("Ratio Vs Filenum --- graph","Ratio Vs Filenum after cut graph",1000,600);
	gPad->SetLogz();
	gr->Draw("ap");
	canv1->SaveAs("ratio vs filenum after cut graph.png");
	/*ofstream outfile1("good_file_ratio");

	ifstream infile_again("count_event");
	for(int i=0;i<2500;i++)
		{
			infile_again>>runnum>>filenum>>count_event>>count_phot>>ratio;
		 // cout<<"filenum"<<filenum<<endl;

		if(ratio>=mean1-n_sigma*sigma1 && ratio<=mean1+n_sigma*sigma1)
		 outfile1<<setw(20)<<setiosflags(ios::left)<<runnum<<setw(20)<<setiosflags(ios::left)<<filenum<<setw(20)
			 <<setiosflags(ios::left)<<count_event<<setw(20)<<setiosflags(ios::left)<<count_phot<<setw(20)<<setiosflags(ios::left)<<ratio<<endl;
		/* if(filenum<10){
			outfile1<<"a1ntp_"<<runnum<<"_pass1.a0"<<filenum<<".rzn.root"<<endl;
		 }
		 else outfile1<<"a1ntp_"<<runnum<<"_pass1.a"<<filenum<<".rzn.root"<<endl;

		}

	outfile1.close();*/

}

void count_after_cut(char *fin="all.lis", char *RootFile="outFile.root") {
	//Set constant variable for mass
	string temp_filename,str_runnum,str_filenum;
	TFile *rootfile;
	TTree *tree;

	if(!infile) {
			cout<<"Can't open infile "<<endl;
			exit(1);
	}

	ofstream outfile("count_event",ios::out);

	while(1) {

		cout<<temp_filename<<endl;
		rootfile=new TFile(temp_filename.c_str());
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

			if ((q_temp > 0.)){
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
					//cout<<"pion"<<endl;
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
		str_runnum=temp_filename.substr(temp_filename.length()-24,5);
		str_filenum=temp_filename.substr(temp_filename.length()-11,2);
		cout<<str_runnum<<"    "<<str_filenum<<"    "<<n_event<<"totalQ="<<totalQ<<endl;
		if(n_event != 0 && totalQ != 0){
			outfile<<setw(20)<<setiosflags(ios::left)<<str_runnum<<setw(20)<<setiosflags(ios::left)<<str_filenum<<setw(20)<<setiosflags(ios::left)<<n_event<<setw(20)<<setiosflags(ios::left)<<totalQ<<setw(20)<<setiosflags(ios::left)<<n_event/totalQ<<endl;
	 	}
	}
	outfile.close();
}
