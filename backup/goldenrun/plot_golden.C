#include <iomanip>

using namespace std;

void plot_golden()
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
  gStyle->SetTitleOffset(0.85,"X");
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

  double runfilenum_array[1985], ratio_array[1985];
  TH2F *hist=new TH2F("Ratio Vs Filenum_hist","Ratio vs Filenum after cut;File num;Event/Charge",50,3650000,3662000,50,10000,20000);
  TH1F *hist1 = new TH1F("Ratio","Ratio after cut;Event/Charge",50,10000,20000);

  for(int i=0;i<1985;i++)
    {
      infile>>runnum>>filenum>>count_event>>count_phot>>ratio;
     //cout<<"filenum="<<filenum<<endl;
      file = runnum*100+filenum;
      //cout<<"file"<<file<<endl;
      hist->Fill(file,ratio);
      hist1->Fill(ratio);
      runfilenum_array[i]=runnum*100+filenum;
      ratio_array[i]=ratio;
      //cout<<"run"<<runfilenum_array[i]<<endl;
    }

  TGraph *gr=new TGraph(1985,runfilenum_array,ratio_array);
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
  for(int i=0;i<1985;i++)
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
