/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DELTA_T_HIST_H_GUARD
#define DELTA_T_HIST_H_GUARD
#include <TGraph.h>
const int N_SIGMA = 3;
//
//Histogram declarations, fills, and write
// j -> type: 0=>Proton,1=>Pip,2=>Electron
// jj -> Fit point
int bins_dt = 500;
int bins_p = 500;
float P_min = 0;
float P_max = 3.5;
float Dt_min = -10;
float Dt_max = 10;

char hname[50];
char htitle[500];
const int num_points = 20;
TH1D *delta_t_hist[3][num_points];
const double bin_width = (P_max - P_min)/num_points;

const int sc_sector_num = 6;
const int sc_paddle_num = 48;
TH2D *delta_t_sec_pad_hist[3][sc_sector_num][sc_paddle_num];

TH2D *delta_t_mass_P = new TH2D("delta_t_mass_P","#Deltat assuming mass of proton",
	bins_p, P_min, P_max, bins_dt, Dt_min, Dt_max);
TH2D *delta_t_mass_P_PID = new TH2D("delta_t_mass_P_PID","#Deltat assuming mass of proton with PID proton",
	bins_p, P_min, P_max, bins_dt, Dt_min, Dt_max);

TH2D *delta_t_mass_PIP = new TH2D("delta_t_mass_PIP","#Deltat assuming mass of #pi^{+}",
	bins_p, P_min, P_max, bins_dt, Dt_min, Dt_max);
TH2D *delta_t_mass_PIP_PID = new TH2D("delta_t_mass_PIP_PID","#Deltat assuming mass of #pi^{+} with PID #pi^{+}",
	bins_p, P_min, P_max, bins_dt, Dt_min, Dt_max);

TH2D *delta_t_mass_electron = new TH2D("delta_t_mass_electron","#Deltat assuming mass of e^{-}",
	bins_p, P_min, P_max, bins_dt, Dt_min, Dt_max);
TH2D *delta_t_mass_electron_PID = new TH2D("delta_t_mass_electron_PID","#Deltat assuming mass of e^{-} with PID e^{-}",
	bins_p, P_min, P_max, bins_dt, Dt_min, Dt_max);

TH2D *delta_t_mass_positron = new TH2D("delta_t_mass_postitron","#Deltat assuming mass of e^{+}",
	bins_p, P_min, P_max, bins_dt, Dt_min, Dt_max);
TH2D *delta_t_mass_positron_PID = new TH2D("delta_t_mass_postitron_PID","#Deltat assuming mass of e^{+} with PID e^{+}",
	bins_p, P_min, P_max, bins_dt, Dt_min, Dt_max);

void makeHists_delta_t(){
	for (int jj = 0; jj < num_points; jj++) {
		sprintf(hname,"delta_t_p_%d",jj);
		sprintf(htitle,"#Deltat P %d",jj);
		delta_t_hist[0][jj] = new TH1D(hname,htitle, bins_dt, Dt_min, Dt_max);

		sprintf(hname,"delta_t_pip_%d",jj);
		sprintf(htitle,"#Deltat #pi^{+} %d",jj);
		delta_t_hist[1][jj] = new TH1D(hname,htitle, bins_dt, Dt_min, Dt_max);

		sprintf(hname,"delta_t_electron_%d",jj);
		sprintf(htitle,"#Deltat electron %d",jj);
		delta_t_hist[2][jj] = new TH1D(hname,htitle, bins_dt, Dt_min, Dt_max);
	}

	for (int jj = 0; jj < sc_sector_num; jj++) {
		for (int jjj = 0; jjj < sc_paddle_num; jjj++) {
			sprintf(hname,"delta_t_p_sec%d_pad%d",jj+1,jjj+1);
			sprintf(htitle,"#Deltat P Sector %d Paddle %d",jj+1,jjj+1);
			delta_t_sec_pad_hist[0][jj][jjj] = new TH2D(hname, htitle,
				bins_p/2, P_min, P_max, bins_dt/2, Dt_min, Dt_max);

			sprintf(hname,"delta_t_pip_sec%d_pad%d",jj+1,jjj+1);
			sprintf(htitle,"#Deltat #pi^{+} Sector %d Paddle %d",jj+1,jjj+1);
			delta_t_sec_pad_hist[1][jj][jjj] = new TH2D(hname, htitle,
				bins_p/2, P_min, P_max, bins_dt/2, Dt_min, Dt_max);

			sprintf(hname,"delta_t_electron_sec%d_pad%d",jj+1,jjj+1);
			sprintf(htitle,"#Deltat electron Sector %d Paddle %d",jj+1,jjj+1);
			delta_t_sec_pad_hist[2][jj][jjj] = new TH2D(hname, htitle,
				bins_p/2, P_min, P_max, bins_dt/2, Dt_min, Dt_max);

		}
	}

}

void Fill_deltat_P(double momentum, double delta_t){
	delta_t_mass_P->Fill(momentum,delta_t);
}

void Fill_deltat_P_PID(double momentum, double delta_t){
	delta_t_mass_P_PID->Fill(momentum,delta_t);
}

void Fill_deltat_PIP(double momentum, double delta_t){
	delta_t_mass_PIP->Fill(momentum,delta_t);
}

void Fill_deltat_PIP_PID(double momentum, double delta_t){
	delta_t_mass_PIP_PID->Fill(momentum,delta_t);
}

void Fill_deltat_electron(double momentum, double delta_t){
	delta_t_mass_electron->Fill(momentum,delta_t);
}

void Fill_deltat_electron_PID(double momentum, double delta_t){
	delta_t_mass_electron_PID->Fill(momentum,delta_t);
}

void Fill_deltat_positron(double momentum, double delta_t){
	delta_t_mass_positron->Fill(momentum,delta_t);
}

void Fill_deltat_positron_PID(double momentum, double delta_t){
	delta_t_mass_positron_PID->Fill(momentum,delta_t);
}

void delta_t_slice_fit(){
	fit_functions.open("../skim/fit_functions.hpp");
	fit_functions << "//Auto Generated fit code from e1d" << endl;
	fit_functions << "#ifndef FIT_FUNCTIONS_H_GUARD\n#define FIT_FUNCTIONS_H_GUARD\n#include \"main.h\"\n" << endl;
	TF1 *peak = new TF1("peak","gaus", -1.5, 1.5);
	//[0]*exp(-[1]*x) + 
	char *func = "[2]*x + [3]";
	delta_t_mass_P->FitSlicesY(peak,0,-1,10,"QRG5");
	TH1D *delta_t_mass_P_0 = (TH1D*)gDirectory->Get("delta_t_mass_P_0");
	TH1D *delta_t_mass_P_1 = (TH1D*)gDirectory->Get("delta_t_mass_P_1");
	TH1D *delta_t_mass_P_2 = (TH1D*)gDirectory->Get("delta_t_mass_P_2");
	double x[500];
	double y_plus[500];
	double y_minus[500];
	int num = 0;
	for (int i = 0; i < 500; i++){
		if(delta_t_mass_P_1->GetBinContent(i) != 0){
			//Get momentum from bin center
			x[num] = (double)delta_t_mass_P_1->GetBinCenter(i);
			//mean + 3sigma
			y_plus[num] = (double)delta_t_mass_P_1->GetBinContent(i) + N_SIGMA * (double)delta_t_mass_P_2->GetBinContent(i);
			//mean - 3simga
			y_minus[num] = (double)delta_t_mass_P_1->GetBinContent(i) - N_SIGMA * (double)delta_t_mass_P_2->GetBinContent(i);
			num++;
		}
	}
	
	TGraph *P = new TGraph(num,x,y_plus);
	TGraph *M = new TGraph(num,x,y_minus);
	TF1 *Proton_Pos_fit = new TF1("Proton_Pos_fit",func);
	TF1 *Proton_Neg_fit = new TF1("Proton_Neg_fit",func);
	P->Fit(Proton_Pos_fit,"Q","",0.2,2);
	P->Write();
	M->Fit(Proton_Neg_fit,"Q","",0.2,2);
	M->Write();
	Proton_Pos_fit->Write();
	Proton_Neg_fit->Write();
	P->Draw("Same");
	M->Draw("Same");
	Proton_Pos_fit->Draw("Same");
	Proton_Neg_fit->Draw("Same");

	fit_functions << "double Proton_Pos_fit(double x){\n\treturn " << Proton_Pos_fit->GetExpFormula("P") << ";\n}"<< endl;
	fit_functions << "double Proton_Neg_fit(double x){\n\treturn " << Proton_Neg_fit->GetExpFormula("P") << ";\n}"<< endl;

	delta_t_mass_PIP->FitSlicesY(peak,0,-1,10,"QRG5");
	TH1D *delta_t_mass_PIP_0 = (TH1D*)gDirectory->Get("delta_t_mass_PIP_0");
	TH1D *delta_t_mass_PIP_1 = (TH1D*)gDirectory->Get("delta_t_mass_PIP_1");
	TH1D *delta_t_mass_PIP_2 = (TH1D*)gDirectory->Get("delta_t_mass_PIP_2");
	double x_pip[500];
	double y_plus_pip[500];
	double y_minus_pip[500];
	num = 0;
	for (int i = 0; i < 500; i++){
		if(delta_t_mass_PIP_1->GetBinContent(i) != 0){
			//Get momentum from bin center
			x_pip[num] = (double)delta_t_mass_PIP_1->GetBinCenter(i);
			//mean + 3sigma
			y_plus_pip[num] = (double)delta_t_mass_PIP_1->GetBinContent(i) + N_SIGMA * (double)delta_t_mass_PIP_2->GetBinContent(i);
			//mean - 3simga
			y_minus_pip[num] = (double)delta_t_mass_PIP_1->GetBinContent(i) - N_SIGMA * (double)delta_t_mass_PIP_2->GetBinContent(i);
			num++;
		}
	}
	
	TGraph *P_pip = new TGraph(num,x_pip,y_plus_pip);
	TGraph *M_pip = new TGraph(num,x_pip,y_minus_pip);
	TF1 *Pip_Pos_fit = new TF1("Pip_Pos_fit",func);
	TF1 *Pip_Neg_fit = new TF1("Pip_Neg_fit",func);
	P_pip->Fit(Pip_Pos_fit,"Q","",0.1,1.75);
	P_pip->Write();
	M_pip->Fit(Pip_Neg_fit,"Q","",0.1,1.75);
	M_pip->Write();
	Pip_Pos_fit->Write();
	Pip_Neg_fit->Write();
	P_pip->Draw("Same");
	M_pip->Draw("Same");
	Pip_Pos_fit->Draw("Same");
	Pip_Neg_fit->Draw("Same");

	fit_functions << "double Pip_Pos_fit(double x){\n\treturn " << Pip_Pos_fit->GetExpFormula("P") << ";\n}"<< endl;
	fit_functions << "double Pip_Neg_fit(double x){\n\treturn " << Pip_Neg_fit->GetExpFormula("P") << ";\n}"<< endl;
	fit_functions << "#endif\n" << endl;
}

void delta_t_Write(){
	delta_t_mass_P->SetXTitle("Momentum (GeV)");
	delta_t_mass_P->SetYTitle("#Deltat");
	delta_t_mass_PIP->SetXTitle("Momentum (GeV)");
	delta_t_mass_PIP->SetYTitle("#Deltat");
	delta_t_mass_P_PID->SetXTitle("Momentum (GeV)");
	delta_t_mass_P_PID->SetYTitle("#Deltat");
	delta_t_mass_PIP_PID->SetXTitle("Momentum (GeV)");
	delta_t_mass_PIP_PID->SetYTitle("#Deltat");
	delta_t_mass_electron->SetXTitle("Momentum (GeV)");
	delta_t_mass_electron->SetYTitle("#Deltat");
	delta_t_mass_electron_PID->SetXTitle("Momentum (GeV)");
	delta_t_mass_electron_PID->SetYTitle("#Deltat");
	delta_t_mass_positron->SetXTitle("Momentum (GeV)");
	delta_t_mass_positron->SetYTitle("#Deltat");
	delta_t_mass_positron_PID->SetXTitle("Momentum (GeV)");
	delta_t_mass_positron_PID->SetYTitle("#Deltat");

	delta_t_slice_fit();

	delta_t_mass_P->Write();
	delta_t_mass_P_PID->Write();



	delta_t_mass_PIP->Write();
	delta_t_mass_PIP_PID->Write();
	delta_t_mass_electron->Write();
	delta_t_mass_electron_PID->Write();
	delta_t_mass_positron->Write();
	delta_t_mass_positron_PID->Write();
}

void delta_t_Fill(double momentum, int charge, double delta_t_proton, double delta_t_pip,double delta_t_electron){
	for (int jj = 0; jj < num_points; jj++) {
		if(momentum > jj * bin_width && momentum <= (jj+1) * bin_width){
			if(charge == 1 && !std::isnan(delta_t_proton) && !std::isnan(delta_t_pip)) {
				delta_t_hist[0][jj]->Fill(delta_t_proton);
				delta_t_hist[1][jj]->Fill(delta_t_pip);
			}
			if(charge == -1) delta_t_hist[2][jj]->Fill(delta_t_electron);
		}
	}
}

void delta_t_slices_Write(){
	Cuts delta_t_cut[3][num_points];
	double fit_dt_min = -1.0;
	double fit_dt_max = 1.0;
	for (int j = 0; j < 3; j++) {
		for (int jj = 0; jj < num_points; jj++) {
			if(j != 2) {
				delta_t_cut[j][num_points].FitGaus(delta_t_hist[j][jj],fit_dt_min,fit_dt_max);
				//cout << j << ',' << jj << ',' << delta_t_cut[j][num_points].mean << ',' << delta_t_cut[j][num_points].sigma << endl;
			}

			delta_t_hist[j][jj]->SetYTitle("#Deltat");
			delta_t_hist[j][jj]->Write();
		}
	}
}

void delta_t_sec_pad(double momentum, int charge,
			double delta_t_proton, double delta_t_pip,double delta_t_electron,
			int sc_sector, int sc_paddle){

	if(charge == 1) {
		delta_t_sec_pad_hist[0][sc_sector-1][sc_paddle-1]->Fill(momentum,delta_t_proton);
		delta_t_sec_pad_hist[1][sc_sector-1][sc_paddle-1]->Fill(momentum,delta_t_pip);
	}
	if(charge == -1) delta_t_sec_pad_hist[2][sc_sector-1][sc_paddle-1]->Fill(momentum,delta_t_electron);

}

void delta_t_sec_pad_Write(){
	for (int j = 0; j < 3; j++) {
		for (int jj = 0; jj < sc_sector_num; jj++) {
			for (int jjj = 0; jjj < sc_paddle_num; jjj++) {
				delta_t_sec_pad_hist[j][jj][jjj]->SetYTitle("#Deltat");
				delta_t_sec_pad_hist[j][jj][jjj]->Write();
			}
		}
	}
}

void delta_T_canvas(){
	TCanvas *can_dt[sc_sector_num][3];
	char can_name[50];
	char *P_PIP_E;
	for (int particle_i = 0; particle_i < 3; particle_i++) {
		for (int sec_i = 0; sec_i < sc_sector_num; sec_i++) {
			if(particle_i == 0) P_PIP_E = "Proton";
			else if(particle_i == 1) P_PIP_E = "Pip";
			else if(particle_i == 2) P_PIP_E = "Electron";

			sprintf(can_name, "Sector %d %s",sec_i+1,P_PIP_E);
			can_dt[sec_i][particle_i] = new TCanvas(can_name,can_name,1200,800);
			can_dt[sec_i][particle_i]->Divide(6, 8);
			for (int pad_i = 0; pad_i < sc_paddle_num; pad_i++) {
				can_dt[sec_i][particle_i]->cd((int)pad_i+1);
				delta_t_sec_pad_hist[particle_i][sec_i][pad_i]->Draw("same""colz");
			}
			can_dt[sec_i][particle_i]->Write();
		}
	}
}

#endif
