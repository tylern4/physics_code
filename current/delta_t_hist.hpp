/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DELTA_T_HIST_H_GUARD
#define DELTA_T_HIST_H_GUARD
//
//Histogram declarations, fills, and write
//
//
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

	delta_t_mass_P->Write();
	delta_t_mass_P_PID->Write();
	delta_t_mass_PIP->Write();
	delta_t_mass_PIP_PID->Write();
	delta_t_mass_electron->Write();
	delta_t_mass_electron_PID->Write();
	delta_t_mass_positron->Write();
	delta_t_mass_positron_PID->Write();
}

void delta_t_slices_Write(){
	Cuts delta_t_cut[3][num_points];
	double fit_dt_min = -1.0;
	double fit_dt_max = 1.0;
	for (int j = 0; j < 3; j++) {
		for (int jj = 0; jj < num_points; jj++) {
			if(j != 2) delta_t_cut[j][num_points].FitGaus(delta_t_hist[j][jj],fit_dt_min,fit_dt_max);
			delta_t_hist[j][jj]->SetYTitle("#Deltat");
			delta_t_hist[j][jj]->Write();
		}
	}
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

void delta_t_Fill(double momentum, int ID, int charge, double delta_t_proton, double delta_t_pip,double delta_t_electron){
	for (int jj = 0; jj < num_points; jj++) {
		if(momentum > jj * bin_width && momentum <= (jj+1) * bin_width){
			if(charge == 1) {
				delta_t_hist[0][jj]->Fill(delta_t_proton);
				delta_t_hist[1][jj]->Fill(delta_t_pip);
			}
			if(charge == -1) delta_t_hist[2][jj]->Fill(delta_t_electron);
		}
	}
}

void delta_t_sec_pad(double momentum, int ID, int charge,
			double delta_t_proton, double delta_t_pip,double delta_t_electron,
			int sc_sector, int sc_paddle){

	if(charge == 1) {
		delta_t_sec_pad_hist[0][sc_sector-1][sc_paddle-1]->Fill(momentum,delta_t_proton);
		delta_t_sec_pad_hist[1][sc_sector-1][sc_paddle-1]->Fill(momentum,delta_t_pip);
	}
	if(charge == -1) delta_t_sec_pad_hist[2][sc_sector-1][sc_paddle-1]->Fill(momentum,delta_t_electron);

}

void delta_T_canvas(){
	TCanvas* can_dt[sc_sector_num][3];
	char can_name[50];
	char * P_PIP_E;
	for (int particle_i = 0; particle_i < 3; particle_i++) {
		for (int sec_i = 0; sec_i < sc_sector_num; sec_i++) {
			if(particle_i == 0) P_PIP_E = "Proton";
			if(particle_i == 1) P_PIP_E = "Pip";
			if(particle_i == 2) P_PIP_E = "Electron";

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
