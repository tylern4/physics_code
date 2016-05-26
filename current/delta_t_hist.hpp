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

TH2D *delta_t_mass_P = new TH2D("delta_t_mass_P","#Deltat assuming mass of proton", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);
TH2D *delta_t_mass_P_PID = new TH2D("delta_t_mass_P_PID","#Deltat assuming mass of proton with PID proton", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);

TH2D *delta_t_mass_PIP = new TH2D("delta_t_mass_PIP","#Deltat assuming mass of #pi^{+}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);
TH2D *delta_t_mass_PIP_PID = new TH2D("delta_t_mass_PIP_PID","#Deltat assuming mass of #pi^{+} with PID #pi^{+}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);

TH2D *delta_t_mass_electron = new TH2D("delta_t_mass_electron","#Deltat assuming mass of e^{-}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);
TH2D *delta_t_mass_electron_PID = new TH2D("delta_t_mass_electron_PID","#Deltat assuming mass of e^{-} with PID e^{-}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);

TH2D *delta_t_mass_positron = new TH2D("delta_t_mass_postitron","#Deltat assuming mass of e^{+}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);
TH2D *delta_t_mass_positron_PID = new TH2D("delta_t_mass_postitron_PID","#Deltat assuming mass of e^{+} with PID e^{+}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);


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

void delta_t_Fill(double momentum, double delta_t, int option){

	switch(option){
		case 1 :	delta_t_mass_P_PID->Fill(momentum,delta_t); 
					break;
		case 2 :	delta_t_mass_PIP_PID->Fill(momentum,delta_t);
					break;
		case 3 :	delta_t_mass_P->Fill(momentum,delta_t); 
					break;
		case 4 :	delta_t_mass_PIP->Fill(momentum,delta_t);
					break;
	}
}

#endif
