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
TH2D *delta_t_mass_PIP = new TH2D("delta_t_mass_PIP","#Deltat assuming mass of #pi^{+}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);
TH2D *delta_t_mass_P_PID = new TH2D("delta_t_mass_P_PID","#Deltat assuming mass of proton with PID proton", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);
TH2D *delta_t_mass_PIP_PID = new TH2D("delta_t_mass_PIP_PID","#Deltat assuming mass of #pi^{+} with PID #pi^{+}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);

TH2D *delta_t_mass_ELECTRON = new TH2D("delta_t_mass_ELECTRON","#Deltat assuming mass of Electron", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);
TH2D *delta_t_mass_ELECTRON_PIP = new TH2D("delta_t_mass_ELECTRON_PIP","#Deltat assuming mass of Electron with PID #pi^{+}", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);
TH2D *delta_t_mass_ELECTRON_P = new TH2D("delta_t_mass_ELECTRON_P","#Deltat assuming mass of Electron with PID Proton", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);

TH2D *delta_t_mass_PIP_E = new TH2D("delta_t_mass_PIP_E","#Deltat assuming mass of PIP with PID Electron", 
	bins_dt, P_min, P_max, bins_p, Dt_min, Dt_max);



inline void delta_t_Fill(double momentum, double delta_t, int option){

	switch(option){
		case 1 :	delta_t_mass_P_PID->Fill(momentum,delta_t);
					break;
		case 2 :	delta_t_mass_PIP_PID->Fill(momentum,delta_t);
					break;
		case 3 :	delta_t_mass_P->Fill(momentum,delta_t);
					break;
		case 4 :	delta_t_mass_PIP->Fill(momentum,delta_t);
					break;
		case 5 :	delta_t_mass_ELECTRON->Fill(momentum,delta_t);
					break;
		case 6 :	delta_t_mass_ELECTRON_PIP->Fill(momentum,delta_t);
					break;
		case 7 :	delta_t_mass_ELECTRON_P->Fill(momentum,delta_t);
					break;
		case 8 :	delta_t_mass_PIP_E->Fill(momentum,delta_t);
					break; 
	}
}

inline void delta_t_Write(){
	delta_t_mass_P->SetXTitle("Momentum (GeV)");
	delta_t_mass_P->SetYTitle("#Deltat");
	delta_t_mass_PIP->SetXTitle("Momentum (GeV)");
	delta_t_mass_PIP->SetYTitle("#Deltat");
	delta_t_mass_P_PID->SetXTitle("Momentum (GeV)");
	delta_t_mass_P_PID->SetYTitle("#Deltat");
	delta_t_mass_PIP_PID->SetXTitle("Momentum (GeV)");
	delta_t_mass_PIP_PID->SetYTitle("#Deltat");
	delta_t_mass_ELECTRON->SetXTitle("Momentum (GeV)");
	delta_t_mass_ELECTRON->SetYTitle("#Deltat");
	delta_t_mass_ELECTRON_PIP->SetXTitle("Momentum (GeV)");
	delta_t_mass_ELECTRON_PIP->SetYTitle("#Deltat");
	delta_t_mass_ELECTRON_P->SetXTitle("Momentum (GeV)");
	delta_t_mass_ELECTRON_P->SetYTitle("#Deltat");
	delta_t_mass_PIP_E->SetXTitle("Momentum (GeV)");
	delta_t_mass_PIP_E->SetYTitle("#Deltat");

	delta_t_mass_P->Write();
	delta_t_mass_PIP->Write();
	delta_t_mass_P_PID->Write();
	delta_t_mass_PIP_PID->Write();
	delta_t_mass_ELECTRON->Write();
	delta_t_mass_ELECTRON_PIP->Write();
	delta_t_mass_ELECTRON_P->Write();
	delta_t_mass_PIP_E->Write(); 
}

#endif
