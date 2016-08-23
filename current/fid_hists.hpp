#ifndef FIDHIST_H_GUARD
#define FIDHIST_H_GUARD

int bins_theta = 500;
int bins_phi = 500;
float Theta_min = -360;
float Theta_max = 360;
float Phi_min = -360;
float Phi_max = 360;

char hname_fid[50];
char htitle_fid[500];
const int sector_num = 6;
TH2D *fid_sec_hist[sector_num];

TH2D *fid_hist = new TH2D("fid","fid",
	bins_phi, Theta_min, Theta_max, bins_theta, Phi_min, Phi_max);

void makeHists_fid(){
	for (int jj = 0; jj < sector_num; jj++) {
		sprintf(hname_fid,"fid_sec%d",jj+1);
		sprintf(htitle_fid,"fid_sec%d",jj+1);
		fid_sec_hist[jj] = new TH2D(hname_fid, htitle_fid,
			bins_phi/2, Theta_min, Theta_max, bins_theta/2, Phi_min, Phi_max);
	}

}

void Fill_fid(double theta, double phi, int sector){
	fid_hist->Fill(theta,phi);
	fid_sec_hist[sector]->Fill(theta,phi);
}

void Fid_Write(){
	fid_hist->SetXTitle("#theta");
	fid_hist->SetYTitle("#phi");

	fid_hist->Write();

	for (int jj = 0; jj < sector_num; jj++) {
		fid_sec_hist[jj]->SetYTitle("#Deltat");
		fid_sec_hist[jj]->Write();
	}
}

/*
void delta_T_canvas(){
	TCanvas* can_dt[sector_num][3];
	char can_name[50];
	char * P_PIP_E;
	for (int particle_i = 0; particle_i < 3; particle_i++) {
		for (int sec_i = 0; sec_i < sector_num; sec_i++) {
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
*/


#endif