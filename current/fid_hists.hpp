#ifndef FIDHIST_H_GUARD
#define FIDHIST_H_GUARD

int bins_theta = 500;
int bins_phi = 500;
float theta_min = 0;
float theta_max = 90;
float phi_min = -360/2.0;
float phi_max = 360/2.0;

char hname_fid[50];
char htitle_fid[500];
const int sector_num = 6;
TH2D *fid_sec_hist[sector_num];

TH2D *fid_hist = new TH2D("fid","fid",
	bins_phi, phi_min, phi_max, bins_theta, theta_min, theta_max);

void makeHists_fid(){
	float min_phi[6] = {0,60,120,-180,-120,-60};
	float max_phi[6] = {60,120,180,-120,-60,0};
	for (int sec = 0; sec < sector_num; sec++) {
		sprintf(hname_fid,"fid_sec%d",sec+1);
		sprintf(htitle_fid,"fid_sec%d",sec+1);
		fid_sec_hist[sec] = new TH2D(hname_fid, htitle_fid,
			bins_phi, min_phi[sec], max_phi[sec], bins_theta, theta_min, theta_max);
	}

}

void Fill_fid(double theta, double phi, int sector){
	fid_hist->Fill(phi,theta);
	fid_sec_hist[sector]->Fill(phi,theta);
}

void Fid_Write(){
	fid_hist->SetYTitle("#theta");
	fid_hist->SetXTitle("#phi");

	fid_hist->Write();

	for (int sec = 0; sec < sector_num; sec++) {
		fid_sec_hist[sec]->SetYTitle("#theta");
		fid_sec_hist[sec]->SetXTitle("#phi");
		fid_sec_hist[sec]->Write();
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