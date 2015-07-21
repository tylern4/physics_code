
using namespace std;

//  dataHandeler
//
//	hopefully this works
//
void dataHandeler(char *fin="all.lis", char *RootFile_output="outFile.root", int MaxEvents=0, int dEvents=10000){
	gROOT->Reset();
	int current_event;
	int num_of_events;
	int total_events = 0;

	TLorentzVector *_e0, *_p0, *_e1;//, *_p1;

	_e0 = new TLorentzVector();
	_p0 = new TLorentzVector();
	_e1 = new TLorentzVector();
	_e0->SetPxPyPzE(0,0,E1D_E0,E1D_E0);
	_p0->SetPxPyPzE(0,0,0,MASS_P);

	TFile *myFile;
	TFile *RootOutputFile;
	TTree *myTree;
	int number_cols=0;
	int number_files = 0;
	char rootFile[500];

	RootOutputFile = new TFile(RootFile_output,"RECREATE");

	cout << "Analyzing file " << fin << endl;


	FILE *input_file = fopen(fin,"r");
	if (input_file == NULL) perror ("Error opening file");

	while (1){

		number_cols = fscanf(input_file,"%s",rootFile);
		if (number_cols<0) break;
		myFile = new TFile(rootFile, "READ");

		myTree = (TTree *)myFile->Get("h10");


		getBranches(myTree);

		num_of_events = (Int_t)myTree->GetEntries();

		current_event = 0;

		while(current_event<num_of_events){

			myTree->GetEntry(current_event);

			////////////if (current_event%10000 == 0)	cout<<current_event<<"/"<<num_of_events<<endl;

			#pragma omp parallel for
			for(int event_number = 0; event_number < gpart; event_number++)
			{

				Px = cx[event_number]*p[event_number];
				Py = cy[event_number]*p[event_number];
				Pz = cz[event_number]*p[event_number];

				x = vx[event_number];
				y = vy[event_number];
				z = vz[event_number];

				ID = id[event_number];


				FillHist();

			}

			current_event++; 		  	// increment event counter
			total_events++; 					// increment total event counter
		}

		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		number_files++; 						// increment file counter

	}
	RootOutputFile->cd();
	WriteHists();
	RootOutputFile->Write();
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout<<total_events<<" events in "<<number_files<< " files."<<endl; // print out stats
}

void PrintEverything(char *fin, char *RootFile_output){
	gROOT->Reset();
	Int_t current_event;
	Int_t num_of_events;
	Int_t total_events = 0;

	TFile *myFile;
	TFile *RootOutputFile;
	TTree *myTree;
	Int_t number_cols=0;
	Int_t number_files = 0;
	char rootFile[500];

	RootOutputFile = new TFile(RootFile_output,"RECREATE");

	FILE *input_file = fopen(fin,"r");

	if (input_file == NULL) perror ("Error opening file");

	while (1){

		number_cols = fscanf(input_file,"%s",rootFile); 

		if (number_cols<0) break;
		myFile = new TFile(rootFile, "READ");

		myTree = (TTree *)myFile->Get("h10");

		getBranches(myTree);

		num_of_events = (Int_t)myTree->GetEntries();

		current_event = 0; 

		while(current_event<1000){

			myTree->GetEntry(current_event);

			///if (current_event%100 == 0)	cout<<current_event<<"/"<<num_of_events<<endl;

			cout << "gpart	" << gpart << endl;

			for(int event_number = 0; event_number < gpart; event_number++){
				cout <<	"id: " <<	id[event_number] << 	endl;
				cout <<	"stat: " <<	stat[event_number] << 	endl;
				cout <<	"dc: " <<	dc[event_number] << 	endl;
				cout <<	"cc: " <<	cc[event_number] << 	endl;
				cout <<	"sc: " <<	sc[event_number] << 	endl;
				cout <<	"ec: " <<	ec[event_number] << 	endl;
				cout <<	"lec: " <<	lec[event_number] << 	endl;
				cout <<	"p: " <<	p[event_number] << 	endl;
				cout <<	"m: " <<	m[event_number] << 	endl;
				cout <<	"q: " <<	q[event_number] << 	endl;
				cout <<	"b: " <<	b[event_number] << 	endl;
				cout <<	"cx: " <<	cx[event_number] << 	endl;
				cout <<	"cy: " <<	cy[event_number] << 	endl;
				cout <<	"cz: " <<	cz[event_number] << 	endl;
				cout <<	"vx: " <<	vx[event_number] << 	endl;
				cout <<	"vy: " <<	vy[event_number] << 	endl;
				cout <<	"vz: " <<	vz[event_number] << 	endl;
			}

			cout << "dc_part:	" << dc_part << endl;

			for(int event_number = 0; event_number < dc_part; event_number++){
				cout << "dc_sect:" <<	dc_sect[event_number] << endl;   
				cout << "dc_trk:" <<	dc_trk[event_number] << endl;   
				cout << "dc_stat:" <<	dc_stat[event_number] << endl;   
				cout << "dc_xcs:" <<	dc_xsc[event_number] << endl;   
				cout << "dc_ysc:" <<	dc_ysc[event_number] << endl;   
				cout << "dc_zsc:" <<	dc_zsc[event_number] << endl;   
				cout << "dc_cxsc:" <<	dc_cxsc[event_number] << endl;   
				cout << "dc_cysc:" <<	dc_cysc[event_number] << endl;   
				cout << "dc_czsc:" <<	dc_czsc[event_number] << endl;   
				cout << "dc_xec:" <<	dc_xec[event_number] << endl;   
				cout << "dc_yec:" <<	dc_yec[event_number] << endl;   
				cout << "dc_zec:" <<	dc_zec[event_number] << endl;   
				cout << "dc_thecc:" <<	dc_thcc[event_number] << endl;   
				cout << "dc_c2:" <<	dc_c2[event_number] << endl;   
			}



			cout << "ec_part:	" << dc_part << endl;

			for(int event_number = 0; event_number < dc_part; event_number++){
				cout << "ec_stat:"	<< ec_stat[event_number] << endl;   
				cout << "ec_sect:"	<< ec_sect[event_number] << endl;   
				cout << "ec_whol:"	<< ec_whol[event_number] << endl;   
				cout << "ec_inst:"	<< ec_inst[event_number] << endl;   
				cout << "ec_oust:"	<< ec_oust[event_number] << endl;   
				cout << "etot:"	<< etot[event_number] << endl;   
				cout << "ec_ei:"	<< ec_ei[event_number] << endl;   
				cout << "ec_eo:"	<< ec_eo[event_number] << endl;   
				cout << "ec_t:"	<< ec_t[event_number] << endl;   
				cout << "ec_r:"	<< ec_r[event_number] << endl;   
				cout << "ec_x:"	<< ech_x[event_number] << endl;   
				cout << "ec_y:"	<< ech_y[event_number] << endl;   
				cout << "ec_z:"	<< ech_z[event_number] << endl;   
				cout << "ec_m2:"	<< ec_m2[event_number] << endl;   
				cout << "ec_m3:"	<< ec_m3[event_number] << endl;   
				cout << "ec_m4:"	<< ec_m4[event_number] << endl;   
				cout << "ec_c2:"	<< ec_c2[event_number] << endl;   
			}


			/*cout << "sc_part:	" << dc_part << endl;
			cout << "cc_part:	" << dc_part << endl;
			cout << "lac_part:	" << dc_part << endl; */
			cout << "	**		" << endl;

			current_event++; 		  	// increment event counter
			total_events++; 					// increment total event counter 
		}

		myTree->Delete(); 						// delete Tree object
		myFile->Close("R"); 					// close input ROOT file.  The R flag deletes TProcessIDs
		number_files++; 						// increment file counter

	}
	RootOutputFile->Close();
	fclose(input_file); 														// close file with input file list
	cout<<total_events<<" events in "<<number_files<< " files."<<endl; // print out stats

}