#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

gROOT->Reset();


	    
void pic_getter(){
	Int_t i,j,k;
	TH1F *hist;
	char *tar;
	char *AAA;
	char *type;
	char *inName;
	inName = "stuff";

	char *rootFile = "AnaTwoPion.root";
	TFile *fd = new TFile(rootFile,"UPDATE");
        TCanvas *can1 = new TCanvas("can1","title",0,0,1280,720); // create the canvas

	// 3 AAA, 4 targets, 3 types

	for(i=0;i<3;i++){

	   switch (i){
	      case 0: 
	         AAA = "Angle";
	         break;
	      case 1:
	         AAA = "Around";
		 break;
	      case 2:
		 AAA = "After";
		 break;
	   }


	   for(j=0;j<4;j++){

	      switch (i){
	         case 0: 
	            tar = "2H";
	            break;
	         case 1:
	            tar = "C";
		    break;
	         case 2:
		    tar = "FeTi";
		    break;
	         case 3:
		    tar = "Pb";
		    break;
	      }


	      for(k=0;k<3;k++){

	         switch (k){
	            case 0: 
	               type = "Chi";
	               break;
	            case 1:
	               type = "Yield";
		       break;
	            case 2:
		       type = "Ratio";
		       break;
	         }

	         sprintf(inName,"%sHist_%s_%s", type,tar,AAA);


                 hist = (TH1F*)fd->Get(inName);
	         hist->Draw();

	         if(k == 1)can1-> SetLogy();
	         else can1-> SetLogy(0);
	         //create the image files
	         char OutCan[200];
	         sprintf(OutCan,"pictures/%s.gif",inName);
	         can1->Print(OutCan);

	      }
	   }
	}

	can1->Close();


}
