#include <iostream>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"

void MakeTimeDeviations(int startrun, int endrun) {

  //Path to stage 0 root;
  // std::string fpath = "/home/cprokop/CJP/DANCE_Analysis/stage1_root";
   std::string fpath = "/home/cprokop/CJP/DANCE_Analysis/stage0_root";
  
  //First read in the time deviations global file used for the stage 0 analysis
  ifstream TimeDevGlobal;
  TimeDevGlobal.open("TimeDeviations.txt");
  
  int detector[160];
  double timedeviation[160];
  
  for(int eye=0; eye<160; eye++) {
    TimeDevGlobal>>detector[eye]>>timedeviation[eye];
    //  cout<<eye<<"  "<<detector[eye]<<"  "<<timedeviation[eye]<<endl;
  } 
  
  const int nruns = endrun-startrun+1;
  
  TFile *fin[nruns];
  TH2D *hTimeDevRel0[nruns];
  TH1D *hTimeDevProj[nruns][162];
  double TimeDeviations[nruns][162];
  ofstream fout[nruns];
  stringstream outfilename[nruns];

  int runcounter=0;
  for(int eye=startrun; eye<endrun+1; eye++) {
    // fin[runcounter] = new TFile(Form("%s/Stage1_Histograms_Run_%d_10ns_CW_2000ns_CBT_500ns_DEBT.root",fpath.c_str(),eye));
    fin[runcounter] = new TFile(Form("%s/Stage0_Histograms_Run_%d.root",fpath.c_str(),eye));
    hTimeDevRel0[runcounter] = (TH2D*)fin[runcounter]->Get("TimeDev");
    
    hTimeDevRel0[runcounter]->GetXaxis()->SetRangeUser(-10.0,10.0);
    
    int tdcounter=0;

    outfilename[runcounter].str();
    outfilename[runcounter] << "TimeDeviations_Run_" << eye << ".txt";
    fout[runcounter].open(outfilename[runcounter].str().c_str());

    for(int jay=0; jay<162; jay++) {
      hTimeDevProj[runcounter][jay] = hTimeDevRel0[runcounter]->ProjectionX(Form("hTimeDevProj_%d_%d",eye,jay),jay+1,jay+1);
       cout<<eye<<"  "<<jay<<"  "<<hTimeDevProj[runcounter][jay]->GetMean()<<endl;
      TimeDeviations[runcounter][jay+1] = hTimeDevProj[runcounter][jay]->GetMean();
      
      if(jay>0) {
	TimeDeviations[runcounter][jay] += TimeDeviations[runcounter][jay-1];
      }
      
      if(jay != 76 && jay != 86) {

	/*
	if(TMath::Abs(TimeDeviations[runcounter][jay]) > 1.0) {
	  std::string text;
	  cin>>text;
	}
	*/
	//	TimeDeviations[runcounter][jay] += timedeviation[tdcounter];
	fout[runcounter]<<jay<<"  \t"<<TimeDeviations[runcounter][jay] + timedeviation[tdcounter]<<"\n";
	tdcounter++;
	
      }
      
    }
    fout[runcounter].close();
    fin[runcounter]->Close();
    runcounter++;
    
  }
  

  
}
