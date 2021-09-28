////////////////////////////////////////////////////////////////////////
//                                                                    //
//   Software Name: DANCE Data Acquisition and Analysis Package       //
//     Subpackage: DANCE_Analysis                                     //
//   Identifying Number: C18105                                       // 
//                                                                    //
////////////////////////////////////////////////////////////////////////
//                                                                    //
//                                                                    //
// Copyright 2019.                                                    //
// Triad National Security, LLC. All rights reserved.                 //
//                                                                    //
//                                                                    //
//                                                                    //
// This program was produced under U.S. Government contract           //
// 89233218CNA000001 for Los Alamos National Laboratory               //
// (LANL), which is operated by Triad National Security, LLC          //
// for the U.S. Department of Energy/National Nuclear Security        //
// Administration. All rights in the program are reserved by          //
// Triad National Security, LLC, and the U.S. Department of           //
// Energy/National Nuclear Security Administration. The Government    //
// is granted for itself and others acting on its behalf a            //
// nonexclusive, paid-up, irrevocable worldwide license in this       //
// material to reproduce, prepare derivative works, distribute        //
// copies to the public, perform publicly and display publicly,       //
// and to permit others to do so.                                     //
//                                                                    //
// This is open source software; you can redistribute it and/or       //
// modify it under the terms of the GPLv2 License. If software        //
// is modified to produce derivative works, such modified             //
// software should be clearly marked, so as not to confuse it         //
// with the version available from LANL. Full text of the GPLv2       //
// License can be found in the License file of the repository         //
// (GPLv2.0_License.txt).                                             //
//                                                                    //
////////////////////////////////////////////////////////////////////////


//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  Cathleen E. Fry        *//
//*  cfry@lanl.gov          *//
//*  analyzer.cpp           *// 
//*  Last Edit: 04/06/21    *//  
//***************************//

//File includes
#include "analyzer.h"

//C/C++ Includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

//ROOT Includes
#include "TRandom.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TCutG.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

stringstream amsg;

//Events
DANCE_Event devent;      //DANCE Event
DANCE_Event devent2;     //Second DANCE Event for processing removed spectra
U235_Event u235event;    //U235 Event
He3_Event he3event;      //He3 Event
Li6_Event li6event;      //Li6 Event
Bkg_Event bkgevent;      //Background Monitor Event

//Histograms
TH2D *hCoinCAEN;                      //coincidence matrix

//Time Diagnostics
TH1D *hEventLength;                   //length in ns of the event (diagnostic)
TH2D *hEventLength_Etot;              //ESum vs Event Legnth
TH2D *hEventLength_MCr;              //Mcr vs Event Legnth
TH2D *hEventTimeDist_Etot;              //ESum vs Event Legnth

TH1D *hTimeBetweenDEvents;            //time between subsequent DANCE events (ns)
TH3D *hTimeBetweenDEvents_ESum_Mcr;   //DANCE ESum vs time between subsequent DANCE events for Mcr==1
TH1D *hTimeBetweenT0s;                //time between T0s (ns)
TH2D *hCrystalIDvsTOF;                //TOF for each crystal 
TH1D *hCrystalTOF;                    //TOF for all crytals
TH2D *hCrystalIDvsTOF_Corr;           //Corrected TOF for each crystal 
TH1D *hCrystalTOF_Corr; //TOF for all crytals
TH1D *hDANCE_Entries_per_T0;
TH1D *hDANCE_Events_per_T0;

//Time Deviations
TH2D *hTimeDev_Rel0;  //Time deviations of all crystals to crystal 0
TH2D *hTimeDev;  //Time deviations relative to adjacent crystals


//Physics Spectra
TH1D *hEn;                   //Neutron Energy from dance events
TH1D *hTOF;                  //TOF for dance events
TH2D *hTOF_Mcl;              //MCl vs TOF from dance events

TH1D *hEn_Corr;              //Neutron Energy from dance events
TH1D *hCrystal_En_Corr;      //Neutron Energy from each crystal dance events
TH2D *hECrystal_En_Corr;      //Neutron Energy from each crystal dance events
TH1D *hTOF_Corr;             //TOF for dance events
TH2D *hTOF_Mcl_Corr;         //MCl vs TOF for dance events


TH2D *hGamma_Mcr1;
TH2D *hGammaCalib_Mcr1;
TH2D* hGammaCalib_Mcr2_late;
TH2D* hGammaCalib_Mcr1_late;

TH2D* ID_TOF;

// QGated spectra
TH3F *En_Ecl_Mcl_QGated[10]; //Max is 10 QGates. 
TH3F *En_Ecr_Mcr_QGated[10]; //Max is 10 QGates. 
TH3F *ID_Ecr_Mcr_QGated[10]; //Max is 10 QGates. 
TH2D *hTOF_Mcl_QGated[10];

//Isomer spectra
TH1D *hIsomer_Prompt[10];   //Prompt singles TOF spectrum
TH1D *hIsomer_Delayed[10];  //Delayed singles TOF spectrum
TH1D *hIsomer_TDiff[10];    //Delayed - Prompt spectrum
vector<double> Isomer_Prompt[10]; //Storage for Prompt TOFs
vector<double> Isomer_Delayed[10]; //Storage for Delayed TOFs

//3D Histograms
TH3F *En_Esum_Mcl;
TH3F *En_Esum_Mcr;
TH3F *En_Esum_Mcr_Pileup;
TH3F *En_Esum_Mcr_NoPileup;

TH3F *hTOF_Esum_Mcl;
TH3F *hTOF_Esum_Mcr;

TH3F *hTOF_Esum_Mcr_Removed[Max_Gamma_Removed];
TH3F *hEn_Esum_Mcr_Removed[Max_Gamma_Removed];

TH3F *En_Ecr1_Ecr2_mcr2;
TH2F *En_Ecr1_mcr1;

TH3F *hEn_TimeBetweenCrystals_Mcr;  // En vs time between subsequent hits of the same crystal (ns) vs mcr

TH3F *hEn_Eg_Mcr;

//Beam Monitors
//U235 Monitor
TH1D *hU235_TOF;  //Raw TOF for U235 Monitor
TH1D *hU235_TOF_Corr; //Corrected TOF for U235 Monitor
TH2D *hU235_PH_TOF;
TH1D *hU235_PulseHeight;  //Energy for U235 Monitor
TH1D *hU235_En;  //Neutron Energy for U235 Monitor 
TH1D *hU235_En_Corr;  //Neutron Energy for U235 Monitor (From Corrected TOF)
TH1D *hU235_Time_Between_Events; //Time between U235 hits

//Li6 Monitor
TH1D *hLi6_TOF;  //Raw TOF for Li6 Monitor
TH1D *hLi6_TOF_Corr; //Corrected TOF for Li6 Monitor
TH1D *hLi6_PulseHeight;  //Energy for Li6 Monitor
TH1D *hLi6_En;  //Neutron Energy for Li6 Monitor 
TH1D *hLi6_En_Corr;  //Neutron Energy for Li6 Monitor (From Corrected TOF)
TH2D *hLi6_PSD; 
TH1D *hLi6_Time_Between_Events; //Time between Li6 hits

//Background Monitor
TH1D *hBkg_TOF;  //Raw TOF for Bkg Monitor
TH1D *hBkg_TOF_Corr; //Corrected TOF for Bkg Monitor
TH1D *hBkg_PulseHeight;  //Energy for Bkg Monitor
TH1D *hBkg_En;  //Neutron Energy for Bkg Monitor 
TH1D *hBkg_En_Corr;  //Neutron Energy for Bkg Monitor (From Corrected TOF)
TH2D *hBkg_PSD; 
TH1D *hBkg_Time_Between_Events; //Time between Bkg hits

//He3 Monitor
TH1D *hHe3_TOF;  //Raw TOF for He3 Monitor
TH1D *hHe3_TOF_Corr; //Corrected TOF for He3 Monitor
TH1D *hHe3_PulseHeight;  //Energy for He3 Monitor
TH1D *hHe3_En;  //Neutron Energy for He3 Monitor 
TH1D *hHe3_En_Corr;  //Neutron Energy for He3 Monitor (From Corrected TOF)
TH1D *hHe3_Time_Between_Events; //Time between He3 hits

// JU Histograms
TH1F *hU235_TOF_gated;
TH1F *hU235_TOF_long_gated;
TH1F *hLi6_TOF_gated;
TH1F *hLi6_TOF_long_gated;
TH1F *hHe3_TOF_gated;
TH1F *hHe3_TOF_long_gated;
TH1F *tof;
TH1F *tof_gated_QM;
TH1F *tof_gated_BM;
TH1F *tof_gated_QM_long;
TH1F *tof_gated_BM_long;
TH1D *esum;	// These esum histos included for convenience
TH1D *esum2;	// they are redundant - could project from 3D
TH1D *esum3;
TH1D *esum4;
TH1D *esum5;

TH2D* hGamma_Mcr2_1stex;
TH2D* hGammaCalib_Mcr2_1stex;

// HighRateDebug histos
TH1F* hID_resonancegated;
TH1F* hID_backgroundgated;
TH2F* ISlow_ID_mcr2;
//TH2F* En_ID;

/* VARIABLES */
double last_timestamp_devent[256];        //This keeps track of the last timestamp valid or not
double last_energy_devent[256];           //This keeps track of the last timestamp valid or not

double last_devent_timestamp;       //This keeps track of the last DANCE event timestamp valid or not
double last_valid_devent_timestamp; //This keeps track of the last DANCE event timestamp that was valid

//TMatrix Things
int reftoindex1[200];
int reftoindex2[200];
int index1[200];
int index2[200];

int totalindex; //Number of TMatrix Pairs read in

//DMatrix Things
int DetMat1[167];		      
int DetMat2[167];
int DetMat3[167];
int DetMat4[167];
int DetMat5[167];
int DetMat6[167];
int DetMat7[167];



//Diagnostics
uint32_t T0_Counter=0;
uint32_t DANCE_Entries_per_T0=0;
uint32_t DANCE_Events_per_T0=0;

vector<int> removed_crystals;


//This function goes through and makes the run-by-run time deviations
int Make_Time_Deviations(int RunNumber) {
  
  DANCE_Info("Analyzer","Making Time Deviations");
  
  //Stringstream for the outpout file name of td_out
  stringstream outfilename;
  outfilename.str();
  outfilename << TIMEDEV_DIR << "/TimeDeviations_Run_" << RunNumber << ".txt";
  
  //initialize the time deviation to 0.  Its cumulative...
  double time_deviation=0;

  //open the time deviations text file
  ofstream td_out;
  td_out.open(outfilename.str().c_str());
  
  //First detector has no offest
  td_out <<"0   \t 0 \n";
    
  //  cout<<"Total Indices: "<<totalindex<<endl;

  //Fill in the rest
  for(int eye=0; eye<totalindex; eye++) {
    // cout<<index1[eye]<<"  "<<index2[eye]<<endl;
    
    //Set the range of the Time Deviations 2D histogram to the proper bin 
    hTimeDev->GetYaxis()->SetRangeUser(index1[eye],index1[eye]+1);
    hTimeDev->GetXaxis()->SetRangeUser(-500,500);

    //Iteratively Close in on the proper range 
    for(int kay=3; kay<103; kay++) {
      hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-(103.0-kay), hTimeDev->GetMean()+(103.0-kay));
    }
    time_deviation += hTimeDev->GetMean();
    td_out <<index2[eye]<<"   \t"<<time_deviation<<"\n";
  }
  
  amsg.str("");
  amsg<<"Made Time Deviations for Run: "<<RunNumber;
  DANCE_Success("Analyzer",amsg.str());
  return 0;
}


int Read_DMatrix() {

  DANCE_Info("Analyzer","Making Detector Matrix");
  
  ifstream matrix_in(DMatrixFile);
  
  if(matrix_in.is_open()) {
    for(int eye=0; eye<167; eye++){
      matrix_in >> DetMat1[eye];
      matrix_in >> DetMat2[eye];
      matrix_in >> DetMat3[eye];
     matrix_in >> DetMat4[eye];
      matrix_in >> DetMat5[eye];
      matrix_in >> DetMat6[eye];
      matrix_in >> DetMat7[eye];
    }
    matrix_in.close();
    amsg.str("");
    amsg<<"Read in Detector Matrix File "<<DMatrixFile;
    DANCE_Success("Analyzer",amsg.str());
    return 0;
  }
  else {
    amsg.str("");
    amsg<<"Failed to Read in Detector Matrix File.  Exiting "<<DMatrixFile;
    DANCE_Error("Analyzer",amsg.str());
    return -1;
  }
  
}


int Read_TMatrix() {

  DANCE_Info("Analyzer","Making Timing Matrix");
  
  ifstream timemat;
  timemat.open(TMatrixFile);
    
  for(int i=0;i<200;i++){
    reftoindex1[i]=-1;
    reftoindex2[i]=-1;
  }
  int counter=0;

  if(timemat.is_open()) {
    while(!timemat.eof()) {
      
      //aaa and bbb are Detector IDs for neighboring crystals
      int aaa,bbb;
      timemat >> aaa >> bbb;

      index1[counter] = aaa;   
      index2[counter] = bbb;
      
      if(timemat.eof()) break;

      //for each pair of neighboring crystals there is a refereence to them.  
      //This way if there is a missing crystal there isnt a gap in the array.  
      reftoindex1[aaa]=counter;   //left entry of pair
      reftoindex2[bbb]=counter;   //right entry of pair
      
      // cout<<counter<<"  "<<aaa<<"  "<<bbb<<"  "<<reftoindex1[aaa]<<"  "<<reftoindex2[bbb]<<endl;
      //number of pairs read in
      counter++;
    }    
    timemat.close();
    amsg.str("");
    amsg<<"Read in Timing Matrix File "<<TMatrixFile;
    DANCE_Success("Analyzer",amsg.str());
  }
  else {
    amsg.str("");
    amsg<<"Failed to Read in Timing Matrix File.  Exiting "<<TMatrixFile;
    DANCE_Error("Analyzer",amsg.str());
    return -1;
  }
  
  return counter;
}


int Create_Analyzer_Histograms(Input_Parameters input_params) {
  
  DANCE_Info("Analyzer","Creating Histograms");

  //Make Histograms
  hCoinCAEN = new TH2D("CoinCAEN","CoinCAEN",162,0,162,162,0,162);  //coincidence matrix


  //Time Deviations
  hTimeDev_Rel0 = new TH2D("TimeDev_Rel0","TimeDev_Rel0",10000,-500,500,162,0,162);  //Time deviations relative to crystal 0
  hTimeDev = new TH2D("TimeDev","TimeDev",10000,-500,500,162,0,162);  //Time deviations

  //Diagnostics
  hEventLength = new TH1D("EventLength","EventLength",10000,0,100);
  hEventLength_Etot = new TH2D("EventLength_Etot","EventLength_Etot",1000,0,10,400,0,20);
  hEventLength_MCr = new TH2D("EventLength_Mcr","EventLength_Mcr",1000,0,10,30,0,30);
  hEventTimeDist_Etot = new TH2D("EventTimeDist_Etot","EventTimeDist_Etot",1000,0,10,400,0,20);

  hTimeBetweenDEvents = new TH1D("TimeBetweenDEvents","TimeBetweenDEvents",1000,0,10000);
#ifdef TurnOffGoFast
  hTimeBetweenDEvents_ESum_Mcr = new TH3D("TimeBetweenDEvents_ESum_Mcr","TimeBetweenDEvents_ESum_Mcr",1000,0,10000,500,0,10,20,0,20);
#endif
  hTimeBetweenT0s = new TH1D("TimeBetweenT0s","TimeBetweenT0s",1000000,0,100000000);  //Time difference between T0 in ns
  
  //RAW TOF
  hCrystalIDvsTOF = new TH2D("ID_CrystalTOF","ID_CrystalTOF",10000,0,1000000,162,0,162); //TOF for each crystal 
  hCrystalTOF = new TH1D("CrystalTOF","CrystalTOF",6000000,0,60000000);
  hTOF = new TH1D("TOF","TOF",6000000,0,60000000);
  //2D out to 1 ms
  hTOF_Mcl = new TH2D("Mcl_TOF","Mcl_TOF",10000,0,1000000,8,0,8);

  //Corrected TOF
  hCrystalIDvsTOF_Corr = new TH2D("ID_CrystalTOF_corr","ID_CrystalTOF_corr",10000,0,10000000,162,0,162); //TOF for each crystal 
  hCrystalTOF_Corr = new TH1D("CrystalTOF_corr","CrystalTOF_corr",600000,0,60000000);
  hTOF_Corr = new TH1D("TOF_corr","TOF_corr",6000000,0,60000000);
  //2D out to 1 ms
  hTOF_Mcl_Corr = new TH2D("Mcl_TOF_corr","Mcl_TOF_corr",10000,0,1000000,8,0,8);



  hDANCE_Entries_per_T0 = new TH1D("DEntriesPerT0","DANCE_Entries_per_T0",100000,0,100000);
  hDANCE_Events_per_T0 = new TH1D("DEventsPerT0","DANCE_Events_per_T0",100000,0,100000);
  
  //Gamma gated on Crystal Mult 1
  hGamma_Mcr1 = new TH2D("ISlow_ID_mcr1","ISlow_ID_mcr1",3500,0,70000,162,0,162);
  hGammaCalib_Mcr1 = new TH2D("ESlow_ID_mcr1","ESlow_ID_mcr1",2000,0.0,20.0,162,0,162);
   hGammaCalib_Mcr1_late = new TH2D("ESlow_ID_mcr1LateTOF","ESlow_ID_mcr1LateTOF",2000,0.0,20.0,162,0,162);
   hGammaCalib_Mcr2_late = new TH2D("ESlow_ID_mcr2LateTOF","ESlow_ID_mcr1LateTOF",2000,0.0,20.0,162,0,162);
    
  //Physics Histograms
  double x[5000];
  int NEbins=0;
  
  for(double lx=log10(NeutronE_From);lx<log10(NeutronE_To);lx=lx+(1./NeutronE_BinsPerDecade)){
    x[NEbins]=pow(10,lx);
    NEbins++;
  }
  NEbins--;
  
  int NoOfEnergyBins=GammaE_NoOfBins;
  double *EtotBins=new double[NoOfEnergyBins+1];
  double *Mbins=new double[21];
  double *Tbins=new double[3001];
  double *IDbins=new double[162];
  
  for(int i=0;i<21;i++){
    Mbins[i]=i;
    //   Mbins[i]=0.5+1.*i;
  };
  
  for(int i=0;i<3001;i++){
    Tbins[i]=i*10.0;
  };
  

  double DEGamma=(GammaE_To-GammaE_From)/GammaE_NoOfBins;
  for(int i=0;i<NoOfEnergyBins+1;i++){
    //EtotBins[i]=i*16./128.;
    EtotBins[i]=GammaE_From+i*DEGamma;
  };

  for (int i=0;i<162;i++){
    IDbins[i]=(double)i;
  }




  // ------ JU physics histograms (tof) ----------------------------------------------- 
  tof=new TH1F("tof","Time of flight [ns]",500000,-50000,49950000);
  tof_gated_QM = new TH1F("tof_gated_QM","TOF gated",100000,0,1000000);
  tof_gated_QM->GetXaxis()->SetTitle("TOF (ns)");
  tof_gated_BM = new TH1F("tof_gated_BM","TOF gated",100000,0,1000000);
  tof_gated_BM->GetXaxis()->SetTitle("TOF (ns)");
  tof_gated_QM_long = new TH1F("tof_gated_QM_long","TOF gated",150000,0,30000000);
  tof_gated_QM_long->GetXaxis()->SetTitle("TOF (ns)");
  tof_gated_BM_long = new TH1F("tof_gated_BM_long","TOF gated",150000,0,30000000);
  tof_gated_BM_long->GetXaxis()->SetTitle("TOF (ns)");
  esum  = new TH1D("Esum", "Esum", 400,0.0,20.0);
  esum2 = new TH1D("Esum2","Esum2",400,0.0,20.0);
  esum3 = new TH1D("Esum3","Esum3",400,0.0,20.0);
  esum4 = new TH1D("Esum4","Esum4",400,0.0,20.0);
  esum5 = new TH1D("Esum5","Esum5",400,0.0,20.0);
  // ----------------------------------------------------------------------------------
  
  //Beam Monitors
  hU235_TOF = new TH1D("U235_TOF","U235_TOF",600000,0,60000000);  //Raw TOF for U235 Monitor
  hU235_TOF_Corr = new TH1D("U235_TOF_corr","U235_TOF_corr",600000,0,60000000); //Corrected TOF for U235 Monitor

  hU235_PH_TOF = new TH2D("U235_ISlow_TOF","U235_ISlow_TOF",6000,0,60000000,250,0,100000);  //Raw PH vs TOF for U235 Monitor

  

  hU235_PulseHeight = new TH1D("U235_ISlow","U235_ISlow",10000,0,100000);  //Energy for U235 Monitor
  hU235_En = new TH1D("U235_En","U235_En",NEbins,x);  //Neutron Energy for U235 Monitor 
  hU235_En_Corr = new TH1D("U235_En_corr","U235_En_corr",NEbins,x);  //Neutron Energy for U235 Monitor (From Corrected TOF)
  hU235_Time_Between_Events = new TH1D("TimeBetweenU235Events","TimeBetweenU235Events",100000,0,10000000);

  hHe3_TOF = new TH1D("He3_TOF","He3_TOF",600000,0,60000000);  //Raw TOF for He3 Monitor
  hHe3_TOF_Corr = new TH1D("He3_TOF_corr","He3_TOF_corr",600000,0,60000000); //Corrected TOF for He3 Monitor
  hHe3_PulseHeight = new TH1D("He3_ISlow","He3_ISlow",10000,0,100000);  //Energy for He3 Monitor
  hHe3_En = new TH1D("He3_En","He3_En",NEbins,x);  //Neutron Energy for He3 Monitor 
  hHe3_En_Corr = new TH1D("He3_En_corr","He3_En_corr",NEbins,x);  //Neutron Energy for He3 Monitor (From Corrected TOF)
  hHe3_Time_Between_Events = new TH1D("TimeBetweenHe3Events","TimeBetweenHe3Events",100000,0,10000000);

  hLi6_TOF = new TH1D("Li6_TOF","Li6_TOF",600000,0,60000000);  //Raw TOF for Li6 Monitor
  hLi6_TOF_Corr = new TH1D("Li6_TOF_corr","Li6_TOF_corr",600000,0,60000000); //Corrected TOF for Li6 Monitor
  hLi6_PulseHeight = new TH1D("Li6_ISlow","Li6_ISlow",10000,0,100000);  //Energy for Li6 Monitor
  hLi6_PSD = new TH2D("Li6_ISlow_IFast","Li6_ISlow_IFast",600,0,60000,600,0,60000);  //Ifast vs late for Li6 Monitor
  hLi6_En = new TH1D("Li6_En","Li6_En",NEbins,x);  //Neutron Energy for Li6 Monitor 
  hLi6_En_Corr = new TH1D("Li6_En_corr","Li6_En_corr",NEbins,x);  //Neutron Energy for Li6 Monitor (From Corrected TOF)
  hLi6_Time_Between_Events = new TH1D("TimeBetweenLi6Events","TimeBetweenLi6Events",100000,0,10000000);

  hBkg_TOF = new TH1D("Bkg_TOF","Bkg_TOF",600000,0,60000000);  //Raw TOF for Bkg Monitor
  hBkg_TOF_Corr = new TH1D("Bkg_TOF_corr","Bkg_TOF_corr",600000,0,60000000); //Corrected TOF for Bkg Monitor
  hBkg_PulseHeight = new TH1D("Bkg_ISlow","Bkg_ISlow",10000,0,100000);  //Energy for Bkg Monitor
  hBkg_En = new TH1D("Bkg_En","Bkg_En",NEbins,x);  //Neutron Energy for Bkg Monitor 
  hBkg_En_Corr = new TH1D("Bkg_En_corr","Bkg_En_corr",NEbins,x);  //Neutron Energy for Bkg Monitor (From Corrected TOF)
  hBkg_Time_Between_Events = new TH1D("TimeBetweenBkgEvents","TimeBetweenBkgEvents",100000,0,10000000);


  // --- JU monitor histogram definitions ----------------------
  //     (keeping same name as previous for use in flux programs) (Also keeping as TH1F for compatability)
  hU235_TOF_gated = new TH1F("h_tof_U235_FF","U235 Gated",10000,0.0,1.0e6);	// gated TOF
  hU235_TOF_long_gated = new TH1F("h_tof_long_U235_FF","U235 Long Gated",45000,0.0,4.5e7);
  hLi6_TOF_gated = new TH1F("h_tof_Li6_Triton","Li6 Gated",10000,0.0,1.0e+6);	// gated TOF
  hLi6_TOF_long_gated = new TH1F("h_tof_long_Li6_Triton","Li6 Long Gated",45000,0.0,4.5e7);
  hHe3_TOF_gated  = new TH1F("h_tof_He3_All","He3 gated",10000,0.0,1.0e+6);	// gated TOF
  hHe3_TOF_long_gated = new TH1F("h_tof_long_He3_All","He3 Long Gated",45000,0.0,4.5e7);
  // ------------------------------------------------------------


  if(input_params.Analysis_Stage==1 || input_params.Read_Simulation==1) {

    hEn = new TH1D("En","En",NEbins,x); //Raw En
    hEn_Corr = new TH1D("En_corr","En_corr",NEbins,x);  //Corrected En 
    hCrystal_En_Corr = new TH1D("CrystalEn_corr","CrystalEn_corr",NEbins,x);  //Corrected En 
    hECrystal_En_Corr = new TH2D("Ecr_En_corr","Ecr_En_corr",NEbins,x,NoOfEnergyBins,EtotBins);  //Corrected En vs ECrystal

#ifdef HighRateDebug
    hID_backgroundgated = new TH1F("ID_background","hID_backgroundgated",162,0,162); 
    hID_resonancegated = new TH1F("ID_resonance","hID_resonancegated",162,0,162);  
    //En_ID=new TH2F("En_ID","En_ID",NEbins,x,162,IDbins);
    //this is broken a little and I don't know why 
    ISlow_ID_mcr2 = new TH2F("ISlow_ID_mcr2","ISlow_ID_mcr2",3500,0.0,70000,162,0,162);
#endif

    En_Esum_Mcl=new TH3F("En_Etot_Mcl","En_Etot_Mcl",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
    En_Esum_Mcr=new TH3F("En_Etot_Mcr","En_Etot_Mcr",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
    En_Esum_Mcr_Pileup=new TH3F("En_Etot_Mcr_PU","En_Etot_Mcr_PU",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
    En_Esum_Mcr_NoPileup=new TH3F("En_Etot_Mcr_noPU","En_Etot_Mcr_noPU",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);

    hTOF_Esum_Mcl=new TH3F("TOF_Etot_Mcl","TOF_Etot_Mcl",2000,0,1000000,NoOfEnergyBins,GammaE_From,GammaE_To,20,0,20);
    hTOF_Esum_Mcr=new TH3F("TOF_Etot_Mcr","TOF_Etot_Mcr",2000,0,1000000,NoOfEnergyBins,GammaE_From,GammaE_To,20,0,20);
#ifdef TurnOffGoFast    
    //Make a plot of energy 1 vs energy 2 for mcr==2 as a function of energy
    En_Ecr1_Ecr2_mcr2 = new TH3F("En_Ecr1_Ecr2_mcr2","En_Ecr1_Ecr2_mcr2",NEbins,x,NoOfEnergyBins,EtotBins,NoOfEnergyBins,EtotBins);

    //cout << En_Ecr1_Ecr2_mcr2->GetSize() << "\t" << (NEbins+2)*(NoOfEnergyBins+2)*(NoOfEnergyBins+2) << endl;

    En_Ecr1_mcr1 = new TH2F("En_Ecr1_mcr1","En_Ecr1_mcr1",NEbins,x,NoOfEnergyBins,EtotBins);

    hEn_TimeBetweenCrystals_Mcr = new TH3F("En_TimeBetweenCrystals_Mcr","En_TimeBetweenCrystals_Mcr",NEbins,x,3000,Tbins,20,Mbins);

    hEn_Eg_Mcr = new TH3F("En_Ecr_Mcr","En_Ecr_Mcr",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
#endif
    
#ifdef Make_Removed_Spectra
    for(int kay=0; kay<Max_Gamma_Removed; kay++) {
      //Same spectra but with one gamma ray at random thrown away
      hTOF_Esum_Mcr_Removed[kay]=new TH3F(Form("TOF_Etot_Mcr_%dremoved",kay),Form("TOF_Etot_Mcr_%dremoved",kay),2000,0,1000000,NoOfEnergyBins,GammaE_From,GammaE_To,20,0,20);
      hEn_Esum_Mcr_Removed[kay]=new TH3F(Form("En_Etot_Mcr_%dremoved",kay),Form("En_Etot_Mcr_%dremoved",kay),NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
    } 
#endif

    if(input_params.QGatedSpectra) {
      
      for (int kay=0; kay<input_params.NQGates; kay++) {
	En_Ecl_Mcl_QGated[kay]=new TH3F(Form("En_Ecl_Mcl_esum%d",kay),
					Form("En_Ecl_Mcl_esum_gated_%2.2f_%2.2f",input_params.QGates[2*kay],input_params.QGates[2*kay+1]),
					NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
      }

      for (int kay=0; kay<input_params.NQGates; kay++) {
	En_Ecr_Mcr_QGated[kay]=new TH3F(Form("En_Ecr_Mcr_esum%d",kay),
					Form("En_Ecr_Mcr_esum_gated_%2.2f_%2.2f",input_params.QGates[2*kay],input_params.QGates[2*kay+1]),
					NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
      }
      
      for (int kay=0; kay<input_params.NQGates; kay++) {
	ID_Ecr_Mcr_QGated[kay]=new TH3F(Form("ID_Ecr_Mcr_esum%d",kay),
					Form("ID_Ecr_Mcr_esum_gated_%2.2f_%2.2f",input_params.QGates[2*kay],input_params.QGates[2*kay+1]),
					162,0,162,400,0,20,20,0,20);
      }
      
      //QGated Mcl_TOF out to 1 ms
      for (int kay=0; kay<input_params.NQGates; kay++) {
	hTOF_Mcl_QGated[kay] = new TH2D(Form("Mcl_TOF_esum%d",kay),
					Form("Mcl_TOF_ESum_Gated_%2.2f_%2.2f",input_params.QGates[2*kay],input_params.QGates[2*kay+1]),
					2000,0,1000000,8,0,8);
      }
    }


    if(input_params.IsomerSpectra) {
      
      for(int isom=0; isom<input_params.NIsomers; isom++) {
	hIsomer_Prompt[isom] = new TH1D(Form("Isomer_prompt%d",isom),
					Form("Isomer_Prompt_QGated_%2.2f_%2.2f_MclGated_%d_%d_TOFGated_%2.2f_%2.2f",
					     input_params.IsomerPromptQGates[isom*2],input_params.IsomerPromptQGates[isom*2+1],
					     input_params.IsomerPromptMclGates[isom*2],input_params.IsomerPromptMclGates[isom*2+1],
					     input_params.IsomerPromptTOFGates[isom*2],input_params.IsomerPromptTOFGates[isom*2+1]),
					input_params.IsomerPromptTOFGates[isom*2+1] - input_params.IsomerPromptTOFGates[isom*2],
					input_params.IsomerPromptTOFGates[isom*2],
					input_params.IsomerPromptTOFGates[isom*2+1]);
	  
	  
	hIsomer_Delayed[isom] = new TH1D(Form("Isomer_delayed%d",isom),
					 Form("Isomer_Delayed_QGated_%2.2f_%2.2f_MclGated_%d_%d_TOFGated_%2.2f_%2.2f",
					      input_params.IsomerDelayedQGates[isom*2],input_params.IsomerDelayedQGates[isom*2+1],
					      input_params.IsomerDelayedMclGates[isom*2],input_params.IsomerDelayedMclGates[isom*2+1],
					      input_params.IsomerDelayedTOFGates[isom*2],input_params.IsomerDelayedTOFGates[isom*2+1]),
					 input_params.IsomerDelayedTOFGates[isom*2+1] - input_params.IsomerDelayedTOFGates[isom*2],
					 input_params.IsomerDelayedTOFGates[isom*2],
					 input_params.IsomerDelayedTOFGates[isom*2+1]);
	  
	hIsomer_TDiff[isom] = new TH1D(Form("Isomer_timeDiff%d",isom),
				       Form("Isomer_TDiff_Delayed_%d_minus_Prompt_%d",isom,isom),
				       (input_params.IsomerDelayedTOFGates[isom*2+1] - input_params.IsomerPromptTOFGates[isom*2+1])-(input_params.IsomerDelayedTOFGates[isom*2] - input_params.IsomerPromptTOFGates[isom*2]),
				       input_params.IsomerDelayedTOFGates[isom*2] - input_params.IsomerPromptTOFGates[isom*2],
				       input_params.IsomerDelayedTOFGates[isom*2+1] - input_params.IsomerPromptTOFGates[isom*2+1]);
	  
      } //End loop over number of isomers
    } //Check on isomers 
    
  }
  DANCE_Success("Analyzer","Created Histograms");
  return 0;
  
}


int Write_Analyzer_Histograms(TFile *fout, Input_Parameters input_params) {
  
  DANCE_Info("Analyzer","Writing Histograms");

  fout->cd();

  hCoinCAEN->Write();
  hTimeDev_Rel0->Write();
  hTimeDev->Write();

  hTimeBetweenT0s->Write();
  

  hDANCE_Entries_per_T0->Write();

  if(input_params.Analysis_Stage==1 || input_params.Read_Simulation == 1) {
    hEventLength->Write();
    hEventLength_Etot->Write();
    hEventLength_MCr->Write();
    hEventTimeDist_Etot->Write();
    hTimeBetweenDEvents->Write();
#ifdef TurnOffGoFast
    hTimeBetweenDEvents_ESum_Mcr->Write();
#endif
    hCrystalIDvsTOF_Corr->Write();
    hCrystalTOF_Corr->Write();
    hCrystalIDvsTOF->Write();
    hCrystalTOF->Write();
    hDANCE_Events_per_T0->Write();

#ifdef HighRateDebug
    hID_backgroundgated->Write();
    hID_resonancegated->Write();
    //En_ID->Write();
    ISlow_ID_mcr2->Write();
#endif

    hGamma_Mcr1->Write();
    hGammaCalib_Mcr1->Write();
    hGammaCalib_Mcr1_late->Write();
    hGammaCalib_Mcr2_late->Write();

    hTOF->Write();
    hTOF_Corr->Write();
    hTOF_Mcl->Write();
    hTOF_Mcl_Corr->Write();

    hEn->Write();
    hEn_Corr->Write();
    hCrystal_En_Corr->Write();
    hECrystal_En_Corr->Write();
    En_Esum_Mcl->Write();
    En_Esum_Mcr->Write();

    En_Esum_Mcr_Pileup->Write();
    En_Esum_Mcr_NoPileup->Write();

    hTOF_Esum_Mcl->Write();
    hTOF_Esum_Mcr->Write();
#ifdef TurnOffGoFast
    En_Ecr1_Ecr2_mcr2->Write();;
    En_Ecr1_mcr1->Write();
      
    hEn_TimeBetweenCrystals_Mcr->Write();

    hEn_Eg_Mcr->Write();
#endif

#ifdef Make_Removed_Spectra
    for(int kay=0; kay<Max_Gamma_Removed; kay++) {
      hTOF_Esum_Mcr_Removed[kay]->Write();
      hEn_Esum_Mcr_Removed[kay]->Write();
    }
#endif

    if(input_params.QGatedSpectra) {
      for (int kay=0; kay<input_params.NQGates; kay++) {
	En_Ecl_Mcl_QGated[kay]->Write();
      }
      for (int kay=0; kay<input_params.NQGates; kay++) {
	En_Ecr_Mcr_QGated[kay]->Write();
      }
      for (int kay=0; kay<input_params.NQGates; kay++) {
	ID_Ecr_Mcr_QGated[kay]->Write();
      }
      for (int kay=0; kay<input_params.NQGates; kay++) {
	hTOF_Mcl_QGated[kay]->Write();
      }
    }

    if(input_params.IsomerSpectra) {
      for (int kay=0; kay<input_params.NIsomers; kay++) {
	hIsomer_Prompt[kay]->Write();
  	hIsomer_Delayed[kay]->Write();
   	hIsomer_TDiff[kay]->Write();
      }
    }
    //JU stage 1 only
    esum  -> Write();
    esum2 -> Write();
    esum3 -> Write();
    esum4 -> Write();
    esum5 -> Write();
  tof_gated_QM -> Write();
  tof_gated_BM -> Write();
  tof_gated_QM_long -> Write();
  tof_gated_BM_long -> Write();




  } //End check for stage 1


  //Beam Monitors
  hU235_TOF->Write();  //Raw TOF for U235 Monitor
  hU235_TOF_Corr->Write(); //Corrected TOF for U235 Monitor
  hU235_PH_TOF->Write();
  hU235_PulseHeight->Write();  //Energy for U235 Monitor
  hU235_En->Write();  //Neutron Energy for U235 Monitor 
  hU235_En_Corr->Write();  //Neutron Energy for U235 Monitor (From Corrected TOF)
  hU235_Time_Between_Events->Write();
  
  hLi6_TOF->Write();  //Raw TOF for Li6 Monitor
  hLi6_TOF_Corr->Write(); //Corrected TOF for Li6 Monitor
  hLi6_PulseHeight->Write();  //Energy for Li6 Monitor
  hLi6_En->Write();  //Neutron Energy for Li6 Monitor 
  hLi6_En_Corr->Write();  //Neutron Energy for Li6 Monitor (From Corrected TOF)
  hLi6_PSD->Write(); 
  hLi6_Time_Between_Events->Write();

  hBkg_TOF->Write();  //Raw TOF for Bkg Monitor
  hBkg_TOF_Corr->Write(); //Corrected TOF for Bkg Monitor
  hBkg_PulseHeight->Write();  //Energy for Bkg Monitor
  hBkg_En->Write();  //Neutron Energy for Bkg Monitor 
  hBkg_En_Corr->Write();  //Neutron Energy for Bkg Monitor (From Corrected TOF)
  hBkg_Time_Between_Events->Write();


  hHe3_TOF->Write();  //Raw TOF for He3 Monitor
  hHe3_TOF_Corr->Write(); //Corrected TOF for He3 Monitor
  hHe3_PulseHeight->Write();  //Energy for He3 Monitor
  hHe3_En->Write();  //Neutron Energy for He3 Monitor 
  hHe3_En_Corr->Write();  //Neutron Energy for He3 Monitor (From Corrected TOF)
  hHe3_Time_Between_Events->Write();

  // JU histos ------
  hU235_TOF_gated-> Write();
  hU235_TOF_long_gated->Write();
  hLi6_TOF_gated -> Write();
  hLi6_TOF_long_gated->Write();
  hHe3_TOF_gated -> Write();
  hHe3_TOF_long_gated->Write();


  DANCE_Success("Analyzer","Wrote Histograms");
  return 0;
}


int Initialize_Analyzer(Input_Parameters input_params) {

  DANCE_Init("Analyzer","Initializing Analyzer");
  
  int func_ret = 0;
  
  totalindex = Read_TMatrix();
  
  if(totalindex <= 0) {
    func_ret = -1;
  }
  
  func_ret = Read_DMatrix();
  func_ret = Create_Analyzer_Histograms(input_params);
  
  for(int eye=0; eye<162; eye++) {
    last_timestamp_devent[eye]=0;
    last_energy_devent[eye]=0;

    last_devent_timestamp=0;
    last_valid_devent_timestamp=0;

  }

  //Diagnostics
  T0_Counter=0;

  DANCE_Entries_per_T0=0;
  DANCE_Events_per_T0=0;

  return func_ret;
}


int Analyze_Data(std::vector<DEVT_BANK> eventvector, Input_Parameters input_params, Analysis_Parameters *analysis_params) {


  int Crystal_Mult=0;
#ifdef Make_Removed_Spectra
  int Crystal_Mult2=0;
#endif

  //Initialize DANCE Event
  devent.Crystal_mult=0;
  devent.ESum=0;

  //No Physics events are valid until created
  devent.Valid=0;
  devent.pileup_detected=0;
  devent2.Valid=0;
  u235event.Valid=0;
  he3event.Valid=0;
  li6event.Valid=0;
  bkgevent.Valid=0;
 
  //Event length
  hEventLength->Fill(eventvector[eventvector.size()-1].timestamp-eventvector[0].timestamp);
  hEventLength_MCr->Fill(eventvector[eventvector.size()-1].timestamp-eventvector[0].timestamp,eventvector.size());

  //Loop over event 
  for(uint32_t eye=0; eye<eventvector.size(); eye++) {

    //Fill ID histogram
    int id_eye=eventvector[eye].ID;
        
    //this is a dance crystal
    if(id_eye<162) {
              
      //Coincidences
      if(eventvector.size() > 1 && eye < (eventvector.size()-1)) {
	
	//start with the next one
	for(uint32_t jay=eye+1; jay<eventvector.size(); jay++) {
	  
	  int id_jay = eventvector[jay].ID;
        
	  
	  if(id_jay<162 && id_jay>=0) {

	    //time difference between crystal jay and eye
	    double ddT = eventvector[jay].timestamp - eventvector[eye].timestamp;

	    //Fill the coincidence matrix
	    hCoinCAEN->Fill(id_eye,id_jay,1);
	    hCoinCAEN->Fill(id_jay,id_eye,1);
	    
	    //Look at time deviations
	    //Check relative to 0 first
	    if(id_eye==0) { 
	      hTimeDev_Rel0->Fill(-ddT,id_jay,1);
	    } 
	    if(id_jay==0) {
	      hTimeDev_Rel0->Fill(ddT,id_eye,1);
	    }
	    
	    //id11 is the index of where this crystal is in TMatrix if it is the first one listed
	    int id11=reftoindex1[eventvector[eye].ID];
	    //same but for if it is the second one listed
	    int id12=reftoindex2[eventvector[eye].ID];
	    
	    //same references to the crystals in the tmatrix as before but for the second crystal 
	    int id21=reftoindex1[eventvector[jay].ID];
	    //same but for if it is the second one listed
	    int id22=reftoindex2[eventvector[jay].ID];
	    
	    //if the leftmost crystal referernce is >=0 
	    //if the left and right are neighbors (reftocrystals are equal) 
	    //eye on left and jay on right
	    if(id11>=0 && id11==id22){
	      hTimeDev->Fill(-ddT,index1[id11]);
	    }
	    //if the leftmost crystal referernce is >=0 
	    //if the left and right are neighbors (reftocrystals are equal)
	    //eye on right and jay on left 
	    if(id12>=0 && id12==id21){
	      hTimeDev->Fill(ddT,index1[id21]);
	    }	  
	  }
	}
      }  //Done with coincidences 

    
      if(eventvector[eye].IsGamma == 1 && eventvector[eye].Valid==1) {

	//Make a DANCE Event
	devent.Crystal_ID[Crystal_Mult] = eventvector[eye].ID;        //Crystal ID
	devent.Cluster_ID[Crystal_Mult] = Crystal_Mult+1;             //Cluster ID
	devent.Islow[Crystal_Mult] = eventvector[eye].Islow;          //Crystal long integral
	devent.Ifast[Crystal_Mult] = eventvector[eye].Ifast;          //Crystal short integral
	devent.timestamp[Crystal_Mult] = eventvector[eye].timestamp;  //timestamp
	devent.tof[Crystal_Mult] = eventvector[eye].TOF;              //time of flight
	devent.tof_corr[Crystal_Mult] = eventvector[eye].TOF_Corr;    //corrected time of flight
	devent.Ecrystal[Crystal_Mult] = eventvector[eye].Eslow;       //Energy if calibrated 
	devent.ESum += eventvector[eye].Eslow; //ESum 
	devent.Crystal_mult++;
	devent.Valid=1;  //event is now valid
	Crystal_Mult++;	
	if(eventvector[eye].pileup_detected==1) {
	  //  cout<<"PU Det 1"<<endl;
	  devent.pileup_detected=1;
	}
	//	cout<<"Gamma  "<<devent.pileup_detected<<endl;
	
	DANCE_Entries_per_T0++;
	analysis_params->DANCE_entries_analyzed++;

      } //end of Valid Check
      
    }  
    
    //this is t0
    if(id_eye==200) {

      analysis_params->events_analyzed++;
      analysis_params->T0_events_analyzed++;
      analysis_params->T0_entries_analyzed++;

      if(input_params.IsomerSpectra) {
	for(int isom=0; isom<input_params.NIsomers; isom++) {
	  
	  for(int xx=0; xx<(int)Isomer_Prompt[isom].size(); xx++) {
	    for(int yy=0; yy<(int)Isomer_Delayed[isom].size(); yy++) {
	      
	      hIsomer_TDiff[isom]->Fill(Isomer_Delayed[isom][yy]-Isomer_Prompt[isom][xx]);	  
	      
	    } //end loop over delayed
	  } //end loop over prompt
	  
	  Isomer_Delayed[isom].clear();
	  Isomer_Prompt[isom].clear();
	  
	} //end loop over isomers
      } //end check on isomers
      
      //Do some T0 diagnostics
      if(analysis_params->last_timestamp[T0_ID] > 0) {
	hTimeBetweenT0s->Fill(eventvector[eye].timestamp-analysis_params->last_last_T0);
      }
      
      //Fill Histos
      hDANCE_Entries_per_T0->Fill(T0_Counter,DANCE_Entries_per_T0);
      hDANCE_Events_per_T0->Fill(T0_Counter,DANCE_Events_per_T0);
      
      //Reset Entries and Events per T0
      DANCE_Entries_per_T0=0;
      DANCE_Events_per_T0=0;
      
      //Incriment T0 counter
      T0_Counter++;
      
    }
   
    //beam monitors
    //He3
    if(id_eye == He3_ID) {
      analysis_params->He3_entries_analyzed++;

      he3event.tof = eventvector[eye].TOF;
      he3event.tof_corr = eventvector[eye].TOF_Corr;
      he3event.Ifast = eventvector[eye].Ifast;
      he3event.Islow = eventvector[eye].Islow;
      he3event.Valid = 1; //he3 event now valid

      //Fill time Diagnostics
      hHe3_Time_Between_Events->Fill(eventvector[eye].timestamp-analysis_params->last_timestamp[id_eye]);
    }
    
    //U235
    if(id_eye == U235_ID) {
      analysis_params->U235_entries_analyzed++;

      u235event.tof = eventvector[eye].TOF;
      u235event.tof_corr = eventvector[eye].TOF_Corr;
      u235event.Ifast = eventvector[eye].Ifast;
      u235event.Islow = eventvector[eye].Islow;
      u235event.Valid = 1; //u235 event now valid

      //Fill time Diagnostics
      hU235_Time_Between_Events->Fill(eventvector[eye].timestamp-analysis_params->last_timestamp[id_eye]);
    }
    
    //Li6
    if(id_eye == Li6_ID) {      
      analysis_params->Li6_entries_analyzed++;
      
      li6event.tof = eventvector[eye].TOF;
      li6event.tof_corr = eventvector[eye].TOF_Corr;
      li6event.Ifast = eventvector[eye].Ifast;
      li6event.Islow = eventvector[eye].Islow;
      li6event.Valid = 1; //li6 event now valid
      
      //Fill time Diagnostics
      hLi6_Time_Between_Events->Fill(eventvector[eye].timestamp-analysis_params->last_timestamp[id_eye]);
    }


    //Bkg
    if(id_eye == Bkg_ID) {      
      analysis_params->Bkg_entries_analyzed++;
      
      bkgevent.tof = eventvector[eye].TOF;
      bkgevent.tof_corr = eventvector[eye].TOF_Corr;
      bkgevent.Ifast = eventvector[eye].Ifast;
      bkgevent.Islow = eventvector[eye].Islow;
      bkgevent.Valid = 1; //bkg event now valid
      
      //Fill time Diagnostics
      hBkg_Time_Between_Events->Fill(eventvector[eye].timestamp-analysis_params->last_timestamp[id_eye]);
    }
    analysis_params->entries_analyzed++;
    
  } //end of eventvector loop


    //Handle various events and do some physics

    //DANCE Events
  if(devent.Valid == 1) {
    
    //   cout<<"processing"<<"  "<<devent.pileup_detected<<endl;
    analysis_params->events_analyzed++;
    analysis_params->DANCE_events_analyzed++;
    
    if(input_params.Analysis_Stage==1 || input_params.Read_Simulation ==1 ) {
      
      
      if(analysis_params->last_timestamp[T0_ID] > 0) {
      
	//   cout<<"DANCE Event: Mult: "<<devent.Crystal_mult<<"  TOF "<<devent.tof[0]<<endl;

	//DANCE Event Diagnostics
	DANCE_Events_per_T0++;

	//DEvent Blocking Time
	if((devent.timestamp[0]-last_valid_devent_timestamp) < input_params.DEvent_Blocking_Time) {
	  devent.Valid=0;
	  cout<<"DEvent Blocked"<<endl;
	}
	//Event not within blocking time
	else {
	  if((devent.timestamp[0]-last_valid_devent_timestamp) < input_params.Coincidence_Window) {
	    cout<<RED<<" Too small of a time difference! "<<devent.timestamp[0]-last_devent_timestamp<<RESET<<endl;
	    // string stuff;
	    // cin>>stuff;
	  }
	  hTimeBetweenDEvents->Fill(devent.timestamp[0]-last_valid_devent_timestamp);
#ifdef TurnOffGoFast	  
          hTimeBetweenDEvents_ESum_Mcr->Fill(devent.timestamp[0]-last_valid_devent_timestamp,devent.ESum,devent.Crystal_mult);
#endif
	  //Update the last valid devent timestamp.  Same non-paralyzable model
	  last_valid_devent_timestamp=devent.timestamp[0];
	}
     
	//Set Last devent timestamp to current one
	last_devent_timestamp=devent.timestamp[0];
      

	//Calculate energies 
	for(int kay=0; kay<devent.Crystal_mult; kay++) {
		
		
	  //Calculate the neutron energy    
	  devent.En[kay] = 0.5*939.565379e6*DANCE_FlightPath*DANCE_FlightPath/((devent.tof[kay])/1e9)/((devent.tof[kay])/1e9)/(2.997924589e8*2.997924589e8); 
	  //Calculate corrrected nuetron energy
	  devent.En_corr[kay] = 0.5*939.565379e6*DANCE_FlightPath*DANCE_FlightPath/((devent.tof_corr[kay])/1e9)/((devent.tof_corr[kay])/1e9)/(2.997924589e8*2.997924589e8); 
#ifdef TurnOffGoFast
	  hEn_Eg_Mcr->Fill(devent.En_corr[kay],devent.Ecrystal[kay],devent.Crystal_mult,1);
#endif
	  //Fill TOF for each crystal
	  hCrystal_En_Corr->Fill(devent.En_corr[kay],1);
	  hECrystal_En_Corr->Fill(devent.En_corr[kay],devent.Ecrystal[kay],1);

	  //Fill the 2D spectrum of TOF vs Crystal
	  hCrystalIDvsTOF->Fill(devent.tof[kay],devent.Crystal_ID[kay],1);
	  hCrystalTOF->Fill(devent.tof[kay],1);

#ifdef HighRateDebug
          if (devent.tof[kay]<92000 && devent.tof[kay]>100000 ) {
            hID_resonancegated->Fill(devent.Crystal_ID[kay]);
          }
          if (devent.tof[kay]>30000000 ) {
            hID_backgroundgated->Fill(devent.Crystal_ID[kay]);
          }
       //   En_ID->Fill(devent.En[kay],devent.Crystal_ID[kay],1);
          if (devent.Crystal_mult==2) {
            ISlow_ID_mcr2->Fill(devent.Islow[kay],devent.Crystal_ID[kay]);
          }
#endif
	
	  //Fill the 2D spectrum of Corrected TOF vs Crystal
	  hCrystalIDvsTOF_Corr->Fill(devent.tof_corr[kay],devent.Crystal_ID[kay],1);
	  hCrystalTOF_Corr->Fill(devent.tof_corr[kay],1);
	
	} //end loop over devent crystals
      } //end check of crystal mult and last_t0_timestamp
      else {
	devent.Valid = 0;
      }
      
      //Handle various events and do some physics
      if(devent.Valid==1) {
      
	//Fill TOF and En
	hEn->Fill(devent.En[0],1);
	hEn_Corr->Fill(devent.En_corr[0],1);
	hTOF->Fill(devent.tof[0]);
	hTOF_Corr->Fill(devent.tof_corr[0]);

	//Only one crystal...
	if(devent.Crystal_mult==1) {
	  devent.Cluster_mult=1;
	  devent.Ecluster[0]=devent.Ecrystal[0];
	}
      
	//Multiple Crystals... Clusterize
	if(devent.Crystal_mult > 1) {
	
	  // Here we will clusterize. Result will be in Cluster_ID[i] where i runs through all the crystals. The minimum of Cluster_ID is 1 so don't forget to lower it by 1. It is somewhat different then Jan Wouters routine but should lead to same results. Label value (after this routine) is the cluster multiplicity
	
	  int jj=0;
	  int Label=0;
	  int intcounter;
	  int Hits;
	  int Internal_ID[163];
	
	  for(int loop=0; loop<devent.Crystal_mult; loop++){
	    devent.Ecluster[loop]=0.;
	  
	    if(loop==0 || devent.Cluster_ID[loop]>Label){
	    
	      Label++;
	      intcounter=0;
	      Hits=1;
	      Internal_ID[0]=devent.Crystal_ID[loop];
	    
	      for(int l=0;l<Hits;l++) {
		for(jj=0;jj<devent.Crystal_mult;jj++) {
		  if(Internal_ID[l]!=devent.Crystal_ID[jj]) { // no need to check the same crystals
		  
		    bool IsNeighbor=false;
		  
		    if(DetMat1[Internal_ID[l]]==DetMat2[devent.Crystal_ID[jj]]) IsNeighbor=true;
		    else if(DetMat1[Internal_ID[l]]==DetMat3[devent.Crystal_ID[jj]]) IsNeighbor=true;
		    else if(DetMat1[Internal_ID[l]]==DetMat4[devent.Crystal_ID[jj]]) IsNeighbor=true;
		    else if(DetMat1[Internal_ID[l]]==DetMat5[devent.Crystal_ID[jj]]) IsNeighbor=true;
		    else if(DetMat1[Internal_ID[l]]==DetMat6[devent.Crystal_ID[jj]]) IsNeighbor=true;
		    else if(DetMat1[Internal_ID[l]]==DetMat7[devent.Crystal_ID[jj]]) IsNeighbor=true;
		  
		    if(IsNeighbor) {	
		      // if neighbor is found label it and add it to internal list (Internal_ID[]) for checking further in l loop (Hits is incremented by one)
		      // this way a new participant of the cluster is being chacked agains the other crystals
		    
		      if((1.*devent.Cluster_ID[jj])>=0 && (1.*devent.Cluster_ID[jj])>(1.*Label)) {
			// here if the crystal is already labeled we skip so that we do not have to repeat already labeled crystals more speed to it
		      
			devent.Cluster_ID[jj]=Label; // this one is important- Adding the crystals to the cluster
		      
			Hits++;
			intcounter++;
			Internal_ID[intcounter]=devent.Crystal_ID[jj];	
		      }
		    }
		  }
		}
	      }
	      if(Hits==1) devent.Cluster_ID[loop]=Label;
	    }
	  }
	
	  devent.Cluster_mult=Label;
	
	  // Fill the cluster energy. Labels go through 1,2,3 ... CrystalMult
	  for(int ii=0;ii<devent.Crystal_mult;ii++){
	    devent.Ecluster[devent.Cluster_ID[ii]-1]+=devent.Ecrystal[ii];
	  };
	}  //Done checking mult > 1
            
	//Mults vs TOF
	hTOF_Mcl->Fill(devent.tof[0],devent.Cluster_mult);
	hTOF_Mcl_Corr->Fill(devent.tof_corr[0],devent.Cluster_mult);

	//Mults vs ESum vs En
	En_Esum_Mcl->Fill(devent.En_corr[0],devent.ESum,devent.Cluster_mult);
	En_Esum_Mcr->Fill(devent.En_corr[0],devent.ESum,devent.Crystal_mult);

	if(devent.pileup_detected==1) {
	  // cout<<"PU Det"<<endl;
	  En_Esum_Mcr_Pileup->Fill(devent.En_corr[0],devent.ESum,devent.Crystal_mult);
	}
	else {
	  En_Esum_Mcr_NoPileup->Fill(devent.En_corr[0],devent.ESum,devent.Crystal_mult);
	}
#ifdef TurnOffGoFast      
	if(devent.Crystal_mult == 1) {
	  En_Ecr1_mcr1->Fill(devent.En_corr[0],devent.Ecrystal[0]);
	}
	if(devent.Crystal_mult == 2) {
	  En_Ecr1_Ecr2_mcr2->Fill(devent.En_corr[0],devent.Ecrystal[0],devent.Ecrystal[1]);
	}
#endif
	//Mults vs ESum vs TOF
	hTOF_Esum_Mcl->Fill(devent.tof_corr[0],devent.ESum,devent.Cluster_mult);
	hTOF_Esum_Mcr->Fill(devent.tof_corr[0],devent.ESum,devent.Crystal_mult);

	//Do Isomer Stuff
	if(input_params.IsomerSpectra) {
	  //Loop over number of isomers
	  for(int isom=0; isom<input_params.NIsomers; isom++) {
	  
	    //prompt
	    if(devent.tof_corr[0] >= input_params.IsomerPromptTOFGates[isom*2] && devent.tof_corr[0] < input_params.IsomerPromptTOFGates[isom*2+1]) {
	      if(devent.Cluster_mult >= input_params.IsomerPromptMclGates[isom*2] && devent.Cluster_mult < input_params.IsomerPromptMclGates[isom*2+1]) {
		if(devent.ESum >= input_params.IsomerPromptQGates[isom*2] && devent.ESum < input_params.IsomerPromptQGates[isom*2+1]) {
		
		  //Push time onto vector
		  Isomer_Prompt[isom].push_back(devent.tof_corr[0]);
		  //Fill Singles Histogram
		  hIsomer_Prompt[isom]->Fill(devent.tof_corr[0]);
		
		} //End check on TOF
	      } //End check on Mcl
	    } //End check on ESum

	    //prompt
	    if(devent.tof_corr[0] >= input_params.IsomerDelayedTOFGates[isom*2] && devent.tof_corr[0] < input_params.IsomerDelayedTOFGates[isom*2+1]) {
	      if(devent.Cluster_mult >= input_params.IsomerDelayedMclGates[isom*2] && devent.Cluster_mult < input_params.IsomerDelayedMclGates[isom*2+1]) {
		if(devent.ESum >= input_params.IsomerDelayedQGates[isom*2] && devent.ESum < input_params.IsomerDelayedQGates[isom*2+1]) {
		
		  //Push time onto vector
		  Isomer_Delayed[isom].push_back(devent.tof_corr[0]);
		  //Fill Singles Histogram
		  hIsomer_Delayed[isom]->Fill(devent.tof_corr[0]);
		
		} //End check on TOF
	      } //End check on Mcl
	    } //End check on ESum
	  } //End loop over Isomers
	} //End check on isomer spectra

	// fill esum histos (for convenience, they are redundant to the 3d histo above --- JU -----)
	esum -> Fill(devent.ESum);
	if(devent.Cluster_mult==2)esum2->Fill(devent.ESum);
	if(devent.Cluster_mult==3)esum3->Fill(devent.ESum);
	if(devent.Cluster_mult==4)esum4->Fill(devent.ESum);
	if(devent.Cluster_mult==5)esum5->Fill(devent.ESum);

	if(input_params.QGatedSpectra) {
	
	  //Loop over the cluster mult
	  for(int jay=0; jay<devent.Cluster_mult; jay++ )  {
	    for(int kay=0; kay<input_params.NQGates; kay++) {
	      if(devent.ESum > input_params.QGates[kay*2] && devent.ESum < input_params.QGates[kay*2+1]) {
		En_Ecl_Mcl_QGated[kay]-> Fill(devent.En_corr[0], devent.Ecluster[jay], devent.Cluster_mult );
		hTOF_Mcl_QGated[kay]->Fill(devent.tof_corr[0],devent.Cluster_mult);
	      } //Done Checking QGates
	    } //Done Looping over QGates
	  } //Done looping over cluster mult
	

	  
	  //Loop over the crystal mult
	  for(int jay=0; jay<devent.Crystal_mult; jay++ )  {
	    for(int kay=0; kay<input_params.NQGates; kay++) {
	      if(devent.ESum >input_params.QGates[kay*2] && devent.ESum < input_params.QGates[kay*2+1]) {
		En_Ecr_Mcr_QGated[kay]-> Fill(devent.En_corr[0], devent.Ecrystal[jay], devent.Crystal_mult );
		ID_Ecr_Mcr_QGated[kay]-> Fill(devent.Crystal_ID[jay], devent.Ecrystal[jay], devent.Crystal_mult );
	      } //Done Checking QGates
	    } //Done Looping over QGates
	  } //Done looping over crystal mult
	}

	//Crystal mult 1
	if(devent.Crystal_mult == 1) {
	  hGamma_Mcr1->Fill(devent.Islow[0], devent.Crystal_ID[0],1);
	  hGammaCalib_Mcr1->Fill(devent.Ecrystal[0],devent.Crystal_ID[0],1);
          if (devent.tof_corr[0] > 14000000){
            hGammaCalib_Mcr1_late->Fill(devent.Ecrystal[0],devent.Crystal_ID[0],1);
          }
	}
        if(devent.Crystal_mult == 2) {
          if (devent.tof_corr[0] > 14000000){
            hGammaCalib_Mcr1_late->Fill(devent.Ecrystal[0],devent.Crystal_ID[0],1);
          }
	}
	
	double largesttimediff=0;

	//Loop over the crystal mult
	if(devent.Crystal_mult>1) {
	  for(int jay=0; jay<devent.Crystal_mult; jay++ )  {
	    
	    if(jay > 0) {
	      if(devent.timestamp[jay]-devent.timestamp[jay-1] > largesttimediff) {
		largesttimediff = devent.timestamp[jay]-devent.timestamp[jay-1];	    
	      }
	    }
	  }
	}
	hEventTimeDist_Etot->Fill(largesttimediff,devent.ESum);

	
	/*
      
	// At this point in the code, we apparently have TOF, Esum, and Mcluster
	// so fill the "gated TOF histograms" - these are higher resolution TOF spectra than hTOF_Esum_Mcl
	if(devent.Cluster_mult>1 && devent.Cluster_mult<7) {	// Edit these limits for now
	  if(devent.ESum > input_params.QGates[0] && devent.ESum < input_params.QGates[1]) {	// " Q-value" gate    ---- JU -------
	    tof_gated_QM -> Fill(devent.tof_corr[0]);
	    tof_gated_QM_long -> Fill(devent.tof_corr[0]);
	  }
	  if(devent.ESum >input_params.QGates[4] && devent.ESum < input_params.QGates[5]) {	// "Background" Q gate
	    tof_gated_BM -> Fill(devent.tof_corr[0]);
	    tof_gated_BM_long -> Fill(devent.tof_corr[0]);
	  }
	}
	*/

	hEventLength_Etot->Fill(eventvector[eventvector.size()-1].TOF-eventvector[0].TOF,devent.ESum);


  
      
#ifdef Make_Removed_Spectra

	//Clear the removed crystals
	removed_crystals.clear();
#ifdef Removed_Verbose
	cout<<endl<<endl;
#endif
	//make the removed spectra
	for(int kay=1; kay<Max_Gamma_Removed; kay++) {

	  //Initialize Second DANCE Event
	  devent2.Crystal_mult=0;
	  devent2.ESum=0;
	
	  //No Physics events are valid until created
	  devent2.Valid=0;
	  Crystal_Mult2=0;
	  	
	  //If it is greater than kay do stuff  
	  if(devent.Crystal_mult>kay) {
#ifdef Removed_Verbose
	    cout<<"kay: "<<kay<<endl;
#endif
	    // //loop over number of removed gammas
	    // for(int el=0; el<kay; el++) {
	    
	    //The first one just gets pushed back
	    if(removed_crystals.size()==0) {
	      int removed_hit = gRandom->Uniform(0,devent.Crystal_mult); //Pick a random crystal
#ifdef Removed_Verbose
	      cout<<"removed hit: "<<removed_hit<<"  vector size: "<<removed_crystals.size()<<endl;
#endif
	      removed_crystals.push_back(removed_hit);
	    }
	      
	    else {
	      //set duplicate false to start
	      //put it in a while loop to remove duplicates
	      while(true) {
		bool duplicate=false;

		//Pick a crystal at random
		int removed_hit = gRandom->Uniform(0,devent.Crystal_mult); //Pick a random crystal
#ifdef Removed_Verbose
		cout<<"removed hit: "<<removed_hit<<"  vector size: "<<removed_crystals.size()<<endl;
#endif
		//Check to see if it is a duplicate
		for(int em=0; em<(int)removed_crystals.size(); em++) {
		  if(removed_hit == removed_crystals[em]) {
		    //duplicate
		    duplicate = true;
		  }
		}
		//Not a duplicate
		if(!duplicate) {
		  removed_crystals.push_back(removed_hit);
		  break;
		} //end check on !duplicate
	      } //end while true
	    } //end check on not first removed
	    //  }
#ifdef Removed_Verbose
	    cout<<"Original: mult "<<devent.Crystal_mult<<"  ";
	    for(int eye=0; eye < devent.Crystal_mult; eye++) {
	      cout<<devent.Crystal_ID[eye]<<"  ";
	    }
	    cout<<endl;
#endif

	    //Throw away last crystal
	    //Loop over event 
	    for(int eye=0; eye < devent.Crystal_mult; eye++) {
	      
	      bool keep_this = true;
	      //Look to see if this hit was removed 
	      for(int jay=0; jay<(int)removed_crystals.size(); jay++) {
		if(eye == removed_crystals[jay]) {
		  keep_this = false;
		}
	      }
	      if(keep_this) {
		
		//Make a DANCE Event 2
		devent2.Crystal_ID[Crystal_Mult2] = devent.Crystal_ID[eye];  //Crystal ID
		devent2.Cluster_ID[Crystal_Mult2] = devent.Cluster_ID[eye];  //??????
		devent2.Ecrystal[Crystal_Mult2] = devent.Ecrystal[eye];   //Energy if calibrated 
		devent2.ESum += devent.Ecrystal[eye];
		devent2.tof[Crystal_Mult2] = devent.tof[Crystal_Mult2];
		devent2.tof_corr[Crystal_Mult2] = devent.tof_corr[Crystal_Mult2];
		devent2.En[Crystal_Mult2] = devent.En[Crystal_Mult2];
		devent2.En_corr[Crystal_Mult2] = devent.En_corr[Crystal_Mult2];   
		devent2.Crystal_mult++;
		devent2.Valid=1;  //event is now valid
		Crystal_Mult2++;	
	      } //end of check on keep_this
	    } //end of loop over original dance crystal mult
	  } //end check on original dance crystal mult > #removed
	  
	
	  if(devent2.Valid == 1 ) {
	    hTOF_Esum_Mcr_Removed[kay]->Fill(devent2.tof_corr[0],devent2.ESum,devent2.Crystal_mult);
	    hEn_Esum_Mcr_Removed[kay]->Fill(devent2.En_corr[0],devent2.ESum,devent2.Crystal_mult);
#ifdef Removed_Verbose
	    cout<<"Removed: mult "<<devent2.Crystal_mult<<"  ";
	    for(int eye=0; eye < devent2.Crystal_mult; eye++) {
	      cout<<devent2.Crystal_ID[eye]<<"  ";
	    }
	    cout<<endl;
#endif
	  }
	
	}     

#endif
      } //end check valid dance event

    } //End of stage and simulation



  } //Done checking for stage 1
  
    //U235 Beam Monitor Events
  if(u235event.Valid == 1) {

    analysis_params->events_analyzed++;
    analysis_params->U235_events_analyzed++;

    if(analysis_params->last_timestamp[T0_ID] > 0) {
      //Fill Pulse Height Spectrum
      hU235_PulseHeight->Fill(u235event.Islow);
     
      //Fill TOF spectra
      hU235_TOF->Fill(u235event.tof);
      hU235_TOF_Corr->Fill(u235event.tof_corr);
      hU235_PH_TOF->Fill(u235event.tof,u235event.Islow);
      
      // JU Histograms ----
      //if(u235event.Islow>5000.0 && u235event.Islow<35000.0)
      if(u235event.Islow>2000.0 && u235event.Islow<20000.0)
	{
	  hU235_TOF_gated -> Fill(u235event.tof);
	  hU235_TOF_long_gated -> Fill(u235event.tof);
	}
      //-------------------
      
    } //end check over last t0 timestamp
    else {
      u235event.Valid=0;
    }
    //Handle various events and do some physics
    if(u235event.Valid==1) {
      //Calculate the neutron energy    
      u235event.En = 0.5*939.565379e6*U235_FlightPath*U235_FlightPath/((u235event.tof)/1e9)/((u235event.tof)/1e9)/(2.997924589e8*2.997924589e8); 
      //Calculate corrrected nuetron energy
      u235event.En_corr = 0.5*939.565379e6*U235_FlightPath*U235_FlightPath/((u235event.tof_corr)/1e9)/((u235event.tof_corr)/1e9)/(2.997924589e8*2.997924589e8); 
      hU235_En->Fill(u235event.En);
      hU235_En_Corr->Fill(u235event.En_corr);
    } 
  }
  
  //He3 Beam Monitor Events
  if(he3event.Valid == 1) {

    analysis_params->events_analyzed++;
    analysis_params->He3_events_analyzed++;

    if(analysis_params->last_timestamp[T0_ID] > 0) {
      //Fill Pulse Height Spectrum
      hHe3_PulseHeight->Fill(he3event.Islow);
   
      //Fill TOF spectra
      hHe3_TOF->Fill(he3event.tof);
      hHe3_TOF_Corr->Fill(he3event.tof_corr);

      // JU Histograms ----
      if(he3event.Islow>0.0 && he3event.Islow<35000.0)
	{
	  hHe3_TOF_gated -> Fill(he3event.tof);
	  hHe3_TOF_long_gated -> Fill(he3event.tof);
	}
      //-------------------

    } //end check over last t0 timestamp
    else {
      he3event.Valid=0;
    }
    //Handle various events and do some physics
    if(he3event.Valid==1) {
      //Calculate the neutron energy    
      he3event.En = 0.5*939.565379e6*He3_FlightPath*He3_FlightPath/((he3event.tof)/1e9)/((he3event.tof)/1e9)/(2.997924589e8*2.997924589e8); 
      //Calculate corrrected nuetron energy
      he3event.En_corr = 0.5*939.565379e6*He3_FlightPath*He3_FlightPath/((he3event.tof_corr)/1e9)/((he3event.tof_corr)/1e9)/(2.997924589e8*2.997924589e8); 
      hHe3_En->Fill(he3event.En);
      hHe3_En_Corr->Fill(he3event.En_corr);
    } 
  }

  //Li6 Beam Monitor Events
  if(li6event.Valid == 1) {
    
    analysis_params->events_analyzed++;
    analysis_params->Li6_events_analyzed++;

    if(analysis_params->last_timestamp[T0_ID] > 0) {
      //Fill Pulse Height Spectrum
      hLi6_PulseHeight->Fill(li6event.Islow);
      hLi6_PSD->Fill(li6event.Islow,li6event.Ifast);
   
      //Fill TOF spectra
      hLi6_TOF->Fill(li6event.tof);
      hLi6_TOF_Corr->Fill(li6event.tof_corr);

      // JU Histograms ----
      //if(li6event.Islow>17000.0 && li6event.Islow<24000.0)
      if(li6event.Islow>17000.0 && li6event.Islow<24000.0)	// Runs 10194 - 
	{
	  hLi6_TOF_gated -> Fill(li6event.tof);
	  hLi6_TOF_long_gated -> Fill(li6event.tof);
	}
      //-------------------
    } //end check over last t0 timestamp
    else {
      li6event.Valid=0;
    }
    //Handle various events and do some physics
    if(li6event.Valid==1) {
      //Calculate the neutron energy    
      li6event.En = 0.5*939.565379e6*Li6_FlightPath*Li6_FlightPath/((li6event.tof)/1e9)/((li6event.tof)/1e9)/(2.997924589e8*2.997924589e8); 
      //Calculate corrrected nuetron energy
      li6event.En_corr = 0.5*939.565379e6*Li6_FlightPath*Li6_FlightPath/((li6event.tof_corr)/1e9)/((li6event.tof_corr)/1e9)/(2.997924589e8*2.997924589e8); 
      hLi6_En->Fill(li6event.En);
      hLi6_En_Corr->Fill(li6event.En_corr);
    } 
  }

  //Bkg Beam Monitor Events
  if(bkgevent.Valid == 1) {
    
    analysis_params->events_analyzed++;
    analysis_params->Bkg_events_analyzed++;

    if(analysis_params->last_timestamp[T0_ID] > 0) {
      //Fill Pulse Height Spectrum
      hBkg_PulseHeight->Fill(bkgevent.Islow);
   
      //Fill TOF spectra
      hBkg_TOF->Fill(bkgevent.tof);
      hBkg_TOF_Corr->Fill(bkgevent.tof_corr);

    } //end check over last t0 timestamp
    else {
      bkgevent.Valid=0;
    }
    //Handle various events and do some physics
    if(bkgevent.Valid==1) {
      //Calculate the neutron energy    
      bkgevent.En = 0.5*939.565379e6*Bkg_FlightPath*Bkg_FlightPath/((bkgevent.tof)/1e9)/((bkgevent.tof)/1e9)/(2.997924589e8*2.997924589e8); 
      //Calculate corrrected nuetron energy
      bkgevent.En_corr = 0.5*939.565379e6*Bkg_FlightPath*Bkg_FlightPath/((bkgevent.tof_corr)/1e9)/((bkgevent.tof_corr)/1e9)/(2.997924589e8*2.997924589e8); 
      hBkg_En->Fill(bkgevent.En);
      hBkg_En_Corr->Fill(bkgevent.En_corr);
    } 
  }


  //Done with Event Processing
  return 0;

}
