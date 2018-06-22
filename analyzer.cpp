//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  analyzer.cpp           *// 
//*  Last Edit: 05/08/18    *//  
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
#include "TCutG.h"
#include "TGraph.h"

using namespace std;

//#define Histogram_DetectorLoad

//Analyzer stats
uint32_t events_analyzed;
uint32_t entries_analyzed;
uint32_t analyzer_counter;

//output binary file
ofstream outputbinfile;

//Structure to write data
DEVT_STAGE1 devt_out;

//Graphs of TOF Corrections
TGraph *gr_DANCE_TOF_Corr;
TGraph *gr_U235_TOF_Corr;
TGraph *gr_Li6_TOF_Corr;
TGraph *gr_He3_TOF_Corr;

//Limits of TOF Corrections
double DANCE_TOF_Corr_Limit[2];   //[0] is lower [1] is upper
double U235_TOF_Corr_Limit[2];   //[0] is lower [1] is upper
double Li6_TOF_Corr_Limit[2];   //[0] is lower [1] is upper
double He3_TOF_Corr_Limit[2];   //[0] is lower [1] is upper

/* HISTOGRAMS */
TH1D *hID;  //ID's present in datastream
TH2D *hCoinCAEN;  //coincidence matrix

//Time Diagnostics
TH1D *hEventLength;  //length in ns of the event (diagnostic)
TH2D *hTimeBetweenCrystals;  //time between subsequent hits of the same crystal (ns)
TH1D *hTimeBetweenDEvents;  //time between subsequent DANCE events (ns)
TH3D *hTimeBetweenDEvents_ESum_Mcr; //DANCE ESum vs time between subsequent DANCE events for Mcr==1
TH1D *hTimeBetweenT0s; //time between T0s (ns)
TH2D *hCrystalIDvsTOF; //TOF for each crystal 
TH1D *hCrystalTOF; //TOF for all crytals
TH2D *hCrystalIDvsTOF_Corr; //TOF for each crystal 
TH1D *hCrystalTOF_Corr; //TOF for all crytals
TH1D *hDetectorLoad; //Average detector load vs TOF

TH1D *hDANCE_Entries_per_T0;
TH1D *hDANCE_Events_per_T0;

//Time Deviations
TH2D *hTimeDev_Rel0;  //Time deviations of all crystals to crystal 0
TH2D *hTimeDev;  //Time deviations relative to adjacent crystals

//Energy Histograms
TH2D *ADC_calib;  //2D PSD Plot (calibrated)
TH2D *ADC_raw;    //2D PSD Plot (uncalibrated)
TH2D *ADC_calib_Invalid;  //2D PSD Plot of invalid events (calibrated)

TH3D *ADC_raw_ID;    //3D Plot of 2D PSD plot (uncalibrated) vs DANCE Crystal ID (z)
TH3D *ADC_calib_ID;  //3D Plot of 2D PSD plot (calibrated) vs DANCE Crystal ID (z)

//Physics Spectra
TH1D *hEn; //Neutron Energy from dance events
TH1D *hTOF; //TOF for dance events
TH1D *hEn_Corr; //Neutron Energy from dance events
TH1D *hTOF_Corr; //TOF for dance events

//Alpha Spectra
TH2D *hAlpha;
TH2D *hAlphaCalib;

//Gamma Spectra
TH2D *hGamma;
TH2D *hGammaCalib;

//QGated 3D Histograms;

//Flag to turn on the QGated spectra
bool QGated_Spectra = true;

//This sets the number of QGates
//const int NQGates = 3;

//This defines the Qgates for the En_Ecl_Mcl QGated spectra
//double QGates[2*NQGates] = {6.5,7.0,   //QGate 1
//			    6.5,7.5,   //QGate 2
//			    6.5,6.7};  //QGate 3
//

TH3F *En_Ecl_Mcl_QGated[10]; //Max is 10 QGates. 


//3D Histograms
TH3F *En_Esum_Mcl;
TH3F *En_Esum_Mcr;

TH3F *hTOF_Esum_Mcl;
TH3F *hTOF_Esum_Mcr;

TH3F *En_Eg_Mcl; // this should have a Qgate on it
TH3F *En_Eg_Mcr; // this should have a Qgate on it

TH3F *Esum_Eg_Mcl; // Eg is Ecluster here
TH3F *Esum_Eg_Mcr; // Eg is Ecrystal here

//Beam Monitors
TH1D *hU235_TOF;  //Raw TOF for U235 Monitor
TH1D *hU235_TOF_Corr; //Corrected TOF for U235 Monitor
TH2D *hU235_PH_TOF;
TH1D *hU235_PulseHeight;  //Energy for U235 Monitor
TH1D *hU235_En;  //Neutron Energy for U235 Monitor 
TH1D *hU235_En_Corr;  //Neutron Energy for U235 Monitor (From Corrected TOF)
TH1D *hU235_Time_Between_Events; //Time between U235 hits

TH1D *hLi6_TOF;  //Raw TOF for Li6 Monitor
TH1D *hLi6_TOF_Corr; //Corrected TOF for Li6 Monitor
TH1D *hLi6_PulseHeight;  //Energy for Li6 Monitor
TH1D *hLi6_En;  //Neutron Energy for Li6 Monitor 
TH1D *hLi6_En_Corr;  //Neutron Energy for Li6 Monitor (From Corrected TOF)
TH2D *hLi6_PSD; 
TH1D *hLi6_Time_Between_Events; //Time between U235 hits

TH1D *hHe3_TOF;  //Raw TOF for He3 Monitor
TH1D *hHe3_TOF_Corr; //Corrected TOF for He3 Monitor
TH1D *hHe3_PulseHeight;  //Energy for He3 Monitor
TH1D *hHe3_En;  //Neutron Energy for He3 Monitor 
TH1D *hHe3_En_Corr;  //Neutron Energy for He3 Monitor (From Corrected TOF)
TH1D *hHe3_Time_Between_Events; //Time between U235 hits

/* CUTS */
TCutG *Gamma_Gate;
TCutG *Alpha_Gate;

//Dance Event
DANCE_Event devent;
U235_Event u235event;
He3_Event he3event;
Li6_Event li6event;

/* VARIABLES */

double last_timestamp[256];        //This keeps track of the last timestamp valid or not
double last_valid_timestamp[256];  //This keeps track of the last timestamp that was valid
double current_timestamp[256];     //current timestamp 

double last_t0_timestamp;
double current_t0_timestamp;

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

//energy calibrations
double slow_offset[200];
double slow_slope[200];
double slow_quad[200];
double fast_slope[200];
double fast_offset[200];

//Diagnostics
uint32_t T0_Counter=0;
double Detector_Load[100000000];
uint32_t DANCE_Entries_per_T0=0;
uint32_t DANCE_Events_per_T0=0;

//This function goes through and makes the run-by-run time deviations
int Make_Time_Deviations(int RunNumber) {
  
  cout<<"Analyzer [INFO]: Making Time Deviations"<<endl;
 
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
    hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-100.0, hTimeDev->GetMean()+100.0);
    hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-50.0, hTimeDev->GetMean()+50.0);
    hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-10.0, hTimeDev->GetMean()+10.0);
    hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-5.0, hTimeDev->GetMean()+5.0);
    hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-5.0, hTimeDev->GetMean()+5.0);
    hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-5.0, hTimeDev->GetMean()+5.0);
    hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-4.0, hTimeDev->GetMean()+4.0);
    hTimeDev->GetXaxis()->SetRangeUser(hTimeDev->GetMean()-4.0, hTimeDev->GetMean()+4.0);
    // cout<<hTimeDev->GetMean()<<"  "<<time_deviation<<endl;
    time_deviation += hTimeDev->GetMean();
    td_out <<index2[eye]<<"   \t"<<time_deviation<<"\n";
  }
  
  cout<<GREEN<<"Analyzer [INFO]: Made Time Deviations for Run: "<<RunNumber<<RESET<<endl;
  return 0;
}




//This function fetches the TOF Correction plots for the moderation time 
int Read_Moderation_Time_Graphs() {
 
  cout<<"Analyzer [INFO]: Reading TOF Corrections"<<endl;

  ifstream tof_corr;
  tof_corr.open("./TOF_Corrections/TOF_Corrections.txt");

  if(tof_corr.is_open()) {
    //Find out how many points there are 
    int npoints;
    tof_corr>>npoints;
    const int N=npoints;
    
    //Time of flight from neutron energy
    double DANCE_TOF[N];
    double U235_TOF[N];
    double He3_TOF[N];
    double Li6_TOF[N];

    //Time of flight plus moderation time
    double DANCE_TOF_Measured[N];
    double U235_TOF_Measured[N];
    double He3_TOF_Measured[N];
    double Li6_TOF_Measured[N];

    for(int eye=0; eye<N; eye++) {
      tof_corr>>DANCE_TOF[eye]>>DANCE_TOF_Measured[eye]>>U235_TOF[eye]>>U235_TOF_Measured[eye]>>Li6_TOF[eye]>>Li6_TOF_Measured[eye]>>He3_TOF[eye]>>He3_TOF_Measured[eye];
    }
    
    //Graphs of TOF Corrections
    gr_DANCE_TOF_Corr = new TGraph(N,DANCE_TOF_Measured,DANCE_TOF);
    gr_U235_TOF_Corr = new TGraph(N,U235_TOF_Measured,U235_TOF);
    gr_He3_TOF_Corr = new TGraph(N,He3_TOF_Measured,He3_TOF);
    gr_Li6_TOF_Corr = new TGraph(N,Li6_TOF_Measured,Li6_TOF);
    
    DANCE_TOF_Corr_Limit[0] = DANCE_TOF_Measured[N-1];
    DANCE_TOF_Corr_Limit[1] = DANCE_TOF_Measured[0];
    U235_TOF_Corr_Limit[0] = U235_TOF_Measured[N-1];
    U235_TOF_Corr_Limit[1] = U235_TOF_Measured[0];
    Li6_TOF_Corr_Limit[0] = Li6_TOF_Measured[N-1];
    Li6_TOF_Corr_Limit[1] = Li6_TOF_Measured[0];
    He3_TOF_Corr_Limit[0] = He3_TOF_Measured[N-1];
    He3_TOF_Corr_Limit[1] = He3_TOF_Measured[0];
       
    cout<<GREEN<<"Analyzer [INFO]: Read In TOF Corrections"<<RESET<<endl;
    return 0;
  }
  else {
    cout<<RED<<"Analyzer [ERROR]: Faild to Read In TOF Corrections"<<RESET<<endl;
    return -1;
  }  
}



//This function makes the output binary file for stage 0 or 1 filled with time-ordered devt_bank structures
int Make_Output_Binfile(int RunNumber, bool read_binary) {
  
  stringstream outfilename;
  outfilename.str();
  
  //stage0 
  if(read_binary==0) {
    outfilename << STAGE0_BIN; 
    outfilename <<"/stage0_run_";
  }
  //stage1
  if(read_binary==1) {
    outfilename << STAGE1_BIN;
    outfilename <<"/stage0_run_";
  }
  outfilename << RunNumber << ".bin";
  
  outputbinfile.open(outfilename.str().c_str(), ios::out | ios::binary);
  
  if(outputbinfile.is_open()) {
    cout<<GREEN<<"Analyzer [INFO]: Succesfully created and opened output binary file: "<<outfilename.str()<<RESET<<endl;
    return 0;
  }
  else {
    cout<<RED<<"Analyzer [ERROR]: Failed to create output binary file: "<<outfilename.str()<<RESET<<endl;
    return -1;
  }
}

int Read_Energy_Calibrations(int RunNumber, bool read_binary, bool read_simulation) {
  
  int id=0;
  double temp[5] = {0,0,0,0,0};

  for(int eye=0; eye<200; eye++) {
    slow_offset[eye]=0;
    slow_slope[eye]=1;
    slow_quad[eye]=0;
    fast_offset[eye]=0;
    fast_slope[eye]=1;  
  }

  cout<<"Analyzer [INFO]: Reading Energy Calibrations"<<endl;

  if(read_simulation == 0) {


    stringstream cal_name;
    cal_name.str();
    cal_name << CALIB_DIR << "/param_out_" << RunNumber << ".txt";
    
    ifstream encal;
    encal.open(cal_name.str().c_str());
  
    bool encalfail=false;
  
    if(read_binary==1) {
      if(encal.is_open()) {
	while(!encal.eof()) {
	  encal>>id>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4];
	  slow_offset[id]=temp[0];
	  slow_slope[id]=temp[1];
	  slow_quad[id]=temp[2];
	  fast_offset[id]=temp[3];
	  fast_slope[id]=temp[4];
	
	  if(encal.eof()) {
	    break;
	  }
	}  
      
	cout<<GREEN<<"Analyzer [INFO]: Read Energy Calibrations File: "<<cal_name.str()<<RESET<<endl;
	return 0;
      }
      else {
	cout<<RED<<"Analyzer [ERROR]: File: "<<cal_name.str()<<" NOT found..."<<RESET<<endl;
	encalfail=true;
      }
    }
    if(read_binary==0 || encalfail==true)
      cout<<"Analyzer [INFO]: Looking for calib_ideal.dat"<<endl;
  
    stringstream idealcal_name;
    idealcal_name.str();
    idealcal_name << CALIB_DIR << "/calib_ideal.dat";
  
    ifstream idealcal;
    idealcal.open(idealcal_name.str().c_str());
  
    if(idealcal.is_open()) {
      while(!idealcal.eof()) {
	idealcal>>id>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4];
	slow_offset[id]=temp[0];
	slow_slope[id]=temp[1];
	slow_quad[id]=temp[2];
	fast_offset[id]=temp[3];
	fast_slope[id]=temp[4];

	  cout<<id<<"  "<<slow_slope[id]<<"  "<<fast_slope[id]<<endl;
      
	if(idealcal.eof()) {
	  break;
	}
      }
      cout<<GREEN<<"Analyzer [INFO]: Opened "<<idealcal_name.str()<<RESET<<endl;
      return 0;
    }
    else {
      cout<<"Analyzer [ERROR]: Can NOT Open Energy Calibration File..."<<RESET<<endl;
      return-1;
    }
  }  
  cout<<GREEN<<"Analyzer [INFO]: Set all calibrations off for simulations"<<RESET<<endl;
  return 0;

}

int Read_PI_Gates() {

  char name[200];
  int Nalphacut=0;
  int Ngammacut=0;
  double x_alphacut[200];
  double y_alphacut[200];
  double x_gammacut[200];
  double y_gammacut[200];

  sprintf(name,"Gates/%s",ALPHAGATE);
  ifstream alphacutin(name);
  if(alphacutin.is_open()) {
    cout << "Analyzer [INFO]: Reading AlphaCut " << name << endl;
    while(!alphacutin.eof()){
      alphacutin  >> x_alphacut[Nalphacut] >> y_alphacut[Nalphacut];
      Nalphacut++;
      if(alphacutin.eof()) break;
    }
    // cout<<Nalphacut<<endl;
    alphacutin.close();
    cout<<GREEN<<"Analyzer [INFO]: Alpha PI cut loaded" <<RESET<<endl;
  }
  else {
    cout<<RED<<"Analyzer [ERROR]: Alpha PI cut "<<name<<" file NOT found." <<RESET<<endl;
  }
  
  
  sprintf(name,"Gates/%s",GAMMAGATE);
  ifstream gammacutin(name);
  if(gammacutin.is_open()) {
    cout<<"Analyzer [INFO]: Reading GammaCut " <<name<<endl;
    while(!gammacutin.eof()){
      gammacutin  >> x_gammacut[Ngammacut] >> y_gammacut[Ngammacut];
      Ngammacut++;
      if(gammacutin.eof()) break;
    }
    // cout<<Ngammacut<<endl;
    gammacutin.close();
    cout<<GREEN<<"Analyzer [INFO]: Gamma PI cut loaded" <<RESET<<endl;
  }
  else {
    cout<<RED<<"Analyzer [INFO]: Gamma PI cut "<<name<<" file NOT found." <<RESET<<endl;
  }

  Alpha_Gate=new TCutG("Alpha_Gate",Nalphacut,x_alphacut,y_alphacut);
  // Alpha_Gate->Print();	
  Gamma_Gate=new TCutG("Gamma_Gate",Ngammacut,x_gammacut,y_gammacut);
  //  Gamma_Gate->Print();
  
  return 0;
}


int Read_DMatrix() {

  ifstream matrix_in(DMatrixFile);
  
  cout <<"Analyzer [INFO]: Reading "<<DMatrixFile<<endl;
  
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
    cout<<GREEN<<"Analyzer [INFO]: Read in Detector Matrix"<<RESET<<endl;
    return 0;
  }
  else {
    cout<<RED<<"Analyzer [ERROR]: Could NOT Read in Detector Matrix"<<RESET<<endl;
    return -1;
  }
  
}


int Read_TMatrix() {

  cout <<"Analyzer [INFO]: Reading "<<TMatrixFile << endl;
  
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
    cout<<GREEN<<"Analyzer [INFO]: Read in TMatrix"<<RESET<<endl;
  }
  else {
    cout<<RED<<"Analyzer [ERROR]: Couldn't open "<<TMatrixFile<<RESET<<endl;
  }
  
  return counter;
}


int Create_Analyzer_Histograms(bool read_binary, bool read_simulation, int NQGates, double QGates[]) {
  
  cout<<"Analyzer [INFO]: Creating Histograms"<<endl;

  //Make Histograms
  hCoinCAEN = new TH2D("CoinCAEN","CoinCAEN",162,0,162,162,0,162);  //coincidence matrix
  hID = new TH1D("hID","hID",256,0,256);
  
  //Time Deviations
  hTimeDev_Rel0 = new TH2D("TimeDev_Rel0","TimeDev_Rel0",10000,-500,500,162,0,162);  //Time deviations relative to crystal 0
  hTimeDev = new TH2D("TimeDev","TimeDev",10000,-500,500,162,0,162);  //Time deviations

  //Diagnostics
  hEventLength = new TH1D("EventLength","EventLength",10000,0,100);
  hTimeBetweenCrystals = new TH2D("TimeBetweenCrystals","TimeBetweenCrystals",10000,0,10000,162,0,162);
  hTimeBetweenDEvents = new TH1D("TimeBetweenDEvents","TimeBetweenDEvents",1000,0,10000);
  hTimeBetweenDEvents_ESum_Mcr = new TH3D("TimeBetweenDEvents_ESum_Mcr","TimeBetweenDEvents_ESum_Mcr",1000,0,10000,500,0,10,20,0,20);
  hTimeBetweenT0s = new TH1D("TimeBetweenT0s","TimeBetweenT0s",1000000,0,100000000);  //Time difference between T0 in ns
  
  //RAW TOF
  hCrystalIDvsTOF = new TH2D("CrystalIDvsTOF","CrystalIDvsTOF",10000,0,1000000,162,0,162); //TOF for each crystal 
  hCrystalTOF = new TH1D("CrystalTOF","CrystalTOF",600000,0,60000000);
  hTOF = new TH1D("TOF","TOF",600000,0,60000000);

  //Corrected TOF
  hCrystalIDvsTOF_Corr = new TH2D("CrystalIDvsTOF_Corrected","CrystalIDvsTOF_Corrected",10000,0,1000000,162,0,162); //TOF for each crystal 
  hCrystalTOF_Corr = new TH1D("CrystalTOF_Corrected","CrystalTOF_Corrected",600000,0,60000000);
  hTOF_Corr = new TH1D("TOF_Corrected","TOF_Corrected",600000,0,60000000);

#ifdef Histogram_DetectorLoad
  if(read_binary==0 && read_simulation==0) {
    hDetectorLoad = new TH1D("DetectorLoad","DetectorLoad",10000000,0,10000000);
  }
#endif

  hDANCE_Entries_per_T0 = new TH1D("DANCE_Entries_per_T0","DANCE_Entries_per_T0",100000,0,100000);
  hDANCE_Events_per_T0 = new TH1D("DANCE_Events_per_T0","DANCE_Events_per_T0",100000,0,100000);
  
  //PSD Histograms
  ADC_calib=new TH2D("ADC_calib","ADC_calib",2400,0.,24.,800,0.,8.);
  ADC_raw=new TH2D("ADC_raw","ADC_raw",1800,0.,72000.,180,0.,7200);
  ADC_calib_Invalid=new TH2D("ADC_calib_Invalid","ADC_calib_Invalid",2400,0.,24.,800,0.,8.);
  
  ADC_raw_ID=new TH3D("ADC_raw_ID","ADC_raw_ID",1800,0.,72000.,180,0.,7200,162,0,162);
  ADC_calib_ID=new TH3D("ADC_calib_ID","ADC_calib_ID",600,0.,24.,200,0.,8.0,162,0,162);


  //Gamma Histograms
  hGamma = new TH2D("hGamma","hGamma",1500,0,30000,162,0,162);
  hGammaCalib = new TH2D("hGammaCalib","hGammaCalib",2000,0.0,20.0,162,0,162);

  //Alpha Histograms
  hAlpha = new TH2D("hAlpha","hAlpha",1500,0,30000,162,0,162);
  hAlphaCalib = new TH2D("hAlphaCalib","hAlphaCalib",500,0.0,5.0,162,0,162);
    
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
  
  for(int i=0;i<21;i++){
    Mbins[i]=0.5+1.*i;
  };
  
  double DEGamma=(GammaE_To-GammaE_From)/GammaE_NoOfBins;
  for(int i=0;i<NoOfEnergyBins+1;i++){
    //EtotBins[i]=i*16./128.;
    EtotBins[i]=GammaE_From+i*DEGamma;
  };
  
  //Beam Monitors
  hU235_TOF = new TH1D("U235_TOF","U235_TOF",600000,0,60000000);  //Raw TOF for U235 Monitor
  hU235_TOF_Corr = new TH1D("U235_TOF_Corrected","U235_TOF_Corrected",600000,0,60000000); //Corrected TOF for U235 Monitor

  hU235_PH_TOF = new TH2D("U235_PH_TOF","U235_PH_TOF",6000,0,60000000,250,0,100000);  //Raw PH vs TOF for U235 Monitor

  

  hU235_PulseHeight = new TH1D("U235_PulseHeight","U235_PulseHeight",10000,0,100000);  //Energy for U235 Monitor
  hU235_En = new TH1D("U235_En","U235_En",NEbins,x);  //Neutron Energy for U235 Monitor 
  hU235_En_Corr = new TH1D("U235_En_Corrected","U235_En_Corrected",NEbins,x);  //Neutron Energy for U235 Monitor (From Corrected TOF)
  hU235_Time_Between_Events = new TH1D("TimeBetweenU235Events","TimeBetweenU235Events",100000,0,10000000);

  hHe3_TOF = new TH1D("He3_TOF","He3_TOF",600000,0,60000000);  //Raw TOF for He3 Monitor
  hHe3_TOF_Corr = new TH1D("He3_TOF_Corrected","He3_TOF_Corrected",600000,0,60000000); //Corrected TOF for He3 Monitor
  hHe3_PulseHeight = new TH1D("He3_PulseHeight","He3_PulseHeight",10000,0,100000);  //Energy for He3 Monitor
  hHe3_En = new TH1D("He3_En","He3_En",NEbins,x);  //Neutron Energy for He3 Monitor 
  hHe3_En_Corr = new TH1D("He3_En_Corrected","He3_En_Corrected",NEbins,x);  //Neutron Energy for He3 Monitor (From Corrected TOF)
  hHe3_Time_Between_Events = new TH1D("TimeBetweenHe3Events","TimeBetweenHe3Events",100000,0,10000000);

  hLi6_TOF = new TH1D("Li6_TOF","Li6_TOF",600000,0,60000000);  //Raw TOF for Li6 Monitor
  hLi6_TOF_Corr = new TH1D("Li6_TOF_Corrected","Li6_TOF_Corrected",600000,0,60000000); //Corrected TOF for Li6 Monitor
  hLi6_PulseHeight = new TH1D("Li6_PulseHeight","Li6_PulseHeight",10000,0,100000);  //Energy for Li6 Monitor
  hLi6_PSD = new TH2D("Li6_PSD","Li6_PSD",600,0,60000,600,0,60000);  //Ifast vs late for Li6 Monitor
  hLi6_En = new TH1D("Li6_En","Li6_En",NEbins,x);  //Neutron Energy for Li6 Monitor 
  hLi6_En_Corr = new TH1D("Li6_En_Corrected","Li6_En_Corrected",NEbins,x);  //Neutron Energy for Li6 Monitor (From Corrected TOF)
  hLi6_Time_Between_Events = new TH1D("TimeBetweenLi6Events","TimeBetweenLi6Events",100000,0,10000000);

  if(read_binary==1 || read_simulation==1) {

    hEn = new TH1D("En","En",NEbins,x); //Raw En
    hEn_Corr = new TH1D("En_Corr","En_Corr",NEbins,x);  //Corrected En 

    En_Esum_Mcl=new TH3F("En_Etot_Mcl","En_Etot_Mcl",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
    En_Esum_Mcr=new TH3F("En_Etot_Mcr","En_Etot_Mcl",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);

    hTOF_Esum_Mcl=new TH3F("TOF_Etot_Mcl","TOF_Etot_Mcl",5100,0,51000000,NoOfEnergyBins,GammaE_From,GammaE_To,20,0,20);
    hTOF_Esum_Mcr=new TH3F("TOF_Etot_Mcr","TOF_Etot_Mcr",5100,0,51000000,NoOfEnergyBins,GammaE_From,GammaE_To,20,0,20);
    
    En_Eg_Mcl=new TH3F("En_Eg_Mcl","En_Eg_Mcl gated on Q",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
    En_Eg_Mcr=new TH3F("En_Eg_Mcr","En_Eg_Mcl gated on Q",NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
    
    Esum_Eg_Mcl=new TH3F("Esum_Eg_Mcl","Esum_Eg_Mcl where Eg is Ecluster",NoOfEnergyBins,EtotBins,NoOfEnergyBins,EtotBins,20,Mbins);
    Esum_Eg_Mcr=new TH3F("Esum_Eg_Mcr","Esum_Eg_Mcr where Eg is Ecrystal",NoOfEnergyBins,EtotBins,NoOfEnergyBins,EtotBins,20,Mbins);


    if(QGated_Spectra) {
      
      for (int kay=0; kay<NQGates; kay++) {
	En_Ecl_Mcl_QGated[kay]=new TH3F(Form("En_Ecl_Mcl_ESum_Gated_%d",kay),
					Form("En_Ecl_Mcl_ESum_Gated_%2.2f_%2.2f",QGates[2*kay],QGates[2*kay+1]),
					NEbins,x,NoOfEnergyBins,EtotBins,20,Mbins);
      }
    }
  }

  cout<<GREEN<<"Analyzer [INFO]: Created Histograms"<<RESET<<endl;
  return 0;

}


int Write_Analyzer_Histograms(TFile *fout, bool read_binary, bool read_simulation, int NQGates, double QGates[]) {
  
  cout<<"Analyzer [INFO]: Writing Histograms"<<endl;
  
  fout->cd();
  hID->Write();
  hCoinCAEN->Write();
  hTimeDev_Rel0->Write();
  hTimeDev->Write();

  hEventLength->Write();
  hTimeBetweenCrystals->Write();
  hTimeBetweenDEvents->Write();
  hTimeBetweenDEvents_ESum_Mcr->Write();
  hTimeBetweenT0s->Write();
  hCrystalIDvsTOF_Corr->Write();
  hCrystalTOF_Corr->Write();
  hCrystalIDvsTOF->Write();
  hCrystalTOF->Write();

  hDANCE_Entries_per_T0->Write();
  hDANCE_Events_per_T0->Write();


#ifdef Histogram_DetectorLoad
  if(read_binary==0 && read_simulation==0) { 
    for(int eye=0; eye<10000000; eye++) {
      hDetectorLoad->SetBinContent(eye+1,Detector_Load[eye]/(160.0*T0_Counter));
    }
    hDetectorLoad->Write();
  }
#endif


  ADC_calib->Write();
  ADC_raw->Write();
  ADC_calib_Invalid->Write();

  ADC_raw_ID->Write();
  ADC_calib_ID->Write();

  hAlpha->Write();
  hAlphaCalib->Write();

  hGamma->Write();
  hGammaCalib->Write();

  Alpha_Gate->Write();
  Gamma_Gate->Write();

  if(read_binary==1 || read_simulation == 1) {
    hTOF->Write();
    hTOF_Corr->Write();

    hEn->Write();
    hEn_Corr->Write();

    En_Esum_Mcl->Write();
    En_Esum_Mcr->Write();

    hTOF_Esum_Mcl->Write();
    hTOF_Esum_Mcr->Write();

    if(QGated_Spectra) {
      for (int kay=0; kay<NQGates; kay++) {
	En_Ecl_Mcl_QGated[kay]->Write();
      }
    }
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

  hHe3_TOF->Write();  //Raw TOF for He3 Monitor
  hHe3_TOF_Corr->Write(); //Corrected TOF for He3 Monitor
  hHe3_PulseHeight->Write();  //Energy for He3 Monitor
  hHe3_En->Write();  //Neutron Energy for He3 Monitor 
  hHe3_En_Corr->Write();  //Neutron Energy for He3 Monitor (From Corrected TOF)
  hHe3_Time_Between_Events->Write();

  cout<<GREEN<<"Analyzer [INFO]: Wrote Histograms"<<RESET<<endl;
  return 0;
}


int Initialize_Analyzer(bool read_binary, bool write_binary, bool read_simulation, int NQGates, double QGates[]) {

  cout<<BLUE<<"Analyzer [INIT]: Initializing Analyzer"<<RESET<<endl;
  
  //Initialize Analysis Stuff

  events_analyzed=0;
  entries_analyzed=0;
  analyzer_counter=1;

  Read_Moderation_Time_Graphs();
  Read_PI_Gates();
  totalindex = Read_TMatrix();
  Read_DMatrix();
  Create_Analyzer_Histograms(read_binary,read_simulation,NQGates,QGates);
  
  for(int eye=0; eye<162; eye++) {
    last_timestamp[eye]=0;
    last_valid_timestamp[eye]=0;
    current_timestamp[eye]=0;

    last_t0_timestamp=0;
    current_t0_timestamp=0;

    last_devent_timestamp=0;
    last_valid_devent_timestamp=0;

  }

  //Diagnostics
  T0_Counter=0;

#ifdef Histogram_DetectorLoad
  for(int eye=0; eye<100000000; eye++) {
    Detector_Load[eye]=0;
  }
#endif

  DANCE_Entries_per_T0=0;
  DANCE_Events_per_T0=0;

  return 0;
}


int Analyze_Data(std::vector<DEVT_BANK> eventvector, bool read_binary, bool write_binary, bool read_simulation, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold, int NQGates, double QGates[]) {


  //Progress indicator
  if( entries_analyzed > analyzer_counter*ProgressInterval) {
    cout<< events_analyzed<<" Events Comprised of "<<entries_analyzed<<" Entries Analyzed. "<<
      "Average Mult: "<<1.0*entries_analyzed/(1.0*events_analyzed)<<endl;
    analyzer_counter++;
  }
  int Crystal_Mult=0;

  //Initialize DANCE Event
  devent.Crystal_mult=0;
  devent.ESum=0;

  //No Physics events are valid until created
  devent.Valid=0;
  u235event.Valid=0;
  he3event.Valid=0;
  li6event.Valid=0;
    
  
  //Fill Event Length (mult 1 events just give 0 so require mult>1 to fill this)
  //  if(eventvector.size() > 1) {
  hEventLength->Fill(eventvector[eventvector.size()-1].TOF-eventvector[0].TOF);
  // }
  // if(eventvector.size() > 10) {
  //  cout<<eventvector.size()<<endl;
  //   cout<<eventvector[eventvector.size()-1].TOF-eventvector[0].TOF<<endl;
  // }
  
  //Loop over event 
  for(uint32_t eye=0; eye<eventvector.size(); eye++) {

    //Fill ID histogram
    hID->Fill(eventvector[eye].ID);
    int id_eye=eventvector[eye].ID;
    
    //Place the current time in the 
    current_timestamp[id_eye] = eventvector[eye].TOF;
    
    //this is a dance crystal
    if(id_eye<162) {

      //Apply Energy Calibration
      double temp_slow = eventvector[eye].Islow + gRandom->Uniform(0,1);
      double temp_fast = eventvector[eye].Ifast + gRandom->Uniform(0,1);
      
      eventvector[eye].Eslow = 0.001*(temp_slow*temp_slow*slow_quad[eventvector[eye].ID] +
				      temp_slow*slow_slope[eventvector[eye].ID] +
				      slow_offset[eventvector[eye].ID]);
      
      eventvector[eye].Efast = 0.001*(temp_fast*fast_slope[eventvector[eye].ID] +
				      fast_offset[eventvector[eye].ID]);
      
      //  cout<<eventvector[eye].Eslow<<"  "<<temp_slow<<endl;

      if(HAVE_Threshold==1) {
	if(eventvector[eye].Eslow < Energy_Threshold) {
	  eventvector[eye].Valid=0;
	}
      }
      
      //Coincidences
      if(eventvector.size() > 1 && eye < (eventvector.size()-1)) {
	
	//start with the next one
	for(uint32_t jay=eye+1; jay<eventvector.size(); jay++) {
	  
	  int id_jay = eventvector[jay].ID;
	  
	  if(id_jay<162) {

	    //time difference between crystal jay and eye
	    double ddT = eventvector[jay].TOF - eventvector[eye].TOF;
	    
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
      
      //Fill time between crystal hits
      hTimeBetweenCrystals->Fill((current_timestamp[id_eye]-last_timestamp[id_eye]),id_eye,1);
      
      //Blocking Time to avoid retrigger problems
      if((current_timestamp[id_eye]-last_valid_timestamp[id_eye]) < Crystal_Blocking_Time) {
	eventvector[eye].Valid=0;
      }
      
      if(eventvector[eye].Valid==1) {
	//Fill ADC Calib
	ADC_calib->Fill(eventvector[eye].Eslow, eventvector[eye].Efast,1);
	ADC_raw->Fill(eventvector[eye].Islow, eventvector[eye].Ifast,1);

	ADC_calib_ID->Fill(eventvector[eye].Eslow, eventvector[eye].Efast, id_eye,1);
	ADC_raw_ID->Fill(eventvector[eye].Islow, eventvector[eye].Ifast,id_eye,1);

	//First check to see if it is in the alpha gate
	if(Alpha_Gate->IsInside(eventvector[eye].Eslow, eventvector[eye].Efast)) {
	  //Fill some Calibration Spectra
	  hAlpha->Fill(eventvector[eye].Islow, eventvector[eye].ID,1);
	  hAlphaCalib->Fill(eventvector[eye].Eslow, eventvector[eye].ID,1);
	}
	
	//If it is not in the alpha gate then check the gamma gate or whether its simulated data
	else {
	  if((Gamma_Gate->IsInside(eventvector[eye].Eslow, eventvector[eye].Efast)) || read_simulation) {
	  
	    //Fill some Calibration Spectra
	    hGamma->Fill(eventvector[eye].Islow, eventvector[eye].ID,1);
	    hGammaCalib->Fill(eventvector[eye].Eslow, eventvector[eye].ID,1);
	    
	    //Make a DANCE Event
	    devent.Crystal_ID[Crystal_Mult] = eventvector[eye].ID;  //Crystal ID
	    devent.Cluster_ID[Crystal_Mult] = Crystal_Mult+1;  //??????
	    devent.Islow[Crystal_Mult] = eventvector[eye].Islow;  //Crystal long integral
	    devent.Ifast[Crystal_Mult] = eventvector[eye].Ifast;  //Crystal short integral
	    devent.tof[Crystal_Mult] = eventvector[eye].TOF;  //time of flight
	    devent.Ecrystal[Crystal_Mult] = eventvector[eye].Eslow;   //Energy if calibrated 
	    devent.ESum += eventvector[eye].Eslow; //ESum 
	    devent.Crystal_mult++;
	    devent.Valid=1;  //event is now valid
	    Crystal_Mult++;	
	  }
	}
	
	DANCE_Entries_per_T0++;
	
      } //end of Valid Check
      else {
	ADC_calib_Invalid->Fill(eventvector[eye].Eslow, eventvector[eye].Efast,1);
      } //end of not valid else
      
      if(read_binary==0 && read_simulation==0) {

#ifdef Histogram_DetectorLoad
	//Fill detecotor load
	if(last_t0_timestamp > 0) {
	  uint32_t temptof = (uint32_t)(eventvector[eye].TOF-last_t0_timestamp);
	  if(temptof>0) {
	    for(int el=(int)temptof; el<((int)(temptof+1000)); el++) {
	      if(el >=0 && el < 100000000) {
		Detector_Load[el]+=1.0;
	      }
	    }
	  }
	}
#endif

      }
    }
      
    
    //this is t0
    if(id_eye==200) {
      
      current_t0_timestamp = eventvector[eye].TOF;
      //Do some T0 diagnostics
      if(last_t0_timestamp > 0) {
	hTimeBetweenT0s->Fill(current_t0_timestamp-last_t0_timestamp);
      }
      
      if(current_t0_timestamp > last_t0_timestamp + 1000000) {

	//Fill Histos
	hDANCE_Entries_per_T0->Fill(T0_Counter,DANCE_Entries_per_T0);
	hDANCE_Events_per_T0->Fill(T0_Counter,DANCE_Events_per_T0);
	
	//Reset Entries and Events per T0
	DANCE_Entries_per_T0=0;
	DANCE_Events_per_T0=0;

	//Incriment T0 counter
	T0_Counter++;
	//Set the old timestamp value once done
	last_t0_timestamp=current_t0_timestamp;
	
      }  
      
    }
   
    //beam monitors
    //He3
    if(id_eye == 241) {
      he3event.tof = eventvector[eye].TOF;
      he3event.Ifast = eventvector[eye].Ifast;
      he3event.Islow = eventvector[eye].Islow;
      he3event.Valid = 1; //he3 event now valid

      //Fill time Diagnostics
      hHe3_Time_Between_Events->Fill(current_timestamp[id_eye]-last_timestamp[id_eye]);
    }
    
    //U235
    if(id_eye == 243) {
      u235event.tof = eventvector[eye].TOF;
      u235event.Ifast = eventvector[eye].Ifast;
      u235event.Islow = eventvector[eye].Islow;
      u235event.Valid = 1; //u235 event now valid

      //Fill time Diagnostics
      hU235_Time_Between_Events->Fill(current_timestamp[id_eye]-last_timestamp[id_eye]);
    }
    
    //Li6
    if(id_eye == 244) {
      li6event.tof = eventvector[eye].TOF;
      li6event.Ifast = eventvector[eye].Ifast;
      li6event.Islow = eventvector[eye].Islow;
      li6event.Valid = 1; //li6 event now valid
      
      //Fill time Diagnostics
      hLi6_Time_Between_Events->Fill(current_timestamp[id_eye]-last_timestamp[id_eye]);
    }
    
    /*Non paralyzable model of dead time.  
    //The valid timestamp is not updated unless the current detector event 
    falls after the dead-time window following the preveious event */
    
    if(eventvector[eye].Valid==1) {
      //Once done the time is the last timestamp
      last_valid_timestamp[id_eye] = current_timestamp[id_eye];
    }
    last_timestamp[id_eye] = current_timestamp[id_eye];

    //Write to Binary
    if(write_binary==1 && outputbinfile.is_open()) {
      devt_out.Ifast = eventvector[eye].Ifast;
      devt_out.Islow = eventvector[eye].Islow;
      devt_out.TOF = eventvector[eye].TOF;
      devt_out.ID = eventvector[eye].ID;
      outputbinfile.write(reinterpret_cast<char*>(&devt_out),sizeof(DEVT_STAGE1));
    }    

    entries_analyzed++;
  } //end of eventvector loop
  
  
  //Handle various events and do some physics
  
  //DANCE Events
  if(devent.Valid == 1) {
    if(last_t0_timestamp > 0) {
      
      //DANCE Event Diagnostics
      DANCE_Events_per_T0++;

      //DEvent Blocking Time
      if((devent.tof[0]-last_valid_devent_timestamp) < DEvent_Blocking_Time) {
	devent.Valid=0;
      }
      else {
	hTimeBetweenDEvents->Fill(devent.tof[0]-last_devent_timestamp);
	hTimeBetweenDEvents_ESum_Mcr->Fill(devent.tof[0]-last_devent_timestamp,devent.ESum,devent.Crystal_mult);

	//Update the last valid devent timestamp.  Same non-paralyzable model
	last_valid_devent_timestamp=devent.tof[0];
      }
     
      //Set Last devent timestamp to current one
      last_devent_timestamp=devent.tof[0];
      
      //TOF now relative to last T0
      for(int kay=0; kay<devent.Crystal_mult; kay++) {
	
	//Make TOF relative to last T0
	devent.tof[kay] -= last_t0_timestamp;
	
	//Correct the TOF for moderation time between ~0 and ~10 MeV
	if(devent.tof[kay] >= DANCE_TOF_Corr_Limit[0] && devent.tof[kay] <= DANCE_TOF_Corr_Limit[1]) {
	  devent.tof_corr[kay] = gr_DANCE_TOF_Corr->Eval(devent.tof[kay]);
	}
	else {
	  devent.tof_corr[kay] = -1;
	  devent.Valid=0;
	}
	
	//Fill the 2D spectrum of TOF vs Crystal
	hCrystalIDvsTOF->Fill(devent.tof[kay],devent.Crystal_ID[kay],1);
	hCrystalTOF->Fill(devent.tof[kay],1);
	
	//Fill the 2D spectrum of Corrected TOF vs Crystal
	hCrystalIDvsTOF_Corr->Fill(devent.tof_corr[kay],devent.Crystal_ID[kay],1);
	hCrystalTOF_Corr->Fill(devent.tof_corr[kay],1);
	
      } //end loop over devent crystals
    } //end check of crystal mult and last_t0_timestamp
    else {
      devent.Valid = 0;
    }
  } //end check of devent.Valid
  
  if(read_binary==1 || read_simulation ==1 ) {
    
    //Handle various events and do some physics
    if(devent.Valid==1) {
      
      //Calculate the neutron energy    
      devent.En = 0.5*939.565379e6*DANCE_FlightPath*DANCE_FlightPath/((devent.tof[0])/1e9)/((devent.tof[0])/1e9)/(2.997924589e8*2.997924589e8); 
      //Calculate corrrected nuetron energy
      devent.En_corr = 0.5*939.565379e6*DANCE_FlightPath*DANCE_FlightPath/((devent.tof_corr[0])/1e9)/((devent.tof_corr[0])/1e9)/(2.997924589e8*2.997924589e8); 
      
      //Fill TOF and En
      hEn->Fill(devent.En,1);
      hEn_Corr->Fill(devent.En_corr,1);
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
            
      //Mults vs ESum vs En
      En_Esum_Mcl->Fill(devent.En_corr,devent.ESum,devent.Cluster_mult);
      En_Esum_Mcr->Fill(devent.En_corr,devent.ESum,devent.Crystal_mult);
      
      //Mults vs ESum vs TOF
      hTOF_Esum_Mcl->Fill(devent.tof_corr[0],devent.ESum,devent.Cluster_mult);
      hTOF_Esum_Mcr->Fill(devent.tof_corr[0],devent.ESum,devent.Crystal_mult);
      
      //Loop over the cluster mult
      for(int jay=0; jay<devent.Cluster_mult; jay++ )  {
	for(int kay=0; kay<NQGates; kay++) {
	  if(devent.ESum > QGates[kay*2] && devent.ESum < QGates[kay*2+1]) {
	    En_Ecl_Mcl_QGated[kay]-> Fill(devent.En, devent.Ecluster[jay], devent.Cluster_mult );
	  } //Done Checking QGates
	} //Done Looping over QGates
      } //Done looping over cluster mult
    }
  } //Done checking for stage 1
  
  //U235 Beam Monitor Events
  if(u235event.Valid == 1) {
    if(last_t0_timestamp > 0) {
      //Fill Pulse Height Spectrum
      hU235_PulseHeight->Fill(u235event.Islow);
      //Make TOF relative to last T0
      u235event.tof -= last_t0_timestamp;
      //Correct the TOF for moderation time between ~0 and ~10 MeV
      if(u235event.tof >= U235_TOF_Corr_Limit[0] && u235event.tof <= U235_TOF_Corr_Limit[1]) {
	u235event.tof_corr = gr_U235_TOF_Corr->Eval(u235event.tof);
      }
      else {
	u235event.tof_corr = -1;
	u235event.Valid=0;
      }
      //Fill TOF spectra
      hU235_TOF->Fill(u235event.tof);
      hU235_TOF_Corr->Fill(u235event.tof_corr);
      hU235_PH_TOF->Fill(u235event.tof,u235event.Islow);
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
    if(last_t0_timestamp > 0) {
      //Fill Pulse Height Spectrum
      hHe3_PulseHeight->Fill(he3event.Islow);
      //Make TOF relative to last T0
      he3event.tof -= last_t0_timestamp;
      //Correct the TOF for moderation time between ~0 and ~10 MeV
      if(he3event.tof >= He3_TOF_Corr_Limit[0] && he3event.tof <= He3_TOF_Corr_Limit[1]) {
	he3event.tof_corr = gr_He3_TOF_Corr->Eval(he3event.tof);
      }
      else {
	he3event.tof_corr = -1;
	he3event.Valid=0;
      }
      //Fill TOF spectra
      hHe3_TOF->Fill(he3event.tof);
      hHe3_TOF_Corr->Fill(he3event.tof_corr);
      
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
    if(last_t0_timestamp > 0) {
      //Fill Pulse Height Spectrum
      hLi6_PulseHeight->Fill(li6event.Islow-li6event.Ifast);
      hLi6_PSD->Fill(li6event.Islow,li6event.Ifast);
      //Make TOF relative to last T0
      li6event.tof -= last_t0_timestamp;
      //Correct the TOF for moderation time between ~0 and ~10 MeV
      if(li6event.tof >= Li6_TOF_Corr_Limit[0] && li6event.tof <= Li6_TOF_Corr_Limit[1]) {
	li6event.tof_corr = gr_Li6_TOF_Corr->Eval(li6event.tof);
      }
      else {
	li6event.tof_corr = -1;
	li6event.Valid=0;
      }
      //Fill TOF spectra
      hLi6_TOF->Fill(li6event.tof);
      hLi6_TOF_Corr->Fill(li6event.tof_corr);
      
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

  //Done with Event Processing
  //increment counter
  events_analyzed++;

  return 0;
}
