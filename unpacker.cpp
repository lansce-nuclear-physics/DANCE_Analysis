
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
//*  unpacker.cpp           *// 
//*  Last Edit: 05/18/22    *//  
//***************************//

//File includes
#include "global.h"
#include "unpacker.h"
#include "unpack_vx725_vx730.h"
#include "sort_functions.h"
#include "eventbuilder.h"
#include "structures.h"
#include "analyzer.h"
#include "validator.h"

//C/C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include <iomanip>
#include <string>
#include <fstream>

//ROOT includes
#include "TH1.h"
#include "TH3.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom.h"

using namespace std;

stringstream umsg;

//output diagnostics file
ofstream outputdiagnosticsfile;

//global unpacker variables
double TimeDeviations[200];

//holds the id of each detector from dancemap
int MapID[20][20];

//function return
int func_ret = 0;

//Histograms
TH3S *hWaveform_ID;
TH3S *hWaveform_ID_NR;
TH2S *hWaveform_Li6;
TH2S *hWaveform_U235;
TH2S *hWaveform_Bkg;
TH2S *hWaveform_He3;
TH2S *hWaveform_T0;

TH2C *hDigital_Probe1_ID;
TH2C *hDigital_Probe2_ID;

TH2F *hID_vs_WFRatio;
TH3F *hID_vs_WFInt_vs_Islow;
TH3F* hID_vs_WFRatio_vs_Islow;
//TH2F *hID_vs_WFRatio_check = new TH2F("WFRatio_ID_check","WFRatio_ID_check",1000,-0.2,0.8,162,0,162); //IK added
TH1I *hScalers;

TH1D* hTimestamps;
TH1D* hTimestampsT0;
TH1D* hTimestampsBM;
TH2D* hTimestampsID; 

TH1S *hWaveforms[20];
int waveform_counter=0;

TH1I *hEventID;



//Histograms for Unpacker Things
int Create_Unpacker_Histograms(Input_Parameters input_params) {
  
  DANCE_Info("Unpacker","Creating Histograms");
  
  hEventID=new TH1I("EventID","EventID",20,0,20);

  //Make a histogram for waveforms
  if(input_params.Read_Binary==0) {
    for(int eye=0; eye<20; eye++){
      hWaveforms[eye] = new TH1S(Form("Waveform%d",eye),Form("Waveform%d",eye),80,0,80);
    }

#ifdef Histogram_Waveforms 
    hWaveform_ID = new TH3S("Waveform_ID","Waveform_ID",40,0,40,2000,0,20000,162,0,162);
    hWaveform_ID_NR = new TH3S("Waveform_ID_NR","Waveform_ID_NR",40,0,40,2000,0,20000,162,0,162);
    hWaveform_T0 = new TH2S("T0_Waveform","T0_Waveform",200,0,200,2000,0,20000);
    hWaveform_Li6 = new TH2S("Li6_Waveform","Li6_Waveform",600,0,600,2000,0,20000);
    hWaveform_U235 = new TH2S("U235_Waveform","U235_Waveform",600,0,600,2000,0,20000);
    hWaveform_Bkg = new TH2S("Bkg_Waveform","Bkg_Waveform",600,0,600,2000,0,20000);
    hWaveform_He3 = new TH2S("He3_Waveform","He3_Waveform_He3",600,0,600,2000,0,20000);
    hID_vs_WFRatio = new TH2F("WFRatio_ID","WFRatio_ID",1000,-0.2,0.8,162,0,162);
    hID_vs_WFInt_vs_Islow = new TH3F("WFInt_Islow_ID","WFInt_Islow_ID",500,-1000,4000,2000,0,40000,162,0,162);
    hID_vs_WFRatio_vs_Islow = new TH3F("WFRatio_Islow_ID","WFRatio_Islow_ID",500,-0.2,0.8,2000,0,40000,162,0,162);
  
#endif

#ifdef MakeTimeStampHistogram
     hTimestamps = new TH1D("Timestamps","Timestamps",1e6,0,4000);
     hTimestampsT0 = new TH1D("T0_TimestampsT0","T0_Timestamps",1e6,0,4000);
     hTimestampsBM = new TH1D("BM_Timestamps","BM_Timestamps",1e6,0,4000);
     hTimestampsID = new TH2D("Timestamps_ID","Timestamps_ID",4e5,0,1000,162,0,162);
#endif

#ifdef Histogram_Digital_Probes
    hDigital_Probe1_ID = new TH2C("ID_DigitalProbe1","ID_DigitalProbe1",600,0,600,256,0,256);
    hDigital_Probe2_ID = new TH2C("ID_DigitalProbe2","ID_DigitalProbe2",600,0,600,256,0,256);
#endif
    
    //Scalers
    hScalers = new TH1I("Scalers","Scalers",35,0,35);
  }


    
  DANCE_Success("Unpacker","Created Unpacker Histograms");
  return 0; 
  
}

int Write_Unpacker_Histograms(TFile *fout, Input_Parameters input_params) {
  
  DANCE_Info("Unpacker","Writing Histograms");
  
  fout->cd();
  hEventID->Write();
  if(input_params.Read_Binary==0) {
#ifdef Histogram_Waveforms
    for(int eye=0; eye<20; eye++){
      hWaveforms[eye]->Write();
    }
    hWaveform_ID->Write();
    hWaveform_ID_NR->Write();
    hWaveform_T0->Write();
    hWaveform_Li6->Write();
    hWaveform_U235->Write();
    hWaveform_Bkg->Write();
    hWaveform_He3->Write();
    hID_vs_WFRatio->Write();
    hID_vs_WFRatio_vs_Islow->Write();
    hID_vs_WFInt_vs_Islow->Write();
#endif

#ifdef Histogram_Digital_Probes
    hDigital_Probe1_ID->Write();
    hDigital_Probe2_ID->Write();
#endif

#ifdef MakeTimeStampHistogram
     hTimestamps->Write();
     hTimestampsT0->Write();
     hTimestampsBM->Write();
     hTimestampsID->Write();
#endif

    hScalers->Write();
  }
  
  DANCE_Success("Unpacker","Wrote Unpacker Histograms");
  return 0;
}

int Write_Root_File(Input_Parameters input_params, Analysis_Parameters *analysis_params){
  //output the rootfile
  //Name of the output root file
  stringstream rootfilename;
  rootfilename.str();
  
  //stage 0
  if(input_params.Analysis_Stage==0 && input_params.Read_Simulation==0) {
    rootfilename << STAGE0_ROOT << "/Stage0_Histograms_Run_";
    rootfilename << input_params.RunNumber;
    if (input_params.SingleSubrun){
      rootfilename << "_";
      rootfilename << input_params.SubRunNumber;
    }
  }
  //stage 1
  else if(input_params.Analysis_Stage==1 && input_params.Read_Simulation==0) {
    rootfilename << STAGE1_ROOT << "/Stage1_Histograms_Run_";
    rootfilename << input_params.RunNumber;
    if (input_params.SingleSubrun){
      rootfilename << "_";
      rootfilename << input_params.SubRunNumber;
    }
  }
  //simulation
  else if(input_params.Analysis_Stage==1 && input_params.Read_Simulation==1) {
    rootfilename << STAGE1_ROOT << "/Stage1_Histograms_";
    rootfilename << input_params.Simulation_File_Name;
  }
  else {
    DANCE_Error("Unpacker","Cannot understand options for making rootfile. Exiting!");
    return -1;
  }
    
  //make the root file  
  TFile *fout;


  //Make the file
  fout = new TFile(Form("%s_%dns_CW_%dns_CBT_%dns_DEBT.root",rootfilename.str().c_str(),
                        (int)input_params.Coincidence_Window,
                        (int)input_params.Crystal_Blocking_Time,
                        (int)input_params.DEvent_Blocking_Time),"RECREATE");
  
  DANCE_Success("Unpacker","Rootfile Created");
  
  //Write histograms
  //hID_vs_WFRatio_check->Write(); //IK added
  Write_Unpacker_Histograms(fout, input_params);
  Write_PI_Gates(fout);
  Write_Eventbuilder_Histograms(fout, input_params, analysis_params);
  Write_Analyzer_Histograms(fout, input_params);

  //Write the root file
  fout->Write();
  DANCE_Success("Unpacker","Rootfile Written");
  return 0;
}

//Functions
int Make_DANCE_Map() {
  
  DANCE_Info("Unpacker","Reading DANCE Map");
  
  ifstream map;
  map.open(DanceMapFile);
  int m1,m2,m3;
  
  if(map.is_open()) {
    while(!map.eof()){
      map>>m1>>m2>>m3;
      MapID[m2][m1]=m3;
      if(map.eof()) {
        break;
      }
    }
    map.close();
    DANCE_Success("Unpacker","Created DANCE Map");
    return 0;
  }
  else {

    umsg.str("");
    umsg<<"Could not Open DANCE Map File: "<<DanceMapFile;
    DANCE_Error("Unpacker",umsg.str());

    return -1;
  }  
}

//This function makes the output binary file for stage 0 or 1 filled with time-ordered devt_bank structures
int Make_Output_Diagnostics_File(int RunNumber) {
  
  stringstream outfilename;
  outfilename.str();
  
  outfilename << DIAGNOSTICS;
  outfilename <<"/diagnostics_run";
  outfilename << RunNumber;
  outfilename << ".txt";
  
  outputdiagnosticsfile.open(outfilename.str().c_str(), ios::out);
  
  if(outputdiagnosticsfile.is_open()) {
    DANCE_Success("Unpacker","Created Output Diagnostics File");
    return 0;
  }
  else {
    stringstream umsg;
    umsg.str("");
    umsg<<"Could not Create Output Diagnostics File: "<<outfilename.str();
    DANCE_Error("Unpacker",umsg.str());

    return -1;
  }
}

int Read_TimeDeviations(Input_Parameters input_params) {

  DANCE_Info("Unpacker","Reading Time Deviations");
  
  for(int eye=0; eye<200; eye++) {
    TimeDeviations[eye]=0;
  }  
  
  //If we want to fit the time deviations then set all of them to zero
  if(input_params.FitTimeDev) {
    DANCE_Info("Unpacker","Time Deviations set to 0");
  }
  
  //if we have time deviations then obtain them from the files
  else if(input_params.Read_Simulation==0) {
    stringstream fname;
    fname.str();
    
    fname<<TIMEDEV_DIR <<"/TimeDeviations_Run_" << input_params.RunNumber << ".txt";
    
    ifstream fdata;
    fdata.open(fname.str().c_str());
    
    if(fdata.is_open()) {
      int tempid;
      double tempoffset;
      while(!fdata.eof()) {
        fdata >> tempid >> tempoffset; 
        TimeDeviations[tempid] = tempoffset;
        if(fdata.eof()) {
          break;
        }
      }
      
      stringstream umsg;
      umsg.str("");
      umsg<<"Time Deviations Read from: "<<fname.str();
      DANCE_Success("Unpacker",umsg.str()); 
    
      return 0;
    }
    else {

      stringstream umsg;
      umsg.str("");
      umsg<<"Failed to Read Time Deviations From: "<<fname.str();
      DANCE_Error("Unpacker",umsg.str()); 
    }
  }
  return 0;
}


double Calculate_Fractional_Time(uint16_t waveform[], uint32_t Ns, uint8_t dual_trace, uint16_t model, Analysis_Parameters *analysis_params) {

  // CALCULATE THE LEADING EDGE using constant fraction "frac"
  uint32_t imin=0;
  double sigmin=1e9;
  double frac=0.04;
  double base=0;
  uint32_t NNN=4; 
  
  //fraction
  frac=0.2;
  
  //loop over the waveform
  for(uint32_t kay=0; kay<Ns; kay++) {
    
    //calculate the baseline 
    if(kay<NNN) {
      base+=(1.*waveform[kay]);
    }                  
   
    //locate the minimum in the waveform (negative polarity so this is the max)
    if((1.*waveform[kay])<sigmin) {
      sigmin=1.*waveform[kay];
      imin=kay;
    }
  }
  
  //Divide baseline by number of samples integrated
  base /= (1.0*NNN);

  
  //Threshold set to a fraction of the pulse minimum 
  double thr=(sigmin-base)*frac+base;
  
  //Fractional time initialized to zero
  double dT=0;
  
  //Leading edge value
  int iLD=0;
  
  //Go backward from the minimum i value
  for(int kay=imin;kay>1;kay--){
    
    //Look for where the waveform crosses the threshold
    if((1.*waveform[kay])<thr && (1.*waveform[kay-1])>thr){
      
      //difference in the signal height about the crossing of the threshold
      double dSig=(1.*waveform[kay-1]-1.*waveform[kay]);
      
      //When there is no dual trace there is 
      if(!dual_trace && model == 730) {
        if(dSig!=0) dT=(1.*waveform[kay-1]-thr)/dSig*2.+(kay-1)*2.;  // this is in ns
        else dT=(kay-1)*2.;
      }
      //When dual trace is on there is 4 ns between samples
      else if((!dual_trace && model == 725) || (dual_trace && model ==730)) {              
        if(dSig!=0) dT=(1.*waveform[kay-1]-thr)/dSig*4.+(kay-1)*4.;  // this is in ns
        else dT=(kay-1)*4.;
      }
      else if(dual_trace && model == 725) {
        if(dSig!=0) dT=(1.*waveform[kay-1]-thr)/dSig*8.+(kay-1)*8.;  // this is in ns
        else dT=(kay-1)*8.;
      }
      else {
        stringstream umsg;
        umsg.str("");
        umsg<<"Not sure what to do with dual trace: "<<dual_trace<<"  and model: "<<model;
        DANCE_Error("Unpacker",umsg.str());

        return -1;
      }
      iLD=kay;
    }                
  }      

  //Use dT to calculate the integral of the end of the waveform
  double wf_counter=0.0;
  double integral=0.0;
  for(int kay=(int)dT/2.0+10; kay<(int)Ns; kay++) {
    integral += base-waveform[kay];    
    wf_counter += 1.0;
  }
  integral /= wf_counter;

  analysis_params->wf_integral=integral;

  return dT;
}

int Unpack_Data(queue<gzFile> &gz_queue, double begin, Input_Parameters input_params, Analysis_Parameters *analysis_params) {

  gzFile gz_in=gz_queue.front();
  ofstream faillog;
  faillog.open("Readout_Status_Failures.txt", ios::app);

  //Boolean control variables
  bool run=true; 
  bool subrun=true;
  //Structures to put data in
  deque<DEVT_BANK> datadeque;                         //Storage container for time sorted data
  
  //CAEN 2015 unpacking
  long devt_padding = 0;                              // padding between banks not divisible by 64 bits
  short imported_peaks[256][16384];                   // this is actually supported channels / supported length of PXXX bank
  unsigned short int wf1[15000];
  CEVT_BANK *evinfo = new CEVT_BANK();                //caen event info
  test_struct_cevt *evaggr = new test_struct_cevt();  //event aggregate
  DEVT_BANK *db_arr = new DEVT_BANK[MaxDEVTArrSize];  //Storage array for entries
  DEVT_STAGE1_WF devt_stage1_wf;                      //Stage1 format for reading binary with WF Integral
  DEVT_STAGE1 devt_stage1;                            //Stage1 format for reading binary without WF Integral

  //CAEN 2018 unpacking
  User_Data_t user_data;                              //This is the fw version and user extra word storage
  Vx725_Vx730_Board_Data_t vx725_vx730_board_data;    //This is the board header data for the Vx725 and Vx730 boards
  Vx725_Vx730_PSD_Data_t vx725_vx730_psd_data;        //This is the PSD data for the Vx725 and Vx730 boards
  Vx725_Vx730_PHA_Data_t vx725_vx730_pha_data;        //This is the PHA data for the Vx725 and Vx730 boards

  //Counters
  uint32_t EVTS=0;              //Total number of entries unpacked since last time sort
  uint64_t BYTES_READ=0;        //Total number of Bytes read since last time sort
  uint64_t TOTAL_BYTES=0;       //Total number of Bytes read 
  uint32_t progresscounter=1;   //Keep track of how many progress statements have been made
  int gzret=1;                  //number of bytes read by gzread
  

  //Some CAEN structures to optimize reads
  V1730_Header_t v1730_header;
  V1730_ChAgg_Header_t v1730_chagg_header;
  uint32_t v1730_chagg_data[Max_ChAgg_Size];

  //MIDAS Bank Stuff
  EventHeader_t head;           //MIDAS event header
  BankHeader_t bhead;           //MIDAS bank header
  Bank32_t bank32;              //MIDAS 32-bit bank
  Bank_t bank;                  //MIDAS 16-bit bank

  uint32_t TotalDataSize=0;     // head.fDataSize;
  uint32_t TotalBankSize=0;     // bhead.fDataSize;
  uint32_t EventBankSize=0;     // bank32fDataSize;

  //Scalers
  Sclr_Totals_t sclr_totals;    //Scaler Totals
  Sclr_Rates_t sclr_rates;      //Scaler Rates
  
  //time profiling for performance
  struct timeval tv;              //Real time  
  double time_elapsed;     //Elapsed time
  double time_elapsed_old; //Previous elapsed time
  
  //Scaler Variables for caen2018 format
  struct timeval timevalue;          //Time at which scalers were recorded
  uint32_t Digitizer_Rates[20];      //Digitizer read rates in bytes per second
  uint16_t ADC_Temp[20][16];         //0x1nA8 ADC Temps in degrees C
  uint32_t Channel_Status[20][16];   //0x1n88 Channel status registers
  uint32_t Acquisition_Status[20];   //0x8104 Acquisition Status
  uint32_t Failure_Status[20];       //0x8178 Board Failure Status
  uint32_t Readout_Status[20];       //0xEF04 Readout Status
  uint32_t Register_0x8504n[20][8];  //0x8500 + 4n (Tells how many buffers are left to readout in each pair)
  uint32_t Register_0x1n2C[20][16];

  //waveform ratio gates
  char gatename[200]; 
  double wf_ratio_low, wf_ratio_high;
  sprintf(gatename,"Gates/%s",PILEUPGATE);
  ifstream pileupcutin(gatename);
  if(pileupcutin.is_open()){
    pileupcutin >> wf_ratio_low >> wf_ratio_high;
  }
  else {
    umsg.str("");
    umsg<<"Failed to Load waveform ratio Cut " << gatename;
    DANCE_Error("Unpacker",umsg.str());
    return -1;
  }
  pileupcutin.close();  


  //channel vector
  vector<int> channels;

  //Start of the unpacking process 
  gettimeofday(&tv,NULL);  
  double unpack_begin = tv.tv_sec+(tv.tv_usec/1000000.0);

  DANCE_Info("Unpacker","Started Unpacking");
  
  while (!gz_queue.empty()) { 
    //Stage 0 unpacking/stage 1 from midas
    if(input_params.Read_Binary==0 && input_params.Read_Simulation==0) {
 
      stringstream umsg;
      umsg.str("");
      umsg<<"Data Format: "<<input_params.DataFormat;
      DANCE_Info("Unpacker",umsg.str());
      
      while(run) {
        while(subrun) {
       
        if(analysis_params->entries_unpacked > EventLimit) {
          run=false;
          while (!gz_queue.empty()){gz_queue.pop();}
        }
        
        //Progess indicator
        if(analysis_params->entries_unpacked > progresscounter*ProgressInterval) {
          progresscounter++;
          cout<<"Processing Run Number: "<<input_params.RunNumber<<endl;
          if(datadeque.size()>0) {
            cout<<"Oldest Time in the Buffer: "<<datadeque[0].TOF<<endl;
            cout<<"Newest Time in the Buffer: "<<datadeque[datadeque.size()-1].TOF<<endl;
          }
          
          cout<<analysis_params->entries_unpacked<<" Entries Unpacked "<<endl;
          cout<<analysis_params->entries_awaiting_timesort<<" Entries Awaiting timesort"<<endl;
          cout<<datadeque.size()<<" Entries Sorted and in the Buffer"<<endl;
#ifdef CheckBufferDepth
          if(analysis_params->max_buffer_utilization < 0.75) {
            cout<<GREEN<<" Max Buffer Utilization: "<<100.0*analysis_params->max_buffer_utilization<<" %"<<RESET<<endl;
          }
          else if(analysis_params->max_buffer_utilization < 0.90) {
            cout<<YELLOW<<" Max Buffer Utilization: "<<100.0*analysis_params->max_buffer_utilization<<" %"<<RESET<<endl;
          }
          else {
            cout<<RED<<" Max Buffer Utilization: "<<100.0*analysis_params->max_buffer_utilization<<" %"<<RESET<<endl;
          }
#endif
          cout<<analysis_params->entries_written_to_binary<<" Entries Written to Binary"<<endl;
          cout<<analysis_params->entries_invalid<<" Entries Invalid ("<<100.0*analysis_params->entries_invalid/analysis_params->entries_processed<<" %)"<<endl;
          cout<<analysis_params->entries_built<<" Entries Built into "<<analysis_params->events_built<<" Events"<<endl;
          cout<<"Analyzed "<<analysis_params->entries_analyzed<<" Entries from "<<analysis_params->events_analyzed<<" Events"<<endl;
 
          cout<<setw(20)<<left<<"Breakdown:"<<setw(12)<<left<<"DANCE"<<setw(12)<<left<<"T0"<<setw(12)<<left<<"Li6"<<setw(12)<<left<<"U235"<<setw(12)<<left<<"He3"<<setw(12)<<left<<"Background"<<setw(12)<<left<<"Unknown"<<setw(12)<<left<<endl;
 
          cout<<setw(20)<<left<<"Entries:"<<setw(12)<<left<<analysis_params->DANCE_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->T0_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->Li6_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->U235_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->He3_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->Bkg_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->Unknown_entries<<endl;
 
          cout<<setw(20)<<left<<"Events:"<<setw(12)<<left<<analysis_params->DANCE_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->T0_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->Li6_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->U235_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->He3_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->Bkg_events_analyzed<<endl;
 
          if(analysis_params->DANCE_events_analyzed > 0) {
            cout<<setw(20)<<left<<"Average Mult:"<<setw(12)<<left<<(1.0*analysis_params->DANCE_entries_analyzed)/(1.0*analysis_params->DANCE_events_analyzed)<<endl;
          }
          cout<<endl;
 
          gettimeofday(&tv,NULL); 
          time_elapsed=tv.tv_sec+(tv.tv_usec/1000000.0);
        
          cout << "Average Entry Processing Rate: "<<(double)analysis_params->entries_unpacked/(time_elapsed-unpack_begin)<<" Entries per second "<<endl;
          cout << "Average Data Read Rate: "<<(double)TOTAL_BYTES/(time_elapsed-unpack_begin)/(1024.0*1024.0)<<" MB/s"<<endl;
          cout << "Instantaneous Data Read Rate: "<<(double)BYTES_READ/(time_elapsed-time_elapsed_old)/(1024.0*1024.0)<<" MB/s"<<endl;
          cout << (double)TOTAL_BYTES/(1024.0*1024.0*1024.0)<<" GiB Read"<<endl<<endl<<endl;
 
          BYTES_READ=0;
          time_elapsed_old = time_elapsed;
        } //end progress indicator
      
 
        //This is the unpacker for the caen2015 data format
        if(strcmp(input_params.DataFormat.c_str(),"caen2015") == 0) {
 
          //Read in the event header
          gzret=gzread(gz_in,&head,sizeof(EventHeader_t));
        
          //As long as gzread reads something start doing unpacking
          if(gzret!=0) {
          
            TotalDataSize = head.fDataSize;
            BYTES_READ += TotalDataSize;
            TOTAL_BYTES += TotalDataSize;
 
#ifdef Unpacker_Verbose
            cout<<"Type: "<<head.fEventId<<endl;
            cout<<"Size: "<<head.fDataSize<<endl;     ///< event size in bytes
            cout<<"TimeStamp "<<head.fTimeStamp<<endl;    ///< event timestamp in nseconds
#endif
           hEventID->Fill(head.fEventId);
           //cout << head.fTimeStamp << endl;
             //Data
            if(head.fEventId==1){
              gzret=gzread(gz_in,&bhead,sizeof(BankHeader_t));        
          
#ifdef Unpacker_Verbose
              cout << "Bank_HEADER " << endl;
              cout << dec <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
              cout << dec << bhead.fFlags << endl;
#endif
          
              TotalBankSize = bhead.fDataSize;
          
              while(TotalBankSize>0) {
                gzret=gzread(gz_in,&bank32,sizeof(Bank32_t));
                TotalBankSize-=sizeof(Bank32_t);
#ifdef Unpacker_Verbose
                cout << "BANK " << endl;
                cout << bank32.fName[0] << bank32.fName[1] << bank32.fName[2]<< bank32.fName[3] << endl;
                cout << dec << bank32.fType << endl;
                cout << dec << bank32.fDataSize << endl;
#endif
                EventBankSize = bank32.fDataSize;
                if (bank32.fName[0]=='C' && bank32.fName[1]=='E') {        // name starts as CE VT_BANK
                  int EVTS_batch=0;
              
                  evaggr->N = 0; // reset how many events we've processed this event
              
                  int number_cevt_events = bank32.fDataSize/sizeof(CEVT_BANK);
              
                  for (int eye = 0; eye < number_cevt_events; ++eye) {
                    gzret=gzread(gz_in,evinfo,sizeof(CEVT_BANK));
                    gzseek(gz_in,devt_padding,SEEK_CUR);
 
#ifdef Unpacker_Verbose 
                    cout<<"cevt event number: "<<eye<<endl;
                    cout<<"position: "<<evinfo->position<<endl;
                    cout<<"extras: "<<evinfo->extras<<endl;
                    cout<<"width: "<<evinfo->width<<endl;
                    cout<<"detector_id: "<<evinfo->detector_id<<endl;
                    cout<<evinfo->integral[0]<<"  "<<evinfo->integral[1]<<endl;
                    cout<<"padding: "<<devt_padding<<endl<<endl;;
#endif
                    TotalBankSize-=sizeof(CEVT_BANK)+devt_padding;
                    EventBankSize-=sizeof(CEVT_BANK)+devt_padding;          
                
                    evaggr->P[evaggr->N] = *evinfo;
                    evaggr->N++;
                    EVTS_batch += 1;
                
#ifdef Unpacker_Verbose 
                    cout << "evaggr->N: " << evaggr->N << endl;
#endif
                  }
                          
                  // snag the trig bank
                  gzret=gzread(gz_in,&bank32,sizeof(Bank32_t));
                  TotalBankSize-=sizeof(Bank32_t);
                  EventBankSize = bank32.fDataSize;
              
                  char *fData;
                  fData=(char*)malloc(bank32.fDataSize);
                  gzret=gzread(gz_in,fData,bank32.fDataSize);
                  TotalBankSize -= EventBankSize;
                  free (fData);        
              
#ifdef Unpacker_Verbose 
                  cout<<"Before trig bank"<<endl;
#endif
 
                  // begin funny place between peaks and cpu
                  while (true) {
                    // the peaks bank should be here
                    gzret=gzread(gz_in,&bank32,sizeof(Bank32_t));
                    TotalBankSize-=sizeof(Bank32_t);
                    EventBankSize = bank32.fDataSize;
 
#ifdef Unpacker_Verbose 
                    cout<<"Peak bank"<<endl;
                    cout<<"Name: "<<bank32.fName[0]<<"  Total Bank Size: "<<TotalBankSize<<"  EventBankSize: "<<EventBankSize<<endl;
#endif
                
                    if(bank32.fName[0]=='p') {
                      int whichpeak = atoi(&bank32.fName[1]);
#ifdef Unpacker_Verbose 
                      cout << "whichpeak: " << whichpeak << endl;
#endif
                      gzret=gzread(gz_in,imported_peaks[whichpeak],bank32.fDataSize);
#ifdef Unpacker_Verbose 
                      cout<<"After Read"<<endl;
#endif
                      //gzread(in,waveform,bank32.fDataSize);
                      TotalBankSize -= EventBankSize;
                    } 
                    else {
#ifdef Unpacker_Verbose 
                      cout<<"CPU Bank"<<endl;
#endif
                      // get the cpu bank information
                      char *fData=(char*)malloc(bank32.fDataSize);
                      gzret=gzread(gz_in,fData,bank32.fDataSize);
                      TotalBankSize -= EventBankSize;
                      free (fData);
                     // gz_queue.pop();
                     // subrun=false;
                     // run=false;
                      break; // you break here because the cpu comes last
                    }
                  }
#ifdef Unpacker_Verbose 
                  cout<<"After trig bank"<<endl;
#endif
 
                  int last_detnum = evaggr->P[0].detector_id;
                  int where_in_peakbank = 0;
                  for (uint32_t evtnum=0;evtnum<evaggr->N;++evtnum) {
                    int current_detnum = evaggr->P[evtnum].detector_id;
                    if (current_detnum != last_detnum) {
                      where_in_peakbank = 0;
                    }
                    uint32_t wflen = evaggr->P[evtnum].width;        // CEVT_BANK variable
                    analysis_params->wf_integral=0; 
                    for (uint wfindex=where_in_peakbank;wfindex<where_in_peakbank+wflen;++wfindex) {
                      // at this point we have reserved only 40 samples in db_arr waveform !!
                      evaggr->wavelets[evtnum][wfindex-where_in_peakbank] = imported_peaks[current_detnum][wfindex];
                    }        
                    where_in_peakbank += wflen;
                    last_detnum = current_detnum;
                
                    uint64_t timestamp_raw = (evaggr->P[evtnum].position & 0x7FFFFFFFFFFF);                // 47 bits for timestamp
 
#ifdef Unpacker_Verbose 
                    cout<<"timestamp_raw: "<<timestamp_raw<<endl;
#endif
                              
                    db_arr[EVTS].timestamp        = (double)(timestamp_raw);                                     //Digitizer timestamp
                    db_arr[EVTS].TOF               = (double)(timestamp_raw);                                     //Time of Flight (Currently in 2ns increments)
                    db_arr[EVTS].Ns                = evaggr->P[evtnum].width;                                     //Number of samples of the waveform
                    db_arr[EVTS].Ifast        = evaggr->P[evtnum].integral[0];                               //Fast integral
                    db_arr[EVTS].Islow        = evaggr->P[evtnum].integral[1]-evaggr->P[evtnum].integral[0]; //Slow integral
                    db_arr[EVTS].board        = (int)((1.*((int)evaggr->P[evtnum].detector_id)-1)/16.);      //Board number
                    db_arr[EVTS].channel        = (1*evaggr->P[evtnum].detector_id-1)-16*db_arr[EVTS].board;   //Channel number
                    db_arr[EVTS].ID             = MapID[db_arr[EVTS].channel][db_arr[EVTS].board];             //ID from DANCE map
                    db_arr[EVTS].Valid = 1;                                                                //Everything starts valid     
                    db_arr[EVTS].InvalidReason = 0;                                                            
 
 
#ifdef Unpacker_Verbose 
                    cout<<(int)db_arr[EVTS].board<<"  "<<(int)db_arr[EVTS].channel<<endl;
#endif
                    // CALCULATE THE LEADING EDGE using constant fraction "frac"
                    int imin=0;
                    double sigmin=1e9;
                    double frac=0.04;
                    double base=0;
                    double secmom=0.;        
                    int NNN=10;
                    // int id=db_arr[EVTS].ID;
 
                    double dT = 0;
                    
                    //Beam monitor waveforms for caen2015 data are ostensibly useless
                    if(db_arr[EVTS].ID <= 200) { 
                      if(db_arr[EVTS].ID<162) frac=0.04;
                      else frac=0.1;
                      
                      frac=0.2;
                      
                      for(int i=0;i<db_arr[EVTS].Ns;i++) {
                        wf1[i]=evaggr->wavelets[evtnum][i]+8192;
                        
                        if(i<NNN) {
                          base+=(1.*wf1[i]);
                          secmom+=(1.*wf1[i]*1.*wf1[i]);
                        }                  
                        if((1.*wf1[i])<sigmin) {
                          sigmin=1.*wf1[i];
                          imin=i;
                        }
                      }
                      
                      double thr=(sigmin-base/(1.*NNN))*frac+base/(1.*NNN);
                      int iLD=0;
                      for(int i=imin;i>1;i--){
                        if((1.*wf1[i])<thr && (1.*wf1[i-1])>thr){
                          double dSig=(1.*wf1[i-1]-1.*wf1[i]);
                          if(dSig!=0) dT=(1.*wf1[i-1]-thr)/dSig*2.+(i-1)*2.;  // this is in ns
                          else dT=(i-1)*2.;
                          iLD=i;
                        }                
                      }
                    } //end loop over ID <= 200
                    
                    //Fill raw IDs
 
#ifdef Histogram_Waveforms

                    //Fill waveform histogram
                    if(input_params.Read_Binary==0) {
                      if(db_arr[EVTS].ID<162) {
                        for(int i=0;i<db_arr[EVTS].Ns;i++) {
                          hWaveform_ID->Fill(i,wf1[i],db_arr[EVTS].ID);
                        }
                      } 
                     
                      if(db_arr[EVTS].ID == He3_ID) {
                        for(int i=0;i<db_arr[EVTS].Ns;i++) {
                          hWaveform_He3->Fill(i,wf1[i]);
                        }                          
                      } 
                      if(db_arr[EVTS].ID == Bkg_ID) {
                        for(int i=0;i<db_arr[EVTS].Ns;i++) {
                          hWaveform_Bkg->Fill(i,wf1[i]);
                        } 
                      }
                      if(db_arr[EVTS].ID == U235_ID) {
                        for(int i=0;i<db_arr[EVTS].Ns;i++) {
                          hWaveform_U235->Fill(i,wf1[i]);
                        }
                      } 
                      if(db_arr[EVTS].ID == Li6_ID) {
                        for(int i=0;i<db_arr[EVTS].Ns;i++) {
                          hWaveform_Li6->Fill(i,wf1[i]);
                        } 
                      }
                      if(db_arr[EVTS].ID == T0_ID) {
                        for(int i=0;i<db_arr[EVTS].Ns;i++) {
                          hWaveform_T0->Fill(i,wf1[i]);
                        } 
                      }
                    }

#endif
 
                    db_arr[EVTS].timestamp *= 2.0;                                                        //timestamp now in ns                 
                    db_arr[EVTS].timestamp += dT;                                                         //Full timestamp in ns
 
                  
                    if(input_params.Analysis_Stage > 0) {
                      //need to add the time deviations before time sorting
                      if(db_arr[EVTS].ID < 200) {
                        db_arr[EVTS].timestamp += TimeDeviations[db_arr[EVTS].ID];
                      }
                      
                      //Add the DANCE delay
                      if(db_arr[EVTS].ID < 162) {
                        db_arr[EVTS].timestamp += DANCE_Delay;
                      }
            
                      //Add the He3 delay
                      if(db_arr[EVTS].ID == He3_ID) {
                        db_arr[EVTS].timestamp += He3_Delay;
                      } 
                      
                      //Add the U235 delay
                      if(db_arr[EVTS].ID == U235_ID) {
                        db_arr[EVTS].timestamp += U235_Delay;
                      } 
                      
                      //Add the Li6 delay
                      if(db_arr[EVTS].ID == Li6_ID) {
                        db_arr[EVTS].timestamp += Li6_Delay;
                      }
                    }
 
                    db_arr[EVTS].TOF = db_arr[EVTS].timestamp;                                            //Start with TOF as Full timestamp in ns + time dev
 
                    
                    //keep track of the smallest timestamp
                    if(db_arr[EVTS].TOF<analysis_params->smallest_timestamp) {
                      analysis_params->smallest_timestamp=db_arr[EVTS].TOF;
                    }    
                    //keep track of the largest timestamp
                    if(db_arr[EVTS].TOF>analysis_params->largest_timestamp) {
                      analysis_params->largest_timestamp=db_arr[EVTS].TOF;
                    }  
                    
                    EVTS++;
            
                    analysis_params->entries_unpacked++;
                    analysis_params->entries_awaiting_timesort++;
 
                    
#ifdef Unpacker_Verbose 
                    cout<<EVTS<<"  "<<analysis_params->entries_unpacked<<endl;
#endif
                  }         //End of loop on eventnum                            
                }  //End of if on CEVT bank
                break;
              } //End of loop on EventBankSize
           
            } //End of event type 1
            
 
            
            //Scalers
            else if(head.fEventId==2) {
                
              TotalBankSize  = head.fDataSize;
 
#ifdef Scaler_Verbose
              cout << dec <<"TotalBankSize (bytes): " << TotalBankSize << endl;
#endif
              gzread(gz_in,&bhead,sizeof(BankHeader_t));
              TotalBankSize -= sizeof( BankHeader_t );
              
#ifdef Scaler_Verbose
              cout << "Bank_HEADER " << endl;
              cout << dec <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
              cout << dec << bhead.fFlags << endl;
#endif
              while (TotalBankSize > 0) {
 
                //MLTM
                gzread( gz_in, &bank32, sizeof( Bank32_t ) );
                TotalBankSize -= sizeof( Bank32_t );
 
#ifdef Scaler_Verbose
                cout << "BANK " << endl;
                cout << bank32.fName[0] << bank32.fName[1] << bank32.fName[2]<< bank32.fName[3] << endl;
                cout << dec << bank32.fType << endl;
                cout << "Size: "<<dec << bank32.fDataSize << endl;
#endif
                //see if data is on an 8-byte boundary
                bool readextra = false;
                if(bank32.fDataSize%8 !=0) {
                  readextra = true; 
                }
                
                if (bank32.fName[0]=='M' && bank32.fName[1]=='L' && bank32.fName[2]== 'T' && bank32.fName[3]=='M') {
                  
                  uint32_t time_seconds;
                  gzret=gzread(gz_in,&time_seconds,sizeof(time_seconds));
                  TotalBankSize-=sizeof(time_seconds);
 
#ifdef Scaler_Verbose
                  cout << time_seconds<<"\n";
#endif
                }
                if (bank32.fName[0]=='S' && bank32.fName[1]=='C' && bank32.fName[2]== 'L' && bank32.fName[3]=='R') {
 
                  gzret=gzread(gz_in,&sclr_totals,sizeof(sclr_totals));
                  TotalBankSize-=sizeof(sclr_totals);
          
                  for(int kay=0; kay<N_SCLR; kay++) {
                    hScalers->SetBinContent(kay+1,sclr_totals.totals[kay]);
#ifdef Scaler_Verbose
                    cout<<kay<<"  "<<sclr_totals.totals[kay]<<endl;
#endif
                  }
                }
                
                if (bank32.fName[0]=='R' && bank32.fName[1]=='A' && bank32.fName[2]== 'T' && bank32.fName[3]=='E') {
                 
                  gzret=gzread(gz_in,&sclr_rates,sizeof(sclr_rates));
                  TotalBankSize-=sizeof(sclr_rates);
 
#ifdef Scaler_Verbose
                  for(int kay=0; kay<N_SCLR; kay++) {
                    cout<<kay<<"  "<<sclr_rates.rates[kay]<<endl;
                  }
#endif
                }
                
                if(readextra) {
                  uint32_t extra = 0;
                  gzret=gzread(gz_in,&extra,sizeof(extra));
                  TotalBankSize -= sizeof(extra);
                }        
                
#ifdef Scaler_Verbose
                cout << dec <<"TotalBankSize (bytes): " << TotalBankSize << endl;
#endif
 
              } //End of while(TotalBankSize > 0)
            
            } //end of scalers
 
            else if(head.fEventId==0x8000 || head.fEventId==0x8001 || head.fEventId==0x8002 ){
 
              //End of Run
              if(head.fEventId==0x8001)  {
                run=false;
                break;
              }
          
              char *fData;
              fData=(char*)malloc(head.fDataSize);
              gzret=gzread(gz_in,fData,head.fDataSize);        
              free (fData);        
            }
            
            else { // other crap
              char *fData;
              fData=(char*)malloc(head.fDataSize);
              gzret=gzread(gz_in,fData,head.fDataSize);
              free (fData);
            } //end of crap
          }  //checking to see if gzret > 0
          /*else {
            gz_queue.pop();
          }*/
        } //end of caen2015
 
        else if(strcmp(input_params.DataFormat.c_str(),"caen2018") == 0) {
 
          //variables
          bool readextra = false;
          
          //counters
          uint32_t dataword =0;
          uint32_t BytesRead =0;
          uint32_t wordstoread = 0;
          uint32_t chaggcounter = 0;
          uint32_t chaggwordstoread = 0;
 
          //Start reading the file
          gzret=gzread(gz_in,&head,sizeof(EventHeader_t));
          
          if(gzret!=0) {
            
            TotalDataSize=head.fDataSize;
            
            BYTES_READ += head.fDataSize;
            TOTAL_BYTES += head.fDataSize;
            
#ifdef Unpacker_Verbose
            cout<<"Type: "<<head.fEventId<<"  TotalDataSize  "<<TotalDataSize<<endl;
#endif
           hEventID->Fill(head.fEventId);                          
            //Data
            if(head.fEventId==1){
              
              TotalBankSize=0;
              EventBankSize=0;
              
              gzret=gzread(gz_in,&bhead,sizeof(BankHeader_t));
#ifdef Unpacker_Verbose
              cout<<"Event Data"<<endl;
              cout << "Bank_HEADER " << endl;
              cout <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
              cout << bhead.fFlags << endl;
#endif
              
              TotalBankSize = bhead.fDataSize;
              
              while(TotalBankSize>0) {
                
                gzret=gzread(gz_in,&bank32,sizeof(Bank32_t));
                TotalBankSize -= sizeof(Bank32_t);
                
#ifdef Unpacker_Verbose
                cout<<"TotalBankSize after Bank Header Read "<<TotalBankSize<<endl;    
                cout << "BANK  " << bank32.fName[0] << bank32.fName[1] << bank32.fName[2]<< bank32.fName[3] << endl;
                cout << dec << bank32.fType << endl;
                cout << dec << bank32.fDataSize << endl;
#endif
          
                EventBankSize = bank32.fDataSize;
                
                //the data lie on 8 byte boundaries so there will be an extra 4 bytes at the end of the data that is "unaccounted" for in the header
                readextra = false;
                if(EventBankSize%8 !=0) {
                  readextra = true; 
                }
                
                //Read the firmware version and board ID
                gzret=gzread(gz_in,&dataword,sizeof(dataword));
                TotalBankSize -= sizeof(dataword);
                EventBankSize -=  sizeof(dataword);
                user_data.fw_majrev = (dataword & MAJREV_MASK);
                user_data.fw_minrev = (dataword & MINREV_MASK) >> 8;
                user_data.modtype = (dataword & MODTYPE_MASK) >> 14;
                user_data.modtype = 730;
                user_data.boardid = (dataword & BOARDID_MASK) >> 26;
                              
                //Read the user extras word
                gzret=gzread(gz_in,&dataword,sizeof(dataword));
                TotalBankSize -= sizeof(dataword);
                EventBankSize -=  sizeof(dataword);
 
                user_data.user_extra = dataword;
                
#ifdef Unpacker_Verbose
                cout<< "board: "<< (int)user_data.boardid <<" is a "<< (int)user_data.modtype <<" with Firmware: "<<(int)user_data.fw_majrev<< "."<<(int)user_data.fw_minrev<<"  "<<user_data.user_extra<<endl;
#endif
                
                while(EventBankSize>0) {
                  
                  BytesRead =0;
                  
                  //Make sure its a Vx725 or Vx730 
                  if(user_data.modtype == 725 || user_data.modtype == 730) {
 
                    //Read in the Vx725_Vx730 header
                    gzret=gzread(gz_in,&v1730_header,sizeof(v1730_header));
                    BytesRead += gzret;
 
                    //Unpack the header information
                    func_ret = unpack_vx725_vx730_board_data(&v1730_header, &vx725_vx730_board_data);
                    
#ifdef Unpacker_Verbose
                    cout<<"header: "<<(int)vx725_vx730_board_data.header<<"  ";
                    cout<<"nwords: "<<vx725_vx730_board_data.boardaggsize<<"  ";
                    cout<<"boardid: "<<(int)vx725_vx730_board_data.boardid<<"  ";
                    cout<<"pattern: "<<vx725_vx730_board_data.pattern<<"  ";
                    cout<<"channelmask: "<<(int)vx725_vx730_board_data.channelmask<<"  ";
                    cout<<"boardaggcounter: "<<vx725_vx730_board_data.boardaggcounter<<"  ";
                    cout<<"boardaggtime: "<<vx725_vx730_board_data.boardaggtime<<endl;
#endif
 
                    //Make sure the board header ID is 10 before proceeding
                    if(vx725_vx730_board_data.header != 10) {
                      cout<<RED<<"Unpacker [ERROR] CAEN Data Header is NOT 10!"<<endl;
                      cout<<"Entries: "<<EVTS<<" Total Entries: "<<analysis_params->entries_unpacked<<endl;
                      cout<<"Unpacker [ERROR] Data beyond this point would be corrupt and thus I am exiting to analysis!"<<RESET<<endl;
                      run = false;
                      return -1;
                      break;
                    }
                    
                    //interpret the channel mask 
                    channels.clear();
                    for(int m=0; m<8; m++) {
#ifdef Unpacker_Verbose
                      cout<<m<<"  "<<(( vx725_vx730_board_data.channelmask >> m) & 0x1)<<endl;
#endif
                      if((( vx725_vx730_board_data.channelmask >> m) & 0x1)) {
                        channels.push_back(2*m);
                      }
                    }
                    
                    //Number of words left in the CAEN data 
                    wordstoread = vx725_vx730_board_data.boardaggsize-4;
                    
                    //Number of channel aggregates unpacked
                    chaggcounter = 0;
                    
                    while (wordstoread>0) {
 
                      //Read the Vx725_Vx730 channel aggregate header
                      gzret=gzread(gz_in,&v1730_chagg_header,sizeof(v1730_chagg_header));
                      BytesRead += gzret;
                      wordstoread -= 2;            
                      
                      
                      //************  PSD ************//
                      if(user_data.fw_majrev == 136) {
                        
                        func_ret = unpack_vx725_vx730_psd_chagg_header(&v1730_chagg_header, &vx725_vx730_psd_data);
                        
#ifdef Unpacker_Verbose
                        cout<<"PSD   DT: "<<(int)vx725_vx730_psd_data.dual_trace<<"  ";
                        cout<<"EQ: "<<(int)vx725_vx730_psd_data.charge_enabled<<"  ";
                        cout<<"ET: "<<(int)vx725_vx730_psd_data.time_enabled<<"  ";
                        cout<<"EE: "<<(int)vx725_vx730_psd_data.extras_enabled<<"  ";
                        cout<<"ES: "<<(int)vx725_vx730_psd_data.waveform_enabled<<"  ";
                        cout<<"EX Opt: "<<(int)vx725_vx730_psd_data.extras_option<<"  ";
                        cout<<"AP: "<<(int)vx725_vx730_psd_data.ap<<"  ";
                        cout<<"DP1: "<<(int)vx725_vx730_psd_data.dp1<<"  ";
                        cout<<"DP2: "<<(int)vx725_vx730_psd_data.dp2<<"  ";
                        cout<<"NSDB8: "<<vx725_vx730_psd_data.nsdb8<<endl;
#endif
                        
                        //Number of words in the channel aggregate left to read
                        chaggwordstoread = vx725_vx730_psd_data.chagg_size - 2;
 
                        //Unpack the channel aggregate
                        while (chaggwordstoread>0) {
 
                          //Read the Vx725_Vx730 PSD channel aggregate
                          gzret=gzread(gz_in,&v1730_chagg_data,vx725_vx730_psd_data.individual_chagg_size*sizeof(uint32_t));
                          BytesRead += gzret;
                          wordstoread -= vx725_vx730_psd_data.individual_chagg_size;
                          chaggwordstoread -= vx725_vx730_psd_data.individual_chagg_size;
			 
                          //unpack the channel agregate
                          func_ret = unpack_vx725_vx730_psd_chagg(v1730_chagg_data, &vx725_vx730_psd_data);
 
                          //Set the remaining analysis variables
                          db_arr[EVTS].Valid = 1;                                                             //Everything starts valid
                          db_arr[EVTS].board = user_data.boardid;                                             //Board ID
                          db_arr[EVTS].channel = vx725_vx730_psd_data.channel + channels[chaggcounter];       //Channel ID
                          db_arr[EVTS].Ifast =  vx725_vx730_psd_data.qshort;                                  //Fast Integral
                          db_arr[EVTS].Islow =  vx725_vx730_psd_data.qlong - vx725_vx730_psd_data.qshort;     //Slow Integral (minus the fast)
                          db_arr[EVTS].InvalidReason = 0;                                                         
			 
                          //Map it
                          db_arr[EVTS].ID = MapID[db_arr[EVTS].channel][db_arr[EVTS].board];  
				
                               //Do waveform analysis and calculate times
                          if(vx725_vx730_psd_data.dual_trace) {
                            db_arr[EVTS].Ns = 4.0*vx725_vx730_psd_data.nsdb8;                                  //Dual trace effectively reduces the sampling frequency
                          }
                          else {
                            db_arr[EVTS].Ns = 8.0*vx725_vx730_psd_data.nsdb8;                                 
                          }
                          
                          double dT=0;
                          
                          analysis_params->wf_integral=0;
 
                          //If the detector is not a DANCE crystal or the use fine time is off
                          // NB: Calculate Frational time ALSO sets analysis_params->wf_integral
			  if ( ! input_params.Use_Firmware_FineTime || db_arr[EVTS].ID >= 162) {
                            dT = Calculate_Fractional_Time(vx725_vx730_psd_data.analog_probe1,                 //Function that calculates the fine time stamp
                                                             db_arr[EVTS].Ns, 
                                                           vx725_vx730_psd_data.dual_trace, 
                                                           user_data.modtype,
                                                           analysis_params);
                          }
                          else {
                            dT = 2.* vx725_vx730_psd_data.fine_time_stamp/1024.;
                          }
                           
                          //Set the timestamps
                          db_arr[EVTS].timestamp = vx725_vx730_psd_data.trigger_time_tag;                       //31-bit time in clock ticks
                          db_arr[EVTS].timestamp += 2147483648*vx725_vx730_psd_data.extended_time_stamp;        //16-bit extended time in clock ticks
                          db_arr[EVTS].timestamp *= 2.0;                                                        //timestamp now in ns                 
                          db_arr[EVTS].timestamp += dT;                                                         //Full timestamp in ns
 
                          if(analysis_params->wf_integral/(1.0*db_arr[EVTS].Islow) < wf_ratio_low || analysis_params->wf_integral/(1.0*db_arr[EVTS].Islow) > wf_ratio_high ) { 
                            db_arr[EVTS].pileup_detected=1;
                          } 
                          else {
                            db_arr[EVTS].pileup_detected=0;
                          }
			  db_arr[EVTS].wfintegral = analysis_params->wf_integral; //IK
                          
                          if(input_params.Analysis_Stage > 0) {
                            //need to add the time deviations before time sorting
                            if(db_arr[EVTS].ID < 200) {
                              db_arr[EVTS].timestamp += TimeDeviations[db_arr[EVTS].ID];
                            }
                            
                            //Add the DANCE delay
                            if(db_arr[EVTS].ID < 162) {
                              db_arr[EVTS].timestamp += DANCE_Delay;
                            }
                
                            //Add the He3 delay
                            if(db_arr[EVTS].ID == He3_ID) {
                              db_arr[EVTS].timestamp += He3_Delay;
                            } 
                            
                            //Add the U235 delay
                            if(db_arr[EVTS].ID == U235_ID) {
                              db_arr[EVTS].timestamp += U235_Delay;
                            } 
                            
                            //Add the Li6 delay
                            if(db_arr[EVTS].ID == Li6_ID) {
                              db_arr[EVTS].timestamp += Li6_Delay;
                            }
                          }
 
                          db_arr[EVTS].TOF = db_arr[EVTS].timestamp;                                            //TOF start as Full timestamp in ns
 
                        
                 
 
#ifdef Unpacker_Verbose
                          cout<<"Valid: "<<(int)db_arr[EVTS].Valid<<"  ";
                          cout<<"Board: "<<(int)db_arr[EVTS].board<<"  ";
                          cout<<"Channel: "<<(int)db_arr[EVTS].channel<<"  ";
                          cout<<"NS: "<<db_arr[EVTS].Ns<<"  ";
                          cout<<"Timestamp: "<<db_arr[EVTS].timestamp<<"  ";
                          cout<<"TOF: "<<db_arr[EVTS].TOF<<"  ";
                          cout<<"Ifast: "<<db_arr[EVTS].Ifast<<"  ";
                          cout<<"ISlow: "<<db_arr[EVTS].Islow<<endl;
#endif
 
#ifdef MakeTimeStampHistogram
                          if (db_arr[EVTS].ID<162){
                            hTimestamps->Fill(db_arr[EVTS].timestamp*1.0e-9);
                            hTimestampsID->Fill(db_arr[EVTS].timestamp*1.0e-9,db_arr[EVTS].ID);
                          }
                          if (db_arr[EVTS].ID==T0_ID){
                            hTimestampsT0->Fill(db_arr[EVTS].timestamp*1.0e-9);
                          }
                          if (db_arr[EVTS].ID==He3_ID || db_arr[EVTS].ID==Li6_ID || db_arr[EVTS].ID==U235_ID || db_arr[EVTS].ID==Bkg_ID ){
                            hTimestampsBM->Fill(db_arr[EVTS].timestamp*1.0e-9);
                          }
#endif                           
                                  
#ifdef Histogram_Digital_Probes
                          //Fill probe histograms
                          if(input_params.Read_Binary==0) {
                            if(db_arr[EVTS].ID<256) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                //digital probes
                                hDigital_Probe1_ID->Fill(kay,db_arr[EVTS].ID,vx725_vx730_psd_data.digital_probe1[kay]);
                                hDigital_Probe2_ID->Fill(kay,db_arr[EVTS].ID,vx725_vx730_psd_data.digital_probe2[kay]);
                              }
                            }
                          }
#endif
                          
#ifdef Histogram_Waveforms
 
                          //Fill waveform histograms
                          if(input_params.Read_Binary==0) {
                            if(db_arr[EVTS].ID<162) {
 
                          hID_vs_WFRatio->Fill(analysis_params->wf_integral/(1.0*db_arr[EVTS].Islow),db_arr[EVTS].ID);
 
                          hID_vs_WFInt_vs_Islow->Fill(analysis_params->wf_integral,db_arr[EVTS].Islow,db_arr[EVTS].ID);
 
			  hID_vs_WFRatio_vs_Islow->Fill(analysis_params->wf_integral/(1.0*db_arr[EVTS].Islow),db_arr[EVTS].Islow,db_arr[EVTS].ID);
                          //        cout<<db_arr[EVTS].ID<<"  "<<analysis_params->wf_integral/(1.0*db_arr[EVTS].Islow)<<endl;
 
                              if(waveform_counter < 20) {
                                if(db_arr[EVTS].Islow > 5000 && db_arr[EVTS].Ifast >500 && db_arr[EVTS].Ifast <1000) {
                                  for(int kay=0; kay<db_arr[EVTS].Ns; kay++) {
                                    hWaveforms[waveform_counter]->Fill(kay,vx725_vx730_psd_data.analog_probe1[kay]);
                                  }
                                  waveform_counter++;
                                }
                              }
                           
                              if(analysis_params->wf_integral<0) {
                                for(int kay=0; kay<db_arr[EVTS].Ns; kay++) {
                                  hWaveform_ID_NR->Fill(kay,vx725_vx730_psd_data.analog_probe1[kay],db_arr[EVTS].ID);
                                }                              
                              }
                              else {
                                for(int kay=0; kay<db_arr[EVTS].Ns; kay++) {
                                  hWaveform_ID->Fill(kay,vx725_vx730_psd_data.analog_probe1[kay],db_arr[EVTS].ID);
                                }
                              }
                              
                            }                        
     
                            if(db_arr[EVTS].ID == He3_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_He3->Fill(kay,vx725_vx730_psd_data.analog_probe1[kay]);
                              } 
                            }
                            if(db_arr[EVTS].ID == Bkg_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_Bkg->Fill(kay,vx725_vx730_psd_data.analog_probe1[kay]);
                              }
                            } 
                            if(db_arr[EVTS].ID == U235_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_U235->Fill(kay,vx725_vx730_psd_data.analog_probe1[kay]);
                              } 
                            }
                            if(db_arr[EVTS].ID == Li6_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_Li6->Fill(kay,vx725_vx730_psd_data.analog_probe1[kay]);
                              } 
                            }
                            
                            if(db_arr[EVTS].ID==T0_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_T0->Fill(kay,vx725_vx730_psd_data.analog_probe1[kay]);
                              }
                            }
                          }
#endif
 
          
                          //keep track of the smallest timestamp
                          if(db_arr[EVTS].TOF<analysis_params->smallest_timestamp) {
                            analysis_params->smallest_timestamp=db_arr[EVTS].TOF;
                          }    
                          //keep track of the largest timestamp
                          if(db_arr[EVTS].TOF>analysis_params->largest_timestamp) {
                            analysis_params->largest_timestamp=db_arr[EVTS].TOF;
                          }  
                      
                          EVTS++;
                  
                          analysis_params->entries_unpacked++;
                          analysis_params->entries_awaiting_timesort++;
 
                       
                          
#ifdef Unpacker_Verbose
                          cout<<"chaggwordstoread: "<<chaggwordstoread<<"  wordstoread: "<<wordstoread<<endl;
#endif
                        } //End of check on chagg words to read
                        
                        //increment the chagg counter
                        chaggcounter++;
 
                      } //End of check on PSD
                      
 
                      //*********** PHA ***********//
                      else if(user_data.fw_majrev == 139) {
                        
                        func_ret = unpack_vx725_vx730_pha_chagg_header(&v1730_chagg_header, &vx725_vx730_pha_data);
                        
#ifdef Unpacker_Verbose
                        cout<<"PHA   DT: "<<(int)vx725_vx730_pha_data.dual_trace<<"  ";
                        cout<<"EE: "<<(int)vx725_vx730_pha_data.energy_enabled<<"  ";
                        cout<<"ET: "<<(int)vx725_vx730_pha_data.time_enabled<<"  ";
                        cout<<"E2: "<<(int)vx725_vx730_pha_data.extras2_enabled<<"  ";
                        cout<<"ES: "<<(int)vx725_vx730_pha_data.waveform_enabled<<"  ";
                        cout<<"EX Opt: "<<(int)vx725_vx730_pha_data.extras2_option<<"  ";
                        cout<<"AP1: "<<(int)vx725_vx730_pha_data.ap1<<"  ";
                        cout<<"AP2: "<<(int)vx725_vx730_pha_data.ap2<<"  ";
                        cout<<"DP: "<<(int)vx725_vx730_pha_data.dp<<"  ";
                        cout<<"NSDB8: "<<vx725_vx730_pha_data.nsdb8<<endl;
#endif
                        
                        
                        //Number of words in the channel aggregate left to read
                        chaggwordstoread = vx725_vx730_pha_data.chagg_size - 2;
 
                        //Unpack the channel aggreate
                        while (chaggwordstoread>0) {
                          
                          //Read the Vx725_Vx730 PSD channel aggregate
                          gzret=gzread(gz_in,&v1730_chagg_data,vx725_vx730_pha_data.individual_chagg_size*sizeof(uint32_t));
                          BytesRead += gzret;
                          wordstoread -= vx725_vx730_pha_data.individual_chagg_size;
                          chaggwordstoread -= vx725_vx730_pha_data.individual_chagg_size;
 
                          //unpack the channel agregate
                          func_ret = unpack_vx725_vx730_pha_chagg(v1730_chagg_data, &vx725_vx730_pha_data);
                          
                          //Set the remaining analysis variables
                          db_arr[EVTS].Valid = 1;
                          db_arr[EVTS].board = user_data.boardid;
                          db_arr[EVTS].channel = vx725_vx730_pha_data.channel + channels[chaggcounter];
                          db_arr[EVTS].Ifast =  vx725_vx730_pha_data.energy;
                          db_arr[EVTS].Islow =  vx725_vx730_pha_data.energy;                
                          db_arr[EVTS].InvalidReason = 0;
 
                          //Map it
                          db_arr[EVTS].ID = MapID[db_arr[EVTS].channel][db_arr[EVTS].board]; 
 
                          
                          //Do waveform analysis and calculate times
                          if(vx725_vx730_pha_data.dual_trace) {
                            db_arr[EVTS].Ns = 4.0*vx725_vx730_pha_data.nsdb8;
                          }
                          else {
                            db_arr[EVTS].Ns = 8.0*vx725_vx730_pha_data.nsdb8;
                          }
                          
                          double dT=0;
                          if ( ! input_params.Use_Firmware_FineTime ) {
                            dT = Calculate_Fractional_Time(vx725_vx730_pha_data.analog_probe1,
                                                           db_arr[EVTS].Ns, 
                                                           vx725_vx730_pha_data.dual_trace, 
                                                           user_data.modtype,
                                                           analysis_params);
                          }
                          else {
                            dT = 2.*vx725_vx730_pha_data.fine_time_stamp/65356.;
                          }
                                                    
                          db_arr[EVTS].timestamp = vx725_vx730_pha_data.trigger_time_tag;                       //31-bit time in clock ticks
                          db_arr[EVTS].timestamp += 2147483648*vx725_vx730_pha_data.extended_time_stamp;        //16-bit extended time in clock ticks
                          db_arr[EVTS].timestamp *= 2.0;                                                        //timestamp now in ns                 
                          db_arr[EVTS].timestamp += dT;                                                         //Full timestamp in ns
                          
 
                          if(input_params.Analysis_Stage > 0) {
                            //need to add the time deviations before time sorting
                            if(db_arr[EVTS].ID < 200) {
                              db_arr[EVTS].timestamp += TimeDeviations[db_arr[EVTS].ID];
                            }
                            
                            //Add the DANCE delay
                            if(db_arr[EVTS].ID < 162) {
                              db_arr[EVTS].timestamp += DANCE_Delay;
                            }
                     
                            //Add the He3 delay
                            if(db_arr[EVTS].ID == He3_ID) {
                              db_arr[EVTS].timestamp += He3_Delay;
                            } 
                            
                            //Add the U235 delay
                            if(db_arr[EVTS].ID == U235_ID) {
                              db_arr[EVTS].timestamp += U235_Delay;
                            } 
                            
                            //Add the Li6 delay
                            if(db_arr[EVTS].ID == Li6_ID) {
                              db_arr[EVTS].timestamp += Li6_Delay;
                            }
                          }
 
                          db_arr[EVTS].TOF = db_arr[EVTS].timestamp;                                            //Start with TOF as Full timestamp in ns
 
#ifdef Unpacker_Verbose
                          cout<<"Valid: "<<(int)db_arr[EVTS].Valid<<"  ";
                          cout<<"Board: "<<(int)db_arr[EVTS].board<<"  ";
                          cout<<"Channel: "<<(int)db_arr[EVTS].channel<<"  ";
                          cout<<"NS: "<<db_arr[EVTS].Ns<<"  ";
                          cout<<"Timestamp: "<<db_arr[EVTS].timestamp<<"  ";
                          cout<<"TOF: "<<db_arr[EVTS].TOF<<"  ";
                          cout<<"Ifast: "<<db_arr[EVTS].Ifast<<"  ";
                          cout<<"ISlow: "<<db_arr[EVTS].Islow<<endl;
#endif
                          
 
          
#ifdef Histogram_Digital_Probes
                          //Fill probe histograms
                          if(input_params.Read_Binary==0) {
                            if(db_arr[EVTS].ID<256) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                //digital probes
                                hDigital_Probe1_ID->Fill(kay,db_arr[EVTS].ID,vx725_vx730_pha_data.digital_probe1[kay]);
                                hDigital_Probe2_ID->Fill(kay,db_arr[EVTS].ID,vx725_vx730_pha_data.digital_probe2[kay]);
                              }
                            }
                          }
#endif
                          
#ifdef Histogram_Waveforms
                          //Fill waveform histograms
                          if(input_params.Read_Binary==0) {
                            if(db_arr[EVTS].ID<162) {
                              for(int kay=0; kay<db_arr[EVTS].Ns; kay++) {
                                hWaveform_ID->Fill(kay,vx725_vx730_pha_data.analog_probe1[kay],db_arr[EVTS].ID);
                              }
                            }                        
                 
                            if(db_arr[EVTS].ID == He3_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_He3->Fill(kay,vx725_vx730_pha_data.analog_probe1[kay]);
                              } 
                            }
                            if(db_arr[EVTS].ID == Bkg_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_Bkg->Fill(kay,vx725_vx730_pha_data.analog_probe1[kay]);
                              }
                            } 
                            if(db_arr[EVTS].ID == U235_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_U235->Fill(kay,vx725_vx730_pha_data.analog_probe1[kay]);
                              } 
                            }
                            if(db_arr[EVTS].ID == Li6_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_Li6->Fill(kay,vx725_vx730_pha_data.analog_probe1[kay]);
                              } 
                            }
                            
                            if(db_arr[EVTS].ID == T0_ID) {
                              for(int kay=0;kay<db_arr[EVTS].Ns;kay++) {
                                hWaveform_T0->Fill(kay,vx725_vx730_pha_data.analog_probe1[kay]);
                              }
                            }
                          }
#endif
 
 
                          //keep track of the smallest timestamp
                          if(db_arr[EVTS].TOF<analysis_params->smallest_timestamp) {
                           analysis_params->smallest_timestamp=db_arr[EVTS].TOF;
                          }    
                          //keep track of the largest timestamp
                          if(db_arr[EVTS].TOF>analysis_params->largest_timestamp) {
                            analysis_params->largest_timestamp=db_arr[EVTS].TOF;
                          }  
                      
                          EVTS++;
            
                          analysis_params->entries_unpacked++;
                          analysis_params->entries_awaiting_timesort++;
 
 
                          
#ifdef Unpacker_Verbose
                          cout<<"chaggwordstoread: "<<chaggwordstoread<<"  wordstoread: "<<wordstoread<<endl;
#endif
                        } //End of check on chagg words to read
                        
                        //increment the chagg counter
                        chaggcounter++;
 
                      } //End of check on PHA
 
 
                    } //End of wordstoread > 0 
                    
#ifdef Unpacker_Verbose
                     cout<<EventBankSize<<"  "<<BytesRead<<endl;
#endif
    
                    EventBankSize -= BytesRead;
                    TotalBankSize -= BytesRead;
                    // BytesRead = 0;
                    
                  } //End of check on 725 or 730
                } //End of check on event bank size
 
                if(readextra) {
                  uint32_t extra = 0;
                  gzret=gzread(gz_in,&extra,sizeof(extra));
                  TotalBankSize -= sizeof(extra);
                } //End of read extra
              } //End of check  on total bank size
            }  //Endf of EventID 1 (Data)
 
            else if(head.fEventId==8){
              
              // this is scaler data
              gzret=gzread(gz_in,&bhead,sizeof(BankHeader_t));        
              
#ifdef Diagnostic_Verbose
              cout << "SCALER " << endl;
              cout << "Bank_HEADER " << endl;
              cout << dec <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
              cout << dec << bhead.fFlags << endl;
#endif    
            
              TotalBankSize = bhead.fDataSize;
 
              //This is the number of active boards
              int nactiveboards = 0;
              
              while(TotalBankSize > 0) {
 
                gzret=gzread(gz_in,&bank,sizeof(BankHeader_t));        
                TotalBankSize-=sizeof(Bank_t);
                
#ifdef Diagnostic_Verbose
                cout<<"TotalBankSize after Bank Header Read "<<TotalBankSize<<endl;
                cout << bank.fName[0] << bank.fName[1] << bank.fName[2]<< bank.fName[3] << endl;
                cout << dec << bank.fType << endl;
                cout << dec << bank.fDataSize << endl;
#endif
                //see if data is on an 8-byte boundary
                bool readextra = false;
                if(bank.fDataSize%8 !=0) {
                  readextra = true; 
                }
                
                //This is the time struct
                if (bank.fName[0]=='T' && bank.fName[1]=='I' && bank.fName[2]== 'M' && bank.fName[3]=='E') {
                  outputdiagnosticsfile << "TIME\n";
#ifdef Diagnostic_Verbose
                  cout<<endl<<"Time"<<endl;
#endif
                  gzret=gzread(gz_in,&timevalue,sizeof(timevalue));
                  TotalBankSize-=sizeof(timevalue);
 
                  outputdiagnosticsfile << timevalue.tv_sec<<"  "<<timevalue.tv_usec<<"\n";
                }
                
                //These are the Digitizer Rates
                if (bank.fName[0]=='S' && bank.fName[1]=='C' && bank.fName[2]== 'L' && bank.fName[3]=='R') {
                  
                  nactiveboards = bank.fDataSize/sizeof(uint32_t);
                  outputdiagnosticsfile << "SCLR  "<<nactiveboards<<"\n";
 
#ifdef Diagnostic_Verbose
                  cout<<endl<<"Digitizer Rates"<<endl;
#endif
                  nactiveboards = bank.fDataSize/sizeof(uint32_t);
                  for(int eye=0; eye<nactiveboards; eye++) {
                    gzret=gzread(gz_in,&Digitizer_Rates[eye],sizeof(uint32_t));
#ifdef Diagnostic_Verbose
                    cout<<eye<<"  "<<Digitizer_Rates[eye]<<endl;
#endif
                    TotalBankSize-=sizeof(uint32_t);
                    outputdiagnosticsfile << Digitizer_Rates[eye]<<"\n";
                  }
                }
 
                //These are the Acquisition Status
                if (bank.fName[0]=='A' && bank.fName[1]=='C' && bank.fName[2]== 'Q' && bank.fName[3]=='S') {
                  nactiveboards = bank.fDataSize/sizeof(uint32_t);
                  outputdiagnosticsfile << "ACQS  "<<nactiveboards<<"\n";
#ifdef Diagnostic_Verbose
                  cout<<endl<<"Acquisition Status"<<endl;
#endif
                  nactiveboards = bank.fDataSize/sizeof(uint32_t);
                  for(int eye=0; eye<nactiveboards; eye++) {
                    gzret=gzread(gz_in,&Acquisition_Status[eye],sizeof(uint32_t));
#ifdef Diagnostic_Verbose
                    cout<<eye<<"  "<<Acquisition_Status[eye]<<endl;
#endif
                    TotalBankSize-=sizeof(uint32_t);
                    outputdiagnosticsfile << Acquisition_Status[eye]<<"\n";
                  }
                }
 
                //These are the Failure Status
                if (bank.fName[0]=='F' && bank.fName[1]=='A' && bank.fName[2]== 'I' && bank.fName[3]=='L') {
                  nactiveboards = bank.fDataSize/sizeof(uint32_t);
                  outputdiagnosticsfile << "FAIL  "<<nactiveboards<<"\n";
#ifdef Diagnostic_Verbose
                  cout<<endl<<"Failure Status"<<endl;
#endif
                  nactiveboards = bank.fDataSize/sizeof(uint32_t);
                  for(int eye=0; eye<nactiveboards; eye++) {
                    gzret=gzread(gz_in,&Failure_Status[eye],sizeof(uint32_t));
 
                    if(Failure_Status[eye] != 0) {
                      faillog<<"Run: "<<input_params.RunNumber<<"  Board: "<<eye<<" Failure_Status: "<<Failure_Status[eye]<<endl;
                    }
#ifdef Diagnostic_Verbose
                    cout<<eye<<"  "<<Failure_Status[eye]<<endl;
#endif
                    TotalBankSize-=sizeof(uint32_t);
                    outputdiagnosticsfile << Failure_Status[eye]<<"\n";
                  }
                }
 
                //These are the Readout Status
                if (bank.fName[0]=='R' && bank.fName[1]=='E' && bank.fName[2]== 'A' && bank.fName[3]=='D') {
                  nactiveboards = bank.fDataSize/sizeof(uint32_t);
                  outputdiagnosticsfile << "READ  "<<nactiveboards<<"\n";
#ifdef Diagnostic_Verbose
                  cout<<endl<<"Readout Status"<<endl;
#endif
                  nactiveboards = bank.fDataSize/sizeof(uint32_t);
                  for(int eye=0; eye<nactiveboards; eye++) {
                    gzret=gzread(gz_in,&Readout_Status[eye],sizeof(uint32_t));
#ifdef Diagnostic_Verbose
                    cout<<eye<<"  "<<Readout_Status[eye]<<endl;
#endif
                    TotalBankSize-=sizeof(uint32_t);
                    outputdiagnosticsfile << Readout_Status[eye]<<"\n";
                  }
                }
 
                //These are the 8500 + 4n register values
                if (bank.fName[0]=='D' && bank.fName[1]=='I' && bank.fName[2]== 'A' && bank.fName[3]=='G') {
                  outputdiagnosticsfile << "DIAG  "<<nactiveboards<<"\n";
 
#ifdef Diagnostic_Verbose
                  cout<<endl<<"0x8500 + 4n Diagnostics"<<endl;
#endif
                  for(int eye=0; eye<nactiveboards; eye++) {
                    for(int jay=0; jay<8; jay++) {
                      gzret=gzread(gz_in,&Register_0x8504n[eye][jay],sizeof(uint32_t));
                      TotalBankSize-=sizeof(uint32_t);
                      outputdiagnosticsfile << Register_0x8504n[eye][jay]<<"  ";
                    }
                    outputdiagnosticsfile <<"\n";
 
                  }
                
#ifdef Diagnostic_Verbose
                  for(int eye=0; eye<8; eye++) {
                    for(int jay=0; jay<nactiveboards; jay++) {
                      cout<<Register_0x8504n[jay][eye]<<"  ";
                    }
                    cout<<endl;
                  }
#endif
                }
 
                //These are the ADC Temps
                if (bank.fName[0]=='T' && bank.fName[1]=='E' && bank.fName[2]== 'M' && bank.fName[3]=='P') {
                  outputdiagnosticsfile << "TEMP  "<<nactiveboards<<"\n";
#ifdef Diagnostic_Verbose
                  cout<<endl<<"ADC Temps"<<endl;
#endif
                  for(int eye=0; eye<nactiveboards; eye++) {
                    for(int jay=0; jay<16; jay++) {
                      gzret=gzread(gz_in,&ADC_Temp[eye][jay],sizeof(uint16_t));
                      TotalBankSize-=sizeof(uint16_t);
                      outputdiagnosticsfile << ADC_Temp[eye][jay]<<"  ";
                    }
                    outputdiagnosticsfile <<"\n";
                  }
                
#ifdef Diagnostic_Verbose
                  for(int eye=0; eye<16; eye++) {
                    for(int jay=0; jay<nactiveboards; jay++) {
                      cout<<ADC_Temp[jay][eye]<<"  ";
                    }
                    cout<<endl;
                  }
#endif
                }
 
                //These are the 0x1n2C values
                if (bank.fName[0]=='1' && bank.fName[1]=='n' && bank.fName[2]== '2' && bank.fName[3]=='C') {
                  outputdiagnosticsfile << "1n2C  "<<nactiveboards<<"\n";
#ifdef Diagnostic_Verbose
                  cout<<endl<<"Register 0x1n2C"<<endl;
#endif
                  for(int eye=0; eye<nactiveboards; eye++) {
                    for(int jay=0; jay<16; jay++) {
                      gzret=gzread(gz_in,&Register_0x1n2C[eye][jay],sizeof(uint32_t));
                      TotalBankSize-=sizeof(uint32_t);
                      outputdiagnosticsfile << Register_0x1n2C[eye][jay]<<"  ";
                    }
                    outputdiagnosticsfile <<"\n";
                  }
                
#ifdef Diagnostic_Verbose
                  for(int eye=0; eye<16; eye++) {
                    for(int jay=0; jay<nactiveboards; jay++) {
                      cout<<Register_0x1n2C[jay][eye]<<"  ";
                    }
                    cout<<endl;
                  }
#endif
                }
                
                //These are the Channel Status
                if (bank.fName[0]=='C' && bank.fName[1]=='H' && bank.fName[2]== 'S' && bank.fName[3]=='T') {
                  outputdiagnosticsfile << "CHST  "<<nactiveboards<<"\n";
#ifdef Diagnostic_Verbose
                  cout<<endl<<"Channel Status"<<endl;
#endif
                  for(int eye=0; eye<nactiveboards; eye++) {
                    for(int jay=0; jay<16; jay++) {
                      gzret=gzread(gz_in,&Channel_Status[eye][jay],sizeof(uint32_t));
                      TotalBankSize-=sizeof(uint32_t);
                      outputdiagnosticsfile << Channel_Status[eye][jay]<<"  ";
                    }
                    outputdiagnosticsfile <<"\n";
                  }
                
#ifdef Diagnostic_Verbose
                  for(int eye=0; eye<16; eye++) {
                    for(int jay=0; jay<nactiveboards; jay++) {
                      cout<<Channel_Status[jay][eye]<<"  ";
                    }
                    cout<<endl;
                  }
#endif
                }
                
                if(readextra) {
                  uint32_t extra = 0;
                  gzret=gzread(gz_in,&extra,sizeof(extra));
                  TotalBankSize -= sizeof(extra);
                }              
              }
#ifdef Diagnostic_Verbose
              cout<<"Done with Scalers. Total Bank Size: "<<TotalBankSize<<endl;
#endif
            }  //End of Event ID 8 (Diagnostics)
 
 
            //Scalers
            else if(head.fEventId==2) {
                
              TotalBankSize  = head.fDataSize;
 
#ifdef Scaler_Verbose
              cout << dec <<"TotalBankSize (bytes): " << TotalBankSize << endl;
#endif
              gzread(gz_in,&bhead,sizeof(BankHeader_t));
              TotalBankSize -= sizeof( BankHeader_t );
              
#ifdef Scaler_Verbose
              cout << "Bank_HEADER " << endl;
              cout << dec <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
              cout << dec << bhead.fFlags << endl;
#endif
              while (TotalBankSize > 0) {
 
                //MLTM
                gzread( gz_in, &bank32, sizeof( Bank32_t ) );
                TotalBankSize -= sizeof( Bank32_t );
 
#ifdef Scaler_Verbose
                cout << "BANK " << endl;
                cout << bank32.fName[0] << bank32.fName[1] << bank32.fName[2]<< bank32.fName[3] << endl;
                cout << dec << bank32.fType << endl;
                cout << "Size: "<<dec << bank32.fDataSize << endl;
#endif
                //see if data is on an 8-byte boundary
                bool readextra = false;
                if(bank32.fDataSize%8 !=0) {
                  readextra = true; 
                }
                
                if (bank32.fName[0]=='M' && bank32.fName[1]=='L' && bank32.fName[2]== 'T' && bank32.fName[3]=='M') {
                  
                  uint32_t time_seconds;
                  gzret=gzread(gz_in,&time_seconds,sizeof(time_seconds));
                  TotalBankSize-=sizeof(time_seconds);
 
#ifdef Scaler_Verbose
                  cout << time_seconds<<"\n";
#endif
                }
                if (bank32.fName[0]=='S' && bank32.fName[1]=='C' && bank32.fName[2]== 'L' && bank32.fName[3]=='R') {
 
                  gzret=gzread(gz_in,&sclr_totals,sizeof(sclr_totals));
                  TotalBankSize-=sizeof(sclr_totals);
          
                  for(int kay=0; kay<N_SCLR; kay++) {
                    hScalers->SetBinContent(kay+1,sclr_totals.totals[kay]);
#ifdef Scaler_Verbose
                    cout<<kay<<"  "<<sclr_totals.totals[kay]<<endl;
#endif
                  }
                }
                
                if (bank32.fName[0]=='R' && bank32.fName[1]=='A' && bank32.fName[2]== 'T' && bank32.fName[3]=='E') {
                 
                  gzret=gzread(gz_in,&sclr_rates,sizeof(sclr_rates));
                  TotalBankSize-=sizeof(sclr_rates);
 
#ifdef Scaler_Verbose
                  for(int kay=0; kay<N_SCLR; kay++) {
                    cout<<kay<<"  "<<sclr_rates.rates[kay]<<endl;
                  }
#endif
                }
                
                if(readextra) {
                  uint32_t extra = 0;
                  gzret=gzread(gz_in,&extra,sizeof(extra));
                  TotalBankSize -= sizeof(extra);
                }        
                
#ifdef Scaler_Verbose
                cout << dec <<"TotalBankSize (bytes): " << TotalBankSize << endl;
#endif
 
              } //End of while(TotalBankSize > 0)
            
            } // End of Event ID 2 (Scalers)
 
            // 0x8000 is a begin of run
            // 0x8001 is an end of run
            // 0x8002 is an ASCII message created by the logger 
          
            else if(head.fEventId==0x8000 || head.fEventId==0x8001 || head.fEventId==0x8002 || head.fEventId == 9 ) {
              //see if its the end of run
              if(head.fEventId==0x8001) {
                subrun=false;
                break;
              }
              
              char *fData;
              fData=(char*)malloc(head.fDataSize);
              gzret=gzread(gz_in,fData,head.fDataSize);
              free (fData);        
            }
                    
            else {
              // cout<<"EventID: "<<head.fEventId<<" Unknown"<<endl;
 
              TotalBankSize=TotalDataSize;
              
              uint32_t garbage;
              while(TotalBankSize > 0) {
                gzret = gzread(gz_in,&garbage, sizeof(garbage));
                TotalBankSize -= sizeof(garbage);
                //cout<<"Garbage Total Bank Size: "<<TotalBankSize<<endl;
              }            
 
            }  //End of catchall   
          }  //End of EventHeader gzret > 0
        } //End of caen2018 format check
        else {
          cout<<RED<<"Unpacker: [ERROR] I dont understand Data Format "<<input_params.DataFormat<<RESET<<endl;
          return -1;
        }
 
        //At this point we need to start ordering and eventbuilding
          if(EVTS >= BlockBufferSize) {
            
            //Sort this block of data
            func_ret = sort_array(db_arr,datadeque,EVTS,input_params,analysis_params);
            if(func_ret) {
              cout<<RED<<"Problem with sort_array in the MIDAS Reader"<<RESET<<endl;
              return -1;
            }
                    
            //Eventbuild
            func_ret = Build_Events(datadeque,input_params,analysis_params);
            if(func_ret) {
              cout<<RED<<"Problem with build_events in the MIDAS Reader"<<RESET<<endl;
              return -1;
            }
  
            //Reset the event counter and smallest timestamp
            EVTS=0;
            analysis_params->entries_awaiting_timesort=0;
            analysis_params->smallest_timestamp=2.814749767e14;
            
          } //end check on block buffer size and eventbuild
        }//end while subrun

        gz_queue.pop();
        if (gz_queue.size()>0) {
          //currently unused
          analysis_params->largest_subrun_timestamp=analysis_params->largest_timestamp;

          //grab the new subrun
          gz_in=gz_queue.front();
          input_params.SubRunNumber++;
          subrun=true;
        }
        else {
          run=false;
          break;
        }
      }  //end of while run


      cout<<"Unpacker [INFO]: Run Length: "<<analysis_params->largest_timestamp/1000000000.0<<" seconds"<<endl;
 
      //Now that we are done sorting we need to empty the buffer
      cout<<GREEN<<"Unpacker [INFO]: Finished Unpacking Data"<<RESET<<endl;
    
      //see if anything is left in the unsorted part
      if(EVTS>0) {
 
        umsg.str("");
        umsg<<"Unpacker [INFO]: There are "<<EVTS<<" Entries left to sort and "<<datadeque.size()<<" Entries left in the Buffer";
        DANCE_Info("Unpacker",umsg.str());
        
        //Sort this block of data
        func_ret = sort_array(db_arr,datadeque,EVTS,input_params,analysis_params);
        if(func_ret) {
          cout<<RED<<"Problem with sort_array in the empty stage of the MIDAS Reader"<<RESET<<endl;
          return -1;
        }
 
        EVTS=0;
        analysis_params->entries_awaiting_timesort=0;
 
      }
    
      if(datadeque.size()>0) {
 
        //need to set the buffer depth to zero
        input_params.Buffer_Depth = 0;
 
        //Eventbuild
        func_ret = Build_Events(datadeque,input_params,analysis_params);
        if(func_ret) {
          cout<<RED<<"Problem with build_events in the empty stage of the MIDAS Reader"<<RESET<<endl;
          return -1;
        }
 
        if(datadeque.size()==0) {
          cout<<GREEN<<"Unpacker [INFO]: Buffer empty, unpacking complete."<<RESET<<endl;
        }
      } //end check on timesort and datadeque
 
    } // END OF MIDAS READER
 
 
    //Stage 1 unpacker
    if(input_params.Read_Binary==1 || input_params.Read_Simulation==1) {
 
      while(run) {
        
        //Event limit control
        if(analysis_params->entries_unpacked > EventLimit) {
          run=false;
          gz_queue.pop();
        }
        
        //Progress indicator
        if(analysis_params->entries_unpacked > progresscounter*ProgressInterval) {
          progresscounter++;
          if(input_params.Read_Simulation == 0) {
            cout<<"Processing Run Number: "<<input_params.RunNumber<<endl;
          }
          else {
            cout<<"Processing Simulated Data"<<endl;
          }
 
          cout<<analysis_params->entries_unpacked<<" Entries Unpacked "<<endl;
          cout<<analysis_params->entries_awaiting_timesort<<" Entries Awaiting timesort"<<endl;
          cout<<datadeque.size()<<" Entries Sorted and in the Buffer"<<endl;
#ifdef CheckBufferDepth
          if(analysis_params->max_buffer_utilization < 0.75) {
            cout<<GREEN<<" Max Buffer Utilization: "<<100.0*analysis_params->max_buffer_utilization<<" %"<<RESET<<endl;
          }
          else if(analysis_params->max_buffer_utilization < 0.90) {
            cout<<YELLOW<<" Max Buffer Utilization: "<<100.0*analysis_params->max_buffer_utilization<<" %"<<RESET<<endl;
          }
          else {
            cout<<RED<<" Max Buffer Utilization: "<<100.0*analysis_params->max_buffer_utilization<<" %"<<RESET<<endl;
          }
#endif
          cout<<analysis_params->entries_written_to_binary<<" Entries Written to Binary"<<endl;
          cout<<analysis_params->entries_invalid<<" Entries Invalid ("<<100.0*analysis_params->entries_invalid/analysis_params->entries_processed<<" %)"<<endl;
          cout<<analysis_params->entries_built<<" Entries Built into "<<analysis_params->events_built<<" Events"<<endl;
          cout<<"Analyzed "<<analysis_params->entries_analyzed<<" Entries from "<<analysis_params->events_analyzed<<" Events"<<endl;
 
          cout<<setw(20)<<left<<"Breakdown:"<<setw(12)<<left<<"DANCE"<<setw(12)<<left<<"T0"<<setw(12)<<left<<"Li6"<<setw(12)<<left<<"U235"<<setw(12)<<left<<"He3"<<setw(12)<<left<<"Background"<<setw(12)<<left<<"Unknown"<<setw(12)<<left<<endl;
 
          cout<<setw(20)<<left<<"Entries:"<<setw(12)<<left<<analysis_params->DANCE_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->T0_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->Li6_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->U235_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->He3_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->Bkg_entries_analyzed;
          cout<<setw(12)<<left<<analysis_params->Unknown_entries<<endl;
 
          cout<<setw(20)<<left<<"Events:"<<setw(12)<<left<<analysis_params->DANCE_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->T0_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->Li6_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->U235_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->He3_events_analyzed;
          cout<<setw(12)<<left<<analysis_params->Bkg_events_analyzed<<endl;
 
          if(analysis_params->DANCE_events_analyzed > 0) {
            cout<<setw(20)<<left<<"Average Mult:"<<setw(12)<<left<<(1.0*analysis_params->DANCE_entries_analyzed)/(1.0*analysis_params->DANCE_events_analyzed)<<endl;
          }
          cout<<endl;
 
          gettimeofday(&tv,NULL); 
          time_elapsed=tv.tv_sec+(tv.tv_usec/1000000.0);
        
          cout << "Average Entry Processing Rate: "<<(double)analysis_params->entries_unpacked/(time_elapsed-unpack_begin)<<" Entries per second "<<endl;
          cout << "Average Data Read Rate: "<<(double)TOTAL_BYTES/(time_elapsed-unpack_begin)/(1024.0*1024.0)<<" MB/s"<<endl;
          cout << "Instantaneous Data Read Rate: "<<(double)BYTES_READ/(time_elapsed-time_elapsed_old)/(1024.0*1024.0)<<" MB/s"<<endl;
          cout << (double)TOTAL_BYTES/(1024.0*1024.0*1024.0)<<" GiB Read"<<endl<<endl<<endl;
          
          BYTES_READ=0;
          time_elapsed_old = time_elapsed;
        }
        
	//Fill the array
	if(input_params.WF_Integral){    //if WF integral read from binaries, fill DEVT_STAGE1_WF struct
       	  gzret=gzread(gz_in,&devt_stage1_wf,sizeof(DEVT_STAGE1_WF));
 	  db_arr[EVTS].timestamp = devt_stage1_wf.timestamp;
          db_arr[EVTS].wfintegral = devt_stage1_wf.wfintegral;
          db_arr[EVTS].Ifast = devt_stage1_wf.Ifast;
          db_arr[EVTS].Islow = devt_stage1_wf.Islow;
          db_arr[EVTS].ID = devt_stage1_wf.ID;
	  // We need to actually check the wf integral to see if pileup
          if ( devt_stage1_wf.wfintegral/(1.0*db_arr[EVTS].Islow) < wf_ratio_low 
                || devt_stage1_wf.wfintegral/(1.0*db_arr[EVTS].Islow) > wf_ratio_high ) 
          { 
            db_arr[EVTS].pileup_detected=1;                                                         
          } 
          else {
            db_arr[EVTS].pileup_detected=0;                                                         
          }
	}
	else {
	  gzret=gzread(gz_in,&devt_stage1,sizeof(DEVT_STAGE1));  //if no WF integral read from binaries, fill DEVT_STAGE1 struct
	  db_arr[EVTS].timestamp = devt_stage1.timestamp;
          db_arr[EVTS].Ifast = devt_stage1.Ifast;
          db_arr[EVTS].Islow = devt_stage1.Islow;
          db_arr[EVTS].ID = devt_stage1.ID;
	}

        // gzseek(gz_in,devt_padding,SEEK_CUR);
        
        if(gzret!=0) {
          
          BYTES_READ += gzret;
          TOTAL_BYTES += gzret;
          
          //Fill the array
          db_arr[EVTS].Valid = 1; //Everything starts valid
          db_arr[EVTS].InvalidReason = 0; //Everything starts valid
 
          //Time Deviations
          if(input_params.Analysis_Stage > 0) {
            //need to add the time deviations before time sorting
            if(db_arr[EVTS].ID < 200) {
              db_arr[EVTS].timestamp += TimeDeviations[db_arr[EVTS].ID];
            }
            
            //Add the DANCE delay
            if(db_arr[EVTS].ID < 162) {
              db_arr[EVTS].timestamp += DANCE_Delay;
            }
 
            //Add the He3 delay
            if(db_arr[EVTS].ID == He3_ID) {
              db_arr[EVTS].timestamp += He3_Delay;
            } 
            
            //Add the U235 delay
            if(db_arr[EVTS].ID == U235_ID) {
              db_arr[EVTS].timestamp += U235_Delay;
            } 
            
            //Add the Li6 delay
            if(db_arr[EVTS].ID == Li6_ID) {
              db_arr[EVTS].timestamp += Li6_Delay;
            }
          }
 
          db_arr[EVTS].TOF = db_arr[EVTS].timestamp;                                            //Start with TOF as Full timestamp in ns
 
 
          //keep track of the smallest timestamp
          if(db_arr[EVTS].TOF<analysis_params->smallest_timestamp) {
            analysis_params->smallest_timestamp=db_arr[EVTS].TOF;
          }    
          //keep track of the largest timestamp
          if(db_arr[EVTS].TOF>analysis_params->largest_timestamp) {
            analysis_params->largest_timestamp=db_arr[EVTS].TOF;
          }  
            
          EVTS++;           
 
          analysis_params->entries_unpacked++;
          analysis_params->entries_awaiting_timesort++;
 
        }
        else {
          run=false;
          gz_queue.pop();
          break;
        }
 
        //At this point we need to start ordering and eventbuilding
        if(EVTS >= BlockBufferSize) {
          
          //sort this block of data
          func_ret = sort_array(db_arr,datadeque,EVTS,input_params,analysis_params);
          if(func_ret) {
            cout<<RED<<"Problem with sort_array in the Binary Reader"<<RESET<<endl;
            return -1;
          }
 
          //Eventbuild
          func_ret = Build_Events(datadeque,input_params, analysis_params);
          if(func_ret) {
            cout<<RED<<"Problem with build_events in the Binary Reader"<<RESET<<endl;
            return -1;
          }
          
          //Reset the event counter and smallest timestamp
          EVTS=0;
          analysis_params->entries_awaiting_timesort=0;
          analysis_params->smallest_timestamp=2.814749767e14;
                  
        }        
      }  //end of while(run)
 
      umsg.str("");
      umsg<<"Run Length: "<<analysis_params->largest_timestamp/1000000000.0<<" seconds";
      DANCE_Info("Unpacker",umsg.str());
 
      //Now that we are done sorting we need to empty the buffer
      DANCE_Info("Unpacker","Finished unpacking data");
    
      //see if anything is left in the unsorted part
      if(EVTS>0) {
        umsg.str("");
        umsg<<"There are "<<EVTS<<" Entries left to sort and "<<datadeque.size()<<" Entries left in the Buffer";
        DANCE_Info("Unpacker",umsg.str());
 
        //Sort this block of data
        func_ret = sort_array(db_arr,datadeque,EVTS,input_params,analysis_params);
        if(func_ret) {
          DANCE_Error("Unpacker","Problem with sort_array in the empty stage of the Binary Reader");
          return -1;
        }
        
        EVTS=0;
        analysis_params->entries_awaiting_timesort=0;
      }
    
      if(datadeque.size()>0) {
 
        //need to set the buffer depth to zero
        input_params.Buffer_Depth = 0;
 
        //Eventbuild
        func_ret = Build_Events(datadeque,input_params, analysis_params);
        if(func_ret) {
          DANCE_Error("Unpacker","Problem with build_events in the empty stage of the Binary Reader");
          return -1;
        }
         
        if(datadeque.size()==0) {
          DANCE_Success("Unpacker","Buffer empty, unpacking complete.");
        }
      }
    } //end of if(input_params.Read_Binary==1)


    func_ret = Write_Root_File(input_params, analysis_params);

    if (func_ret) {
      return -1;
    }
    if (gz_queue.size()==1){gz_queue.pop();} 
  } 
  //Make the time deviations if needed (Likely only a stage 0 thing)
  if(input_params.FitTimeDev) {
    Make_Time_Deviations(input_params.RunNumber);
  }



 
  //Unpacking is finished.  Return the total number of events unpacked and analyzed
  return analysis_params->entries_unpacked;
 
 
}

int Initialize_Unpacker(Input_Parameters input_params) {  

  DANCE_Init("Unpacker","Initializing");

  int func_ret = 0;

  cout<<"Buffer Depth: "<<input_params.Buffer_Depth<<" seconds"<<endl;
  cout<<"Coincidence Window: "<<input_params.Coincidence_Window<<" ns"<<endl;
  cout<<"Crystal Blocking Time: "<<input_params.Crystal_Blocking_Time<<" ns"<<endl;
  cout<<"DANCE Event Blocking Time: "<<input_params.DEvent_Blocking_Time<<" ns"<<endl;
  
  if(input_params.HAVE_Threshold) {
    cout<<"Energy Thresholds are ON and set to: "<<input_params.Energy_Threshold<<" MeV"<<endl;
  }
  else {
    cout<<"Energy Thresholds are OFF"<<endl;
  }
  if(input_params.FitTimeDev) {
    cout<<"Time Deviations will be determined following analysis"<<endl;
  }
  cout<<endl;
 
  //initialize histograms
  func_ret += Create_Unpacker_Histograms(input_params);
  
  //initialize DANCE Map
  func_ret +=  Make_DANCE_Map();

  //initiliaze the time deviations
  func_ret += Read_TimeDeviations(input_params);

  //Make the ouput diagnostics file
  if(input_params.Analysis_Stage == 0 && (strcmp(input_params.DataFormat.c_str(),"caen2018") == 0) && input_params.Read_Simulation==0) {
    func_ret += Make_Output_Diagnostics_File(input_params.RunNumber);
  }
 
  if(func_ret==0) {
    DANCE_Success("Unpacker","Initialized");
  }
  else {
    DANCE_Error("Unpacker","Initialization Failed. Exiting!");
  }
  return func_ret;

}
