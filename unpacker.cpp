//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  unpacker.cpp           *// 
//*  Last Edit: 03/22/18    *//  
//***************************//

//File includes
#include "global.h"
#include "unpacker.h"
#include "structures.h"
#include "sort_functions.h"
#include "analyzer.h"

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

//ROOT includes
#include "TH1.h"
#include "TH3.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom.h"

using namespace std;

//Some Global Unpacker Variables (DONT CHANGE UNLESS NEEDED)

/*size of the DEVT array in the unpacker*/
#define MaxDEVTArrSize 100000000  //this should be a big number like 1e8 but pay attention

/*size of the block buffer. The unpacker waits this many
  events before looking at the deque and time sorting.  This 
  number has huge efficiency implications so dont change it if
  you dont understand it */
#define BlockBufferSize 350000

/*Time depth (in seconds) of the buffer.  While you cant know if you 
  will fail the program will tell you if you did... */
#define BufferDepth 10000  //Read the whole thing in.  DANCE doesnt always flush till the end anyway. 


//Verbosity and Error Checking
//#define CheckTheDeque
//#define Eventbuilder_Verbose
//#define Unpacker_Verbose
#define Scaler_Verbose

//output diagnostics file
ofstream outputdiagnosticsfile;

//global unpacker variables
double TimeDeviations[200];

//holds the id of each detector from dancemap
int MapID[20][20];

//Histograms
TH3S *hWaveform_ID;
TH1D *hID_Raw;
//TH2S *hWaveform_T0;

//Histograms for Unpacker Things
int Create_Unpacker_Histograms(bool read_binary) {
  
  cout<<"Unpacker [INFO]: Creating Histograms"<<endl;
  
  //Make a histogram for waveforms
  if(read_binary==0) {
    hWaveform_ID = new TH3S("Waveform_ID","Waveform_ID",80,0,80,2000,0,20000,256,0,256);
    hID_Raw = new TH1D("hID_Raw","hID_Raw",256,0,256);
    //  hWaveform_T0 = new TH2S("Waveform_T0","Waveform_T0",80,0,80,2000,0,20000);
  }
  
  cout<<GREEN<<"Unpacker [INFO]: Created Histograms"<<RESET<<endl;
  return 0; 
  
}

int Write_Unpacker_Histograms(TFile *fout, bool read_binary) {
  
  cout<<"Unpacker [INFO]: Writing Histograms"<<endl;
  
  fout->cd();

  if(read_binary==0) {
    hWaveform_ID->Write();
    hID_Raw->Write();
    // hWaveform_T0->Write();
  }
  cout<<GREEN<<"Unpacker [INFO]: Wrote Histograms"<<RESET<<endl;
  return 0;
}

//Functions
int Make_DANCE_Map() {
  cout << "Unpacker [INFO]: Reading "<<DanceMapFile<<endl;
  
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
    cout<<GREEN<<"Unpacker [INFO]: Made DANCE Map"<<RESET<<endl;
    return 0;
  }
  else {
    cout<<RED<<"Unpacker [ERROR]: Could Not Open DANCE Map File"<<endl;
    return -1;
  }  
}

//This function makes the output binary file for stage 0 or 1 filled with time-ordered devt_bank structures
int Make_Output_Diagnostics_File(int RunNumber) {
  
  stringstream outfilename;
  outfilename.str();
  
  outfilename << DIAGNOSTICS;
  outfilename <<"/diagnostics_run";
  outfilename << RunNumber << ".txt";
  
  outputdiagnosticsfile.open(outfilename.str().c_str(), ios::out);
  
  if(outputdiagnosticsfile.is_open()) {
    cout<<GREEN<<"Analyzer [INFO]: Succesfully created and opened output diagnostics file: "<<outfilename.str()<<RESET<<endl;
    return 0;
  }
  else {
    cout<<RED<<"Analyzer [ERROR]: Failed to create output diagnostics file: "<<outfilename.str()<<RESET<<endl;
    return -1;
  }
}

int Read_TimeDeviations(int runnum, bool FitTimeDev) {
  cout << "Unpacker [INFO]: Reading Time Deviations" << endl;
  
  for(int eye=0; eye<200; eye++) {
    TimeDeviations[eye]=0;
  }  
  
  //If we want to fit the time deviations then set all of them to zero
  if(FitTimeDev) {
    cout<<GREEN<<"Unpacker [INFO]: Set all time deviations to 0 "<<RESET<<endl;
  }
  
  //if we have time deviations then obtain them from the files
  else {
    stringstream fname;
    fname.str();
    
    fname<<TIMEDEV_DIR"/TimeDeviations_Run_" << std::setfill('0') << std::setw(5) << runnum << ".txt";
    cout<<"Unpacker [INFO]: Looking for TimeDeviations File: "<<fname.str()<<endl;
    
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
      cout<<GREEN<<"Unpacker [INFO]: Set all time deviations to values in "<<fname.str()<<RESET<<endl;
      return 0;
    }
    else {
      cout<<RED<<"Unpacker [ERROR]: Couldnt find "<<fname.str()<<" Setting all time deviations to 0"<<RESET<<endl;
    }
  }
  return 0;
}


int Unpack_Data(gzFile &gz_in, double begin, int runnum, bool read_binary, bool write_binary, double CoincidenceWindow, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold, bool FitTimeDev, string DataFormat, int NQGates, double QGates[]) {

  ofstream faillog;
  faillog.open("Readout_Status_Failures.txt", ios::app);

  cout<<BLUE<<"Unpacker [INIT]: Initializing Unpacker"<<RESET<<endl<<endl;

  cout<<"Buffer Depth: "<<BufferDepth<<" seconds"<<endl;
  cout<<"Coincidence Window: "<<CoincidenceWindow<<" ns"<<endl;
  cout<<"Crystal Blocking Time: "<<Crystal_Blocking_Time<<" ns"<<endl;
  cout<<"DANCE Event Blocking Time: "<<DEvent_Blocking_Time<<" ns"<<endl;
  
  if(HAVE_Threshold) {
    cout<<"Energy Thresholds are ON and set to: "<<Energy_Threshold<<" MeV"<<endl;
  }
  else {
    cout<<"Energy Thresholds are OFF"<<endl;
  }
  if(FitTimeDev) {
    cout<<"Time Deviations will be determined following analysis"<<endl;
  }
  cout<<endl;

  //initialize histograms
  Create_Unpacker_Histograms(read_binary);
  
  //initialize DANCE Map
  Make_DANCE_Map();

  //initiliaze the time deviations
  Read_TimeDeviations(runnum,FitTimeDev);

  //Make the ouput diagnostics file
  if(read_binary == 0 && (strcmp(DataFormat.c_str(),"caen2018") == 0)) {
    Make_Output_Diagnostics_File(runnum);
  }
  
  //Make the output binary file
  if(write_binary == 1) {
    Make_Output_Binfile(runnum, read_binary);
  }
  
  //keep track of timestamps
  double last_timestamp[2000];
  for(int eye=0; eye<2000; eye++) {
    last_timestamp[eye]=0;
  }
  
  //define some things for the unpacker
  double buffer_length = (double)1000000000.0*BufferDepth; //clock ticks

  bool run=true; 
  bool first_sort=true;
  bool event_building_active=false;  //this says whether or not we are event building yet
  bool emptythedeque=true; //once done clean out the deque
  bool timesort=true;

  //Structures to put data in
  deque<DEVT_BANK> datadeque;                         //Storage container for time sorted data
  vector<DEVT_BANK> eventvector;                      //Vector to store events for analysis
  CEVT_BANK *evinfo = new CEVT_BANK();                //caen event info
  test_struct_cevt *evaggr = new test_struct_cevt();  //event aggregate
  DEVT_BANK *db_arr = new DEVT_BANK[MaxDEVTArrSize];  //Storage array for entries
  DEVT_STAGE1 devt_stage1;                            //Stage1 format for reading binary
  
  //Counters
  uint32_t TOTAL_EVTS=0;        //Total number of entries unpacked
  uint32_t EVTS=0;              //Total number of entries unpacked since last time sort
  uint64_t BYTES_READ=0;        //Total number of Bytes read since last time sort
  uint64_t TOTAL_BYTES=0;       //Total number of Bytes read 
  uint32_t Events_Sent=0;       //Total number of entries sent to the analyzer
  uint32_t progresscounter=1;   //Keep track of how many progress statements have been made
  int gzret=1;                  //number of bytes read by gzread
  
  double smallest_timestamp=400000000000000.0;  //This just needs to be a large number
  
  //Some CAEN structures to optimize reads
  V1730_Header_t v1730_header;
  V1730_ChAgg_Header_t v1730_chagg_header;

  //MIDAS Bank Stuff
  EventHeader_t head;           //MIDAS event header
  BankHeader_t bhead;           //MIDAS bank header
  Bank32_t bank32;              //MIDAS 32-bit bank
  Bank_t bank;                  //MIDAS 16-bit bank

  uint32_t TotalDataSize=0;     // head.fDataSize;
  uint32_t TotalBankSize=0;     // bhead.fDataSize;
  uint32_t EventBankSize=0;     // bank32fDataSize;

  long devt_padding = 0;        // padding between banks not divisible by 64 bits
  
  short imported_peaks[256][16384]; // this is actually supported channels / supported length of PXXX bank
  unsigned short int wf1[15000];

  //time profiling for performance
  struct timeval tv;   	   //Real time  
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

  //waveform
  // vector<uint16_t> waveform;
  uint16_t waveform[1000];

  //channel vector
  vector<int> channels;

  //Start of the unpacking process 
  gettimeofday(&tv,NULL);  
  double unpack_begin = tv.tv_sec+(tv.tv_usec/1000000.0);

  cout<<GREEN<<"Unpacker [INFO]: Started Unpacking: "<<RESET<<endl;
  cout<<"Unpacker [INFO]: Data Format: "<<DataFormat<<endl;
  
  //Stage 0 unpacking
  if(read_binary==0) {
    while(run) {
     
      if(TOTAL_EVTS > EventLimit) {
	run=false;
      }
      
      //Progess indicator
      if(TOTAL_EVTS > progresscounter*ProgressInterval) {
	progresscounter++;
	cout<<"Processing Run Number: "<<runnum<<endl;
	if(timesort && datadeque.size()>0) {
	  cout<<"Oldest Time in the Buffer: "<<datadeque[0].TOF<<endl;
	  cout<<"Newest Time in the Buffer: "<<datadeque[datadeque.size()-1].TOF<<endl;
	}	
	cout<<TOTAL_EVTS<<" Events Unpacked "<<endl;
	cout<<Events_Sent<<" Events Analyzed"<<endl;
	cout<<datadeque.size()<<" Events in the Buffer"<<endl;
	gettimeofday(&tv,NULL); 
	time_elapsed=tv.tv_sec+(tv.tv_usec/1000000.0);
      
	cout << "Average Event Processing Rate: "<<(double)TOTAL_EVTS/(time_elapsed-unpack_begin)<<" Events per second "<<endl;
	cout << "Average Data Read Rate: "<<(double)TOTAL_BYTES/(time_elapsed-unpack_begin)/(1024.0*1024.0)<<" MB/s"<<endl;
	cout << "Instantaneous Data Read Rate: "<<(double)BYTES_READ/(time_elapsed-time_elapsed_old)/(1024.0*1024.0)<<" MB/s"<<endl;
	cout << (double)TOTAL_BYTES/(1024.0*1024.0*1024.0)<<" GiB Read"<<endl<<endl;
      
	BYTES_READ=0;
	time_elapsed_old = time_elapsed;
      }
    

      //This is the unpacker for the caen2015 data format
      if(strcmp(DataFormat.c_str(),"caen2015") == 0) {

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
	  cout<<"TimeStamp "<<head.fTimeStamp<<endl;    ///< event timestamp in seconds
#endif
	
	  if(head.fEventId==0x8000 || head.fEventId==0x8001 || head.fEventId==0x8002 ){

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
      
	  //Data
	  else if(head.fEventId==1){
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
	      if (bank32.fName[0]=='C' && bank32.fName[1]=='E') {	// name starts as CE VT_BANK
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
		  uint32_t wflen = evaggr->P[evtnum].width;	// CEVT_BANK variable
              
		  for (uint wfindex=where_in_peakbank;wfindex<where_in_peakbank+wflen;++wfindex) {
		    // at this point we have reserved onkly 40 samples in db_arr waveform !!
		    evaggr->wavelets[evtnum][wfindex-where_in_peakbank] = imported_peaks[current_detnum][wfindex];
		  }        
		  where_in_peakbank += wflen;
		  last_detnum = current_detnum;
	      
		  uint64_t timestamp_raw = (evaggr->P[evtnum].position & 0x7FFFFFFFFFFF);                // 47 bits for timestamp

#ifdef Unpacker_Verbose 
		  cout<<"timestamp_raw: "<<timestamp_raw<<endl;
#endif
	      	      
		  db_arr[EVTS].timestamp	= (double)(timestamp_raw);                                     //Digitizer timestamp
		  db_arr[EVTS].TOF       	= (double)(timestamp_raw);                                     //Time of Flight (Currently in 2ns incriments)
		  db_arr[EVTS].Ns		= evaggr->P[evtnum].width;                                     //Number of samples of the waveform
		  db_arr[EVTS].Ifast	= evaggr->P[evtnum].integral[0];                               //Fast integral
		  db_arr[EVTS].Islow	= evaggr->P[evtnum].integral[1]-evaggr->P[evtnum].integral[0]; //Slow integral
		  db_arr[EVTS].board	= (int)((1.*((int)evaggr->P[evtnum].detector_id)-1)/16.);      //Board number
		  db_arr[EVTS].channel	= (1*evaggr->P[evtnum].detector_id-1)-16*db_arr[EVTS].board;   //Channel number
		  db_arr[EVTS].ID     	= MapID[db_arr[EVTS].channel][db_arr[EVTS].board];             //ID from DANCE map
		  db_arr[EVTS].Valid = 1;                                                                //Everything starts valid     
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
		
		  if(db_arr[EVTS].ID<162) frac=0.04;
		  else frac=0.1;
		
		  frac=0.2;
		
		  for(int i=0;i<db_arr[EVTS].Ns;i++) {
		  
		    wf1[i]=evaggr->wavelets[evtnum][i]+8192;

		    //Fill waveform histogram
		    if(read_binary==0) {
		      if(db_arr[EVTS].ID<256) {
			hWaveform_ID->Fill(i,wf1[i],db_arr[EVTS].ID);
			hID_Raw->Fill(db_arr[EVTS].channel+(db_arr[EVTS].board*16));  //Channel + (Board *16)
		      }
		      //  if(db_arr[EVTS].ID==200) {
		      //    hWaveform_T0->Fill(i,wf1[i],1);
		      //  }
		    }
		  
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
		  double dT=0;
		  int iLD=0;
		  for(int i=imin;i>1;i--){
		    if((1.*wf1[i])<thr && (1.*wf1[i-1])>thr){
		      double dSig=(1.*wf1[i-1]-1.*wf1[i]);
		      if(dSig!=0) dT=(1.*wf1[i-1]-thr)/dSig*2.+(i-1)*2.;  // this is in ns
		      else dT=(i-1)*2.;
		      iLD=i;
		    }		
		  }
		
		  //Put the TOF into 1ns units including the CFD time
		  db_arr[EVTS].TOF=dT+2.*db_arr[EVTS].timestamp;
   
		
		  //need to add the time deviations (if any) before time sorting
		  if(db_arr[EVTS].ID < 200) {
		    db_arr[EVTS].TOF += TimeDeviations[db_arr[EVTS].ID];
		  }

		  //keep track of the smallest timestamp
		  if(db_arr[EVTS].TOF<smallest_timestamp) {
		    smallest_timestamp=db_arr[EVTS].TOF;
		  }    
		
		  //increment counters
		  EVTS++;
		  TOTAL_EVTS++;

#ifdef Unpacker_Verbose 
		  cout<<EVTS<<"  "<<TOTAL_EVTS<<endl;
#endif
		}	 //End of loop on eventnum			    
	      }  //End of if on CEVT bank
	      break;
	    } //End of loop on EventBankSize
	  }

	  else { //scalers and other crap
	    char *fData;
	    fData=(char*)malloc(head.fDataSize);
	    gzret=gzread(gz_in,fData,head.fDataSize);
	    free (fData);
	  } //end of scalers and stuff
	}  //checking to see if gzret > 0
      }

      else if(strcmp(DataFormat.c_str(),"caen2018") == 0) {
	
	//Start reading the file
	gzret=gzread(gz_in,&head,sizeof(EventHeader_t));
	
	if(gzret!=0) {
	  
	  TotalDataSize=head.fDataSize;
	  
	  BYTES_READ += head.fDataSize;
	  TOTAL_BYTES += head.fDataSize;
	  
#ifdef Unpacker_Verbose
	  cout<<"Type: "<<head.fEventId<<"  TotalDataSize  "<<TotalDataSize<<endl;
#endif
	
	  // 0x8000 is a begin of run
	  // 0x8001 is an end of run
	  // 0x8002 is an ASCII message created by the logger 
	
	  if(head.fEventId==0x8000 || head.fEventId==0x8001 || head.fEventId==0x8002 ) {
	    //see if its the end of run
	    if(head.fEventId==0x8001) {
	      run=false;
	      break;
	    }
	    
	    char *fData;
	    fData=(char*)malloc(head.fDataSize);
	    gzret=gzread(gz_in,fData,head.fDataSize);
	    free (fData);	
	  }
	  
	  else if(head.fEventId==8){
	    
	    // this is scaler data
	    gzret=gzread(gz_in,&bhead,sizeof(BankHeader_t));	
	    
#ifdef Scaler_Verbose
	    cout << "SCALER " << endl;
	    cout << "Bank_HEADER " << endl;
	    cout << dec <<"TotalBankSize (bytes): " << bhead.fDataSize << endl;
	    cout << dec << bhead.fFlags << endl;
#endif    
	  
	    TotalBankSize = bhead.fDataSize;

	    //This is the numbef of active boards
	    int nactiveboards = 0;
	    
	    while(TotalBankSize > 0) {

	      gzret=gzread(gz_in,&bank,sizeof(BankHeader_t));	
	      TotalBankSize-=sizeof(Bank_t);
	      
#ifdef Scaler_Verbose
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
#ifdef Scaler_Verbose
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

#ifdef Scaler_Verbose
		cout<<endl<<"Digitizer Rates"<<endl;
#endif
		nactiveboards = bank.fDataSize/sizeof(uint32_t);
		for(int eye=0; eye<nactiveboards; eye++) {
		  gzret=gzread(gz_in,&Digitizer_Rates[eye],sizeof(uint32_t));
#ifdef Scaler_Verbose
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
#ifdef Scaler_Verbose
		cout<<endl<<"Acquisition Status"<<endl;
#endif
		nactiveboards = bank.fDataSize/sizeof(uint32_t);
		for(int eye=0; eye<nactiveboards; eye++) {
		  gzret=gzread(gz_in,&Acquisition_Status[eye],sizeof(uint32_t));
#ifdef Scaler_Verbose
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
#ifdef Scaler_Verbose
		cout<<endl<<"Failure Status"<<endl;
#endif
		nactiveboards = bank.fDataSize/sizeof(uint32_t);
		for(int eye=0; eye<nactiveboards; eye++) {
		  gzret=gzread(gz_in,&Failure_Status[eye],sizeof(uint32_t));

		  if(Failure_Status[eye] != 0) {
		    faillog<<"Run: "<<runnum<<"  Board: "<<eye<<" Failure_Status: "<<Failure_Status[eye]<<endl;
		  }
#ifdef Scaler_Verbose
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
#ifdef Scaler_Verbose
		cout<<endl<<"Readout Status"<<endl;
#endif
		nactiveboards = bank.fDataSize/sizeof(uint32_t);
		for(int eye=0; eye<nactiveboards; eye++) {
		  gzret=gzread(gz_in,&Readout_Status[eye],sizeof(uint32_t));
#ifdef Scaler_Verbose
		  cout<<eye<<"  "<<Readout_Status[eye]<<endl;
#endif
		  TotalBankSize-=sizeof(uint32_t);
		  outputdiagnosticsfile << Readout_Status[eye]<<"\n";
		}
	      }

	      //These are the 8500 + 4n register values
	      if (bank.fName[0]=='D' && bank.fName[1]=='I' && bank.fName[2]== 'A' && bank.fName[3]=='G') {
		outputdiagnosticsfile << "DIAG  "<<nactiveboards<<"\n";

#ifdef Scaler_Verbose
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
	      
#ifdef Scaler_Verbose
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
#ifdef Scaler_Verbose
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
	      
#ifdef Scaler_Verbose
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
#ifdef Scaler_Verbose
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
	      
#ifdef Scaler_Verbose
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
#ifdef Scaler_Verbose
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
	      
#ifdef Scaler_Verbose
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
#ifdef Scaler_Verbose
	    cout<<"Done with Scalers. Total Bank Size: "<<TotalBankSize<<endl;
#endif
	  }  //End of Event ID 2
    
	  //Data
	  else if(head.fEventId==1){
	    
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
	      bool readextra = false;
	      if(EventBankSize%8 !=0) {
		readextra = true; 
	      }
	      
	      //Read the firmware version and board ID
	      uint32_t firmware_version;
	      gzret=gzread(gz_in,&firmware_version,sizeof(firmware_version));
	      TotalBankSize -= sizeof(firmware_version);
	      EventBankSize -=  sizeof(firmware_version);
	           
#ifdef Unpacker_Verbose
	      cout<< "board: "<<((firmware_version & 0xFC000000) >> 26)<<" Firmware: "<<(firmware_version & 0xFF)<< "."<<((firmware_version & 0x3F00) >> 8)<<endl;
#endif
	      
	      while(EventBankSize>0) {
	    
		uint32_t dataword =0;
		uint32_t BytesRead =0;
		
		//Unpack digitizer header (4 32 bit words)
		gzret=gzread(gz_in,&v1730_header,sizeof(v1730_header));
		BytesRead += sizeof(v1730_header);
		
		//extract stuff
		uint32_t nwords = v1730_header.dataword_1 & 0xFFFFFFF;
		int header = (v1730_header.dataword_1 & 0xF0000000) >> 28;
		if(header != 10) {
		  cout<<RED<<"Unpacker [ERROR] CAEN Data Header is NOT 10!"<<RESET<<endl;
		  cout<<RED<<"Events: "<<EVTS<<" Total Events: "<<TOTAL_EVTS<<RESET<<endl;
		  cout<<RED<<"Unpacker [ERROR] Data beyond this point would be corrupt and thus I am exiting to analysis!"<<RESET<<endl;
		  run = false;
		  break;
		}
		else {
		  
		  //Check for reported board failures
		  int failbit = (v1730_header.dataword_2 & 0x04000000) >> 26;
		  if(failbit) {
		    cout<<RED<<"Unpacker [WARNING] CAEN Fail Bit on Board "<<((firmware_version & 0xFC000000) >> 26)<<" is Active!"<<RESET<<endl;
		  }

		  int channelmask = (v1730_header.dataword_2 & 0xFF);
		  
#ifdef Unpacker_Verbose
		  cout<< "header: "<<header<<"  nwords: "<<nwords<<" channelmask: "<<channelmask<<endl;
#endif
		  
			  //interpret the channel mask 
		  channels.clear();
		  for(int m=0; m<8; m++) {
#ifdef Unpacker_Verbose
		    cout<<m<<"  "<<((channelmask >> m) & 0x1)<<endl;
#endif
		    if(((channelmask >> m) & 0x1)) {
		      channels.push_back(2*m);
		    }
		  }
		  
		  int wordstoread = nwords-4;
		  int chaggcounter = 0;
		  
		  while (wordstoread>0) {
		    
		    gzret=gzread(gz_in,&v1730_chagg_header,sizeof(v1730_chagg_header));
		    BytesRead += sizeof(v1730_chagg_header);
		    wordstoread -= sizeof(v1730_chagg_header);
		    
		    int chaggsize = (v1730_chagg_header.dataword_1 & 0x3FFFFF);
		    int chagghead = (v1730_chagg_header.dataword_1 & 0x80000000) >> 31;
		    if(chagghead !=1) {
		      cout<<"Unpacker [ERROR] CAEN Channel Aggregate Header is NOT 1"<<endl;
		      return -1;
		    }
		    
		    int nsdb8 = (v1730_chagg_header.dataword_2 & 0xFFFF);
#ifdef Unpacker_Verbose
		    cout<<"Waveform Size: "<<db_arr[EVTS].Ns<<endl;
#endif
		    int extras_format = (v1730_chagg_header.dataword_2 & 0x7000000) >> 24;
#ifdef Unpacker_Verbose 
		    cout<<"Extras Format: "<<extras_format<<endl;
#endif
		    int extras_enabled = (v1730_chagg_header.dataword_2 & 0x10000000) >> 28;
#ifdef Unpacker_Verbose
		    cout<<"Extras Enabled: "<<extras_enabled<<endl;
#endif
		    int chaggwordstoread = chaggsize-2;
		    
		    //Unpack the channel aggreate
		    while (chaggwordstoread>0) {

		      //Need to set the board and Ns here
		      db_arr[EVTS].Valid = 1;
		      db_arr[EVTS].board = ((firmware_version & 0xFC000000) >> 26);
		      db_arr[EVTS].Ns = 8*nsdb8;
		      
		      //TTT and Ch
		      gzret=gzread(gz_in,&dataword,sizeof(dataword));
		      db_arr[EVTS].timestamp = (dataword & 0x7FFFFFFF);
#ifdef Unpacker_Verbose
		      cout<<"TTT: "<<db_arr[EVTS].timestamp<<endl;
#endif	      
		      // int ch = (dataword & 0x80000000) >> 31;
		      db_arr[EVTS].channel = ((dataword & 0x80000000) >> 31) + channels[chaggcounter];

		      //Map it
		      db_arr[EVTS].ID = MapID[db_arr[EVTS].channel][db_arr[EVTS].board];         

		      BytesRead += sizeof(dataword);
		      wordstoread--;
		      chaggwordstoread--;
		      
		      //unpack waveform
		      gzret=gzread(gz_in,&waveform[0],db_arr[EVTS].Ns*sizeof(uint16_t));
		      
		      BytesRead += gzret;
		      wordstoread -= nsdb8*4;
		      chaggwordstoread -= nsdb8*4;
		      
		      for(int eye=0; eye<db_arr[EVTS].Ns; eye++) {
			waveform[eye] = waveform[eye] & 0x3FFF;
		      }
		      
 		      
		      //Extras
		      if(extras_enabled) {
			gzret=gzread(gz_in,&dataword,sizeof(dataword));
			BytesRead += sizeof(dataword);
			wordstoread--;
			chaggwordstoread--;
			if(extras_format <= 2) {
			  
			  //Add the upper bits to the 47-bit timestamp
			  uint16_t timehigh = (dataword & 0xFFFF0000) >> 16;
			  db_arr[EVTS].timestamp += timehigh*2147483648; 
			  
#ifdef Unpacker_Verbose
			  cout<<"time high: "<<timehigh<<" timestamp: "<<db_arr[EVTS].timestamp<<endl;;
#endif	 
			}
		      }
			
		      //Do the waveform analysis
		      // CALCULATE THE LEADING EDGE using constant fraction "frac"
		      int imin=0;
		      double sigmin=1e9;
		      double frac=0.04;
		      double base=0;
		      double secmom=0.;	
		      int NNN=10;
		      // int id=db_arr[EVTS].ID;
		      
		      //  if(db_arr[EVTS].ID<162) frac=0.04;
		      // else frac=0.1;
		
		      frac=0.2;
		
		      for(int i=0;i<db_arr[EVTS].Ns;i++) {

			//Fill waveform histogram
			if(read_binary==0) {
			  if(db_arr[EVTS].ID<256) {
			    hWaveform_ID->Fill(i,waveform[i],db_arr[EVTS].ID,1);
			    hID_Raw->Fill(db_arr[EVTS].channel+(db_arr[EVTS].board*16));  //Channel + (Board *16)
			  }
			  // if(db_arr[EVTS].ID==200) {
			  //   hWaveform_T0->Fill(i,wf1[i],1);
			  // }
			}
			
			if(i<NNN) {
			  base+=(1.*waveform[i]);
			  secmom+=(1.*waveform[i]*1.*waveform[i]);
			}		  
			if((1.*waveform[i])<sigmin) {
			  sigmin=1.*waveform[i];
			  imin=i;
			}
		      }
		
		      double thr=(sigmin-base/(1.*NNN))*frac+base/(1.*NNN);
		      double dT=0;
		      int iLD=0;
		      for(int i=imin;i>1;i--){
			if((1.*waveform[i])<thr && (1.*waveform[i-1])>thr){
			  double dSig=(1.*waveform[i-1]-1.*waveform[i]);
			  if(dSig!=0) dT=(1.*waveform[i-1]-thr)/dSig*2.+(i-1)*2.;  // this is in ns
			  else dT=(i-1)*2.;
			  iLD=i;
			}		
		      }      

		      //TOF
		      db_arr[EVTS].TOF=dT+2.*db_arr[EVTS].timestamp;

		      //need to add the time deviations (if any) before time sorting
		      if(db_arr[EVTS].ID < 200) {
			db_arr[EVTS].TOF += TimeDeviations[db_arr[EVTS].ID];
		      }
		      
		      //keep track of the smallest timestamp
		      if(db_arr[EVTS].TOF<smallest_timestamp) {
			smallest_timestamp=db_arr[EVTS].TOF;
		      }    
		      
		      //Energies
		      gzret=gzread(gz_in,&dataword,sizeof(dataword));
		      db_arr[EVTS].Ifast = (dataword & 0x7FFF);
		      db_arr[EVTS].Islow = (dataword & 0xFFFF0000) >> 16;

#ifdef Unpacker_Verbose	      
		      cout<<"Ifast: "<<db_arr[EVTS].Ifast<<"  ISlow: "<<db_arr[EVTS].Islow<<endl;
#endif
      
		      BytesRead += sizeof(dataword);
		      wordstoread--;	    
		      chaggwordstoread--;
		      
		      //incriment counters
		      EVTS++;
		      TOTAL_EVTS++;

#ifdef Unpacker_Verbose
		      cout<<"chaggwordstoread: "<<chaggwordstoread<<"  wordstoread: "<<wordstoread<<endl;
#endif
		    }
		    chaggcounter++;
		  }
		  
		  EventBankSize -= BytesRead;
		  TotalBankSize -= BytesRead;
		  
		}
	      }
	      
	      if(readextra) {
		uint32_t extra = 0;
		gzret=gzread(gz_in,&extra,sizeof(extra));
		TotalBankSize -= sizeof(extra);
	      }
	    }
	  }
	  else {
	    cout<<"EventID: "<<head.fEventId<<" Unknown"<<endl;
	  }
	}	
      }
      else {
	cout<<RED<<"Unpacker: [ERROR] I dont understand Data Format "<<DataFormat<<RESET<<endl;
	return -1;
      }

      //Do the time sorting
      if(timesort) {

	//At this point we need to start ordering and eventbuilding
	if(EVTS >= BlockBufferSize) {
		
#ifdef Unpacker_Verbose
	  cout<<endl<<"total events at start: "<<EVTS+datadeque.size()<<" deque size: "<<datadeque.size()<<endl;
#endif
	  int EVT_SORT=EVTS;
	  int this_send=0;
	  int deque_size=datadeque.size();
		
	  //the first time through we want to sort and push everything from the first "block" onto the buffer
	  if(!first_sort) {
		  
	    //check to see if there are any timestamps that originate before the first one in the buffer.  If so event building is broken and the analysis is wrong
	    if(event_building_active) {
	      if(smallest_timestamp < datadeque[0].TOF) {
		cout<<"WARNING THE SMALLEST TIMESTAMP IS LOWER THAN THE SMALLEST ONE IN THE DEQUE!!"<<endl;
		cout<<"smallest: "<<smallest_timestamp<<"  smallest in deque: "<<datadeque[0].TOF<<" largest in deque: "<<datadeque[datadeque.size()-1].TOF<<endl;
		cout<<"Deque Depth: "<<(datadeque[datadeque.size()-1].TOF - datadeque[0].TOF)/(1.0e9)<<" seconds"<<endl;
		cout<<"Make the deque: "<<(datadeque[0].TOF-smallest_timestamp)/(1.0e9)<<" seconds deeper"<<endl;
		cout<<"Exiting"<<endl;
		ofstream failfile;
		failfile.open("Failed_Analysis.txt", ios::out | ios::app);
		failfile << "Run: "<<runnum<<" Failed due to insufficient buffer depth...  Add: "<<(datadeque[0].TOF-smallest_timestamp)/(1.0e9)<<" seconds\n";
		failfile.close();
		return -1;
	      }   
	    }
		  
	    //location where the smallest timestamp sits in the sorted data
	    int first_index=0;
		  
	    //find where the smallest time stamp sits in the already time sorted data
	    for(int k=datadeque.size()-1; k>=0; k--) {
	      if(smallest_timestamp >= datadeque[k].TOF) {
		first_index=k-1;
		break;
	      }
	    }
#ifdef Unpacker_Verbose
	    cout<<"start at "<<first_index<<" of "<<datadeque.size()<<endl;
#endif

	    //place everything after that onto the unsorted array
	    for(uint k=first_index; k<datadeque.size(); k++) {
#ifdef Unpacker_Verbose
	      cout<<"EVT_SORT: "<<EVT_SORT<<"  k: "<<k<<endl;
#endif
	      db_arr[EVT_SORT]=datadeque[k];
	      EVT_SORT++;
	    }
		      
	    //remove the ones put onto the array so we dont double things
	    for(int k=0; k<(deque_size-first_index); k++) {
	      datadeque.pop_back();
	    }
		    
#ifdef Unpacker_Verbose
	    cout<<"Deque Size after reduction: "<<datadeque.size()<<"  About to sort "<<EVT_SORT<<" events"<<endl;
#endif
	  } //end loop over !first_sort
		
	  //the first sort is over
	  first_sort=false;
		
	  //sort the unsorted data
	  heapSort(db_arr, EVT_SORT);

	  //push the now sorted data onto the sorted buffer
	  for(int j=0; j<EVT_SORT; j++) {
	    //Update and push it onto the deque if valid
	    last_timestamp[db_arr[j].ID] = db_arr[j].TOF;
	    datadeque.push_back(db_arr[j]);
	  }
		    
#ifdef CheckTheDeque
	  cout<<"Checking deque"<<endl;
	  for(int k=0; k<datadeque.size()-1; k++) {
	    if(datadeque[k+1].TOF < datadeque[k].TOF) {
	      cout<<"problem with entry "<<k<<endl;
	    }
#ifdef Unpacker_Verbose
	    cout<<k<<"  "<<datadeque[k].TOF<<endl;
#endif
	  }
#endif
		
#ifdef Eventbuilder_Verbose
	  cout<<"About to event build deque size: "<<datadeque.size()<<endl;
#endif
	  //Eventbuild
	  while(true) {
		  
	    //check to see if the buffer is longer than the length specificed in global.h
	    if((datadeque[datadeque.size()-1].TOF - datadeque[0].TOF) > buffer_length) {
		    
	      //we have started to build
	      event_building_active=true;
		    
	      //clear the event vector
	      eventvector.clear();
		    
	      //first timestamp in the deque
	      double first_entry_time = datadeque[0].TOF;  //start of the event in time
	      eventvector.push_back(datadeque[0]); //put the first event in the events vector
	      datadeque.pop_front();  //remove the first entry in the deque
	      Events_Sent++;
	      this_send++;
		    
	      bool event_build =true;  //bool to do eventbuilding
		    
	      while(event_build) {
		if(datadeque[0].TOF < (first_entry_time + CoincidenceWindow)) {
		  eventvector.push_back(datadeque[0]); //put the first event in the events vector
		  datadeque.pop_front();  //remove the first entry in the deque     
		  Events_Sent++;
		  this_send++;
		}
		else {
		  event_build = false;
		  break;
		}
	      }
			
	      if(eventvector.size()>0) {
#ifdef Eventbuilder_Verbose
		cout<<"Processing Event with Size: "<<eventvector.size()<<"  " <<datadeque.size()<<" Entries in the deque"<<endl;
#endif
		//Send it to the analyzer
		Analyze_Data(eventvector, read_binary, write_binary,Crystal_Blocking_Time,DEvent_Blocking_Time, HAVE_Threshold, Energy_Threshold,NQGates,QGates);
	      }
	    }
	    else {
#ifdef Eventbuilder_Verbose
	      cout<<"event build complete: "<<datadeque.size()<<endl;
#endif
	      break;
	    }
	  }
		
	  //Reset the event counter and smallest timestamp
	  EVTS=0;
	  smallest_timestamp=400000000000000.0;
		
#ifdef Unpacker_Verbose
	  cout<<"total events at end: "<<datadeque.size()+this_send<<"  Total Sent: "<<Events_Sent<<endl;
#endif
	}	
      } //end of loop on if timesort
      else {
	EVTS=0;
      } //end else


    }  //end of while run

  
    //Now that we are done sorting we need to empty the buffer
    cout<<GREEN<<"Unpacker [INFO]: Finsihed Unpacking Data"<<RESET<<endl;
  
    //see if anything is left in the unsorted part
    if(EVTS>0) {
      cout<<GREEN<<"Unpacker [INFO]: There are "<<EVTS<<" Events left to sort and "<<datadeque.size()<<" Events left in the Buffer" RESET<<endl;
   
#ifdef Eventbuilder_Verbose
      cout<<"Sorting remaing data"<<first_sort<<" "<<EVTS<<endl;
      cout<<endl<<"total events at start: "<<EVTS+datadeque.size()<<" deque size: "<<datadeque.size()<<endl;
#endif
      int EVT_SORT=EVTS;
      int deque_size=datadeque.size();
    
      //the first time through we want to sort and push everything from the first "block" onto the buffer
      if(!first_sort) {
      
	//check to see if there are any timestamps that originate before the first one in the buffer.  If so event building is broken and the analysis is wrong
	if(event_building_active) {
	  if(smallest_timestamp < datadeque[0].TOF) {
	    cout<<RED<<"WARNING THE SMALLEST TIMESTAMP IS LOWER THAN THE SMALLEST ONE IN THE DEQUE!!"<<RESET<<endl;
	    cout<<RED<<"smallest: "<<smallest_timestamp<<"  smallest in deque: "<<datadeque[0].TOF<<RESET<<endl;
	    cout<<RED<<"Exiting"<<RESET<<endl;
	    ofstream failfile;
	    failfile.open("Failed_Analysis.txt", ios::out | ios::app);
	    failfile << "Run: "<<runnum<<" Failed due to insufficient buffer depth...  Add: "<<(datadeque[0].TOF-smallest_timestamp)/(1.0e9)<<" seconds\n";
	    failfile.close();
	    return -1;
	  }   
	}
      
	//location where the smallest timestamp sits in the sorted data
	int first_index=0;
      
	//find where the smallest time stamp sits in the already time sorted data
	for(int k=datadeque.size()-1; k>=0; k--) {
	  if(smallest_timestamp >= datadeque[k].TOF) {
	    first_index=k-1;
	    break;
	  }
	}
#ifdef Eventbuilder_Verbose
	cout<<"start at "<<first_index<<" of "<<datadeque.size()<<endl;
#endif
      	//place everything after that onto the unsorted array
	for(uint k=first_index; k<datadeque.size(); k++) {
#ifdef Eventbuilder_Verbose
	  cout<<"EVT_SORT: "<<EVT_SORT<<"  k: "<<k<<endl;
#endif
	  db_arr[EVT_SORT]=datadeque[k];
	  EVT_SORT++;
	}
      
	//remove the ones put onto the array so we dont double things
	for(int k=0; k<(deque_size-first_index); k++) {
	  datadeque.pop_back();
	}
#ifdef Eventbuilder_Verbose
	cout<<"Deque Size after reduction: "<<datadeque.size()<<"  About to sort "<<EVT_SORT<<" events"<<endl;
#endif
      } //end loop over !first_sort
    
      //the first sort is over
      first_sort=false;
    
      //sort the unsorted data
      heapSort(db_arr, EVT_SORT);
      
      //push the now sorted data onto the sorted buffer
      for(int j=0; j<EVT_SORT; j++) {
	//Update and push it onto the deque if valid
	last_timestamp[db_arr[j].ID] = db_arr[j].TOF;
	datadeque.push_back(db_arr[j]);
      }
#ifdef Eventbuilder_Verbose
      cout<<"About to event build deque size: "<<datadeque.size()<<endl;
#endif
      EVTS=0;
    }
  
    if(timesort && datadeque.size()>0) {
      // uint64_t entry_counter=0;

      while(emptythedeque) {
	eventvector.clear();
	double first_entry_time = datadeque[0].TOF;  //start of the event in time
	eventvector.push_back(datadeque[0]); //put the first event in the events vector
 	datadeque.pop_front();  //remove the first entry in the deque
     	Events_Sent++;
    
	bool event_build =true;  //initilize bool to do eventbuilding
    
	while(event_build && datadeque.size()>0) {
	  if(datadeque[0].TOF < (first_entry_time + CoincidenceWindow)) {
	    eventvector.push_back(datadeque[0]); //put the first event in the events vector
    	    datadeque.pop_front();  //remove the first entry in the deque     
	    Events_Sent++;
	  }
	  else {
	    event_build = false;
	    break;
	  }
	}
    
	if(eventvector.size()>100) {
#ifdef Eventbuilder_Verbose
	  cout<<"Processing Event with Size: "<<eventvector.size()<<"  " <<datadeque.size()<<" Entries in the deque"<<endl;
#endif
	  }
    
	if(eventvector.size()>0) {
	  Analyze_Data(eventvector, read_binary, write_binary,Crystal_Blocking_Time,DEvent_Blocking_Time, HAVE_Threshold,Energy_Threshold,NQGates,QGates);
	}
    
	if(datadeque.size()==0) {
	  cout<<GREEN<<"Unpacker [INFO]: Buffer empty, unpacking complete."<<RESET<<endl;
	  emptythedeque=false;
	}
      } 
    }

  } // END OF MIDAS READER


  //Stage 1 unpacker
  if(read_binary==1) {

    while(run) {
      
      //Event limit control
      if(TOTAL_EVTS > EventLimit) {
	run=false;
      }
      
      //Progress indicator
      if(TOTAL_EVTS > progresscounter*ProgressInterval) {
	progresscounter++;
	cout<<"Processing Run Number: "<<runnum<<endl;
	if(timesort && datadeque.size()>0) {
	  cout<<"Oldest Time in the Buffer: "<<datadeque[0].TOF<<endl;
	  cout<<"Newest Time in the Buffer: "<<datadeque[datadeque.size()-1].TOF<<endl;
	}	
	cout<<TOTAL_EVTS<<" Events Unpacked "<<endl;
	cout<<Events_Sent<<" Events Analyzed"<<endl;
	cout<<datadeque.size()<<" Events in the Buffer"<<endl;
	gettimeofday(&tv,NULL); 
	time_elapsed=tv.tv_sec+(tv.tv_usec/1000000.0);
      
	cout << "Average Event Processing Rate: "<<(double)TOTAL_EVTS/(time_elapsed-unpack_begin)<<" Events per second "<<endl;
	cout << "Average Data Read Rate: "<<(double)TOTAL_BYTES/(time_elapsed-unpack_begin)/(1024.0*1024.0)<<" MB/s"<<endl;
	cout << "Instantaneous Data Read Rate: "<<(double)BYTES_READ/(time_elapsed-time_elapsed_old)/(1024.0*1024.0)<<" MB/s"<<endl;
	cout << (double)TOTAL_BYTES/(1024.0*1024.0*1024.0)<<" GiB Read"<<endl<<endl;
      
	BYTES_READ=0;
	time_elapsed_old = time_elapsed;
      }
      
      gzret=gzread(gz_in,&devt_stage1,sizeof(DEVT_STAGE1));
      // gzseek(gz_in,devt_padding,SEEK_CUR);
      
      if(gzret!=0) {
	
 	BYTES_READ += gzret;
	TOTAL_BYTES += gzret;

	//Fill the array
	db_arr[TOTAL_EVTS].TOF = devt_stage1.TOF;
	db_arr[TOTAL_EVTS].Ifast = devt_stage1.Ifast;
	db_arr[TOTAL_EVTS].Islow = devt_stage1.Islow;
	db_arr[TOTAL_EVTS].ID = devt_stage1.ID;
	db_arr[TOTAL_EVTS].Valid = 1; //Everything starts valid

	
	//need to add the time deviations before time sorting
	if(db_arr[TOTAL_EVTS].ID < 200) {
	  db_arr[TOTAL_EVTS].TOF += TimeDeviations[db_arr[TOTAL_EVTS].ID];
	}
	
	//Add the DANCE delay
	if(db_arr[TOTAL_EVTS].ID < 162) {
	  db_arr[TOTAL_EVTS].TOF += DANCE_Delay;
	}
	
	//Add the He3 delay
	if(db_arr[TOTAL_EVTS].ID == 241) {
	  db_arr[TOTAL_EVTS].TOF += He3_Delay;
	} 
	
	//Add the U235 delay
	if(db_arr[TOTAL_EVTS].ID == 243) {
	  db_arr[TOTAL_EVTS].TOF += U235_Delay;
	} 
	
	//Add the Li6 delay
	if(db_arr[TOTAL_EVTS].ID == 244) {
	  db_arr[TOTAL_EVTS].TOF += Li6_Delay;
	}
       	
	TOTAL_EVTS++;

      }
      else {
	run=false;
	break;
      }
    }  //end of while(run)
    
    cout<<GREEN<<"Unpacker [INFO]: Unpacking Complete"<<RESET<<endl;

    //Since it is stage 1 it is very likely that things are already closely ordered in time but the offsets can play with things so resort 
    cout<<"Unpacker [INFO]: Time Sorting"<<endl;
    heapSort(db_arr, TOTAL_EVTS);
    cout<<GREEN<<"Unpacker [INFO]: Time Sort Complete"<<RESET<<endl;

    cout<<"Unpacker [INFO]: Starting Analysis"<<endl;

    uint64_t stage1_counter=0;  //keep track of entries sorted into events     
    bool event_build =true;  //initilize bool to do eventbuilding
    
    while(event_build) {
      
      double first_entry_time = db_arr[stage1_counter].TOF;  //start of the event in time
      eventvector.push_back(db_arr[stage1_counter]); //put the first event in the events vector
      Events_Sent++;
      stage1_counter++;
      
      bool add_to_current_event = true; 
      
      while(add_to_current_event && stage1_counter<=TOTAL_EVTS) {
	if(db_arr[stage1_counter].TOF < (first_entry_time + CoincidenceWindow)) {
	  eventvector.push_back(db_arr[stage1_counter]); //put the first event in the events vector
	  Events_Sent++;
	  stage1_counter++;
	}
	else {
	  //If the timestamp is beyond the current start of event + coincidence window this event is over
	  add_to_current_event = false;
	  
	  //Send the current event to analyzer
	  if(eventvector.size()>0) {
	    Analyze_Data(eventvector, read_binary, write_binary,Crystal_Blocking_Time,DEvent_Blocking_Time, HAVE_Threshold, Energy_Threshold,NQGates,QGates);
	    eventvector.clear();
	  }
	  
	  break;
	}
      } //end of while(add_to_current_event && eye<=TOTAL_EVTS) 
      
      //Check to see if the next event in the array is beyond the number of events unpacked
      if(stage1_counter >= TOTAL_EVTS + 1) {
	event_build = false;
	break;
      } //end of if(stage1_counter == TOTAL_EVENTS + 1)
    } //end of while(event_build)
  } //end of if(read_binary==1)

  //Make the time deviations if needed (Likely only a stage 0 thing)
  if(FitTimeDev) {
    Make_Time_Deviations(runnum);
  }
  
  //Unpacking is finished.  Return the total number of events unpacked and analyzed
  return TOTAL_EVTS;
  
  
}

  

