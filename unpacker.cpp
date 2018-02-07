//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  unpacker.cpp           *// 
//*  Last Edit: 02/07/18    *//  
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
#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"

using namespace std;

//Some Global Unpacker Variables (DONT CHANGE UNLESS NEEDED)

/*size of the DEVT array in the unpacker*/
#define MaxDEVTArrSize 10000000  //this should be a big number 

/*size of the block buffer. The unpacker waits this many
  events before looking at the deque and time sorting.  This 
  number has huge efficiency implications so dont change it if
  you dont understand it */
#define BlockBufferSize 350000

/*Time depth (in seconds) of the buffer.  While you cant know if you 
  will fail the program will tell you if you did... */
#define BufferDepth 10000  //Read the whole thing in.  DANCE doesnt always flush till the end anyway. 


//Verbosity and Error Checking
//define CheckTheDeque
//#define Eventbuilder_Verbose
//#define Unpacker_Verbose


//global unpacker variables
double TimeDeviations[200];

//holds the id of each detector from dancemap
int MapID[20][20];

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
    
    fname<<"./TimeDeviations/TimeDeviations_Run_" << std::setfill('0') << std::setw(5) << runnum << ".txt";
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


int Unpack_Data(gzFile &gz_in, double begin, int runnum, bool read_binary, bool write_binary, double CoincidenceWindow, double Crystal_Blocking_Time, double DEvent_Blocking_Time, bool HAVE_Threshold, double Energy_Threshold, bool FitTimeDev) {

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
  
  //initialize DANCE Map
  Make_DANCE_Map();

  //initiliaze the time deviations
  Read_TimeDeviations(runnum,FitTimeDev);

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
  
  //MIDAS Bank Stuff
  EventHeader_t head;           //MIDAS event header
  BankHeader_t bhead;           //MIDAS bank header
  Bank32_t bank32;              //MIDAS 32-bit bank

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
  
  //Start of the unpacking process 
  gettimeofday(&tv,NULL);  
  double unpack_begin = tv.tv_sec+(tv.tv_usec/1000000.0);

  cout<<GREEN<<"Unpacker [INFO]: Started Unpacking: "<<RESET<<endl;
  
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
			  Analyze_Data(eventvector, read_binary, write_binary,Crystal_Blocking_Time,DEvent_Blocking_Time, HAVE_Threshold, Energy_Threshold);
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
		}
		else {
		  EVTS=0;
		}
	      }				    
	    }
	    break;
	  }		
	}
	else { //scalers and other crap
	  char *fData;
	  fData=(char*)malloc(head.fDataSize);
	  gzret=gzread(gz_in,fData,head.fDataSize);
	  free (fData);
	}
      }
    }

  
    //Now that we are done sorting we need to empty the buffer
    cout<<GREEN<<"Unpacker [INFO]: Finsihed Unpacking Data"<<RESET<<endl;
  
    //see if anything is left in the unsorted part
    if(EVTS>0) {
      cout<<GREEN<<"Unpacker [INFO]: There are "<<EVTS<<" Events left to sort and "<<datadeque.size()<<" Events left in the Buffer" RESET<<endl;
   
#ifdef Unpacker_Verbose
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
#ifdef Unpacker_Verbose
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
    
	if(eventvector.size()>0) {
	  //	cout<<"Processing Event with Size: "<<eventvector.size()<<"  " <<datadeque.size()<<" Entries in the deque"<<endl;
	  Analyze_Data(eventvector, read_binary, write_binary,Crystal_Blocking_Time,DEvent_Blocking_Time, HAVE_Threshold,Energy_Threshold);
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
	    Analyze_Data(eventvector, read_binary, write_binary,Crystal_Blocking_Time,DEvent_Blocking_Time, HAVE_Threshold, Energy_Threshold);
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

  

