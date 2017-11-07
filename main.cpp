#include "global.h"
#include "main.h"
#include "unpacker.h"
#include "analyzer.h"

//#define VERBOSE

#include <sys/time.h>

using namespace std;


double CalculateCFD_500MSPS(std::vector<Double_t> trace,  //waveform data
			    double CFDSF,              //CFD fraction 
			    int CFDDelay,           //CFD delay
			    int CFDThreshold,       //CFD threshold
			    int length,             //length of smoothing
			    int gap,                //gap between smoothing sums 
			    int boardid,            //boardid
			    int chanid) {           //chanid 
  
  CFDOut.clear();
  
  Double_t cfdvalue = -1;
  Int_t cfdtrigpoint = -1;
  Int_t A=0;
  Int_t m=0;

  Double_t S0;
  Double_t S1;
  Double_t S2;
  Double_t S3;

  //The 500 MSPS modules have a fixed length and gap of 5 and 1 clock ticks respectively
  
  double CFDWeight = CFDSF;
  double CFDTriggerPoint =0;
  double CFDTime=0;
  
  if(trace.size() > 0) {
    
    //Calculate the CFD Response
    for(m=0; m<(CFDDelay+2*length+gap); m++) {
      CFDOut.push_back(0);
    }
    
    m=CFDDelay+2*length+gap;
    
    for(A=(m-CFDDelay-2*length-gap); A<(m-CFDDelay-gap-length); A++){
      S3+=trace[A];
    }
    for(A=(m-CFDDelay-length); A<(m-CFDDelay); A++) {
      S2+=trace[A];
    }  
    for(A=(m-2*length-gap); A<(m-gap-length); A++){
      S1+=trace[A];
    }
    for(A=m-length; A<m; A++){
      S0+=trace[A];
    }
   
    for(m=(CFDDelay+2*length+gap+1); m<(Int_t)(trace.size()); m++) {
     
      S3=S3-trace[m-CFDDelay-2*length-gap-1]+trace[m-CFDDelay-gap-length-1];
      S2=S2-trace[m-CFDDelay-length-1]+trace[m-CFDDelay-1];
      S1=S1-trace[m-2*length-gap-1]+trace[m-gap-length-1];
      S0=S0-trace[m-length-1]+trace[m-1];
      
      cfdvalue = CFDWeight*(S0-S1)-(S2-S3);
      CFDOut.push_back(cfdvalue);     
      
      //Fill histogram
      hCFD[boardid][chanid]->Fill(m,cfdvalue,1);
    }

    //Find the CFD over threshold point
    for(m=0; m<(Int_t)trace.size(); m++) {
      if(CFDOut[m-1]<=CFDThreshold && CFDOut[m]>CFDThreshold) {
	cfdtrigpoint = m-1;
	break;
      }
    }
    
    //Find the zero crossing after threshold and calculate time stamp
    for(m=cfdtrigpoint; m<(Int_t)trace.size(); m++) {
      if(CFDOut[m-1]>=0 && CFDOut[m]<0) {
       	CFDTime = CFDOut[m-1]/(CFDOut[m-1]+fabs(CFDOut[m]));
      	CFDTriggerPoint = m-1;
	break;
      }
    }    
  }
  else {
    cout<<"Trace is length 0!"<<endl;
    return -105;
  }

  return CFDTriggerPoint + CFDTime;  //this is the location of the waveform where the CFD triggers
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Analyze_Data(): Where all the magic happens.  
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


/*

int Analyze_Data(vector<RawData> events) {
  
  double ring_time=0;
  double wedge_time=0;

  for(int i=0; i<MAXNBOARDS; i++) {
    Board_Mult[i]=0;
    for(int j=0; j<CHPERBOARD; j++) {
      channel_present[i][j]=false;
      time_stamps[i][j]=-1;
      last_time_stamp[i][j]=0;
      channel_energy[i][j]=0;
    }
  }
  
  // int boardofinterest=0;

  for (int i=0; i<events.size(); i++) {
    
    //  cout<<"board: "<<events[i].boardid<<"  channel: "<<events[i].chanid<<"  energy: "<<events[i].qlong<<endl;

    hEnergy[events[i].boardid][events[i].chanid]->Fill(events[i].qlong);
    hHit[events[i].boardid]->Fill(events[i].timefull,events[i].chanid,1);
      
    //Fill the waveform
    WaveForm.clear();
    Processed_WaveForm.clear();

    for(int k=0; k<events[i].Ns; k++) {
      WaveForm.push_back(events[i].waveform[k]);
      hWaveform[events[i].boardid][events[i].chanid]->Fill(k,WaveForm[k],1);
    }
    
    //Calculate a Baseline
    double Baseline=0;
    for(int k=0; k<120; k++) {
      Baseline += WaveForm[k];
    }
    Baseline /= 120.0;

    //negative polarity signals 
    if(Polarity[events[i].boardid] == 1) {
      for (int k=0; k<WaveForm.size(); k++) {
	Processed_WaveForm.push_back((double)(-1*(WaveForm[k]-Baseline)+2000.0));
	hScaledWaveform[events[i].boardid][events[i].chanid]->Fill(k,Processed_WaveForm[k],1);
      }
    }
    //positive polarity signals 
    else {
      for (int k=0; k<WaveForm.size(); k++) {
	Processed_WaveForm.push_back(WaveForm[k]-(Baseline-2000.0));
	hScaledWaveform[events[i].boardid][events[i].chanid]->Fill(k,Processed_WaveForm[k],1);
      }
    }
    

    double CFDTimeStamp = CalculateCFD_500MSPS(Processed_WaveForm,CFD_scale,CFD_delay,CFD_threshold,CFD_length,CFD_gap,events[i].boardid,events[i].chanid);
    hCFDTimestamp[events[i].boardid]->Fill(CFDTimeStamp,events[i].chanid,1);

    double FineTimeStamp = events[i].timefull - 200.0 + CFDTimeStamp;    
    
    channel_present[events[i].boardid][events[i].chanid]=true;
    time_stamps[events[i].boardid][events[i].chanid] = FineTimeStamp;
    
    //  cout.precision(15);
    //  cout<<"board: "<<events[i].boardid<<"  channel: "<<events[i].chanid<<"  energy: "<<events[i].qlong<<"  time: "<< FineTimeStamp<<endl;
    
    last_time_stamp[events[i].boardid][events[i].chanid] = FineTimeStamp;
        
    for(int k=CFDTimeStamp; k<CFDTimeStamp+400; k++) {
      channel_energy[events[i].boardid][events[i].chanid] += (Processed_WaveForm[k]-2000.0);
    }
    
    channel_energy[events[i].boardid][events[i].chanid] /= 40.0;
    
    if(events[i].boardid==0 && events[i].chanid == 4) {
      if( channel_energy[events[i].boardid][events[i].chanid] > 30000.0 && channel_energy[events[i].boardid][events[i].chanid] < 40000.0) {
	ring_time=FineTimeStamp;
	//	cout<<event_counter<<" Ring: "<<ring_time<<endl;
      }
      else {
	ring_time=-1;
      }
    }
    
    if(events[i].boardid==1 && events[i].chanid == 4) {
      if( channel_energy[events[i].boardid][events[i].chanid] > 30000.0 &&  channel_energy[events[i].boardid][events[i].chanid] < 40000.0) {
	wedge_time=FineTimeStamp;
	//	cout<<event_counter<<" Wedge: "<<wedge_time<<endl;
      }
      else {
	wedge_time=-1;
      }
    }
    
    hEnergy2D[events[i].boardid]->Fill(channel_energy[events[i].boardid][events[i].chanid],events[i].chanid,1);
    
    Board_Mult[events[i].boardid]++;
    Events_Analyzed++;
    
  }

  //Front back timing
  if(Board_Mult[0]==1 && Board_Mult[1]==1) {
    if(ring_time > 0 && wedge_time>0) {
      hRing_m_Wedge->Fill(ring_time-wedge_time);
      // cout<<event_counter<<" Diff: "<<ring_time-wedge_time<<endl;
    }   
  }

  for(int k=0; k<MAXNBOARDS; k++) {
    for(int l=0; l<CHPERBOARD-1; l++) {
    
      if(channel_present[k][l] && channel_present[k][l+1]) {
	if(channel_energy[k][l] > 5000 && channel_energy[k][l] < 40000 && channel_energy[k][l+1] > 5000 && channel_energy[k][l+1] < 40000) {
	  hTDiff[k]->Fill(time_stamps[k][l+1]-time_stamps[k][l],l,1);
	  //cout<<k<<"  "<<l<<"  "<<time_stamps[k][l+1]-time_stamps[k][l]<<endl;
	  hCoinc[k]->Fill(l);
	}	
      }
    }      
  }

  event_counter++;
}
  
*/

int main(int argc, char *argv[]) {

  //Set CFD parameters
  CFD_scale = 1.0;
  CFD_delay = 40;
  CFD_threshold = 8000;
  CFD_length = 40;
  CFD_gap = 10;

  /*
  //Make Histograms
    
  for(int i=0; i<MAXNBOARDS; i++) {
    hTimeStampDiff[i] = new TH2D(Form("hTimestampDiff_Board_%d",i),Form("hTimestampDiff_Board_%d",i),4000,0,4000,16,0,16);
    hTimeStampDiff2[i] = new TH2D(Form("hTimestampDiff2_Board_%d",i),Form("hTimestampDiff2_Board_%d",i),4000,0,4000,16,0,16);
    hCFDTimestamp[i] = new TH2D(Form("hCFDTimestamp_Board_%d",i),Form("hCFDTimestamp_Board_%d",i),4000,0,800,16,0,16);
    hTDiff[i] = new TH2D(Form("hTDiff_Board_%d",i),Form("hTDiff_Board_%d",i),10000,-1000,1000,16,0,16);
    hEnergy2D[i] = new TH2D(Form("hEnergy_Board_%d",i),Form("hEnergy_Board_%d",i),1000,0,100000,16,0,16);
    hCoinc[i] = new TH1D(Form("hCoinc_Board_%d",i),Form("hCoinc_Board_%d",i),16,0,16);
    hHit[i] = new TH2D(Form("hHit_Board_%d",i),Form("hHit_Board_%d",i),10000,0,0.5e11,16,0,16);
    for(int j=0; j<CHPERBOARD; j++) {
      hEnergy[i][j] = new TH1D(Form("hEnergy_%d_%d",i,j),Form("hEnergy_%d_%d",i,j),100000,0,100000);
      hSWLongGate[i][j] = new TH1D(Form("hSWLongGate_%d_%d",i,j),Form("hSWLongGate_%d_%d",i,j),100000,0,100000000);
      hWaveform[i][j] = new TH2D(Form("hWaveform_%d_%d",i,j),Form("hWaveform_%d_%d",i,j),800,0,800,1800,0,18000);
      hScaledWaveform[i][j] = new TH2D(Form("hScaledWaveform_%d_%d",i,j),Form("hScaledWaveform_%d_%d",i,j),800,0,800,1800,-2000,16000);
      hScaledWaveform2[i][j] = new TH2D(Form("hScaledWaveform2_%d_%d",i,j),Form("hScaledWaveform2_%d_%d",i,j),800,0,800,1800,-2000,16000);
      hCFD[i][j] = new TH2D(Form("hCFD_%d_%d",i,j),Form("hCFD_%d_%d",i,j),800,0,800,2000,-100000,100000);
    }
  }
  
  */
  string pathtomidas = "./MIDAS_Data/";
  string pathtoroot = "./RootFiles/";

  bool midas_format = false;
  bool digites_format = false;

  //takes one argument the run number;
  if(argc>2) {
    cout<<"Too many arguments provided.  I just need a run number"<<endl;
    return -1;
  }
  
  //get the file number from command line args
  RunNum = atoi(argv[1]);

  //make the file
  gzFile gz_in;
  
  //Figure out what we are reading in.
  //MIDAS Files have a .mid or .mid.gz ending and digites are .bin or .bin.gz

  stringstream runname;
  runname.str();

  stringstream midasrunname;
  midasrunname.str();
  midasrunname << pathtomidas << "run" << std::setfill('0') << std::setw(6) << RunNum << ".mid";
  cout<<"Checking for: "<<midasrunname.str()<<endl;

  //Look for uncompressed .mid files
  gz_in=gzopen(midasrunname.str().c_str(),"rb");
  
  //check to see if its open
  if(gz_in) {
    cout<<"File "<<midasrunname.str().c_str()<<" Found"<<endl;
    midas_format = true;
    digites_format = false;
    runname << midasrunname.str();
  }
  else {
    //look for compressed .mid.gz files
    midasrunname << ".gz";
    cout<<"Checking for: "<<midasrunname.str()<<endl;
    gz_in=gzopen(midasrunname.str().c_str(),"rb");
    if(gz_in) {
      cout<<"File "<<midasrunname.str().c_str()<<" Found"<<endl;
      midas_format = true;
      digites_format = false;
      runname << midasrunname.str();
    }
  }
  
  
  //Make sure after all of that we have a file
  cout<<"Run Name: "<<runname.str()<<endl;

  gz_in=gzopen(runname.str().c_str(),"rb");
  
  //check to see if its open
  if(!gz_in) {
    cout<<"No Files for run "<<RunNum<< " Found. Exiting."<<endl;
    return -1;
  }
  
  //Initialize Things
  Initialize_Analyzer();

 
  
  fout = new TFile(Form("%sHistograms_Run_%d.root",pathtoroot.c_str(),RunNum),"RECREATE");

  /*
  for(int i=0; i<MAXNBOARDS; i++) {
    for(int j=0; j<CHPERBOARD; j++) {
      hEnergy[i][j] = new TH1D(Form("Energy_%d_%d",i,j),Form("Energy_%d_%d",i,j),100000,0,100000);
    }
  }
  */
  
  struct timeval tv;  						// real time  
  double begin, end, time_elapsed;			// start,stop, elapsed time

  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  if(midas_format) {

    int events_analyzed=  Unpack_Data(gz_in, begin, RunNum);
    
    cout<<"analyzed: "<<events_analyzed<<" events"<<endl;
  }
      
  WriteHistograms(fout);
  
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin); ;
  printf("\nTime elapsed:		%7.2lf\n", time_elapsed);
  
}
