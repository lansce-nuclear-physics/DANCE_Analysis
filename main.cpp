#include "global.h"
#include "main.h"
#include "unpacker.h"
#include "analyzer.h"

//#define VERBOSE

#include <sys/time.h>

using namespace std;

int main(int argc, char *argv[]) {

//Control Variables
  bool Read_Binary=false;
  bool Write_Binary=false;
  double Coincidence_Window=0;

  string pathtobinary = "./stage0_bin/";
  string pathtoroot = "./RootFiles/";

  // ./DANCE_Analysis  PathToData  RunNunber  .cfgFile
  if(argc>4) {
    cout<<RED<<"Main [ERROR]: Too many arguments provided.  I just need a path, a run number, and a .cfg file"<<RESET<<endl;
    return -1;
  }
  
  //get the file number from command line args
  int RunNum = atoi(argv[2]);
  string pathtodata = argv[1];
  string cfgfile = argv[3];

  ifstream cfgf;
  cfgf.open(cfgfile.c_str());

  if(cfgf.is_open()) {
    string item;
    while(!cfgf.eof()) {
      cfgf>>item;
      if(item.compare("Coincidence_Window") == 0) {
      	cfgf>>Coincidence_Window;
      }
      if(item.compare("Read_Binary") == 0) {
      	cfgf>>Read_Binary;
      }
      if(item.compare("Write_Binary") == 0) {
	cfgf>>Write_Binary;
      }      
    }
    
    cout<<GREEN<<"Main [INFO]: Read Configuration File: "<<cfgfile<<RESET<<endl;
    cout<<"Coincidence Window: "<<Coincidence_Window<<"  "<<Read_Binary<<"  "<<Write_Binary<<endl;
  }
  else {
    cout<<RED<<"Main [ERROR]: Failed to Read Configuration File: "<<cfgfile<<RESET<<endl;
    return -1;
  }

  
   
  //make the file
  gzFile gz_in;
  
  //Figure out what we are reading in.
  //MIDAS Files have a .mid or .mid.gz ending and stage0 are .bin or .bin.gz

    stringstream runname;
    runname.str();
    
    if(Read_Binary == 0) {

    stringstream midasrunname;
    midasrunname.str();
    midasrunname << pathtodata << "/run" << std::setfill('0') << std::setw(6) << RunNum << ".mid";
    cout<<"Main [INFO]: Checking for: "<<midasrunname.str()<<endl;
    
    //Look for uncompressed .mid files
    gz_in=gzopen(midasrunname.str().c_str(),"rb");
    
    //check to see if its open
    if(gz_in) {
      cout<<GREEN<<"Main [INFO]: File "<<midasrunname.str().c_str()<<" Found"<<RESET<<endl;
      runname << midasrunname.str();
    }
    else {
      //look for compressed .mid.gz files
      midasrunname << ".gz";
      cout<<"Main [INFO]: Checking for: "<<midasrunname.str()<<endl;
      gz_in=gzopen(midasrunname.str().c_str(),"rb");
      if(gz_in) {
	cout<<GREEN<<"Main [INFO]: File "<<midasrunname.str().c_str()<<" Found"<<RESET<<endl;
	runname << midasrunname.str();
      }
    }
  }
  
    if(Read_Binary == 1) {
    
    stringstream binaryrunname;
    binaryrunname.str();
    binaryrunname << pathtobinary << "stage0_run_" << RunNum << ".bin";
    cout<<"Main [INFO]: Checking for: "<<binaryrunname.str()<<endl;
    
    //Look for uncompressed .bin files
    gz_in=gzopen(binaryrunname.str().c_str(),"rb");
    
    //check to see if its open
    if(gz_in) {
      cout<<GREEN<<"Main [INFO]: File "<<binaryrunname.str().c_str()<<" Found"<<RESET<<endl;
      runname << binaryrunname.str();
    }
    else {
      //look for compressed .bin.gz files
      binaryrunname << ".gz";
      cout<<"Main [INFO]: Checking for: "<<binaryrunname.str()<<endl;
      gz_in=gzopen(binaryrunname.str().c_str(),"rb");
      if(gz_in) {
	cout<<GREEN<<"Main [INFO]: File "<<binaryrunname.str().c_str()<<" Found"<<RESET<<endl;
	runname << binaryrunname.str();
      }
    }
  }
  
  
  //Make sure after all of that we have a file
  gz_in=gzopen(runname.str().c_str(),"rb");

   //check to see if its open
  if(!gz_in) {
    cout<<RED<<"Main [ERROR]: No Files for run "<<RunNum<< " Found. Exiting."<<RESET<<endl;
    return -1;
  }
  
  //Initialize Things
  Initialize_Analyzer(Read_Binary, Write_Binary);

  stringstream rootfilename;
  rootfilename.str();
  
  //stage 0
  if(Read_Binary==0) {
    rootfilename << "./stage0_root/Stage0_Histograms_Run_";
    rootfilename << RunNum << ".root";
  }
 //stage 1
  if(Read_Binary==1) {
    rootfilename << "./stage1_root/Stage1_Histograms_Run_";
    rootfilename << RunNum << ".root";
  }

  //make the root file  
  TFile *fout = new TFile(Form("%s",rootfilename.str().c_str()),"RECREATE");
  
  struct timeval tv;  						// real time  
  double begin, end, time_elapsed;			// start,stop, elapsed time


  //Initialize some stuff
  Read_Energy_Calibrations(RunNum, Read_Binary);

  //start time
  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  int events_analyzed=  Unpack_Data(gz_in, begin, RunNum, Read_Binary, Write_Binary, Coincidence_Window);
  cout<<GREEN<<"Main [INFO]: Analysis Complete. Analyzed: "<<events_analyzed<<" Events"<<RESET<<endl;
  
  //Write histograms
  Write_Analyzer_Histograms(fout, Read_Binary);

  fout->Write();
  cout<<GREEN<<"Main [INFO]: Rootfile Written"<<RESET<<endl;
  
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin); ;
  cout<<GREEN<<"Main [INFO]: Time Elapsed: "<<time_elapsed<<" Seconds"<<RESET<<endl;
}
