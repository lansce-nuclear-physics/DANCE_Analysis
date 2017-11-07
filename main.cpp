#include "global.h"
#include "main.h"
#include "unpacker.h"
#include "analyzer.h"

//#define VERBOSE

#include <sys/time.h>

using namespace std;

int main(int argc, char *argv[]) {


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
  int RunNum = atoi(argv[1]);

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

  //make the root file  
  TFile *fout = new TFile(Form("%sHistograms_Run_%d.root",pathtoroot.c_str(),RunNum),"RECREATE");
  
  struct timeval tv;  						// real time  
  double begin, end, time_elapsed;			// start,stop, elapsed time


  //Initialize some stuff
  Read_Energy_Calibrations(RunNum);

  //start time
  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  if(midas_format) {
    int events_analyzed=  Unpack_Data(gz_in, begin, RunNum);
    cout<<"analyzed: "<<events_analyzed<<" events"<<endl;
  }
  
  //Write histograms
  Write_Analyzer_Histograms(fout);

  fout->Write();
  cout<<"Rootfile Written"<<endl;
  
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin); ;
  printf("\n Analysis Complete:  Time elapsed: %7.2lf\n", time_elapsed);
  
}
