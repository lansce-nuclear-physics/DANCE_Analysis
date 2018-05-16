//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  main.cpp               *// 
//*  Last Edit: 05/16/18    *//  
//***************************//

//File includes
#include "global.h"
#include "main.h"
#include "unpacker.h"
#include "analyzer.h"

//C/C++ includes
#include <sys/time.h>

using namespace std;

int main(int argc, char *argv[]) {

  //Read in the version information
  ifstream version_file;
  version_file.open(".version");
  
  string filler, version;
  version_file >> filler >> version;
  
  //Output the version information
  cout<<"You are Running DANCE Analysis version "<<version<<endl;

  //Control Variables that get read in from .cfg files passed along to various functions
  double Crystal_Blocking_Time=0;
  double DEvent_Blocking_Time=0;
  bool Read_Binary=false;
  bool Write_Binary=false;
  double Coincidence_Window=0;
  bool HAVE_Threshold=false;
  double Energy_Threshold=0.15; //MeV
  bool FitTimeDev=0;
  string DataFormat;
  string Simulation_File_Name;
  int NQGates;
  double QGates[20];  //This gives 10 pairs
  bool Read_Simulation=false;


  //Control things
  int RunNum=0;
  string pathtodata;
  string cfgfile;

  if(argc<2) {
    cout<<RED<<"Main [ERROR]: Too few arguments provided.  See README file"<<RESET<<endl;
    cout<<RED<<"Main [ERROR]: for Stage0: \"./DANCE_Analysis pathtodata runnumber cfgfile.cfg \""<<RESET<<endl;
    cout<<RED<<"Main [ERROR]: for Stage1: \"./DANCE_Analysis runnumber cfgfile.cfg \""<<RESET<<endl;
    cout<<RED<<"Main [ERROR]: for Simulations: \"./DANCE_Analysis cfgfile.cfg \""<<RESET<<endl;
    return -1;
  }
  
  else if(argc==2) {
    cfgfile = argv[1];
  }
  else if(argc==3) {
    RunNum = atoi(argv[1]);
    cfgfile = argv[2];
  }
  else if(argc==4) {
    pathtodata = argv[1];
    RunNum = atoi(argv[2]);
    cfgfile = argv[3];
  }
  else {
    cout<<RED<<"Main [ERROR]: Too many arguments provided.  See README file"<<RESET<<endl;
    cout<<RED<<"Main [ERROR]: for Stage0: \"./DANCE_Analysis pathtodata runnumber cfgfile.cfg \""<<RESET<<endl;
    cout<<RED<<"Main [ERROR]: for Stage1: \"./DANCE_Analysis runnumber cfgfile.cfg \""<<RESET<<endl;
    cout<<RED<<"Main [ERROR]: for Simulations: \"./DANCE_Analysis cfgfile.cfg \""<<RESET<<endl;
    return -1;
  }
  


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
      if(item.compare("Read_Simulation") == 0) {
	cfgf>>Read_Simulation;
      }      
      if(item.compare("Simulation_File_Name") == 0) {
	cfgf>>Simulation_File_Name;
      } 
      if(item.compare("Crystal_Blocking_Time") == 0) {
	cfgf>>Crystal_Blocking_Time;
      }    
      if(item.compare("DEvent_Blocking_Time") == 0) {
	cfgf>>DEvent_Blocking_Time;
      } 
      if(item.compare("HAVE_Threshold") == 0) {
	cfgf>>HAVE_Threshold;
      } 
      if(item.compare("Energy_Threshold") == 0) {
	cfgf>>Energy_Threshold;
      }  
      if(item.compare("FitTimeDev") == 0) {
	cfgf>>FitTimeDev;
      } 
      if(item.compare("DataFormat") == 0) {
	cfgf>>DataFormat;
      }
      if(item.compare("NQGates") == 0) {
	cfgf>>NQGates;
	for(int eye=0; eye<2*NQGates; eye++) {
	  cfgf >> QGates[eye];
	}
      }
    }
    
    cout<<GREEN<<"Main [INFO]: Read Configuration File: "<<cfgfile<<RESET<<endl;
    cout<<"Coincidence Window: "<<Coincidence_Window<<endl;
    cout<<"Read Binary: "<<Read_Binary<<endl;
    cout<<"Write Binary: "<<Write_Binary<<endl;
    cout<<"Read Simulation: "<<Read_Simulation<<endl;
    if(Read_Simulation) {
      cout<<"Simulated File Name: "<<Simulation_File_Name<<endl;
    }
    cout<<"Crystal Blocking Time: "<<Crystal_Blocking_Time<<endl;
    cout<<"DANCE Event Blocking Time: "<<DEvent_Blocking_Time<<endl;
    cout<<"Have Threshold: "<<HAVE_Threshold<<endl;
    cout<<"Energy Threshold: "<<Energy_Threshold<<endl;
    cout<<"Fit Time Deviations: "<<FitTimeDev<<endl;
    cout<<"Data Format: "<<DataFormat<<endl;
    cout<<"Number of Q-Value Gates: "<<NQGates<<endl;
    for(int eye=0; eye<NQGates; eye++) {
      cout<<"Q-Value Gate "<<eye<<": "<<QGates[2*eye]<<" to "<<QGates[2*eye+1]<<" MeV"<<endl;
    }
  }
  else {
    cout<<RED<<"Main [ERROR]: Failed to Read Configuration File: "<<cfgfile<<RESET<<endl;
    return -1;
  }
  
   
  //make the file handle
  gzFile gz_in;
  
  //Figure out what we are reading in.
  //MIDAS Files have a .mid or .mid.gz ending and staged analysis files are .bin or .bin.gz

  //Name of the midas or bin file
  stringstream runname;
  runname.str();
    
  //Stage 0
  if(Read_Binary == 0 && Read_Simulation == 0) {

    stringstream midasrunname;
    midasrunname.str();
    if(strcmp(DataFormat.c_str(),"caen2015") == 0) {
      midasrunname << pathtodata << "/run" << std::setfill('0') << std::setw(6) << RunNum << ".mid";
    }
    else if(strcmp(DataFormat.c_str(),"caen2018") == 0) {
      midasrunname << pathtodata << "/run" << std::setfill('0') << std::setw(6) << RunNum << ".mid";
    }
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
  
  //Stage 1
  else if(Read_Binary == 1 && Read_Simulation == 0) {
    
    stringstream binaryrunname;
    binaryrunname.str();
    binaryrunname << STAGE0_BIN << "/stage0_run_" << RunNum << ".bin";
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

  //Simulated Data from GEANT4
  else if(Read_Binary == 0 && Read_Simulation == 1) {
    
    stringstream simulationrunname;
    simulationrunname.str();
    simulationrunname << STAGE0_SIM << "/"<< Simulation_File_Name <<".bin";
    cout<<"Main [INFO]: Checking for: "<<simulationrunname.str()<<endl;
    
    //Look for uncompressed .bin files
    gz_in=gzopen(simulationrunname.str().c_str(),"rb");
    
    //check to see if its open
    if(gz_in) {
      cout<<GREEN<<"Main [INFO]: File "<<simulationrunname.str().c_str()<<" Found"<<RESET<<endl;
      runname << simulationrunname.str();
    }
    else {
      //look for compressed .bin.gz files
      simulationrunname << ".gz";
      cout<<"Main [INFO]: Checking for: "<<simulationrunname.str()<<endl;
      gz_in=gzopen(simulationrunname.str().c_str(),"rb");
      if(gz_in) {
	cout<<GREEN<<"Main [INFO]: File "<<simulationrunname.str().c_str()<<" Found"<<RESET<<endl;
	runname << simulationrunname.str();
      }
    }
  }

  else {
    cout<<RED<<"Main [ERROR]: Conflict Between Read_Binary (set to "<<Read_Binary<<") and Read_Simulation (set to ("<<Read_Simulation<<"). Exiting"<<RESET<<endl;
    return -1;
  }
  
  //Make sure after all of that we have a file
  gz_in=gzopen(runname.str().c_str(),"rb");

  //check to see if its open
  if(!gz_in) {
    cout<<RED<<"Main [ERROR]: No Files for run "<<RunNum<< " Found. Exiting."<<RESET<<endl;
    return -1;
  }
  
  //Initialize Things
  Initialize_Analyzer(Read_Binary, Write_Binary,Read_Simulation,NQGates,QGates);

  //Name of the output root file
  stringstream rootfilename;
  rootfilename.str();
  
  //stage 0
  if(Read_Binary==0 && Read_Simulation==0) {
    rootfilename << "./stage0_root/Stage0_Histograms_Run_";
    rootfilename << RunNum;
  }
  //stage 1
  else if(Read_Binary==1 && Read_Simulation==0) {
    rootfilename << "./stage1_root/Stage1_Histograms_Run_";
    rootfilename << RunNum;
  }
  
  else if(Read_Binary==0 && Read_Simulation==1) {
    rootfilename << "./stage1_root/Stage1_Histograms_";
    rootfilename << Simulation_File_Name;
  }
  else {
    cout<<RED<<"Main [ERROR]: Cannot understand options for making rootfile. Exiting."<<RESET<<endl;
    return -1;
  }
    
  //make the root file  
  TFile *fout;
  
  //Time profiling stuff
  struct timeval tv;  						// real time  
  double begin, end, time_elapsed;			// start,stop, elapsed time


  //Initialize some stuff
  Read_Energy_Calibrations(RunNum, Read_Binary, Read_Simulation);

  //start time
  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  int events_analyzed=  Unpack_Data(gz_in, begin, RunNum, Read_Binary, Write_Binary, Read_Simulation, Coincidence_Window,Crystal_Blocking_Time,DEvent_Blocking_Time, HAVE_Threshold, Energy_Threshold,FitTimeDev,DataFormat,NQGates,QGates);
  cout<<GREEN<<"Main [INFO]: Analysis Complete. Analyzed: "<<events_analyzed<<" Events"<<RESET<<endl;
  
  //Make the file
  fout = new TFile(Form("%s_%dns_CW_%dns_CBT_%dns_DEBT.root",rootfilename.str().c_str(),(int)Coincidence_Window,(int)Crystal_Blocking_Time,(int)DEvent_Blocking_Time),"RECREATE");
  
  cout<<GREEN<<"Main [INFO]: Rootfile Created"<<RESET<<endl;
  
  //Write histograms
  Write_Unpacker_Histograms(fout, Read_Binary);
  Write_Analyzer_Histograms(fout, Read_Binary, Read_Simulation,NQGates,QGates);

  //Write the root file
  fout->Write();
  cout<<GREEN<<"Main [INFO]: Rootfile Written"<<RESET<<endl;
  
  //Print the time elapsed for the entire process
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin); ;
  cout<<GREEN<<"Main [INFO]: Time Elapsed: "<<time_elapsed<<" Seconds"<<RESET<<endl;
}
