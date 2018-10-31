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

  Input_Parameters input_params;
  
  //Initialize input parameters
  input_params.Crystal_Blocking_Time=0;
  input_params.DEvent_Blocking_Time=0;
  input_params.Read_Binary=false;
  input_params.Write_Binary=false;
  input_params.Read_Simulation=false;
  input_params.Coincidence_Window=10;
  input_params.HAVE_Threshold=false;
  input_params.Energy_Threshold=0.15; //MeV
  input_params.FitTimeDev=0;
  input_params.QGatedSpectra=false;
  input_params.NQGates=0;
  input_params.IsomerSpectra=false;
  input_params.NIsomers=0;
  input_params.Evaluate_DeadTime=false;
  input_params.Artificial_TOF=0;
  input_params.Long_Gate=1000;
  input_params.Use_Firmware_FineTime=false;

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
  //  else if(argc==3) {
  //   RunNum = atoi(argv[1]);
  //   cfgfile = argv[2];
  // }
  else if(argc==4) {
    pathtodata = argv[1];
    RunNum = atoi(argv[2]);
    cfgfile = argv[3];
  }
  else {
    cout<<RED<<"Main [ERROR]: Too many or too few arguments provided.  See README file"<<RESET<<endl;
    cout<<RED<<"Main [ERROR]: for Stage0 and Stage1: \"./DANCE_Analysis pathtodata runnumber cfgfile.cfg \""<<RESET<<endl;
    // cout<<RED<<"Main [ERROR]: for Stage1: \"./DANCE_Analysis runnumber cfgfile.cfg \""<<RESET<<endl;
    cout<<RED<<"Main [ERROR]: for Simulations: \"./DANCE_Analysis cfgfile.cfg \""<<RESET<<endl;
    return -1;
  }
  
  input_params.RunNumber = RunNum;

  cout<<"Main [INFO]: Run Number: "<<input_params.RunNumber<<endl;

  ifstream cfgf;
  cfgf.open(cfgfile.c_str());

  if(cfgf.is_open()) {

    cout<<"Main [INFO]: Configuration File "<<cfgfile.c_str()<<" is Open"<<endl;

    string item;
    while(!cfgf.eof()) {
      cfgf>>item;
      if(item.compare("Coincidence_Window") == 0) {
      	cfgf>>input_params.Coincidence_Window;
      }
      if(item.compare("Read_Binary") == 0) {
      	cfgf>>input_params.Read_Binary;
      }
      if(item.compare("Write_Binary") == 0) {
	cfgf>>input_params.Write_Binary;
      }      
      if(item.compare("Read_Simulation") == 0) {
	cfgf>>input_params.Read_Simulation;
      }      
      if(item.compare("Simulation_File_Name") == 0) {
	cfgf>>input_params.Simulation_File_Name;
      } 
      if(item.compare("Crystal_Blocking_Time") == 0) {
	cfgf>>input_params.Crystal_Blocking_Time;
      }    
      if(item.compare("DEvent_Blocking_Time") == 0) {
	cfgf>>input_params.DEvent_Blocking_Time;
      } 
      if(item.compare("HAVE_Threshold") == 0) {
	cfgf>>input_params.HAVE_Threshold;
      } 
      if(item.compare("Energy_Threshold") == 0) {
	cfgf>>input_params.Energy_Threshold;
      }  
      if(item.compare("FitTimeDev") == 0) {
	cfgf>>input_params.FitTimeDev;
      } 
      if(item.compare("DataFormat") == 0) {
	cfgf>>input_params.DataFormat;
      }
      if(item.compare("NQGates") == 0) {
	cfgf>>input_params.NQGates;
	for(int eye=0; eye<2*input_params.NQGates; eye++) {
	  cfgf >> input_params.QGates[eye];
	}
      }
      if(item.compare("NIsomers") == 0) {
	cfgf>>input_params.NIsomers;
	for(int eye=0; eye<input_params.NIsomers; eye++) {
	  cfgf >> input_params.IsomerPromptQGates[2*eye] >> input_params.IsomerPromptQGates[2*eye+1];
	  cfgf >> input_params.IsomerPromptMclGates[2*eye] >> input_params.IsomerPromptMclGates[2*eye+1];
	  cfgf >> input_params.IsomerPromptTOFGates[2*eye] >> input_params.IsomerPromptTOFGates[2*eye+1];
	  cfgf >> input_params.IsomerDelayedQGates[2*eye] >> input_params.IsomerDelayedQGates[2*eye+1];
	  cfgf >> input_params.IsomerDelayedMclGates[2*eye] >> input_params.IsomerDelayedMclGates[2*eye+1];
	  cfgf >> input_params.IsomerDelayedTOFGates[2*eye] >> input_params.IsomerDelayedTOFGates[2*eye+1];
	}
      }

      if(item.compare("JMOD_Background") == 0) {
	cfgf>>input_params.JMOD_Background;
      }

      if(item.compare("Evaluate_DeadTime") == 0) {
	cfgf>>input_params.Evaluate_DeadTime;
      }   

      if(item.compare("Artificial_TOF") == 0) {
	cfgf>>input_params.Artificial_TOF;
      }   
      if(item.compare("DetectorLoad_FileName") == 0) {
	cfgf>>input_params.DetectorLoad_FileName;
      }   
      if(item.compare("DetectorLoad_HistName") == 0) {
	cfgf>>input_params.DetectorLoad_HistName;
      }   
      if(item.compare("Long_Gate") == 0) {
	cfgf>>input_params.Long_Gate;
      }  
      if(item.compare("Use_Firmware_FineTime") == 0) {
	cfgf>>input_params.Use_Firmware_FineTime;
      }  
    }

    //Set the bool for QGates
    if(input_params.NQGates>0) {
      input_params.QGatedSpectra = true;
    }
    else {
      input_params.QGatedSpectra = false;
    }
    
    //Set the bool for Isomers
    if(input_params.NIsomers>0) {
      input_params.IsomerSpectra = true;
    }
    else {
      input_params.IsomerSpectra = false;
    }
    
    
    cout<<GREEN<<"Main [INFO]: Read Configuration File: "<<cfgfile<<RESET<<endl;
    cout<<"Coincidence Window: "<<input_params.Coincidence_Window<<endl;
    cout<<"Read Binary: "<<input_params.Read_Binary<<endl;
    cout<<"Write Binary: "<<input_params.Write_Binary<<endl;
    cout<<"Read Simulation: "<<input_params.Read_Simulation<<endl;
    if(input_params.Read_Simulation) {
      cout<<"Simulated File Name: "<<input_params.Simulation_File_Name<<endl;
    }
    cout<<"Crystal Blocking Time: "<<input_params.Crystal_Blocking_Time<<endl;
    cout<<"DANCE Event Blocking Time: "<<input_params.DEvent_Blocking_Time<<endl;
    cout<<"Have Threshold: "<<input_params.HAVE_Threshold<<endl;
    cout<<"Energy Threshold: "<<input_params.Energy_Threshold<<endl;
    cout<<"Fit Time Deviations: "<<input_params.FitTimeDev<<endl;
    cout<<"Data Format: "<<input_params.DataFormat<<endl;
    cout<<"Number of Q-Value Gates: "<<input_params.NQGates<<endl;

    for(int eye=0; eye<input_params.NQGates; eye++) {
      cout<<"Q-Value Gate "<<eye<<": "<<input_params.QGates[2*eye]<<" to "<<input_params.QGates[2*eye+1]<<" MeV"<<endl;
    }

    for(int eye=0; eye<input_params.NIsomers; eye++) {
      cout<<"Isomer Gate "<<eye<<endl;
      cout<<"Prompt :";
      cout<<"  ESum from: "<<input_params.IsomerPromptQGates[2*eye]<<" to "<<input_params.IsomerPromptQGates[2*eye+1]<<" MeV";
      cout<<"  Mcl from: "<<input_params.IsomerPromptMclGates[2*eye]<<" to "<<input_params.IsomerPromptMclGates[2*eye+1];
      cout<<"  TOF from: "<<input_params.IsomerPromptTOFGates[2*eye]<<" to "<<input_params.IsomerPromptTOFGates[2*eye+1]<<" ns"<<endl;
      cout<<"Delayed :";
      cout<<"  ESum from: "<<input_params.IsomerDelayedQGates[2*eye]<<" to "<<input_params.IsomerDelayedQGates[2*eye+1]<<" MeV";
      cout<<"  Mcl from: "<<input_params.IsomerDelayedMclGates[2*eye]<<" to "<<input_params.IsomerDelayedMclGates[2*eye+1];
      cout<<"  TOF from: "<<input_params.IsomerDelayedTOFGates[2*eye]<<" to "<<input_params.IsomerDelayedTOFGates[2*eye+1]<<" ns"<<endl;	 
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
  if(input_params.Read_Binary == 0 && input_params.Read_Simulation == 0) {

    stringstream midasrunname;
    midasrunname.str();
    if(strcmp(input_params.DataFormat.c_str(),"caen2015") == 0) {
      midasrunname << pathtodata << "/run" << std::setfill('0') << std::setw(6) << RunNum << ".mid";
    }
    else if(strcmp(input_params.DataFormat.c_str(),"caen2018") == 0) {
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
  else if(input_params.Read_Binary == 1 && input_params.Read_Simulation == 0) {
    
    stringstream binaryrunname;
    binaryrunname.str();
    binaryrunname << pathtodata << "/stage0_run_" << RunNum << ".bin";
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
  else if(input_params.Read_Binary == 0 && input_params.Read_Simulation == 1) {
    
    stringstream simulationrunname;
    simulationrunname.str();
    simulationrunname << STAGE0_SIM << "/"<< input_params.Simulation_File_Name <<".bin";
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
    cout<<RED<<"Main [ERROR]: Conflict Between Read_Binary (set to "<<input_params.Read_Binary<<") and Read_Simulation (set to ("<<input_params.Read_Simulation<<"). Exiting"<<RESET<<endl;
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
  Initialize_Analyzer(input_params);
  //  Initialize_Analyzer(Read_Binary, Write_Binary,Read_Simulation,NQGates,QGates);

  //Name of the output root file
  stringstream rootfilename;
  rootfilename.str();
  
  //stage 0
  if(input_params.Read_Binary==0 && input_params.Read_Simulation==0) {
    rootfilename << "./stage0_root/Stage0_Histograms_Run_";
    rootfilename << RunNum;
  }
  //stage 1
  else if(input_params.Read_Binary==1 && input_params.Read_Simulation==0) {
    rootfilename << "./stage1_root/Stage1_Histograms_Run_";
    rootfilename << RunNum;
  }
  //simulation
  else if(input_params.Read_Binary==0 && input_params.Read_Simulation==1) {
    rootfilename << "./stage1_root/Stage1_Histograms_";
    rootfilename << input_params.Simulation_File_Name;
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
  Read_Energy_Calibrations(input_params);

  //start time
  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  int events_analyzed=  Unpack_Data(gz_in, begin, input_params);
  cout<<GREEN<<"Main [INFO]: Analysis Complete. Analyzed: "<<events_analyzed<<" Events"<<RESET<<endl;
  
  //Make the file
  fout = new TFile(Form("%s_%dns_CW_%dns_CBT_%dns_DEBT.root",rootfilename.str().c_str(),
			(int)input_params.Coincidence_Window,
			(int)input_params.Crystal_Blocking_Time,
			(int)input_params.DEvent_Blocking_Time),"RECREATE");
  
  cout<<GREEN<<"Main [INFO]: Rootfile Created"<<RESET<<endl;
  
  //Write histograms
  Write_Unpacker_Histograms(fout, input_params);
  Write_Analyzer_Histograms(fout, input_params);

  //Write the root file
  fout->Write();
  cout<<GREEN<<"Main [INFO]: Rootfile Written"<<RESET<<endl;
  
  //Print the time elapsed for the entire process
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin); ;
  cout<<GREEN<<"Main [INFO]: Time Elapsed: "<<time_elapsed<<" Seconds"<<RESET<<endl;
}
