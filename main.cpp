//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  main.cpp               *// 
//*  Last Edit: 09/08/19    *//  
//***************************//

//File includes
#include "global.h"
#include "main.h"
#include "unpacker.h"
#include "analyzer.h"
#include "calibrator.h"
#include "validator.h"
#include "eventbuilder.h"

//Root include
#include "TROOT.h"

//C/C++ includes
#include <sys/time.h>
#include <queue>

using namespace std;

stringstream mmsg;

int main(int argc, char *argv[]) {

  int func_ret=0;

  const char* rootver=gROOT->GetVersion();
  if (rootver[0]!='5' || rootver[2]=='2' || rootver[3]=='6'){
       DANCE_Error("Main","Wrong root version, go to 5.34/36, otherwise time deviations won't make correctly");
    return -1;
  }

  //Read in the version information
  ifstream version_file;
  version_file.open(".version");
  
  string filler, version;
  version_file >> filler >> version;
  

  //Output the version information
  mmsg.str("");
  mmsg<<"DANCE Analyzer Version "<<version<<" Initializing";
  DANCE_Init("Main",mmsg.str());

  Analysis_Parameters analysis_params;
  for(int eye=0; eye<256; eye++) {
    analysis_params.last_timestamp[eye]=0;
    analysis_params.last_Islow[eye]=65535;
    analysis_params.last_Eslow[eye]=20;
    analysis_params.last_Efast[eye]=20;

    analysis_params.last_valid_timestamp[eye]=0;
    analysis_params.last_valid_Islow[eye]=65535;
    analysis_params.last_valid_Eslow[eye]=20;

    analysis_params.entries_unpacked=0;
    analysis_params.entries_awaiting_timesort=0;
    analysis_params.entries_written_to_binary=0;
    analysis_params.entries_processed=0;
    analysis_params.entries_invalid=0;
    analysis_params.entries_built=0;
    analysis_params.events_built=0;

    analysis_params.entries_analyzed=0;
    analysis_params.DANCE_entries_analyzed=0;
    analysis_params.T0_entries_analyzed=0;
    analysis_params.He3_entries_analyzed=0;
    analysis_params.Li6_entries_analyzed=0;
    analysis_params.Bkg_entries_analyzed=0;
    analysis_params.U235_entries_analyzed=0;
    analysis_params.Unknown_entries=0;

    analysis_params.events_analyzed=0;
    analysis_params.DANCE_events_analyzed=0;
    analysis_params.T0_events_analyzed=0;
    analysis_params.He3_events_analyzed=0;
    analysis_params.Li6_events_analyzed=0;
    analysis_params.Bkg_events_analyzed=0;
    analysis_params.U235_events_analyzed=0;

    analysis_params.max_buffer_utilization=0;
    
    analysis_params.first_sort=true;
    analysis_params.event_building_active=false; 
    analysis_params.smallest_timestamp=2.814749767e14; //this is the largest the clock can be on CAEN boards
    analysis_params.largest_timestamp=0;
    analysis_params.largest_subrun_timestamp=0;
    analysis_params.last_subrun_timestamp=0;
  }
  
  Input_Parameters input_params;

  //Initialize input parameters
  input_params.Crystal_Blocking_Time=0;
  input_params.DEvent_Blocking_Time=0;
  input_params.Read_Binary=false;
  input_params.Write_Binary=false;
  input_params.Read_Simulation=false;
  input_params.SingleSubrun=false;
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
  input_params.Analysis_Stage = 0;
  input_params.Buffer_Depth = 10;
      
  //Control things
  int RunNum=0;
  int SubRunNum=0;
  string pathtodata;
  string cfgfile;
  
  //Make sure the number of arguments is reasonable
  if(argc==2) {
    cfgfile = argv[1];
  }
  else if(argc==4) {
    pathtodata = argv[1];
    RunNum = atoi(argv[2]);
    cfgfile = argv[3];
  }
  else if(argc==5) {
    pathtodata = argv[1];
    RunNum = atoi(argv[2]);
    cfgfile = argv[4];
    SubRunNum = atoi(argv[3]);
    input_params.SingleSubrun = true;
  }

  else {
    DANCE_Error("Main","Too many or too few arguments provided.  See README file");
    DANCE_Error("Main","for Stage 0 and Stage 1: \"./DANCE_Analysis pathtodata runnumber (optional subrunnumber) cfgfile.cfg");
    DANCE_Error("Main","for Simulations: \"./DANCE_Analysis cfgfile.cfg");
    return -1;
  }

  //Set the run number
  input_params.RunNumber = RunNum;
  input_params.NumSubRun = 0;

  //Display the run number
  mmsg.str("");
  mmsg<<"Run Number: "<<input_params.RunNumber;
  DANCE_Info("Main",mmsg.str());

  DANCE_Info("Main","Opening Configuration File");
  
  ifstream cfgf;
  cfgf.open(cfgfile.c_str());

  if(cfgf.is_open()) {
    
    mmsg.str("");
    mmsg<<"Configuration File "<<cfgfile.c_str()<<" is Open";
    DANCE_Success("Main",mmsg.str());

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
      if(item.compare("Analysis_Stage") == 0) {
	cfgf>>input_params.Analysis_Stage;
      }
      if(item.compare("Buffer_Depth") == 0) {
	cfgf>>input_params.Buffer_Depth;
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
    
    mmsg.str("");
    mmsg<<"Read Configuration File: "<<cfgfile;
    DANCE_Success("Main",mmsg.str());

    cout<<"Analysis Stage: "<<input_params.Analysis_Stage<<endl;
    cout<<"Coincidence Window: "<<input_params.Coincidence_Window<<endl;
    cout<<"Read Binary: "<<input_params.Read_Binary<<endl;
    cout<<"Write Binary: "<<input_params.Write_Binary<<endl;
    cout<<"Read Simulation: "<<input_params.Read_Simulation<<endl;
    if(input_params.Read_Simulation) {
      cout<<"Simulated File Name: "<<input_params.Simulation_File_Name<<endl;
    }
 
    cout<<"Buffer Depth: "<<input_params.Buffer_Depth<<" seconds"<<endl;
     
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
    cout<<"Long Gate (\"Minimum Time Between Crystal Hits\"): "<< input_params.Long_Gate<<endl;
    cout<<"Use Firmware Fine Time: "<<input_params.Use_Firmware_FineTime<<endl;
  }
  
  //If no configuration file then exit
  else {
    mmsg.str("");
    mmsg<<"Failed to Read Configuration File: "<<cfgfile;
    DANCE_Error("Main",mmsg.str());
    return -1;
  }
 
  //make the file handle
  gzFile gz_in;
  queue<gzFile> gz_queue;

  
  //Figure out what we are reading in.
  //MIDAS Files have a .mid or .mid.gz ending and staged analysis files are .bin or .bin.gz

  //Name of the midas or bin file
  stringstream runname;
  runname.str();
    
  //MIDAS input
  if(input_params.Read_Binary == 0 && input_params.Read_Simulation == 0) {

    stringstream midasrunname;
    midasrunname.str();
    stringstream midassubrunname;
    midassubrunname.str();

    midasrunname << pathtodata << "/run" << std::setfill('0') << std::setw(6) << RunNum << ".mid";
    midassubrunname << pathtodata << "/run" << std::setfill('0') << std::setw(6) << RunNum << "_";
    if (input_params.SingleSubrun){
      midassubrunname << std::setw(3) << SubRunNum;
    }
    else {
      midassubrunname << std::setw(3) << 0;
    } 
    midassubrunname << ".mid";

    mmsg.str("");
    mmsg<<"Checking for: "<<midassubrunname.str();
    DANCE_Info("Main",mmsg.str());
    
    //Look for uncompressed .mid subrun files
    gz_in=gzopen(midassubrunname.str().c_str(),"rb");
    
    //check to see if its open
    if (input_params.SingleSubrun){
       if(gz_in) {
         mmsg.str("");
         mmsg<<"File "<<midassubrunname.str().c_str()<<" Found";
         DANCE_Success("Main",mmsg.str());
          
         runname << midassubrunname.str();
         gz_queue.push(gz_in);
      }
      else {
        midassubrunname << ".gz";         
        mmsg<<"Checking for: "<<midassubrunname.str();
        DANCE_Info("Main",mmsg.str());
        input_params.SubRunNumber=0;
        
        gz_in=gzopen(midassubrunname.str().c_str(),"rb");
        if (gz_in) {
          mmsg.str("");
          mmsg<<"File "<<midassubrunname.str().c_str()<<" Found";
          DANCE_Success("Main",mmsg.str());
          
          runname << midassubrunname.str();
          gz_queue.push(gz_in);
        }
      }
    }
    else { 
        if(gz_in) {
          mmsg.str("");
          mmsg<<"File "<<midassubrunname.str().c_str()<<" Found";
          DANCE_Success("Main",mmsg.str());
          
          runname << midassubrunname.str();
          input_params.SubRunNumber=0;
 
          while(gz_in) {
            input_params.NumSubRun++;
            gz_queue.push(gz_in);
            midassubrunname.str("");
            midassubrunname << pathtodata << "/run" << std::setfill('0') << std::setw(6) << RunNum << "_" << std::setw(3) << input_params.NumSubRun << ".mid";
            gz_in=gzopen(midassubrunname.str().c_str(),"rb");
          }
        }
      else {
        //look for compressed .mid.gz subrun files
        midassubrunname << ".gz";
        
        mmsg.str("");
        mmsg<<"Checking for: "<<midassubrunname.str();
        DANCE_Info("Main",mmsg.str());
        input_params.SubRunNumber=0;
        
        gz_in=gzopen(midassubrunname.str().c_str(),"rb");
        if(gz_in) {
          mmsg.str("");
          mmsg<<"File "<<midassubrunname.str().c_str()<<" Found";
          DANCE_Success("Main",mmsg.str());

          runname << midassubrunname.str();
          input_params.SubRunNumber=0;
          while(gz_in) {
            input_params.NumSubRun++;
            gz_queue.push(gz_in);
            midassubrunname.str("");
            midassubrunname << pathtodata << "/run" << std::setfill('0') << std::setw(6) << RunNum << "_" << std::setw(3) << input_params.NumSubRun << ".mid.gz";
            gz_in=gzopen(midassubrunname.str().c_str(),"rb");
          }
        }
        else { //look for uncompressed .mid files (no subrun)
          gz_in=gzopen(midasrunname.str().c_str(),"rb");
        
          if(gz_in) {
            mmsg.str("");
            mmsg<<"File "<<midasrunname.str().c_str()<<" Found";
            DANCE_Success("Main",mmsg.str());
            runname << midasrunname.str();
            input_params.SubRunNumber=-1;
            gz_queue.push(gz_in);
          }
          else { //look for .mid.gz files (no subrun)
            midasrunname << ".gz";
            gz_in=gzopen(midasrunname.str().c_str(),"rb");
          
            if(gz_in) {
              mmsg.str("");
              mmsg<<"File "<<midasrunname.str().c_str()<<" Found";
              DANCE_Success("Main",mmsg.str());
              runname << midasrunname.str();
              input_params.SubRunNumber=-1;
              gz_queue.push(gz_in);
            }
          }
        }
      } //end .mid.gz
    } //end not single subrun
  } //end read midas 
 
  //Binary input that is not simulation
  else if(input_params.Read_Binary == 1 && input_params.Read_Simulation == 0) {
    
    stringstream binaryrunname;
    binaryrunname.str();
    binaryrunname << pathtodata << "/stage0_run_" << RunNum << ".bin";
    stringstream binarysubrunname;
    binarysubrunname.str();
    binarysubrunname << pathtodata << "/stage0_run_" << RunNum << "_" <<input_params.NumSubRun<< ".bin";
    mmsg.str("");
    mmsg<<"Checking for: "<<binarysubrunname.str();
    DANCE_Info("Main",mmsg.str());
    
    //Look for uncompressed .bin subrun files
    gz_in=gzopen(binarysubrunname.str().c_str(),"rb");
    
    //check to see if its open
    if(gz_in) {
      mmsg.str("");
      mmsg<<"File "<<binarysubrunname.str().c_str()<<" Found";
      DANCE_Success("Main",mmsg.str());
      runname << binarysubrunname.str();
      input_params.SubRunNumber=0;

      while(gz_in) {
        input_params.NumSubRun++;
        gz_queue.push(gz_in);
        binarysubrunname.str("");
        binarysubrunname << pathtodata << "/stage0_run_" << RunNum << "_" <<input_params.NumSubRun<< ".bin";
        gz_in=gzopen(binarysubrunname.str().c_str(),"rb");
      } 
    }
    else {
      //look for compressed .bin.gz subrun files
      binarysubrunname << ".gz";
      mmsg.str("");
      mmsg<<"Checking for: "<<binarysubrunname.str();
      DANCE_Info("Main",mmsg.str());
      gz_in=gzopen(binarysubrunname.str().c_str(),"rb");

      if(gz_in) {
	mmsg.str("");
	mmsg<<"File "<<binarysubrunname.str().c_str()<<" Found";
	DANCE_Success("Main",mmsg.str());
	runname << binarysubrunname.str();
        input_params.SubRunNumber=0;

        while(gz_in) {
          input_params.NumSubRun++;
          gz_queue.push(gz_in);
          binarysubrunname.str("");
          binarysubrunname << pathtodata << "/stage0_run_" << RunNum << "_" <<input_params.NumSubRun<< ".bin.gz";
          gz_in=gzopen(binarysubrunname.str().c_str(),"rb");
        } 

      }
      else {
        mmsg.str("");
        mmsg<<"Checking for: "<<binaryrunname.str().c_str();
        DANCE_Info("Main",mmsg.str());
        gz_in=gzopen(binaryrunname.str().c_str(),"rb");
        if (gz_in) {
          mmsg.str("");
          mmsg<<"File "<<binaryrunname.str().c_str()<<" Found";
          DANCE_Success("Main",mmsg.str());
          runname << binaryrunname.str();
          input_params.SubRunNumber=-1;
          gz_queue.push(gz_in);
        }
        else {
          binaryrunname << ".gz";
          mmsg.str("");
          mmsg<<"Checking for: "<<binaryrunname.str().c_str();
          gz_in=gzopen(binaryrunname.str().c_str(),"rb");
          if(gz_in) {
            DANCE_Info("Main",mmsg.str());
            DANCE_Success("Main",mmsg.str());
            runname << binaryrunname.str();
            input_params.SubRunNumber=-1;
            gz_queue.push(gz_in); 
          }
        }
      }
    }
  }

  //Simulated Data from GEANT4
  else if(input_params.Read_Binary == 0 && input_params.Read_Simulation == 1) {
    
    stringstream simulationrunname;
    simulationrunname.str();
    simulationrunname << STAGE0_SIM << "/"<< input_params.Simulation_File_Name <<".bin";

    mmsg.str("");
    mmsg<<"Checking for: "<<simulationrunname.str();
    DANCE_Info("Main",mmsg.str());
    
    //Look for uncompressed .bin files
    gz_in=gzopen(simulationrunname.str().c_str(),"rb");
    
    //check to see if its open
    if(gz_in) {
      mmsg.str("");
      mmsg<<"File "<<simulationrunname.str().c_str()<<" Found";
      DANCE_Success("Main",mmsg.str());
      runname << simulationrunname.str();
      gz_queue.push(gz_in);
      input_params.SubRunNumber=-1;
    }
    else {
      //look for compressed .bin.gz files
      simulationrunname << ".gz";
      mmsg.str("");
      mmsg<<"Checking for: "<<simulationrunname.str();
      DANCE_Info("Main",mmsg.str());
      gz_in=gzopen(simulationrunname.str().c_str(),"rb");
      if(gz_in) {

	mmsg.str("");
	mmsg<<"File "<<simulationrunname.str().c_str()<<" Found";
	DANCE_Success("Main",mmsg.str());
	runname << simulationrunname.str();
        gz_queue.push(gz_in);
        input_params.SubRunNumber=-1;
      }
    }
  }

  else {
    mmsg.str("");
    mmsg<<"Conflict Between Read_Binary (set to "<<input_params.Read_Binary<<") and Read_Simulation set to ("<<input_params.Read_Simulation<<"). Exiting";
    DANCE_Error("Main",mmsg.str());
    return -1;
  }
  
  if(gz_queue.empty()) {
    mmsg.str("");
    mmsg<<"File queue empty, this is bad.  Exiting.";
    DANCE_Error("Main",mmsg.str());
    return -1;
  }

  //Time profiling stuff
  struct timeval tv;  						// real time  
  double begin, end, time_elapsed;			// start,stop, elapsed time

  //start time
  gettimeofday(&tv,NULL); 
  begin=tv.tv_sec+(tv.tv_usec/1000000.0);

  //Initialize Unpacker
  func_ret = Initialize_Unpacker(input_params);
  if(func_ret) {
    return func_ret;
  }
  //initialize eventbuilder
  func_ret = Initialize_Eventbuilder(input_params);
  if(func_ret) {
    return func_ret;
  }
  //initialize calibrator
  func_ret = Initialize_Calibrator(input_params);
  if(func_ret) {
    return func_ret;
  }
  //initialize validator
  func_ret = Initialize_Validator(input_params);
  if(func_ret) {
    return func_ret;
  }
  //Initialize analyzer
  func_ret = Initialize_Analyzer(input_params);
  if(func_ret) {
    return func_ret;
  }

  //Launch the unpacker
  int events_analyzed=  Unpack_Data(gz_queue, begin, input_params, &analysis_params);

  mmsg.str("");
  mmsg<<"Analysis Complete. Analyzed: "<<events_analyzed<<" Entries";
  DANCE_Success("Main",mmsg.str());
  
  //Print the time elapsed for the entire process
  gettimeofday(&tv,NULL);  
  end=tv.tv_sec+(tv.tv_usec/1000000.0);
  time_elapsed = (double)(end-begin); ;

  mmsg.str("");
  mmsg<<"Time Elapsed: "<<time_elapsed<<" Seconds";
  DANCE_Info("Main",mmsg.str());

}
