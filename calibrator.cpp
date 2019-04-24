//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  calibrator.cpp         *// 
//*  Last Edit: 04/22/19    *//  
//***************************//

#include "calibrator.h"
#include "global.h"
#include "message.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "TRandom.h"

using namespace std;

//energy calibrations
double slow_offset[200];
double slow_slope[200];
double slow_quad[200];
double fast_slope[200];
double fast_offset[200];
double fast_quad[200];

stringstream cmsg;


int Read_Energy_Calibrations(Input_Parameters input_params) {
  
  int id=0;
  double temp[6] = {0,0,0,0,0,0};

  for(int eye=0; eye<200; eye++) {
    slow_offset[eye]=0;
    slow_slope[eye]=1;
    slow_quad[eye]=0;
    fast_offset[eye]=0;
    fast_slope[eye]=1;  
    fast_quad[eye]=0;  
  }

  DANCE_Info("Calibrator","Reading Energy Calibrations");

  if(input_params.Read_Simulation == 0) {


    stringstream cal_name;
    cal_name.str();
    cal_name << CALIB_DIR << "/param_out_" << input_params.RunNumber << ".txt";
    
    ifstream encal;
    encal.open(cal_name.str().c_str());
  
    bool encalfail=false;
  
    if(input_params.Analysis_Stage==1) {
      if(encal.is_open()) {
	while(!encal.eof()) {
	  //encal>>id>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4]>>temp[5];
	  encal>>id>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4];
	  slow_offset[id]=temp[0];
	  slow_slope[id]=temp[1];
	  slow_quad[id]=temp[2];
	  fast_offset[id]=temp[3];
	  fast_slope[id]=temp[4];
	  //fast_quad[id]=temp[5];

	  if(encal.eof()) {
	    break;
	  }
	}  
      
	DANCE_Success("Calibrator","Read Energy Calibrations");
	return 0;
      }
      else {

	cmsg.str("");
	cmsg<<"File: "<<cal_name.str()<<" NOT found...";
	
	DANCE_Error("Calibrator",cmsg.str());

	encalfail=true;
      }
    }
    if(input_params.Analysis_Stage==0 || encalfail==true)

      DANCE_Info("Calibrator","Looking for calib_ideal.dat");
  
    stringstream idealcal_name;
    idealcal_name.str();
    idealcal_name << CALIB_DIR << "/calib_ideal.dat";
  
    ifstream idealcal;
    idealcal.open(idealcal_name.str().c_str());
  
    if(idealcal.is_open()) {
      while(!idealcal.eof()) {
	
	//	idealcal>>id>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4]>>temp[5];
	idealcal>>id>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4];
	slow_offset[id]=temp[0];
	slow_slope[id]=temp[1];
	slow_quad[id]=temp[2];
	fast_offset[id]=temp[3];
	fast_slope[id]=temp[4];
	//	fast_quad[id]=temp[5];

	cout<<id<<"  "<<slow_slope[id]<<"  "<<fast_slope[id]<<endl;
      
	if(idealcal.eof()) {
	  break;
	}
      }

      cmsg.str("");
      cmsg<<"Opened "<<idealcal_name.str();
      DANCE_Success("Calibrator", cmsg.str());

      return 0;
    }
    else {
      DANCE_Error("Calibrator","Can NOT Open Energy Calibration File");
      return-1;
    }
  }  
  DANCE_Info("Calibrator", "Set all calibrations off for simulations");
  return 0;

}


int Calibrate_DANCE(DEVT_BANK *devt_bank) {

  //Apply Energy Calibration
#ifdef Calibrator_Verbose
  cout<<"Calibrator:  ID: "<<devt_bank->ID<<endl;
#endif

  //DANCE Ball
  if(devt_bank->ID < 162) {
    
    double temp_slow = devt_bank->Islow + gRandom->Uniform(0,1);
    double temp_fast = devt_bank->Ifast + gRandom->Uniform(0,1);
    
    devt_bank->Eslow = 0.001*(temp_slow*temp_slow*slow_quad[devt_bank->ID] +
			      temp_slow*slow_slope[devt_bank->ID] +
			      slow_offset[devt_bank->ID]);
    
    devt_bank->Efast = 0.001*(temp_fast*temp_fast*fast_quad[devt_bank->ID] +
			      temp_fast*fast_slope[devt_bank->ID] +
			      fast_offset[devt_bank->ID]);
  
#ifdef Calibrator_Verbose
    cout<<"Calibrator: fast: "<< devt_bank->Efast<<"  slow: "<< devt_bank->Eslow<<endl;
#endif
    
  }
  
  return 0;

}


int Initialize_Calibrator(Input_Parameters input_params) {

  DANCE_Init("Calibrator","Initializing");

  Read_Energy_Calibrations(input_params);

  DANCE_Success("Calibrator","Initialized");

  return 0;
  
}
