
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
//*  validator.cpp          *// 
//*  Last Edit: 04/15/19    *//  
//***************************//

#include "validator.h"
#include "message.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std; 

stringstream msg;

//PI Cuts
TCutG *Gamma_Gate;
TCutG *Alpha_Gate;
TCutG *Retrigger_Gate;

//PI Gates for Alphas, Gammas, and Retrigger
int Read_PI_Gates() {

  DANCE_Info("Validator","Reading PI Gates");

  char name[200];
  int Nalphacut=0;
  int Ngammacut=0;
  int Nretriggercut=0;

  double x_alphacut[200];
  double y_alphacut[200];
  double x_gammacut[200];
  double y_gammacut[200];
  double x_retriggercut[200];
  double y_retriggercut[200];


  sprintf(name,"Gates/%s",ALPHAGATE);
  ifstream alphacutin(name);
  if(alphacutin.is_open()) {
    
    DANCE_Info("Validator","Reading Alpha PI Cut");

    while(!alphacutin.eof()){
      alphacutin  >> x_alphacut[Nalphacut] >> y_alphacut[Nalphacut];
      Nalphacut++;
      if(alphacutin.eof()) break;
    }
    // cout<<Nalphacut<<endl;
    alphacutin.close();
    
    
    msg.str("");
    msg<<"Alpha PI Cut " << name<< " Loaded";
    DANCE_Success("Validator",msg.str());
  }
  else {
    
    msg.str("");
    msg<<"Failed to Load Alpha PI Cut " << name;
    DANCE_Error("Validator",msg.str());
    return -1;
  }
  
  
  sprintf(name,"Gates/%s",GAMMAGATE);
  ifstream gammacutin(name);
  if(gammacutin.is_open()) {

    DANCE_Info("Validator","Reading Gamma PI Cut");

    while(!gammacutin.eof()){
      gammacutin  >> x_gammacut[Ngammacut] >> y_gammacut[Ngammacut];
      Ngammacut++;
      if(gammacutin.eof()) break;
    }
    // cout<<Ngammacut<<endl;
    gammacutin.close();
   
    msg.str("");
    msg<<"Gamma PI Cut " << name<< " Loaded";
    DANCE_Success("Validator",msg.str());
  }
  else {
   
    msg.str("");
    msg<<"Failed to Load Gamma PI Cut " << name;
    DANCE_Error("Validator",msg.str());
    return -1;
  }


  sprintf(name,"Gates/%s",RETRIGGERGATE);
  ifstream retriggercutin(name);
  if(retriggercutin.is_open()) {

    DANCE_Info("Validator","Reading Retrigger PI Cut");

    while(!retriggercutin.eof()){
      retriggercutin  >> x_retriggercut[Nretriggercut] >> y_retriggercut[Nretriggercut];
      Nretriggercut++;
      if(retriggercutin.eof()) break;
    }
    // cout<<Nretriggercut<<endl;
    retriggercutin.close();
    msg.str("");
    msg<<"Retrigger PI Cut " << name<< " Loaded";
    DANCE_Success("Validator",msg.str());
  }
  else {
    
    msg.str("");
    msg<<"Failed to Load Retrigger PI Cut " << name;
    DANCE_Error("Validator",msg.str());
    return -1;
  }



  Alpha_Gate=new TCutG("Alpha_Gate",Nalphacut,x_alphacut,y_alphacut);
  // Alpha_Gate->Print();	
  Gamma_Gate=new TCutG("Gamma_Gate",Ngammacut,x_gammacut,y_gammacut);
  //  Gamma_Gate->Print();
  Retrigger_Gate=new TCutG("Retrigger_Gate",Nretriggercut,x_retriggercut,y_retriggercut);


  DANCE_Success("Validator","Read PI Gates");

  return 0;
}

int Write_PI_Gates(TFile *fout){
  
  DANCE_Info("Validator","Writing PI Gates");

  fout->cd();
  
  Gamma_Gate->Write();
  Alpha_Gate->Write();
  Retrigger_Gate->Write();

  DANCE_Success("Validator","Wrote PI Gates");
  return 0;
}


//Upper level discriminator
int Check_ULD(DEVT_BANK *devt_bank) {
  
  //Upper Level Discriminator to reject junk that appears huge because of invalid conversion between signed and unsigned integers.  Current max value is 65535

  if(devt_bank->Islow > 62500 || devt_bank->Ifast > 31000) {
    devt_bank->Valid=0;
    devt_bank->InvalidReason += 32;
#ifdef Validator_Verbose
    cout<<RED<<"Validator: Event Invalid from ULD"<<RESET<<endl;
#endif
  }
  return 0;
}


//Energy Threshold
int Check_Threshold(DEVT_BANK *devt_bank, Input_Parameters input_params) {
  
  if(input_params.HAVE_Threshold==1) {
    if(devt_bank->Eslow < input_params.Energy_Threshold) {
      devt_bank->Valid=0;
      devt_bank->InvalidReason |= Invalid::Threshold;
#ifdef Validator_Verbose
      cout<<RED<<"Validator: Invalid from Threshold"<<RESET<<endl;
#endif
    }
  }
  return 0;
}


//Crystal Blocking
int Check_Crystal_Blocking(DEVT_BANK *devt_bank, Analysis_Parameters *analysis_params, Input_Parameters input_params) {
  
  int ID = devt_bank->ID;
  
  double timediff = devt_bank->timestamp - analysis_params->last_timestamp[ID];

  //Blocking Time to avoid retrigger problems
  if(timediff < input_params.Crystal_Blocking_Time) {
    devt_bank->Valid=0;
    devt_bank->InvalidReason |= Invalid::CrystalBlocking;
#ifdef Validator_Verbose
    cout<<RED<<"Validator: Invalid from Blocking Time"<<RESET<<endl;
#endif
  }
  return 0;
}


//PI Gate to avoid retrigger problems
int Check_Retrigger(DEVT_BANK *devt_bank, Analysis_Parameters *analysis_params) {
  
  //Crystal ID
  int ID = devt_bank->ID;
  
  //Time between this crystal hit and the last crystal hit
  double timediff = devt_bank->timestamp - analysis_params->last_timestamp[ID];

  //Ratio of this crystal hit and the last crystal hit
  double slowratio = devt_bank->Islow / analysis_params->last_Islow[ID];
  double fastratio = devt_bank->Efast / analysis_params->last_Efast[ID];
  
  //  cout<<"timediff: "<<timediff<<" ratio: "<<ratio<<"  "<<Retrigger_Gate<<endl;

  if(Retrigger_Gate->IsInside(timediff,slowratio)) {
    devt_bank->Valid=0;
    devt_bank->InvalidReason += 8;
#ifdef Validator_Verbose
    cout<<RED<<"Validator: Invalid from Retrigger"<<RESET<<endl;
#endif
  }

  if(Retrigger_Gate->IsInside(timediff,fastratio)) {
    devt_bank->Valid=0;
    devt_bank->InvalidReason |= Invalid::RetriggerFast;
#ifdef Validator_Verbose
    cout<<RED<<"Validator: Invalid from fast Retrigger"<<RESET<<endl;
#endif
  }



  return 0;
}


int Check_Alpha(DEVT_BANK *devt_bank) {
  
  if(Alpha_Gate->IsInside(devt_bank->Eslow, devt_bank->Efast)) {
    devt_bank->IsGamma = 0;
    devt_bank->IsAlpha = 1;    
#ifdef Validator_Verbose
    cout<<RED<<"Validator: Alpha"<<RESET<<endl;
#endif
  }
  return 0;
}


int Check_Gamma(DEVT_BANK *devt_bank) {
  
  if(Gamma_Gate->IsInside(devt_bank->Eslow, devt_bank->Efast)) {
    devt_bank->IsGamma = 1;
    devt_bank->IsAlpha = 0;    
#ifdef Validator_Verbose
    cout<<GREEN<<"Validator: Gamma"<<RESET<<endl;
#endif
  }
  return 0;
}
  

int Initialize_Validator(Input_Parameters input_params) {

  DANCE_Init("Validator","Initializing");

  int func_ret=0;
  func_ret = Read_PI_Gates();

  if(func_ret==0) {
    DANCE_Success("Validator","Initialized");
  }
  else {
    DANCE_Error("Validator","Initialization Failed. Exiting!");
  }
  return func_ret;

}
