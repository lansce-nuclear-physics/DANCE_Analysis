
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
//*  validator.h            *// 
//*  Last Edit: 04/15/19    *//  
//***************************//

#ifndef VALIDATOR_H
#define VALIDATOR_H

//File Includes
#include "structures.h"

#include "TCutG.h"
#include "TFile.h"

//Function Prototypes
int Read_PI_Gates();
int Write_PI_Gates(TFile *fout);

int Check_ULD(DEVT_BANK *devt_bank);  //Upper level discriminator
int Check_Threshold(DEVT_BANK *devt_bank);  //Threshold
int Check_Alpha(DEVT_BANK *devt_bank);
int Check_Gamma(DEVT_BANK *devt_bank);
int Check_Crystal_Blocking(DEVT_BANK *devt_bank, Analysis_Parameters *analysis_params, Input_Parameters input_params);
int Check_Retrigger(DEVT_BANK *devt_bank, Analysis_Parameters *analysis_params);  //Retrigger Gate
int Check_Threshold(DEVT_BANK *devt_bank, Input_Parameters input_params);  //Energy Threshold
int Initialize_Validator(Input_Parameters input_params);

#endif
