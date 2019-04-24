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
