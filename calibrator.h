//***************************//
//*  Christopher J. Prokop  *//
//*  cprokop@lanl.gov       *//
//*  calibrator.h           *// 
//*  Last Edit: 05/08/18    *//  
//***************************//

#ifndef CALIBRATOR_H
#define CALIBRATOR_H

//File Includes
#include "structures.h"

//Function Prototypes
int Read_Energy_Calibrations( Input_Parameters input_params);
int Calibrate_DANCE(DEVT_BANK *devt_bank);
int Initialize_Calibrator(Input_Parameters input_params);

#endif
